#include "pgbam/bam_io.hpp"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <fstream>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <sstream>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <vector>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "pgbam/error.hpp"
#include "pgbam/fingerprint.hpp"
#include "pgbam/gaf.hpp"
#include "pgbam/gbwt_backend.hpp"
#include "pgbam/sidecar.hpp"

namespace pgbam {
namespace {

// Natural-sort comparison matching htslib/samtools strnum_cmp used by
// `samtools sort -n`.  Digit runs are compared numerically so that e.g.
// "read/9/ccs" < "read/10/ccs" rather than the other way around.
bool qname_less(std::string_view a, std::string_view b) {
  std::size_t i = 0, j = 0;
  while (i < a.size() && j < b.size()) {
    if (std::isdigit(static_cast<unsigned char>(a[i])) &&
        std::isdigit(static_cast<unsigned char>(b[j]))) {
      while (i < a.size() && a[i] == '0') ++i;
      while (j < b.size() && b[j] == '0') ++j;
      std::size_t ia = i, jb = j;
      while (ia < a.size() && std::isdigit(static_cast<unsigned char>(a[ia]))) ++ia;
      while (jb < b.size() && std::isdigit(static_cast<unsigned char>(b[jb]))) ++jb;
      std::size_t alen = ia - i, blen = jb - j;
      if (alen != blen) return alen < blen;
      while (i < ia) {
        if (a[i] != b[j]) return a[i] < b[j];
        ++i; ++j;
      }
      i = ia; j = jb;
    } else {
      if (a[i] != b[j]) return static_cast<unsigned char>(a[i]) < static_cast<unsigned char>(b[j]);
      ++i; ++j;
    }
  }
  return i == a.size() && j < b.size();
}

using SamFilePtr = std::unique_ptr<samFile, decltype(&hts_close)>;
using SamHeaderPtr = std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)>;

struct VectorHash {
  std::size_t operator()(const std::vector<std::uint64_t>& values) const noexcept {
    std::size_t seed = 0;
    for (std::uint64_t value : values) {
      seed ^= std::hash<std::uint64_t>{}(value) + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
    }
    return seed;
  }
};

struct WorkItem {
  std::size_t ordinal = 0;
  std::vector<bam1_t*> records;
  std::vector<std::vector<OrientedNode>> walks;
};

struct ResultItem {
  std::size_t ordinal = 0;
  std::vector<bam1_t*> records;
  std::vector<SubpathResult> subpaths;
};

class BamRecord {
public:
  BamRecord() : record_(bam_init1()) {
    if (record_ == nullptr) {
      throw Error("unable to allocate BAM record");
    }
  }

  ~BamRecord() {
    if (record_ != nullptr) {
      bam_destroy1(record_);
    }
  }

  bam1_t* get() const { return record_; }

private:
  bam1_t* record_ = nullptr;
};

class SidecarInterner {
public:
  explicit SidecarInterner(SidecarWriter& writer) : writer_(writer) {}

  std::uint32_t intern(const std::vector<std::uint64_t>& thread_ids) {
    auto found = set_to_id_.find(thread_ids);
    if (found != set_to_id_.end()) {
      return found->second;
    }

    const std::uint32_t assigned = next_set_id_++;
    set_to_id_.emplace(thread_ids, assigned);
    writer_.write_set(SidecarSetRecord{assigned, thread_ids});
    return assigned;
  }

private:
  SidecarWriter& writer_;
  std::unordered_map<std::vector<std::uint64_t>, std::uint32_t, VectorHash> set_to_id_;
  std::uint32_t next_set_id_ = 0;
};

void destroy_records(std::vector<bam1_t*>& records) {
  for (bam1_t*& record : records) {
    if (record != nullptr) {
      bam_destroy1(record);
      record = nullptr;
    }
  }
  records.clear();
}

bam1_t* duplicate_bam_record(const bam1_t* record) {
  bam1_t* duplicate = bam_dup1(record);
  if (duplicate == nullptr) {
    throw Error("unable to duplicate BAM record for " + std::string(bam_get_qname(record)));
  }
  return duplicate;
}

bool subpath_less(const SubpathResult& left, const SubpathResult& right) {
  if (left.begin_offset != right.begin_offset) {
    return left.begin_offset < right.begin_offset;
  }
  if (left.end_offset != right.end_offset) {
    return left.end_offset < right.end_offset;
  }
  return std::lexicographical_compare(left.thread_ids.begin(), left.thread_ids.end(),
                                      right.thread_ids.begin(), right.thread_ids.end());
}

bool subpath_equal(const SubpathResult& left, const SubpathResult& right) {
  return left.begin_offset == right.begin_offset &&
         left.end_offset == right.end_offset &&
         left.thread_ids == right.thread_ids;
}

void deduplicate_subpaths(std::vector<SubpathResult>& subpaths) {
  std::sort(subpaths.begin(), subpaths.end(), subpath_less);
  subpaths.erase(std::unique(subpaths.begin(), subpaths.end(), subpath_equal), subpaths.end());
}

std::string build_command_line(const AnnotateOptions& options) {
  std::ostringstream out;
  out << "pgbam annotate --bam " << options.bam_path
      << " --gaf " << options.gaf_path
      << (options.use_gbz() ? " --gbz " + options.gbz_path : " --gbwt " + options.gbwt_path)
      << " --out-bam " << options.out_bam_path
      << " --out-sets " << options.out_sets_path
      << " --threads " << options.threads;
  if (options.use_r_index()) {
    out << " --r-index " << options.r_index_path;
  }
  if (options.primary_only) {
    out << " --primary-only";
  }
  return out.str();
}

void set_or_remove_array_tag(bam1_t* record, const char tag[2], const std::vector<std::uint32_t>& values) {
  if (values.empty()) {
    if (uint8_t* existing = bam_aux_get(record, tag)) {
      bam_aux_del(record, existing);
    }
    return;
  }

  if (bam_aux_update_array(record, tag, 'I', static_cast<std::uint32_t>(values.size()),
                           const_cast<std::uint32_t*>(values.data())) != 0) {
    throw Error("failed to update BAM aux tag " + std::string(tag, 2));
  }
}

sam_hdr_t* prepare_header(const sam_hdr_t* input_header, const AnnotateOptions& options) {
  sam_hdr_t* header = sam_hdr_dup(input_header);
  if (header == nullptr) {
    throw Error("unable to duplicate BAM header");
  }

  const std::string command_line = build_command_line(options);
  if (sam_hdr_add_pg(header, "pgbam", "VN", "0.1.0", "CL", command_line.c_str(), nullptr) != 0) {
    sam_hdr_destroy(header);
    throw Error("unable to add @PG line to BAM header");
  }

  const std::string comment = "@CO\tpgbam tags: hs:B:I=set ids, hb:B:I=begin node offsets, he:B:I=end-exclusive node offsets\n";
  if (sam_hdr_add_lines(header, comment.c_str(), comment.size()) != 0) {
    sam_hdr_destroy(header);
    throw Error("unable to add @CO line to BAM header");
  }

  return header;
}

SamFilePtr open_sam_file(const std::string& path, const char* mode) {
  samFile* handle = sam_open(path.c_str(), mode);
  if (handle == nullptr) {
    throw Error("cannot open BAM file: " + path);
  }
  return SamFilePtr(handle, hts_close);
}

SamHeaderPtr read_header(samFile* input, const std::string& path) {
  sam_hdr_t* header = sam_hdr_read(input);
  if (header == nullptr) {
    throw Error("cannot read BAM header: " + path);
  }
  return SamHeaderPtr(header, sam_hdr_destroy);
}

int read_bam_record(samFile* input, sam_hdr_t* header, bam1_t* record, const std::string& path) {
  const int status = sam_read1(input, header, record);
  if (status < -1) {
    throw Error("failed to read BAM record from " + path);
  }
  return status;
}

std::string make_gaf_qname_order_error(const std::string& path,
                                       const std::string& previous_qname,
                                       std::size_t previous_line,
                                       const std::string_view current_qname,
                                       std::size_t current_line) {
  std::ostringstream out;
  out << "GAF must be sorted by qname in column 1: " << path
      << " (line " << current_line << " qname '" << current_qname
      << "' appears after line " << previous_line << " qname '" << previous_qname
      << "'; use: sort -k1,1 input.gaf > sorted.gaf)";
  return out.str();
}

std::string make_bam_qname_order_error(const std::string& path,
                                       const std::string& previous_qname,
                                       std::size_t previous_record,
                                       const std::string_view current_qname,
                                       std::size_t current_record) {
  std::ostringstream out;
  out << "BAM must be qname-sorted: " << path
      << " (record " << current_record << " qname '" << current_qname
      << "' appears after record " << previous_record << " qname '" << previous_qname
      << "'; use: samtools sort -n -o sorted.bam input.bam)";
  return out.str();
}

std::string make_gaf_mapq_order_error(const std::string& path,
                                      const std::string& qname,
                                      std::size_t previous_line,
                                      std::uint32_t previous_mapq,
                                      std::size_t current_line,
                                      std::uint32_t current_mapq) {
  std::ostringstream out;
  out << "GAF must be sorted by qname and then non-increasing MAPQ when --primary-only is used: " << path
      << " (line " << current_line << " qname '" << qname << "' mapq " << current_mapq
      << " appears after line " << previous_line << " mapq " << previous_mapq
      << "; use: sort -k1,1 -k12,12nr input.gaf > sorted.gaf)";
  return out.str();
}

void ensure_gaf_order(const AnnotateOptions& options) {
  const std::string& path = options.gaf_path;
  std::ifstream input(path);
  if (!input) {
    throw Error("cannot open GAF file: " + path);
  }

  std::string previous_qname;
  std::string line;
  std::size_t previous_line = 0;
  std::size_t line_number = 0;
  std::uint32_t previous_mapq = 0;
  while (std::getline(input, line)) {
    ++line_number;
    auto record = parse_gaf_line(line);
    if (!record) {
      continue;
    }
    if (!previous_qname.empty() && qname_less(record->qname, previous_qname)) {
      throw Error(make_gaf_qname_order_error(path, previous_qname, previous_line, record->qname, line_number));
    }
    if (options.primary_only && previous_qname == record->qname && previous_mapq < record->mapq) {
      throw Error(make_gaf_mapq_order_error(path, record->qname, previous_line, previous_mapq, line_number, record->mapq));
    }
    previous_qname = record->qname;
    previous_line = line_number;
    previous_mapq = record->mapq;
  }
}

void ensure_bam_qname_sorted(const std::string& path, std::size_t threads) {
  auto input = open_sam_file(path, "rb");
  if (threads > 1) {
    hts_set_threads(input.get(), static_cast<int>(threads));
  }
  auto header = read_header(input.get(), path);
  BamRecord record;
  std::string previous;
  std::size_t previous_record = 0;
  std::size_t record_number = 0;
  while (true) {
    const int status = read_bam_record(input.get(), header.get(), record.get(), path);
    if (status == -1) {
      break;
    }
    ++record_number;
    const std::string_view qname = bam_get_qname(record.get());
    if (!previous.empty() && qname_less(qname, previous)) {
      throw Error(make_bam_qname_order_error(path, previous, previous_record, qname, record_number));
    }
    previous.assign(qname.data(), qname.size());
    previous_record = record_number;
  }
}

void ensure_qname_sorted_inputs(const AnnotateOptions& options) {
  ensure_bam_qname_sorted(options.bam_path, options.threads);
  ensure_gaf_order(options);
}

std::optional<GafRecord> read_next_gaf_record(std::istream& input) {
  std::string line;
  while (std::getline(input, line)) {
    if (auto record = parse_gaf_line(line)) {
      return record;
    }
  }
  return std::nullopt;
}

template<class PrepareItem>
void annotate_stream(samFile* input,
                     sam_hdr_t* input_header,
                     const std::string& input_path,
                     samFile* output,
                     sam_hdr_t* output_header,
                     const GraphIndex& graph,
                     std::size_t threads,
                     SidecarInterner& interner,
                     PrepareItem prepare_item) {
  if (threads > 1) {
    hts_set_threads(input, static_cast<int>(threads));
    hts_set_threads(output, static_cast<int>(threads));
  }

  if (sam_hdr_write(output, output_header) != 0) {
    throw Error("cannot write BAM header");
  }

  std::deque<WorkItem> work_queue;
  std::deque<ResultItem> result_queue;
  std::mutex work_mutex;
  std::mutex result_mutex;
  std::mutex error_mutex;
  std::condition_variable work_ready;
  std::condition_variable work_space;
  std::condition_variable result_ready;
  bool input_done = false;
  bool workers_done = false;
  std::atomic<bool> has_failure{false};
  std::exception_ptr failure;
  const std::size_t max_work_queue = std::max<std::size_t>(1, threads * 8);

  auto fail = [&](std::exception_ptr error) {
    std::lock_guard<std::mutex> lock(error_mutex);
    if (!failure) {
      failure = error;
      has_failure.store(true);
    }
    work_ready.notify_all();
    work_space.notify_all();
    result_ready.notify_all();
  };

  auto destroy_queued_records = [&]() {
    {
      std::lock_guard<std::mutex> lock(work_mutex);
      for (WorkItem& item : work_queue) {
        destroy_records(item.records);
      }
      work_queue.clear();
    }
    {
      std::lock_guard<std::mutex> lock(result_mutex);
      for (ResultItem& item : result_queue) {
        destroy_records(item.records);
      }
      result_queue.clear();
    }
  };

  std::vector<std::thread> workers;
  workers.reserve(threads);
  for (std::size_t thread_id = 0; thread_id < threads; ++thread_id) {
    workers.emplace_back([&]() {
      try {
        while (true) {
          WorkItem item;
          {
            std::unique_lock<std::mutex> lock(work_mutex);
            work_ready.wait(lock, [&]() { return has_failure.load() || !work_queue.empty() || input_done; });
            if (has_failure.load()) {
              return;
            }
            if (work_queue.empty()) {
              if (input_done) {
                return;
              }
              continue;
            }
            item = std::move(work_queue.front());
            work_queue.pop_front();
          }
          work_space.notify_one();

          ResultItem result;
          result.ordinal = item.ordinal;
          result.records = std::move(item.records);
          for (const std::vector<OrientedNode>& walk : item.walks) {
            std::vector<SubpathResult> walk_subpaths = graph.find_subpaths(walk);
            result.subpaths.insert(result.subpaths.end(),
                                   std::make_move_iterator(walk_subpaths.begin()),
                                   std::make_move_iterator(walk_subpaths.end()));
          }
          if (result.subpaths.size() > 1) {
            deduplicate_subpaths(result.subpaths);
          }

          {
            std::lock_guard<std::mutex> lock(result_mutex);
            result_queue.push_back(std::move(result));
          }
          result_ready.notify_one();
        }
      } catch (...) {
        fail(std::current_exception());
      }
    });
  }

  std::thread writer([&]() {
    try {
      std::map<std::size_t, ResultItem> pending;
      std::size_t next_ordinal = 0;

      while (true) {
        {
          std::unique_lock<std::mutex> lock(result_mutex);
          result_ready.wait(lock, [&]() { return has_failure.load() || !result_queue.empty() || workers_done; });
          while (!result_queue.empty()) {
            ResultItem item = std::move(result_queue.front());
            result_queue.pop_front();
            pending.emplace(item.ordinal, std::move(item));
          }
        }

        while (true) {
          auto next = pending.find(next_ordinal);
          if (next == pending.end()) {
            break;
          }

          ResultItem item = std::move(next->second);
          pending.erase(next);

          ReadAnnotation annotation;
          annotation.hs.reserve(item.subpaths.size());
          annotation.hb.reserve(item.subpaths.size());
          annotation.he.reserve(item.subpaths.size());

          for (const SubpathResult& subpath : item.subpaths) {
            if (subpath.thread_ids.empty()) {
              continue;
            }
            annotation.hs.push_back(interner.intern(subpath.thread_ids));
            annotation.hb.push_back(subpath.begin_offset);
            annotation.he.push_back(subpath.end_offset);
          }

          for (bam1_t* record : item.records) {
            set_or_remove_array_tag(record, "hs", annotation.hs);
            set_or_remove_array_tag(record, "hb", annotation.hb);
            set_or_remove_array_tag(record, "he", annotation.he);

            if (sam_write1(output, output_header, record) < 0) {
              throw Error("failed to write BAM record");
            }
          }

          destroy_records(item.records);
          ++next_ordinal;
        }

        if (has_failure.load()) {
          for (auto& [_, pending_item] : pending) {
            destroy_records(pending_item.records);
          }
          return;
        }

        if (workers_done && pending.empty()) {
          return;
        }
      }
    } catch (...) {
      fail(std::current_exception());
    }
  });

  try {
    BamRecord record;
    std::size_t ordinal = 0;
    std::size_t record_number = 0;
    bool have_record = false;
    auto read_next_record = [&]() {
      const int read_status = read_bam_record(input, input_header, record.get(), input_path);
      if (read_status == -1) {
        have_record = false;
        return;
      }
      have_record = true;
      ++record_number;
    };

    read_next_record();
    while (have_record) {
      WorkItem item;
      item.ordinal = ordinal++;

      const std::string qname = bam_get_qname(record.get());
      const std::size_t first_record_number = record_number;
      while (true) {
        item.records.push_back(duplicate_bam_record(record.get()));
        read_next_record();
        if (!have_record || bam_get_qname(record.get()) != qname) {
          break;
        }
      }

      try {
        prepare_item(qname, first_record_number, item);
      } catch (...) {
        destroy_records(item.records);
        throw;
      }

      {
        std::unique_lock<std::mutex> lock(work_mutex);
        work_space.wait(lock, [&]() { return has_failure.load() || work_queue.size() < max_work_queue; });
        if (has_failure.load()) {
          destroy_records(item.records);
          std::lock_guard<std::mutex> error_lock(error_mutex);
          if (failure) {
            std::rethrow_exception(failure);
          }
          throw Error("annotation failed");
        }
        work_queue.push_back(std::move(item));
      }
      work_ready.notify_one();
    }
  } catch (...) {
    fail(std::current_exception());
  }

  {
    std::lock_guard<std::mutex> lock(work_mutex);
    input_done = true;
  }
  work_ready.notify_all();

  for (std::thread& worker : workers) {
    if (worker.joinable()) {
      worker.join();
    }
  }

  {
    std::lock_guard<std::mutex> lock(result_mutex);
    workers_done = true;
  }
  result_ready.notify_all();

  if (writer.joinable()) {
    writer.join();
  }

  {
    std::lock_guard<std::mutex> lock(error_mutex);
    if (failure) {
      destroy_queued_records();
      std::rethrow_exception(failure);
    }
  }
}

void run_sorted_join(const AnnotateOptions& options,
                     const GraphIndex& graph,
                     SidecarInterner& interner) {
  auto input = open_sam_file(options.bam_path, "rb");
  auto input_header = read_header(input.get(), options.bam_path);

  auto output = open_sam_file(options.out_bam_path, "wb");
  SamHeaderPtr output_header(prepare_header(input_header.get(), options), sam_hdr_destroy);

  std::ifstream gaf_input(options.gaf_path);
  if (!gaf_input) {
    throw Error("cannot open GAF file: " + options.gaf_path);
  }

  std::string previous_bam_qname;
  std::size_t previous_bam_record = 0;
  std::size_t gaf_record_number = 0;
  std::string previous_gaf_qname;
  std::size_t previous_gaf_record = 0;
  std::optional<GafRecord> current_gaf;

  auto advance_gaf = [&]() {
    current_gaf = read_next_gaf_record(gaf_input);
    if (!current_gaf) {
      return;
    }
    ++gaf_record_number;
    if (!previous_gaf_qname.empty() && qname_less(current_gaf->qname, previous_gaf_qname)) {
      throw Error(make_gaf_qname_order_error(options.gaf_path,
                                             previous_gaf_qname,
                                             previous_gaf_record,
                                             current_gaf->qname,
                                             gaf_record_number));
    }
    previous_gaf_qname = current_gaf->qname;
    previous_gaf_record = gaf_record_number;
  };

  auto collect_gaf_walks_for_current_qname = [&]() {
    std::vector<std::vector<OrientedNode>> walks;
    if (!current_gaf) {
      return walks;
    }

    const std::string qname = current_gaf->qname;
    walks.push_back(current_gaf->nodes);
    std::uint32_t previous_mapq = current_gaf->mapq;
    std::size_t previous_mapq_record = gaf_record_number;
    advance_gaf();

    while (current_gaf && current_gaf->qname == qname) {
      if (options.primary_only && previous_mapq < current_gaf->mapq) {
        throw Error(make_gaf_mapq_order_error(options.gaf_path,
                                              qname,
                                              previous_mapq_record,
                                              previous_mapq,
                                              gaf_record_number,
                                              current_gaf->mapq));
      }
      if (!options.primary_only) {
        walks.push_back(current_gaf->nodes);
      }
      previous_mapq = current_gaf->mapq;
      previous_mapq_record = gaf_record_number;
      advance_gaf();
    }

    return walks;
  };

  advance_gaf();
  annotate_stream(input.get(), input_header.get(), options.bam_path, output.get(), output_header.get(),
                  graph, options.threads, interner,
                  [&](std::string_view qname, std::size_t first_record_number, WorkItem& item) {
                    if (!previous_bam_qname.empty() && qname_less(qname, previous_bam_qname)) {
                      throw Error(make_bam_qname_order_error(options.bam_path,
                                                             previous_bam_qname,
                                                             previous_bam_record,
                                                             qname,
                                                             first_record_number));
                    }
                    previous_bam_qname.assign(qname.data(), qname.size());
                    previous_bam_record = first_record_number;

                    while (current_gaf && qname_less(current_gaf->qname, qname)) {
                      advance_gaf();
                    }

                    if (current_gaf && current_gaf->qname == qname) {
                      item.walks = collect_gaf_walks_for_current_qname();
                    }
                  });
}

}  // namespace

int run_annotate(const AnnotateOptions& options) {
  ensure_qname_sorted_inputs(options);

  const std::unique_ptr<GraphIndex> graph = make_graph_index(GraphIndex::Config{
    options.gbz_path,
    options.gbwt_path,
    options.r_index_path,
  });

  const std::string fingerprint = sha256_file(options.use_gbz() ? options.gbz_path : options.gbwt_path);
  SidecarWriter sidecar(options.out_sets_path);
  sidecar.write_header(SidecarHeader{1, fingerprint, options.use_r_index()});
  SidecarInterner interner(sidecar);

  run_sorted_join(options, *graph, interner);
  return 0;
}

}  // namespace pgbam
