#include "pgbam/bam_io.hpp"

#include <algorithm>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <fstream>
#include <map>
#include <memory>
#include <mutex>
#include <atomic>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <utility>
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
  bam1_t* record = nullptr;
  const std::vector<OrientedNode>* walk = nullptr;
};

struct ResultItem {
  std::size_t ordinal = 0;
  bam1_t* record = nullptr;
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

std::string build_command_line(const AnnotateOptions& options) {
  std::ostringstream out;
  out << "pgbam annotate --bam " << options.bam_path
      << " --gaf " << options.gaf_path
      << (options.use_gbz() ? " --gbz " + options.gbz_path : " --gbwt " + options.gbwt_path)
      << " --out-bam " << options.out_bam_path
      << " --out-sets " << options.out_sets_path;
  if (options.use_r_index()) {
    out << " --r-index " << options.r_index_path;
  }
  out << " --threads " << options.threads;
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

std::unordered_map<std::string, std::vector<OrientedNode>> load_gaf_lookup(const std::string& path) {
  std::ifstream input(path);
  if (!input) {
    throw Error("cannot open GAF file: " + path);
  }
  return read_gaf_lookup(input);
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

}  // namespace

int run_annotate(const AnnotateOptions& options) {
  const auto gaf_lookup = load_gaf_lookup(options.gaf_path);
  const std::unique_ptr<GraphIndex> graph = make_graph_index(GraphIndex::Config{
    options.gbz_path,
    options.gbwt_path,
    options.r_index_path,
  });

  const std::string fingerprint = sha256_file(options.use_gbz() ? options.gbz_path : options.gbwt_path);
  SidecarWriter sidecar(options.out_sets_path);
  sidecar.write_header(SidecarHeader{1, fingerprint, options.use_r_index()});

  samFile* input = sam_open(options.bam_path.c_str(), "rb");
  if (input == nullptr) {
    throw Error("cannot open BAM file: " + options.bam_path);
  }
  std::unique_ptr<samFile, decltype(&hts_close)> input_guard(input, hts_close);
  if (options.threads > 1) {
    hts_set_threads(input, static_cast<int>(options.threads));
  }

  sam_hdr_t* input_header = sam_hdr_read(input);
  if (input_header == nullptr) {
    throw Error("cannot read BAM header: " + options.bam_path);
  }
  std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)> input_header_guard(input_header, sam_hdr_destroy);

  sam_hdr_t* output_header = prepare_header(input_header, options);
  std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)> output_header_guard(output_header, sam_hdr_destroy);

  samFile* output = sam_open(options.out_bam_path.c_str(), "wb");
  if (output == nullptr) {
    throw Error("cannot open output BAM file: " + options.out_bam_path);
  }
  std::unique_ptr<samFile, decltype(&hts_close)> output_guard(output, hts_close);
  if (options.threads > 1) {
    hts_set_threads(output, static_cast<int>(options.threads));
  }

  if (sam_hdr_write(output, output_header) != 0) {
    throw Error("cannot write BAM header: " + options.out_bam_path);
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
  const std::size_t max_work_queue = std::max<std::size_t>(1, options.threads * 8);

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
        if (item.record != nullptr) {
          bam_destroy1(item.record);
          item.record = nullptr;
        }
      }
      work_queue.clear();
    }
    {
      std::lock_guard<std::mutex> lock(result_mutex);
      for (ResultItem& item : result_queue) {
        if (item.record != nullptr) {
          bam_destroy1(item.record);
          item.record = nullptr;
        }
      }
      result_queue.clear();
    }
  };

  std::vector<std::thread> workers;
  workers.reserve(options.threads);
  for (std::size_t thread_id = 0; thread_id < options.threads; ++thread_id) {
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
          result.record = item.record;
          if (item.walk != nullptr) {
            result.subpaths = graph->find_subpaths(*item.walk);
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
      std::unordered_map<std::vector<std::uint64_t>, std::uint32_t, VectorHash> set_to_id;
      std::map<std::size_t, ResultItem> pending;
      std::uint32_t next_set_id = 0;
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

            auto found = set_to_id.find(subpath.thread_ids);
            if (found == set_to_id.end()) {
              const std::uint32_t assigned = next_set_id++;
              set_to_id.emplace(subpath.thread_ids, assigned);
              sidecar.write_set(SidecarSetRecord{assigned, subpath.thread_ids});
              found = set_to_id.find(subpath.thread_ids);
            }

            annotation.hs.push_back(found->second);
            annotation.hb.push_back(subpath.begin_offset);
            annotation.he.push_back(subpath.end_offset);
          }

          set_or_remove_array_tag(item.record, "hs", annotation.hs);
          set_or_remove_array_tag(item.record, "hb", annotation.hb);
          set_or_remove_array_tag(item.record, "he", annotation.he);

          if (sam_write1(output, output_header, item.record) < 0) {
            bam_destroy1(item.record);
            throw Error("failed to write BAM record");
          }

          bam_destroy1(item.record);
          ++next_ordinal;
        }

        if (has_failure.load()) {
          for (auto& [_, pending_item] : pending) {
            if (pending_item.record != nullptr) {
              bam_destroy1(pending_item.record);
            }
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
    while (true) {
      const int read_status = sam_read1(input, input_header, record.get());
      if (read_status < -1) {
        throw Error("failed to read BAM record from " + options.bam_path);
      }
      if (read_status == -1) {
        break;
      }

      const std::string qname = bam_get_qname(record.get());
      const auto found = gaf_lookup.find(qname);

      bam1_t* duplicate = bam_dup1(record.get());
      if (duplicate == nullptr) {
        throw Error("unable to duplicate BAM record for " + qname);
      }

      WorkItem item;
      item.ordinal = ordinal++;
      item.record = duplicate;
      item.walk = (found == gaf_lookup.end() ? nullptr : &found->second);

      {
        std::unique_lock<std::mutex> lock(work_mutex);
        work_space.wait(lock, [&]() { return has_failure.load() || work_queue.size() < max_work_queue; });
        if (has_failure.load()) {
          bam_destroy1(duplicate);
          duplicate = nullptr;
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

  return 0;
}

}  // namespace pgbam
