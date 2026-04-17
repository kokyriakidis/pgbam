#include "pgbam/gbwt_backend.hpp"

#include <algorithm>
#include <fstream>
#include <limits>
#include <memory>
#include <vector>

#include <gbwt/fast_locate.h>
#include <gbwt/gbwt.h>
#include <gbwtgraph/gbz.h>
#include <gbwtgraph/utils.h>
#include <handlegraph/path_metadata.hpp>

#include "pgbam/error.hpp"

namespace pgbam {
namespace {

class GbwtGraphIndex final : public GraphIndex {
public:
  explicit GbwtGraphIndex(const Config& config) {
    if (!config.gbz_path.empty()) {
      auto gbz = std::make_unique<gbwtgraph::GBZ>();
      std::ifstream input(config.gbz_path, std::ios::binary);
      if (!input) {
        throw Error("cannot open GBZ file: " + config.gbz_path);
      }
      gbz->simple_sds_load(input);
      index_ = &gbz->index;
      gbz_ = std::move(gbz);
    } else {
      auto index = std::make_unique<gbwt::GBWT>();
      std::ifstream input(config.gbwt_path, std::ios::binary);
      if (!input) {
        throw Error("cannot open GBWT file: " + config.gbwt_path);
      }
      index->simple_sds_load(input);
      index_ = index.get();
      owned_index_ = std::move(index);
    }

    reference_samples_ = gbwtgraph::parse_reference_samples_tag(*index_);

    if (!config.r_index_path.empty()) {
      auto fast_locate = std::make_unique<gbwt::FastLocate>();
      std::ifstream input(config.r_index_path, std::ios::binary);
      if (!input) {
        throw Error("cannot open r-index file: " + config.r_index_path);
      }
      fast_locate->load(input);
      fast_locate->setGBWT(*index_);
      fast_locate_ = std::move(fast_locate);
    }
  }

  std::vector<SubpathResult> find_subpaths(const std::vector<OrientedNode>& walk) const override {
    std::vector<gbwt::node_type> encoded;
    encoded.reserve(walk.size());
    for (const OrientedNode& node : walk) {
      encoded.push_back(gbwt::Node::encode(node.id, node.reverse));
    }

    std::vector<SubpathResult> result;
    std::size_t start = 0;
    while (start < encoded.size()) {
      gbwt::SearchState state;
      gbwt::FastLocate::size_type first_occurrence = gbwt::FastLocate::NO_POSITION;
      if (fast_locate_) {
        state = fast_locate_->find(encoded[start], first_occurrence);
      } else {
        state = index_->find(encoded[start]);
      }

      if (state.empty()) {
        ++start;
        continue;
      }

      std::size_t end = start + 1;
      while (end < encoded.size()) {
        gbwt::SearchState next_state;
        gbwt::FastLocate::size_type next_first = first_occurrence;
        if (fast_locate_) {
          next_state = fast_locate_->extend(state, encoded[end], next_first);
        } else {
          next_state = index_->extend(state, encoded[end]);
        }
        if (next_state.empty()) {
          break;
        }
        state = next_state;
        first_occurrence = next_first;
        ++end;
      }

      std::vector<gbwt::size_type> located;
      if (fast_locate_) {
        located = fast_locate_->locate(state, first_occurrence);
      } else {
        located = index_->locate(state);
      }

      std::vector<std::uint64_t> path_ids;
      path_ids.reserve(located.size());
      for (gbwt::size_type sequence_id : located) {
        path_ids.push_back(static_cast<std::uint64_t>(gbwt::Path::id(sequence_id)));
      }
      std::sort(path_ids.begin(), path_ids.end());
      path_ids.erase(std::unique(path_ids.begin(), path_ids.end()), path_ids.end());

      result.push_back(SubpathResult{
        static_cast<std::uint32_t>(start),
        static_cast<std::uint32_t>(end),
        std::move(path_ids),
      });
      start = end;
    }

    return result;
  }

  ThreadMetadata decode_thread(std::uint64_t thread_id) const override {
    const gbwt::size_type path_id = static_cast<gbwt::size_type>(thread_id);
    const gbwtgraph::PathSense sense = gbwtgraph::get_path_sense(*index_, path_id, reference_samples_);

    ThreadMetadata metadata;
    metadata.thread_id = thread_id;
    metadata.path_id = thread_id;
    metadata.sample = gbwtgraph::get_path_sample_name(*index_, path_id, sense);
    const std::size_t haplotype = gbwtgraph::get_path_haplotype(*index_, path_id, sense);
    metadata.haplotype_known = (haplotype != handlegraph::PathMetadata::NO_HAPLOTYPE);
    metadata.haplotype = static_cast<std::uint64_t>(haplotype);
    metadata.locus = gbwtgraph::get_path_locus_name(*index_, path_id, sense);
    metadata.path_name = gbwtgraph::compose_path_name(*index_, path_id, sense);
    return metadata;
  }

private:
  std::unique_ptr<gbwtgraph::GBZ> gbz_;
  std::unique_ptr<gbwt::GBWT> owned_index_;
  const gbwt::GBWT* index_ = nullptr;
  std::unique_ptr<gbwt::FastLocate> fast_locate_;
  gbwtgraph::sample_name_set reference_samples_;
};

}  // namespace

std::unique_ptr<GraphIndex> make_graph_index(const GraphIndex::Config& config) {
  return std::make_unique<GbwtGraphIndex>(config);
}

}  // namespace pgbam
