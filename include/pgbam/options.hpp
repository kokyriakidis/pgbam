#pragma once

#include <cstddef>
#include <optional>
#include <string>

namespace pgbam {

struct AnnotateOptions {
  std::string bam_path;
  std::string gaf_path;
  std::string gbz_path;
  std::string gbwt_path;
  std::string r_index_path;
  std::string out_bam_path;
  std::string out_sets_path;
  std::size_t threads = 1;

  bool use_gbz() const { return !gbz_path.empty(); }
  bool use_r_index() const { return !r_index_path.empty(); }
};

struct DecodeOptions {
  std::string sets_path;
  std::string gbz_path;
  std::string gbwt_path;
  std::string out_path;

  bool use_gbz() const { return !gbz_path.empty(); }
};

struct CommandLine {
  enum class Command { Annotate, Decode };

  Command command = Command::Annotate;
  AnnotateOptions annotate;
  DecodeOptions decode;
};

}  // namespace pgbam
