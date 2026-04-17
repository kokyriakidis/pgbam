#include "pgbam/cli.hpp"

#include <charconv>
#include <limits>
#include <string>
#include <string_view>

#include "pgbam/error.hpp"

namespace pgbam {
namespace {

std::string require_value(int& index, int argc, char** argv, std::string_view option) {
  if (index + 1 >= argc) {
    throw Error("missing value for option " + std::string(option));
  }
  ++index;
  return argv[index];
}

std::size_t parse_threads(const std::string& value) {
  std::size_t parsed = 0;
  const char* begin = value.data();
  const char* end = value.data() + value.size();
  const auto result = std::from_chars(begin, end, parsed);
  if (result.ec != std::errc{} || result.ptr != end || parsed == 0) {
    throw Error("--threads must be at least 1");
  }
  return parsed;
}

void validate_annotate(const AnnotateOptions& options) {
  if (options.bam_path.empty()) {
    throw Error("--bam is required");
  }
  if (options.gaf_path.empty()) {
    throw Error("--gaf is required");
  }
  if (options.out_bam_path.empty()) {
    throw Error("--out-bam is required");
  }
  if (options.out_sets_path.empty()) {
    throw Error("--out-sets is required");
  }
  if (options.gbz_path.empty() == options.gbwt_path.empty()) {
    throw Error("exactly one of --gbz or --gbwt is required");
  }
}

void validate_decode(const DecodeOptions& options) {
  if (options.sets_path.empty()) {
    throw Error("--sets is required");
  }
  if (options.out_path.empty()) {
    throw Error("--out is required");
  }
  if (options.gbz_path.empty() == options.gbwt_path.empty()) {
    throw Error("exactly one of --gbz or --gbwt is required");
  }
}

}  // namespace

CommandLine parse_command_line(int argc, char** argv) {
  if (argc < 2) {
    throw Error("expected a subcommand: annotate or decode");
  }

  CommandLine result;
  const std::string command = argv[1];
  if (command == "annotate") {
    result.command = CommandLine::Command::Annotate;
    for (int index = 2; index < argc; ++index) {
      const std::string_view option = argv[index];
      if (option == "--bam") {
        result.annotate.bam_path = require_value(index, argc, argv, option);
      } else if (option == "--gaf") {
        result.annotate.gaf_path = require_value(index, argc, argv, option);
      } else if (option == "--gbz") {
        result.annotate.gbz_path = require_value(index, argc, argv, option);
      } else if (option == "--gbwt") {
        result.annotate.gbwt_path = require_value(index, argc, argv, option);
      } else if (option == "--r-index") {
        result.annotate.r_index_path = require_value(index, argc, argv, option);
      } else if (option == "--out-bam") {
        result.annotate.out_bam_path = require_value(index, argc, argv, option);
      } else if (option == "--out-sets") {
        result.annotate.out_sets_path = require_value(index, argc, argv, option);
      } else if (option == "--threads") {
        result.annotate.threads = parse_threads(require_value(index, argc, argv, option));
      } else {
        throw Error("unknown annotate option: " + std::string(option));
      }
    }
    validate_annotate(result.annotate);
    return result;
  }

  if (command == "decode") {
    result.command = CommandLine::Command::Decode;
    for (int index = 2; index < argc; ++index) {
      const std::string_view option = argv[index];
      if (option == "--sets") {
        result.decode.sets_path = require_value(index, argc, argv, option);
      } else if (option == "--gbz") {
        result.decode.gbz_path = require_value(index, argc, argv, option);
      } else if (option == "--gbwt") {
        result.decode.gbwt_path = require_value(index, argc, argv, option);
      } else if (option == "--out") {
        result.decode.out_path = require_value(index, argc, argv, option);
      } else {
        throw Error("unknown decode option: " + std::string(option));
      }
    }
    validate_decode(result.decode);
    return result;
  }

  throw Error("unknown subcommand: " + command);
}

}  // namespace pgbam
