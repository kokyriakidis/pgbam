#include <cassert>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "pgbam/cli.hpp"
#include "pgbam/error.hpp"
#include "pgbam/gaf.hpp"
#include "pgbam/sidecar.hpp"

int main() {
  static_assert(sizeof(pgbam::OrientedNode) == sizeof(std::uint64_t));

  {
    const auto nodes = pgbam::parse_target_walk(">2>3<4");
    const pgbam::OrientedNode first{2, false};
    const pgbam::OrientedNode second{3, false};
    const pgbam::OrientedNode third{4, true};
    assert(nodes.size() == 3);
    assert(nodes[0] == first);
    assert(nodes[1] == second);
    assert(nodes[2] == third);
  }

  {
    const auto parsed = pgbam::parse_gaf_line("read1\t6\t0\t6\t+\t>2>3>4\t12\t2\t8\t6\t6\t60\tcs:Z::6");
    const pgbam::OrientedNode first{2, false};
    const pgbam::OrientedNode last{4, false};
    assert(parsed.has_value());
    assert(parsed->qname == "read1");
    assert(parsed->query_length == 6);
    assert(parsed->nodes.size() == 3);
    assert(parsed->nodes[0] == first);
    assert(parsed->nodes[2] == last);
  }

  {
    std::istringstream gaf(
      "read1\t6\t0\t6\t+\t>2>3>4\t12\t2\t8\t6\t6\t60\tcs:Z::6\n"
      "read2\t6\t0\t6\t+\t>5<6\t12\t2\t8\t6\t6\t60\tcs:Z::6\n");
    const auto lookup = pgbam::read_gaf_lookup(gaf);
    assert(lookup.size() == 2);
    assert(lookup.at("read1").size() == 3);
    assert(lookup.at("read2").size() == 2);
  }

  {
    std::istringstream gaf(
      "read1\t6\t0\t6\t+\t>2>3>4\t12\t2\t8\t6\t6\t60\tcs:Z::6\n"
      "read1\t6\t0\t6\t+\t>5<6\t12\t2\t8\t6\t6\t60\tcs:Z::6\n");
    bool threw = false;
    try {
      (void) pgbam::read_gaf_lookup(gaf);
    } catch (const pgbam::Error&) {
      threw = true;
    }
    assert(threw);
  }

  {
    bool threw = false;
    try {
      (void) pgbam::parse_target_walk(">2>");
    } catch (const pgbam::Error&) {
      threw = true;
    }
    assert(threw);
  }

  {
    const std::filesystem::path path = std::filesystem::temp_directory_path() / "pgbam_test_sidecar.pgs";
    {
      pgbam::SidecarWriter writer(path.string());
      writer.write_header(pgbam::SidecarHeader{1, "fingerprint", true});
      writer.write_set(pgbam::SidecarSetRecord{7, {1, 2, 3}});
    }

    const pgbam::LoadedSidecar loaded = pgbam::read_sidecar(path.string());
    assert(loaded.header.version == 1);
    assert(loaded.header.fingerprint == "fingerprint");
    assert(loaded.header.used_r_index);
    assert(loaded.sets.at(7) == std::vector<std::uint64_t>({1, 2, 3}));
    std::filesystem::remove(path);
  }

  {
    const std::filesystem::path path = std::filesystem::temp_directory_path() / "pgbam_test_sidecar_duplicate.pgs";
    {
      std::ofstream out(path, std::ios::binary);
      const char magic[] = {'P', 'G', 'S', '1'};
      const std::uint32_t version = 1;
      const std::uint8_t used_r_index = 0;
      const std::uint32_t fingerprint_size = 2;
      const char fingerprint[] = {'o', 'k'};
      const std::uint32_t set_id = 4;
      const std::uint32_t count = 1;
      const std::uint64_t thread = 9;
      out.write(magic, sizeof(magic));
      out.write(reinterpret_cast<const char*>(&version), sizeof(version));
      out.write(reinterpret_cast<const char*>(&used_r_index), sizeof(used_r_index));
      out.write(reinterpret_cast<const char*>(&fingerprint_size), sizeof(fingerprint_size));
      out.write(fingerprint, sizeof(fingerprint));
      out.write(reinterpret_cast<const char*>(&set_id), sizeof(set_id));
      out.write(reinterpret_cast<const char*>(&count), sizeof(count));
      out.write(reinterpret_cast<const char*>(&thread), sizeof(thread));
      out.write(reinterpret_cast<const char*>(&set_id), sizeof(set_id));
      out.write(reinterpret_cast<const char*>(&count), sizeof(count));
      out.write(reinterpret_cast<const char*>(&thread), sizeof(thread));
    }

    bool threw = false;
    try {
      (void) pgbam::read_sidecar(path.string());
    } catch (const pgbam::Error&) {
      threw = true;
    }
    assert(threw);
    std::filesystem::remove(path);
  }

  {
    std::vector<std::string> args = {
      "pgbam",
      "annotate",
      "--bam", "in.bam",
      "--gaf", "in.gaf",
      "--gbwt", "graph.gbwt",
      "--out-bam", "out.bam",
      "--out-sets", "out.pgs",
      "--threads", "4",
    };
    std::vector<char*> argv;
    argv.reserve(args.size());
    for (std::string& arg : args) {
      argv.push_back(arg.data());
    }

    const pgbam::CommandLine command = pgbam::parse_command_line(static_cast<int>(argv.size()), argv.data());
    assert(command.command == pgbam::CommandLine::Command::Annotate);
    assert(command.annotate.gbwt_path == "graph.gbwt");
    assert(command.annotate.threads == 4);
    assert(!command.annotate.primary_only);
  }

  {
    std::vector<std::string> args = {
      "pgbam",
      "annotate",
      "--bam", "in.bam",
      "--gaf", "in.gaf",
      "--gbwt", "graph.gbwt",
      "--out-bam", "out.bam",
      "--out-sets", "out.pgs",
      "--primary-only",
    };
    std::vector<char*> argv;
    argv.reserve(args.size());
    for (std::string& arg : args) {
      argv.push_back(arg.data());
    }

    const pgbam::CommandLine command = pgbam::parse_command_line(static_cast<int>(argv.size()), argv.data());
    assert(command.annotate.primary_only);
  }

  {
    std::vector<std::string> args = {
      "pgbam",
      "annotate",
      "--bam", "in.bam",
      "--gaf", "in.gaf",
      "--gbwt", "graph.gbwt",
      "--out-bam", "out.bam",
      "--out-sets", "out.pgs",
      "--threads", "4x",
    };
    std::vector<char*> argv;
    argv.reserve(args.size());
    for (std::string& arg : args) {
      argv.push_back(arg.data());
    }

    bool threw = false;
    try {
      (void) pgbam::parse_command_line(static_cast<int>(argv.size()), argv.data());
    } catch (const pgbam::Error&) {
      threw = true;
    }
    assert(threw);
  }

  {
    std::vector<std::string> args = {
      "pgbam",
      "annotate",
      "--bam", "in.bam",
      "--gaf", "in.gaf",
      "--gbwt", "graph.gbwt",
      "--out-bam", "out.bam",
      "--out-sets", "out.pgs",
      "--invalid-option",
    };
    std::vector<char*> argv;
    argv.reserve(args.size());
    for (std::string& arg : args) {
      argv.push_back(arg.data());
    }

    bool threw = false;
    try {
      (void) pgbam::parse_command_line(static_cast<int>(argv.size()), argv.data());
    } catch (const pgbam::Error&) {
      threw = true;
    }
    assert(threw);
  }

  return 0;
}
