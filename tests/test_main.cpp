#include <cassert>
#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "pgbam/bam_io.hpp"
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
      "annotate-bam",
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
    assert(command.command == pgbam::CommandLine::Command::AnnotateBam);
    assert(command.annotate_bam.gbwt_path == "graph.gbwt");
    assert(command.annotate_bam.threads == 4);
    assert(!command.annotate_bam.primary_only);
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
    };
    std::vector<char*> argv;
    argv.reserve(args.size());
    for (std::string& arg : args) {
      argv.push_back(arg.data());
    }

    const pgbam::CommandLine command = pgbam::parse_command_line(static_cast<int>(argv.size()), argv.data());
    assert(command.command == pgbam::CommandLine::Command::AnnotateBam);
  }

  {
    std::vector<std::string> args = {
      "pgbam",
      "annotate-bam",
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
    assert(command.annotate_bam.primary_only);
  }

  {
    std::vector<std::string> args = {
      "pgbam",
      "annotate-bam",
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
      "annotate-bam",
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

  // qname_less: basic ordering
  {
    assert(pgbam::qname_less("read/1/ccs", "read/2/ccs"));
    assert(pgbam::qname_less("read/9/ccs", "read/10/ccs"));
    assert(!pgbam::qname_less("read/10/ccs", "read/9/ccs"));
    assert(!pgbam::qname_less("read/1/ccs", "read/1/ccs"));
    assert(pgbam::qname_less("a", "b"));
    assert(!pgbam::qname_less("b", "a"));
  }

  // qname_less: the exact digit-length collision from HG002 HiFi data.
  // Without LC_ALL=C, locale-aware sort places 100270982 before 10027098
  // because it ignores '/' separators, making "10027098ccs" > "100270982ccs"
  // at position 8 ('c' > '2').  qname_less must put 10027098 first.
  {
    const std::string short_zmw = "m84031_231217_034919_s2/10027098/ccs";
    const std::string long_zmw_a = "m84031_231217_034919_s2/100270982/ccs";
    const std::string long_zmw_b = "m84031_231217_034919_s2/100270986/ccs";
    const std::string long_zmw_c = "m84031_231217_034919_s2/100270992/ccs";

    assert(pgbam::qname_less(short_zmw, long_zmw_a));
    assert(pgbam::qname_less(short_zmw, long_zmw_b));
    assert(pgbam::qname_less(short_zmw, long_zmw_c));
    assert(!pgbam::qname_less(long_zmw_a, short_zmw));
    assert(!pgbam::qname_less(long_zmw_b, short_zmw));
    assert(!pgbam::qname_less(long_zmw_c, short_zmw));
    assert(pgbam::qname_less(long_zmw_a, long_zmw_b));
    assert(pgbam::qname_less(long_zmw_b, long_zmw_c));
  }

  // sort-gaf CLI parsing
  {
    std::vector<std::string> args = {
      "pgbam", "sort-gaf",
      "--in", "input.gaf",
      "--out", "sorted.gaf",
    };
    std::vector<char*> argv;
    argv.reserve(args.size());
    for (std::string& arg : args) {
      argv.push_back(arg.data());
    }
    const pgbam::CommandLine command = pgbam::parse_command_line(static_cast<int>(argv.size()), argv.data());
    assert(command.command == pgbam::CommandLine::Command::SortGaf);
    assert(command.sort_gaf.in_path == "input.gaf");
    assert(command.sort_gaf.out_path == "sorted.gaf");
  }

  // sort-gaf: missing --in or --out must throw
  {
    std::vector<std::string> args = {"pgbam", "sort-gaf", "--out", "sorted.gaf"};
    std::vector<char*> argv;
    for (std::string& arg : args) argv.push_back(arg.data());
    bool threw = false;
    try { (void) pgbam::parse_command_line(static_cast<int>(argv.size()), argv.data()); }
    catch (const pgbam::Error&) { threw = true; }
    assert(threw);
  }

  // sort-gaf integration: locale-naive order in → correct order out
  {
    const std::filesystem::path in_path  = std::filesystem::temp_directory_path() / "pgbam_test_sort_gaf_in.gaf";
    const std::filesystem::path out_path = std::filesystem::temp_directory_path() / "pgbam_test_sort_gaf_out.gaf";

    // Write qnames in locale-naive (wrong) order: long ZMW before short ZMW.
    // "LC_ALL=C sort -k1,1" must reorder them so 10027098 comes first.
    {
      std::ofstream f(in_path);
      f << "m84031_231217_034919_s2/100270982/ccs\t10\t0\t10\t+\t>1\t10\t0\t10\t10\t10\t60\n";
      f << "m84031_231217_034919_s2/100270986/ccs\t10\t0\t10\t+\t>1\t10\t0\t10\t10\t10\t55\n";
      f << "m84031_231217_034919_s2/10027098/ccs\t10\t0\t10\t+\t>1\t10\t0\t10\t10\t10\t50\n";
      f << "m84031_231217_034919_s2/100270992/ccs\t10\t0\t10\t+\t>1\t10\t0\t10\t10\t10\t45\n";
    }

    pgbam::SortGafOptions opts;
    opts.in_path  = in_path.string();
    opts.out_path = out_path.string();
    assert(pgbam::run_sort_gaf(opts) == 0);

    // First qname in output must be the 8-digit ZMW (10027098).
    std::ifstream result(out_path);
    std::string first_line;
    std::getline(result, first_line);
    assert(first_line.find("m84031_231217_034919_s2/10027098/ccs") == 0);

    std::filesystem::remove(in_path);
    std::filesystem::remove(out_path);
  }

  return 0;
}
