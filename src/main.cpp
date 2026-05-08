#include <exception>
#include <iostream>

#include "pgbam/bam_io.hpp"
#include "pgbam/cli.hpp"

int main(int argc, char** argv) {
  try {
    const pgbam::CommandLine command_line = pgbam::parse_command_line(argc, argv);
    switch (command_line.command) {
      case pgbam::CommandLine::Command::AnnotateBam:
        return pgbam::run_annotate_bam(command_line.annotate_bam);
      case pgbam::CommandLine::Command::Decode:
        return pgbam::run_decode(command_line.decode);
      case pgbam::CommandLine::Command::SortGaf:
        return pgbam::run_sort_gaf(command_line.sort_gaf);
    }
  } catch (const std::exception& error) {
    std::cerr << "pgbam: " << error.what() << '\n';
  }

  return 1;
}
