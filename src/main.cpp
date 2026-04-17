#include <exception>
#include <iostream>

#include "pgbam/bam_io.hpp"
#include "pgbam/cli.hpp"

int main(int argc, char** argv) {
  try {
    const pgbam::CommandLine command_line = pgbam::parse_command_line(argc, argv);
    switch (command_line.command) {
      case pgbam::CommandLine::Command::Annotate:
        return pgbam::run_annotate(command_line.annotate);
      case pgbam::CommandLine::Command::Decode:
        return pgbam::run_decode(command_line.decode);
    }
  } catch (const std::exception& error) {
    std::cerr << "pgbam: " << error.what() << '\n';
  }

  return 1;
}
