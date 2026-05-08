#include "pgbam/bam_io.hpp"

#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

#include "pgbam/error.hpp"

namespace pgbam {

int run_sort_gaf(const SortGafOptions& options) {
  const int out_fd = open(options.out_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
  if (out_fd < 0) {
    throw Error("cannot open output file for sort-gaf: " + options.out_path);
  }

  const pid_t pid = fork();
  if (pid < 0) {
    close(out_fd);
    throw Error("fork failed for sort-gaf");
  }

  if (pid == 0) {
    dup2(out_fd, STDOUT_FILENO);
    close(out_fd);
    // LC_ALL=C: prevents locale from treating '/' and '_' as non-significant.
    // -k1,1V:   version sort — digit spans compared numerically, matching qname_less.
    // -k12,12nr: within each qname group, highest MAPQ first.
    setenv("LC_ALL", "C", 1);
    execlp("sort", "sort", "-k1,1V", "-k12,12nr",
           options.in_path.c_str(), static_cast<char*>(nullptr));
    _exit(127);
  }

  close(out_fd);

  int status = 0;
  waitpid(pid, &status, 0);

  if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
    const int code = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
    throw Error("sort-gaf: sort exited with code " + std::to_string(code) +
                " (GNU sort is required; on macOS: brew install coreutils)");
  }

  return 0;
}

}  // namespace pgbam
