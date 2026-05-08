#pragma once

#include <stdexcept>
#include <string>

namespace pgbam {

class Error : public std::runtime_error {
public:
  explicit Error(const std::string& message) : std::runtime_error(message) {}
};

class MemoryBudgetExceeded : public Error {
public:
  MemoryBudgetExceeded(const std::string& message, std::size_t estimated_bytes, std::size_t budget_bytes)
    : Error(message), estimated_bytes_(estimated_bytes), budget_bytes_(budget_bytes) {}

  std::size_t estimated_bytes() const { return estimated_bytes_; }
  std::size_t budget_bytes() const { return budget_bytes_; }

private:
  std::size_t estimated_bytes_ = 0;
  std::size_t budget_bytes_ = 0;
};

}  // namespace pgbam
