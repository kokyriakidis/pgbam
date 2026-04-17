#pragma once

#include <stdexcept>
#include <string>

namespace pgbam {

class Error : public std::runtime_error {
public:
  explicit Error(const std::string& message) : std::runtime_error(message) {}
};

}  // namespace pgbam
