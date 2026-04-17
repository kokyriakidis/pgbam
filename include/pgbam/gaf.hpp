#pragma once

#include <istream>
#include <optional>
#include <string>
#include <vector>

#include "pgbam/types.hpp"

namespace pgbam {

std::vector<OrientedNode> parse_target_walk(const std::string& walk);
std::optional<GafRecord> parse_gaf_line(const std::string& line);
std::vector<GafRecord> read_gaf_records(std::istream& in);

}  // namespace pgbam
