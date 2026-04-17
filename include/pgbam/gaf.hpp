#pragma once

#include <istream>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "pgbam/types.hpp"

namespace pgbam {

std::vector<OrientedNode> parse_target_walk(std::string_view walk);
std::optional<GafRecord> parse_gaf_line(const std::string& line);
std::vector<GafRecord> read_gaf_records(std::istream& in);
std::unordered_map<std::string, std::vector<OrientedNode>> read_gaf_lookup(std::istream& in);

}  // namespace pgbam
