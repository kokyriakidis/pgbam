#pragma once

#include <istream>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "pgbam/types.hpp"

namespace pgbam {

// Natural-sort comparison: digit runs are compared numerically so that
// e.g. "read/9/ccs" < "read/10/ccs". Matches the ordering used by
// `samtools sort -n` and expected by pgbam's GAF sort validation.
bool qname_less(std::string_view a, std::string_view b);

std::vector<OrientedNode> parse_target_walk(std::string_view walk);
std::optional<GafRecord> parse_gaf_line(const std::string& line);
std::vector<GafRecord> read_gaf_records(std::istream& in);
std::unordered_map<std::string, std::vector<OrientedNode>> read_gaf_lookup(std::istream& in,
                                                                           std::size_t memory_budget_bytes = 0);

}  // namespace pgbam
