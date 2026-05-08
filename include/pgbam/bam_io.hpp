#pragma once

#include "pgbam/options.hpp"

namespace pgbam {

int run_annotate_bam(const AnnotateOptions& options);
int run_decode(const DecodeOptions& options);
int run_sort_gaf(const SortGafOptions& options);

}  // namespace pgbam
