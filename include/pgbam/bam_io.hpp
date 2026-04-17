#pragma once

#include "pgbam/options.hpp"

namespace pgbam {

int run_annotate(const AnnotateOptions& options);
int run_decode(const DecodeOptions& options);

}  // namespace pgbam
