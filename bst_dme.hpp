#pragma once

#include "cts_common.hpp"

#include <vector>

namespace cts {

Topology build_bst_dme(const std::vector<Sink>& sinks, double skew_bound_um);

}  // namespace cts
