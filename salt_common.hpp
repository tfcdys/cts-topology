#pragma once

#include "cts_common.hpp"

#include <vector>

namespace cts {

struct SaltOptions {
    double eps = 0.0;
    int refine_level = 3;
};

struct SlltOptions {
    double skew_bound_um = 100.0;
    double eps = 0.0;
    int refine_level = 3;
};

Topology build_salt_tree(const std::vector<Sink>& sinks,
                         const SaltOptions& options = {});

Topology build_sllt_tree(const std::vector<Sink>& sinks,
                         const SlltOptions& options = {});

}  // namespace cts
