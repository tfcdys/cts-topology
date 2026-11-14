// SLLT: bounded-skew tree fused with SALT shallow-light refinement.

#include "salt_common.hpp"

int main(int argc, char** argv) {
    return cts::run_topology_cli(
        argc, argv, "sllt",
        [](const std::vector<cts::Sink>& sinks,
           const cts::BuildOptions&) {
            return cts::build_sllt_tree(sinks);
        },
        false, false, false);
}
