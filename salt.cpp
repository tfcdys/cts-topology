// Standalone SALT topology executable.

#include "salt_common.hpp"

int main(int argc, char** argv) {
    return cts::run_topology_cli(
        argc, argv, "salt",
        [](const std::vector<cts::Sink>& sinks,
           const cts::BuildOptions&) {
            return cts::build_salt_tree(sinks);
        },
        false, false, false);
}
