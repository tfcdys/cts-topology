// BST-DME: bounded-skew deferred-merge embedding.
//
// Algorithm concept reference:
//   https://github.com/abk-openroad/BST-DME
//
// This is an independent, simplified implementation using a linear delay model
// and a nearest-neighbor topology.  The global skew bound is converted to a
// per-merge budget; a deterministic right-child bias is used.

#include "bst_dme.hpp"
#include "cts_common.hpp"
#include "dme_common.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace cts {
namespace {

Topology build_bst_dme_impl(const std::vector<Sink>& sinks, double skew_bound_um) {
    DmeTree tree = build_nearest_neighbor_tree(sinks);
    const std::size_t n = sinks.size();
    const int depth = dme_tree_depth(tree);
    const double local_skew =
        depth > 0 ? skew_bound_um / static_cast<double>(depth) : 0.0;

    // Bottom-up: zero-skew radii, then reallocate a bounded amount of delay
    // from the left child to the right child while preserving eL + eR = L.
    for (std::size_t i = n; i < tree.nodes.size(); ++i) {
        DmeNode& node = tree.nodes[i];
        DmeNode& left = tree.nodes[static_cast<std::size_t>(node.left)];
        DmeNode& right = tree.nodes[static_cast<std::size_t>(node.right)];

        const double distance = dme_segment_distance(left.core, right.core);
        double e_left = (right.delay - left.delay + distance) / 2.0;
        double e_right = (left.delay - right.delay + distance) / 2.0;
        if (e_left < 0.0) {
            e_left = 0.0;
            e_right = std::max(0.0, left.delay - right.delay);
        } else if (e_right < 0.0) {
            e_right = 0.0;
            e_left = std::max(0.0, right.delay - left.delay);
        }

        if (local_skew > 0.0) {
            double delta = -0.5 * local_skew;
            if (e_left + delta < 0.0) {
                delta = -e_left;
            }
            if (e_right - delta < 0.0) {
                delta = e_right;
            }
            e_left += delta;
            e_right -= delta;
        }

        left.radius_to_parent = e_left;
        right.radius_to_parent = e_right;
        node.core = dme_merge_segment(left.core, e_left, right.core, e_right);
        node.delay = std::max(left.delay + e_left, right.delay + e_right);
        node.capacitance =
            left.capacitance + right.capacitance + e_left + e_right;
    }

    // Top-down: root is the midpoint of the root merging segment; all
    // non-sink children are embedded exactly on their merging segments.
    tree.nodes[static_cast<std::size_t>(tree.root)].pos =
        dme_midpoint(tree.nodes[static_cast<std::size_t>(tree.root)].core.a,
                     tree.nodes[static_cast<std::size_t>(tree.root)].core.b);
    std::vector<int> stack{tree.root};
    while (!stack.empty()) {
        const int current = stack.back();
        stack.pop_back();
        const DmeNode& parent = tree.nodes[static_cast<std::size_t>(current)];
        for (const int child_id : {parent.left, parent.right}) {
            if (child_id < 0) {
                continue;
            }
            DmeNode& child = tree.nodes[static_cast<std::size_t>(child_id)];
            if (child.sink_index >= 0) {
                const Sink& sink =
                    sinks[static_cast<std::size_t>(child.sink_index)];
                child.pos = dme_point(sink.x, sink.y);
            } else {
                child.pos = dme_choose_exact(parent.pos, child.core,
                                             child.radius_to_parent);
            }
            stack.push_back(child_id);
        }
    }
    return emit_dme_topology(sinks, tree);
}

}  // namespace

Topology build_bst_dme(const std::vector<Sink>& sinks, double skew_bound_um) {
    return build_bst_dme_impl(sinks, skew_bound_um);
}

}  // namespace cts

#ifndef BST_DME_NO_MAIN
int main(int argc, char** argv) {
    return cts::run_topology_cli(
        argc, argv, "bst_dme",
        [](const std::vector<cts::Sink>& sinks,
           const cts::BuildOptions& options) {
            return cts::build_bst_dme(sinks, options.skew_bound_um);
        },
        false, true, true);
}
#endif
