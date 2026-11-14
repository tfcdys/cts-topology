// ZST-DME: zero-skew deferred-merge embedding.
//
// Algorithm concept reference:
//   https://github.com/geziangfinn/ZSTDME
//
// This is an independent, simplified implementation using a linear delay model
// and a nearest-neighbor topology.  See docs/dme_algorithms.md.

#include "cts_common.hpp"
#include "dme_common.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace cts {
namespace {

// Snaking post-pass: when a leaf path is shorter than the zero-skew target,
// insert a small U-shaped detour on its parent edge.  This is the missing
// top-down snaking step for cases where a clamped merging radius cannot be met
// exactly on the child's merging segment.
Topology apply_zero_skew_snaking(Topology topology) {
    const int n = topology.sink_nodes;
    if (n <= 1 || topology.root < 0) {
        return topology;
    }
    const std::size_t node_count = topology.nodes.size();
    std::vector<std::vector<std::pair<int, int>>> children(node_count);
    std::vector<int> parent_edge(node_count, -1);
    std::vector<int> parent(node_count, -1);
    for (std::size_t i = 0; i < topology.edges.size(); ++i) {
        const Edge& edge = topology.edges[i];
        children[static_cast<std::size_t>(edge.source)].push_back(
            {edge.target, static_cast<int>(i)});
        parent_edge[static_cast<std::size_t>(edge.target)] =
            static_cast<int>(i);
        parent[static_cast<std::size_t>(edge.target)] = edge.source;
    }

    std::vector<double> path(node_count, 0.0);
    std::vector<int> stack{topology.root};
    while (!stack.empty()) {
        const int current = stack.back();
        stack.pop_back();
        for (const auto& [child, edge_index] : children[static_cast<std::size_t>(current)]) {
            const Edge& edge = topology.edges[static_cast<std::size_t>(edge_index)];
            const Node& a = topology.nodes[static_cast<std::size_t>(edge.source)];
            const Node& b = topology.nodes[static_cast<std::size_t>(edge.target)];
            path[static_cast<std::size_t>(child)] =
                path[static_cast<std::size_t>(current)] +
                std::abs(a.x - b.x) + std::abs(a.y - b.y);
            stack.push_back(child);
        }
    }

    double target = 0.0;
    for (const Node& node : topology.nodes) {
        if (node.kind == "sink") {
            target = std::max(target, path[static_cast<std::size_t>(node.id)]);
        }
    }

    struct Adjustment {
        int edge_index;
        int sink_node;
        double delta;
    };
    std::vector<Adjustment> adjustments;
    double added = 0.0;
    for (const Node& node : topology.nodes) {
        if (node.kind != "sink") {
            continue;
        }
        const double delta =
            target - path[static_cast<std::size_t>(node.id)];
        if (delta <= 1e-6) {
            continue;
        }
        const int edge_index = parent_edge[static_cast<std::size_t>(node.id)];
        if (edge_index < 0) {
            continue;
        }
        adjustments.push_back(
            Adjustment{edge_index, node.id, delta});
        added += delta;
    }

    if (adjustments.empty()) {
        return topology;
    }

    for (const Adjustment& adjustment : adjustments) {
        const Edge& edge =
            topology.edges[static_cast<std::size_t>(adjustment.edge_index)];
        const Node& parent_node =
            topology.nodes[static_cast<std::size_t>(edge.source)];
        const Node& sink_node =
            topology.nodes[static_cast<std::size_t>(adjustment.sink_node)];
        const double dx = sink_node.x - parent_node.x;
        const double dy = sink_node.y - parent_node.y;
        const double half = 0.5 * adjustment.delta;
        double qx = sink_node.x;
        double qy = sink_node.y;
        if (std::abs(dx) >= std::abs(dy) && std::abs(dx) > 1e-12) {
            qx += (dx >= 0.0 ? 1.0 : -1.0) * half;
        } else {
            qy += (dy >= 0.0 ? 1.0 : -1.0) * half;
        }
        Node bend;
        bend.id = static_cast<int>(topology.nodes.size());
        bend.x = qx;
        bend.y = qy;
        bend.kind = "bend";
        topology.nodes.push_back(bend);
        ++topology.bend_nodes;
        const int bend_id = bend.id;
        topology.edges[static_cast<std::size_t>(adjustment.edge_index)].target =
            bend_id;
        topology.edges.push_back(
            Edge{bend_id, adjustment.sink_node, "HV"});
    }

    std::ostringstream metadata;
    metadata << "{\"snaking_adjustments\":" << adjustments.size()
             << ",\"added_wirelength_um\":" << format_number(added) << "}";
    topology.metadata_json = metadata.str();
    return topology;
}

Topology build_zst_dme(const std::vector<Sink>& sinks) {
    DmeTree tree = build_nearest_neighbor_tree(sinks);
    const std::size_t n = sinks.size();

    // Bottom-up: compute child radii, parent merging segment, and delay.
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

        left.radius_to_parent = e_left;
        right.radius_to_parent = e_right;
        node.core = dme_merge_segment(left.core, e_left, right.core, e_right);
        node.delay = std::max(left.delay + e_left, right.delay + e_right);
        node.capacitance =
            left.capacitance + right.capacitance + e_left + e_right;
    }

    // Top-down: choose the root and exact child locations on merging segments.
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
    return apply_zero_skew_snaking(emit_dme_topology(sinks, tree));
}

}  // namespace
}  // namespace cts

int main(int argc, char** argv) {
    return cts::run_topology_cli(
        argc, argv, "zst_dme",
        [](const std::vector<cts::Sink>& sinks,
           const cts::BuildOptions&) {
            return cts::build_zst_dme(sinks);
        },
        false, true, false);
}
