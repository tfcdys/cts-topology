// SALT and SLLT (BST + SALT) topology builders.
//
// SALT reference:
//   https://github.com/OSCC-Project/iEDA/tree/master/src/operation/iCTS
//   https://github.com/chengengjie/salt
//
// SLLT is inspired by iCTS `bstSaltTree`: build a bounded-skew tree, convert it
// to a SALT tree, then apply the SALT/CBS shallow-light refinement.

#include "salt_common.hpp"

#include "bst_dme.hpp"
#include "salt/base/rsa.h"
#include "salt/refine/refine.h"
#include "salt/salt.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cts {
namespace {

constexpr int kSaltScale = 1000;  // 1 nm integer grid
constexpr double kSaltEps = 1e-9;

struct ScaledInput {
    std::vector<salt::Point> sink_points;
    salt::Point source;
};

int scaled_round(double value) {
    return static_cast<int>(std::llround(value * kSaltScale));
}


bool all_collinear(const std::vector<Sink>& sinks) {
    if (sinks.size() <= 2) return true;
    const Sink& a = sinks[0];
    std::size_t second = 1;
    while (second < sinks.size() &&
           std::abs(sinks[second].x - a.x) <= kSaltEps &&
           std::abs(sinks[second].y - a.y) <= kSaltEps) {
        ++second;
    }
    if (second >= sinks.size()) return true;
    const double dx = sinks[second].x - a.x;
    const double dy = sinks[second].y - a.y;
    for (std::size_t i = second + 1; i < sinks.size(); ++i) {
        const double cross =
            dx * (sinks[i].y - a.y) - dy * (sinks[i].x - a.x);
        if (std::abs(cross) > 1e-7 * std::max(1.0, std::abs(dx) + std::abs(dy))) {
            return false;
        }
    }
    return true;
}

ScaledInput make_scaled_input(const std::vector<Sink>& sinks,
                              const salt::Point& desired_source) {
    ScaledInput input;
    input.sink_points.reserve(sinks.size());
    std::set<std::pair<int, int>> occupied;
    auto unique_point = [&](int x, int y) {
        while (occupied.count({x, y}) != 0) {
            ++x;
        }
        occupied.insert({x, y});
        return salt::Point(x, y);
    };
    for (const Sink& sink : sinks) {
        input.sink_points.push_back(
            unique_point(scaled_round(sink.x), scaled_round(sink.y)));
    }
    input.source = unique_point(desired_source.x, desired_source.y);
    return input;
}

salt::Point bbox_center(const std::vector<Sink>& sinks) {
    if (sinks.empty()) {
        return salt::Point(0, 0);
    }
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymin = std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();
    for (const Sink& sink : sinks) {
        xmin = std::min(xmin, sink.x);
        xmax = std::max(xmax, sink.x);
        ymin = std::min(ymin, sink.y);
        ymax = std::max(ymax, sink.y);
    }
    return salt::Point(scaled_round(0.5 * (xmin + xmax)),
                       scaled_round(0.5 * (ymin + ymax)));
}

std::vector<std::shared_ptr<salt::Pin>> make_pins(
    const std::vector<Sink>& sinks, const ScaledInput& input) {
    std::vector<std::shared_ptr<salt::Pin>> pins;
    pins.reserve(sinks.size() + 1);
    pins.push_back(std::make_shared<salt::Pin>(input.source, 0, 1.0));
    for (std::size_t i = 0; i < sinks.size(); ++i) {
        pins.push_back(std::make_shared<salt::Pin>(
            input.sink_points[i], static_cast<int>(i) + 1,
            std::max(1e-9, sinks[i].cap_ff)));
    }
    return pins;
}

std::string salt_metadata(const std::string& variant, double eps,
                          int refine_level, const salt::Point& source) {
    std::ostringstream out;
    out << "{\"variant\":\"" << variant << "\","
        << "\"library\":\"SALT (MIT) + FLUTE3 (BSD)\","
        << "\"source\":\"bbox-center-or-bst-root\","
        << "\"eps\":" << format_number(eps) << ","
        << "\"refine_level\":" << refine_level << ","
        << "\"scale_nm_per_um\":" << kSaltScale << ","
        << "\"source_x_um\":" << format_number(
               static_cast<double>(source.x) / kSaltScale)
        << ",\"source_y_um\":" << format_number(
               static_cast<double>(source.y) / kSaltScale)
        << "}";
    return out.str();
}

Topology convert_salt_tree(const salt::Tree& tree,
                           const std::vector<Sink>& sinks,
                           const std::string& metadata) {
    Topology topology;
    if (!tree.source) {
        if (!sinks.empty()) {
            fail("SALT returned an empty tree");
        }
        topology.metadata_json = metadata;
        return topology;
    }

    auto node_for_sink = [&](const salt::TreeNode* node, int& sink_index) {
        sink_index = -1;
        if (node->pin && node->pin->id >= 1 &&
            static_cast<std::size_t>(node->pin->id) <= sinks.size()) {
            sink_index = node->pin->id - 1;
            return true;
        }
        return false;
    };

    auto add_root = [&](const salt::TreeNode* node) {
        Node result;
        result.id = static_cast<int>(topology.nodes.size());
        result.kind = "root";
        result.x = static_cast<double>(node->loc.x) / kSaltScale;
        result.y = static_cast<double>(node->loc.y) / kSaltScale;
        topology.root = result.id;
        ++topology.internal_nodes;
        topology.nodes.push_back(result);
        return result.id;
    };

    auto add_sink = [&](int sink_index) {
        const Sink& sink = sinks[static_cast<std::size_t>(sink_index)];
        Node result;
        result.id = static_cast<int>(topology.nodes.size());
        result.kind = "sink";
        result.sink_id = sink.id;
        result.x = sink.x;
        result.y = sink.y;
        result.cap_ff = sink.cap_ff;
        ++topology.sink_nodes;
        topology.nodes.push_back(result);
        return result.id;
    };

    auto add_internal = [&](const salt::TreeNode* node) {
        Node result;
        result.id = static_cast<int>(topology.nodes.size());
        result.kind = "internal";
        result.x = static_cast<double>(node->loc.x) / kSaltScale;
        result.y = static_cast<double>(node->loc.y) / kSaltScale;
        ++topology.internal_nodes;
        topology.nodes.push_back(result);
        return result.id;
    };

    auto add_edge = [&](int parent, int child) {
        topology.edges.push_back(Edge{parent, child, "HV"});
    };

    std::function<void(const salt::TreeNode*, int)> visit =
        [&](const salt::TreeNode* node, int parent) {
            int sink_index = -1;
            if (node_for_sink(node, sink_index)) {
                const int sink_id = add_sink(sink_index);
                add_edge(parent, sink_id);
                if (!node->children.empty()) {
                    // A SALT tree may route through a sink pin. Keep the sink
                    // as a leaf and continue the routing through a coincident
                    // internal clone.
                    const int clone = add_internal(node);
                    add_edge(parent, clone);
                    for (const auto& child : node->children) {
                        visit(child.get(), clone);
                    }
                }
                return;
            }
            const int current = add_internal(node);
            add_edge(parent, current);
            for (const auto& child : node->children) {
                visit(child.get(), current);
            }
        };

    const int root = add_root(tree.source.get());
    for (const auto& child : tree.source->children) {
        visit(child.get(), root);
    }

    if (static_cast<int>(topology.edges.size()) !=
        static_cast<int>(topology.nodes.size()) - 1) {
        fail("SALT output is not a tree");
    }
    // SALT/FLUTE may collapse coincident pins (duplicate sink coordinates).
    // Reattach missing sink IDs to the nearest non-sink routing node.
    std::set<int> present_sink_ids;
    for (const auto& node : topology.nodes) {
        if (node.kind == "sink") {
            present_sink_ids.insert(node.sink_id);
        }
    }
    for (std::size_t sink_index = 0; sink_index < sinks.size(); ++sink_index) {
        const Sink& sink = sinks[sink_index];
        if (present_sink_ids.count(sink.id) != 0) {
            continue;
        }
        int best = -1;
        double best_distance = std::numeric_limits<double>::infinity();
        for (const auto& node : topology.nodes) {
            if (node.kind == "sink") {
                continue;
            }
            const double distance =
                std::abs(node.x - sink.x) + std::abs(node.y - sink.y);
            if (distance < best_distance - kSaltEps) {
                best_distance = distance;
                best = node.id;
            }
        }
        if (best < 0) {
            fail("SALT could not attach a dropped sink pin");
        }
        const int sink_id = add_sink(static_cast<int>(sink_index));
        add_edge(best, sink_id);
    }
    if (topology.sink_nodes != static_cast<int>(sinks.size())) {
        fail("SALT output lost sink pins");
    }
    topology.metadata_json = metadata;
    return topology;
}

// Adapted from iCTS TreeSaltBuilder/CBS (Mulan PSL v2), using only the SALT
// public API. It fuses a pre-built bounded-skew tree with SALT refinement.
class TreeSaltBuilder {
public:
    void run(const salt::Net& net, salt::Tree& input_tree, double eps,
             int refine_level = 3) {
        auto tree = input_tree;
        if (refine_level >= 1) {
            salt::Refine::flip(tree);
            salt::Refine::uShift(tree);
        }
        init(tree, net.source());
        dfs(tree.source, source_, eps);
        finalize(net, input_tree);
        input_tree.RemoveTopoRedundantSteiner();
        salt::RsaBuilder rsa_builder;
        rsa_builder.ReplaceRootChildren(input_tree);
        input_tree.UpdateId();
        if (refine_level >= 1) {
            salt::Refine::cancelIntersect(input_tree);
            salt::Refine::flip(input_tree);
            salt::Refine::uShift(input_tree);
            if (refine_level >= 2) {
                salt::Refine::substitute(input_tree, eps, refine_level == 3);
            }
        }
    }

private:
    void init(salt::Tree& min_tree, std::shared_ptr<salt::Pin> source_pin) {
        min_tree.UpdateId();
        const auto nodes = min_tree.ObtainNodes();
        nodes_.resize(nodes.size());
        shortest_.resize(nodes.size());
        current_.resize(nodes.size());
        for (const auto& node : nodes) {
            nodes_[static_cast<std::size_t>(node->id)] =
                std::make_shared<salt::TreeNode>(node->loc, node->pin,
                                                 node->id);
            shortest_[static_cast<std::size_t>(node->id)] =
                utils::Dist(source_pin->loc, node->loc);
            current_[static_cast<std::size_t>(node->id)] =
                std::numeric_limits<double>::max();
        }
        current_[static_cast<std::size_t>(source_pin->id)] = 0.0;
        source_ = nodes_[static_cast<std::size_t>(source_pin->id)];
    }

    void finalize(const salt::Net& net, salt::Tree& tree) {
        for (const auto& node : nodes_) {
            if (node->parent) {
                nodes_[static_cast<std::size_t>(node->parent->id)]
                    ->children.push_back(node);
            }
        }
        tree.source = source_;
        tree.net = &net;
    }

    bool relax(const std::shared_ptr<salt::TreeNode>& u,
               const std::shared_ptr<salt::TreeNode>& v) {
        const double new_latency =
            current_[static_cast<std::size_t>(u->id)] +
            utils::Dist(u->loc, v->loc);
        if (current_[static_cast<std::size_t>(v->id)] > new_latency) {
            current_[static_cast<std::size_t>(v->id)] = new_latency;
            v->parent = u;
            return true;
        }
        if (std::abs(current_[static_cast<std::size_t>(v->id)] -
                     new_latency) <= kSaltEps &&
            utils::Dist(u->loc, v->loc) < v->WireToParentChecked()) {
            v->parent = u;
            return true;
        }
        return false;
    }

    void dfs(const std::shared_ptr<salt::TreeNode>& tree_node,
             const std::shared_ptr<salt::TreeNode>& cbs_node, double eps) {
        if (tree_node->pin &&
            current_[static_cast<std::size_t>(cbs_node->id)] >
                (1.0 + eps) *
                    shortest_[static_cast<std::size_t>(cbs_node->id)]) {
            cbs_node->parent = source_;
            current_[static_cast<std::size_t>(cbs_node->id)] =
                shortest_[static_cast<std::size_t>(cbs_node->id)];
        }
        for (const auto& child : tree_node->children) {
            relax(cbs_node, nodes_[static_cast<std::size_t>(child->id)]);
            dfs(child, nodes_[static_cast<std::size_t>(child->id)], eps);
            relax(nodes_[static_cast<std::size_t>(child->id)], cbs_node);
        }
    }

    std::vector<std::shared_ptr<salt::TreeNode>> nodes_;
    std::vector<double> shortest_;
    std::vector<double> current_;
    std::shared_ptr<salt::TreeNode> source_;
};

void build_net_and_pins(const std::vector<Sink>& sinks,
                        const ScaledInput& input, salt::Net& net) {
    net.init(0, "cts", make_pins(sinks, input));
}

salt::Tree build_sllt_base_tree(const std::vector<Sink>& sinks,
                                double skew_bound_um,
                                const ScaledInput& input, salt::Net& net) {
    const Topology base = build_bst_dme(sinks, skew_bound_um);
    std::vector<std::shared_ptr<salt::Pin>> pins =
        make_pins(sinks, input);
    if (base.root >= 0 &&
        static_cast<std::size_t>(base.root) < base.nodes.size()) {
        const Node& root_node =
            base.nodes[static_cast<std::size_t>(base.root)];
        pins[0]->loc = salt::Point(scaled_round(root_node.x),
                                   scaled_round(root_node.y));
    }
    net.init(0, "sllt", pins);
    std::unordered_map<int, std::shared_ptr<salt::TreeNode>> nodes;
    for (const auto& node : base.nodes) {
        salt::Point loc(scaled_round(node.x), scaled_round(node.y));
        std::shared_ptr<salt::Pin> pin;
        if (node.kind == "sink") {
            for (std::size_t i = 0; i < sinks.size(); ++i) {
                if (sinks[i].id == node.sink_id) {
                    pin = pins[i + 1];
                    break;
                }
            }
            if (!pin) {
                fail("SLLT: sink mapping failed");
            }
        } else if (node.id == base.root) {
            pin = pins[0];
        }
        nodes[node.id] = std::make_shared<salt::TreeNode>(loc, pin);
    }
    std::vector<std::vector<int>> children(base.nodes.size());
    for (const auto& edge : base.edges) {
        children[static_cast<std::size_t>(edge.source)].push_back(edge.target);
    }
    std::vector<int> stack{base.root};
    while (!stack.empty()) {
        const int id = stack.back();
        stack.pop_back();
        for (const int child : children[static_cast<std::size_t>(id)]) {
            nodes[static_cast<std::size_t>(child)]->parent =
                nodes[static_cast<std::size_t>(id)];
            nodes[static_cast<std::size_t>(id)]->children.push_back(
                nodes[static_cast<std::size_t>(child)]);
            stack.push_back(child);
        }
    }
    salt::Tree tree(nodes[static_cast<std::size_t>(base.root)], &net);
    tree.UpdateId();
    return tree;
}

}  // namespace

Topology build_salt_tree(const std::vector<Sink>& sinks,
                         const SaltOptions& options) {
    if (sinks.empty()) {
        Topology topology;
        topology.metadata_json =
            salt_metadata("SALT", options.eps, options.refine_level,
                          salt::Point(0, 0));
        return topology;
    }
    const ScaledInput input = make_scaled_input(sinks, bbox_center(sinks));
    salt::Net net;
    build_net_and_pins(sinks, input, net);
    salt::Tree tree;
    salt::SaltBuilder builder;
    builder.Run(net, tree, options.eps, options.refine_level);
    return convert_salt_tree(
        tree, sinks,
        salt_metadata("SALT", options.eps, options.refine_level,
                      input.source));
}

Topology build_sllt_tree(const std::vector<Sink>& sinks,
                         const SlltOptions& options) {
    if (sinks.empty()) {
        Topology topology;
        topology.metadata_json =
            salt_metadata("SLLT (BST + SALT)", options.eps,
                          options.refine_level, salt::Point(0, 0));
        return topology;
    }
    if (sinks.size() <= 2 || all_collinear(sinks)) {
        // SALT's tree refinement expects a non-degenerate multi-level tree.
        return build_bst_dme(sinks, options.skew_bound_um);
    }
    salt::Net net;
    const ScaledInput input =
        make_scaled_input(sinks, bbox_center(sinks));
    salt::Tree tree = build_sllt_base_tree(sinks, options.skew_bound_um,
                                           input, net);
    TreeSaltBuilder builder;
    builder.run(net, tree, options.eps, options.refine_level);
    return convert_salt_tree(
        tree, sinks,
        salt_metadata("SLLT (BST + SALT)", options.eps,
                      options.refine_level, input.source));
}

}  // namespace cts
