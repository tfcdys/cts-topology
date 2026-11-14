#include "cts_common.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace cts {

Bounds bounds_of(const std::vector<Sink>& sinks, const std::vector<int>& ids) {
    Bounds b;
    b.xmin = b.ymin = std::numeric_limits<double>::infinity();
    b.xmax = b.ymax = -std::numeric_limits<double>::infinity();
    for (const int index : ids) {
        const Sink& s = sinks[static_cast<std::size_t>(index)];
        b.xmin = std::min(b.xmin, s.x);
        b.xmax = std::max(b.xmax, s.x);
        b.ymin = std::min(b.ymin, s.y);
        b.ymax = std::max(b.ymax, s.y);
    }
    return b;
}

class HTreeBuilder {
public:
    explicit HTreeBuilder(const std::vector<Sink>& sinks,
                          bool capacitance_aware = false)
        : sinks_(sinks), capacitance_aware_(capacitance_aware) {
        sink_node_ids_.assign(sinks_.size(), -1);
    }

    Topology build() {
        std::vector<int> ids(sinks_.size());
        std::iota(ids.begin(), ids.end(), 0);
        topology_.root = build_group(ids, true);

        if (topology_.sink_nodes != static_cast<int>(sinks_.size())) {
            fail("internal error: not every sink was inserted into the H-tree");
        }

        // Tree invariant sanity check: a rooted tree with n nodes has n-1 edges.
        const int expected_edges =
            static_cast<int>(topology_.nodes.size()) - 1;
        if (static_cast<int>(topology_.edges.size()) != expected_edges) {
            fail("internal error: H-tree is not a tree");
        }
        return std::move(topology_);
    }

private:
    const std::vector<Sink>& sinks_;
    Topology topology_;
    std::vector<int> sink_node_ids_;
    bool capacitance_aware_ = false;

    // In uniform mode every sink has the same effective capacitance weight,
    // so the partition degenerates to the count-balanced geometric split.
    double cap_weight(int sink_index) const {
        if (!capacitance_aware_) {
            return 1.0;
        }
        return sinks_[static_cast<std::size_t>(sink_index)].cap_ff;
    }

    // Choose how many leading elements of a coordinate-sorted list form the
    // low side.  In cap-aware mode, stay within 40%-60% of the count and pick
    // the split with the smallest total-capacitance imbalance.
    std::size_t choose_split(const std::vector<int>& sorted) const {
        const std::size_t n = sorted.size();
        if (!capacitance_aware_ || n < 5) {
            return n / 2;
        }
        const std::size_t low = std::max<std::size_t>(
            1, static_cast<std::size_t>(std::ceil(0.40 * n)));
        const std::size_t high = std::min<std::size_t>(
            n - 1, static_cast<std::size_t>(std::floor(0.60 * n)));
        if (low > high) {
            return n / 2;
        }

        double total_cap = 0.0;
        for (const int index : sorted) {
            total_cap += cap_weight(index);
        }
        if (total_cap <= 0.0) {
            return n / 2;
        }

        double running_cap = 0.0;
        double best_score = std::numeric_limits<double>::infinity();
        std::size_t best_split = n / 2;
        for (std::size_t i = 0; i < n; ++i) {
            running_cap += cap_weight(sorted[i]);
            const std::size_t split = i + 1;
            if (split < low || split > high) {
                continue;
            }
            const double cap_imbalance =
                std::abs(2.0 * running_cap - total_cap) / total_cap;
            const double count_deviation =
                std::abs(static_cast<double>(split) -
                         static_cast<double>(n) / 2.0) /
                (static_cast<double>(n) / 2.0);
            const double score = cap_imbalance + 0.10 * count_deviation;
            if (score < best_score) {
                best_score = score;
                best_split = split;
            }
        }
        return best_split;
    }

    int add_node(double x, double y, const std::string& kind,
                 int sink_index = -1) {
        Node node;
        node.id = static_cast<int>(topology_.nodes.size());
        node.x = x;
        node.y = y;
        node.kind = kind;
        if (sink_index >= 0) {
            const Sink& sink = sinks_[static_cast<std::size_t>(sink_index)];
            node.sink_id = sink.id;
            node.cap_ff = sink.cap_ff;
            ++topology_.sink_nodes;
            sink_node_ids_[static_cast<std::size_t>(sink_index)] = node.id;
        } else if (kind == "bend") {
            ++topology_.bend_nodes;
        } else if (kind == "internal" || kind == "root") {
            ++topology_.internal_nodes;
        }
        topology_.nodes.push_back(node);
        return node.id;
    }

    void add_edge(int source, int target, const std::string& route) {
        if (source == target) {
            fail("internal error: self edge in H-tree");
        }
        topology_.edges.push_back(Edge{source, target, route});
    }

    static bool compare_x(const Sink& a, const Sink& b) {
        if (a.x != b.x) return a.x < b.x;
        if (a.y != b.y) return a.y < b.y;
        return a.id < b.id;
    }

    static bool compare_y(const Sink& a, const Sink& b) {
        if (a.y != b.y) return a.y < b.y;
        if (a.x != b.x) return a.x < b.x;
        return a.id < b.id;
    }

    // Split a group into four roughly equal groups:
    //   [0] left/bottom, [1] left/top, [2] right/bottom, [3] right/top.
    std::array<std::vector<int>, 4> quadrisect(
        const std::vector<int>& ids) const {
        std::array<std::vector<int>, 4> groups;

        std::vector<int> by_x = ids;
        std::sort(by_x.begin(), by_x.end(), [this](int a, int b) {
            return compare_x(sinks_[static_cast<std::size_t>(a)],
                             sinks_[static_cast<std::size_t>(b)]);
        });

        const std::size_t half = choose_split(by_x);
        std::vector<int> left(by_x.begin(), by_x.begin() + half);
        std::vector<int> right(by_x.begin() + half, by_x.end());

        auto split_y = [this](std::vector<int>& side, std::vector<int>& bottom,
                              std::vector<int>& top) {
            std::sort(side.begin(), side.end(), [this](int a, int b) {
                return compare_y(sinks_[static_cast<std::size_t>(a)],
                                 sinks_[static_cast<std::size_t>(b)]);
            });
            const std::size_t mid = choose_split(side);
            bottom.assign(side.begin(), side.begin() + mid);
            top.assign(side.begin() + mid, side.end());
        };

        std::vector<int> left_bottom;
        std::vector<int> left_top;
        std::vector<int> right_bottom;
        std::vector<int> right_top;
        split_y(left, left_bottom, left_top);
        split_y(right, right_bottom, right_top);

        groups[0] = std::move(left_bottom);
        groups[1] = std::move(left_top);
        groups[2] = std::move(right_bottom);
        groups[3] = std::move(right_top);
        return groups;
    }

    int build_group(const std::vector<int>& ids, bool is_root) {
        if (ids.size() == 1) {
            const int sink_index = ids.front();
            const Sink& sink = sinks_[static_cast<std::size_t>(sink_index)];
            return add_node(sink.x, sink.y, "sink", sink_index);
        }

        const Bounds b = bounds_of(sinks_, ids);
        const double center_x = 0.5 * (b.xmin + b.xmax);
        const double center_y = 0.5 * (b.ymin + b.ymax);
        const int parent = add_node(center_x, center_y,
                                    is_root ? "root" : "internal");

        const auto groups = quadrisect(ids);

        std::vector<int> top_children;
        std::vector<int> bottom_children;
        std::array<int, 4> child_node{-1, -1, -1, -1};
        for (int g = 0; g < 4; ++g) {
            if (groups[static_cast<std::size_t>(g)].empty()) {
                continue;
            }
            const int child = build_group(groups[static_cast<std::size_t>(g)],
                                          false);
            child_node[static_cast<std::size_t>(g)] = child;
            if (g == 1 || g == 3) {
                top_children.push_back(child);
            } else {
                bottom_children.push_back(child);
            }
        }

        connect_side(parent, center_x, center_y, top_children, true);
        connect_side(parent, center_x, center_y, bottom_children, false);
        return parent;
    }

    void connect_side(int parent, double parent_x, double parent_y,
                      const std::vector<int>& children, bool is_top) {
        if (children.empty()) {
            return;
        }

        double average_y = 0.0;
        for (const int child : children) {
            average_y += topology_.nodes[static_cast<std::size_t>(child)].y;
        }
        average_y /= static_cast<double>(children.size());

        // Keep the bar on the same side of the parent as the child group.
        const double bar_y = is_top ? std::max(parent_y, average_y)
                                    : std::min(parent_y, average_y);

        if (std::abs(bar_y - parent_y) < 1e-12) {
            for (const int child : children) {
                add_edge(parent, child, "HV");
            }
            return;
        }

        const int bar = add_node(parent_x, bar_y, "bend");
        add_edge(parent, bar, "V");
        for (const int child : children) {
            add_edge(bar, child, "HV");
        }
    }
};

Topology build_htree(const std::vector<Sink>& sinks, bool capacitance_aware) {
    HTreeBuilder builder(sinks, capacitance_aware);
    return builder.build();
}

}  // namespace cts

int main(int argc, char** argv) {
    return cts::run_topology_cli(
        argc, argv, "htree",
        [](const std::vector<cts::Sink>& sinks, const cts::BuildOptions& options) {
            return cts::build_htree(sinks, options.capacitance_aware);
        },
        true, false, false);
}
