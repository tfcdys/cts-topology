// GH-tree: generalized H-tree topology inspired by TritonCTS.
//
// Reference:
//   https://github.com/The-OpenROAD-Project/TritonCTS
//   K. Han, A. B. Kahng, J. Li, "Optimal Generalized H-Tree Topology and
//   Buffering for High-Performance and Low-Power Clock Distribution".
//
// This is a simplified topology-only implementation: deterministic capacitated
// K-means clustering into four child regions, recursive H-shaped connections,
// no buffer insertion / RC optimization / blockage handling.

#include "cts_common.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace cts {
namespace {

constexpr int kKMeansIterations = 8;
constexpr double kGhEps = 1e-9;

struct Point {
    double x = 0.0;
    double y = 0.0;
};

struct Cluster {
    std::vector<int> members;
    double x = 0.0;
    double y = 0.0;
    double cap = 0.0;
};

double l1(const Point& a, const Point& b) {
    return std::abs(a.x - b.x) + std::abs(a.y - b.y);
}

struct Bounds {
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymin = std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();
};

Bounds bounds_of(const std::vector<Sink>& sinks, const std::vector<int>& ids) {
    Bounds b;
    for (const int i : ids) {
        const Sink& s = sinks[static_cast<std::size_t>(i)];
        b.xmin = std::min(b.xmin, s.x);
        b.xmax = std::max(b.xmax, s.x);
        b.ymin = std::min(b.ymin, s.y);
        b.ymax = std::max(b.ymax, s.y);
    }
    return b;
}

class GhTreeBuilder {
public:
    GhTreeBuilder(const std::vector<Sink>& sinks, bool cap_aware)
        : sinks_(sinks), cap_aware_(cap_aware) {}

    Topology build() {
        topology_.nodes.clear();
        topology_.edges.clear();
        if (!sinks_.empty()) {
            std::vector<int> ids(sinks_.size());
            std::iota(ids.begin(), ids.end(), 0);
            topology_.root = build_group(ids, true);
        }
        if (!sinks_.empty() &&
            static_cast<int>(topology_.edges.size()) !=
                static_cast<int>(topology_.nodes.size()) - 1) {
            fail("GH-tree internal error: output is not a tree");
        }
        topology_.metadata_json =
            "{\"variant\":\"GH-tree\",\"clustering\":\"capacitated K-means "
            "(4-way, deterministic)\",\"kmeans_iterations\":" +
            std::to_string(kKMeansIterations) + ",\"capacitance_aware\":" +
            (cap_aware_ ? "true" : "false") + "}";
        return std::move(topology_);
    }

private:
    const std::vector<Sink>& sinks_;
    bool cap_aware_ = false;
    Topology topology_;

    double weight(int index) const {
        return cap_aware_ ? std::max(1e-9, sinks_[static_cast<std::size_t>(index)].cap_ff)
                          : 1.0;
    }

    Point weighted_centroid(const std::vector<int>& ids) const {
        double total = 0.0;
        double x = 0.0;
        double y = 0.0;
        for (const int i : ids) {
            const Sink& s = sinks_[static_cast<std::size_t>(i)];
            const double w = weight(i);
            total += w;
            x += w * s.x;
            y += w * s.y;
        }
        if (total <= 0.0) {
            const Bounds b = bounds_of(sinks_, ids);
            return Point{0.5 * (b.xmin + b.xmax), 0.5 * (b.ymin + b.ymax)};
        }
        return Point{x / total, y / total};
    }

    int add_node(const Point& pos, const std::string& kind,
                 int sink_index = -1) {
        Node node;
        node.id = static_cast<int>(topology_.nodes.size());
        node.x = pos.x;
        node.y = pos.y;
        node.kind = kind;
        if (sink_index >= 0) {
            const Sink& sink = sinks_[static_cast<std::size_t>(sink_index)];
            node.sink_id = sink.id;
            node.cap_ff = sink.cap_ff;
            ++topology_.sink_nodes;
        } else if (kind == "bend") {
            ++topology_.bend_nodes;
        } else {
            ++topology_.internal_nodes;
        }
        topology_.nodes.push_back(node);
        return node.id;
    }

    void add_edge(int source, int target, const std::string& route) {
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

    std::array<std::vector<int>, 4> initial_quadrants(
        const std::vector<int>& ids) const {
        std::array<std::vector<int>, 4> groups;
        std::vector<int> by_x = ids;
        std::sort(by_x.begin(), by_x.end(), [&](int a, int b) {
            return compare_x(sinks_[static_cast<std::size_t>(a)],
                             sinks_[static_cast<std::size_t>(b)]);
        });
        const std::size_t half = by_x.size() / 2;
        std::vector<int> left(by_x.begin(), by_x.begin() + half);
        std::vector<int> right(by_x.begin() + half, by_x.end());
        auto split = [&](std::vector<int>& side, std::vector<int>& bottom,
                         std::vector<int>& top) {
            std::sort(side.begin(), side.end(), [&](int a, int b) {
                return compare_y(sinks_[static_cast<std::size_t>(a)],
                                 sinks_[static_cast<std::size_t>(b)]);
            });
            const std::size_t mid = side.size() / 2;
            bottom.assign(side.begin(), side.begin() + mid);
            top.assign(side.begin() + mid, side.end());
        };
        split(left, groups[0], groups[1]);
        split(right, groups[2], groups[3]);
        return groups;
    }

    std::array<std::vector<int>, 4> kmeans_quadrants(
        const std::vector<int>& ids) {
        auto groups = initial_quadrants(ids);
        std::array<Point, 4> centers;
        for (int g = 0; g < 4; ++g) {
            if (groups[static_cast<std::size_t>(g)].empty()) {
                return groups;
            }
            centers[static_cast<std::size_t>(g)] =
                weighted_centroid(groups[static_cast<std::size_t>(g)]);
        }

        std::array<std::vector<int>, 4> previous;
        for (int iteration = 0; iteration < kKMeansIterations; ++iteration) {
            std::array<std::vector<int>, 4> assignment;
            for (const int id : ids) {
                const Point p{sinks_[static_cast<std::size_t>(id)].x,
                              sinks_[static_cast<std::size_t>(id)].y};
                int best = 0;
                double best_distance = l1(p, centers[0]);
                for (int g = 1; g < 4; ++g) {
                    const double distance = l1(p, centers[static_cast<std::size_t>(g)]);
                    if (distance < best_distance - kGhEps) {
                        best_distance = distance;
                        best = g;
                    }
                }
                assignment[static_cast<std::size_t>(best)].push_back(id);
            }
            bool empty = false;
            for (const auto& group : assignment) {
                if (group.empty()) empty = true;
            }
            if (empty) {
                return groups;
            }
            bool unchanged = true;
            for (int g = 0; g < 4; ++g) {
                std::vector<int> a = previous[static_cast<std::size_t>(g)];
                std::vector<int> b = assignment[static_cast<std::size_t>(g)];
                std::sort(a.begin(), a.end());
                std::sort(b.begin(), b.end());
                if (a != b) unchanged = false;
            }
            groups = assignment;
            for (int g = 0; g < 4; ++g) {
                centers[static_cast<std::size_t>(g)] =
                    weighted_centroid(groups[static_cast<std::size_t>(g)]);
            }
            if (unchanged) {
                break;
            }
            previous = groups;
        }

        std::size_t largest = 0;
        for (const auto& group : groups) {
            largest = std::max(largest, group.size());
        }
        // Refined K-means can collapse on degenerate inputs; keep the balanced
        // quadrisection fallback when it does.
        if (largest * 4 > ids.size() * 3) {
            return initial_quadrants(ids);
        }
        return groups;
    }

    void connect_side(int parent, const Point& parent_pos,
                      const std::vector<int>& children, bool is_top) {
        if (children.empty()) {
            return;
        }
        double average_y = 0.0;
        for (const int child : children) {
            average_y += topology_.nodes[static_cast<std::size_t>(child)].y;
        }
        average_y /= static_cast<double>(children.size());
        const double bar_y = is_top ? std::max(parent_pos.y, average_y)
                                    : std::min(parent_pos.y, average_y);
        if (std::abs(bar_y - parent_pos.y) < 1e-12) {
            for (const int child : children) {
                add_edge(parent, child, "HV");
            }
            return;
        }
        const int bend =
            add_node(Point{parent_pos.x, bar_y}, "bend");
        add_edge(parent, bend, "V");
        for (const int child : children) {
            add_edge(bend, child, "HV");
        }
    }

    int build_group(const std::vector<int>& ids, bool is_root) {
        if (ids.size() == 1) {
            const int index = ids.front();
            const Sink& sink = sinks_[static_cast<std::size_t>(index)];
            return add_node(Point{sink.x, sink.y}, "sink", index);
        }

        const Bounds bounds = bounds_of(sinks_, ids);
        const Point root_pos{0.5 * (bounds.xmin + bounds.xmax),
                             0.5 * (bounds.ymin + bounds.ymax)};
        const int parent =
            add_node(root_pos, is_root ? "root" : "internal");

        auto groups = kmeans_quadrants(ids);
        std::vector<int> children;
        for (auto& group : groups) {
            if (group.empty()) {
                continue;
            }
            children.push_back(build_group(group, false));
        }

        std::vector<int> top;
        std::vector<int> bottom;
        for (const int child : children) {
            if (topology_.nodes[static_cast<std::size_t>(child)].y >=
                root_pos.y) {
                top.push_back(child);
            } else {
                bottom.push_back(child);
            }
        }
        connect_side(parent, root_pos, top, true);
        connect_side(parent, root_pos, bottom, false);
        return parent;
    }
};

}  // namespace

Topology build_gh_tree(const std::vector<Sink>& sinks, bool cap_aware) {
    GhTreeBuilder builder(sinks, cap_aware);
    return builder.build();
}

}  // namespace cts

#ifndef GH_TREE_NO_MAIN
int main(int argc, char** argv) {
    return cts::run_topology_cli(
        argc, argv, "gh_tree",
        [](const std::vector<cts::Sink>& sinks,
           const cts::BuildOptions& options) {
            return cts::build_gh_tree(sinks, options.capacitance_aware);
        },
        true, false, false);
}
#endif
