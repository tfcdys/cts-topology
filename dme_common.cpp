#include "dme_common.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <tuple>
#include <utility>

namespace cts {

DmePoint dme_point(double x, double y) {
    return DmePoint{x, y};
}

DmePoint dme_midpoint(const DmePoint& a, const DmePoint& b) {
    return DmePoint{0.5 * (a.x + b.x), 0.5 * (a.y + b.y)};
}

double dme_l1(const DmePoint& a, const DmePoint& b) {
    return std::abs(a.x - b.x) + std::abs(a.y - b.y);
}

bool dme_same_point(const DmePoint& a, const DmePoint& b) {
    return dme_l1(a, b) <= kDmeEps;
}

DmeSegment dme_make_segment(const DmePoint& a, const DmePoint& b) {
    DmeSegment s;
    s.a = a;
    s.b = b;
    s.valid = true;
    return s;
}

DmeSegment dme_make_point_segment(const DmePoint& p) {
    return dme_make_segment(p, p);
}

bool dme_segment_is_point(const DmeSegment& s) {
    return dme_same_point(s.a, s.b);
}

double dme_cross(const DmePoint& a, const DmePoint& b, const DmePoint& c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

DmePoint dme_lerp(const DmePoint& a, const DmePoint& b, double t) {
    return DmePoint{a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t};
}

DmePoint dme_closest_point_on_segment(const DmePoint& p,
                                      const DmeSegment& s) {
    if (!s.valid) {
        return p;
    }
    const double dx = s.b.x - s.a.x;
    const double dy = s.b.y - s.a.y;
    double best_distance = std::numeric_limits<double>::infinity();
    DmePoint best = s.a;

    auto consider = [&](double t) {
        t = std::clamp(t, 0.0, 1.0);
        const DmePoint q = dme_lerp(s.a, s.b, t);
        const double distance = dme_l1(p, q);
        if (distance < best_distance) {
            best_distance = distance;
            best = q;
        }
    };

    consider(0.0);
    consider(1.0);
    if (std::abs(dx) > kDmeEps) {
        consider((p.x - s.a.x) / dx);
    }
    if (std::abs(dy) > kDmeEps) {
        consider((p.y - s.a.y) / dy);
    }
    return best;
}

double dme_point_segment_distance(const DmePoint& p, const DmeSegment& s) {
    return dme_l1(p, dme_closest_point_on_segment(p, s));
}

bool dme_point_on_segment(const DmePoint& p, const DmeSegment& s) {
    return s.valid && dme_point_segment_distance(p, s) <= 1e-6;
}

bool dme_segment_intersection(const DmeSegment& s1, const DmeSegment& s2,
                              DmeSegment& out) {
    out.valid = false;
    if (!s1.valid || !s2.valid) {
        return false;
    }
    if (dme_segment_is_point(s1) && dme_segment_is_point(s2)) {
        if (dme_same_point(s1.a, s2.a)) {
            out = dme_make_point_segment(s1.a);
            return true;
        }
        return false;
    }
    if (dme_segment_is_point(s1)) {
        if (dme_point_on_segment(s1.a, s2)) {
            out = dme_make_point_segment(s1.a);
            return true;
        }
        return false;
    }
    if (dme_segment_is_point(s2)) {
        return dme_segment_intersection(s2, s1, out);
    }

    const DmePoint p = s1.a;
    const DmePoint r{s1.b.x - s1.a.x, s1.b.y - s1.a.y};
    const DmePoint q = s2.a;
    const DmePoint s{s2.b.x - s2.a.x, s2.b.y - s2.a.y};
    const double rxs = r.x * s.y - r.y * s.x;
    const double qpxr = (q.x - p.x) * r.y - (q.y - p.y) * r.x;

    if (std::abs(rxs) <= kDmeEps) {
        if (std::abs(qpxr) > 1e-6) {
            return false;
        }
        // Collinear: project onto the dominant axis and intersect intervals.
        const bool use_x = std::abs(r.x) >= std::abs(r.y);
        const double a0 = use_x ? s1.a.x : s1.a.y;
        const double a1 = use_x ? s1.b.x : s1.b.y;
        const double b0 = use_x ? s2.a.x : s2.a.y;
        const double b1 = use_x ? s2.b.x : s2.b.y;
        const double lo = std::max(std::min(a0, a1), std::min(b0, b1));
        const double hi = std::min(std::max(a0, a1), std::max(b0, b1));
        if (lo > hi + 1e-7) {
            return false;
        }
        auto at_axis = [&](double value) {
            if (std::abs(a1 - a0) <= kDmeEps) {
                return s1.a;
            }
            return dme_lerp(s1.a, s1.b, (value - a0) / (a1 - a0));
        };
        out.a = at_axis(lo);
        out.b = at_axis(hi);
        out.valid = true;
        return true;
    }

    const double t = ((q.x - p.x) * s.y - (q.y - p.y) * s.x) / rxs;
    const double u = ((q.x - p.x) * r.y - (q.y - p.y) * r.x) / rxs;
    if (t < -1e-7 || t > 1.0 + 1e-7 || u < -1e-7 || u > 1.0 + 1e-7) {
        return false;
    }
    const DmePoint hit = dme_lerp(s1.a, s1.b, std::clamp(t, 0.0, 1.0));
    out = dme_make_point_segment(hit);
    return true;
}

double dme_segment_distance(const DmeSegment& a, const DmeSegment& b) {
    DmeSegment overlap;
    if (dme_segment_intersection(a, b, overlap)) {
        return 0.0;
    }
    double best = std::numeric_limits<double>::infinity();
    for (const DmePoint& p : {a.a, a.b}) {
        best = std::min(best, dme_point_segment_distance(p, b));
    }
    for (const DmePoint& p : {b.a, b.b}) {
        best = std::min(best, dme_point_segment_distance(p, a));
    }
    return best;
}

std::vector<DmePoint> dme_convex_hull(std::vector<DmePoint> points) {
    std::sort(points.begin(), points.end(), [](const DmePoint& a,
                                               const DmePoint& b) {
        if (a.x != b.x) return a.x < b.x;
        return a.y < b.y;
    });
    std::vector<DmePoint> unique_points;
    for (const DmePoint& p : points) {
        if (unique_points.empty() ||
            !dme_same_point(unique_points.back(), p)) {
            unique_points.push_back(p);
        }
    }
    if (unique_points.size() <= 2) {
        return unique_points;
    }

    std::vector<DmePoint> hull;
    for (const DmePoint& p : unique_points) {
        while (hull.size() >= 2 &&
               dme_cross(hull[hull.size() - 2], hull.back(), p) <= 1e-9) {
            hull.pop_back();
        }
        hull.push_back(p);
    }
    const std::size_t lower_size = hull.size();
    for (std::size_t i = unique_points.size() - 1; i-- > 0;) {
        const DmePoint& p = unique_points[i];
        while (hull.size() > lower_size &&
               dme_cross(hull[hull.size() - 2], hull.back(), p) <= 1e-9) {
            hull.pop_back();
        }
        hull.push_back(p);
    }
    hull.pop_back();
    return hull;
}

std::vector<DmePoint> dme_trr_polygon(const DmeSegment& core, double radius) {
    if (!core.valid) {
        return {};
    }
    if (radius <= kDmeEps) {
        if (dme_segment_is_point(core)) {
            return {core.a};
        }
        return {core.a, core.b};
    }
    std::vector<DmePoint> points;
    points.reserve(8);
    for (const DmePoint& p : {core.a, core.b}) {
        points.push_back(dme_point(p.x + radius, p.y));
        points.push_back(dme_point(p.x - radius, p.y));
        points.push_back(dme_point(p.x, p.y + radius));
        points.push_back(dme_point(p.x, p.y - radius));
    }
    auto hull = dme_convex_hull(std::move(points));
    if (hull.size() < 2) {
        return {core.a};
    }
    return hull;
}

DmeSegment dme_clip_segment_to_polygon(const DmeSegment& segment,
                                       std::vector<DmePoint> polygon) {
    if (!segment.valid) {
        return DmeSegment{};
    }
    if (polygon.size() < 3) {
        if (polygon.empty()) {
            return DmeSegment{};
        }
        if (polygon.size() == 1) {
            if (dme_point_segment_distance(polygon.front(), segment) <= 1e-6) {
                return dme_make_point_segment(polygon.front());
            }
            return DmeSegment{};
        }
        DmeSegment poly_segment = dme_make_segment(polygon[0], polygon[1]);
        DmeSegment result;
        if (dme_segment_intersection(segment, poly_segment, result)) {
            return result;
        }
        return DmeSegment{};
    }

    double area2 = 0.0;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const DmePoint& a = polygon[i];
        const DmePoint& b = polygon[(i + 1) % polygon.size()];
        area2 += a.x * b.y - a.y * b.x;
    }
    if (area2 < 0.0) {
        std::reverse(polygon.begin(), polygon.end());
    }

    double t0 = 0.0;
    double t1 = 1.0;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const DmePoint& a = polygon[i];
        const DmePoint& b = polygon[(i + 1) % polygon.size()];
        const double f0 = dme_cross(a, b, segment.a);
        const double f1 = dme_cross(a, b, segment.b);
        if (f0 < -1e-7 && f1 < -1e-7) {
            return DmeSegment{};
        }
        if (std::abs(f0 - f1) <= 1e-12) {
            continue;
        }
        const double t = f0 / (f0 - f1);
        if (f0 < -1e-7) {
            t0 = std::max(t0, t);
        } else if (f1 < -1e-7) {
            t1 = std::min(t1, t);
        }
    }
    if (t1 < t0 - 1e-7) {
        return DmeSegment{};
    }
    return dme_make_segment(dme_lerp(segment.a, segment.b, std::clamp(t0, 0.0, 1.0)),
                            dme_lerp(segment.a, segment.b, std::clamp(t1, 0.0, 1.0)));
}

std::vector<DmeSegment> dme_boundary_segments(const DmeSegment& core,
                                             double radius) {
    std::vector<DmeSegment> result;
    if (!core.valid) {
        return result;
    }
    if (radius <= kDmeEps) {
        result.push_back(core);
        return result;
    }
    const auto polygon = dme_trr_polygon(core, radius);
    if (polygon.size() < 3) {
        result.push_back(core);
        return result;
    }
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        result.push_back(
            dme_make_segment(polygon[i], polygon[(i + 1) % polygon.size()]));
    }
    return result;
}

DmeSegment dme_merge_segment(const DmeSegment& a, double radius_a,
                             const DmeSegment& b, double radius_b) {
    DmeSegment best;
    double best_length = -1.0;
    auto consider = [&](const DmeSegment& piece) {
        if (!piece.valid) {
            return;
        }
        const double length = dme_l1(piece.a, piece.b);
        if (best_length < 0.0 || length > best_length + 1e-7) {
            best_length = length;
            best = piece;
        }
    };

    const auto boundary_a = dme_boundary_segments(a, radius_a);
    const auto boundary_b = dme_boundary_segments(b, radius_b);

    for (const DmeSegment& edge : boundary_a) {
        if (radius_b <= kDmeEps) {
            DmeSegment hit;
            if (dme_segment_intersection(edge, b, hit)) {
                consider(hit);
            }
        } else {
            consider(dme_clip_segment_to_polygon(
                edge, dme_trr_polygon(b, radius_b)));
        }
    }
    for (const DmeSegment& edge : boundary_b) {
        if (radius_a <= kDmeEps) {
            DmeSegment hit;
            if (dme_segment_intersection(edge, a, hit)) {
                consider(hit);
            }
        } else {
            consider(dme_clip_segment_to_polygon(
                edge, dme_trr_polygon(a, radius_a)));
        }
    }

    if (best_length >= 0.0) {
        return best;
    }

    // Geometric fallback for numerically degenerate cases.
    const DmePoint pa = dme_closest_point_on_segment(b.a, a);
    const DmePoint pb = dme_closest_point_on_segment(a.a, b);
    return dme_make_point_segment(dme_midpoint(pa, pb));
}

std::vector<DmePoint> dme_circle_intersections(const DmePoint& center,
                                               double radius,
                                               const DmeSegment& segment) {
    std::vector<DmePoint> result;
    if (!segment.valid) {
        return result;
    }
    if (radius <= kDmeEps) {
        const DmePoint q = dme_closest_point_on_segment(center, segment);
        if (dme_l1(center, q) <= 1e-6) {
            result.push_back(q);
        }
        return result;
    }

    const std::array<DmePoint, 4> diamond{
        dme_point(center.x - radius, center.y),
        dme_point(center.x, center.y + radius),
        dme_point(center.x + radius, center.y),
        dme_point(center.x, center.y - radius)};
    auto add = [&](const DmePoint& p) {
        for (const DmePoint& existing : result) {
            if (dme_same_point(existing, p)) {
                return;
            }
        }
        result.push_back(p);
    };
    for (int i = 0; i < 4; ++i) {
        DmeSegment hit;
        if (dme_segment_intersection(
                segment, dme_make_segment(diamond[static_cast<std::size_t>(i)],
                                          diamond[static_cast<std::size_t>((i + 1) % 4)]),
                hit)) {
            add(hit.a);
            add(hit.b);
        }
    }
    if (std::abs(dme_l1(center, segment.a) - radius) <= 1e-6) {
        add(segment.a);
    }
    if (std::abs(dme_l1(center, segment.b) - radius) <= 1e-6) {
        add(segment.b);
    }
    return result;
}

DmePoint dme_choose_exact(const DmePoint& parent, const DmeSegment& core,
                          double radius) {
    const auto candidates = dme_circle_intersections(parent, radius, core);
    if (!candidates.empty()) {
        const DmePoint middle = dme_midpoint(core.a, core.b);
        DmePoint best = candidates.front();
        double best_distance = dme_l1(best, middle);
        for (const DmePoint& candidate : candidates) {
            const double distance = dme_l1(candidate, middle);
            if (distance < best_distance) {
                best_distance = distance;
                best = candidate;
            }
        }
        return best;
    }
    return dme_closest_point_on_segment(parent, core);
}

DmeTree build_nearest_neighbor_tree(const std::vector<Sink>& sinks) {
    if (sinks.empty()) {
        fail("DME requires at least one sink");
    }
    if (sinks.size() > 2000) {
        fail("DME nearest-neighbor topology is limited to 2000 sinks "
             "in this implementation; use htree for larger cases");
    }

    struct Cluster {
        int node = -1;
        DmePoint representative;
    };

    const std::size_t n = sinks.size();
    DmeTree tree;
    tree.nodes.assign(2 * n - 1, DmeNode{});
    std::vector<Cluster> active;
    active.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        DmeNode& leaf = tree.nodes[i];
        leaf.sink_index = static_cast<int>(i);
        leaf.pos = dme_point(sinks[i].x, sinks[i].y);
        leaf.core = dme_make_point_segment(leaf.pos);
        leaf.capacitance = 1.0;
        active.push_back(
            Cluster{static_cast<int>(i), dme_point(sinks[i].x, sinks[i].y)});
    }

    int next_id = static_cast<int>(n);
    while (active.size() > 1) {
        int best_i = 0;
        int best_j = 1;
        double best_distance = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < active.size(); ++i) {
            for (std::size_t j = i + 1; j < active.size(); ++j) {
                const double distance =
                    dme_l1(active[i].representative, active[j].representative);
                if (distance < best_distance - 1e-12 ||
                    (std::abs(distance - best_distance) <= 1e-12 &&
                     std::tie(active[i].node, active[j].node) <
                         std::tie(active[best_i].node, active[best_j].node))) {
                    best_distance = distance;
                    best_i = static_cast<int>(i);
                    best_j = static_cast<int>(j);
                }
            }
        }

        const int left = active[static_cast<std::size_t>(best_i)].node;
        const int right = active[static_cast<std::size_t>(best_j)].node;
        const DmePoint left_rep =
            active[static_cast<std::size_t>(best_i)].representative;
        const DmePoint right_rep =
            active[static_cast<std::size_t>(best_j)].representative;
        const int parent = next_id++;
        tree.nodes[static_cast<std::size_t>(parent)].left = left;
        tree.nodes[static_cast<std::size_t>(parent)].right = right;
        tree.nodes[static_cast<std::size_t>(left)].parent = parent;
        tree.nodes[static_cast<std::size_t>(right)].parent = parent;

        active.erase(active.begin() + best_j);
        active.erase(active.begin() + best_i);
        active.push_back(
            Cluster{parent, dme_midpoint(left_rep, right_rep)});
    }
    tree.root = active.front().node;
    if (tree.root != static_cast<int>(2 * n - 2)) {
        fail("internal error: unexpected DME root id");
    }
    return tree;
}

int dme_tree_depth(const DmeTree& tree) {
    if (tree.root < 0) {
        return 0;
    }
    std::vector<int> depth(tree.nodes.size(), 0);
    std::vector<int> stack{tree.root};
    int max_depth = 0;
    while (!stack.empty()) {
        const int current = stack.back();
        stack.pop_back();
        for (const int child :
             {tree.nodes[static_cast<std::size_t>(current)].left,
              tree.nodes[static_cast<std::size_t>(current)].right}) {
            if (child >= 0) {
                depth[static_cast<std::size_t>(child)] =
                    depth[static_cast<std::size_t>(current)] + 1;
                max_depth = std::max(
                    max_depth, depth[static_cast<std::size_t>(child)]);
                stack.push_back(child);
            }
        }
    }
    return max_depth;
}

Topology emit_dme_topology(const std::vector<Sink>& sinks,
                           const DmeTree& tree) {
    Topology topology;
    topology.nodes.reserve(tree.nodes.size());
    for (std::size_t i = 0; i < tree.nodes.size(); ++i) {
        const DmeNode& source = tree.nodes[i];
        Node node;
        node.id = static_cast<int>(i);
        node.x = source.pos.x;
        node.y = source.pos.y;
        if (static_cast<int>(i) == tree.root && source.sink_index < 0) {
            node.kind = "root";
            ++topology.internal_nodes;
        } else if (source.sink_index >= 0) {
            node.kind = "sink";
            const Sink& sink =
                sinks[static_cast<std::size_t>(source.sink_index)];
            node.sink_id = sink.id;
            node.cap_ff = sink.cap_ff;
            ++topology.sink_nodes;
        } else {
            node.kind = "internal";
            ++topology.internal_nodes;
        }
        topology.nodes.push_back(node);
    }
    for (std::size_t i = 0; i < tree.nodes.size(); ++i) {
        const DmeNode& source = tree.nodes[i];
        const int parent_id = static_cast<int>(i);
        if (source.left >= 0) {
            topology.edges.push_back(Edge{parent_id, source.left, "HV"});
        }
        if (source.right >= 0) {
            topology.edges.push_back(Edge{parent_id, source.right, "HV"});
        }
    }
    topology.root = tree.root;
    if (static_cast<int>(topology.edges.size()) !=
        static_cast<int>(topology.nodes.size()) - 1) {
        fail("internal error: DME output is not a tree");
    }
    return topology;
}

}  // namespace cts
