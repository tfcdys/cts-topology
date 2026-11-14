// HP-DME: path-aware heterogeneous parallel DME.
//
// Extends PT-DME with:
//   - cumulative path_range propagation (path_min/path_max),
//   - a path-aware ternary merge cost J,
//   - geometric outlier detection and matching,
//   - four-terminal shared-bar H-tree candidates,
//   - deterministic parallel candidate evaluation.

#include "hp_dme.hpp"

#include "cts_common.hpp"
#include "dme_common.hpp"
#include "third_party/delaunator.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

namespace cts {
namespace {

constexpr double kHpEps = 1e-9;
constexpr int kTernarySamples = 32;

struct Point {
    double x = 0.0;
    double y = 0.0;
};

struct Segment {
    Point a;
    Point b;
};

struct HpActiveNode {
    int id = -1;
    Point pos;
    double capacitance = 0.0;
    double path_min = 0.0;
    double path_max = 0.0;
    int tree_root = -1;
    int wait_rounds = 0;
};

struct LocalNode {
    Point pos;
    int parent = -1;
};

// A candidate merge carries its actual local topology: nodes[0] is the new
// root; other local nodes are junctions, each with a parent local index.
struct MergeCandidate {
    Point root;
    std::vector<LocalNode> local_nodes;
    std::vector<int> children;       // active-node indices
    std::vector<int> attachments;    // local_nodes index for each child
    std::vector<double> branch_lengths;
    double path_min = 0.0;
    double path_max = 0.0;
    double wire_length = 0.0;
    double skew = 0.0;
    double max_path = 0.0;
    double cost = std::numeric_limits<double>::infinity();
    int serial = 0;
    bool valid = false;
};

struct CostRefs {
    double h_ref = 1.0;
    double s_ref = 1.0;
    double w_ref = 1.0;
};

struct Triangle {
    std::array<int, 3> v{-1, -1, -1};
    MergeCandidate ternary;
    double cost = std::numeric_limits<double>::infinity();
};

double l1(const Point& a, const Point& b) {
    return std::abs(a.x - b.x) + std::abs(a.y - b.y);
}

bool same_point(const Point& a, const Point& b) {
    return l1(a, b) <= 1e-8;
}

Point l1_path_midpoint(const Point& a, const Point& b) {
    const double lh = std::abs(b.x - a.x);
    const double lv = std::abs(b.y - a.y);
    const double length = lh + lv;
    if (length <= kHpEps) {
        return a;
    }
    const double half = 0.5 * length;
    if (half <= lh) {
        return Point{a.x + (b.x > a.x ? 1.0 : (b.x < a.x ? -1.0 : 0.0)) * half,
                     a.y};
    }
    return Point{b.x,
                 a.y + (b.y > a.y ? 1.0 : (b.y < a.y ? -1.0 : 0.0)) *
                           (half - lh)};
}

void add_segment(std::vector<Segment>& segments, const Point& a,
                 const Point& b) {
    if (l1(a, b) <= kHpEps) {
        return;
    }
    for (const Segment& existing : segments) {
        if ((same_point(existing.a, a) && same_point(existing.b, b)) ||
            (same_point(existing.a, b) && same_point(existing.b, a))) {
            return;
        }
    }
    segments.push_back(Segment{a, b});
}

std::vector<Segment> ternary_candidate_segments(const Point& m,
                                                const Point& c) {
    std::vector<Segment> segments;
    if (same_point(m, c)) {
        return segments;
    }
    add_segment(segments, m, Point{c.x, m.y});
    add_segment(segments, Point{c.x, m.y}, c);
    add_segment(segments, m, Point{m.x, c.y});
    add_segment(segments, Point{m.x, c.y}, c);
    return segments;
}

Point point_on_segment(const Segment& s, double t) {
    return Point{s.a.x + (s.b.x - s.a.x) * t,
                 s.a.y + (s.b.y - s.a.y) * t};
}

std::vector<double> segment_parameters(const Segment& segment,
                                       const std::vector<Point>& points) {
    std::vector<double> values{0.0, 1.0};
    const double dx = segment.b.x - segment.a.x;
    const double dy = segment.b.y - segment.a.y;
    for (int i = 1; i < kTernarySamples; ++i) {
        values.push_back(static_cast<double>(i) / kTernarySamples);
    }
    for (const Point& p : points) {
        if (std::abs(dx) > kHpEps) {
            values.push_back((p.x - segment.a.x) / dx);
        }
        if (std::abs(dy) > kHpEps) {
            values.push_back((p.y - segment.a.y) / dy);
        }
    }
    for (double& t : values) {
        t = std::clamp(t, 0.0, 1.0);
    }
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end(),
                             [](double a, double b) {
                                 return std::abs(a - b) <= 1e-12;
                             }),
                 values.end());
    return values;
}

bool candidate_better(const MergeCandidate& lhs, const MergeCandidate& rhs) {
    if (!rhs.valid) return true;
    if (lhs.cost < rhs.cost - kHpEps) return true;
    if (rhs.cost < lhs.cost - kHpEps) return false;
    if (lhs.max_path < rhs.max_path - kHpEps) return true;
    if (rhs.max_path < lhs.max_path - kHpEps) return false;
    if (lhs.skew < rhs.skew - kHpEps) return true;
    if (rhs.skew < lhs.skew - kHpEps) return false;
    if (lhs.wire_length < rhs.wire_length - kHpEps) return true;
    if (rhs.wire_length < lhs.wire_length - kHpEps) return false;
    if (lhs.root.x < rhs.root.x - kHpEps) return true;
    if (rhs.root.x < lhs.root.x - kHpEps) return false;
    if (lhs.root.y < rhs.root.y - kHpEps) return true;
    if (rhs.root.y < lhs.root.y - kHpEps) return false;
    return lhs.serial < rhs.serial;
}

CostRefs make_refs(const std::vector<const HpActiveNode*>& children) {
    CostRefs refs;
    double max_path = 0.0;
    double min_path = std::numeric_limits<double>::infinity();
    double max_skew = 0.0;
    double min_x = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    for (const HpActiveNode* child : children) {
        max_path = std::max(max_path, child->path_max);
        min_path = std::min(min_path, child->path_min);
        max_skew = std::max(max_skew, child->path_max - child->path_min);
        min_x = std::min(min_x, child->pos.x);
        max_x = std::max(max_x, child->pos.x);
        min_y = std::min(min_y, child->pos.y);
        max_y = std::max(max_y, child->pos.y);
    }
    refs.h_ref = std::max(1.0, max_path);
    refs.s_ref = std::max(1.0, max_skew);
    refs.w_ref = std::max(1.0, (max_x - min_x) + (max_y - min_y));
    return refs;
}

double evaluate_merge_cost(const MergeCandidate& candidate,
                           const CostRefs& refs, const HpDmeOptions& options) {
    return options.alpha * (candidate.path_max / refs.h_ref) +
           options.beta * (candidate.skew / refs.s_ref) +
           options.gamma * (candidate.wire_length / refs.w_ref);
}

MergeCandidate make_candidate(const Point& root,
                              std::vector<LocalNode> local_nodes,
                              const std::vector<int>& children,
                              const std::vector<int>& attachments,
                              const std::vector<const HpActiveNode*>& child_ptrs,
                              const CostRefs& refs,
                              const HpDmeOptions& options, int serial) {
    MergeCandidate candidate;
    candidate.root = root;
    candidate.local_nodes = std::move(local_nodes);
    candidate.children = children;
    candidate.attachments = attachments;
    candidate.serial = serial;

    std::vector<double> local_distance(candidate.local_nodes.size(), 0.0);
    double wire = 0.0;
    for (std::size_t i = 1; i < candidate.local_nodes.size(); ++i) {
        const int parent = candidate.local_nodes[i].parent;
        if (parent < 0 || static_cast<std::size_t>(parent) >= i) {
            return candidate;
        }
        const double length =
            l1(candidate.local_nodes[i].pos,
               candidate.local_nodes[static_cast<std::size_t>(parent)].pos);
        local_distance[i] = local_distance[static_cast<std::size_t>(parent)] + length;
        wire += length;
    }

    candidate.branch_lengths.resize(children.size(), 0.0);
    double high = -std::numeric_limits<double>::infinity();
    double low = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < children.size(); ++i) {
        const int attachment = attachments[i];
        if (attachment < 0 ||
            static_cast<std::size_t>(attachment) >= candidate.local_nodes.size()) {
            return candidate;
        }
        const double branch =
            local_distance[static_cast<std::size_t>(attachment)] +
            l1(candidate.local_nodes[static_cast<std::size_t>(attachment)].pos,
               child_ptrs[i]->pos);
        candidate.branch_lengths[i] = branch;
        wire += l1(candidate.local_nodes[static_cast<std::size_t>(attachment)].pos,
                   child_ptrs[i]->pos);
        high = std::max(high, child_ptrs[i]->path_max + branch);
        low = std::min(low, child_ptrs[i]->path_min + branch);
    }

    candidate.path_max = high;
    candidate.path_min = low;
    candidate.skew = high - low;
    candidate.max_path = high;
    candidate.wire_length = wire;
    candidate.cost = evaluate_merge_cost(candidate, refs, options);
    candidate.valid = true;
    return candidate;
}

std::vector<const HpActiveNode*> ptrs(
    const std::vector<HpActiveNode>& active, const std::vector<int>& indices) {
    std::vector<const HpActiveNode*> result;
    result.reserve(indices.size());
    for (const int index : indices) {
        result.push_back(&active[static_cast<std::size_t>(index)]);
    }
    return result;
}

std::vector<int> ternary_pair_order(const std::vector<HpActiveNode>& active,
                                    const std::array<int, 3>& v) {
    struct PairChoice {
        double distance = std::numeric_limits<double>::infinity();
        std::array<int, 2> ids{0, 0};
        int a = 0;
        int b = 1;
        int c = 2;
    } best;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            const int k = 3 - i - j;
            const double distance =
                l1(active[static_cast<std::size_t>(v[static_cast<std::size_t>(i)])].pos,
                   active[static_cast<std::size_t>(v[static_cast<std::size_t>(j)])].pos);
            std::array<int, 2> ids{
                active[static_cast<std::size_t>(v[static_cast<std::size_t>(i)])].id,
                active[static_cast<std::size_t>(v[static_cast<std::size_t>(j)])].id};
            std::sort(ids.begin(), ids.end());
            if (distance < best.distance - kHpEps ||
                (std::abs(distance - best.distance) <= kHpEps &&
                 ids < best.ids)) {
                best = PairChoice{distance, ids, i, j, k};
            }
        }
    }
    return {v[static_cast<std::size_t>(best.a)],
            v[static_cast<std::size_t>(best.b)],
            v[static_cast<std::size_t>(best.c)]};
}

MergeCandidate optimize_ternary_merge(const std::vector<HpActiveNode>& active,
                                      const std::array<int, 3>& vertices,
                                      const HpDmeOptions& options,
                                      int serial_base = 0) {
    const std::vector<int> ordered = ternary_pair_order(active, vertices);
    const std::vector<int> children{vertices[0], vertices[1], vertices[2]};
    const auto child_ptrs = ptrs(active, children);
    const CostRefs refs = make_refs(child_ptrs);

    const Point pa = active[static_cast<std::size_t>(ordered[0])].pos;
    const Point pb = active[static_cast<std::size_t>(ordered[1])].pos;
    const Point pc = active[static_cast<std::size_t>(ordered[2])].pos;
    const Point m = l1_path_midpoint(pa, pb);
    const auto segments = ternary_candidate_segments(m, pc);

    MergeCandidate best;
    int serial = serial_base;
    auto consider = [&](const Point& root) {
        MergeCandidate candidate = make_candidate(
            root, {LocalNode{root, -1}}, children, {0, 0, 0}, child_ptrs,
            refs, options, serial++);
        if (candidate_better(candidate, best)) {
            best = std::move(candidate);
        }
    };

    if (segments.empty()) {
        consider(m);
        return best;
    }
    const std::vector<Point> points{pa, pb, pc, m};
    for (const Segment& segment : segments) {
        for (const double t : segment_parameters(segment, points)) {
            consider(point_on_segment(segment, t));
        }
    }
    return best;
}

double point_to_segment_distance(const Point& p, const Point& a,
                                 const Point& b) {
    const DmePoint dp = dme_point(p.x, p.y);
    const DmeSegment s = dme_make_segment(dme_point(a.x, a.y),
                                          dme_point(b.x, b.y));
    return dme_point_segment_distance(dp, s);
}

bool point_in_triangle(const Point& p, const Point& a, const Point& b,
                       const Point& c) {
    auto cross = [](const Point& u, const Point& v, const Point& w) {
        return (v.x - u.x) * (w.y - u.y) - (v.y - u.y) * (w.x - u.x);
    };
    const double c1 = cross(a, b, p);
    const double c2 = cross(b, c, p);
    const double c3 = cross(c, a, p);
    const double tolerance = 1e-8;
    return (c1 >= -tolerance && c2 >= -tolerance && c3 >= -tolerance) ||
           (c1 <= tolerance && c2 <= tolerance && c3 <= tolerance);
}

double point_to_triangle_distance(const Point& p, const Point& a,
                                  const Point& b, const Point& c) {
    if (point_in_triangle(p, a, b, c)) {
        return 0.0;
    }
    return std::min({point_to_segment_distance(p, a, b),
                     point_to_segment_distance(p, b, c),
                     point_to_segment_distance(p, c, a)});
}

std::vector<int> pairing_members(int pairing, int side) {
    // side 0 = first pair, side 1 = second pair.
    static const int members[3][2][2] = {
        {{0, 1}, {2, 3}},
        {{0, 2}, {1, 3}},
        {{0, 3}, {1, 2}},
    };
    return {members[pairing][side][0], members[pairing][side][1]};
}

MergeCandidate optimize_four_terminal_htree(
    const std::vector<HpActiveNode>& active, int i0, int i1, int i2, int ip,
    const HpDmeOptions& options, int serial_base = 0) {
    const std::array<int, 4> child_index{i0, i1, i2, ip};
    const std::vector<int> children{i0, i1, i2, ip};
    const auto child_ptrs = ptrs(active, children);
    const CostRefs refs = make_refs(child_ptrs);

    double min_x = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    std::vector<double> xs;
    std::vector<double> ys;
    for (const int index : child_index) {
        const Point p = active[static_cast<std::size_t>(index)].pos;
        xs.push_back(p.x);
        ys.push_back(p.y);
        min_x = std::min(min_x, p.x);
        max_x = std::max(max_x, p.x);
        min_y = std::min(min_y, p.y);
        max_y = std::max(max_y, p.y);
    }
    std::sort(xs.begin(), xs.end());
    std::sort(ys.begin(), ys.end());

    std::vector<double> x0_candidates{
        0.5 * (min_x + max_x),
        std::accumulate(xs.begin(), xs.end(), 0.0) / 4.0,
        0.5 * (xs[1] + xs[2]),
        xs[0], xs[1], xs[2], xs[3]};
    std::vector<double> y0_candidates{
        0.5 * (min_y + max_y),
        std::accumulate(ys.begin(), ys.end(), 0.0) / 4.0,
        0.5 * (ys[1] + ys[2])};
    auto normalize = [](std::vector<double>& values) {
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end(),
                                 [](double a, double b) {
                                     return std::abs(a - b) <= kHpEps;
                                 }),
                     values.end());
    };
    normalize(x0_candidates);
    normalize(y0_candidates);

    MergeCandidate best;
    int serial = serial_base;
    for (int pairing = 0; pairing < 3; ++pairing) {
        for (int orientation = 0; orientation < 2; ++orientation) {
            const int top_side = orientation == 0 ? 0 : 1;
            const int bottom_side = 1 - top_side;
            const std::vector<int> top_pos = pairing_members(pairing, top_side);
            const std::vector<int> bottom_pos =
                pairing_members(pairing, bottom_side);
            const Point top_a = active[static_cast<std::size_t>(
                child_index[static_cast<std::size_t>(top_pos[0])])].pos;
            const Point top_b = active[static_cast<std::size_t>(
                child_index[static_cast<std::size_t>(top_pos[1])])].pos;
            const Point bottom_a = active[static_cast<std::size_t>(
                child_index[static_cast<std::size_t>(bottom_pos[0])])].pos;
            const Point bottom_b = active[static_cast<std::size_t>(
                child_index[static_cast<std::size_t>(bottom_pos[1])])].pos;
            const double y_top = 0.5 * (top_a.y + top_b.y);
            const double y_bottom = 0.5 * (bottom_a.y + bottom_b.y);
            const double top_lo = std::min(top_a.x, top_b.x);
            const double top_hi = std::max(top_a.x, top_b.x);
            const double bottom_lo = std::min(bottom_a.x, bottom_b.x);
            const double bottom_hi = std::max(bottom_a.x, bottom_b.x);

            std::vector<double> y_roots = y0_candidates;
            y_roots.push_back(y_top);
            y_roots.push_back(y_bottom);
            y_roots.push_back(0.5 * (y_top + y_bottom));
            normalize(y_roots);

            for (const double x0 : x0_candidates) {
                const double x_top = std::clamp(x0, top_lo, top_hi);
                const double x_bottom = std::clamp(x0, bottom_lo, bottom_hi);
                const Point top_junction{x_top, y_top};
                const Point bottom_junction{x_bottom, y_bottom};
                for (const double y_root : y_roots) {
                    const Point root{x0, y_root};
                    std::vector<LocalNode> local_nodes;
                    local_nodes.push_back(LocalNode{root, -1});
                    auto attachment_for = [&](const Point& p) {
                        if (same_point(p, root)) {
                            return 0;
                        }
                        for (std::size_t i = 1; i < local_nodes.size(); ++i) {
                            if (same_point(local_nodes[i].pos, p)) {
                                return static_cast<int>(i);
                            }
                        }
                        local_nodes.push_back(
                            LocalNode{p, 0});
                        return static_cast<int>(local_nodes.size() - 1);
                    };
                    std::vector<int> attachments(4, 0);
                    for (int i = 0; i < 4; ++i) {
                        bool top = false;
                        for (const int p : top_pos) {
                            if (p == i) top = true;
                        }
                        attachments[static_cast<std::size_t>(i)] =
                            attachment_for(top ? top_junction : bottom_junction);
                    }
                    MergeCandidate candidate = make_candidate(
                        root, std::move(local_nodes), children, attachments,
                        child_ptrs, refs, options, serial++);
                    if (candidate_better(candidate, best)) {
                        best = std::move(candidate);
                    }
                }
            }
        }
    }
    return best;
}

MergeCandidate optimize_binary_merge(const std::vector<HpActiveNode>& active,
                                     int ia, int ib,
                                     const HpDmeOptions& options,
                                     int serial_base = 0) {
    const std::vector<int> children{ia, ib};
    const auto child_ptrs = ptrs(active, children);
    const CostRefs refs = make_refs(child_ptrs);
    const Point a = active[static_cast<std::size_t>(ia)].pos;
    const Point b = active[static_cast<std::size_t>(ib)].pos;

    std::vector<Point> path;
    path.push_back(a);
    path.push_back(Point{b.x, a.y});
    path.push_back(b);
    const double length = l1(a, b);
    MergeCandidate best;
    int serial = serial_base;
    auto consider = [&](const Point& root) {
        MergeCandidate candidate = make_candidate(
            root, {LocalNode{root, -1}}, children, {0, 0}, child_ptrs, refs,
            options, serial++);
        if (candidate_better(candidate, best)) {
            best = std::move(candidate);
        }
    };
    if (length <= kHpEps) {
        consider(a);
        return best;
    }
    const int samples = 64;
    for (int i = 0; i <= samples; ++i) {
        const double t = static_cast<double>(i) / samples;
        const double walked = t * length;
        const double lh = std::abs(b.x - a.x);
        Point root;
        if (walked <= lh) {
            root = Point{a.x + (b.x >= a.x ? 1.0 : -1.0) * walked, a.y};
        } else {
            root = Point{b.x,
                         a.y + (b.y >= a.y ? 1.0 : -1.0) * (walked - lh)};
        }
        consider(root);
    }
    return best;
}

std::vector<Triangle> build_delaunay_triangles(
    const std::vector<HpActiveNode>& active) {
    if (active.size() < 3) {
        return {};
    }
    std::vector<double> coordinates;
    coordinates.reserve(active.size() * 2);
    for (const HpActiveNode& node : active) {
        coordinates.push_back(node.pos.x);
        coordinates.push_back(node.pos.y);
    }
    std::vector<std::array<int, 3>> raw;
    try {
        delaunator::Delaunator delaunay(coordinates);
        raw.reserve(delaunay.triangles.size() / 3);
        for (std::size_t i = 0; i + 2 < delaunay.triangles.size(); i += 3) {
            const std::size_t i0 = delaunay.triangles[i];
            const std::size_t i1 = delaunay.triangles[i + 1];
            const std::size_t i2 = delaunay.triangles[i + 2];
            if (i0 >= active.size() || i1 >= active.size() ||
                i2 >= active.size()) {
                continue;
            }
            std::array<int, 3> indices{static_cast<int>(i0),
                                       static_cast<int>(i1),
                                       static_cast<int>(i2)};
            const Point a = active[i0].pos;
            const Point b = active[i1].pos;
            const Point c = active[i2].pos;
            const double area2 =
                (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
            if (std::abs(area2) <= 1e-10) {
                continue;
            }
            std::sort(indices.begin(), indices.end(), [&](int lhs, int rhs) {
                return active[static_cast<std::size_t>(lhs)].id <
                       active[static_cast<std::size_t>(rhs)].id;
            });
            raw.push_back(indices);
        }
    } catch (const std::exception&) {
        return {};
    }
    std::sort(raw.begin(), raw.end(), [&](const std::array<int, 3>& lhs,
                                          const std::array<int, 3>& rhs) {
        std::array<int, 3> lid{active[static_cast<std::size_t>(lhs[0])].id,
                               active[static_cast<std::size_t>(lhs[1])].id,
                               active[static_cast<std::size_t>(lhs[2])].id};
        std::array<int, 3> rid{active[static_cast<std::size_t>(rhs[0])].id,
                               active[static_cast<std::size_t>(rhs[1])].id,
                               active[static_cast<std::size_t>(rhs[2])].id};
        return lid < rid;
    });
    raw.erase(std::unique(raw.begin(), raw.end()), raw.end());
    std::vector<Triangle> triangles;
    triangles.reserve(raw.size());
    for (const auto& t : raw) {
        Triangle triangle;
        triangle.v = t;
        triangles.push_back(triangle);
    }
    return triangles;
}

bool all_collinear(const std::vector<HpActiveNode>& active) {
    if (active.size() <= 2) {
        return true;
    }
    const Point p0 = active[0].pos;
    std::size_t second = 1;
    while (second < active.size() && same_point(p0, active[second].pos)) {
        ++second;
    }
    if (second >= active.size()) {
        return true;
    }
    const Point p1 = active[second].pos;
    const double dx = p1.x - p0.x;
    const double dy = p1.y - p0.y;
    const double span = std::max(std::abs(dx), std::abs(dy));
    if (span <= kHpEps) {
        return true;
    }
    const double tolerance = 1e-8 * std::max(1.0, span);
    for (const HpActiveNode& node : active) {
        const double cross =
            dx * (node.pos.y - p0.y) - dy * (node.pos.x - p0.x);
        if (std::abs(cross) > tolerance) {
            return false;
        }
    }
    return true;
}

double median(std::vector<double> values) {
    if (values.empty()) {
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    const std::size_t middle = values.size() / 2;
    if (values.size() % 2 == 1) {
        return values[middle];
    }
    return 0.5 * (values[middle - 1] + values[middle]);
}

struct MatchCell {
    double distance = 0.0;
    MergeCandidate candidate;
};

struct MatchedOutlier {
    int outlier_index = -1;
    int selected_index = -1;
    MergeCandidate candidate;
    double match_cost = 0.0;
    int wait_rounds = 0;
};


std::string hp_metadata(const HpDmeResult& result,
                        const HpDmeOptions& options) {
    std::ostringstream out;
    out << "{\"variant\":\"HP-DME\","
        << "\"delaunay\":\"delaunator-cpp (MIT)\","
        << "\"delay_model\":\"path-aware-linear\","
        << "\"alpha\":" << format_number(options.alpha) << ","
        << "\"beta\":" << format_number(options.beta) << ","
        << "\"gamma\":" << format_number(options.gamma) << ","
        << "\"lambda\":" << format_number(options.lambda) << ","
        << "\"threads\":" << result.thread_count << ","
        << "\"runtime_ms\":" << format_number(result.runtime_ms) << ","
        << "\"ternary_merges\":" << result.ternary_merges << ","
        << "\"four_terminal_merges\":" << result.four_terminal_merges << ","
        << "\"outliers_detected\":" << result.outliers_detected << ","
        << "\"matched_outliers\":" << result.matched_outliers << ","
        << "\"match_rate\":"
        << format_number(result.outliers_detected == 0
                             ? 0.0
                             : static_cast<double>(result.matched_outliers) /
                                   static_cast<double>(result.outliers_detected))
        << ","
        << "\"fallback_binary\":" << result.fallback_binary << ","
        << "\"max_wait_rounds\":" << result.max_wait_rounds << ","
        << "\"average_matched_outlier_wait\":"
        << format_number(result.average_matched_outlier_wait) << ",\"rounds\":[";
    for (std::size_t i = 0; i < result.rounds.size(); ++i) {
        const HpDmeRoundStat& r = result.rounds[i];
        if (i) out << ',';
        out << "{\"round\":" << r.round
            << ",\"active_before\":" << r.active_before
            << ",\"duplicate_merges\":" << r.duplicate_merges
            << ",\"triangles\":" << r.triangles
            << ",\"selected_triangles\":" << r.selected_triangles
            << ",\"ternary_merges\":" << r.ternary_merges
            << ",\"four_terminal_merges\":" << r.four_terminal_merges
            << ",\"outliers\":" << r.outliers
            << ",\"matched_outliers\":" << r.matched_outliers
            << ",\"fallback_binary\":" << r.fallback_binary
            << ",\"carry_over\":" << r.carry_over
            << ",\"active_after\":" << r.active_after
            << ",\"average_outlier_wait\":"
            << format_number(r.average_outlier_wait)
            << ",\"max_wait_rounds\":" << r.max_wait_rounds << '}';
    }
    out << "]}";
    return out.str();
}

class HpDmeBuilder {
public:
    HpDmeBuilder(const std::vector<Sink>& sinks, HpDmeOptions options)
        : sinks_(sinks), options_(options) {}

    HpDmeResult run() {
        const auto start = std::chrono::steady_clock::now();
        initialize();
        if (!active_.empty()) {
            while (active_.size() > 2) {
                run_round();
            }
            finalize();
        }
        const auto end = std::chrono::steady_clock::now();
        result_.runtime_ms =
            std::chrono::duration<double, std::milli>(end - start).count();
        result_.max_wait_rounds = std::max(result_.max_wait_rounds, max_wait_);
        double wait_sum = 0.0;
        int wait_count = 0;
        for (const HpDmeRoundStat& round : result_.rounds) {
            if (round.matched_outliers > 0) {
                wait_sum += round.average_outlier_wait *
                            static_cast<double>(round.matched_outliers);
                wait_count += round.matched_outliers;
            }
        }
        result_.average_matched_outlier_wait =
            wait_count > 0 ? wait_sum / static_cast<double>(wait_count) : 0.0;
        result_.topology = std::move(topology_);
        result_.topology.metadata_json =
            hp_metadata(result_, options_);
        return std::move(result_);
    }

private:
    const std::vector<Sink>& sinks_;
    HpDmeOptions options_;
    Topology topology_;
    std::vector<HpActiveNode> active_;
    HpDmeResult result_;
    int next_active_id_ = 0;
    int next_tree_id_ = 0;
    int max_wait_ = 0;

    void initialize() {
        topology_.nodes.clear();
        topology_.edges.clear();
        active_.clear();
        result_.rounds.clear();
        result_.thread_count = 1;
        if (sinks_.empty()) {
            topology_.root = -1;
            return;
        }
        const std::size_t n = sinks_.size();
        active_.reserve(n);
        topology_.nodes.reserve(2 * n);
        for (std::size_t i = 0; i < n; ++i) {
            const Sink& sink = sinks_[i];
            Node node;
            node.id = next_tree_id_++;
            node.x = sink.x;
            node.y = sink.y;
            node.kind = "sink";
            node.sink_id = sink.id;
            node.cap_ff = sink.cap_ff;
            topology_.nodes.push_back(node);
            ++topology_.sink_nodes;
            active_.push_back(HpActiveNode{
                static_cast<int>(i), Point{sink.x, sink.y}, sink.cap_ff,
                0.0, 0.0, node.id, 0});
            ++next_active_id_;
        }
    }

    int make_tree_node(const Point& pos, const std::string& kind) {
        Node node;
        node.id = next_tree_id_++;
        node.x = pos.x;
        node.y = pos.y;
        node.kind = kind;
        topology_.nodes.push_back(node);
        if (kind == "internal") {
            ++topology_.internal_nodes;
        } else if (kind == "bend") {
            ++topology_.bend_nodes;
        }
        return node.id;
    }

    HpActiveNode commit_candidate(const MergeCandidate& candidate,
                                  const std::vector<int>& child_indices) {
        if (!candidate.valid) {
            fail("HP-DME internal error: committing invalid candidate");
        }
        std::vector<int> local_ids(candidate.local_nodes.size(), -1);
        for (std::size_t i = 0; i < candidate.local_nodes.size(); ++i) {
            local_ids[i] = make_tree_node(
                candidate.local_nodes[i].pos, i == 0 ? "internal" : "bend");
        }
        for (std::size_t i = 1; i < candidate.local_nodes.size(); ++i) {
            const int parent = candidate.local_nodes[i].parent;
            topology_.edges.push_back(
                Edge{local_ids[static_cast<std::size_t>(parent)],
                     local_ids[i], "HV"});
        }
        double capacitance = 0.0;
        for (std::size_t i = 0; i < child_indices.size(); ++i) {
            const HpActiveNode& child =
                active_[static_cast<std::size_t>(child_indices[i])];
            capacitance += child.capacitance;
            topology_.edges.push_back(
                Edge{local_ids[static_cast<std::size_t>(
                          candidate.attachments[i])],
                     child.tree_root, "HV"});
        }
        HpActiveNode node;
        node.id = next_active_id_++;
        node.pos = candidate.root;
        node.capacitance = capacitance;
        node.path_min = candidate.path_min;
        node.path_max = candidate.path_max;
        node.tree_root = local_ids.front();
        node.wait_rounds = 0;
        return node;
    }

    HpActiveNode commit_binary(const MergeCandidate& candidate, int ia, int ib) {
        return commit_candidate(candidate, {ia, ib});
    }

    int handle_duplicate_coordinates() {
        if (active_.size() < 2) {
            return 0;
        }
        std::vector<std::size_t> order(active_.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](std::size_t lhs,
                                                  std::size_t rhs) {
            const HpActiveNode& a = active_[lhs];
            const HpActiveNode& b = active_[rhs];
            if (std::abs(a.pos.x - b.pos.x) > kHpEps) {
                return a.pos.x < b.pos.x;
            }
            if (std::abs(a.pos.y - b.pos.y) > kHpEps) {
                return a.pos.y < b.pos.y;
            }
            return a.id < b.id;
        });
        bool duplicate = false;
        for (std::size_t i = 1; i < order.size(); ++i) {
            if (same_point(active_[order[i]].pos,
                           active_[order[i - 1]].pos)) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            return 0;
        }

        std::vector<HpActiveNode> next;
        next.reserve(active_.size());
        int merges = 0;
        for (std::size_t i = 0; i < order.size();) {
            std::size_t j = i + 1;
            while (j < order.size() &&
                   same_point(active_[order[i]].pos,
                              active_[order[j]].pos)) {
                ++j;
            }
            std::vector<std::size_t> group(order.begin() + i, order.begin() + j);
            std::sort(group.begin(), group.end(), [&](std::size_t lhs,
                                                       std::size_t rhs) {
                return active_[lhs].id < active_[rhs].id;
            });
            int current_index = static_cast<int>(group.front());
            for (std::size_t k = 1; k < group.size(); ++k) {
                const int right_index = static_cast<int>(group[k]);
                const MergeCandidate candidate =
                    optimize_binary_merge(active_, current_index, right_index,
                                          options_, 0);
                HpActiveNode merged =
                    commit_binary(candidate, current_index, right_index);
                active_.push_back(merged);
                current_index = static_cast<int>(active_.size()) - 1;
                ++merges;
            }
            next.push_back(active_[static_cast<std::size_t>(current_index)]);
            i = j;
        }
        active_ = std::move(next);
        return merges;
    }

    int merge_collinear_round() {
        double min_x = std::numeric_limits<double>::infinity();
        double max_x = -std::numeric_limits<double>::infinity();
        double min_y = std::numeric_limits<double>::infinity();
        double max_y = -std::numeric_limits<double>::infinity();
        for (const HpActiveNode& node : active_) {
            min_x = std::min(min_x, node.pos.x);
            max_x = std::max(max_x, node.pos.x);
            min_y = std::min(min_y, node.pos.y);
            max_y = std::max(max_y, node.pos.y);
        }
        const bool primary_x = (max_x - min_x) >= (max_y - min_y);
        std::vector<std::size_t> order(active_.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](std::size_t lhs,
                                                  std::size_t rhs) {
            const HpActiveNode& a = active_[lhs];
            const HpActiveNode& b = active_[rhs];
            const double a1 = primary_x ? a.pos.x : a.pos.y;
            const double b1 = primary_x ? b.pos.x : b.pos.y;
            const double a2 = primary_x ? a.pos.y : a.pos.x;
            const double b2 = primary_x ? b.pos.y : b.pos.x;
            if (std::abs(a1 - b1) > kHpEps) return a1 < b1;
            if (std::abs(a2 - b2) > kHpEps) return a2 < b2;
            return a.id < b.id;
        });

        std::vector<HpActiveNode> next;
        next.reserve(active_.size());
        int merges = 0;
        std::size_t i = 0;
        for (; i + 2 < order.size(); i += 3) {
            const int ia = static_cast<int>(order[i]);
            const int ib = static_cast<int>(order[i + 1]);
            const int ic = static_cast<int>(order[i + 2]);
            const std::array<int, 3> v{ia, ib, ic};
            const MergeCandidate candidate =
                optimize_ternary_merge(active_, v, options_, merges * 10000);
            next.push_back(commit_candidate(candidate, {ia, ib, ic}));
            ++merges;
        }
        for (; i < order.size(); ++i) {
            HpActiveNode node = active_[order[i]];
            ++node.wait_rounds;
            max_wait_ = std::max(max_wait_, node.wait_rounds);
            next.push_back(node);
        }
        active_ = std::move(next);
        return merges;
    }

    std::vector<int> uncovered_nodes(const std::vector<Triangle>& selected) const {
        std::vector<bool> used(active_.size(), false);
        for (const Triangle& triangle : selected) {
            for (const int v : triangle.v) {
                used[static_cast<std::size_t>(v)] = true;
            }
        }
        std::vector<int> uncovered;
        for (std::size_t i = 0; i < active_.size(); ++i) {
            if (!used[i]) uncovered.push_back(static_cast<int>(i));
        }
        return uncovered;
    }

    std::vector<int> detect_outliers(const std::vector<int>& uncovered) const {
        if (active_.size() < 2 || uncovered.empty()) {
            return {};
        }
        std::vector<double> nearest(active_.size(),
                                    std::numeric_limits<double>::infinity());
        for (std::size_t i = 0; i < active_.size(); ++i) {
            for (std::size_t j = i + 1; j < active_.size(); ++j) {
                const double distance = l1(active_[i].pos, active_[j].pos);
                nearest[i] = std::min(nearest[i], distance);
                nearest[j] = std::min(nearest[j], distance);
            }
        }
        const double median_distance = median(nearest);
        if (median_distance <= kHpEps) {
            return {};
        }
        std::vector<int> outliers;
        for (const int index : uncovered) {
            const double ratio = nearest[static_cast<std::size_t>(index)] /
                                 (median_distance + kHpEps);
            if (ratio > options_.outlier_ratio_threshold) {
                outliers.push_back(index);
            }
        }
        std::sort(outliers.begin(), outliers.end(), [&](int lhs, int rhs) {
            const double dl =
                nearest[static_cast<std::size_t>(lhs)] / median_distance;
            const double dr =
                nearest[static_cast<std::size_t>(rhs)] / median_distance;
            if (std::abs(dl - dr) > kHpEps) return dl > dr;
            if (active_[static_cast<std::size_t>(lhs)].wait_rounds !=
                active_[static_cast<std::size_t>(rhs)].wait_rounds) {
                return active_[static_cast<std::size_t>(lhs)].wait_rounds >
                       active_[static_cast<std::size_t>(rhs)].wait_rounds;
            }
            return active_[static_cast<std::size_t>(lhs)].id <
                   active_[static_cast<std::size_t>(rhs)].id;
        });
        return outliers;
    }

    std::array<int, 3> sorted_vertex_ids(const Triangle& triangle) const {
        std::array<int, 3> ids{
            active_[static_cast<std::size_t>(triangle.v[0])].id,
            active_[static_cast<std::size_t>(triangle.v[1])].id,
            active_[static_cast<std::size_t>(triangle.v[2])].id};
        std::sort(ids.begin(), ids.end());
        return ids;
    }

    void sort_triangles(std::vector<Triangle>& triangles) const {
        std::sort(triangles.begin(), triangles.end(),
                  [&](const Triangle& lhs, const Triangle& rhs) {
                      if (std::abs(lhs.cost - rhs.cost) > kHpEps) {
                          return lhs.cost < rhs.cost;
                      }
                      const int lw = std::max(
                          {active_[static_cast<std::size_t>(lhs.v[0])].wait_rounds,
                           active_[static_cast<std::size_t>(lhs.v[1])].wait_rounds,
                           active_[static_cast<std::size_t>(lhs.v[2])].wait_rounds});
                      const int rw = std::max(
                          {active_[static_cast<std::size_t>(rhs.v[0])].wait_rounds,
                           active_[static_cast<std::size_t>(rhs.v[1])].wait_rounds,
                           active_[static_cast<std::size_t>(rhs.v[2])].wait_rounds});
                      if (lw != rw) return lw > rw;
                      return sorted_vertex_ids(lhs) < sorted_vertex_ids(rhs);
                  });
    }

    std::vector<Triangle> greedy_mis(const std::vector<Triangle>& triangles) const {
        std::vector<bool> used(active_.size(), false);
        std::vector<Triangle> selected;
        for (const Triangle& triangle : triangles) {
            if (!used[static_cast<std::size_t>(triangle.v[0])] &&
                !used[static_cast<std::size_t>(triangle.v[1])] &&
                !used[static_cast<std::size_t>(triangle.v[2])]) {
                selected.push_back(triangle);
                used[static_cast<std::size_t>(triangle.v[0])] = true;
                used[static_cast<std::size_t>(triangle.v[1])] = true;
                used[static_cast<std::size_t>(triangle.v[2])] = true;
            }
        }
        for (std::size_t i = 0; i < selected.size(); ++i) {
            for (std::size_t j = i + 1; j < selected.size(); ++j) {
                for (const int a : selected[i].v) {
                    for (const int b : selected[j].v) {
                        if (a == b) {
                            fail("HP-DME internal error: selected triangles share a vertex");
                        }
                    }
                }
            }
        }
        for (const Triangle& triangle : triangles) {
            bool disjoint = true;
            for (const int v : triangle.v) {
                if (used[static_cast<std::size_t>(v)]) {
                    disjoint = false;
                    break;
                }
            }
            if (disjoint) {
                fail("HP-DME internal error: independent set is not maximal");
            }
        }
        return selected;
    }

    void compute_ternary_candidates(std::vector<Triangle>& triangles) {
        if (triangles.empty()) return;
        unsigned threads = options_.threads;
        if (threads == 0) threads = std::thread::hardware_concurrency();
        if (threads == 0) threads = 1;
        threads = std::min<unsigned>(threads,
                                     static_cast<unsigned>(triangles.size()));
        result_.thread_count = std::max(result_.thread_count,
                                        static_cast<int>(threads));
        std::atomic<std::size_t> next{0};
        auto worker = [&]() {
            while (true) {
                const std::size_t i = next.fetch_add(1);
                if (i >= triangles.size()) break;
                Triangle& triangle = triangles[i];
                triangle.ternary = optimize_ternary_merge(
                    active_, triangle.v, options_, 0);
                triangle.cost = triangle.ternary.cost;
            }
        };
        std::vector<std::thread> workers;
        workers.reserve(threads);
        for (unsigned i = 0; i < threads; ++i) workers.emplace_back(worker);
        for (std::thread& thread : workers) thread.join();
    }

    std::vector<MatchCell> compute_match_cells(
        const std::vector<int>& outliers, const std::vector<Triangle>& selected,
        int& threads_used) {
        const std::size_t o_count = outliers.size();
        const std::size_t s_count = selected.size();
        std::vector<MatchCell> cells(o_count * s_count);
        threads_used = 1;
        if (cells.empty()) return cells;
        unsigned threads = options_.threads;
        if (threads == 0) threads = std::thread::hardware_concurrency();
        if (threads == 0) threads = 1;
        threads = std::min<unsigned>(threads,
                                     static_cast<unsigned>(cells.size()));
        threads_used = static_cast<int>(threads);
        std::atomic<std::size_t> next{0};
        auto worker = [&]() {
            while (true) {
                const std::size_t i = next.fetch_add(1);
                if (i >= cells.size()) break;
                const std::size_t oi = i / s_count;
                const std::size_t si = i % s_count;
                const int outlier = outliers[oi];
                const Triangle& triangle = selected[si];
                const Point p = active_[static_cast<std::size_t>(outlier)].pos;
                const Point a = active_[static_cast<std::size_t>(triangle.v[0])].pos;
                const Point b = active_[static_cast<std::size_t>(triangle.v[1])].pos;
                const Point c = active_[static_cast<std::size_t>(triangle.v[2])].pos;
                MatchCell& cell = cells[i];
                cell.distance = point_to_triangle_distance(p, a, b, c);
                cell.candidate = optimize_four_terminal_htree(
                    active_, triangle.v[0], triangle.v[1], triangle.v[2],
                    outlier, options_, 0);
            }
        };
        std::vector<std::thread> workers;
        workers.reserve(threads);
        for (unsigned i = 0; i < threads; ++i) workers.emplace_back(worker);
        for (std::thread& thread : workers) thread.join();
        return cells;
    }

    void binary_fallback(std::vector<HpActiveNode>& next, int& used_count,
                         int& fallback_count, int serial_base) {
        int best_i = 0;
        int best_j = 1;
        double best_distance = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < active_.size(); ++i) {
            for (std::size_t j = i + 1; j < active_.size(); ++j) {
                const double distance = l1(active_[i].pos, active_[j].pos);
                if (distance < best_distance - kHpEps ||
                    (std::abs(distance - best_distance) <= kHpEps &&
                     std::tie(active_[i].id, active_[j].id) <
                         std::tie(active_[best_i].id, active_[best_j].id))) {
                    best_distance = distance;
                    best_i = static_cast<int>(i);
                    best_j = static_cast<int>(j);
                }
            }
        }
        const MergeCandidate candidate = optimize_binary_merge(
            active_, best_i, best_j, options_, serial_base);
        next.push_back(commit_binary(candidate, best_i, best_j));
        used_count = 2;
        ++fallback_count;
        for (std::size_t i = 0; i < active_.size(); ++i) {
            if (static_cast<int>(i) == best_i ||
                static_cast<int>(i) == best_j) {
                continue;
            }
            HpActiveNode node = active_[i];
            ++node.wait_rounds;
            max_wait_ = std::max(max_wait_, node.wait_rounds);
            next.push_back(node);
        }
    }

    void finalize() {
        if (active_.empty()) {
            fail("HP-DME produced an empty active set");
        }
        if (active_.size() == 2) {
            const MergeCandidate candidate =
                optimize_binary_merge(active_, 0, 1, options_, 900000000);
            active_ = {commit_binary(candidate, 0, 1)};
        }
        if (active_.size() != 1) {
            fail("HP-DME did not converge to one root");
        }
        const int root_id = active_[0].tree_root;
        if (root_id < 0 ||
            static_cast<std::size_t>(root_id) >= topology_.nodes.size()) {
            fail("HP-DME internal error: invalid final root");
        }
        if (topology_.nodes[static_cast<std::size_t>(root_id)].kind == "internal") {
            topology_.nodes[static_cast<std::size_t>(root_id)].kind = "root";
        }
        topology_.root = root_id;
        if (static_cast<int>(topology_.edges.size()) !=
            static_cast<int>(topology_.nodes.size()) - 1) {
            fail("HP-DME internal error: output is not a tree");
        }
        validate_final_paths();
    }

    void validate_final_paths() {
        const std::size_t n = topology_.nodes.size();
        std::vector<std::vector<std::pair<int, double>>> adjacency(n);
        for (const Edge& edge : topology_.edges) {
            const Node& a = topology_.nodes[static_cast<std::size_t>(edge.source)];
            const Node& b = topology_.nodes[static_cast<std::size_t>(edge.target)];
            adjacency[static_cast<std::size_t>(edge.source)].push_back(
                {edge.target, std::abs(a.x - b.x) + std::abs(a.y - b.y)});
        }
        double max_path = -1.0;
        double min_path = std::numeric_limits<double>::infinity();
        std::vector<std::pair<int, double>> stack{{topology_.root, 0.0}};
        while (!stack.empty()) {
            const auto [node, distance] = stack.back();
            stack.pop_back();
            if (topology_.nodes[static_cast<std::size_t>(node)].kind == "sink") {
                max_path = std::max(max_path, distance);
                min_path = std::min(min_path, distance);
            }
            for (const auto& [child, length] :
                 adjacency[static_cast<std::size_t>(node)]) {
                stack.push_back({child, distance + length});
            }
        }
        if (max_path < 0.0) {
            fail("HP-DME internal error: no sink reached in final validation");
        }
        if (std::abs(max_path - active_[0].path_max) > 1e-6 ||
            std::abs(min_path - active_[0].path_min) > 1e-6) {
            fail("HP-DME internal error: path ranges disagree with final DFS");
        }
    }

    void run_round() {
        HpDmeRoundStat stat;
        stat.round = static_cast<int>(result_.rounds.size());
        stat.active_before = static_cast<int>(active_.size());

        stat.duplicate_merges = handle_duplicate_coordinates();
        if (static_cast<int>(active_.size()) !=
            stat.active_before - stat.duplicate_merges) {
            fail("HP-DME internal error: duplicate count invariant");
        }
        const int base_active = static_cast<int>(active_.size());
        if (active_.size() <= 2) {
            stat.active_after = static_cast<int>(active_.size());
            stat.carry_over = stat.active_after;
            result_.rounds.push_back(stat);
            return;
        }

        if (all_collinear(active_)) {
            const int before = static_cast<int>(active_.size());
            stat.ternary_merges = merge_collinear_round();
            stat.active_after = static_cast<int>(active_.size());
            stat.carry_over = stat.active_after - stat.ternary_merges;
            stat.max_wait_rounds = max_wait_;
            if (stat.active_after != before - 2 * stat.ternary_merges) {
                fail("HP-DME internal error: collinear invariant");
            }
            result_.rounds.push_back(stat);
            return;
        }

        std::vector<Triangle> triangles = build_delaunay_triangles(active_);
        if (triangles.empty()) {
            std::vector<HpActiveNode> next;
            int used = 0;
            int fallback = 0;
            binary_fallback(next, used, fallback, 1000000);
            stat.fallback_binary = fallback;
            stat.active_after = static_cast<int>(next.size());
            stat.carry_over = static_cast<int>(next.size()) - 1;
            active_ = std::move(next);
            result_.fallback_binary += fallback;
            if (stat.active_after != base_active - 1) {
                fail("HP-DME internal error: binary fallback invariant");
            }
            result_.rounds.push_back(stat);
            return;
        }

        compute_ternary_candidates(triangles);
        stat.triangles = static_cast<int>(triangles.size());
        sort_triangles(triangles);
        const std::vector<Triangle> selected = greedy_mis(triangles);
        stat.selected_triangles = static_cast<int>(selected.size());
        if (selected.empty()) {
            std::vector<HpActiveNode> next;
            int used = 0;
            int fallback = 0;
            binary_fallback(next, used, fallback, 2000000);
            stat.fallback_binary = fallback;
            stat.active_after = static_cast<int>(next.size());
            stat.carry_over = static_cast<int>(next.size()) - 1;
            active_ = std::move(next);
            result_.fallback_binary += fallback;
            if (stat.active_after != base_active - 1) {
                fail("HP-DME internal error: binary fallback invariant");
            }
            result_.rounds.push_back(stat);
            return;
        }

        const std::vector<int> uncovered = uncovered_nodes(selected);
        const std::vector<int> outliers = detect_outliers(uncovered);
        stat.outliers = static_cast<int>(outliers.size());
        result_.outliers_detected += stat.outliers;

        std::vector<MatchedOutlier> matches;
        std::vector<std::vector<MatchCell>> cell_rows;
        std::vector<bool> matched_outlier(outliers.size(), false);
        std::vector<bool> triangle_matched(selected.size(), false);
        if (!outliers.empty()) {
            int match_threads = 1;
            const std::vector<MatchCell> cells =
                compute_match_cells(outliers, selected, match_threads);
            result_.thread_count = std::max(result_.thread_count, match_threads);
            const std::size_t s_count = selected.size();
            std::vector<double> distances;
            distances.reserve(cells.size());
            for (const MatchCell& cell : cells) distances.push_back(cell.distance);
            const double d_ref = std::max(1e-9, median(distances));
            for (std::size_t oi = 0; oi < outliers.size(); ++oi) {
                const int outlier = outliers[oi];
                int best_triangle = -1;
                double best_match_cost = std::numeric_limits<double>::infinity();
                std::array<int, 3> best_ids{0, 0, 0};
                for (std::size_t si = 0; si < selected.size(); ++si) {
                    if (triangle_matched[si]) continue;
                    const MatchCell& cell = cells[oi * s_count + si];
                    if (!cell.candidate.valid) continue;
                    const double match_cost =
                        cell.candidate.cost +
                        options_.lambda * cell.distance / d_ref;
                    const auto ids = sorted_vertex_ids(selected[si]);
                    if (match_cost < best_match_cost - kHpEps ||
                        (std::abs(match_cost - best_match_cost) <= kHpEps &&
                         (best_triangle < 0 || ids < best_ids))) {
                        best_match_cost = match_cost;
                        best_triangle = static_cast<int>(si);
                        best_ids = ids;
                    }
                }
                if (best_triangle < 0) continue;
                triangle_matched[static_cast<std::size_t>(best_triangle)] = true;
                matched_outlier[oi] = true;
                MatchedOutlier matched;
                matched.outlier_index = outlier;
                matched.selected_index = best_triangle;
                matched.candidate =
                    cells[oi * s_count + static_cast<std::size_t>(best_triangle)]
                        .candidate;
                matched.match_cost = best_match_cost;
                matched.wait_rounds = active_[static_cast<std::size_t>(outlier)]
                                          .wait_rounds;
                matches.push_back(std::move(matched));
            }
        }
        stat.matched_outliers = static_cast<int>(matches.size());
        result_.matched_outliers += stat.matched_outliers;

        double outlier_wait_sum = 0.0;
        for (const MatchedOutlier& match : matches) {
            outlier_wait_sum += match.wait_rounds;
        }
        stat.average_outlier_wait =
            matches.empty() ? 0.0
                            : outlier_wait_sum /
                                  static_cast<double>(matches.size());

        std::vector<bool> used(active_.size(), false);
        for (std::size_t si = 0; si < selected.size(); ++si) {
            for (const int v : selected[si].v) {
                used[static_cast<std::size_t>(v)] = true;
            }
        }
        for (const MatchedOutlier& match : matches) {
            used[static_cast<std::size_t>(match.outlier_index)] = true;
        }

        std::vector<HpActiveNode> next;
        next.reserve(active_.size());
        for (std::size_t si = 0; si < selected.size(); ++si) {
            if (triangle_matched[si]) {
                continue;
            }
            const Triangle& triangle = selected[si];
            next.push_back(commit_candidate(
                triangle.ternary, {triangle.v[0], triangle.v[1], triangle.v[2]}));
            ++stat.ternary_merges;
            ++result_.ternary_merges;
        }
        for (const MatchedOutlier& match : matches) {
            const Triangle& triangle = selected[match.selected_index];
            next.push_back(commit_candidate(
                match.candidate,
                {triangle.v[0], triangle.v[1], triangle.v[2],
                 match.outlier_index}));
            ++stat.four_terminal_merges;
            ++result_.four_terminal_merges;
        }

        int fallback_count = 0;
        // Long-wait unmatched outlier fallback: only unused nodes participate.
        int chosen_outlier = -1;
        std::vector<int> unused;
        for (std::size_t i = 0; i < active_.size(); ++i) {
            if (!used[i]) unused.push_back(static_cast<int>(i));
        }
        for (std::size_t oi = 0; oi < outliers.size(); ++oi) {
            if (matched_outlier[oi]) continue;
            const int index = outliers[oi];
            if (active_[static_cast<std::size_t>(index)].wait_rounds >=
                options_.max_wait_rounds) {
                chosen_outlier = index;
                break;
            }
        }
        if (chosen_outlier >= 0 && unused.size() >= 2) {
            int best_other = -1;
            double best_distance = std::numeric_limits<double>::infinity();
            for (const int candidate : unused) {
                if (candidate == chosen_outlier) continue;
                const double distance =
                    l1(active_[static_cast<std::size_t>(chosen_outlier)].pos,
                       active_[static_cast<std::size_t>(candidate)].pos);
                if (distance < best_distance - kHpEps ||
                    (std::abs(distance - best_distance) <= kHpEps &&
                     active_[static_cast<std::size_t>(candidate)].id <
                         (best_other < 0
                              ? std::numeric_limits<int>::max()
                              : active_[static_cast<std::size_t>(best_other)].id))) {
                    best_distance = distance;
                    best_other = candidate;
                }
            }
            if (best_other >= 0) {
                const MergeCandidate candidate = optimize_binary_merge(
                    active_, chosen_outlier, best_other, options_, 3000000);
                next.push_back(commit_binary(candidate, chosen_outlier,
                                             best_other));
                used[static_cast<std::size_t>(chosen_outlier)] = true;
                used[static_cast<std::size_t>(best_other)] = true;
                ++fallback_count;
            }
        }

        for (std::size_t i = 0; i < active_.size(); ++i) {
            if (used[i]) continue;
            HpActiveNode node = active_[i];
            ++node.wait_rounds;
            max_wait_ = std::max(max_wait_, node.wait_rounds);
            next.push_back(node);
            ++stat.carry_over;
        }

        const int expected = base_active - 2 * stat.ternary_merges -
                             3 * stat.four_terminal_merges - fallback_count;
        if (static_cast<int>(next.size()) != expected) {
            fail("HP-DME internal error: active-count invariant");
        }
        stat.fallback_binary = fallback_count;
        stat.active_after = static_cast<int>(next.size());
        stat.max_wait_rounds = max_wait_;
        result_.fallback_binary += fallback_count;
        result_.max_wait_rounds = std::max(result_.max_wait_rounds, max_wait_);
        active_ = std::move(next);
        result_.rounds.push_back(stat);
    }
};

}  // namespace

HpDmeResult build_hp_dme(const std::vector<Sink>& sinks,
                         const HpDmeOptions& options) {
    HpDmeBuilder builder(sinks, options);
    return builder.run();
}

}  // namespace cts

#ifndef HP_DME_NO_MAIN
int main(int argc, char** argv) {
    return cts::run_topology_cli(
        argc, argv, "hp_dme",
        [](const std::vector<cts::Sink>& sinks,
           const cts::BuildOptions&) {
            return cts::build_hp_dme(sinks).topology;
        },
        false, false, false);
}
#endif
