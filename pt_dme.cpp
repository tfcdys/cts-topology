// PT-DME: Parallel Triangular Deferred-Merge Embedding.
//
// This is the implementation of docs/pt-dme.md.  It uses Delaunator-C++ for
// the Delaunay triangulation and a deterministic greedy maximal independent
// set of vertex-disjoint triangles.  Ternary roots are evaluated in parallel;
// topology IDs are assigned deterministically by the main thread afterwards.

#include "pt_dme.hpp"

#include "cts_common.hpp"
#include "dme_common.hpp"
#include "third_party/delaunator.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

namespace cts {
namespace {

constexpr double kPtEps = 1e-9;
constexpr int kTernarySearchIterations = 80;

struct Point {
    double x = 0.0;
    double y = 0.0;
};

struct Segment {
    Point a;
    Point b;
};

struct ActiveNode {
    int id = -1;
    Point pos;
    double capacitance = 0.0;
    double delay_min = 0.0;
    double delay_max = 0.0;
    int tree_root = -1;
    int wait_rounds = 0;
};

struct Triangle {
    std::array<int, 3> v{-1, -1, -1};  // indices into active_nodes
    Point root;
    double cost = 0.0;
};

struct TernarySolution {
    Point root;
    double cost = 0.0;
};

struct LocalMergeResult {
    Point root;
    double capacitance = 0.0;
    double delay_min = 0.0;
    double delay_max = 0.0;
};

double l1(const Point& a, const Point& b) {
    return std::abs(a.x - b.x) + std::abs(a.y - b.y);
}

bool same_point(const Point& a, const Point& b) {
    return l1(a, b) <= 1e-8;
}

bool point_equal(const Point& a, const Point& b) {
    return std::abs(a.x - b.x) <= kPtEps &&
           std::abs(a.y - b.y) <= kPtEps;
}

Point l1_path_midpoint(const Point& a, const Point& b) {
    const double lh = std::abs(b.x - a.x);
    const double lv = std::abs(b.y - a.y);
    const double length = lh + lv;
    if (length <= kPtEps) {
        return a;
    }
    const double half = 0.5 * length;
    if (half <= lh) {
        return Point{
            a.x + (b.x > a.x ? 1.0 : (b.x < a.x ? -1.0 : 0.0)) * half,
            a.y};
    }
    return Point{
        b.x,
        a.y + (b.y > a.y ? 1.0 : (b.y < a.y ? -1.0 : 0.0)) *
                  (half - lh)};
}

void add_segment(std::vector<Segment>& segments, const Point& a,
                 const Point& b) {
    if (l1(a, b) <= kPtEps) {
        return;
    }
    for (const Segment& existing : segments) {
        if ((point_equal(existing.a, a) && point_equal(existing.b, b)) ||
            (point_equal(existing.a, b) && point_equal(existing.b, a))) {
            return;
        }
    }
    segments.push_back(Segment{a, b});
}

std::vector<Segment> candidate_segments(const Point& m, const Point& c) {
    std::vector<Segment> segments;
    if (same_point(m, c)) {
        return segments;
    }
    // HV from M to C.
    add_segment(segments, m, Point{c.x, m.y});
    add_segment(segments, Point{c.x, m.y}, c);
    // VH from M to C.
    add_segment(segments, m, Point{m.x, c.y});
    add_segment(segments, Point{m.x, c.y}, c);
    return segments;
}

double objective(const Point& r, const Point& a, const Point& b,
                 const Point& c) {
    return std::max({l1(r, a), l1(r, b), l1(r, c)});
}

Point point_on_segment(const Segment& s, double t) {
    return Point{s.a.x + (s.b.x - s.a.x) * t,
                 s.a.y + (s.b.y - s.a.y) * t};
}

Point best_on_segment(const Segment& segment, const Point& a, const Point& b,
                      const Point& c, const Point& m) {
    auto better = [&](const Point& lhs, double lhs_cost, const Point& rhs,
                      double rhs_cost) {
        if (lhs_cost < rhs_cost - kPtEps) return true;
        if (rhs_cost < lhs_cost - kPtEps) return false;
        const double dl = l1(lhs, m);
        const double dr = l1(rhs, m);
        if (dl < dr - kPtEps) return true;
        if (dr < dl - kPtEps) return false;
        if (std::abs(lhs.x - rhs.x) > kPtEps) return lhs.x < rhs.x;
        return lhs.y < rhs.y;
    };

    double lo = 0.0;
    double hi = 1.0;
    for (int iteration = 0; iteration < kTernarySearchIterations;
         ++iteration) {
        const double m1 = lo + (hi - lo) / 3.0;
        const double m2 = hi - (hi - lo) / 3.0;
        const double f1 = objective(point_on_segment(segment, m1), a, b, c);
        const double f2 = objective(point_on_segment(segment, m2), a, b, c);
        if (f1 <= f2) {
            hi = m2;
        } else {
            lo = m1;
        }
    }

    Point best = point_on_segment(segment, lo);
    double best_cost = objective(best, a, b, c);
    const std::array<double, 4> candidates{
        0.0, 1.0, lo, hi};
    for (const double t : candidates) {
        const Point p = point_on_segment(segment, t);
        const double value = objective(p, a, b, c);
        if (better(p, value, best, best_cost)) {
            best = p;
            best_cost = value;
        }
    }
    return best;
}

std::array<int, 3> sorted_ids(const std::array<int, 3>& v,
                              const std::vector<ActiveNode>& active) {
    std::array<int, 3> ids{active[static_cast<std::size_t>(v[0])].id,
                           active[static_cast<std::size_t>(v[1])].id,
                           active[static_cast<std::size_t>(v[2])].id};
    std::sort(ids.begin(), ids.end());
    return ids;
}

// Implements section 5 of docs/pt-dme.md.
TernarySolution compute_ternary_root(const ActiveNode& a,
                                     const ActiveNode& b,
                                     const ActiveNode& c) {
    struct PairChoice {
        double distance;
        std::array<int, 2> ids;
        int a;
        int b;
        int c;
    };

    const ActiveNode* nodes[3] = {&a, &b, &c};
    PairChoice best{std::numeric_limits<double>::infinity(), {0, 0}, 0, 0, 0};
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            const int k = 3 - i - j;
            const double distance = l1(nodes[i]->pos, nodes[j]->pos);
            std::array<int, 2> ids{nodes[i]->id, nodes[j]->id};
            std::sort(ids.begin(), ids.end());
            const bool better =
                distance < best.distance - kPtEps ||
                (std::abs(distance - best.distance) <= kPtEps &&
                 ids < best.ids);
            if (better) {
                best = PairChoice{distance, ids, i, j, k};
            }
        }
    }

    const Point pa = nodes[best.a]->pos;
    const Point pb = nodes[best.b]->pos;
    const Point pc = nodes[best.c]->pos;
    const Point m = l1_path_midpoint(pa, pb);

    TernarySolution solution;
    const auto segments = candidate_segments(m, pc);
    if (segments.empty()) {
        solution.root = m;
        solution.cost = objective(m, pa, pb, pc);
        return solution;
    }

    bool have_best = false;
    Point best_root;
    double best_cost = std::numeric_limits<double>::infinity();
    auto consider = [&](const Point& p, double cost) {
        if (!have_best || cost < best_cost - kPtEps ||
            (std::abs(cost - best_cost) <= kPtEps &&
             (l1(p, m) < l1(best_root, m) - kPtEps ||
              (std::abs(l1(p, m) - l1(best_root, m)) <= kPtEps &&
               (p.x < best_root.x - kPtEps ||
                (std::abs(p.x - best_root.x) <= kPtEps &&
                 p.y < best_root.y - kPtEps)))))) {
            have_best = true;
            best_root = p;
            best_cost = cost;
        }
    };

    for (const Segment& segment : segments) {
        const Point p = best_on_segment(segment, pa, pb, pc, m);
        consider(p, objective(p, pa, pb, pc));
    }
    solution.root = best_root;
    solution.cost = best_cost;
    return solution;
}

bool triangle_area_nonzero(const std::array<int, 3>& t,
                           const std::vector<ActiveNode>& active) {
    const Point a = active[static_cast<std::size_t>(t[0])].pos;
    const Point b = active[static_cast<std::size_t>(t[1])].pos;
    const Point c = active[static_cast<std::size_t>(t[2])].pos;
    const double area2 =
        (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    return std::abs(area2) > 1e-10;
}

std::vector<Triangle> build_delaunay_triangles(
    const std::vector<ActiveNode>& active) {
    if (active.size() < 3) {
        return {};
    }
    std::vector<double> coordinates;
    coordinates.reserve(active.size() * 2);
    for (const ActiveNode& node : active) {
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
            if (!triangle_area_nonzero(indices, active)) {
                continue;
            }
            // Canonical vertex order by ActiveNode id.
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
        const auto lid = sorted_ids(lhs, active);
        const auto rid = sorted_ids(rhs, active);
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

bool all_collinear(const std::vector<ActiveNode>& active) {
    if (active.size() <= 2) {
        return true;
    }
    const Point p0 = active[0].pos;
    std::size_t second = 1;
    while (second < active.size() &&
           same_point(p0, active[second].pos)) {
        ++second;
    }
    if (second >= active.size()) {
        return true;
    }
    const Point p1 = active[second].pos;
    const double dx = p1.x - p0.x;
    const double dy = p1.y - p0.y;
    const double span = std::max(std::abs(dx), std::abs(dy));
    if (span <= kPtEps) {
        return true;
    }
    const double tolerance = 1e-8 * std::max(1.0, span);
    for (const ActiveNode& node : active) {
        const double cross =
            dx * (node.pos.y - p0.y) - dy * (node.pos.x - p0.x);
        if (std::abs(cross) > tolerance) {
            return false;
        }
    }
    return true;
}

void compute_triangle_costs(std::vector<Triangle>& triangles,
                            const std::vector<ActiveNode>& active,
                            unsigned requested_threads, int& used_threads) {
    used_threads = 1;
    if (triangles.empty()) {
        return;
    }
    unsigned threads = requested_threads;
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
    }
    if (threads == 0) {
        threads = 1;
    }
    threads = std::min<unsigned>(threads,
                                 static_cast<unsigned>(triangles.size()));
    used_threads = static_cast<int>(threads);
    if (threads <= 1) {
        for (Triangle& triangle : triangles) {
            const ActiveNode& a =
                active[static_cast<std::size_t>(triangle.v[0])];
            const ActiveNode& b =
                active[static_cast<std::size_t>(triangle.v[1])];
            const ActiveNode& c =
                active[static_cast<std::size_t>(triangle.v[2])];
            const TernarySolution solution = compute_ternary_root(a, b, c);
            triangle.root = solution.root;
            triangle.cost = solution.cost;
        }
        return;
    }

    std::atomic<std::size_t> next{0};
    auto worker = [&]() {
        while (true) {
            const std::size_t i = next.fetch_add(1);
            if (i >= triangles.size()) {
                break;
            }
            Triangle& triangle = triangles[i];
            const ActiveNode& a =
                active[static_cast<std::size_t>(triangle.v[0])];
            const ActiveNode& b =
                active[static_cast<std::size_t>(triangle.v[1])];
            const ActiveNode& c =
                active[static_cast<std::size_t>(triangle.v[2])];
            const TernarySolution solution = compute_ternary_root(a, b, c);
            triangle.root = solution.root;
            triangle.cost = solution.cost;
        }
    };
    std::vector<std::thread> workers;
    workers.reserve(threads);
    for (unsigned i = 0; i < threads; ++i) {
        workers.emplace_back(worker);
    }
    for (std::thread& thread : workers) {
        thread.join();
    }
}

int max_wait_of(const Triangle& triangle,
                const std::vector<ActiveNode>& active) {
    return std::max({active[static_cast<std::size_t>(triangle.v[0])].wait_rounds,
                     active[static_cast<std::size_t>(triangle.v[1])].wait_rounds,
                     active[static_cast<std::size_t>(triangle.v[2])].wait_rounds});
}

void sort_triangles(std::vector<Triangle>& triangles,
                    const std::vector<ActiveNode>& active) {
    std::sort(triangles.begin(), triangles.end(),
              [&](const Triangle& lhs, const Triangle& rhs) {
                  if (std::abs(lhs.cost - rhs.cost) > kPtEps) {
                      return lhs.cost < rhs.cost;
                  }
                  const int lwait = max_wait_of(lhs, active);
                  const int rwait = max_wait_of(rhs, active);
                  if (lwait != rwait) {
                      return lwait > rwait;
                  }
                  return sorted_ids(lhs.v, active) <
                         sorted_ids(rhs.v, active);
              });
}

std::vector<Triangle> greedy_maximal_independent_set(
    const std::vector<Triangle>& triangles,
    std::size_t active_count) {
    std::vector<bool> used(active_count, false);
    std::vector<Triangle> selected;
    for (const Triangle& triangle : triangles) {
        const int a = triangle.v[0];
        const int b = triangle.v[1];
        const int c = triangle.v[2];
        if (!used[static_cast<std::size_t>(a)] &&
            !used[static_cast<std::size_t>(b)] &&
            !used[static_cast<std::size_t>(c)]) {
            selected.push_back(triangle);
            used[static_cast<std::size_t>(a)] = true;
            used[static_cast<std::size_t>(b)] = true;
            used[static_cast<std::size_t>(c)] = true;
        }
    }

    // Verify vertex-disjointness and maximality.
    for (std::size_t i = 0; i < selected.size(); ++i) {
        for (std::size_t j = i + 1; j < selected.size(); ++j) {
            for (int a : selected[i].v) {
                for (int b : selected[j].v) {
                    if (a == b) {
                        fail("PT-DME internal error: selected triangles share a vertex");
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
            fail("PT-DME internal error: independent set is not maximal");
        }
    }
    return selected;
}

std::string make_metadata(const std::vector<PtDmeRoundStat>& rounds,
                          int threads, double runtime_ms, int max_wait) {
    std::ostringstream out;
    out << "{\"variant\":\"PT-DME\","
        << "\"delaunay\":\"delaunator-cpp (MIT)\","
        << "\"delay_model\":\"linear-proxy\","
        << "\"threads\":" << threads << ","
        << "\"runtime_ms\":" << format_number(runtime_ms) << ","
        << "\"max_wait_rounds\":" << max_wait << ","
        << "\"rounds\":[";
    for (std::size_t i = 0; i < rounds.size(); ++i) {
        const PtDmeRoundStat& r = rounds[i];
        if (i) out << ',';
        out << "{\"round\":" << r.round
            << ",\"active_before\":" << r.active_before
            << ",\"duplicate_merges\":" << r.duplicate_merges
            << ",\"collinear\":" << (r.collinear ? "true" : "false")
            << ",\"triangles\":" << r.triangles
            << ",\"selected_triangles\":" << r.selected_triangles
            << ",\"ternary_merges\":" << r.ternary_merges
            << ",\"binary_fallback\":" << r.binary_fallback
            << ",\"carry_over\":" << r.carry_over
            << ",\"active_after\":" << r.active_after
            << ",\"max_wait_rounds\":" << r.max_wait_rounds << '}';
    }
    out << "]}";
    return out.str();
}

std::vector<LocalMergeResult> compute_selected_merges(
    const std::vector<Triangle>& selected,
    const std::vector<ActiveNode>& active, unsigned requested_threads,
    int& used_threads) {
    std::vector<LocalMergeResult> results(selected.size());
    used_threads = 1;
    if (selected.empty()) {
        return results;
    }
    unsigned threads = requested_threads;
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
    }
    if (threads == 0) {
        threads = 1;
    }
    threads = std::min<unsigned>(threads,
                                 static_cast<unsigned>(selected.size()));
    used_threads = static_cast<int>(threads);

    auto compute_one = [&](std::size_t i) {
        const Triangle& triangle = selected[i];
        const ActiveNode& a =
            active[static_cast<std::size_t>(triangle.v[0])];
        const ActiveNode& b =
            active[static_cast<std::size_t>(triangle.v[1])];
        const ActiveNode& c =
            active[static_cast<std::size_t>(triangle.v[2])];
        const double da = l1(triangle.root, a.pos);
        const double db = l1(triangle.root, b.pos);
        const double dc = l1(triangle.root, c.pos);
        LocalMergeResult& result = results[i];
        result.root = triangle.root;
        result.capacitance =
            a.capacitance + b.capacitance + c.capacitance;
        result.delay_min =
            std::min({a.delay_min + da, b.delay_min + db, c.delay_min + dc});
        result.delay_max =
            std::max({a.delay_max + da, b.delay_max + db, c.delay_max + dc});
    };

    if (threads <= 1) {
        for (std::size_t i = 0; i < selected.size(); ++i) {
            compute_one(i);
        }
        return results;
    }
    std::atomic<std::size_t> next{0};
    auto worker = [&]() {
        while (true) {
            const std::size_t i = next.fetch_add(1);
            if (i >= selected.size()) {
                break;
            }
            compute_one(i);
        }
    };
    std::vector<std::thread> workers;
    workers.reserve(threads);
    for (unsigned i = 0; i < threads; ++i) {
        workers.emplace_back(worker);
    }
    for (std::thread& thread : workers) {
        thread.join();
    }
    return results;
}

class PtDmeBuilder {
public:
    PtDmeBuilder(const std::vector<Sink>& sinks, PtDmeOptions options)
        : sinks_(sinks), options_(options) {}

    PtDmeResult run() {
        const auto start = std::chrono::steady_clock::now();
        initialize();

        if (!active_.empty()) {
            while (active_.size() > 2) {
                run_round();
            }
            finalize();
        }

        const auto end = std::chrono::steady_clock::now();
        const double runtime_ms =
            std::chrono::duration<double, std::milli>(end - start).count();
        result_.runtime_ms = runtime_ms;
        result_.topology = std::move(topology_);
        result_.max_wait_rounds = max_wait_;
        result_.topology.metadata_json =
            make_metadata(result_.rounds, result_.thread_count, runtime_ms,
                          max_wait_);
        return std::move(result_);
    }

private:
    const std::vector<Sink>& sinks_;
    PtDmeOptions options_;
    Topology topology_;
    std::vector<ActiveNode> active_;
    PtDmeResult result_;
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

        active_.reserve(sinks_.size());
        topology_.nodes.reserve(2 * sinks_.size());
        for (std::size_t i = 0; i < sinks_.size(); ++i) {
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
            active_.push_back(ActiveNode{
                static_cast<int>(i), Point{sink.x, sink.y}, sink.cap_ff,
                0.0, 0.0, node.id, 0});
            ++next_active_id_;
        }
    }

    Point binary_root(const ActiveNode& a, const ActiveNode& b) const {
        if (same_point(a.pos, b.pos)) {
            return a.pos;
        }
        const DmePoint da = dme_point(a.pos.x, a.pos.y);
        const DmePoint db = dme_point(b.pos.x, b.pos.y);
        const double distance = dme_l1(da, db);
        const DmeSegment core_a = dme_make_point_segment(da);
        const DmeSegment core_b = dme_make_point_segment(db);
        const DmeSegment merge = dme_merge_segment(
            core_a, 0.5 * distance, core_b, 0.5 * distance);
        const DmePoint root = dme_midpoint(merge.a, merge.b);
        return Point{root.x, root.y};
    }

    int make_tree_node(const Point& pos) {
        Node node;
        node.id = next_tree_id_++;
        node.x = pos.x;
        node.y = pos.y;
        node.kind = "internal";
        topology_.nodes.push_back(node);
        ++topology_.internal_nodes;
        return node.id;
    }

    void add_branch(int parent, int child) {
        topology_.edges.push_back(Edge{parent, child, "HV"});
    }

    ActiveNode merge_ternary(const ActiveNode& a, const ActiveNode& b,
                             const ActiveNode& c, const Point& root,
                             const LocalMergeResult* local = nullptr) {
        const int tree_root = make_tree_node(root);
        add_branch(tree_root, a.tree_root);
        add_branch(tree_root, b.tree_root);
        add_branch(tree_root, c.tree_root);
        ActiveNode node;
        node.id = next_active_id_++;
        node.pos = root;
        if (local != nullptr) {
            node.capacitance = local->capacitance;
            node.delay_min = local->delay_min;
            node.delay_max = local->delay_max;
        } else {
            const double da = l1(root, a.pos);
            const double db = l1(root, b.pos);
            const double dc = l1(root, c.pos);
            node.capacitance =
                a.capacitance + b.capacitance + c.capacitance;
            node.delay_min = std::min(
                {a.delay_min + da, b.delay_min + db, c.delay_min + dc});
            node.delay_max = std::max(
                {a.delay_max + da, b.delay_max + db, c.delay_max + dc});
        }
        node.tree_root = tree_root;
        node.wait_rounds = 0;
        max_wait_ = std::max(max_wait_, node.wait_rounds);
        return node;
    }

    ActiveNode merge_binary(const ActiveNode& a, const ActiveNode& b,
                            bool zero_length = false) {
        const Point root = zero_length ? a.pos : binary_root(a, b);
        const int tree_root = make_tree_node(root);
        add_branch(tree_root, a.tree_root);
        add_branch(tree_root, b.tree_root);
        const double da = l1(root, a.pos);
        const double db = l1(root, b.pos);
        ActiveNode node;
        node.id = next_active_id_++;
        node.pos = root;
        node.capacitance = a.capacitance + b.capacitance;
        node.delay_min = std::min(a.delay_min + da, b.delay_min + db);
        node.delay_max = std::max(a.delay_max + da, b.delay_max + db);
        node.tree_root = tree_root;
        node.wait_rounds = 0;
        max_wait_ = std::max(max_wait_, node.wait_rounds);
        return node;
    }

    void finalize() {
        if (active_.empty()) {
            fail("PT-DME produced an empty active set");
        }
        if (active_.size() == 2) {
            const ActiveNode merged = merge_binary(active_[0], active_[1]);
            active_ = {merged};
        }
        if (active_.size() != 1) {
            fail("PT-DME did not converge to one root");
        }
        const int root_id = active_[0].tree_root;
        if (root_id < 0 || static_cast<std::size_t>(root_id) >= topology_.nodes.size()) {
            fail("PT-DME internal error: invalid final root");
        }
        if (topology_.nodes[static_cast<std::size_t>(root_id)].kind == "internal") {
            topology_.nodes[static_cast<std::size_t>(root_id)].kind = "root";
        }
        topology_.root = root_id;
        if (static_cast<int>(topology_.edges.size()) !=
            static_cast<int>(topology_.nodes.size()) - 1) {
            fail("PT-DME internal error: output is not a tree");
        }
    }

    int handle_duplicate_coordinates() {
        if (active_.size() < 2) {
            return 0;
        }
        std::vector<std::size_t> order(active_.size());
        for (std::size_t i = 0; i < order.size(); ++i) {
            order[i] = i;
        }
        std::sort(order.begin(), order.end(), [&](std::size_t lhs,
                                                  std::size_t rhs) {
            const ActiveNode& a = active_[lhs];
            const ActiveNode& b = active_[rhs];
            if (std::abs(a.pos.x - b.pos.x) > kPtEps) {
                return a.pos.x < b.pos.x;
            }
            if (std::abs(a.pos.y - b.pos.y) > kPtEps) {
                return a.pos.y < b.pos.y;
            }
            return a.id < b.id;
        });

        bool any_duplicate = false;
        for (std::size_t i = 1; i < order.size(); ++i) {
            if (same_point(active_[order[i]].pos,
                           active_[order[i - 1]].pos)) {
                any_duplicate = true;
                break;
            }
        }
        if (!any_duplicate) {
            return 0;
        }

        std::vector<ActiveNode> next;
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
            ActiveNode current = active_[group.front()];
            for (std::size_t k = 1; k < group.size(); ++k) {
                current = merge_binary(current, active_[group[k]], true);
                ++merges;
            }
            next.push_back(current);
            i = j;
        }
        active_ = std::move(next);
        return merges;
    }

    int merge_collinear_round() {
        const double min_x =
            std::min_element(active_.begin(), active_.end(),
                             [](const ActiveNode& a, const ActiveNode& b) {
                                 return a.pos.x < b.pos.x;
                             })->pos.x;
        const double max_x =
            std::max_element(active_.begin(), active_.end(),
                             [](const ActiveNode& a, const ActiveNode& b) {
                                 return a.pos.x < b.pos.x;
                             })->pos.x;
        const double min_y =
            std::min_element(active_.begin(), active_.end(),
                             [](const ActiveNode& a, const ActiveNode& b) {
                                 return a.pos.y < b.pos.y;
                             })->pos.y;
        const double max_y =
            std::max_element(active_.begin(), active_.end(),
                             [](const ActiveNode& a, const ActiveNode& b) {
                                 return a.pos.y < b.pos.y;
                             })->pos.y;

        std::vector<std::size_t> order(active_.size());
        for (std::size_t i = 0; i < order.size(); ++i) {
            order[i] = i;
        }
        const bool primary_x = (max_x - min_x) >= (max_y - min_y);
        std::sort(order.begin(), order.end(), [&](std::size_t lhs,
                                                  std::size_t rhs) {
            const ActiveNode& a = active_[lhs];
            const ActiveNode& b = active_[rhs];
            const double a1 = primary_x ? a.pos.x : a.pos.y;
            const double b1 = primary_x ? b.pos.x : b.pos.y;
            const double a2 = primary_x ? a.pos.y : a.pos.x;
            const double b2 = primary_x ? b.pos.y : b.pos.x;
            if (std::abs(a1 - b1) > kPtEps) return a1 < b1;
            if (std::abs(a2 - b2) > kPtEps) return a2 < b2;
            return a.id < b.id;
        });

        std::vector<ActiveNode> next;
        next.reserve(active_.size());
        int merges = 0;
        std::size_t i = 0;
        for (; i + 2 < order.size(); i += 3) {
            const ActiveNode& a = active_[order[i]];
            const ActiveNode& b = active_[order[i + 1]];
            const ActiveNode& c = active_[order[i + 2]];
            const TernarySolution solution = compute_ternary_root(a, b, c);
            next.push_back(merge_ternary(a, b, c, solution.root));
            ++merges;
        }
        for (; i < order.size(); ++i) {
            ActiveNode node = active_[order[i]];
            ++node.wait_rounds;
            max_wait_ = std::max(max_wait_, node.wait_rounds);
            next.push_back(node);
        }
        active_ = std::move(next);
        return merges;
    }

    void binary_fallback(std::vector<ActiveNode>& next, int& used_count) {
        int best_i = 0;
        int best_j = 1;
        double best_distance = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < active_.size(); ++i) {
            for (std::size_t j = i + 1; j < active_.size(); ++j) {
                const double distance = l1(active_[i].pos, active_[j].pos);
                if (distance < best_distance - kPtEps ||
                    (std::abs(distance - best_distance) <= kPtEps &&
                     std::tie(active_[i].id, active_[j].id) <
                         std::tie(active_[best_i].id, active_[best_j].id))) {
                    best_distance = distance;
                    best_i = static_cast<int>(i);
                    best_j = static_cast<int>(j);
                }
            }
        }
        ActiveNode merged = merge_binary(active_[static_cast<std::size_t>(best_i)],
                                         active_[static_cast<std::size_t>(best_j)]);
        next.push_back(merged);
        used_count = 2;
        for (std::size_t i = 0; i < active_.size(); ++i) {
            if (static_cast<int>(i) == best_i ||
                static_cast<int>(i) == best_j) {
                continue;
            }
            ActiveNode node = active_[i];
            ++node.wait_rounds;
            max_wait_ = std::max(max_wait_, node.wait_rounds);
            next.push_back(node);
        }
    }

    void run_round() {
        PtDmeRoundStat stat;
        stat.round = static_cast<int>(result_.rounds.size());
        stat.active_before = static_cast<int>(active_.size());

        stat.duplicate_merges = handle_duplicate_coordinates();
        if (static_cast<int>(active_.size()) !=
            stat.active_before - stat.duplicate_merges) {
            fail("PT-DME internal error: duplicate active-count invariant");
        }
        const int base_active = static_cast<int>(active_.size());
        if (active_.size() <= 2) {
            stat.active_after = static_cast<int>(active_.size());
            stat.carry_over = stat.active_after;
            result_.rounds.push_back(stat);
            return;
        }

        if (all_collinear(active_)) {
            stat.collinear = true;
            const int before = static_cast<int>(active_.size());
            stat.ternary_merges = merge_collinear_round();
            stat.active_after = static_cast<int>(active_.size());
            stat.carry_over = stat.active_after - stat.ternary_merges;
            stat.max_wait_rounds = max_wait_;
            if (stat.active_after != before - 2 * stat.ternary_merges) {
                fail("PT-DME internal error: collinear active-count invariant");
            }
            result_.rounds.push_back(stat);
            return;
        }

        std::vector<Triangle> triangles = build_delaunay_triangles(active_);
        stat.triangles = static_cast<int>(triangles.size());
        if (triangles.empty()) {
            std::vector<ActiveNode> next;
            int used = 0;
            binary_fallback(next, used);
            stat.binary_fallback = 1;
            stat.carry_over = static_cast<int>(next.size()) - 1;
            stat.active_after = static_cast<int>(next.size());
            if (stat.active_after != base_active - 1) {
                fail("PT-DME internal error: binary fallback active-count invariant");
            }
            active_ = std::move(next);
            stat.max_wait_rounds = max_wait_;
            result_.rounds.push_back(stat);
            return;
        }

        int threads = 1;
        compute_triangle_costs(triangles, active_, options_.threads, threads);
        result_.thread_count = std::max(result_.thread_count, threads);
        sort_triangles(triangles, active_);

        const std::vector<Triangle> selected =
            greedy_maximal_independent_set(triangles, active_.size());
        stat.selected_triangles = static_cast<int>(selected.size());
        if (selected.empty()) {
            std::vector<ActiveNode> next;
            int used = 0;
            binary_fallback(next, used);
            stat.binary_fallback = 1;
            stat.carry_over = static_cast<int>(next.size()) - 1;
            stat.active_after = static_cast<int>(next.size());
            if (stat.active_after != base_active - 1) {
                fail("PT-DME internal error: binary fallback active-count invariant");
            }
            active_ = std::move(next);
            stat.max_wait_rounds = max_wait_;
            result_.rounds.push_back(stat);
            return;
        }

        int selected_threads = 1;
        const std::vector<LocalMergeResult> local_results =
            compute_selected_merges(selected, active_, options_.threads,
                                    selected_threads);
        result_.thread_count = std::max(result_.thread_count, selected_threads);

        std::vector<bool> used(active_.size(), false);
        std::vector<ActiveNode> next;
        next.reserve(active_.size() - 2 * selected.size());
        for (std::size_t i = 0; i < selected.size(); ++i) {
            const Triangle& triangle = selected[i];
            used[static_cast<std::size_t>(triangle.v[0])] = true;
            used[static_cast<std::size_t>(triangle.v[1])] = true;
            used[static_cast<std::size_t>(triangle.v[2])] = true;
            next.push_back(merge_ternary(
                active_[static_cast<std::size_t>(triangle.v[0])],
                active_[static_cast<std::size_t>(triangle.v[1])],
                active_[static_cast<std::size_t>(triangle.v[2])],
                triangle.root, &local_results[i]));
            ++stat.ternary_merges;
        }
        for (std::size_t i = 0; i < active_.size(); ++i) {
            if (used[i]) {
                continue;
            }
            ActiveNode node = active_[i];
            ++node.wait_rounds;
            max_wait_ = std::max(max_wait_, node.wait_rounds);
            next.push_back(node);
            ++stat.carry_over;
        }
        const int expected =
            base_active - 2 * static_cast<int>(selected.size());
        if (static_cast<int>(next.size()) != expected) {
            fail("PT-DME internal error: ternary active-count invariant");
        }
        active_ = std::move(next);
        stat.active_after = static_cast<int>(active_.size());
        stat.max_wait_rounds = max_wait_;
        result_.rounds.push_back(stat);
    }
};

}  // namespace

PtDmeResult build_pt_dme(const std::vector<Sink>& sinks,
                         const PtDmeOptions& options) {
    PtDmeBuilder builder(sinks, options);
    return builder.run();
}

}  // namespace cts

#ifndef PT_DME_NO_MAIN
int main(int argc, char** argv) {
    return cts::run_topology_cli(
        argc, argv, "pt_dme",
        [](const std::vector<cts::Sink>& sinks,
           const cts::BuildOptions&) {
            return cts::build_pt_dme(sinks).topology;
        },
        false, false, false);
}
#endif
