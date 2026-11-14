// Deterministic PT-DME regression tests for docs/pt-dme.md section 11.

#include "pt_dme.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <utility>
#include <set>
#include <string>
#include <vector>

namespace {

using cts::Sink;
using cts::Topology;

struct Metrics {
    double wirelength = 0.0;
    double max_path = 0.0;
    double min_path = 0.0;
    double skew = 0.0;
};

Metrics metrics_of(const Topology& topology) {
    Metrics metrics;
    if (topology.nodes.empty() || topology.root < 0) {
        return metrics;
    }
    std::vector<std::vector<std::pair<int, double>>> adjacency(
        topology.nodes.size());
    for (const auto& edge : topology.edges) {
        const auto& a = topology.nodes[static_cast<std::size_t>(edge.source)];
        const auto& b = topology.nodes[static_cast<std::size_t>(edge.target)];
        const double length = std::abs(a.x - b.x) + std::abs(a.y - b.y);
        metrics.wirelength += length;
        adjacency[static_cast<std::size_t>(edge.source)].push_back(
            {edge.target, length});
    }
    double max_path = -1.0;
    double min_path = std::numeric_limits<double>::infinity();
    std::vector<std::pair<int, double>> stack{{topology.root, 0.0}};
    while (!stack.empty()) {
        const auto [node, distance] = stack.back();
        stack.pop_back();
        const auto& current = topology.nodes[static_cast<std::size_t>(node)];
        if (current.kind == "sink") {
            max_path = std::max(max_path, distance);
            min_path = std::min(min_path, distance);
        }
        for (const auto& child : adjacency[static_cast<std::size_t>(node)]) {
            stack.push_back({child.first, distance + child.second});
        }
    }
    if (max_path < 0.0) {
        return metrics;
    }
    metrics.max_path = max_path;
    metrics.min_path = min_path;
    metrics.skew = max_path - min_path;
    return metrics;
}

bool validate(const std::vector<Sink>& sinks, const cts::PtDmeResult& result,
              std::string& error) {
    const Topology& topology = result.topology;
    if (sinks.empty()) {
        if (!topology.nodes.empty() || !topology.edges.empty()) {
            error = "empty input produced a non-empty topology";
            return false;
        }
        return true;
    }
    if (topology.root < 0 ||
        static_cast<std::size_t>(topology.root) >= topology.nodes.size()) {
        error = "invalid root";
        return false;
    }
    if (static_cast<int>(topology.edges.size()) !=
        static_cast<int>(topology.nodes.size()) - 1) {
        error = "edge count is not node count - 1";
        return false;
    }
    if (topology.sink_nodes != static_cast<int>(sinks.size())) {
        error = "sink node count mismatch";
        return false;
    }

    std::vector<int> indegree(topology.nodes.size(), 0);
    std::vector<std::vector<int>> children(topology.nodes.size());
    for (const auto& edge : topology.edges) {
        if (edge.source < 0 || edge.target < 0 ||
            static_cast<std::size_t>(edge.source) >= topology.nodes.size() ||
            static_cast<std::size_t>(edge.target) >= topology.nodes.size()) {
            error = "edge endpoint out of range";
            return false;
        }
        ++indegree[static_cast<std::size_t>(edge.target)];
        children[static_cast<std::size_t>(edge.source)].push_back(edge.target);
    }
    if (indegree[static_cast<std::size_t>(topology.root)] != 0) {
        error = "root has an incoming edge";
        return false;
    }
    for (std::size_t i = 0; i < indegree.size(); ++i) {
        if (static_cast<int>(i) == topology.root) {
            continue;
        }
        if (indegree[i] != 1) {
            error = "non-root indegree is not 1";
            return false;
        }
    }

    std::vector<bool> visited(topology.nodes.size(), false);
    std::vector<int> stack{topology.root};
    while (!stack.empty()) {
        const int current = stack.back();
        stack.pop_back();
        if (visited[static_cast<std::size_t>(current)]) {
            error = "cycle or duplicate traversal";
            return false;
        }
        visited[static_cast<std::size_t>(current)] = true;
        for (const int child : children[static_cast<std::size_t>(current)]) {
            stack.push_back(child);
        }
    }
    if (std::count(visited.begin(), visited.end(), true) !=
        static_cast<long>(topology.nodes.size())) {
        error = "not all nodes are reachable";
        return false;
    }

    std::set<int> sink_ids;
    std::size_t sink_count = 0;
    for (const auto& node : topology.nodes) {
        if (node.kind == "sink") {
            ++sink_count;
            if (!sink_ids.insert(node.sink_id).second) {
                error = "duplicate sink id";
                return false;
            }
            if (!children[static_cast<std::size_t>(node.id)].empty()) {
                error = "sink node is not a leaf";
                return false;
            }
        }
    }
    if (sink_count != sinks.size()) {
        error = "sink node count does not match input";
        return false;
    }
    for (const Sink& sink : sinks) {
        if (sink_ids.count(sink.id) == 0) {
            error = "missing sink " + std::to_string(sink.id);
            return false;
        }
    }
    return true;
}

void print_rounds(const cts::PtDmeResult& result) {
    std::cout << "  rounds=" << result.rounds.size()
              << " threads=" << result.thread_count
              << " runtime_ms=" << result.runtime_ms
              << " max_wait=" << result.max_wait_rounds << "\n";
    for (const auto& round : result.rounds) {
        std::cout << "    r" << round.round << ": active "
                  << round.active_before << " -> " << round.active_after
                  << " dup=" << round.duplicate_merges
                  << " collinear=" << (round.collinear ? 1 : 0)
                  << " tri=" << round.triangles
                  << " sel=" << round.selected_triangles
                  << " ternary=" << round.ternary_merges
                  << " binary=" << round.binary_fallback
                  << " carry=" << round.carry_over << "\n";
    }
}

bool run_case(const std::string& name, const std::vector<Sink>& sinks) {
    const auto result = cts::build_pt_dme(sinks);
    std::string error;
    if (!validate(sinks, result, error)) {
        std::cerr << "FAIL " << name << ": " << error << "\n";
        return false;
    }
    std::cout << "PASS " << name << ": sinks=" << sinks.size()
              << " nodes=" << result.topology.nodes.size()
              << " edges=" << result.topology.edges.size();
    if (!sinks.empty()) {
        const Metrics metrics = metrics_of(result.topology);
        std::cout << " WL=" << metrics.wirelength
                  << " max=" << metrics.max_path
                  << " min=" << metrics.min_path
                  << " skew=" << metrics.skew;
    }
    std::cout << "\n";
    if (!sinks.empty()) {
        print_rounds(result);
    }
    return true;
}

std::vector<Sink> random_sinks(int count, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> coordinate(0.0, 1000.0);
    std::vector<Sink> sinks;
    sinks.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
        sinks.push_back(Sink{
            i + 1, coordinate(rng), coordinate(rng),
            10.0 + static_cast<double>(i % 20)});
    }
    return sinks;
}

std::vector<Sink> skewed_sinks(int count) {
    std::mt19937_64 rng(20261114);
    std::normal_distribution<double> dense(500.0, 8.0);
    std::uniform_real_distribution<double> far(800.0, 1000.0);
    std::vector<Sink> sinks;
    sinks.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
        double x;
        double y;
        if (i < 900) {
            x = std::clamp(dense(rng), 450.0, 550.0);
            y = std::clamp(dense(rng), 450.0, 550.0);
        } else {
            x = far(rng);
            y = far(rng);
        }
        sinks.push_back(Sink{i + 1, x, y, 20.0});
    }
    return sinks;
}

bool test_small_counts() {
    bool ok = true;
    for (int count = 0; count <= 5; ++count) {
        std::vector<Sink> sinks;
        for (int i = 0; i < count; ++i) {
            sinks.push_back(Sink{i + 1, 10.0 * i, 5.0 * i * i, 20.0});
        }
        ok = run_case("small_" + std::to_string(count), sinks) && ok;
    }
    return ok;
}

bool test_parallel_determinism() {
    const std::vector<Sink> sinks = random_sinks(300, 77);
    cts::PtDmeOptions serial_options;
    serial_options.threads = 1;
    cts::PtDmeOptions parallel_options;
    parallel_options.threads = 4;
    const cts::PtDmeResult serial =
        cts::build_pt_dme(sinks, serial_options);
    const cts::PtDmeResult parallel =
        cts::build_pt_dme(sinks, parallel_options);
    std::string error;
    if (!validate(sinks, serial, error) || !validate(sinks, parallel, error)) {
        std::cerr << "FAIL parallel_determinism: " << error << "\n";
        return false;
    }
    if (serial.topology.nodes.size() != parallel.topology.nodes.size() ||
        serial.topology.edges.size() != parallel.topology.edges.size()) {
        std::cerr << "FAIL parallel_determinism: size mismatch\n";
        return false;
    }
    for (std::size_t i = 0; i < serial.topology.nodes.size(); ++i) {
        const auto& a = serial.topology.nodes[i];
        const auto& b = parallel.topology.nodes[i];
        if (a.id != b.id || a.kind != b.kind ||
            std::abs(a.x - b.x) > 1e-8 || std::abs(a.y - b.y) > 1e-8 ||
            a.sink_id != b.sink_id) {
            std::cerr << "FAIL parallel_determinism: node mismatch\n";
            return false;
        }
    }
    for (std::size_t i = 0; i < serial.topology.edges.size(); ++i) {
        const auto& a = serial.topology.edges[i];
        const auto& b = parallel.topology.edges[i];
        if (a.source != b.source || a.target != b.target) {
            std::cerr << "FAIL parallel_determinism: edge mismatch\n";
            return false;
        }
    }
    std::cout << "PASS parallel_determinism: threads 1 vs 4 identical topology\n";
    return true;
}

}  // namespace

int main() {
    bool ok = true;
    ok = test_small_counts() && ok;
    ok = run_case("A_random_1000", random_sinks(1000, 20261114)) && ok;
    ok = test_parallel_determinism() && ok;
    ok = run_case("B_skewed_1000", skewed_sinks(1000)) && ok;

    std::vector<Sink> collinear;
    for (int i = 0; i < 1000; ++i) {
        collinear.push_back(Sink{i + 1, static_cast<double>(i), 500.0, 20.0});
    }
    ok = run_case("C_collinear_1000", collinear) && ok;

    std::vector<Sink> duplicates;
    for (int i = 0; i < 1000; ++i) {
        const int group = i % 8;
        duplicates.push_back(Sink{
            i + 1, 100.0 + 90.0 * (group % 4), 150.0 + 90.0 * (group / 4),
            20.0});
    }
    ok = run_case("D_duplicates_1000", duplicates) && ok;

    std::cout << (ok ? "ALL PT-DME TESTS PASS" : "PT-DME TESTS FAILED")
              << "\n";
    return ok ? 0 : 1;
}
