// Deterministic GH-tree regression tests.

#include "cts_common.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

namespace cts {
Topology build_gh_tree(const std::vector<Sink>& sinks, bool cap_aware);
}

namespace {

using cts::Sink;
using cts::Topology;

bool validate(const std::vector<Sink>& sinks, const Topology& topology,
              std::string& error) {
    if (sinks.empty()) {
        if (!topology.nodes.empty() || !topology.edges.empty()) {
            error = "empty input produced non-empty topology";
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
        error = "edge count is not nodes-1";
        return false;
    }
    if (topology.sink_nodes != static_cast<int>(sinks.size())) {
        error = "sink count mismatch";
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
        error = "root has incoming edge";
        return false;
    }
    for (std::size_t i = 0; i < indegree.size(); ++i) {
        if (static_cast<int>(i) == topology.root) continue;
        if (indegree[i] != 1) {
            error = "non-root indegree != 1";
            return false;
        }
    }
    std::vector<bool> visited(topology.nodes.size(), false);
    std::vector<int> stack{topology.root};
    while (!stack.empty()) {
        const int current = stack.back();
        stack.pop_back();
        if (visited[static_cast<std::size_t>(current)]) {
            error = "cycle";
            return false;
        }
        visited[static_cast<std::size_t>(current)] = true;
        for (const int child : children[static_cast<std::size_t>(current)]) {
            stack.push_back(child);
        }
    }
    if (std::count(visited.begin(), visited.end(), true) !=
        static_cast<long>(topology.nodes.size())) {
        error = "not all nodes reachable";
        return false;
    }
    std::set<int> sink_ids;
    for (const auto& node : topology.nodes) {
        if (node.kind == "sink") {
            if (!sink_ids.insert(node.sink_id).second) {
                error = "duplicate sink";
                return false;
            }
            if (!children[static_cast<std::size_t>(node.id)].empty()) {
                error = "sink is not leaf";
                return false;
            }
        }
    }
    if (sink_ids.size() != sinks.size()) {
        error = "sink coverage mismatch";
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

bool run_case(const std::string& name, const std::vector<Sink>& sinks) {
    const Topology topology = cts::build_gh_tree(sinks, false);
    std::string error;
    if (!validate(sinks, topology, error)) {
        std::cerr << "FAIL " << name << ": " << error << "\n";
        return false;
    }
    std::cout << "PASS " << name << ": sinks=" << sinks.size()
              << " nodes=" << topology.nodes.size()
              << " edges=" << topology.edges.size() << "\n";
    return true;
}

std::vector<Sink> random_sinks(int count, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> coordinate(0.0, 1000.0);
    std::vector<Sink> sinks;
    sinks.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
        sinks.push_back(Sink{i + 1, coordinate(rng), coordinate(rng), 20.0});
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

bool test_determinism() {
    const std::vector<Sink> sinks = random_sinks(500, 12345);
    const Topology a = cts::build_gh_tree(sinks, false);
    const Topology b = cts::build_gh_tree(sinks, false);
    if (a.nodes.size() != b.nodes.size() || a.edges.size() != b.edges.size()) {
        std::cerr << "FAIL determinism: size\n";
        return false;
    }
    for (std::size_t i = 0; i < a.nodes.size(); ++i) {
        if (a.nodes[i].id != b.nodes[i].id ||
            std::abs(a.nodes[i].x - b.nodes[i].x) > 1e-8 ||
            std::abs(a.nodes[i].y - b.nodes[i].y) > 1e-8 ||
            a.nodes[i].kind != b.nodes[i].kind) {
            std::cerr << "FAIL determinism: node\n";
            return false;
        }
    }
    for (std::size_t i = 0; i < a.edges.size(); ++i) {
        if (a.edges[i].source != b.edges[i].source ||
            a.edges[i].target != b.edges[i].target) {
            std::cerr << "FAIL determinism: edge\n";
            return false;
        }
    }
    std::cout << "PASS determinism: repeated builds identical\n";
    return true;
}

}  // namespace

int main() {
    bool ok = true;
    ok = test_small_counts() && ok;
    ok = run_case("A_random_1000", random_sinks(1000, 20261114)) && ok;

    std::vector<Sink> collinear;
    for (int i = 0; i < 1000; ++i) {
        collinear.push_back(Sink{i + 1, static_cast<double>(i), 500.0, 20.0});
    }
    ok = run_case("B_collinear_1000", collinear) && ok;

    std::vector<Sink> duplicates;
    for (int i = 0; i < 1000; ++i) {
        const int group = i % 8;
        duplicates.push_back(Sink{
            i + 1, 100.0 + 90.0 * (group % 4), 150.0 + 90.0 * (group / 4),
            20.0});
    }
    ok = run_case("C_duplicates_1000", duplicates) && ok;

    ok = test_determinism() && ok;
    std::cout << (ok ? "ALL GH-TREE TESTS PASS" : "GH-TREE TESTS FAILED")
              << "\n";
    return ok ? 0 : 1;
}
