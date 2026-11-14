#pragma once

#include "cts_common.hpp"

#include <vector>

namespace cts {

struct PtDmeOptions {
    // 0 means auto-detect from std::thread::hardware_concurrency().
    unsigned threads = 0;
};

struct PtDmeRoundStat {
    int round = 0;
    int active_before = 0;
    int duplicate_merges = 0;
    int triangles = 0;
    int selected_triangles = 0;
    int ternary_merges = 0;
    int binary_fallback = 0;
    bool collinear = false;
    int carry_over = 0;
    int active_after = 0;
    int max_wait_rounds = 0;
};

struct PtDmeResult {
    Topology topology;
    std::vector<PtDmeRoundStat> rounds;
    double runtime_ms = 0.0;
    int thread_count = 0;
    int max_wait_rounds = 0;
};

PtDmeResult build_pt_dme(const std::vector<Sink>& sinks,
                         const PtDmeOptions& options = {});

}  // namespace cts
