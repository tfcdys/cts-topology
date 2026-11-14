#pragma once

#include "cts_common.hpp"

#include <vector>

namespace cts {

struct HpDmeOptions {
    double alpha = 1.0;
    double beta = 1.0;
    double gamma = 0.2;
    double lambda = 0.2;
    double outlier_ratio_threshold = 2.5;
    int max_wait_rounds = 4;
    unsigned threads = 0;  // 0 = auto
};

struct HpDmeRoundStat {
    int round = 0;
    int active_before = 0;
    int duplicate_merges = 0;
    int triangles = 0;
    int selected_triangles = 0;
    int ternary_merges = 0;
    int four_terminal_merges = 0;
    int outliers = 0;
    int matched_outliers = 0;
    int fallback_binary = 0;
    int carry_over = 0;
    int active_after = 0;
    double average_outlier_wait = 0.0;
    int max_wait_rounds = 0;
};

struct HpDmeResult {
    Topology topology;
    std::vector<HpDmeRoundStat> rounds;
    double runtime_ms = 0.0;
    int thread_count = 0;
    int ternary_merges = 0;
    int four_terminal_merges = 0;
    int outliers_detected = 0;
    int matched_outliers = 0;
    int fallback_binary = 0;
    int max_wait_rounds = 0;
    double average_matched_outlier_wait = 0.0;
};

HpDmeResult build_hp_dme(const std::vector<Sink>& sinks,
                         const HpDmeOptions& options = {});

}  // namespace cts
