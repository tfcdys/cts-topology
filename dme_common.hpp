#pragma once

#include "cts_common.hpp"

#include <vector>

namespace cts {

constexpr double kDmeEps = 1e-8;

struct DmePoint {
    double x = 0.0;
    double y = 0.0;
};

struct DmeSegment {
    DmePoint a;
    DmePoint b;
    bool valid = false;
};

struct DmeNode {
    int left = -1;
    int right = -1;
    int parent = -1;
    int sink_index = -1;
    DmePoint pos;
    DmeSegment core;
    double delay = 0.0;
    double capacitance = 0.0;
    double radius_to_parent = 0.0;
};

struct DmeTree {
    std::vector<DmeNode> nodes;
    int root = -1;
};

DmePoint dme_point(double x, double y);
DmePoint dme_midpoint(const DmePoint& a, const DmePoint& b);
double dme_l1(const DmePoint& a, const DmePoint& b);
bool dme_same_point(const DmePoint& a, const DmePoint& b);
DmeSegment dme_make_segment(const DmePoint& a, const DmePoint& b);
DmeSegment dme_make_point_segment(const DmePoint& p);
DmePoint dme_closest_point_on_segment(const DmePoint& p,
                                      const DmeSegment& s);
double dme_point_segment_distance(const DmePoint& p, const DmeSegment& s);
double dme_segment_distance(const DmeSegment& a, const DmeSegment& b);
DmeSegment dme_merge_segment(const DmeSegment& a, double radius_a,
                             const DmeSegment& b, double radius_b);
DmePoint dme_choose_exact(const DmePoint& parent, const DmeSegment& core,
                          double radius);

DmeTree build_nearest_neighbor_tree(const std::vector<Sink>& sinks);
int dme_tree_depth(const DmeTree& tree);
Topology emit_dme_topology(const std::vector<Sink>& sinks,
                           const DmeTree& tree);

}  // namespace cts
