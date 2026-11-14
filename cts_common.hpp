#pragma once

#include <cstdint>
#include <filesystem>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace cts {

namespace fs = std::filesystem;

constexpr int kFormatVersion = 1;

struct Sink {
    int id = 0;
    double x = 0.0;
    double y = 0.0;
    double cap_ff = 0.0;
};

struct Node {
    int id = 0;
    double x = 0.0;
    double y = 0.0;
    std::string kind;  // root | internal | bend | sink
    int sink_id = -1;
    double cap_ff = 0.0;
};

struct Edge {
    int source = 0;
    int target = 0;
    std::string route = "HV";  // Manhattan route: HV or VH
};

struct Topology {
    std::vector<Node> nodes;
    std::vector<Edge> edges;
    int root = -1;
    int sink_nodes = 0;
    int internal_nodes = 0;
    int bend_nodes = 0;
    double build_runtime_ms = 0.0;
    // Optional JSON object emitted under "algorithm_metadata".
    std::string metadata_json;
};

struct Bounds {
    double xmin = 0.0;
    double xmax = 0.0;
    double ymin = 0.0;
    double ymax = 0.0;
};

struct BuildOptions {
    std::string case_name;
    std::string case_file;
    std::string algorithm = "htree";
    bool capacitance_aware = false;
    bool dme = false;
    bool bounded_skew = false;
    double skew_bound_um = 0.0;
};

[[noreturn]] void fail(const std::string& message);

std::string to_lower(std::string value);
std::string trim(const std::string& value);
std::vector<std::string> split(const std::string& line, char delimiter);
std::string json_escape(const std::string& value);
double parse_double(const std::string& text, const std::string& what);
int parse_int(const std::string& text, const std::string& what);
std::uint64_t parse_u64(const std::string& text, const std::string& what);
std::string format_number(double value, int precision = 9);

std::string case_stem_from_path(const fs::path& path);
fs::path default_case_path(const std::string& name);
fs::path resolve_case_path(const std::string& value);
std::vector<Sink> read_case_csv(const fs::path& path);
void write_topology_json(const fs::path& path, const Topology& topology,
                         const BuildOptions& options);

struct Arguments {
    std::map<std::string, std::string> values;
    std::set<std::string> flags;
    std::vector<std::string> positional;

    bool has(const std::string& key) const;
    std::string get(const std::string& key,
                    const std::string& fallback = "") const;
};

Arguments parse_arguments(int argc, char** argv, int start);
void require_no_positional(const Arguments& args);
int count_from_case_name(const std::string& name);
Bounds parse_bounds(const std::string& text);

struct GenerationConfig {
    std::string name;
    int count = 0;
    std::uint64_t seed = 0;
    Bounds bounds{0.0, 1000.0, 0.0, 1000.0};
    double cap_min_ff = 5.0;
    double cap_max_ff = 50.0;
};

std::uint64_t default_seed(int count, const std::string& name);
std::vector<Sink> generate_sinks(const GenerationConfig& config);
void write_case_csv(const fs::path& path, const std::vector<Sink>& sinks);
void write_case_meta(const fs::path& case_path, const GenerationConfig& config,
                     std::size_t written);
int generate_one(const GenerationConfig& config, const fs::path& out_path);
int run_case_gen_cli(int argc, char** argv);

using TopologyBuilder =
    std::function<Topology(const std::vector<Sink>&, const BuildOptions&)>;

int run_topology_cli(int argc, char** argv, const std::string& algorithm,
                     const TopologyBuilder& builder, bool supports_cap_aware,
                     bool dme, bool bounded_skew);

}  // namespace cts
