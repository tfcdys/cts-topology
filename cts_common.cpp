#include "cts_common.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>

namespace cts {

void fail(const std::string& message) {
    throw std::runtime_error(message);
}

std::string to_lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

std::string trim(const std::string& value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::vector<std::string> split(const std::string& line, char delimiter) {
    std::vector<std::string> fields;
    std::string field;
    std::istringstream input(line);
    while (std::getline(input, field, delimiter)) {
        fields.push_back(trim(field));
    }
    if (!line.empty() && line.back() == delimiter) {
        fields.emplace_back();
    }
    return fields;
}

std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (const char c : value) {
        switch (c) {
            case '\\': out << "\\\\"; break;
            case '"': out << "\\\""; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                        << static_cast<int>(static_cast<unsigned char>(c))
                        << std::dec << std::setfill(' ');
                } else {
                    out << c;
                }
        }
    }
    return out.str();
}

double parse_double(const std::string& text, const std::string& what) {
    if (text.empty()) {
        fail("empty numeric value for " + what);
    }
    try {
        std::size_t used = 0;
        const double value = std::stod(text, &used);
        if (used != text.size() || !std::isfinite(value)) {
            fail("invalid numeric value for " + what + ": " + text);
        }
        return value;
    } catch (const std::exception&) {
        fail("invalid numeric value for " + what + ": " + text);
    }
}

int parse_int(const std::string& text, const std::string& what) {
    try {
        std::size_t used = 0;
        const long long value = std::stoll(text, &used);
        if (used != text.size() || value < std::numeric_limits<int>::min() ||
            value > std::numeric_limits<int>::max()) {
            fail("invalid integer value for " + what + ": " + text);
        }
        return static_cast<int>(value);
    } catch (const std::exception&) {
        fail("invalid integer value for " + what + ": " + text);
    }
}

std::uint64_t parse_u64(const std::string& text, const std::string& what) {
    try {
        std::size_t used = 0;
        const unsigned long long value = std::stoull(text, &used);
        if (used != text.size()) {
            fail("invalid unsigned integer value for " + what + ": " + text);
        }
        return static_cast<std::uint64_t>(value);
    } catch (const std::exception&) {
        fail("invalid unsigned integer value for " + what + ": " + text);
    }
}

std::string format_number(double value, int precision) {
    if (std::abs(value) < 5e-13) {
        value = 0.0;
    }
    std::ostringstream out;
    out << std::setprecision(precision) << std::fixed << value;
    std::string text = out.str();
    const auto dot = text.find('.');
    if (dot != std::string::npos) {
        while (text.size() > dot + 2 && text.back() == '0') {
            text.pop_back();
        }
    }
    return text;
}

std::string case_stem_from_path(const fs::path& path) {
    std::string stem = path.stem().string();
    return stem.empty() ? "case" : stem;
}

fs::path default_case_path(const std::string& name) {
    return fs::path("data") / "cases" / (name + ".csv");
}

fs::path resolve_case_path(const std::string& value) {
    const fs::path supplied(value);
    if (supplied.has_extension() || supplied.has_parent_path()) {
        return supplied;
    }
    return default_case_path(value);
}

std::vector<Sink> read_case_csv(const fs::path& path) {
    std::ifstream input(path);
    if (!input) {
        fail("cannot open case file: " + path.string());
    }

    std::vector<Sink> sinks;
    std::string line;
    std::size_t line_number = 0;
    bool header_checked = false;
    int id_col = 0;
    int x_col = 1;
    int y_col = 2;
    int cap_col = 3;

    while (std::getline(input, line)) {
        ++line_number;
        line = trim(line);
        if (line.empty() || line[0] == '#') {
            continue;
        }
        const auto fields = split(line, ',');
        if (!header_checked) {
            header_checked = true;
            bool is_header = false;
            for (const auto& field : fields) {
                const std::string lowered = to_lower(field);
                if (lowered.find("sink") != std::string::npos ||
                    lowered == "id" || lowered.find("x") != std::string::npos) {
                    is_header = true;
                    break;
                }
            }
            if (is_header) {
                for (std::size_t i = 0; i < fields.size(); ++i) {
                    const std::string lowered = to_lower(fields[i]);
                    if (lowered == "sink_id" || lowered == "id" ||
                        lowered == "sink") {
                        id_col = static_cast<int>(i);
                    } else if (lowered == "x_um" || lowered == "x") {
                        x_col = static_cast<int>(i);
                    } else if (lowered == "y_um" || lowered == "y") {
                        y_col = static_cast<int>(i);
                    } else if (lowered == "cap_ff" || lowered == "cap" ||
                               lowered == "capacitance" ||
                               lowered == "capacitance_ff") {
                        cap_col = static_cast<int>(i);
                    }
                }
                if (line_number == 1) {
                    continue;
                }
            }
        }

        const int needed = std::max({id_col, x_col, y_col, cap_col});
        if (static_cast<int>(fields.size()) <= needed) {
            fail("case file " + path.string() + " line " +
                 std::to_string(line_number) + " has too few columns");
        }
        Sink s;
        s.id = parse_int(fields[static_cast<std::size_t>(id_col)], "sink id");
        s.x = parse_double(fields[static_cast<std::size_t>(x_col)], "sink x");
        s.y = parse_double(fields[static_cast<std::size_t>(y_col)], "sink y");
        s.cap_ff = parse_double(fields[static_cast<std::size_t>(cap_col)],
                                "sink capacitance");
        if (s.cap_ff < 0.0) {
            fail("negative sink capacitance in case file line " +
                 std::to_string(line_number));
        }
        sinks.push_back(s);
    }

    if (sinks.empty()) {
        fail("case file contains no sink rows: " + path.string());
    }

    std::set<int> ids;
    for (const auto& sink : sinks) {
        if (!ids.insert(sink.id).second) {
            fail("duplicate sink id " + std::to_string(sink.id) +
                 " in case file");
        }
    }
    return sinks;
}

void write_topology_json(const fs::path& path, const Topology& topology,
                         const BuildOptions& options) {
    if (!path.parent_path().empty()) {
        fs::create_directories(path.parent_path());
    }
    std::ofstream out(path);
    if (!out) {
        fail("cannot create topology file: " + path.string());
    }
    out << std::setprecision(10);
    out << "{\n";
    out << "  \"format_version\": " << kFormatVersion << ",\n";
    out << "  \"algorithm\": \"" << json_escape(options.algorithm) << "\",\n";
    out << "  \"case\": {\n";
    out << "    \"name\": \"" << json_escape(options.case_name) << "\",\n";
    out << "    \"file\": \"" << json_escape(options.case_file) << "\",\n";
    out << "    \"sink_count\": " << topology.sink_nodes << "\n";
    out << "  },\n";
    out << "  \"units\": {\"coordinate\": \"um\", \"capacitance\": \"fF\", "
           "\"length\": \"um\"},\n";
    out << "  \"capacitance_aware\": "
        << (options.capacitance_aware ? "true" : "false") << ",\n";
    out << "  \"capacitance_mode\": \""
        << (options.capacitance_aware ? "cap-weighted" : "uniform") << "\",\n";
    if (options.dme) {
        out << "  \"dme\": {\"delay_model\": \"linear\", "
               "\"topology\": \"nearest-neighbor\", \"skew_bound_um\": "
            << format_number(options.skew_bound_um) << "},\n";
    }
    if (!topology.metadata_json.empty()) {
        out << "  \"algorithm_metadata\": " << topology.metadata_json
            << ",\n";
    }
    out << "  \"root\": " << topology.root << ",\n";
    out << "  \"node_count\": " << topology.nodes.size() << ",\n";
    out << "  \"edge_count\": " << topology.edges.size() << ",\n";
    out << "  \"build_runtime_ms\": "
        << format_number(topology.build_runtime_ms) << ",\n";
    out << "  \"statistics\": {\"internal_nodes\": " << topology.internal_nodes
        << ", \"bend_nodes\": " << topology.bend_nodes
        << ", \"sink_nodes\": " << topology.sink_nodes << "},\n";

    out << "  \"nodes\": [\n";
    for (std::size_t i = 0; i < topology.nodes.size(); ++i) {
        const Node& node = topology.nodes[i];
        out << "    {\"id\": " << node.id
            << ", \"x\": " << format_number(node.x)
            << ", \"y\": " << format_number(node.y)
            << ", \"kind\": \"" << node.kind << "\"";
        if (node.kind == "sink") {
            out << ", \"sink_id\": " << node.sink_id
                << ", \"cap_ff\": " << format_number(node.cap_ff);
        }
        out << "}" << (i + 1 == topology.nodes.size() ? "\n" : ",\n");
    }
    out << "  ],\n";
    out << "  \"edges\": [\n";
    for (std::size_t i = 0; i < topology.edges.size(); ++i) {
        const Edge& edge = topology.edges[i];
        out << "    {\"source\": " << edge.source
            << ", \"target\": " << edge.target
            << ", \"route\": \"" << edge.route << "\"}"
            << (i + 1 == topology.edges.size() ? "\n" : ",\n");
    }
    out << "  ]\n";
    out << "}\n";
}

bool Arguments::has(const std::string& key) const {
    return values.count(key) != 0 || flags.count(key) != 0;
}

std::string Arguments::get(const std::string& key,
                           const std::string& fallback) const {
    const auto it = values.find(key);
    return it == values.end() ? fallback : it->second;
}

Arguments parse_arguments(int argc, char** argv, int start) {
    Arguments args;
    for (int i = start; i < argc; ++i) {
        std::string token = argv[i];
        if (token.rfind("--", 0) == 0) {
            token = token.substr(2);
            std::string value;
            const auto equals = token.find('=');
            if (equals != std::string::npos) {
                value = token.substr(equals + 1);
                token = token.substr(0, equals);
                args.values["--" + token] = value;
            } else if (i + 1 < argc &&
                       std::string(argv[i + 1]).rfind("--", 0) != 0) {
                args.values["--" + token] = argv[++i];
            } else {
                args.flags.insert("--" + token);
            }
        } else {
            args.positional.push_back(token);
        }
    }
    return args;
}

void require_no_positional(const Arguments& args) {
    if (!args.positional.empty()) {
        fail("unexpected positional argument: " + args.positional.front());
    }
}

int count_from_case_name(const std::string& name) {
    const auto pos = name.find_last_of('_');
    if (pos == std::string::npos) {
        return 0;
    }
    const std::string suffix = name.substr(pos + 1);
    if (suffix.empty()) {
        return 0;
    }
    for (const char c : suffix) {
        if (!std::isdigit(static_cast<unsigned char>(c))) {
            return 0;
        }
    }
    try {
        return std::stoi(suffix);
    } catch (const std::exception&) {
        return 0;
    }
}

Bounds parse_bounds(const std::string& text) {
    const auto fields = split(text, ',');
    if (fields.size() != 4) {
        fail("--bounds requires four comma-separated values: "
             "xmin,xmax,ymin,ymax");
    }
    Bounds bounds;
    bounds.xmin = parse_double(fields[0], "bounds xmin");
    bounds.xmax = parse_double(fields[1], "bounds xmax");
    bounds.ymin = parse_double(fields[2], "bounds ymin");
    bounds.ymax = parse_double(fields[3], "bounds ymax");
    if (!(bounds.xmin < bounds.xmax) || !(bounds.ymin < bounds.ymax)) {
        fail("--bounds requires xmin<xmax and ymin<ymax");
    }
    return bounds;
}

std::uint64_t default_seed(int count, const std::string& name) {
    std::uint64_t value = 1469598103934665603ULL;
    for (const char c : name) {
        value ^= static_cast<unsigned char>(c);
        value *= 1099511628211ULL;
    }
    value ^= static_cast<std::uint64_t>(count) * 0x9E3779B97F4A7C15ULL;
    value ^= value >> 30;
    value *= 0xBF58476D1CE4E5B9ULL;
    value ^= value >> 27;
    value *= 0x94D049BB133111EBULL;
    value ^= value >> 31;
    return value;
}

std::vector<Sink> generate_sinks(const GenerationConfig& config) {
    if (config.count <= 0) {
        fail("--count must be a positive integer");
    }
    if (!(config.bounds.xmin < config.bounds.xmax) ||
        !(config.bounds.ymin < config.bounds.ymax)) {
        fail("invalid bounds: require xmin<xmax and ymin<ymax");
    }
    if (config.cap_min_ff < 0.0 || config.cap_max_ff < config.cap_min_ff) {
        fail("invalid capacitance range");
    }

    std::mt19937_64 rng(config.seed);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::normal_distribution<double> normal(0.0, 1.0);

    const double width = config.bounds.xmax - config.bounds.xmin;
    const double height = config.bounds.ymax - config.bounds.ymin;
    const double margin = 0.04 * std::min(width, height);

    int cluster_count =
        std::clamp(static_cast<int>(std::lround(std::sqrt(config.count) / 2.0)),
                   3, 14);
    if (config.count < 100) {
        cluster_count = std::max(3, cluster_count);
    }

    struct Cluster {
        double x;
        double y;
        double sx;
        double sy;
    };
    std::vector<Cluster> clusters;
    clusters.reserve(static_cast<std::size_t>(cluster_count));
    for (int i = 0; i < cluster_count; ++i) {
        Cluster c;
        c.x = config.bounds.xmin + margin + unit(rng) * (width - 2.0 * margin);
        c.y = config.bounds.ymin + margin + unit(rng) * (height - 2.0 * margin);
        const double scale = 0.025 + 0.045 * unit(rng);
        c.sx = width * scale;
        c.sy = height * scale;
        clusters.push_back(c);
    }

    const double clustered_fraction = 0.72;
    std::vector<Sink> sinks;
    sinks.reserve(static_cast<std::size_t>(config.count));
    for (int i = 0; i < config.count; ++i) {
        double x;
        double y;
        if (unit(rng) < clustered_fraction) {
            const Cluster& c =
                clusters[static_cast<std::size_t>(rng() % clusters.size())];
            x = c.x + normal(rng) * c.sx;
            y = c.y + normal(rng) * c.sy;
        } else {
            x = config.bounds.xmin + margin + unit(rng) * (width - 2.0 * margin);
            y = config.bounds.ymin + margin + unit(rng) * (height - 2.0 * margin);
        }
        x = std::clamp(x, config.bounds.xmin, config.bounds.xmax);
        y = std::clamp(y, config.bounds.ymin, config.bounds.ymax);
        const double u = unit(rng);
        const double shaped = 0.35 + 0.65 * u;
        double cap = config.cap_min_ff +
                     (config.cap_max_ff - config.cap_min_ff) * shaped;
        if (unit(rng) < 0.08) {
            cap *= 1.35;
        }
        cap = std::clamp(cap, config.cap_min_ff, config.cap_max_ff);
        sinks.push_back(Sink{i + 1, x, y, cap});
    }
    return sinks;
}

void write_case_csv(const fs::path& path, const std::vector<Sink>& sinks) {
    if (!path.parent_path().empty()) {
        fs::create_directories(path.parent_path());
    }
    std::ofstream out(path);
    if (!out) {
        fail("cannot create case file: " + path.string());
    }
    out << "sink_id,x_um,y_um,cap_ff\n";
    out << std::setprecision(10);
    for (const auto& sink : sinks) {
        out << sink.id << ',' << format_number(sink.x) << ','
            << format_number(sink.y) << ',' << format_number(sink.cap_ff)
            << '\n';
    }
}

void write_case_meta(const fs::path& case_path, const GenerationConfig& config,
                     std::size_t written) {
    fs::path meta_path = case_path;
    meta_path.replace_extension(".meta.json");
    std::ofstream out(meta_path);
    if (!out) {
        fail("cannot create case metadata file: " + meta_path.string());
    }
    out << "{\n";
    out << "  \"format_version\": " << kFormatVersion << ",\n";
    out << "  \"name\": \"" << json_escape(config.name) << "\",\n";
    out << "  \"count\": " << written << ",\n";
    out << "  \"seed\": " << config.seed << ",\n";
    out << "  \"bounds_um\": {\"xmin\": " << format_number(config.bounds.xmin)
        << ", \"xmax\": " << format_number(config.bounds.xmax)
        << ", \"ymin\": " << format_number(config.bounds.ymin)
        << ", \"ymax\": " << format_number(config.bounds.ymax) << "},\n";
    out << "  \"capacitance_ff\": {\"min\": " << format_number(config.cap_min_ff)
        << ", \"max\": " << format_number(config.cap_max_ff) << "},\n";
    out << "  \"distribution\": \"72% gaussian clusters + 28% uniform "
           "background\",\n";
    out << "  \"file\": \"" << json_escape(case_path.generic_string()) << "\"\n";
    out << "}\n";
}

int generate_one(const GenerationConfig& config, const fs::path& out_path) {
    const auto sinks = generate_sinks(config);
    write_case_csv(out_path, sinks);
    write_case_meta(out_path, config, sinks.size());
    std::cout << "Generated case: " << out_path.string() << "\n";
    std::cout << "  sinks=" << sinks.size() << " seed=" << config.seed
              << " bounds=[" << config.bounds.xmin << ','
              << config.bounds.xmax << "]x[" << config.bounds.ymin << ','
              << config.bounds.ymax << "] cap=[" << config.cap_min_ff << ','
              << config.cap_max_ff << "] fF\n";
    return 0;
}

namespace {

void print_case_gen_usage() {
    std::cout
        << "CTS testcase generator\n\n"
        << "Usage:\n"
        << "  case_gen generate --case <name> --count <N> [options]\n"
        << "  case_gen generate-all [--outdir data/cases] [--force]\n\n"
        << "Generate options:\n"
        << "  --seed <u64>          deterministic RNG seed\n"
        << "  --bounds xmin,xmax,ymin,ymax   coordinate domain\n"
        << "  --cap-min <fF>        minimum sink capacitance (default 5)\n"
        << "  --cap-max <fF>        maximum sink capacitance (default 50)\n"
        << "  --out <file.csv>      output case file\n"
        << "  --force               overwrite an existing case\n";
}

int cmd_generate(const Arguments& args) {
    require_no_positional(args);
    std::string name = args.get("--case", args.get("--name"));
    const std::string count_text = args.get("--count");
    int count = 0;
    if (!count_text.empty()) {
        count = parse_int(count_text, "--count");
    }
    if (name.empty() && count > 0) {
        name = "sinks_" + std::to_string(count);
    }
    if (name.empty()) {
        fail("generate requires --case <name> or --count <N>");
    }
    if (count <= 0) {
        count = count_from_case_name(name);
    }
    if (count <= 0) {
        fail("generate requires a positive --count");
    }

    GenerationConfig config;
    config.name = name;
    config.count = count;
    config.seed = args.has("--seed") ? parse_u64(args.get("--seed"), "--seed")
                                     : default_seed(count, name);
    if (args.has("--bounds")) {
        config.bounds = parse_bounds(args.get("--bounds"));
    }
    if (args.has("--cap-min")) {
        config.cap_min_ff = parse_double(args.get("--cap-min"), "--cap-min");
    }
    if (args.has("--cap-max")) {
        config.cap_max_ff = parse_double(args.get("--cap-max"), "--cap-max");
    }

    const fs::path out_path = args.has("--out")
                                  ? fs::path(args.get("--out"))
                                  : default_case_path(name);
    if (fs::exists(out_path) && !args.has("--force")) {
        std::cout << "Case already exists (use --force to regenerate): "
                  << out_path.string() << "\n";
        return 0;
    }
    return generate_one(config, out_path);
}

int cmd_generate_all(const Arguments& args) {
    require_no_positional(args);
    const fs::path outdir = args.has("--outdir")
                                ? fs::path(args.get("--outdir"))
                                : fs::path("data/cases");
    const std::array<int, 6> counts{20, 50, 100, 1000, 10000, 100000};
    for (const int count : counts) {
        const std::string name = "sinks_" + std::to_string(count);
        const fs::path out_path = outdir / (name + ".csv");
        if (fs::exists(out_path) && !args.has("--force")) {
            std::cout << "Case already exists (use --force to regenerate): "
                      << out_path.string() << "\n";
            continue;
        }
        GenerationConfig config;
        config.name = name;
        config.count = count;
        config.seed = default_seed(count, name);
        generate_one(config, out_path);
    }
    return 0;
}

void print_build_usage(const std::string& algorithm, bool supports_cap_aware,
                       bool bounded_skew) {
    std::cout << "CTS " << algorithm << " topology builder\n\n"
              << "Usage:\n"
              << "  " << algorithm
              << " --case <name|file.csv> [--out results/" << algorithm
              << "_<case>.json]\n\n"
              << "Options:\n"
              << "  --case <name|file.csv>  case name under data/cases, or CSV path\n"
              << "  --case-file <file>      explicit case CSV path\n"
              << "  --out <file.json>       output topology JSON\n";
    if (supports_cap_aware) {
        std::cout << "  --cap-aware             enable capacitance-weighted "
                     "quadrisection (default off)\n";
    }
    if (bounded_skew) {
        std::cout << "  --skew-bound <um>       global skew bound (default 100)\n";
    }
}

}  // namespace

int run_topology_cli(int argc, char** argv, const std::string& algorithm,
                     const TopologyBuilder& builder, bool supports_cap_aware,
                     bool dme, bool bounded_skew) {
    try {
        if (argc < 2) {
            print_build_usage(algorithm, supports_cap_aware, bounded_skew);
            return 0;
        }
        const Arguments args = parse_arguments(argc, argv, 1);
        if (args.has("--help") || args.has("-h")) {
            print_build_usage(algorithm, supports_cap_aware, bounded_skew);
            return 0;
        }
        require_no_positional(args);

        const std::string case_value = args.get("--case");
        if (case_value.empty()) {
            fail("--case is required");
        }
        const fs::path case_path =
            args.has("--case-file") ? fs::path(args.get("--case-file"))
                                    : resolve_case_path(case_value);
        if (!fs::exists(case_path)) {
            fail("case file does not exist: " + case_path.string());
        }
        const bool capacitance_aware = args.has("--cap-aware");
        if (capacitance_aware && !supports_cap_aware) {
            fail("--cap-aware is not supported by " + algorithm +
                 "; capacitance awareness is off");
        }

        double skew_bound = 0.0;
        if (bounded_skew) {
            skew_bound =
                args.has("--skew-bound")
                    ? parse_double(args.get("--skew-bound"), "--skew-bound")
                    : 100.0;
            if (skew_bound < 0.0) {
                fail("--skew-bound must be non-negative");
            }
        } else if (args.has("--skew-bound")) {
            fail("--skew-bound is only valid for bst_dme");
        }

        const auto sinks = read_case_csv(case_path);
        const std::string case_name = case_stem_from_path(case_path);

        BuildOptions options;
        options.case_name = case_name;
        options.case_file = case_path.lexically_normal().generic_string();
        options.algorithm = algorithm;
        options.capacitance_aware = capacitance_aware;
        options.dme = dme;
        options.bounded_skew = bounded_skew;
        options.skew_bound_um = skew_bound;

        const auto build_start = std::chrono::steady_clock::now();
        Topology topology = builder(sinks, options);
        const auto build_end = std::chrono::steady_clock::now();
        topology.build_runtime_ms =
            std::chrono::duration<double, std::milli>(build_end - build_start)
                .count();
        const fs::path out_path =
            args.has("--out") ? fs::path(args.get("--out"))
                              : fs::path("results") /
                                    (algorithm + "_" + case_name + ".json");
        write_topology_json(out_path, topology, options);

        std::cout << "Built topology: " << out_path.string() << "\n";
        std::cout << "  algorithm=" << algorithm << " case=" << case_name
                  << " sinks=" << topology.sink_nodes
                  << " internal=" << topology.internal_nodes
                  << " bends=" << topology.bend_nodes
                  << " nodes=" << topology.nodes.size()
                  << " edges=" << topology.edges.size()
                  << " cap_aware=" << (capacitance_aware ? "true" : "false");
        if (bounded_skew) {
            std::cout << " skew_bound_um=" << format_number(skew_bound);
        }
        std::cout << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return 2;
    }
}

int run_case_gen_cli(int argc, char** argv) {
    try {
        if (argc < 2) {
            print_case_gen_usage();
            return 0;
        }
        const std::string command = argv[1];
        if (command == "-h" || command == "--help" || command == "help") {
            print_case_gen_usage();
            return 0;
        }
        const Arguments args = parse_arguments(argc, argv, 2);
        if (args.has("--help")) {
            print_case_gen_usage();
            return 0;
        }
        if (command == "generate" || command == "gen") {
            return cmd_generate(args);
        }
        if (command == "generate-all" || command == "gen-all") {
            return cmd_generate_all(args);
        }
        fail("unknown case_gen command '" + command + "'. Try --help");
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return 2;
    }
}

}  // namespace cts
