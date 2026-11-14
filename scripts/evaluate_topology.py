#!/usr/bin/env python3
"""Evaluate a CTS clock-tree topology result.

The evaluator is deliberately independent of the topology generator.  It reads
the JSON topology, reconstructs the directed graph, and verifies:

  * every expected sink is present exactly once and is a reachable leaf,
  * the graph is a single rooted tree containing every node,
  * total Manhattan wirelength,
  * maximum/minimum root-to-sink Manhattan path length and their difference.

Only the Python standard library is required.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import deque
from pathlib import Path
from typing import Any


def finite_number(value: Any, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{what} is not numeric: {value!r}")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{what} is not finite: {value!r}")
    return result


def read_case_csv(path: Path) -> dict[int, dict[str, float]]:
    sinks: dict[int, dict[str, float]] = {}
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"case file has no header: {path}")
        fields = {name.strip().lower(): name for name in reader.fieldnames}
        id_key = fields.get("sink_id") or fields.get("id")
        x_key = fields.get("x_um") or fields.get("x")
        y_key = fields.get("y_um") or fields.get("y")
        cap_key = (
            fields.get("cap_ff")
            or fields.get("capacitance_ff")
            or fields.get("cap")
            or fields.get("capacitance")
        )
        if not all((id_key, x_key, y_key)):
            raise ValueError(
                f"case file is missing sink_id/x_um/y_um columns: {path}"
            )
        for row_number, row in enumerate(reader, start=2):
            try:
                sink_id = int(str(row[id_key]).strip())
                x = float(str(row[x_key]).strip())
                y = float(str(row[y_key]).strip())
                cap = float(str(row[cap_key]).strip()) if cap_key else 0.0
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"invalid numeric value in {path} line {row_number}: {error}"
                ) from error
            if sink_id in sinks:
                raise ValueError(f"duplicate sink id {sink_id} in {path}")
            if not (math.isfinite(x) and math.isfinite(y) and math.isfinite(cap)):
                raise ValueError(f"non-finite sink value in {path} line {row_number}")
            sinks[sink_id] = {"x": x, "y": y, "cap_ff": cap}
    if not sinks:
        raise ValueError(f"case file has no sink rows: {path}")
    return sinks


def resolve_case_file(
    topology: dict[str, Any], topology_path: Path, explicit: str | None
) -> tuple[Path | None, list[str]]:
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit))
    case_info = topology.get("case") or {}
    recorded = case_info.get("file")
    if isinstance(recorded, str) and recorded:
        candidates.append(Path(recorded))
    name = case_info.get("name")
    if isinstance(name, str) and name:
        candidates.extend(
            [
                Path("data") / "cases" / f"{name}.csv",
                topology_path.parent / f"{name}.csv",
                Path(__file__).resolve().parent.parent
                / "data"
                / "cases"
                / f"{name}.csv",
            ]
        )

    tried: list[str] = []
    seen: set[str] = set()
    for candidate in candidates:
        text = str(candidate)
        if text in seen:
            continue
        seen.add(text)
        tried.append(text)
        if candidate.is_file():
            return candidate.resolve(), tried
    return None, tried


def load_topology(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as stream:
            data = json.load(stream)
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid JSON in {path}: {error}") from error
    if not isinstance(data, dict):
        raise ValueError(f"topology root must be a JSON object: {path}")
    if not isinstance(data.get("nodes"), list) or not isinstance(
        data.get("edges"), list
    ):
        raise ValueError("topology must contain arrays named 'nodes' and 'edges'")
    return data


def node_kind(node: dict[str, Any]) -> str:
    return str(node.get("kind", node.get("type", ""))).lower()


def evaluate(
    topology: dict[str, Any],
    topology_path: Path,
    case_file: Path | None,
    expected_sinks: dict[int, dict[str, float]] | None,
) -> dict[str, Any]:
    nodes_raw = topology.get("nodes", [])
    edges_raw = topology.get("edges", [])
    checks: dict[str, Any] = {}
    errors: list[str] = []
    warnings: list[str] = []

    nodes: dict[int, dict[str, Any]] = {}
    sink_id_nodes: dict[int, dict[str, Any]] = {}
    duplicate_sink_ids: list[int] = []
    for index, raw in enumerate(nodes_raw):
        if not isinstance(raw, dict):
            raise ValueError(f"node #{index} is not an object")
        try:
            node_id = int(raw["id"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"node #{index} has an invalid id") from error
        if node_id in nodes:
            raise ValueError(f"duplicate node id {node_id}")
        x = finite_number(raw.get("x"), f"node {node_id} x")
        y = finite_number(raw.get("y"), f"node {node_id} y")
        kind = node_kind(raw)
        if not kind:
            raise ValueError(f"node {node_id} has no kind/type")
        entry: dict[str, Any] = {"id": node_id, "x": x, "y": y, "kind": kind}
        if kind == "sink" or "sink_id" in raw:
            if "sink_id" not in raw:
                errors.append(f"node {node_id} looks like a sink but has no sink_id")
            else:
                try:
                    entry["sink_id"] = int(raw["sink_id"])
                except (TypeError, ValueError) as error:
                    raise ValueError(
                        f"node {node_id} has an invalid sink_id"
                    ) from error
                if entry["sink_id"] in sink_id_nodes:
                    duplicate_sink_ids.append(entry["sink_id"])
                else:
                    sink_id_nodes[entry["sink_id"]] = entry
            if "cap_ff" in raw:
                entry["cap_ff"] = finite_number(
                    raw["cap_ff"], f"node {node_id} cap_ff"
                )
        nodes[node_id] = entry

    if not nodes:
        raise ValueError("topology has no nodes")

    root_value = topology.get("root")
    if root_value is None:
        roots = [node_id for node_id, node in nodes.items() if node["kind"] == "root"]
        if len(roots) != 1:
            raise ValueError("topology.root is missing and no unique root node exists")
        root = roots[0]
    else:
        try:
            root = int(root_value)
        except (TypeError, ValueError) as error:
            raise ValueError("topology.root is not an integer") from error
    if root not in nodes:
        raise ValueError(f"root node {root} does not exist")

    children: dict[int, list[int]] = {node_id: [] for node_id in nodes}
    parent: dict[int, int] = {}
    indegree: dict[int, int] = {node_id: 0 for node_id in nodes}
    edge_lengths: dict[tuple[int, int], float] = {}
    self_edges = 0
    duplicate_directed_edges = 0
    seen_edges: set[tuple[int, int]] = set()
    total_wirelength = 0.0

    for index, raw in enumerate(edges_raw):
        if not isinstance(raw, dict):
            raise ValueError(f"edge #{index} is not an object")
        source_raw = raw.get("source", raw.get("s"))
        target_raw = raw.get("target", raw.get("t"))
        try:
            source = int(source_raw)
            target = int(target_raw)
        except (TypeError, ValueError) as error:
            raise ValueError(f"edge #{index} has an invalid endpoint") from error
        if source not in nodes or target not in nodes:
            errors.append(
                f"edge #{index} references missing node(s): {source}->{target}"
            )
            continue
        if source == target:
            self_edges += 1
            continue
        key = (source, target)
        if key in seen_edges:
            duplicate_directed_edges += 1
            continue
        seen_edges.add(key)
        length = abs(nodes[source]["x"] - nodes[target]["x"]) + abs(
            nodes[source]["y"] - nodes[target]["y"]
        )
        edge_lengths[key] = length
        total_wirelength += length
        children[source].append(target)
        indegree[target] += 1
        parent[target] = source

    non_root_zero_indegree = [
        node_id
        for node_id, degree in indegree.items()
        if node_id != root and degree == 0
    ]
    multiple_parent = [
        node_id for node_id, degree in indegree.items() if degree > 1
    ]
    single_root = indegree[root] == 0 and not non_root_zero_indegree
    edge_count_ok = len(seen_edges) == len(nodes) - 1

    # Kahn's algorithm detects directed cycles even in a disconnected graph.
    kahn_indegree = dict(indegree)
    zero_queue: deque[int] = deque(
        node_id for node_id, degree in kahn_indegree.items() if degree == 0
    )
    processed = 0
    while zero_queue:
        current = zero_queue.popleft()
        processed += 1
        for child in children.get(current, []):
            kahn_indegree[child] -= 1
            if kahn_indegree[child] == 0:
                zero_queue.append(child)
    no_cycles = (
        processed == len(nodes)
        and not multiple_parent
        and not self_edges
        and not duplicate_directed_edges
    )

    reachable: set[int] = set()
    queue: deque[int] = deque([root])
    while queue:
        current = queue.popleft()
        if current in reachable:
            continue
        reachable.add(current)
        queue.extend(children[current])
    all_nodes_reachable = len(reachable) == len(nodes)

    sink_nodes = {
        node["sink_id"]: node for node in nodes.values() if "sink_id" in node
    }
    sink_node_ids = {node["id"] for node in sink_nodes.values()}
    sink_leaf_violations = [
        node["sink_id"]
        for node in sink_nodes.values()
        if children.get(node["id"])
    ]
    unreachable_sinks = [
        sink_id
        for sink_id, node in sink_nodes.items()
        if node["id"] not in reachable
    ]

    if expected_sinks is not None:
        expected_ids = set(expected_sinks)
        actual_ids = set(sink_nodes)
        missing_ids = sorted(expected_ids - actual_ids)
        extra_ids = sorted(actual_ids - expected_ids)
        coordinate_mismatches: list[int] = []
        capacitance_mismatches: list[int] = []
        for sink_id in sorted(expected_ids & actual_ids):
            expected = expected_sinks[sink_id]
            actual = sink_nodes[sink_id]
            if (
                abs(expected["x"] - actual["x"]) > 1e-6
                or abs(expected["y"] - actual["y"]) > 1e-6
            ):
                coordinate_mismatches.append(sink_id)
            actual_cap = actual.get("cap_ff")
            if actual_cap is None:
                capacitance_mismatches.append(sink_id)
            elif abs(expected["cap_ff"] - actual_cap) > 1e-6:
                capacitance_mismatches.append(sink_id)
    else:
        missing_ids = []
        extra_ids = []
        coordinate_mismatches = []
        capacitance_mismatches = []

    declared_sink_count: int | None = None
    case_info = topology.get("case")
    if isinstance(case_info, dict):
        try:
            if "sink_count" in case_info:
                declared_sink_count = int(case_info["sink_count"])
        except (TypeError, ValueError):
            declared_sink_count = None
    declared_count_ok = declared_sink_count is None or (
        declared_sink_count == len(sink_nodes)
    )

    all_sinks_covered = (
        not missing_ids
        and not extra_ids
        and not coordinate_mismatches
        and not capacitance_mismatches
        and not sink_leaf_violations
        and not unreachable_sinks
        and not duplicate_sink_ids
        and declared_count_ok
        and len(sink_nodes) > 0
    )

    tree_complete = (
        single_root
        and no_cycles
        and edge_count_ok
        and all_nodes_reachable
        and not errors
        and not sink_leaf_violations
        and not duplicate_sink_ids
        and not duplicate_directed_edges
        and not self_edges
        and len(sink_nodes) > 0
    )

    metrics: dict[str, Any] = {}
    if tree_complete:
        max_length = -math.inf
        min_length = math.inf
        max_sink: int | None = None
        min_sink: int | None = None
        path_sum = 0.0
        sink_paths = 0
        stack: list[tuple[int, float]] = [(root, 0.0)]
        while stack:
            current, distance = stack.pop()
            if current in sink_node_ids or nodes[current]["kind"] == "sink":
                sink_paths += 1
                path_sum += distance
                sink_id = nodes[current].get("sink_id", current)
                if distance > max_length + 1e-12:
                    max_length = distance
                    max_sink = int(sink_id)
                if distance < min_length - 1e-12:
                    min_length = distance
                    min_sink = int(sink_id)
            for child in children.get(current, []):
                edge = edge_lengths.get((current, child), 0.0)
                stack.append((child, distance + edge))
        if sink_paths:
            metrics = {
                "total_manhattan_wirelength_um": total_wirelength,
                "max_root_to_sink_manhattan_um": max_length,
                "max_path_sink_id": max_sink,
                "min_root_to_sink_manhattan_um": min_length,
                "min_path_sink_id": min_sink,
                "max_minus_min_um": max_length - min_length,
                "average_root_to_sink_manhattan_um": path_sum / sink_paths,
                "sink_path_count": sink_paths,
            }

    dme_info = topology.get("dme")
    skew_bound: float | None = None
    within_skew_bound: bool | None = None
    if isinstance(dme_info, dict):
        raw_bound = dme_info.get("skew_bound_um")
        if isinstance(raw_bound, (int, float)):
            skew_bound = float(raw_bound)
            if metrics:
                within_skew_bound = (
                    metrics["max_minus_min_um"] <= skew_bound + 1e-6
                )

    checks.update(
        {
            "case_file_found": expected_sinks is not None,
            "expected_sinks": len(expected_sinks) if expected_sinks is not None else None,
            "declared_sink_count": declared_sink_count,
            "sink_nodes": len(sink_nodes),
            "covered_sinks": len(set(expected_sinks or {}) & set(sink_nodes)),
            "missing_sink_ids": missing_ids[:20],
            "missing_sink_count": len(missing_ids),
            "extra_sink_ids": extra_ids[:20],
            "extra_sink_count": len(extra_ids),
            "coordinate_mismatch_ids": coordinate_mismatches[:20],
            "coordinate_mismatch_count": len(coordinate_mismatches),
            "capacitance_mismatch_ids": capacitance_mismatches[:20],
            "capacitance_mismatch_count": len(capacitance_mismatches),
            "duplicate_sink_ids": sorted(set(duplicate_sink_ids))[:20],
            "duplicate_sink_count": len(set(duplicate_sink_ids)),
            "unreachable_sink_ids": unreachable_sinks[:20],
            "unreachable_sink_count": len(unreachable_sinks),
            "sink_leaf_violations": sink_leaf_violations[:20],
            "sink_leaf_violation_count": len(sink_leaf_violations),
            "all_sinks_covered": all_sinks_covered,
            "root": root,
            "nodes": len(nodes),
            "edges": len(seen_edges),
            "expected_tree_edges": len(nodes) - 1,
            "rooted_tree": tree_complete,
            "tree_complete": tree_complete,
            "single_root": single_root,
            "no_cycles": no_cycles,
            "edge_count_ok": edge_count_ok,
            "all_nodes_reachable": all_nodes_reachable,
            "all_sinks_reachable": not unreachable_sinks,
            "self_edges": self_edges,
            "duplicate_directed_edges": duplicate_directed_edges,
            "errors": errors[:20],
            "dme_skew_bound_um": skew_bound,
            "within_dme_skew_bound": within_skew_bound,
        }
    )

    if expected_sinks is None:
        warnings.append(
            "source case CSV was not found; coverage was checked against the "
            "topology's declared sink count and sink nodes only"
        )

    passed = bool(
        all_sinks_covered
        and tree_complete
        and (within_skew_bound is not False)
    )
    return {
        "evaluator_version": 1,
        "topology_file": str(topology_path),
        "case_file": str(case_file) if case_file else None,
        "case": (topology.get("case") or {}).get("name")
        if isinstance(topology.get("case"), dict)
        else topology_path.stem,
        "topology": topology.get("algorithm", "unknown"),
        "capacitance_aware": topology.get("capacitance_aware"),
        "capacitance_mode": topology.get("capacitance_mode"),
        "dme": topology.get("dme") if isinstance(topology.get("dme"), dict) else None,
        "algorithm_runtime_ms": topology.get("build_runtime_ms"),
        "checks": checks,
        "metrics": metrics,
        "warnings": warnings,
        "pass": passed,
    }


def capability_label(result: dict[str, Any]) -> str:
    aware = result.get("capacitance_aware")
    if aware is True:
        return "on (cap-weighted)"
    if aware is False:
        return "off (uniform weight)"
    return "not specified"


def runtime_label(result: dict[str, Any]) -> str:
    runtime = result.get("algorithm_runtime_ms")
    if isinstance(runtime, (int, float)):
        return f"{float(runtime):.3f} ms"
    return "not specified"


def dme_label(result: dict[str, Any]) -> str:
    dme = result.get("dme")
    if not isinstance(dme, dict):
        return "n/a"
    model = str(dme.get("delay_model", "linear"))
    bound = dme.get("skew_bound_um")
    if isinstance(bound, (int, float)) and bound > 0:
        return f"{model}, topology={dme.get('topology', 'n/a')}, bound={bound:g} um"
    return f"{model}, topology={dme.get('topology', 'n/a')}, zero-skew"


def human_report(result: dict[str, Any]) -> str:
    checks = result["checks"]
    metrics = result["metrics"]
    status = "PASS" if result["pass"] else "FAIL"
    coverage = (
        f"{checks['covered_sinks']}/{checks['expected_sinks']}"
        if checks["expected_sinks"] is not None
        else f"{checks['sink_nodes']} sink nodes"
    )
    lines = [
        "=== CTS Topology Evaluation ===",
        f"Case                 : {result['case']}",
        f"Topology             : {result['topology']}",
        f"Capacitance aware    : {capability_label(result)}",
        f"DME embedding        : {dme_label(result)}",
        f"Topology file        : {result['topology_file']}",
        f"Case file            : {result['case_file'] or '(not found)'}",
        f"Algorithm runtime    : {runtime_label(result)}",
        f"Sink coverage        : {'PASS' if checks['all_sinks_covered'] else 'FAIL'} ({coverage})",
        f"Complete clock tree  : {'PASS' if checks['tree_complete'] else 'FAIL'}",
        f"Rooted / acyclic     : {checks['single_root']} / {checks['no_cycles']}",
        f"All nodes reachable  : {checks['all_nodes_reachable']}",
        f"Nodes / edges        : {checks['nodes']} / {checks['edges']} (tree requires {checks['expected_tree_edges']})",
    ]
    if checks.get("dme_skew_bound_um") is not None:
        within = checks.get("within_dme_skew_bound")
        lines.append(
            f"DME skew bound       : {checks['dme_skew_bound_um']:g} um "
            f"(within bound: {within})"
        )
    if checks["missing_sink_count"]:
        lines.append(
            f"Missing sinks        : {checks['missing_sink_count']} (e.g. {checks['missing_sink_ids']})"
        )
    if checks["extra_sink_count"]:
        lines.append(
            f"Extra sinks          : {checks['extra_sink_count']} (e.g. {checks['extra_sink_ids']})"
        )
    if checks["coordinate_mismatch_count"]:
        lines.append(
            f"Coordinate mismatch  : {checks['coordinate_mismatch_count']}"
        )
    if checks["capacitance_mismatch_count"]:
        lines.append(
            f"Capacitance mismatch : {checks['capacitance_mismatch_count']}"
        )
    if checks["errors"]:
        lines.append(f"Graph errors         : {checks['errors']}")

    if metrics:
        lines.extend(
            [
                "",
                f"Total Manhattan wirelength      : {metrics['total_manhattan_wirelength_um']:.6f} um",
                f"Max root-to-sink path (Manhattan): {metrics['max_root_to_sink_manhattan_um']:.6f} um (sink {metrics['max_path_sink_id']})",
                f"Min root-to-sink path (Manhattan): {metrics['min_root_to_sink_manhattan_um']:.6f} um (sink {metrics['min_path_sink_id']})",
                f"Max - min path difference        : {metrics['max_minus_min_um']:.6f} um",
                f"Average root-to-sink path        : {metrics['average_root_to_sink_manhattan_um']:.6f} um",
            ]
        )
    else:
        lines.append("")
        lines.append("Path metrics unavailable: root is not connected to sink nodes.")
    for warning in result["warnings"]:
        lines.append(f"Warning              : {warning}")
    lines.append("")
    lines.append(f"Result               : {status}")
    return "\n".join(lines) + "\n"


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Validate a CTS clock-tree topology and report path metrics."
    )
    parser.add_argument("--topology", required=True, help="topology JSON file")
    parser.add_argument(
        "--case-file",
        default=None,
        help="source case CSV (default: discover from topology metadata)",
    )
    parser.add_argument("--json-out", default=None, help="write machine-readable result")
    parser.add_argument("--report-out", default=None, help="write human-readable report")
    parser.add_argument("--quiet", action="store_true", help="do not print the report")
    args = parser.parse_args(argv)

    topology_path = Path(args.topology)
    try:
        topology = load_topology(topology_path)
        case_path, tried = resolve_case_file(topology, topology_path, args.case_file)
        expected_sinks = read_case_csv(case_path) if case_path else None
        result = evaluate(topology, topology_path, case_path, expected_sinks)
        if case_path is None:
            result["warnings"].append("case file candidates tried: " + ", ".join(tried))
        report = human_report(result)
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    if args.json_out:
        output = Path(args.json_out)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(result, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
    if args.report_out:
        output = Path(args.report_out)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(report, encoding="utf-8")
    if not args.quiet:
        sys.stdout.write(report)
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
