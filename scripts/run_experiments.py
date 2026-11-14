#!/usr/bin/env python3
"""Run the common algorithm x case matrix and write a summary table.

Metrics and build runtime come from each topology binary / evaluator.  Build
runtime is measured in-process around the topology builder and is repeated five
times; the table reports the median.
"""

from __future__ import annotations

import csv
import json
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
CASES = ["sinks_grid16", "sinks_100", "sinks_1000"]
ALGORITHMS = {
    "htree": ("build/h_tree", []),
    "zst_dme": ("build/zst_dme", []),
    "bst_dme": ("build/bst_dme", ["--skew-bound", "100"]),
    "pt_dme": ("build/pt_dme", []),
    "hp_dme": ("build/hp_dme", []),
    "gh_tree": ("build/gh_tree", []),
    "sllt": ("build/sllt", []),
    "salt": ("build/salt", []),
}
REPEATS = 5


def run(command: list[str], check: bool = True) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        command,
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if check and completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"{completed.stderr}"
        )
    return completed


def main() -> int:
    rows: list[dict[str, object]] = []
    for algorithm, (binary, extra_args) in ALGORITHMS.items():
        for case in CASES:
            output = RESULTS / f"{algorithm}_{case}.json"
            case_file = f"data/cases/{case}.csv"
            timings: list[float] = []
            topology: dict[str, object] | None = None
            for _ in range(REPEATS):
                run([str(ROOT / binary), "--case", case_file,
                     "--out", str(output), *extra_args])
                topology = json.loads(output.read_text(encoding="utf-8"))
                timings.append(float(topology["build_runtime_ms"]))

            report_only = algorithm in {"gh_tree", "sllt", "salt"}
            if report_only:
                eval_json = ROOT / ".cache" / f"{algorithm}_{case}_eval.json"
                eval_txt = ROOT / ".cache" / f"{algorithm}_{case}_eval.txt"
                eval_json.parent.mkdir(parents=True, exist_ok=True)
            else:
                eval_json = RESULTS / f"{algorithm}_{case}.eval.json"
                eval_txt = RESULTS / f"{algorithm}_{case}.eval.txt"
            run([
                sys.executable,
                str(ROOT / "scripts/evaluate_topology.py"),
                "--topology", str(output),
                "--case-file", str(ROOT / case_file),
                "--json-out", str(eval_json),
                "--report-out", str(eval_txt),
                "--quiet",
            ], check=False)
            evaluation = json.loads(eval_json.read_text(encoding="utf-8"))
            checks = evaluation["checks"]
            metrics = evaluation.get("metrics") or {}
            if not metrics:
                raise RuntimeError(f"no metrics for: {algorithm} {case}")
            if not checks["all_sinks_covered"] or not checks["tree_complete"]:
                raise RuntimeError(
                    f"coverage/tree failure: {algorithm} {case}"
                )
            bound_ok = checks.get("within_dme_skew_bound")
            if evaluation["pass"]:
                status = "PASS"
            elif bound_ok is False:
                status = "SKEW_BOUND_FAIL"
            else:
                status = "FAIL"
            assert topology is not None
            row = {
                "algorithm": algorithm,
                "case": case,
                "sinks": evaluation["checks"]["sink_nodes"],
                "nodes": topology["node_count"],
                "edges": topology["edge_count"],
                "wirelength_um": metrics["total_manhattan_wirelength_um"],
                "max_path_um": metrics["max_root_to_sink_manhattan_um"],
                "min_path_um": metrics["min_root_to_sink_manhattan_um"],
                "skew_um": metrics["max_minus_min_um"],
                "runtime_ms": statistics.median(timings),
                "runtime_min_ms": min(timings),
                "runtime_max_ms": max(timings),
                "repeats": REPEATS,
                "status": status,
            }
            rows.append(row)
            if report_only:
                eval_json.unlink(missing_ok=True)
                eval_txt.unlink(missing_ok=True)
            print(
                f"{algorithm:8s} {case:14s} "
                f"WL={row['wirelength_um']:.6f} "
                f"max={row['max_path_um']:.6f} "
                f"min={row['min_path_um']:.6f} "
                f"skew={row['skew_um']:.6f} "
                f"runtime={row['runtime_ms']:.3f} ms status={row['status']}"
            )

    csv_path = RESULTS / "experiment_summary.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    label = {
        "htree": "H-tree",
        "zst_dme": "ZST-DME",
        "bst_dme": "BST-DME (B=100)",
        "pt_dme": "PT-DME",
        "hp_dme": "HP-DME",
        "gh_tree": "GH-tree",
        "sllt": "SLLT (BST+SALT)",
        "salt": "SALT",
    }
    lines = [
        "# 算法 × Case 实验汇总",
        "",
        "所有算法均关闭电容感知；BST-DME 的全局 skew 上界为 `100um`。"
        "运行时间为拓扑构建器内部耗时（不含 JSON 写出、评估和可视化），"
        f"重复 {REPEATS} 次取中位数。",
        "",
        "| 算法 | Case | Sinks | Nodes/Edges | 曼哈顿线长 (um) | 最大路径 (um) | 最小路径 (um) | 差值 (um) | 运行时间 (ms) | 状态 |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            f"| {label[str(row['algorithm'])]} | {row['case']} | {row['sinks']} | "
            f"{row['nodes']}/{row['edges']} | "
            f"{row['wirelength_um']:.3f} | {row['max_path_um']:.3f} | "
            f"{row['min_path_um']:.3f} | {row['skew_um']:.3f} | "
            f"{row['runtime_ms']:.3f} | {row['status']} |"
        )
    lines.extend([
        "",
        "说明：状态由评估器给出；`SKEW_BOUND_FAIL` 表示覆盖率/树结构通过，"
        "但实测 skew 超过该 DME 拓扑声明的上界。",
        "",
    ])
    md_path = RESULTS / "experiment_summary.md"
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")

    baseline_csv = RESULTS / "baseline_pt_dme" / "experiment_summary_before_hp.csv"
    if baseline_csv.is_file():
        baseline_rows = list(csv.DictReader(baseline_csv.open(encoding="utf-8")))
        baseline = {(r["algorithm"], r["case"]): r for r in baseline_rows}
        current_pt = {r["case"]: r for r in rows if r["algorithm"] == "pt_dme"}
        hp_rows = {r["case"]: r for r in rows if r["algorithm"] == "hp_dme"}
        comparison = []
        for case in CASES:
            pt = baseline.get(("pt_dme", case))
            hp = hp_rows.get(case)
            pt_now = current_pt.get(case)
            if pt is None or hp is None or pt_now is None:
                continue
            hp_json = json.loads(
                (RESULTS / f"hp_dme_{case}.json").read_text(encoding="utf-8")
            )
            meta = hp_json.get("algorithm_metadata", {})
            comparison.append({
                "case": case,
                "sinks": hp["sinks"],
                "pt_wirelength_um": float(pt["wirelength_um"]),
                "hp_wirelength_um": float(hp["wirelength_um"]),
                "wirelength_change_pct": 100.0 * (
                    float(hp["wirelength_um"]) / float(pt["wirelength_um"]) - 1.0
                ),
                "pt_max_path_um": float(pt["max_path_um"]),
                "hp_max_path_um": float(hp["max_path_um"]),
                "max_path_change_pct": 100.0 * (
                    float(hp["max_path_um"]) / float(pt["max_path_um"]) - 1.0
                ),
                "pt_min_path_um": float(pt["min_path_um"]),
                "hp_min_path_um": float(hp["min_path_um"]),
                "pt_skew_um": float(pt["skew_um"]),
                "hp_skew_um": float(hp["skew_um"]),
                "skew_change_pct": 100.0 * (
                    float(hp["skew_um"]) / float(pt["skew_um"]) - 1.0
                ),
                "pt_runtime_ms": float(pt_now["runtime_ms"]),
                "hp_runtime_ms": float(hp["runtime_ms"]),
                "pt_nodes": int(pt["nodes"]),
                "hp_nodes": int(hp["nodes"]),
                "ternary_merges": int(meta.get("ternary_merges", 0)),
                "four_terminal_merges": int(meta.get("four_terminal_merges", 0)),
                "outliers": int(meta.get("outliers_detected", 0)),
                "matched_outliers": int(meta.get("matched_outliers", 0)),
                "match_rate": float(meta.get("match_rate", 0.0)),
                "fallback_binary": int(meta.get("fallback_binary", 0)),
            })
        comparison_csv = RESULTS / "hp_vs_pt_comparison.csv"
        with comparison_csv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(comparison[0]))
            writer.writeheader()
            writer.writerows(comparison)
        lines = [
            "# HP-DME vs 原 PT-DME",
            "",
            "基线 PT-DME 结果保留在 `results/baseline_pt_dme/`。"
            "线长/最大路径/路径差取优化前 PT-DME 基线，运行时间取同一次最终实验中的 "
            "PT-DME 与 HP-DME，以避免不同批次计时噪声。",
            "",
            "| Case | PT WL | HP WL | WL 变化 | PT Max | HP Max | Max 变化 | PT Skew | HP Skew | Skew 变化 | PT Runtime | HP Runtime |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for row in comparison:
            lines.append(
                f"| {row['case']} | {row['pt_wirelength_um']:.3f} | "
                f"{row['hp_wirelength_um']:.3f} | "
                f"{row['wirelength_change_pct']:+.2f}% | "
                f"{row['pt_max_path_um']:.3f} | {row['hp_max_path_um']:.3f} | "
                f"{row['max_path_change_pct']:+.2f}% | "
                f"{row['pt_skew_um']:.3f} | {row['hp_skew_um']:.3f} | "
                f"{row['skew_change_pct']:+.2f}% | "
                f"{row['pt_runtime_ms']:.3f} | {row['hp_runtime_ms']:.3f} |"
            )
        comparison_md = RESULTS / "hp_vs_pt_comparison.md"
        comparison_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"Wrote {comparison_csv}")
        print(f"Wrote {comparison_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
