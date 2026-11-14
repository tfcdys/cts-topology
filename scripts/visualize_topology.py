#!/usr/bin/env python3
"""Visualize a CTS clock-tree topology.

The script is topology-agnostic: any result following the documented JSON
schema can be rendered.  Two connection modes are supported:

  euclidean : draw the abstract tree edge as a straight segment
  manhattan : draw the abstract tree edge as an axis-aligned L/Z route

All sink points are always drawn; no point is down-sampled.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path
from typing import Any

# Keep Matplotlib's cache inside the project when the home cache is read-only.
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(__file__).resolve().parent.parent / ".cache" / "mpl")
)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D
except ImportError as error:  # pragma: no cover - depends on environment
    print(
        "error: matplotlib is required for visualization.\n"
        "       Install it with: python3 -m pip install matplotlib\n"
        f"       Import failure: {error}",
        file=sys.stderr,
    )
    raise SystemExit(2)


BG = "#0b1020"
PANEL = "#0f1730"
GRID = "#233252"
EDGE = "#68c2ff"
SINK_CMAP = "plasma"
ROOT = "#ffd166"
TEXT = "#e8eefc"
MUTED = "#8ea3c8"


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Visualize every sink and every edge of a CTS clock tree."
    )
    parser.add_argument("--input", required=True, help="topology JSON file")
    parser.add_argument(
        "--eval",
        default=None,
        help="optional evaluation JSON (default: auto-discover <stem>.eval.json)",
    )
    parser.add_argument(
        "--mode",
        choices=("euclidean", "manhattan", "both"),
        default="both",
        help="edge connection style (default: both)",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="output directory (default: directory of the topology file)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="explicit output path (only when --mode is not 'both')",
    )
    parser.add_argument(
        "--format",
        choices=("png", "svg", "pdf"),
        default="png",
        help="output image format (default: png)",
    )
    parser.add_argument("--dpi", type=int, default=170, help="raster DPI")
    parser.add_argument("--figsize", type=float, default=10.0, help="figure inches")
    parser.add_argument("--title", default=None, help="override plot title")
    parser.add_argument(
        "--no-colorbar", action="store_true", help="hide the sink-capacitance colorbar"
    )
    return parser.parse_args(argv)


def load_json(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as stream:
            data = json.load(stream)
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid JSON in {path}: {error}") from error
    if not isinstance(data, dict):
        raise ValueError("topology root must be a JSON object")
    if not isinstance(data.get("nodes"), list) or not isinstance(
        data.get("edges"), list
    ):
        raise ValueError("topology must contain 'nodes' and 'edges' arrays")
    return data


def node_kind(node: dict[str, Any]) -> str:
    return str(node.get("kind", node.get("type", ""))).lower()


class TopologyView:
    def __init__(self, data: dict[str, Any]) -> None:
        self.data = data
        self.algorithm = str(data.get("algorithm", "unknown"))
        case = data.get("case")
        self.case_name = (
            str(case.get("name", "case")) if isinstance(case, dict) else "case"
        )

        self.nodes: dict[int, dict[str, Any]] = {}
        for index, raw in enumerate(data["nodes"]):
            if not isinstance(raw, dict):
                raise ValueError(f"node #{index} is not an object")
            node_id = int(raw["id"])
            if node_id in self.nodes:
                raise ValueError(f"duplicate node id {node_id}")
            x = float(raw["x"])
            y = float(raw["y"])
            if not (math.isfinite(x) and math.isfinite(y)):
                raise ValueError(f"node {node_id} has non-finite coordinates")
            self.nodes[node_id] = {
                "id": node_id,
                "x": x,
                "y": y,
                "kind": node_kind(raw),
                "raw": raw,
            }
        if not self.nodes:
            raise ValueError("topology has no nodes")

        root_raw = data.get("root")
        if root_raw is None:
            roots = [
                node["id"] for node in self.nodes.values() if node["kind"] == "root"
            ]
            if len(roots) != 1:
                raise ValueError("cannot identify a unique root node")
            self.root = roots[0]
        else:
            self.root = int(root_raw)
        if self.root not in self.nodes:
            raise ValueError(f"root node {self.root} is missing")

        self.edges: list[tuple[int, int, str]] = []
        for index, raw in enumerate(data["edges"]):
            if not isinstance(raw, dict):
                raise ValueError(f"edge #{index} is not an object")
            source = int(raw.get("source", raw.get("s")))
            target = int(raw.get("target", raw.get("t")))
            if source not in self.nodes or target not in self.nodes:
                raise ValueError(f"edge #{index} references a missing node")
            route = str(raw.get("route", "HV")).upper()
            if route not in ("HV", "VH", "V", "H"):
                route = "HV"
            self.edges.append((source, target, route))

        self.sinks = [
            node
            for node in self.nodes.values()
            if node.get("kind") == "sink" or "sink_id" in node["raw"]
        ]
        self.sink_x = np.asarray([node["x"] for node in self.sinks], dtype=float)
        self.sink_y = np.asarray([node["y"] for node in self.sinks], dtype=float)
        self.sink_caps = np.asarray(
            [
                float(node["raw"].get("cap_ff", float("nan")))
                for node in self.sinks
            ],
            dtype=float,
        )
        self.has_caps = bool(
            len(self.sink_caps) and np.all(np.isfinite(self.sink_caps))
        )
        cap_flag = data.get("capacitance_aware")
        self.capacitance_aware = cap_flag if isinstance(cap_flag, bool) else None
        if self.has_caps:
            self.cap_min = float(np.min(self.sink_caps))
            self.cap_max = float(np.max(self.sink_caps))
            self.raw_caps_uniform = abs(self.cap_max - self.cap_min) <= 1e-9
        else:
            self.cap_min = 0.0
            self.cap_max = 0.0
            self.raw_caps_uniform = False
        # With cap-awareness off every sink has the same effective weight, so
        # the visualization deliberately uses one sink color.
        self.caps_uniform = (
            self.raw_caps_uniform or self.capacitance_aware is False
        )

    def segments(self, mode: str) -> np.ndarray:
        segments: list[list[tuple[float, float]]] = []
        for source, target, route in self.edges:
            a = self.nodes[source]
            b = self.nodes[target]
            x1, y1, x2, y2 = a["x"], a["y"], b["x"], b["y"]
            if mode == "euclidean":
                segments.append([(x1, y1), (x2, y2)])
                continue
            if route in ("V", "VH"):
                if abs(y1 - y2) > 1e-12:
                    segments.append([(x1, y1), (x1, y2)])
                if abs(x1 - x2) > 1e-12:
                    segments.append([(x1, y2), (x2, y2)])
            else:
                if abs(x1 - x2) > 1e-12:
                    segments.append([(x1, y1), (x2, y1)])
                if abs(y1 - y2) > 1e-12:
                    segments.append([(x2, y1), (x2, y2)])
            if abs(x1 - x2) <= 1e-12 and abs(y1 - y2) <= 1e-12:
                segments.append([(x1, y1), (x2, y2)])
        if not segments:
            return np.empty((0, 2, 2), dtype=float)
        return np.asarray(segments, dtype=float)


def discover_eval(input_path: Path, explicit: str | None) -> Path | None:
    if explicit:
        path = Path(explicit)
        return path if path.is_file() else None
    candidate = input_path.with_suffix(".eval.json")
    if candidate.is_file():
        return candidate
    candidate = input_path.parent / f"{input_path.stem}_eval.json"
    return candidate if candidate.is_file() else None


def load_metrics(path: Path | None) -> dict[str, Any]:
    if not path:
        return {}
    try:
        with path.open("r", encoding="utf-8") as stream:
            data = json.load(stream)
        return data if isinstance(data, dict) else {}
    except (OSError, json.JSONDecodeError):
        return {}


def render(
    view: TopologyView,
    mode: str,
    output: Path,
    args: argparse.Namespace,
    metrics: dict[str, Any],
) -> None:
    sink_count = len(view.sinks)
    density = max(1, sink_count)

    fig = plt.figure(figsize=(args.figsize, args.figsize), facecolor=BG)
    if view.has_caps and not view.caps_uniform and not args.no_colorbar:
        ax = fig.add_axes([0.065, 0.075, 0.765, 0.79])
    else:
        ax = fig.add_axes([0.075, 0.075, 0.86, 0.79])
    ax.set_facecolor(PANEL)

    edge_lw = float(np.clip(1.0 - 0.12 * math.log10(density), 0.34, 0.9))
    edge_alpha = float(np.clip(1.0 - 0.13 * math.log10(density), 0.32, 0.9))
    segments = view.segments(mode)
    if len(segments):
        collection = LineCollection(
            segments,
            colors=EDGE,
            linewidths=edge_lw,
            alpha=edge_alpha,
            zorder=2,
            antialiased=True,
            rasterized=(args.format == "png"),
        )
        ax.add_collection(collection)

    # Internal / bend points are only lightly indicated; sinks are the focus.
    internal = [
        node
        for node in view.nodes.values()
        if node["kind"] in ("internal", "bend", "root")
    ]
    if internal and sink_count <= 5000:
        ax.scatter(
            [node["x"] for node in internal],
            [node["y"] for node in internal],
            s=4.0,
            c="#9ec9ff",
            alpha=0.42,
            linewidths=0,
            zorder=3,
        )

    marker_size = (
        28.0
        if sink_count <= 100
        else 14.0
        if sink_count <= 1000
        else 5.0
        if sink_count <= 10000
        else 1.7
        if sink_count <= 50000
        else 1.0
    )
    if view.has_caps and not view.caps_uniform:
        scatter = ax.scatter(
            view.sink_x,
            view.sink_y,
            c=view.sink_caps,
            cmap=SINK_CMAP,
            s=marker_size,
            alpha=0.9 if sink_count <= 10000 else 0.82,
            linewidths=0,
            zorder=5,
            rasterized=(args.format == "png"),
        )
        if not args.no_colorbar:
            cax = fig.add_axes([0.858, 0.115, 0.022, 0.72])
            colorbar = fig.colorbar(scatter, cax=cax)
            colorbar.set_label("Sink capacitance (fF)", color=TEXT, fontsize=9)
            colorbar.ax.tick_params(colors=MUTED, labelsize=8)
            colorbar.outline.set_edgecolor(GRID)
    else:
        ax.scatter(
            view.sink_x,
            view.sink_y,
            c="#ff8c42",
            s=marker_size,
            alpha=0.88,
            linewidths=0,
            zorder=5,
            rasterized=(args.format == "png"),
        )

    root = view.nodes[view.root]
    ax.scatter(
        [root["x"]],
        [root["y"]],
        marker="*",
        s=230,
        c=ROOT,
        edgecolors="#fff3c4",
        linewidths=0.8,
        zorder=7,
    )

    xs = [node["x"] for node in view.nodes.values()]
    ys = [node["y"] for node in view.nodes.values()]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    xpad = max((xmax - xmin) * 0.035, 1.0)
    ypad = max((ymax - ymin) * 0.035, 1.0)
    ax.set_xlim(xmin - xpad, xmax + xpad)
    ax.set_ylim(ymin - ypad, ymax + ypad)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color=GRID, linewidth=0.55, alpha=0.55, zorder=0)
    ax.tick_params(colors=MUTED, labelsize=8)
    for spine in ax.spines.values():
        spine.set_color(GRID)
    ax.set_xlabel("X (um)", color=MUTED, fontsize=9)
    ax.set_ylabel("Y (um)", color=MUTED, fontsize=9)

    mode_label = "Euclidean" if mode == "euclidean" else "Manhattan (L-route)"
    title = args.title or (
        f"{view.case_name}  |  {view.algorithm}  |  {sink_count:,} sinks"
    )
    fig.text(
        0.065,
        0.968,
        title,
        color=TEXT,
        fontsize=15 if sink_count < 10000 else 13,
        fontweight="bold",
        ha="left",
        va="top",
    )

    info_bits = [f"connection: {mode_label}"]
    if view.capacitance_aware is False:
        info_bits.append("cap-aware: off (uniform weight)")
    elif view.capacitance_aware is True:
        info_bits.append("cap-aware: on")
    if view.has_caps and view.raw_caps_uniform:
        info_bits.append(f"caps: {view.cap_min:g} fF (uniform)")

    metric_bits: list[str] = []
    metrics_data = metrics.get("metrics") or {}
    if metrics_data:
        wl = metrics_data.get("total_manhattan_wirelength_um")
        skew = metrics_data.get("max_minus_min_um")
        if isinstance(wl, (int, float)):
            metric_bits.append(f"Manhattan WL: {wl / 1000.0:,.3f} mm")
        if isinstance(skew, (int, float)):
            metric_bits.append(f"max-min path: {skew:,.2f} um")
    eval_pass = metrics.get("pass")
    if isinstance(eval_pass, bool):
        metric_bits.append("evaluation: PASS" if eval_pass else "evaluation: FAIL")

    fig.text(
        0.065,
        0.928,
        "   •   ".join(info_bits),
        color=MUTED,
        fontsize=9,
        ha="left",
        va="top",
    )
    if metric_bits:
        fig.text(
            0.065,
            0.898,
            "   •   ".join(metric_bits),
            color=MUTED,
            fontsize=9,
            ha="left",
            va="top",
        )

    legend_handles = [
        Line2D([0], [0], color=EDGE, lw=1.1, label="clock tree"),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#ff8c42",
            markeredgecolor="none",
            markersize=5,
            label="sink",
        ),
        Line2D(
            [0],
            [0],
            marker="*",
            color="none",
            markerfacecolor=ROOT,
            markeredgecolor="none",
            markersize=10,
            label="clock root",
        ),
    ]
    legend = ax.legend(
        handles=legend_handles,
        loc="lower right",
        frameon=True,
        fontsize=8,
        facecolor=PANEL,
        edgecolor=GRID,
        labelcolor=TEXT,
    )
    legend.get_frame().set_alpha(0.9)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output,
        dpi=args.dpi,
        facecolor=fig.get_facecolor(),
    )
    plt.close(fig)
    print(f"Saved visualization: {output}")


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    input_path = Path(args.input)
    try:
        data = load_json(input_path)
        view = TopologyView(data)
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    eval_path = discover_eval(input_path, args.eval)
    metrics = load_metrics(eval_path)

    if args.mode == "both":
        if args.output:
            print("error: --output cannot be used with --mode both", file=sys.stderr)
            return 2
        modes = ("euclidean", "manhattan")
    else:
        modes = (args.mode,)

    outdir = Path(args.outdir) if args.outdir else input_path.parent
    for mode in modes:
        if args.output and len(modes) == 1:
            output = Path(args.output)
        else:
            output = outdir / f"{input_path.stem}_{mode}.{args.format}"
        render(view, mode, output, args, metrics)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
