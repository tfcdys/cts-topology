#!/usr/bin/env bash
# Run one CTS case end-to-end: generate (if needed), build, evaluate, visualize.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

GEN_BIN="build/case_gen"
ALG_BIN=""
CASE_INPUT=""
TOPOLOGY="htree"
MODE="both"
FORMAT="png"
RUN_ALL=0
DO_EVAL=1
DO_VIZ=1
REGEN_CASE=0
COUNT=""
SEED=""
DPI="170"
CAP_AWARE=0
SKEW_BOUND=""

usage() {
    cat <<'EOF'
用法:
  ./run_cts.sh --case <case名称|case.csv> [--topology NAME] [选项]
  ./run_cts.sh <case名称> [topology]
  ./run_cts.sh --all [--topology NAME]

选项:
  -c, --case NAME       指定 case 名称（对应 data/cases/NAME.csv）
      --case-file FILE  直接指定 case CSV 路径
  -t, --topology NAME   拓扑算法: htree|zst_dme|bst_dme|pt_dme|hp_dme|gh_tree|sllt|salt
  -m, --mode MODE       可视化连线模式: euclidean|manhattan|both（默认 both）
      --format FMT      图片格式: png|svg|pdf（默认 png）
      --count N         case 不存在时自动生成的 sink 数量
      --seed N          自动生成 case 时使用的随机种子
      --all             依次运行 20/50/100/1000/10000/100000 全部标准 case
      --cap-aware       启用电容感知四叉划分（默认关闭，所有 sink 电容等权）
      --skew-bound N    BST-DME 的全局 skew 上界，单位 um（默认 100）
      --regen-case      强制重新生成 case（会覆盖原来的 CSV）
      --no-eval         跳过评估
      --no-viz          跳过可视化
      --dpi N           可视化 DPI（默认 170）
  -h, --help            显示帮助

示例:
  ./run_cts.sh --case sinks_1000 --topology htree
  ./run_cts.sh --case sinks_grid16 --topology htree
  ./run_cts.sh --case sinks_100 --topology zst_dme
  ./run_cts.sh --case sinks_100 --topology bst_dme --skew-bound 100
  ./run_cts.sh --case sinks_grid16 --topology pt_dme
  ./run_cts.sh --case sinks_grid16 --topology hp_dme
  ./run_cts.sh --case sinks_1000 --topology hp_dme
  ./run_cts.sh --case sinks_1000 --topology gh_tree
  ./run_cts.sh --case sinks_1000 --topology sllt
  ./run_cts.sh --case sinks_1000 --topology salt
  ./run_cts.sh sinks_20 htree --mode both
  ./run_cts.sh --all --topology htree --mode manhattan
EOF
}

die() {
    printf 'error: %s\n' "$*" >&2
    exit 1
}

# Prefer the project-local virtualenv when it exists.
if [[ -n "${PYTHON:-}" ]]; then
    PYTHON_BIN="$PYTHON"
elif [[ -x "$ROOT/.venv/bin/python" ]]; then
    PYTHON_BIN="$ROOT/.venv/bin/python"
else
    PYTHON_BIN="python3"
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        -c|--case)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            CASE_INPUT="$2"
            shift 2
            ;;
        --case-file)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            CASE_INPUT="$2"
            shift 2
            ;;
        -t|--topology)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            TOPOLOGY="$2"
            shift 2
            ;;
        -m|--mode)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            MODE="$2"
            shift 2
            ;;
        --format)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            FORMAT="$2"
            shift 2
            ;;
        --count)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            COUNT="$2"
            shift 2
            ;;
        --seed)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            SEED="$2"
            shift 2
            ;;
        --dpi)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            DPI="$2"
            shift 2
            ;;
        --all)
            RUN_ALL=1
            shift
            ;;
        --cap-aware)
            CAP_AWARE=1
            shift
            ;;
        --skew-bound)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            SKEW_BOUND="$2"
            shift 2
            ;;
        --regen-case)
            REGEN_CASE=1
            shift
            ;;
        --no-eval)
            DO_EVAL=0
            shift
            ;;
        --no-viz)
            DO_VIZ=0
            shift
            ;;
        --)
            shift
            break
            ;;
        -*)
            die "unknown option: $1 (use --help)"
            ;;
        *)
            if [[ -z "$CASE_INPUT" ]]; then
                CASE_INPUT="$1"
            elif [[ "$TOPOLOGY" == "htree" ]]; then
                TOPOLOGY="$1"
            else
                die "unexpected argument: $1"
            fi
            shift
            ;;
    esac
done

[[ -z "${1:-}" ]] || die "unexpected argument: $1"

case "$MODE" in
    euclidean|manhattan|both) ;;
    *) die "--mode must be euclidean, manhattan, or both" ;;
esac
case "$FORMAT" in
    png|svg|pdf) ;;
    *) die "--format must be png, svg, or pdf" ;;
esac

if [[ "$RUN_ALL" -eq 0 && -z "$CASE_INPUT" ]]; then
    usage
    exit 1
fi

case "$TOPOLOGY" in
    htree)   ALG_BIN="build/h_tree" ;;
    zst_dme) ALG_BIN="build/zst_dme" ;;
    bst_dme) ALG_BIN="build/bst_dme" ;;
    pt_dme)  ALG_BIN="build/pt_dme" ;;
    hp_dme)  ALG_BIN="build/hp_dme" ;;
    gh_tree) ALG_BIN="build/gh_tree" ;;
    sllt)    ALG_BIN="build/sllt" ;;
    salt)    ALG_BIN="build/salt" ;;
    *) die "unsupported topology '$TOPOLOGY' (use htree, zst_dme, bst_dme, pt_dme, hp_dme, gh_tree, sllt, or salt)" ;;
esac

REPORT_ONLY=0
case "$TOPOLOGY" in
    gh_tree|sllt|salt) REPORT_ONLY=1 ;;
esac

if [[ "$CAP_AWARE" -eq 1 && "$TOPOLOGY" != "htree" && "$TOPOLOGY" != "gh_tree" ]]; then
    die "--cap-aware is only supported by htree and gh_tree"
fi
if [[ "$REPORT_ONLY" -eq 1 && "$MODE" != "manhattan" ]]; then
    echo "[visualize] $TOPOLOGY evaluation policy: manhattan only"
    MODE="manhattan"
fi
if [[ -n "$SKEW_BOUND" && "$TOPOLOGY" != "bst_dme" ]]; then
    die "--skew-bound is only supported by bst_dme"
fi

echo "[build] compiling topology binaries"
make -s

if [[ "$DO_VIZ" -eq 1 ]] && ! "$PYTHON_BIN" -c 'import matplotlib' >/dev/null 2>&1; then
    cat >&2 <<EOF
error: matplotlib is not available for $PYTHON_BIN
       python3 -m venv .venv
       .venv/bin/pip install -r requirements.txt
       or set PYTHON=/path/to/python with matplotlib
EOF
    exit 2
fi

path_for_case() {
    local case_value="$1"
    if [[ "$case_value" == */* || "$case_value" == *.csv ]]; then
        printf '%s\n' "$case_value"
    else
        printf 'data/cases/%s.csv\n' "$case_value"
    fi
}

case_name_for() {
    local case_file="$1"
    local base
    base="$(basename "$case_file")"
    printf '%s\n' "${base%.csv}"
}

infer_count_from_name() {
    local name="$1"
    if [[ "$name" =~ ^sinks_([0-9]+)$ ]]; then
        printf '%s\n' "${BASH_REMATCH[1]}"
    elif [[ "$name" =~ ^case_([0-9]+)$ ]]; then
        printf '%s\n' "${BASH_REMATCH[1]}"
    else
        printf '0\n'
    fi
}

run_one() {
    local case_arg="$1"
    local case_file case_name topo_file eval_json eval_txt

    case_file="$(path_for_case "$case_arg")"
    case_name="$(case_name_for "$case_file")"

    if [[ "$REGEN_CASE" -eq 1 || ! -f "$case_file" ]]; then
        local n="$COUNT"
        if [[ -z "$n" ]]; then
            n="$(infer_count_from_name "$case_name")"
        fi
        [[ "$n" -gt 0 ]] || die \
            "case file not found: $case_file (provide --count N to generate it)"
        echo "[case] generating $case_name ($n sinks)"
        local gen_args=(generate --case "$case_name" --count "$n" --out "$case_file" --force)
        if [[ -n "$SEED" ]]; then
            gen_args+=(--seed "$SEED")
        fi
        "$GEN_BIN" "${gen_args[@]}"
    fi

    topo_file="results/${TOPOLOGY}_${case_name}.json"
    echo
    echo "============================================================"
    echo "[topology] case=$case_name algorithm=$TOPOLOGY cap_aware=$CAP_AWARE"
    echo "============================================================"
    local build_args=(--case "$case_file" --out "$topo_file")
    if [[ "$CAP_AWARE" -eq 1 ]]; then
        build_args+=(--cap-aware)
    fi
    if [[ -n "$SKEW_BOUND" ]]; then
        build_args+=(--skew-bound "$SKEW_BOUND")
    fi
    "$ALG_BIN" "${build_args[@]}"

    if [[ "$DO_EVAL" -eq 1 ]]; then
        eval_json=""
        eval_txt="results/${TOPOLOGY}_${case_name}.eval.txt"
        local eval_args=(
            scripts/evaluate_topology.py
            --topology "$topo_file"
            --case-file "$case_file"
            --report-out "$eval_txt"
        )
        if [[ "$REPORT_ONLY" -eq 1 ]]; then
            # Report-only topologies: no eval JSON, report + Manhattan image.
            eval_txt="results/${TOPOLOGY}_${case_name}_report.txt"
            eval_args=(
                scripts/evaluate_topology.py
                --topology "$topo_file"
                --case-file "$case_file"
                --report-out "$eval_txt"
            )
        else
            eval_json="results/${TOPOLOGY}_${case_name}.eval.json"
            eval_args+=(
                --json-out "$eval_json"
            )
        fi
        echo
        echo "[evaluate] $topo_file"
        set +e
        "$PYTHON_BIN" "${eval_args[@]}"
        local eval_status=$?
        set -e
        echo "[evaluate] report saved: $eval_txt"
        [[ "$eval_status" -eq 0 ]] || die "evaluation failed for $case_name (see $eval_txt)"
    fi

    if [[ "$DO_VIZ" -eq 1 ]]; then
        echo
        echo "[visualize] mode=$MODE format=$FORMAT"
        local viz_args=(
            scripts/visualize_topology.py
            --input "$topo_file"
            --mode "$MODE"
            --format "$FORMAT"
            --dpi "$DPI"
            --outdir results
        )
        if [[ "$DO_EVAL" -eq 1 && -n "$eval_json" && -f "$eval_json" ]]; then
            viz_args+=(--eval "$eval_json")
        fi
        MPLCONFIGDIR="$ROOT/.cache/mpl" "$PYTHON_BIN" "${viz_args[@]}"
    fi

    echo
    echo "[done] case=$case_name topology=$TOPOLOGY"
    echo "       topology : $topo_file"
    if [[ "$DO_EVAL" -eq 1 ]]; then
        echo "       report   : $eval_txt"
        grep -E 'Algorithm runtime|Total Manhattan wirelength|Max root-to-sink|Min root-to-sink|Max - min path' \
            "$eval_txt" | sed 's/^/       /' || true
    fi
    if [[ "$DO_VIZ" -eq 1 ]]; then
        if [[ "$MODE" == "both" ]]; then
            echo "       figures  : results/${TOPOLOGY}_${case_name}_euclidean.${FORMAT}"
            echo "                  results/${TOPOLOGY}_${case_name}_manhattan.${FORMAT}"
        else
            echo "       figure   : results/${TOPOLOGY}_${case_name}_${MODE}.${FORMAT}"
        fi
    fi
}

if [[ "$RUN_ALL" -eq 1 ]]; then
    for n in 20 50 100 1000 10000 100000; do
        run_one "sinks_$n"
    done
else
    run_one "$CASE_INPUT"
fi
