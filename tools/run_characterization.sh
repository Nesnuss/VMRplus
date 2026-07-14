#!/usr/bin/env bash
#
# Characterization-test harness for the VMR+ pipeline (safety net, option A).
# Run this ON THE MACHINE WHERE THE PIPELINE ACTUALLY WORKS (external binaries
# installed + NCBI reachable). It runs the pipeline IN THE FOREGROUND and waits,
# capturing every output, stdout, stderr and the exit code.
#
# Three modes:
#   baseline  Run the CURRENT VMR+_1.7.14.py TWICE and compare A vs B. This tells
#             you what is stable vs what is noise/NCBI drift BEFORE you trust any
#             golden files. Do this first, with the ORIGINAL checked out.
#   golden    Run once and save the result as the GOLDEN reference. Do this with
#             the ORIGINAL (pre-refactor) checked out.
#   check     Run once and compare against the saved golden. Do this after
#             checking out the REFACTORED version (same filename, same flags).
#
# Configure the variables below, then:
#   tools/run_characterization.sh baseline
#   tools/run_characterization.sh golden
#   tools/run_characterization.sh check
#
set -euo pipefail

# ── Configure these for your machine ────────────────────────────────────────
SCRIPT="${SCRIPT:-./VMR+_1.7.14.py}"          # the entry point (do NOT rename)
CONFIG="${CONFIG:-VMR_config_template.ini}"    # the -config .ini
THR="${THR:-10}"                               # -thr value
WORK="${WORK:-testes/VMR_char}"                # where this harness stores results
# ────────────────────────────────────────────────────────────────────────────

HERE="$(cd "$(dirname "$0")" && pwd)"
COMPARE="python3 $HERE/compare_runs.py"

# Read `output = ...` from the [general] section of the .ini.
OUTPUT_DIR="$(python3 - "$CONFIG" <<'PY'
import configparser, sys
c = configparser.ConfigParser()
c.read(sys.argv[1])
print(c["general"]["output"].strip())
PY
)"

run_once() {  # $1 = destination dir for this run's artifacts
    local dest="$1"
    rm -rf "$dest"
    mkdir -p "$dest"

    if [ -e "$OUTPUT_DIR" ]; then
        echo "ERROR: pipeline output dir already exists: $OUTPUT_DIR" >&2
        echo "Move or delete it so unique_dir() does not append a suffix." >&2
        exit 1
    fi

    echo ">> Running (foreground): python3 $SCRIPT -config $CONFIG -thr $THR"
    set +e
    python3 "$SCRIPT" -config "$CONFIG" -thr "$THR" \
        1>"$dest/stdout" 2>"$dest/stderr"
    local code=$?
    set -e
    echo "$code" > "$dest/exit_code"
    echo ">> exit code: $code (saved to $dest/exit_code)"

    # Move the pipeline output into this run's folder under a fixed name so the
    # comparator sees identical relative paths between runs.
    mv "$OUTPUT_DIR" "$dest/output"
    echo ">> artifacts in: $dest/  (output/, stdout, stderr, exit_code)"
}

compare_runs() {  # $1, $2 = two run dirs
    echo; echo "==== comparing $1 vs $2 ===="
    echo "-- exit codes:"; cat "$1/exit_code" "$2/exit_code"
    echo "-- output tree:"; $COMPARE "$1/output" "$2/output" || true
    echo "-- stdout:";      $COMPARE "$1/stdout" "$2/stdout" || true
    echo "-- stderr:";      $COMPARE "$1/stderr" "$2/stderr" || true
}

case "${1:-}" in
    baseline)
        run_once "$WORK/A"
        run_once "$WORK/B"
        compare_runs "$WORK/A" "$WORK/B"
        echo; echo "Review anything flagged above: it is NCBI drift or unhandled noise."
        ;;
    golden)
        run_once "$WORK/golden"
        echo; echo "Golden reference saved in $WORK/golden (keep it safe)."
        ;;
    check)
        run_once "$WORK/candidate"
        compare_runs "$WORK/golden" "$WORK/candidate"
        ;;
    *)
        echo "usage: $0 {baseline|golden|check}" >&2
        exit 2
        ;;
esac
