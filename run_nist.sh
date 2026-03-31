# Full pipeline to execute NIST SP800-90B entropy assessment (IID or Non-IID)
# over all sequences in the data directory and save results to a file.
#
# Usage:
#   ./run_nist.sh              # Non-IID assessment (default)
#   ./run_nist.sh --iid        # IID assessment
#   ./run_nist.sh --iid --bits 4   # IID with 4 bits-per-symbol (for gen_4bit sequences)
#   ./run_nist.sh --iid --bits 8   # IID with 8 bits-per-symbol (for gen_8bit sequences)
#
# Environment overrides:
#   NIST_TOOL=/path/to/ea_non_iid ./run_nist.sh
#   NIST_IID_TOOL=/path/to/ea_iid ./run_nist.sh --iid

set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"

# --- Defaults ---
MODE="non_iid"
BITS_PER_SYMBOL=1

# --- Argument parsing ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --iid)
            MODE="iid"
            shift
            ;;
        --bits)
            BITS_PER_SYMBOL="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: $0 [--iid] [--bits <1|4|8>]"
            exit 1
            ;;
    esac
done

# --- Tool paths ---
NIST_NON_IID_DEFAULT="$ROOT/nist_entropy/cpp/ea_non_iid"
NIST_IID_DEFAULT="$ROOT/nist_entropy/cpp/ea_iid"

NIST_NON_IID_TOOL="${NIST_TOOL:-$NIST_NON_IID_DEFAULT}"
NIST_IID_TOOL="${NIST_IID_TOOL:-$NIST_IID_DEFAULT}"

DATA_DIR="$ROOT/data"

# --- Select tool and output file based on mode ---
if [[ "$MODE" == "iid" ]]; then
    ACTIVE_TOOL="$NIST_IID_TOOL"
    RESULTS="$ROOT/results_iid.txt"
    echo "Mode: IID assessment (bits-per-symbol=$BITS_PER_SYMBOL)"
else
    ACTIVE_TOOL="$NIST_NON_IID_TOOL"
    RESULTS="$ROOT/results.txt"
    echo "Mode: Non-IID assessment"
fi

# --- Ensure data directory exists ---
if [[ ! -d "$DATA_DIR" ]]; then
    echo "Error: '$DATA_DIR' not found. Please run the C++ data generator first."
    exit 1
fi

# --- Build tool if missing ---
if [[ ! -x "$ACTIVE_TOOL" ]]; then
    echo "Binary not found at: $ACTIVE_TOOL"
    if [[ -f "$ROOT/nist_entropy/cpp/Makefile" ]]; then
        echo "Attempting to build in nist_entropy/cpp ..."
        if [[ "$MODE" == "iid" ]]; then
            make -C "$ROOT/nist_entropy/cpp" iid
        else
            make -C "$ROOT/nist_entropy/cpp" non_iid
        fi
    fi
fi

# --- Re-check after build attempt ---
if [[ ! -x "$ACTIVE_TOOL" ]]; then
    echo "Error: tool not found at: $ACTIVE_TOOL"
    if [[ "$MODE" == "iid" ]]; then
        echo "Run 'make iid' inside nist_entropy/cpp"
    else
        echo "Run 'make non_iid' inside nist_entropy/cpp"
    fi
    exit 1
fi

# --- Empty results file ---
: > "$RESULTS"

echo "Running NIST assessment with: $ACTIVE_TOOL"
echo "Results will be written to: $RESULTS"
echo ""

shopt -s nullglob
files=( "$DATA_DIR"/sequence_*.bin )
if (( ${#files[@]} == 0 )); then
    echo "No files matching $DATA_DIR/sequence_*.bin"
    exit 1
fi

# --- Process each sequence ---
for file in "${files[@]}"; do
    echo "Processing $file..."
    {
        echo "--- START ${file} ---"
        if [[ "$MODE" == "iid" ]]; then
            # -i: IID test  -a: treat file as binary (1 sample per byte)
            "$ACTIVE_TOOL" -i -a "$file" "$BITS_PER_SYMBOL"
        else
            "$ACTIVE_TOOL" -v -a "$file"
        fi
        echo "--- END ${file} ---"
    } >> "$RESULTS"
done

echo ""
echo "Done. Results -> $RESULTS"
