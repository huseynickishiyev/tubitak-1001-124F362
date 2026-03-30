# Full pipeline to execute the entire data feed from output dir to the NIST root so that we can
# easily feed finite number of sequences to test suite.
#!/usr/bin/env bash
set -euo pipefail

# Repo root = where this script lives
ROOT="$(cd "$(dirname "$0")" && pwd)"

# Prefer the binary built in cpp/
NIST_TOOL_DEFAULT="$ROOT/nist_entropy/cpp/ea_non_iid"
NIST_TOOL="${NIST_TOOL:-$NIST_TOOL_DEFAULT}"

RESULTS="$ROOT/results.txt"
DATA_DIR="$ROOT/data"

# Ensure data exists
if [[ ! -d "$DATA_DIR" ]]; then
  echo "Error: '$DATA_DIR' not found. Please run the C++ data generator first."
  exit 1
fi

# Build if missing
if [[ ! -x "$NIST_TOOL" ]]; then
  echo "Binary not found at: $NIST_TOOL"
  if [[ -f "$ROOT/nist_entropy/cpp/Makefile" ]]; then
    echo "Building ea_non_iid in nist_entropy/cpp ..."
    make -C "$ROOT/nist_entropy/cpp" non_iid
  fi
fi

# Re-check after build
if [[ ! -x "$NIST_TOOL" ]]; then
  echo "Error: still no ea_non_iid at: $NIST_TOOL"
  echo "Tip: run 'make non_iid' inside nist_entropy/cpp"
  exit 1
fi

# Empty results file
: > "$RESULTS"

echo "Running NIST non-IID assessment with: $NIST_TOOL"
shopt -s nullglob
files=( "$DATA_DIR"/sequence_*.bin )
if (( ${#files[@]} == 0 )); then
  echo "No files like $DATA_DIR/sequence_*.bin"
  exit 1
fi

for file in "${files[@]}"; do
  echo "Processing $file..."
  {
    echo "--- START ${file} ---"
    "$NIST_TOOL" -v -a "$file"
    echo "--- END ${file} ---"
  } >> "$RESULTS"
done

echo "Done. Results -> $RESULTS"
