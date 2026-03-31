# Author: Huseyn Kishiyev
# Helper script to analyze output from the NIST SP800-90B IID test suite (ea_iid).
# Parses results_iid.txt produced by: ./run_nist.sh --iid

import re
import pandas as pd
from scipy.stats import pearsonr, spearmanr

# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def parse_iid_results(filepath):
    """
    Parses the output of the NIST ea_iid tool from results_iid.txt.

    The IID suite reports two categories of results per sequence:
      1. IID Assumption Tests — a pass/fail result for each of the 19 permutation
         tests that collectively assess whether the data is consistent with IID.
      2. Entropy Estimates — numerical min-entropy estimates from 6 estimators,
         identical in name to the non-IID suite but run under the IID assumption.

    Returns two DataFrames:
      - df_entropy : one row per sequence, columns = estimator IDs (E1..E6)
      - df_iid     : one row per sequence, columns = test names, values = pass/fail (1/0)
    """

    with open(filepath, 'r') as f:
        text = f.read()

    # --- Entropy estimator name -> short ID mapping ---
    estimator_map = {
        "Most Common Value Estimate":   "E1",
        "Collision Test Estimate":      "E2",
        "Markov Test Estimate":         "E3",
        "Compression Test Estimate":    "E4",
        "t-Tuple Test Estimate":        "E5",
        "LRS Test Estimate":            "E6",
    }

    # --- IID permutation test names (19 total per NIST SP800-90B Table 5) ---
    iid_tests = [
        "Excursion Test Statistic",
        "Number of Directional Runs",
        "Length of Directional Runs",
        "Number of Increases and Decreases",
        "Number of Runs Based on Median",
        "Length of Runs Based on Median",
        "Average Collision Test Statistic",
        "Maximum Collision Test Statistic",
        "Periodicity Test",
        "Covariance Test",
        "Compression Test",
        "LRS Test",
        "L1-Score",
        "Chi-square Independence Test",
        "LZ78Y Prediction Test",
        "MultiMMC Prediction Test",
        "Lag Prediction Test",
        "MultiMCW Prediction Test",
        "Compression Test (Bypass)",
    ]

    entropy_rows = []
    iid_rows = []

    sequences = text.split('--- START ')
    for block in sequences:
        if not block.strip():
            continue

        entropy_row = {}
        iid_row = {}

        for line in block.split('\n'):
            # --- Parse entropy estimates ---
            for name, eid in estimator_map.items():
                if name in line:
                    match = re.search(r'=\s*(\d+\.\d+)', line)
                    if match:
                        entropy_row[eid] = float(match.group(1))

            # --- Parse IID pass/fail results ---
            for test in iid_tests:
                if test in line:
                    # ea_iid outputs "passed" or "failed" on the same line
                    if re.search(r'\bpassed\b', line, re.IGNORECASE):
                        iid_row[test] = 1
                    elif re.search(r'\bfailed\b', line, re.IGNORECASE):
                        iid_row[test] = 0

        if len(entropy_row) == 6:
            entropy_rows.append(entropy_row)
        if iid_row:
            iid_rows.append(iid_row)

    if not entropy_rows:
        print("Error: Could not parse any entropy estimates from the file.")
        print("Check that results_iid.txt is not empty and contains valid ea_iid output.")
        return pd.DataFrame(), pd.DataFrame()

    df_entropy = pd.DataFrame(entropy_rows)
    df_entropy = df_entropy[[f"E{i}" for i in range(1, 7) if f"E{i}" in df_entropy.columns]]

    df_iid = pd.DataFrame(iid_rows) if iid_rows else pd.DataFrame()

    return df_entropy, df_iid


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

def print_correlation_matrix(df, method='pearson', threshold=0.4):
    """Prints a correlation matrix with bold highlighting above threshold."""
    print(f"\n## {method.capitalize()} Correlation Matrix (Entropy Estimates)")
    print("-" * (9 * len(df.columns) + 6))

    corr = df.corr(method=method)
    header = "      " + " ".join([f"{col:>7}" for col in corr.columns])
    print(header)

    for idx, row in corr.iterrows():
        row_str = f"{idx:<5} "
        for col, val in row.items():
            if abs(val) > threshold and idx != col:
                row_str += f"\033[1m{val:>7.3f}\033[0m "
            else:
                row_str += f"{val:>7.3f} "
        print(row_str)


def print_iid_summary(df_iid):
    """Prints pass rate per IID permutation test across all sequences."""
    if df_iid.empty:
        print("\nNo IID test results found — check that ea_iid was run with -i flag.")
        return

    print(f"\n## IID Permutation Test Pass Rates  (n={len(df_iid)} sequences)")
    print("-" * 60)
    print(f"{'Test':<45} {'Pass Rate':>10}  {'Result':>8}")
    print("-" * 60)

    all_passed = True
    for col in df_iid.columns:
        rate = df_iid[col].mean()
        verdict = "OK" if rate >= 0.99 else ("WARN" if rate >= 0.95 else "FAIL")
        if verdict != "OK":
            all_passed = False
        # Highlight failures and warnings
        if verdict == "FAIL":
            line = f"\033[1;31m{col:<45} {rate:>10.1%}  {verdict:>8}\033[0m"
        elif verdict == "WARN":
            line = f"\033[1;33m{col:<45} {rate:>10.1%}  {verdict:>8}\033[0m"
        else:
            line = f"{col:<45} {rate:>10.1%}  {verdict:>8}"
        print(line)

    print("-" * 60)
    if all_passed:
        print("Overall IID verdict: \033[1;32mPASS\033[0m — data is consistent with IID assumption.")
    else:
        print("Overall IID verdict: \033[1;31mFAIL\033[0m — one or more tests indicate non-IID behaviour.")


def print_entropy_summary(df):
    """Prints min/mean/max for each estimator and the final min-entropy (min of all)."""
    print(f"\n## Entropy Estimate Summary  (n={len(df)} sequences)")
    print("-" * 55)
    print(f"{'Estimator':<8} {'Min':>8} {'Mean':>8} {'Max':>8} {'Std':>8}")
    print("-" * 55)
    for col in df.columns:
        print(f"{col:<8} {df[col].min():>8.4f} {df[col].mean():>8.4f} "
              f"{df[col].max():>8.4f} {df[col].std():>8.4f}")
    print("-" * 55)

    # Final entropy = minimum across all estimators and sequences (conservative bound)
    final_h = df.min().min()
    limiting = df.min().idxmin()
    print(f"\nConservative min-entropy bound : \033[1m{final_h:.4f} bits\033[0m  (limited by {limiting})")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

df_entropy, df_iid = parse_iid_results('results_iid.txt')

if not df_entropy.empty:
    print(f"Successfully parsed IID results for {len(df_entropy)} sequences.\n")

    print_entropy_summary(df_entropy)
    print_iid_summary(df_iid)
    print_correlation_matrix(df_entropy, method='pearson',  threshold=0.4)
    print_correlation_matrix(df_entropy, method='spearman', threshold=0.4)
