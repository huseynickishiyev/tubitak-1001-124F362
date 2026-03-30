# Author: Huseyn Kishiyev
# Helper script to analyze complex output data from NIST Non-IID test results.
import pandas as pd
import re
from scipy.stats import pearsonr, spearmanr

def parse_nist_results(filepath):
    """Parses the output from the ea_non_iid tool."""
    with open(filepath, 'r') as f:
        text = f.read()

    estimator_map = {
        # Name from NIST tool output -> Official Test ID
        "Most Common Value Estimate (bit string)": "E1",
        "Collision Test Estimate (bit string)": "E2",
        "Markov Test Estimate (bit string)": "E3",
        "Compression Test Estimate (bit string)": "E4",
        "T-Tuple Test Estimate (bit string)": "E5",
        "LRS Test Estimate (bit string)": "E6",
        "Multi Most Common in Window (MultiMCW)": "E7",
        "Lag Prediction Test Estimate (bit string)": "E8",
        "Multi Markov Model with Counting (MultiMMC)": "E9",
        "LZ78Y Prediction Test Estimate": "E10"
    }
    
    all_results = []
    sequences = text.split('--- START ')
    
    for seq_block in sequences:
        if not seq_block.strip():
            continue
            
        current_sequence_results = {}
        lines = seq_block.split('\n')
        
        for line in lines:
            for est_name, est_id in estimator_map.items():
                if est_name in line and ("=" in line or ":" in line):
                    match = re.search(r'=\s*(\d+\.\d+)', line)
                    if match:
                        h_hat_value = float(match.group(1))
                        current_sequence_results[est_id] = h_hat_value
        
        if len(current_sequence_results) == 10:
             all_results.append(current_sequence_results)

    if not all_results:
        print("Error: Could not parse any complete data sets from 'results.txt'.")
        print("Please check that the file is not empty and contains valid NIST tool output.")
        return pd.DataFrame()

    df = pd.DataFrame(all_results)
    # Ensure columns are ordered E1, E2, ..., E10
    df = df[[f"E{i}" for i in range(1, 11)]]
    return df

def print_correlation_table(df, method='pearson', threshold=0.4):
    """Calculates and prints a styled correlation matrix."""
    print(f"\n## {method.capitalize()} Correlation Matrix")
    print("-" * (8 * len(df.columns) + 6))
    
    corr_matrix = df.corr(method=method)
    
    header = "      " + " ".join([f"{col:>6}" for col in corr_matrix.columns])
    print(header)
    
    for index, row in corr_matrix.iterrows():
        row_str = f"{index:<5}"
        for col_name, value in row.items():
            if abs(value) > threshold and index != col_name:
                row_str += f"\033[1m{value:>6.3f}\033[0m "
            else:
                row_str += f"{value:>6.3f} "
        print(row_str)

# --- Main execution ---
results_df = parse_nist_results('results.txt')

if not results_df.empty:
    print(f"Successfully parsed results for {len(results_df)} sequences.")
    
    print_correlation_table(results_df, method='pearson', threshold=0.4)
    print_correlation_table(results_df, method='spearman', threshold=0.4)
