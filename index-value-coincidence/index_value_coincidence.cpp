#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <string>

/**
 * Index-Value Coincidence Estimator
 * Based on: "Estimating Entropy via Index-Value Coincidence" (Aslan et al.)
 * Authors: Huseyn Kishiyev & Osman Furkan Atakul
 *
 * This class implements the logic to estimate min-entropy by analyzing 
 * the coincidence of sequence values with a cyclic alphabet assignment.
 */
class IndexValueCoincidence {
public:
    // Structure to hold the estimation results
    struct Result {
        double min_entropy;
        std::vector<double> probabilities; // p_i values
        std::vector<double> q_values;      // Intermediate q_i estimates
        long long total_terminal_points;   // t
        long long exhausting_points;       // t_EP
    };

    /**
     * Estimates the entropy of a given sequence.
     * * @param sequence The input data (integers representing alphabet indices).
     * @param alphabet_size The size of the alphabet 'k'.
     * For binary, k=2. For 8-bit, k=256.
     * @return Result struct containing entropy and stats.
     */
    static Result estimate(const std::vector<int>& sequence, int alphabet_size) {
        
        
        
        // Counts for each alphabet symbol (t_1 ... t_k)
        // Note: Using 0-based indexing for implementation, so counts[0] maps to paper's t_1 (or t_x1)
        std::vector<long long> counts(alphabet_size, 0);
        long long exhausting_points = 0;

        // The current index in the alphabet we are checking against (0 to k-1)
        int current_alphabet_idx = 0;

        // Iterate through the entire sequence
        // Note: This is a linear pass O(N), very efficient for stream processing.
        for (int val : sequence) {
            
            // Sanity check for input data integrity
            if (val < 0 || val >= alphabet_size) {
                std::cerr << "Error: Sequence value " << val << " out of alphabet range [0, " << alphabet_size - 1 << "]" << std::endl;
                return {0.0, {}, {}, 0, 0};
            }

            // Check for Coincidence: Does sequence value match our current alphabet cycle value?
            if (val == current_alphabet_idx) {
                // MATCH found (Index-Value Coincidence)
                counts[val]++; 
                
                // Per paper Section 2: "After an i-index-value coincidence point, labeling restarts from x1"
                current_alphabet_idx = 0; 
            } else {
                // NO MATCH
                // Move to the next value in the alphabet
                current_alphabet_idx++;

                // Check if we ran out of alphabet options (Exhausting Point)
                if (current_alphabet_idx >= alphabet_size) {
                    // We tried x_1 through x_k and none matched the sequence positions.
                    exhausting_points++;
                    
                    // Per paper: "In that case, the assignment starts over from x1"
                    current_alphabet_idx = 0;
                }
            }
        }

        
        
        // Total terminal points t = sum(t_i) + t_EP
        long long sum_coincidences = std::accumulate(counts.begin(), counts.end(), 0LL);
        long long total_t = sum_coincidences + exhausting_points;

        if (total_t == 0) {
            // Edge case: empty sequence or no termination logic triggered (highly unlikely)
            return {0.0, {}, {}, 0, 0};
        }

        std::vector<double> q_values(alphabet_size);
        
        // We use a simplified iterative calculation for q_i to avoid floating point explosion
        // with the product series in the denominator.
        // Formula derived from paper: 
        // q_1 = t_1 / t
        // q_2 = t_2 / (t * (1-q_1))  -> Denom is essentially "remaining opportunities after t_1 removed"
        // q_i = t_i / (t - sum(t_1..t_{i-1}))
        
        double current_denominator = static_cast<double>(total_t);

        for (int i = 0; i < alphabet_size; ++i) {
            if (current_denominator <= 0) {
                 // Should ideally not happen unless counts exceed total somehow
                q_values[i] = 0.0;
            } else {
                q_values[i] = counts[i] / current_denominator;
            }
            
            // For the next iteration, we remove the counts we just accounted for.
            // This is mathematically equivalent to multiplying by (1-q_i) in the denominator 
            // of the original formula, but numerically more stable.
            current_denominator -= counts[i];
        }


        // p_i = q_i / sum(q_k)
        double sum_q = std::accumulate(q_values.begin(), q_values.end(), 0.0);
        std::vector<double> p_values(alphabet_size);

        // Track max probability for min-entropy
        double max_p = 0.0;

        for (int i = 0; i < alphabet_size; ++i) {
            if (sum_q > 0) {
                p_values[i] = q_values[i] / sum_q;
            } else {
                p_values[i] = 0.0;
            }

            if (p_values[i] > max_p) {
                max_p = p_values[i];
            }
        }

        // --- Min-Entropy ---
        
        // H_inf = -log2(max_p)
        double min_entropy = 0.0;
        if (max_p > 0) {
            min_entropy = -std::log2(max_p);
        }

        return {min_entropy, p_values, q_values, total_t, exhausting_points};
    }
};

// --- Helper function for printing verification tables ---
void print_test_case(std::string name, const std::vector<int>& data, int k) {
    std::cout << "========================================" << std::endl;
    std::cout << "TEST CASE: " << name << std::endl;
    std::cout << "Sequence Length: " << data.size() << std::endl;
    std::cout << "Alphabet Size: " << k << std::endl;
    
    // Run the estimator
    auto result = IndexValueCoincidence::estimate(data, k);

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "RESULTS:" << std::endl;
    std::cout << "Total Terminal Points (t): " << result.total_terminal_points << std::endl;
    std::cout << "Exhausting Points (EP):    " << result.exhausting_points << std::endl;
    
    std::cout << "\nIntermediate Probabilities (q):" << std::endl;
    for(size_t i=0; i<result.q_values.size(); i++) {
        std::cout << "  q_" << i << ": " << std::fixed << std::setprecision(4) << result.q_values[i];
    }
    
    std::cout << "\n\nEstimated Distribution (p):" << std::endl;
    for(size_t i=0; i<result.probabilities.size(); i++) {
        std::cout << "  p_" << i << ": " << std::fixed << std::setprecision(4) << result.probabilities[i];
    }
    
    std::cout << "\n\nMin-Entropy Estimate:" << std::endl;
    std::cout << "  H_inf: " << std::fixed << std::setprecision(4) << result.min_entropy << std::endl;
    std::cout << "========================================" << std::endl << std::endl;
}

int main() {
    // ---------------------------------------------------------
    // Verification against Paper Example 2.1 (Table 3)
    // Sequence: 0 1 1 0 0 1 0 1 0 0 1 0 0 0 1 1
    // Alphabet: {0, 1}
    // Expected Output from Paper:
    //  t = 11
    //  EP = 3
    //  q_0 = 0.5454, q_1 = 0.3999
    //  p_0 = 0.5770, p_1 = 0.4230
    //  H_inf = 0.7935
    // ---------------------------------------------------------

    std::vector<int> paper_example_seq = {0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1};
    print_test_case("Paper Example 2.1 (Binary)", paper_example_seq, 2);


    // ---------------------------------------------------------
    // Additional Test: Random-ish 4-symbol sequence
    // Just to ensure generic k logic holds up without crashing
    // ---------------------------------------------------------
    std::vector<int> symbol_seq = {
        0, 1, 2, 3, 0, 0, 1, 2, // matches expected often
        3, 2, 1, 0, 3, 2, 1, 0  // mismatches often
    };
    print_test_case("Synthetic 4-Symbol Sequence", symbol_seq, 4);

    return 0;
}
