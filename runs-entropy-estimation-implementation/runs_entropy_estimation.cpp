/**
 * Entropy Estimation for Random Number Generators
 * Author: Huseyn Kishiyev
 *
 * Implements the statistics proposed in:
 * "New Statistics for Entropy Estimation of Random Number Generators"
 * by Ali Doğanaksoy, Serhat Sağdıçoğlu, and Zülfükar Saygı
 *
 * Implements:
 *   1. Runs statistic (Algorithm 1)
 *   2. t-Runs statistic (Algorithm 2)
 *   3. Min-entropy estimation via runs statistic
 *   4. Peculiar collision statistic (for 4-symbol alphabet)
 */

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <iomanip>
#include <functional>
#include <limits>
#include <stdexcept>

// ============================================================
// Section 3 — Runs Statistic
// ============================================================

/**
 * Algorithm 1: Count the number of runs in a sequence.
 * A run is a maximal contiguous subsequence of identical symbols.
 * Complexity: O(n)
 */
template <typename T>
int countRuns(const std::vector<T>& sigma) {
    int n = static_cast<int>(sigma.size());
    if (n == 0) return 0;
    int runs = 1;
    for (int i = 0; i < n - 1; ++i) {
        if (sigma[i] != sigma[i + 1])
            ++runs;
    }
    return runs;
}

/**
 * Algorithm 2: Count the number of t-runs for all t in a sequence.
 * R[t] gives the count of runs of length exactly t (1-indexed).
 * Complexity: O(n)
 *
 * Returns a vector of size n+1 (R[1]..R[n]; R[0] unused).
 */
template <typename T>
std::vector<int> countTRuns(const std::vector<T>& sigma) {
    int n = static_cast<int>(sigma.size());
    std::vector<int> R(n + 1, 0);
    if (n == 0) return R;
    if (n == 1) { R[1] = 1; return R; }
    int j = 1; // current run length accumulator
    for (int i = 0; i < n - 1; ++i) {
        if (sigma[i] != sigma[i + 1]) {
            R[j] += 1;
            j = 1;
        } else {
            j += 1;
        }
    }
    R[j] += 1; // close the final run
    return R;
}

// ============================================================
// Proposition 1 — Expected number of runs
// E_p[R] = (n-1)(1 - phi) + 1,  phi = sum(p_i^2)
// ============================================================
double expectedRuns(int n, const std::vector<double>& probs) {
    double phi = 0.0;
    for (double p : probs) phi += p * p;
    return (n - 1) * (1.0 - phi) + 1.0;
}

// ============================================================
// Proposition 2 — Variance of the number of runs
// Var_p[R] = (n-1)(1-phi)*phi
// ============================================================
double varianceRuns(int n, const std::vector<double>& probs) {
    double phi = 0.0;
    for (double p : probs) phi += p * p;
    return (n - 1) * (1.0 - phi) * phi;
}

// ============================================================
// Proposition 3 — Expected number of t-runs (binary case, Remark 2)
// E[R_t] = (n-t-1)(p^2*q^t + q^2*p^t) + 2*(p^t*q + q^t*p)
// ============================================================
double expectedTRunsBinary(int n, int t, double p) {
    double q = 1.0 - p;
    double pt = std::pow(p, t);
    double qt = std::pow(q, t);
    return (n - t - 1) * (p * p * qt + q * q * pt) + 2.0 * (pt * q + qt * p);
}

// ============================================================
// Section 3.1 — Min-entropy estimation via Runs Statistic
// Binary case: E[R] = 2(n-1)*p*(1-p) + 1
// Solve for p given observed R, then H_inf = -log2(p_max)
// ============================================================
double binaryRunsMinEntropy(int n, int R_obs) {
    // The maximum of E[R] for binary sequences is (n+1)/2 at p=0.5.
    // If the observed R exceeds this (possible in finite samples),
    // the best estimate is p=0.5, giving max entropy of 1 bit.
    double max_expected = (n + 1.0) / 2.0;
    if (static_cast<double>(R_obs) >= max_expected) return 1.0;

    // Solve: 2(n-1)*p^2 - 2(n-1)*p + (R_obs - 1) = 0
    double a = 2.0 * (n - 1);
    double b = -2.0 * (n - 1);
    double c = static_cast<double>(R_obs - 1);
    double disc = b * b - 4.0 * a * c;
    if (disc < 0.0) return 1.0; // numerical safety — clamp to max entropy

    double sqrtDisc = std::sqrt(disc);
    double p1 = (-b + sqrtDisc) / (2.0 * a);
    double p2 = (-b - sqrtDisc) / (2.0 * a);

    // Step 4: p = min(max(p1,p2), 1)
    double p = std::min(std::max(p1, p2), 1.0);
    if (p <= 0.0) return std::numeric_limits<double>::infinity();
    return -std::log2(p);
}

/**
 * General (k-symbol) runs min-entropy estimation using near-uniform family.
 * phi(theta) = theta^2 + (1-theta)^2/(k-1)
 * E[R] = (n-1)*(1-phi(theta)) + 1  => solve for theta in [1/k, 1]
 * Uses bisection (E[R] is monotone decreasing in theta on [1/k,1]).
 */
double generalRunsMinEntropy(int n, int k, int R_obs) {
    auto phi = [&](double theta) -> double {
        double other = (1.0 - theta) / (k - 1);
        return theta * theta + (k - 1) * other * other;
    };
    // f(theta) = E[R](theta) - R_obs; f is decreasing in theta
    auto f = [&](double theta) -> double {
        return (n - 1) * (1.0 - phi(theta)) + 1.0 - R_obs;
    };

    double lo = 1.0 / k, hi = 1.0;
    double flo = f(lo), fhi = f(hi);

    // If R_obs is above or below the achievable range, clamp
    if (flo < 0) return -std::log2(1.0 / k); // min entropy for uniform
    if (fhi > 0) return 0.0;                  // max bias

    for (int iter = 0; iter < 300; ++iter) {
        double mid = (lo + hi) / 2.0;
        if (f(mid) > 0) lo = mid;
        else            hi = mid;
    }
    double theta = (lo + hi) / 2.0;
    return -std::log2(theta);
}

// ============================================================
// Section 4 — Collision Statistic (Definition 4)
// ============================================================

template <typename T>
struct CollisionResult {
    std::vector<int> times;    // T_1, T_2, ..., T_m (1-indexed positions in sigma)
    std::vector<T>   symbols;  // symbol at each collision
};

/**
 * Compute collision times per Definition 4.
 * T_j = min{l > T_{j-1} : exists r in (T_{j-1}, l) with sigma[l] == sigma[r]}
 */
template <typename T>
CollisionResult<T> computeCollisions(const std::vector<T>& sigma) {
    int n = static_cast<int>(sigma.size());
    CollisionResult<T> result;
    int prev = 0; // T_{j-1} in 1-indexed (T_0 = 0)

    while (true) {
        bool found = false;
        for (int l = prev; l < n && !found; ++l) {
            for (int r = prev; r < l; ++r) {
                if (sigma[r] == sigma[l]) {
                    result.times.push_back(l + 1); // 1-indexed
                    result.symbols.push_back(sigma[l]);
                    prev = l + 1;  // prev stores T_{j-1} in 1-indexed
                    found = true;
                    break;
                }
            }
        }
        if (!found) break;
    }
    return result;
}

/**
 * Standard collision statistic: S_m = T_m / m
 */
template <typename T>
double collisionStatistic(const CollisionResult<T>& cr) {
    if (cr.times.empty()) return 0.0;
    int m = static_cast<int>(cr.times.size());
    return static_cast<double>(cr.times.back()) / m;
}

/**
 * Approximate min-entropy from collision statistic.
 * Per NIST SP 800-90B: H ~= log2(S_m - 1)
 */
double collisionMinEntropy(double S) {
    if (S <= 1.0) return 0.0;
    return std::log2(S - 1.0);
}

// ============================================================
// Section 4.1 — Peculiar Collision Statistic (4-symbol alphabet)
// ============================================================

/**
 * Count collisions per symbol.
 */
template <typename T>
std::map<T, int> peculiarCollisionCounts(const CollisionResult<T>& cr) {
    std::map<T, int> counts;
    for (const T& sym : cr.symbols) counts[sym]++;
    return counts;
}

/**
 * pc(0) under near-uniform distribution, 4-symbol alphabet (Equation 7):
 *   pc(0) = theta^2/9 * (53 - 78*theta + 42*theta^2 - 8*theta^3)
 * where theta = pg(0) in [1/4, 1].
 */
double pcNearUniform4(double theta) {
    return (theta * theta / 9.0)
        * (53.0 - 78.0*theta + 42.0*theta*theta - 8.0*theta*theta*theta);
}

/**
 * pc(0) under inverted near-uniform distribution, 4-symbol alphabet.
 * Three sub-cases (Equations 8, 9, 10):
 *   [1/4, 1/3): pg(0)=pg(1)=pg(2)=theta, pg(3)=1-3*theta  (Eq. 8)
 *   [1/3, 1/2): pg(0)=pg(1)=theta, pg(2)=1-2*theta          (Eq. 9)
 *   [1/2, 1]:   pg(0)=theta, pg(1)=1-theta                   (Eq. 10)
 */
double pcInvertedNearUniform4(double theta) {
    if (theta < 0.25) return std::numeric_limits<double>::quiet_NaN();
    if (theta < 1.0 / 3.0)
        return theta*theta * (3.0 + 10.0*theta - 6.0*theta*theta - 72.0*theta*theta*theta);
    else if (theta < 0.5)
        return theta*theta * (3.0 + 4.0*theta - 12.0*theta*theta);
    else
        return theta*theta * (3.0 - 2.0*theta);
}

/**
 * Bisection solver: find theta in [lo,hi] such that f(theta) == target.
 * f must be monotone on the interval.
 * Returns NaN if no root found.
 */
double bisectSolve(std::function<double(double)> f,
                   double target, double lo, double hi, int iters = 300)
{
    double flo = f(lo) - target;
    double fhi = f(hi) - target;
    if (flo * fhi > 0) return std::numeric_limits<double>::quiet_NaN();
    for (int i = 0; i < iters; ++i) {
        double mid = (lo + hi) / 2.0;
        double fmid = f(mid) - target;
        if (flo * fmid <= 0) { hi = mid; fhi = fmid; }
        else                  { lo = mid; flo = fmid; }
    }
    return (lo + hi) / 2.0;
}

/**
 * Peculiar collision entropy bounds for 4-symbol alphabet.
 *
 * Given pc_obs = nc(max_symbol) / C (collision rate for most-colliding symbol),
 * and C = total number of collisions:
 *   - Near-uniform  => lower bound on pg(0)  (pl)
 *   - Inv. near-uniform => upper bound on pg(0)  (pu)
 *   - 99% CI: z = 2.576
 *
 * Returns bounds: pl <= pg(0) <= pu, and entropy bounds -log2(pu) <= H <= -log2(pl).
 */
struct PeculiarCollisionResult {
    double pc_obs;
    double pg0_hat_nu;        // near-uniform estimate (lower bound estimate)
    double pg0_hat_inu;       // inv. near-uniform estimate (upper bound estimate)
    double pg0_lower;         // pl = pg0_hat_nu  - CI
    double pg0_upper;         // pu = pg0_hat_inu + CI
    double entropy_lower;     // -log2(pu)
    double entropy_upper;     // -log2(pl)
};

PeculiarCollisionResult peculiarCollisionEntropy4(double pc_obs, int C) {
    // --- Near-uniform: solve Eq. 7 for theta in [0.25, 1] ---
    // pcNearUniform4 is increasing in theta on [0.25, 1]
    double theta_nu = bisectSolve(pcNearUniform4, pc_obs, 0.25, 1.0);
    if (std::isnan(theta_nu)) theta_nu = 0.25; // fallback to uniform

    // --- Inverted near-uniform: try all three sub-ranges ---
    double theta_inu = std::numeric_limits<double>::quiet_NaN();
    double t10 = bisectSolve(pcInvertedNearUniform4, pc_obs, 0.5,       1.0      );
    double t9  = bisectSolve(pcInvertedNearUniform4, pc_obs, 1.0/3.0,   0.5      );
    double t8  = bisectSolve(pcInvertedNearUniform4, pc_obs, 0.25,      1.0/3.0  );

    // Pick the solution consistent with its sub-range
    if (!std::isnan(t10) && t10 >= 0.5 && t10 <= 1.0)
        theta_inu = t10;
    else if (!std::isnan(t9) && t9 >= 1.0/3.0 && t9 < 0.5)
        theta_inu = t9;
    else if (!std::isnan(t8) && t8 >= 0.25 && t8 < 1.0/3.0)
        theta_inu = t8;

    if (std::isnan(theta_inu)) theta_inu = 1.0; // fallback

    const double z = 2.576; // 99% confidence

    double ci_nu  = z * std::sqrt(theta_nu  * (1.0 - theta_nu)  / C);
    double ci_inu = z * std::sqrt(theta_inu * (1.0 - theta_inu) / C);

    double pl = std::max(theta_nu  - ci_nu,  1e-12);
    double pu = std::min(theta_inu + ci_inu, 1.0  );

    PeculiarCollisionResult res;
    res.pc_obs        = pc_obs;
    res.pg0_hat_nu    = theta_nu;
    res.pg0_hat_inu   = theta_inu;
    res.pg0_lower     = pl;
    res.pg0_upper     = pu;
    res.entropy_lower = -std::log2(pu); // -log2(pu) <= H
    res.entropy_upper = -std::log2(pl); // H <= -log2(pl)
    return res;
}

// ============================================================
// Utility
// ============================================================
void printSep(const std::string& title) {
    std::cout << "\n" << std::string(62, '=') << "\n"
              << "  " << title << "\n"
              << std::string(62, '=') << "\n";
}

template <typename T>
void printSeq(const std::vector<T>& seq) {
    for (const auto& v : seq) std::cout << v << " ";
    std::cout << "\n";
}

// ============================================================
// Main: reproduce all examples from the paper
// ============================================================
int main() {
    std::cout << std::fixed << std::setprecision(6);

    // ----------------------------------------------------------
    // Algorithm 1 — countRuns
    // ----------------------------------------------------------
    printSep("Algorithm 1: Runs Statistic (Eq. 1)");

    std::vector<int> seq_ex3 = {2,2,0,1,0,2,0,1,2,1,2,0,1,2,1,0,0,1,0,0,0};
    std::cout << "Sequence (paper Example 3):\n  ";
    printSeq(seq_ex3);
    std::cout << "Number of runs: " << countRuns(seq_ex3) << "\n";

    // ----------------------------------------------------------
    // Algorithm 2 — countTRuns
    // ----------------------------------------------------------
    printSep("Algorithm 2: t-Runs Statistic");

    auto tR = countTRuns(seq_ex3);
    std::cout << "t-run counts for Example 3 sequence:\n";
    for (int t = 1; t <= (int)seq_ex3.size(); ++t)
        if (tR[t] > 0) std::cout << "  t=" << t << ": " << tR[t] << "\n";

    // ----------------------------------------------------------
    // Proposition 1 — Expected runs + Example 1
    // ----------------------------------------------------------
    printSep("Proposition 1: Expected Runs & Min-Entropy Estimation");

    {
        int n = 100, R = 42;
        std::cout << "Example 1 (binary, n=" << n << ", R=" << R << "):\n";
        // Solve quadratic manually to show both roots
        double a = 2.0*(n-1), b = -2.0*(n-1), c = R - 1.0;
        double disc = b*b - 4*a*c;
        double p1 = (-b + std::sqrt(disc))/(2*a);
        double p2 = (-b - std::sqrt(disc))/(2*a);
        std::cout << "  Roots: p1=" << p1 << "  p2=" << p2
                  << "  (paper: 0.71, 0.29)\n";
        std::cout << "  Min-entropy = -log2(" << p1 << ") = "
                  << binaryRunsMinEntropy(n, R) << " bits\n";
    }

    // ----------------------------------------------------------
    // Example 2 — Quaternary
    // ----------------------------------------------------------
    {
        int n = 21, R = 13, k = 4;
        std::cout << "\nExample 2 (quaternary k=4, n=" << n << ", R=" << R << "):\n";
        double H = generalRunsMinEntropy(n, k, R);
        double theta = std::pow(2.0, -H);
        double other = (1.0 - theta) / (k - 1);
        std::cout << "  Min-entropy = " << H << " bits\n";
        std::cout << "  p1=" << theta << "  p2=p3=p4=" << other
                  << "  (paper: 0.5854, 0.1382)\n";
    }

    // ----------------------------------------------------------
    // Proposition 2 — Variance
    // ----------------------------------------------------------
    printSep("Proposition 2: Variance of Runs");

    {
        std::vector<double> p_half = {0.5, 0.5};
        int n = 100;
        std::cout << "Binary p=0.5, n=100:\n";
        std::cout << "  E[R]   = " << expectedRuns(n, p_half)
                  << "  (expected " << (n+1)/2.0 << ")\n";
        std::cout << "  Var[R] = " << varianceRuns(n, p_half) << "\n";
    }

    // ----------------------------------------------------------
    // Proposition 3 — Expected t-runs
    // ----------------------------------------------------------
    printSep("Proposition 3: Expected t-Runs (binary, Remark 2)");

    {
        int n = 100; double p = 0.5;
        std::cout << "Binary p=0.5, n=100:\n";
        for (int t = 1; t <= 5; ++t)
            std::cout << "  E[R_" << t << "] = " << expectedTRunsBinary(n, t, p) << "\n";
        std::cout << "  (Remark 2 formula: (n-t+3)/2^{t+1})\n";
        for (int t = 1; t <= 5; ++t)
            std::cout << "  Remark2 check t=" << t << ": " << (n-t+3.0)/std::pow(2,t+1) << "\n";
    }

    // ----------------------------------------------------------
    // Section 4 — Standard Collision Statistic (Example 3)
    // ----------------------------------------------------------
    printSep("Section 4: Collision Statistic (Example 3)");

    {
        // Paper uses 3-symbol alphabet {0,1,2} here
        auto cr = computeCollisions(seq_ex3);
        std::cout << "Sequence: "; printSeq(seq_ex3);
        std::cout << "Collision times (1-indexed): "; printSeq(cr.times);
        std::cout << "Collision symbols:            "; printSeq(cr.symbols);
        int m = static_cast<int>(cr.times.size());
        double S = collisionStatistic(cr);
        std::cout << "m=" << m << "  S_m=" << S
                  << "  H_inf~=" << collisionMinEntropy(S)
                  << "  (paper: 6 collisions, H=0.5016)\n";
    }

    // ----------------------------------------------------------
    // Section 4.1 — Peculiar Collision Statistic (Example 4)
    // ----------------------------------------------------------
    printSep("Section 4.1: Peculiar Collision Statistic (Example 4)");

    {
        // 4-symbol alphabet {0,1,2,3}
        std::vector<int> s1 = {0,0,0,2,0,0,2,1,0,0,1,3,0,0,3,2,0,1,0,1,3};
        std::vector<int> s2 = {2,2,0,2,0,0,2,1,0,0,1,3,0,0,3,2,0,1,0,1,3};

        auto cr1 = computeCollisions(s1);
        auto cr2 = computeCollisions(s2);

        std::cout << "sigma1: "; printSeq(s1);
        std::cout << "sigma2: "; printSeq(s2);

        auto cnt1 = peculiarCollisionCounts(cr1);
        auto cnt2 = peculiarCollisionCounts(cr2);
        int C1 = static_cast<int>(cr1.times.size());
        int C2 = static_cast<int>(cr2.times.size());

        std::cout << "\nCollision times sigma1: "; printSeq(cr1.times);
        std::cout << "Collision times sigma2: "; printSeq(cr2.times);
        std::cout << "(Both have same positions per the paper)\n";

        std::cout << "\nsigma1 peculiar collision counts (total C=" << C1 << "):\n";
        for (auto& [sym, c] : cnt1) std::cout << "  symbol " << sym << ": " << c << "\n";

        std::cout << "\nsigma2 peculiar collision counts (total C=" << C2 << "):\n";
        for (auto& [sym, c] : cnt2) std::cout << "  symbol " << sym << ": " << c << "\n";

        // Max-collision symbol is 0 in both; pc_obs = nc(0)/C
        double pc1 = static_cast<double>(cnt1.count(0) ? cnt1.at(0) : 0) / C1;
        double pc2 = static_cast<double>(cnt2.count(0) ? cnt2.at(0) : 0) / C2;

        std::cout << "\npc_obs(0) sigma1 = " << cnt1[0] << "/" << C1
                  << " = " << pc1 << "  (paper: 5/6)\n";
        std::cout << "pc_obs(0) sigma2 = " << cnt2[0] << "/" << C2
                  << " = " << pc2 << "  (paper: 4/6)\n";

        auto r1 = peculiarCollisionEntropy4(pc1, C1);
        auto r2 = peculiarCollisionEntropy4(pc2, C2);

        std::cout << "\n--- Entropy bounds (99% CI) ---\n";
        std::cout << "\nsigma1:\n"
                  << "  pg(0) in [" << r1.pg0_lower << ", " << r1.pg0_upper << "]\n"
                  << "  H_inf in [" << r1.entropy_lower << ", " << r1.entropy_upper << "] bits\n";

        std::cout << "\nsigma2:\n"
                  << "  pg(0) in [" << r2.pg0_lower << ", " << r2.pg0_upper << "]\n"
                  << "  H_inf in [" << r2.entropy_lower << ", " << r2.entropy_upper << "] bits\n";

        std::cout << "\nConclusion: Both sequences have identical collision positions,\n"
                  << "so the standard collision statistic gives the SAME min-entropy estimate.\n"
                  << "The peculiar collision statistic yields DIFFERENT bounds for sigma1 vs sigma2,\n"
                  << "demonstrating its greater sensitivity (as claimed in the paper).\n";

        // Verify with standard collision (should be identical for both)
        double S1 = collisionStatistic(cr1);
        double S2 = collisionStatistic(cr2);
        std::cout << "\nStandard collision S_m: sigma1=" << S1
                  << "  sigma2=" << S2 << "  (identical, as expected)\n";
    }

    return 0;
}
