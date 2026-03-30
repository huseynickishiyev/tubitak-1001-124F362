# Fisher-Yates vs Sattolo: A Benchmark on AES-Generated Binary Sequences

## Introduction

A practical comparison of two in-place shuffle algorithms run over 100 binary sequences, each 1,000,000 bits long, generated via OpenSSL's `RAND_bytes` (AES-CTR-DRBG). The sequences are stored as NIST-ready bit-packed `.bin` files (125,000 bytes each). Both per-sequence and total shuffle times are reported separately.

Compiled with `g++ -O2 -std=c++17` on Linux Ubuntu 24.04 x86-64.

---

## The Algorithms

### In-Place Fisher-Yates (Knuth shuffle)

The standard. Described by Fisher and Yates in 1938, formalized for computers by Knuth in 1969, and popularized as the go-to shuffle ever since.

```
for i from n-1 down to 1:
    j = random integer in [0, i]   // inclusive
    swap(A[i], A[j])
```

Every permutation of length `n` has exactly `1/n!` probability of being produced. That is provable by induction: at step `i`, there are `i+1` equally likely choices, so every arrangement of the remaining elements is equally reachable. The result is a uniform draw over all `n!` permutations. Fixed points are allowed. Multi-cycle permutations are allowed. Everything is allowed.

Time complexity: O(n). Space: O(1). Swaps: exactly n-1.

### Sattolo

Published by Sandra Sattolo in 1986. Looks almost identical to Fisher-Yates:

```
for i from n-1 down to 1:
    j = random integer in [0, i-1]  // i is excluded
    swap(A[i], A[j])
```

The only difference is `i-1` instead of `i` as the upper bound. That single change means `j` can never equal `i`, so no element ever stays in its own position during any step. The result is that every output permutation forms exactly one cycle of length `n`. The algorithm samples uniformly over the `(n-1)!` single-cycle permutations, not over all `n!` permutations.

Time complexity: O(n). Space: O(1). Swaps: exactly n-1.

### Why the distinction matters

For `n = 1,000,000`, Fisher-Yates samples from `n!` permutations. Sattolo samples from `(n-1)!`, which is roughly `1/n` of Fisher-Yates's space. Everything Sattolo cannot produce: any permutation with at least one fixed point, any permutation with more than one cycle. That is a lot of excluded ground.

If you are running a NIST SP 800-22 permutation test or any test that assumes a uniform draw over permutation space, Sattolo is the wrong tool. The math is not ambiguous on this.

---

## Experiment Setup

**Data generation.** 100 sequences, each 1,000,000 bits, produced by `RAND_bytes`. That gives 125,000 bytes per sequence in bit-packed form. Each bit is one binary symbol. Bit balance across the first three Fisher-Yates output files: 1,499,950 ones out of 3,000,000 bits, ratio 0.5000.

**Why bit-packed matters.** The previous experiment used one `uint8_t` per symbol (1,000,000 bytes per sequence). NIST's test tools expect one bit per symbol. Feeding byte-packed data to `sts` or `ea_iid` produces wrong results or outright rejection. So the pipeline here is: generate raw bytes via `RAND_bytes`, unpack each bit into a `uint8_t` symbol array, shuffle that array, repack into bit-packed `.bin`, write to disk.

**Timing.** Two separate measurements per algorithm:

- **Per-sequence shuffle time**: only the shuffle step, measured with `high_resolution_clock`. Unpack, repack, and disk I/O are excluded.
- **Total shuffle time (100 sequences)**: the cumulative sum of all 100 per-sequence shuffle times.
- **Total wall time (100 sequences)**: unpack + shuffle + repack + disk write, all included.

5 warmup rounds before any measurement.

---

## Results

### Per-sequence shuffle time (ms)

| Algorithm | Mean | Median | Std Dev | Min | Max | Q3 |
|---|---|---|---|---|---|---|
| Fisher-Yates | 4.66 | 4.31 | 1.31 | 4.05 | 11.44 | 4.72 |
| Sattolo | 4.61 | 4.28 | 1.08 | 4.15 | 11.30 | 4.74 |
| Delta (Sa/FY) | -1.07% | -0.70% | -17.6% | +2.47% | -1.22% | +0.42% |

The medians are 0.03 ms apart. Both algorithms do exactly n-1 swaps and call `uniform_int_distribution` once per iteration. The speed difference is noise.

One thing that stands out: Sattolo's standard deviation is about 17.6% lower than Fisher-Yates's. This is a real pattern across runs, not a one-off. The likely cause is that Sattolo's slightly narrower sampling range reduces the frequency of rejection-sampling edge cases inside the distribution, leading to marginally more predictable iteration times. This has no practical consequence but it is worth noting.

Outliers (above 2x median) appeared at random sequence indices for both algorithms, not clustered at the start or end. That is the signature of OS scheduler preemption, not anything algorithmic.

### Total time over 100 sequences

| | Fisher-Yates | Sattolo | Delta |
|---|---|---|---|
| Shuffle only (ms) | 465.54 | 461.04 | -0.97% |
| Full pipeline including I/O (ms) | 1,206.98 | 1,197.84 | -0.76% |

Total shuffle time lines up with per-sequence averages as expected: `465.54 ms ≈ 100 × 4.66 ms`. The unpack and repack steps add about 741 ms on top of the shuffle time across 100 sequences. That is more expensive than the shuffling itself, which is worth keeping in mind if you are optimizing a pipeline.

### Comparison with the byte-packed experiment

| Format | Working set | FY median (ms) | Sa median (ms) | FY std dev (ms) |
|---|---|---|---|---|
| Byte-packed (1 MB/seq) | ~1.0 MB | 4.93 | 4.91 | 1.51 |
| Bit-packed (125 KB/seq) | ~0.25 MB | 4.31 | 4.28 | 1.31 |

The working set shrinks 8x in the bit-packed setup but total time does not drop proportionally because of the extra conversion steps. The shuffle median is slightly lower, probably because 125 KB fits more comfortably in L3 cache than 1 MB.

---

## Discussion

Speed-wise there is nothing interesting here. Both algorithms are O(n), both do n-1 swaps, both are bottlenecked by the same `mt19937` calls. Any difference in timing is within measurement noise for practical purposes.

The interesting part is what the algorithms actually produce.

Fisher-Yates gives you a uniform random permutation. Every arrangement is equally possible. That is what statistical tests assume when they describe a "random permutation" of a sequence. NIST SP 800-22's permutation-based tests, for example, are built on this assumption.

Sattolo gives you a uniform random single-cycle permutation. That is a strict subset. You never get fixed points. You never get multi-cycle arrangements. For n = 1,000,000, Sattolo is working over approximately 1/1,000,000 of the permutation space that Fisher-Yates covers. If a test compares the original sequence against shuffled versions to assess independence, and the shuffled versions are all structurally biased in the same direction (every element moved, single cycle), the test is no longer measuring what it thinks it is measuring.

This is not a subtle bias you can dismiss. It is a hard structural constraint built into the algorithm.

So: if you need permutation-based IID testing, use Fisher-Yates. If you specifically need single-cycle permutations (derangements, certain combinatorial applications), use Sattolo. Do not swap them.

---

## Output Files

The benchmark writes shuffled sequences to disk in NIST-compatible bit-packed format:

```
nist_sequences/fisher_yates/seq_000.bin  ...  seq_099.bin   # 100 files, 125,000 bytes each
nist_sequences/sattolo/seq_000.bin       ...  seq_099.bin   # 100 files, 125,000 bytes each
```

Fisher-Yates outputs are ready for direct use with `sts` or `ea_iid`. Sattolo outputs can be processed by the same tools but the results should not be interpreted as IID permutation test evidence.

---

## References

1. D. E. Knuth, *The Art of Computer Programming, Vol. 2: Seminumerical Algorithms*, 3rd ed. Addison-Wesley, 1997, pp. 145-146.
2. R. A. Fisher and F. Yates, *Statistical Tables for Biological, Agricultural and Medical Research*. Oliver and Boyd, 1938.
3. S. Sattolo, "An algorithm to generate a random cyclic permutation," *Information Processing Letters*, vol. 22, no. 6, pp. 315-317, 1986.
4. R. Durstenfeld, "Algorithm 235: Random permutation," *Communications of the ACM*, vol. 7, no. 7, p. 420, 1964.
5. OpenSSL Project, "RAND_bytes(3)," 2024. https://www.openssl.org/docs/man3.0/man3/RAND_bytes.html
6. E. Barker and J. Kelsey, "Recommendation for Random Number Generation Using Deterministic Random Bit Generators," NIST SP 800-90A Rev. 1, 2015.
7. A. Rukhin et al., "A Statistical Test Suite for Random and Pseudorandom Number Generators for Cryptographic Applications," NIST SP 800-22 Rev. 1a, 2010.
8. J. Axelsson, "Correctness proofs of Fisher-Yates shuffle," arXiv:1805.00371, 2018.
