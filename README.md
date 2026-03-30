# Entropy Estimation for RNGs — Research Repository

**Authors:** Huseyn Kishiyev, Osman Furkan Atakul, Özge Baran Yelim

This repository contains C++ implementations of novel entropy estimators for random number generators (RNGs), a full suite of synthetic dataset generators for NIST SP 800-90B evaluation, and Python tooling for post-analysis. The work accompanies and validates the algorithms described in the following papers:

- Doğanaksoy, Sağdıçoğlu, Saygı — *"New Statistics for Entropy Estimation of Random Number Generators"*
- Aslan et al. — *"Estimating Entropy via Index-Value Coincidence"*
- Aslan et al. — *"Observations on NIST Entropy Suite"*

Disclaimer: The "Authors" named throughout this repository (Huseyn Kishiyev, Osman Furkan Atakul, Özge Baran Yelim) refer solely to the individuals who implemented the software and built the evaluation pipeline. The underlying mathematical methods, estimators, and statistical frameworks are the original work of the respective researchers cited in each source file and in the references above. No claim is made over any of the original intellectual contributions. This repository exists purely as a research implementation effort. The authors declare no conflict of interest.

---

## Table of Contents

1. [Repository Structure](#repository-structure)
2. [Dependencies](#dependencies)
   - [Ubuntu 24.04 LTS](#ubuntu-2404-lts)
   - [macOS (Apple Silicon)](#macos-apple-silicon)
3. [Source Files](#source-files)
   - [runs\_entropy\_estimation.cpp](#runs_entropy_estimationcpp)
   - [index\_value\_coincidence.cpp](#index_value_coincidencecpp)
   - [generate\_data.cpp](#generate_datacpp)
   - [aes\_full.cpp](#aes_fullcpp)
   - [gen\_biased\_binary.cpp](#gen_biased_binarycpp)
   - [gen\_4bit.cpp](#gen_4bitcpp)
   - [gen\_8bit.cpp](#gen_8bitcopp)
   - [gen\_noniid.cpp](#gen_noniidcpp)
4. [Pipeline Scripts](#pipeline-scripts)
   - [run\_nist.sh](#run_nistsh)
   - [analyze\_results.py](#analyze_resultspy)
5. [End-to-End Workflow](#end-to-end-workflow)
6. [NIST SP 800-90B Tool Reference](#nist-sp-800-90b-tool-reference)
7. [Dataset Summary](#dataset-summary)

---

## Repository Structure

```
.
├── data-generation/
│   ├── aes_full.cpp
│   ├── generate_data.cpp
│   ├── gen_4bit.cpp
│   ├── gen_8bit.cpp
│   ├── gen_biased_binary.cpp
│   ├── gen_noniid.cpp
│   ├── data/                          # AES-256-CBC sequences (sequence_0.bin … sequence_199.bin)
│   ├── output_aes_cbc_unpacked/       # AES-128-CBC sequences (sequence_0000.bin … sequence_0199.bin)
│   ├── output_4bit_unpacked/
│   ├── output_8bit/
│   ├── output_biased_binary_unpacked/
│   └── output_noniid_unpacked/
├── index-value-coincidence/
│   └── index_value_coincidence.cpp
├── runs-entropy-estimation-implementation/
│   └── runs_entropy_estimation.cpp
├── shuffle-algorithms-comparison/
│   ├── shuffle_algorithms_comparison.cpp
│   └── explanation.md
├── analyze_results.py
├── run_nist.sh
├── results.txt                        # NIST tool output (created by run_nist.sh)
└── LICENSE
```

---

## Dependencies

### Ubuntu 24.04 LTS

**System packages:**

```bash
sudo apt update
sudo apt install -y \
    build-essential \
    g++ \
    libssl-dev \
    python3 \
    python3-pip \
    python3-venv
```

**Python packages:**

```bash
pip3 install pandas scipy
```

**NIST SP 800-90B Entropy Assessment Tool:**

```bash
git clone https://github.com/usnistgov/SP800-90B_EntropyAssessment.git nist_entropy
cd nist_entropy/cpp
make non_iid
```

> The `ea_non_iid` binary will be built at `nist_entropy/cpp/ea_non_iid`. `run_nist.sh` automatically resolves this path.

---

### macOS (Apple Silicon)

**Homebrew and Xcode Command Line Tools:**

```bash
xcode-select --install
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

**System packages:**

```bash
brew install openssl@3 python
```

Because macOS ships a stub `openssl` that points to LibreSSL, you must explicitly pass the Homebrew OpenSSL path when compiling any file that uses `<openssl/evp.h>`:

```bash
export OPENSSL_PREFIX=$(brew --prefix openssl@3)
# e.g., /opt/homebrew/opt/openssl@3
```

You will use `-I$OPENSSL_PREFIX/include` and `-L$OPENSSL_PREFIX/lib` in every compile command that links OpenSSL (see per-file instructions below).

**Python packages:**

```bash
pip3 install pandas scipy
```

**NIST SP 800-90B Entropy Assessment Tool:**

```bash
git clone https://github.com/usnistgov/SP800-90B_EntropyAssessment.git nist_entropy
cd nist_entropy/cpp
# Apple Silicon may need these flags:
CXXFLAGS="-I$(brew --prefix openssl@3)/include" \
LDFLAGS="-L$(brew --prefix openssl@3)/lib" \
make non_iid
```

---

## Source Files

---

### `runs_entropy_estimation.cpp`

**Purpose:** Implements four entropy estimation methods from Doğanaksoy et al.:

| Method | Description |
|---|---|
| `countRuns` | Algorithm 1 — counts total runs in a sequence |
| `countTRuns` | Algorithm 2 — counts runs of exact length *t* |
| `binaryRunsMinEntropy` | Solves the quadratic for binary min-entropy from observed run count |
| `generalRunsMinEntropy` | Bisection solver for *k*-symbol near-uniform min-entropy |
| `computeCollisions` | Computes collision times per Definition 4 |
| `collisionStatistic` | Standard collision statistic *S_m* |
| `peculiarCollisionEntropy4` | Peculiar collision estimator for 4-symbol alphabet with 99% CI bounds |

The `main()` function reproduces all paper examples (Examples 1–4, Propositions 1–3) with expected values annotated in comments for verification.

**No external dependencies** — uses only the C++ standard library.

**Compile:**

```bash
# Ubuntu / macOS
g++ -O2 -std=c++17 -o runs_entropy runs_entropy_estimation.cpp
```

**Run:**

```bash
./runs_entropy
```

**Expected output:** Reproduced paper examples with run counts, expected/variance statistics, min-entropy estimates, and peculiar collision bounds.

---

### `index_value_coincidence.cpp`

**Purpose:** Implements the Index-Value Coincidence (IVC) estimator from Aslan et al. The algorithm performs a single O(N) pass over a sequence, counting *coincidence points* (where a sequence value matches the current cyclic alphabet index) and *exhausting points* (where no match was found across the full alphabet). Min-entropy is derived from the resulting probability distribution.

**Key formulas:**

```
q_i = counts[i] / running_denominator
p_i = q_i / sum(q)
H_inf = -log2(max(p_i))
```

The `main()` function verifies against **paper Example 2.1 (Table 3)**:

| Quantity | Paper | Implementation |
|---|---|---|
| t (terminal points) | 11 | 11 |
| EP (exhausting points) | 3 | 3 |
| q_0 | 0.5454 | 0.5454 |
| p_0 | 0.5770 | 0.5770 |
| H_inf | 0.7935 | 0.7935 |

**No external dependencies.**

**Compile:**

```bash
# Ubuntu / macOS
g++ -O2 -std=c++17 -o ivc index_value_coincidence.cpp
```

**Run:**

```bash
./ivc
```

---

### `generate_data.cpp`

**Purpose:** Primary dataset generator. Produces **200 full-entropy IID binary sequences** of 1,000,000 bits each, stored as 1 byte per bit (`0x00` / `0x01`). Uses **AES-256-CBC** with a random key and IV per sequence, encrypting an all-zero plaintext block and expanding the ciphertext bits. This is the data generation methodology described in Aslan et al.

**Output directory:** `data/`  
**Output format:** `data/sequence_N.bin` (N = 0 … 199), each file is exactly **1,000,000 bytes**.

**Dependencies:** OpenSSL (`libssl`, `libcrypto`)

**Compile:**

```bash
# Ubuntu
g++ -O2 -std=c++17 -o generate_data generate_data.cpp -lssl -lcrypto

# macOS (Apple Silicon)
OPENSSL_PREFIX=$(brew --prefix openssl@3)
g++ -O2 -std=c++17 \
    -I$OPENSSL_PREFIX/include \
    -L$OPENSSL_PREFIX/lib \
    -o generate_data generate_data.cpp \
    -lssl -lcrypto
```

**Run:**

```bash
mkdir -p data
./generate_data
```

---

### `aes_full.cpp`

**Purpose:** Alternative AES dataset generator using **AES-128-CBC** via OpenSSL's EVP API. Functionally equivalent to `generate_data.cpp` but uses a 128-bit key instead of 256-bit. Produces 200 sequences of 1,000,000 binary samples each.

**Output directory:** `output_aes_cbc_unpacked/`  
**Output format:** `output_aes_cbc_unpacked/sequence_XXXX.bin` (zero-padded 4-digit index)

**Dependencies:** OpenSSL (`libssl`, `libcrypto`)

**Compile:**

```bash
# Ubuntu
g++ -O2 -o gen_aes_cbc_unpacked aes_full.cpp -lssl -lcrypto

# macOS (Apple Silicon)
OPENSSL_PREFIX=$(brew --prefix openssl@3)
g++ -O2 \
    -I$OPENSSL_PREFIX/include \
    -L$OPENSSL_PREFIX/lib \
    -o gen_aes_cbc_unpacked aes_full.cpp \
    -lssl -lcrypto
```

**Run:**

```bash
mkdir -p output_aes_cbc_unpacked
./gen_aes_cbc_unpacked
```

**NIST invocation:**

```bash
./nist_entropy/cpp/ea_iid -i -a output_aes_cbc_unpacked/sequence_0000.bin 1
```

---

### `gen_biased_binary.cpp`

**Purpose:** Dataset 2 — generates **biased binary sequences** with P(1) = 0.7, P(0) = 0.3. Intended to test how NIST estimators respond to known, moderate bias. Theoretical min-entropy: −log₂(0.7) ≈ **0.515 bits**.

**Output directory:** `output_biased_binary_unpacked/`  
**Samples:** 200 sequences × 1,000,000 samples, 1 byte per bit  
**Dependencies:** C++ standard library only (`<random>`)

**Compile:**

```bash
g++ -O2 -std=c++17 -o gen_biased_binary gen_biased_binary.cpp
```

**Run:**

```bash
mkdir -p output_biased_binary_unpacked
./gen_biased_binary
```

**NIST invocation:**

```bash
./nist_entropy/cpp/ea_iid -i -a output_biased_binary_unpacked/sequence_0000.bin 1
```

---

### `gen_4bit.cpp`

**Purpose:** Dataset 3 — generates **4-bit near-uniform biased sequences** over a 16-symbol alphabet (0x0–0xF). Symbol 0 has P = 0.25; all other 15 symbols share the remaining probability equally (P = 0.05 each). Each sample is stored as one byte in `[0, 15]`.

**Output directory:** `output_4bit_unpacked/`  
**Samples:** 200 sequences × 250,000 samples  
**Dependencies:** C++ standard library only (`<random>`)

**Compile:**

```bash
g++ -O2 -std=c++17 -o gen_4bit gen_4bit.cpp
```

**Run:**

```bash
mkdir -p output_4bit_unpacked
./gen_4bit
```

**NIST invocation** (pass bits-per-symbol = 4):

```bash
./nist_entropy/cpp/ea_iid -i -a output_4bit_unpacked/sequence_0000.bin 4
```

---

### `gen_8bit.cpp`

**Purpose:** Dataset 4 — generates **8-bit near-uniform biased sequences** over a 256-symbol alphabet. Symbol 0x00 has P = 0.06; all other 255 symbols share 0.94 equally (~0.00369 each). Uses binary search over the CDF for efficient sampling.

**Output directory:** `output_8bit/`  
**Samples:** 200 sequences × 125,000 samples  
**Dependencies:** C++ standard library only (`<random>`)

**Compile:**

```bash
g++ -O2 -std=c++17 -o gen_8bit gen_8bit.cpp
```

**Run:**

```bash
mkdir -p output_8bit
./gen_8bit
```

**NIST invocation:**

```bash
./nist_entropy/cpp/ea_iid -i -a output_8bit/sequence_0000.bin 8
```

---

### `gen_noniid.cpp`

**Purpose:** Dataset 5 — generates **non-IID binary sequences** with a deterministic parity dependency. Every 8th bit (positions 7, 15, 23, …, 0-indexed) is set to the XOR (parity) of the preceding 7 independently generated bits. All other bits are i.i.d. uniform. This introduces a structured, detectable dependence without changing the marginal distribution.

**Generation rule:**

```
s_{8k} = s_{8k-7} XOR s_{8k-6} XOR ... XOR s_{8k-1}
```

**Output directory:** `output_noniid_unpacked/`  
**Samples:** 200 sequences × 1,000,000 samples, 1 byte per bit  
**Dependencies:** C++ standard library only (`<random>`)

**Compile:**

```bash
g++ -O2 -std=c++17 -o gen_noniid gen_noniid.cpp
```

**Run:**

```bash
mkdir -p output_noniid_unpacked
./gen_noniid
```

**NIST invocation** — use the non-IID assessment (this dataset is intentionally not IID):

```bash
./nist_entropy/cpp/ea_non_iid output_noniid_unpacked/sequence_0000.bin 1
```

To observe IID test *failure* behavior, you may also run it through the IID suite:

```bash
./nist_entropy/cpp/ea_iid -i -a output_noniid_unpacked/sequence_0000.bin 1
```

---

## Pipeline Scripts

### `run_nist.sh`

**Purpose:** Automates the full NIST SP 800-90B non-IID assessment across all sequences in the `data/` directory. Iterates over every `sequence_*.bin` file, calls `ea_non_iid -v -a` on each, and appends all output to `results.txt` with `--- START ... ---` / `--- END ... ---` delimiters for downstream parsing.

**Auto-build:** If `ea_non_iid` is missing and a `Makefile` is present in `nist_entropy/cpp/`, the script will attempt to build the tool automatically.

**Usage:**

```bash
# First generate data (example using generate_data.cpp output):
mkdir -p data && ./generate_data

# Then run the full NIST pipeline:
chmod +x run_nist.sh
./run_nist.sh
```

**Options:** Override the tool path at runtime:

```bash
NIST_TOOL=/path/to/custom/ea_non_iid ./run_nist.sh
```

**Output:** `results.txt` — concatenation of all `ea_non_iid` outputs, one block per sequence.

---

### `analyze_results.py`

**Purpose:** Parses the structured output in `results.txt` and produces a correlation analysis across the 10 NIST Non-IID estimators (E1–E10). Estimator labels map to:

| ID | NIST Estimator Name |
|---|---|
| E1 | Most Common Value Estimate (bit string) |
| E2 | Collision Test Estimate (bit string) |
| E3 | Markov Test Estimate (bit string) |
| E4 | Compression Test Estimate (bit string) |
| E5 | T-Tuple Test Estimate (bit string) |
| E6 | LRS Test Estimate (bit string) |
| E7 | Multi Most Common in Window (MultiMCW) |
| E8 | Lag Prediction Test Estimate (bit string) |
| E9 | Multi Markov Model with Counting (MultiMMC) |
| E10 | LZ78Y Prediction Test Estimate |

The script computes and prints both **Pearson** and **Spearman** correlation matrices. Correlations above the threshold (default 0.4) are **bold-highlighted** in the terminal output.

**Dependencies:** `pandas`, `scipy`

**Run:**

```bash
python3 analyze_results.py
```

> `results.txt` must be present in the same directory. Run `run_nist.sh` first.

**Example output:**

```
Successfully parsed results for 200 sequences.

## Pearson Correlation Matrix
--------------------------------------------------------------
         E1    E2    E3    E4    E5    E6    E7    E8    E9   E10
E1    1.000 0.312 0.284 ...
...
```

---

## End-to-End Workflow

Below is the recommended full pipeline from dataset generation to analysis.

**Step 1 — Build dependencies:**

```bash
# Ubuntu
sudo apt install -y build-essential libssl-dev python3 python3-pip
pip3 install pandas scipy

# macOS
brew install openssl@3 python
pip3 install pandas scipy
```

**Step 2 — Clone and build the NIST tool:**

```bash
git clone https://github.com/usnistgov/SP800-90B_EntropyAssessment.git nist_entropy
cd nist_entropy/cpp && make non_iid && cd ../..
```

**Step 3 — Compile the data generator:**

```bash
# Ubuntu
g++ -O2 -std=c++17 -o generate_data generate_data.cpp -lssl -lcrypto

# macOS (Apple Silicon)
OPENSSL_PREFIX=$(brew --prefix openssl@3)
g++ -O2 -std=c++17 -I$OPENSSL_PREFIX/include -L$OPENSSL_PREFIX/lib \
    -o generate_data generate_data.cpp -lssl -lcrypto
```

**Step 4 — Generate sequences:**

```bash
mkdir -p data
./generate_data
```

**Step 5 — Run NIST assessment:**

```bash
./run_nist.sh
# Output: results.txt
```

**Step 6 — Analyze:**

```bash
python3 analyze_results.py
```

**Step 7 — Run standalone estimator demos:**

```bash
# Runs / t-Runs / Collision estimators
g++ -O2 -std=c++17 -o runs_entropy runs_entropy_estimation.cpp
./runs_entropy

# Index-Value Coincidence estimator
g++ -O2 -std=c++17 -o ivc index_value_coincidence.cpp
./ivc
```

---

## NIST SP 800-90B Tool Reference

The `run_nist.sh` script calls `ea_non_iid` with flags `-v -a`:

| Flag | Meaning |
|---|---|
| (no flag) | Run the non-IID entropy assessment |
| `-v` | Verbose output — prints all 10 estimator results |
| `-a` | Binary data format — input is raw bytes representing 1-bit symbols |

For multi-bit symbol datasets (4-bit, 8-bit), call `ea_non_iid` directly and pass the bits-per-symbol as a positional argument:

```bash
# 4-bit symbols
./nist_entropy/cpp/ea_non_iid output_4bit_unpacked/sequence_0000.bin 4

# 8-bit symbols
./nist_entropy/cpp/ea_non_iid output_8bit/sequence_0000.bin 8
```

For IID testing (Datasets 1, 2), use `ea_iid`:

```bash
./nist_entropy/cpp/ea_iid -i -a <file> <bits_per_symbol>
```

Full NIST documentation: [https://github.com/usnistgov/SP800-90B_EntropyAssessment](https://github.com/usnistgov/SP800-90B_EntropyAssessment)

---

## Dataset Summary

| # | File | Generator | Alphabet | Samples | Distribution | Theoretical H_inf |
|---|---|---|---|---|---|---|
| 1 | `generate_data.cpp` | AES-256-CBC | Binary | 1,000,000 | Uniform | 1.000 bit |
| 1b | `aes_full.cpp` | AES-128-CBC | Binary | 1,000,000 | Uniform | 1.000 bit |
| 2 | `gen_biased_binary.cpp` | MT19937 | Binary | 1,000,000 | P(1)=0.70 | ~0.515 bits |
| 3 | `gen_4bit.cpp` | MT19937 | 4-bit (k=16) | 250,000 | P(0)=0.25, rest=0.05 | ~2.0 bits |
| 4 | `gen_8bit.cpp` | MT19937 | 8-bit (k=256) | 125,000 | P(0)=0.06, rest≈0.0037 | ~4.06 bits |
| 5 | `gen_noniid.cpp` | MT19937 + parity | Binary | 1,000,000 | Uniform marginal, non-IID | ~1.0 bits (marginal) |
