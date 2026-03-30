/*
 * Author: Huseyn Kishiyev
 * Dataset 3: 4-bit Near-Uniform Biased Sequence Generator (H = 0.5 bits/bit)
 *
 * Alphabet: 4-bit symbols 0x0 – 0xF (16 possible values).
 * Distribution:
 *   P(0x0) = 0.25
 *   P(all other 15 symbols) = 0.05 each
 *
 * Generates 200 sequences of length 250,000 samples.
 *
 * IMPORTANT FOR NIST:
 *   Each 4-bit sample is written as ONE BYTE whose value is in {0,1,...,15}.
 *   We do NOT pack two nibbles into one byte.
 *
 * Output: output_4bit_unpacked/sequence_XXXX.bin
 *
 * Compile:
 *   g++ -O2 -std=c++17 -o gen_4bit_unpacked gen_4bit_unpacked.cpp
 *
 * Run:
 *   mkdir -p output_4bit_unpacked && ./gen_4bit_unpacked
 *
 * NIST:
 *   ./ea_iid -i -a output_4bit_unpacked/sequence_0000.bin 4
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <random>
#include <vector>
#include <array>

static const int NUM_SEQUENCES   = 200;
static const int SAMPLES_PER_SEQ = 250000;

// Symbol 0 has probability 0.25, symbols 1..15 each have probability 0.05.
static std::array<double, 16> make_cdf()
{
    std::array<double, 16> cdf{};
    cdf[0] = 0.25;
    for (int i = 1; i < 16; i++) {
        cdf[i] = cdf[i - 1] + 0.05;
    }
    cdf[15] = 1.0;
    return cdf;
}

static inline unsigned char sample_nibble(const std::array<double, 16>& cdf, double u)
{
    for (unsigned char s = 0; s < 16; s++) {
        if (u < cdf[s]) return s;
    }
    return 15;
}

int main()
{
    std::random_device rd;
    auto cdf = make_cdf();
    std::vector<unsigned char> buf(SAMPLES_PER_SEQ);

    for (int seq = 0; seq < NUM_SEQUENCES; seq++) {
        uint64_t seed = (static_cast<uint64_t>(rd()) << 32) | static_cast<uint64_t>(rd());
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> uni(0.0, 1.0);

        for (int i = 0; i < SAMPLES_PER_SEQ; i++) {
            buf[i] = sample_nibble(cdf, uni(rng)); // one 4-bit symbol stored in one byte
        }

        char filename[128];
        std::snprintf(filename, sizeof(filename),
                      "output_4bit_unpacked/sequence_%04d.bin", seq);

        FILE *fp = std::fopen(filename, "wb");
        if (!fp) {
            std::perror(filename);
            return 1;
        }

        size_t written = std::fwrite(buf.data(), 1, buf.size(), fp);
        std::fclose(fp);

        if (written != buf.size()) {
            std::fprintf(stderr, "Failed to write complete file for sequence %d\n", seq);
            return 1;
        }

        if ((seq + 1) % 20 == 0) {
            std::fprintf(stdout, "4-bit unpacked: generated %d / %d sequences\n",
                         seq + 1, NUM_SEQUENCES);
        }
    }

    std::fprintf(stdout, "Done. %d sequences written to output_4bit_unpacked/\n",
                 NUM_SEQUENCES);
    return 0;
}
