/*
 * Author: Huseyn Kishiyev
 * Dataset 4: 8-bit Near-Uniform Biased Sequence Generator
 *
 * Alphabet: 8-bit bytes (0x00 – 0xFF, 256 possible values).
 * Distribution:
 *   P(0x00) = 0.06
 *   P(all other 255 symbols) = 0.94 / 255
 *
 * Generates 200 sequences of length 125,000 samples.
 * One sample = one byte.
 *
 * Output: output_8bit/sequence_XXXX.bin
 *
 * Compile:
 *   g++ -O2 -std=c++17 -o gen_8bit_fixed gen_8bit_fixed.cpp
 *
 * Run:
 *   mkdir -p output_8bit && ./gen_8bit_fixed
 *
 * NIST:
 *   ./ea_iid -i -a output_8bit/sequence_0000.bin 8
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <random>
#include <vector>
#include <array>

static const int NUM_SEQUENCES   = 200;
static const int SAMPLES_PER_SEQ = 125000;

static const double P_ZERO  = 0.06;
static const double P_OTHER = (1.0 - P_ZERO) / 255.0;

static std::array<double, 256> make_cdf()
{
    std::array<double, 256> cdf{};
    cdf[0] = P_ZERO;
    for (int i = 1; i < 256; i++) {
        cdf[i] = cdf[i - 1] + P_OTHER;
    }
    cdf[255] = 1.0;
    return cdf;
}

static inline unsigned char sample_byte(const std::array<double, 256>& cdf, double u)
{
    int lo = 0, hi = 255;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (u < cdf[mid]) hi = mid;
        else lo = mid + 1;
    }
    return static_cast<unsigned char>(lo);
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
            buf[i] = sample_byte(cdf, uni(rng));
        }

        char filename[128];
        std::snprintf(filename, sizeof(filename),
                      "output_8bit/sequence_%04d.bin", seq);

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
            std::fprintf(stdout, "8-bit: generated %d / %d sequences\n",
                         seq + 1, NUM_SEQUENCES);
        }
    }

    std::fprintf(stdout, "Done. %d sequences written to output_8bit/\n",
                 NUM_SEQUENCES);
    return 0;
}
