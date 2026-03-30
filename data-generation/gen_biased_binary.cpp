/* 
 * Author: Huseyn Kishiyev
 * Dataset 2: Biased Binary Sequence Generator
 *
 * Distribution:
 *   P(1) = 0.7
 *   P(0) = 0.3
 *
 * Generates 200 sequences of length 1,000,000 binary samples.
 *
 * IMPORTANT FOR NIST:
 *   Each binary sample is stored as ONE BYTE with value 0x00 or 0x01.
 *   We do NOT pack 8 bits into one byte.
 *
 * Output: output_biased_binary_unpacked/sequence_XXXX.bin
 *
 * Compile:
 *   g++ -O2 -std=c++17 -o gen_biased_binary_unpacked gen_biased_binary_unpacked.cpp
 *
 * Run:
 *   mkdir -p output_biased_binary_unpacked && ./gen_biased_binary_unpacked
 *
 * NIST:
 *   ./ea_iid -i -a output_biased_binary_unpacked/sequence_0000.bin 1
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <random>
#include <vector>

static const int NUM_SEQUENCES  = 200;
static const int SEQUENCE_BITS  = 1000000;
static const double P_ONE       = 0.7;

int main()
{
    std::random_device rd;
    std::vector<unsigned char> buf(SEQUENCE_BITS);

    for (int seq = 0; seq < NUM_SEQUENCES; seq++) {
        uint64_t seed = (static_cast<uint64_t>(rd()) << 32) | static_cast<uint64_t>(rd());
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> uni(0.0, 1.0);

        for (int i = 0; i < SEQUENCE_BITS; i++) {
            buf[i] = (uni(rng) < P_ONE) ? 1 : 0;
        }

        char filename[128];
        std::snprintf(filename, sizeof(filename),
                      "output_biased_binary_unpacked/sequence_%04d.bin", seq);

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
            std::fprintf(stdout, "Biased binary unpacked: generated %d / %d sequences\n",
                         seq + 1, NUM_SEQUENCES);
        }
    }

    std::fprintf(stdout, "Done. %d sequences written to output_biased_binary_unpacked/\n",
                 NUM_SEQUENCES);
    return 0;
}
