/*
 * Author:Huseyn Kishiyev
 * Dataset 5: Non-IID Binary Sequence Generator
 *
 * Generation rule:
 *   Let S = (s_1, s_2, ..., s_L), L = 1,000,000.
 *   All elements are generated independently except:
 *     s_{8k} = (s_{8k-7} + s_{8k-6} + ... + s_{8k-1}) mod 2
 *
 * IMPORTANT FOR NIST:
 *   Each binary sample is stored as ONE BYTE with value 0x00 or 0x01.
 *   We do NOT pack bits into bytes.
 *
 * Output: output_noniid_unpacked/sequence_XXXX.bin
 *
 * Compile:
 *   g++ -O2 -std=c++17 -o gen_noniid_unpacked gen_noniid_unpacked.cpp
 *
 * Run:
 *   mkdir -p output_noniid_unpacked && ./gen_noniid_unpacked
 *
 * NIST:
 *   ./ea_non_iid output_noniid_unpacked/sequence_0000.bin 1
 * or
 *   ./ea_iid -i -a output_noniid_unpacked/sequence_0000.bin 1
 * if you intentionally want to show IID-test failure behavior.
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <random>
#include <vector>

static const int NUM_SEQUENCES = 200;
static const int SEQUENCE_BITS = 1000000;

int main()
{
    std::random_device rd;
    std::vector<unsigned char> bits(SEQUENCE_BITS);

    for (int seq = 0; seq < NUM_SEQUENCES; seq++) {
        uint64_t seed = (static_cast<uint64_t>(rd()) << 32) | static_cast<uint64_t>(rd());
        std::mt19937 rng(static_cast<uint32_t>(seed));
        std::uniform_int_distribution<int> bit_dist(0, 1);

        for (int i = 0; i < SEQUENCE_BITS; i++) {
            bits[i] = static_cast<unsigned char>(bit_dist(rng));
        }

        for (int k = 1; k * 8 <= SEQUENCE_BITS; k++) {
            int parity_pos = k * 8 - 1;  // 0-indexed positions 7, 15, 23, ...
            unsigned char parity = 0;
            for (int j = 1; j <= 7; j++) {
                parity ^= bits[parity_pos - j];
            }
            bits[parity_pos] = parity;
        }

        char filename[128];
        std::snprintf(filename, sizeof(filename),
                      "output_noniid_unpacked/sequence_%04d.bin", seq);

        FILE *fp = std::fopen(filename, "wb");
        if (!fp) {
            std::perror(filename);
            return 1;
        }

        size_t written = std::fwrite(bits.data(), 1, bits.size(), fp);
        std::fclose(fp);

        if (written != bits.size()) {
            std::fprintf(stderr, "Failed to write complete file for sequence %d\n", seq);
            return 1;
        }

        if ((seq + 1) % 20 == 0) {
            std::fprintf(stdout, "Non-IID unpacked: generated %d / %d sequences\n",
                         seq + 1, NUM_SEQUENCES);
        }
    }

    std::fprintf(stdout, "Done. %d sequences written to output_noniid_unpacked/\n",
                 NUM_SEQUENCES);
    return 0;
}
