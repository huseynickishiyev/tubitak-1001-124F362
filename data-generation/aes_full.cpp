/*
 * Dataset 1: AES-CBC Uniform Binary Sequence Generator for NIST SP800-90B IID Tests
 *
 * Generates 200 binary sequences of length 1,000,000 samples each.
 * Uses AES-128 in CBC mode via OpenSSL to produce pseudorandom bytes.
 * Each output bit is unpacked into one byte with value 0x00 or 0x01.
 *
 * Output format for NIST SP800-90B ea_iid:
 *   Each sequence is written as a raw binary file: sequence_XXXX.bin
 *   where XXXX is the zero-padded sequence index (0000–0199).
 *
 * Each output file contains exactly 1,000,000 bytes:
 *   each byte is one binary sample in {0,1}.
 *
 * Compile:
 *   g++ -O2 -o gen_aes_cbc_unpacked gen_aes_cbc_unpacked.cpp -lssl -lcrypto
 *
 * Run:
 *   mkdir -p output_aes_cbc_unpacked && ./gen_aes_cbc_unpacked
 */

#include <openssl/evp.h>
#include <openssl/rand.h>
#include <cstdio>
#include <cstdlib>
#include <vector>

static const int NUM_SEQUENCES    = 200;
static const int SEQUENCE_BITS    = 1000000;                  // desired number of binary samples
static const int PACKED_BYTES     = (SEQUENCE_BITS + 7) / 8;  // 125000 bytes
static const int AES_BLOCK_SIZE_  = 16;
static const int PADDED_BYTES =
    ((PACKED_BYTES + AES_BLOCK_SIZE_ - 1) / AES_BLOCK_SIZE_) * AES_BLOCK_SIZE_; // 125008

/*
 * Encrypt 'out_len' zero bytes with AES-128-CBC using the given key and IV.
 * Since CBC with padding disabled requires a multiple-of-16 input length,
 * out_len must be divisible by 16.
 */
static bool aes_cbc_keystream(const unsigned char *key,
                              const unsigned char *iv,
                              unsigned char *out,
                              int out_len)
{
    EVP_CIPHER_CTX *ctx = EVP_CIPHER_CTX_new();
    if (!ctx) return false;

    if (1 != EVP_EncryptInit_ex(ctx, EVP_aes_128_cbc(), nullptr, key, iv)) {
        EVP_CIPHER_CTX_free(ctx);
        return false;
    }

    EVP_CIPHER_CTX_set_padding(ctx, 0);

    std::vector<unsigned char> plaintext(out_len, 0x00);

    int len1 = 0;
    int len2 = 0;

    if (1 != EVP_EncryptUpdate(ctx, out, &len1, plaintext.data(), out_len)) {
        EVP_CIPHER_CTX_free(ctx);
        return false;
    }

    if (1 != EVP_EncryptFinal_ex(ctx, out + len1, &len2)) {
        EVP_CIPHER_CTX_free(ctx);
        return false;
    }

    EVP_CIPHER_CTX_free(ctx);
    return (len1 + len2 == out_len);
}

int main()
{
    std::vector<unsigned char> packed(PADDED_BYTES);
    std::vector<unsigned char> unpacked(SEQUENCE_BITS); // 1 byte per bit/sample

    for (int seq = 0; seq < NUM_SEQUENCES; ++seq) {
        unsigned char key[16];
        unsigned char iv[16];

        if (1 != RAND_bytes(key, sizeof(key)) ||
            1 != RAND_bytes(iv, sizeof(iv))) {
            std::fprintf(stderr, "RAND_bytes failed\n");
            return 1;
        }

        if (!aes_cbc_keystream(key, iv, packed.data(), PADDED_BYTES)) {
            std::fprintf(stderr, "AES-CBC encryption failed at sequence %d\n", seq);
            return 1;
        }

        // Unpack first 1,000,000 bits into bytes {0,1}
        int out_idx = 0;
        for (int i = 0; i < PACKED_BYTES && out_idx < SEQUENCE_BITS; ++i) {
            unsigned char byte = packed[i];
            for (int b = 7; b >= 0 && out_idx < SEQUENCE_BITS; --b) {
                unpacked[out_idx++] = (byte >> b) & 0x01;
            }
        }

        char filename[128];
        std::snprintf(filename, sizeof(filename),
                      "output_aes_cbc_unpacked/sequence_%04d.bin", seq);

        FILE *fp = std::fopen(filename, "wb");
        if (!fp) {
            std::perror(filename);
            return 1;
        }

        size_t written = std::fwrite(unpacked.data(), 1, unpacked.size(), fp);
        std::fclose(fp);

        if (written != unpacked.size()) {
            std::fprintf(stderr, "Failed to write complete file for sequence %d\n", seq);
            return 1;
        }

        if ((seq + 1) % 20 == 0) {
            std::fprintf(stdout, "AES-CBC unpacked: generated %d / %d sequences\n",
                         seq + 1, NUM_SEQUENCES);
        }
    }

    std::fprintf(stdout, "Done. %d sequences written to output_aes_cbc_unpacked/\n",
                 NUM_SEQUENCES);

    return 0;
}
