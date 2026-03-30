/*
    Author: Huseyn Kishiyev
    Program is designed to serve as a data generator for 
    "Observations on NIST Entropy Suite, (Aslan et al.)"
    and implements the data generation method stated in the paper.
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sys/stat.h> // mkdir

// OpenSSL
#include <openssl/evp.h>
#include <openssl/rand.h>
#include <openssl/err.h>

static void handleErrors() {
    ERR_print_errors_fp(stderr);
    abort();
}

/*
  Generate 1,000,000 binary symbols (bits), stored as 1 byte per bit (0x00/0x01).
  Source: AES-256-CBC encryption of an all-zero plaintext with random key+IV.
  We disable padding and produce enough ciphertext bytes, then expand to bits.

  Bit extraction order: MSB first in each ciphertext byte (bit 7 down to 0).
  This is conventional and deterministic.
*/
static void generate_cbc_binary_bits_1B_per_bit(
    const std::string& filename,
    const int num_bits /* e.g., 1'000'000 */
) {
    // --- 1) Key & IV ---
    unsigned char key[32];  // 256-bit key
    unsigned char iv[16];   // 128-bit IV
    if (RAND_bytes(key, sizeof key) != 1) handleErrors();
    if (RAND_bytes(iv,  sizeof iv)  != 1) handleErrors();

    // --- 2) Cipher context (AES-256-CBC, no padding) ---
    EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
    if (!ctx) handleErrors();

    if (EVP_EncryptInit_ex(ctx, EVP_aes_256_cbc(), nullptr, key, iv) != 1) handleErrors();
    if (EVP_CIPHER_CTX_set_padding(ctx, 0) != 1) handleErrors();

    // --- 3) Zero plaintext sized to a whole number of blocks >= num_bits ---
    const int block_bytes = 16;          // AES block = 128 bits = 16 bytes
    const int block_bits  = block_bytes * 8;
    const int needed_bits_rounded = ((num_bits + block_bits - 1) / block_bits) * block_bits;
    const int needed_bytes        = needed_bits_rounded / 8;

    std::vector<unsigned char> zero_pt(needed_bytes, 0x00);
    std::vector<unsigned char> ct(needed_bytes + block_bytes); // out buffer

    int outlen = 0, tot = 0;
    if (EVP_EncryptUpdate(ctx, ct.data(), &outlen, zero_pt.data(), (int)zero_pt.size()) != 1) handleErrors();
    tot += outlen;

    int fin = 0;
    if (EVP_EncryptFinal_ex(ctx, ct.data() + tot, &fin) != 1) handleErrors();
    tot += fin;

    EVP_CIPHER_CTX_free(ctx);

    if (tot < needed_bytes) {
        std::cerr << "Unexpected: ciphertext shorter than requested\n";
        handleErrors();
    }

    // --- 4) Expand to 1 byte per bit (0x00 or 0x01), MSB-first per byte ---
    std::vector<unsigned char> bit_bytes(num_bits);
    int bit_idx = 0;
    for (int i = 0; i < needed_bytes && bit_idx < num_bits; ++i) {
        unsigned char b = ct[i];
        for (int k = 7; k >= 0 && bit_idx < num_bits; --k) {
            unsigned char bit = (b >> k) & 0x01;
            bit_bytes[bit_idx++] = bit;          // 0x00 or 0x01
        }
    }

    // --- 5) Exactly num_bits bytes ---
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        std::cerr << "Error opening file for writing: " << filename << "\n";
        return;
    }
    out.write(reinterpret_cast<const char*>(bit_bytes.data()), bit_bytes.size());
    out.close();
}

int main() {
    const int NUM_SEQUENCES = 200;
    const int SEQ_BITS = 1'000'000; // one million *bits* -> one million *bytes* in output

    mkdir("data", 0755);

    std::cout << "Generating " << NUM_SEQUENCES
              << " IID full-entropy binary sequences (AES-256-CBC)...\n";

    for (int i = 0; i < NUM_SEQUENCES; ++i) {
        std::string filename = "data/sequence_" + std::to_string(i) + ".bin";
        generate_cbc_binary_bits_1B_per_bit(filename, SEQ_BITS);
        std::cout << "Generated " << filename << " (" << SEQ_BITS << " bytes)\n";
    }

    std::cout << "Done.\n";
    return 0;
}
