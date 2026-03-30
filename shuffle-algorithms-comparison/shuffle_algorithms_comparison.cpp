/*
 * Huseyn Kishiyev
 * Comparison of In-Place Fisher Yates and
 * Sattolo Shuffling Algorithms 
*/
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <openssl/rand.h>
#include <openssl/err.h>

struct BitSequence {
    std::vector<uint8_t> packed;
    int n_bits;
};

BitSequence generate_bitpacked(int n_bits) {
    int n_bytes = (n_bits + 7) / 8;
    BitSequence bs; bs.n_bits = n_bits; bs.packed.resize(n_bytes);
    if (RAND_bytes(bs.packed.data(), n_bytes) != 1) {
        ERR_print_errors_fp(stderr);
        throw std::runtime_error("RAND_bytes failed");
    }
    int rem = n_bits % 8;
    if (rem) bs.packed.back() &= (uint8_t)(0xFF << (8 - rem));
    return bs;
}

std::vector<uint8_t> unpack_bits(const BitSequence& bs) {
    std::vector<uint8_t> syms(bs.n_bits);
    for (int i = 0; i < bs.n_bits; ++i)
        syms[i] = (bs.packed[i/8] >> (7-(i%8))) & 1;
    return syms;
}

BitSequence pack_bits(const std::vector<uint8_t>& syms) {
    BitSequence bs; bs.n_bits = (int)syms.size();
    bs.packed.assign((bs.n_bits+7)/8, 0);
    for (int i = 0; i < bs.n_bits; ++i)
        if (syms[i]) bs.packed[i/8] |= (uint8_t)(1 << (7-(i%8)));
    return bs;
}

template<typename T>
void fisher_yates(std::vector<T>& a, std::mt19937& rng) {
    for (int i = (int)a.size()-1; i > 0; --i) {
        std::uniform_int_distribution<int> d(0, i);
        std::swap(a[i], a[d(rng)]);
    }
}

template<typename T>
void sattolo(std::vector<T>& a, std::mt19937& rng) {
    for (int i = (int)a.size()-1; i > 0; --i) {
        std::uniform_int_distribution<int> d(0, i-1);
        std::swap(a[i], a[d(rng)]);
    }
}

struct Stats {
    double mean_ms, stddev_ms, min_ms, max_ms, q1_ms, median_ms, q3_ms;
};

Stats compute_stats(std::vector<double> ns) {
    std::sort(ns.begin(), ns.end());
    int n = ns.size();
    double sum = 0; for (auto v : ns) sum += v;
    double mean = sum/n, sq = 0;
    for (auto v : ns) sq += (v-mean)*(v-mean);
    auto q = [&](double p) {
        double pos = (n-1)*p; int lo = (int)pos, hi = std::min(lo+1,n-1);
        return ns[lo]+(ns[hi]-ns[lo])*(pos-lo);
    };
    return { mean/1e6, std::sqrt(sq/n)/1e6,
             ns.front()/1e6, ns.back()/1e6,
             q(0.25)/1e6, q(0.50)/1e6, q(0.75)/1e6 };
}

int main() {
    constexpr int N_SEQ  = 100;
    constexpr int N_BITS = 1'000'000;
    constexpr int WARMUP = 5;

    std::filesystem::create_directories("/home/claude/nist_sequences/fisher_yates");
    std::filesystem::create_directories("/home/claude/nist_sequences/sattolo");

    std::mt19937 rng(std::random_device{}());

    std::cerr << "Generating " << N_SEQ << " bit-packed AES sequences...\n";
    std::vector<BitSequence> seqs(N_SEQ);
    for (auto& s : seqs) s = generate_bitpacked(N_BITS);
    std::cerr << "Done.\n";

    // Warmup
    for (int w = 0; w < WARMUP; ++w) {
        auto tmp = unpack_bits(seqs[0]); fisher_yates(tmp, rng);
        tmp = unpack_bits(seqs[0]);      sattolo(tmp, rng);
    }

    std::vector<double> fy_ns, sa_ns;
    fy_ns.reserve(N_SEQ); sa_ns.reserve(N_SEQ);

    // ── Fisher-Yates: toplam süreyi de ölç ───────────────────────────────────
    auto fy_total_t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N_SEQ; ++i) {
        auto syms = unpack_bits(seqs[i]);
        auto t0 = std::chrono::high_resolution_clock::now();
        fisher_yates(syms, rng);
        auto t1 = std::chrono::high_resolution_clock::now();
        fy_ns.push_back(std::chrono::duration<double,std::nano>(t1-t0).count());
        auto out = pack_bits(syms);
        char fname[128];
        std::snprintf(fname,sizeof(fname),
            "/home/claude/nist_sequences/fisher_yates/seq_%03d.bin",i);
        std::ofstream f(fname,std::ios::binary);
        f.write(reinterpret_cast<const char*>(out.packed.data()),out.packed.size());
    }
    auto fy_total_t1 = std::chrono::high_resolution_clock::now();
    double fy_total_ms = std::chrono::duration<double,std::milli>(fy_total_t1-fy_total_t0).count();

    // ── Sattolo: toplam süreyi de ölç ────────────────────────────────────────
    auto sa_total_t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N_SEQ; ++i) {
        auto syms = unpack_bits(seqs[i]);
        auto t0 = std::chrono::high_resolution_clock::now();
        sattolo(syms, rng);
        auto t1 = std::chrono::high_resolution_clock::now();
        sa_ns.push_back(std::chrono::duration<double,std::nano>(t1-t0).count());
        auto out = pack_bits(syms);
        char fname[128];
        std::snprintf(fname,sizeof(fname),
            "/home/claude/nist_sequences/sattolo/seq_%03d.bin",i);
        std::ofstream f(fname,std::ios::binary);
        f.write(reinterpret_cast<const char*>(out.packed.data()),out.packed.size());
    }
    auto sa_total_t1 = std::chrono::high_resolution_clock::now();
    double sa_total_ms = std::chrono::duration<double,std::milli>(sa_total_t1-sa_total_t0).count();

    // Toplam sadece shuffle sürelerinin toplamı (unpack/repack/IO hariç)
    double fy_shuffle_sum_ms = 0, sa_shuffle_sum_ms = 0;
    for (auto v : fy_ns) fy_shuffle_sum_ms += v/1e6;
    for (auto v : sa_ns) sa_shuffle_sum_ms += v/1e6;

    auto fy = compute_stats(fy_ns);
    auto sa = compute_stats(sa_ns);

    std::cout << "\n================================================================\n";
    std::cout << "  AES Bit-Packed Shuffle Benchmark  |  NIST-ready output\n";
    std::cout << "  Sequences: " << N_SEQ
              << "  |  Bits/seq: " << N_BITS
              << "  |  Bytes/seq: " << N_BITS/8 << "\n";
    std::cout << "================================================================\n\n";

    std::cout << "── Tek sekans istatistikleri (ms) ──────────────────────────────\n";
    std::cout << std::left
              << std::setw(14) << "Algoritma"
              << std::setw(10) << "Ortalama"
              << std::setw(10) << "Ortanca"
              << std::setw(10) << "Std.Sap."
              << std::setw(8)  << "Min"
              << std::setw(8)  << "Maks"
              << "\n";
    std::cout << std::string(60,'-') << "\n";

    auto pr = [](const char* name, const Stats& s){
        std::cout << std::left << std::fixed << std::setprecision(2)
                  << std::setw(14) << name
                  << std::setw(10) << s.mean_ms
                  << std::setw(10) << s.median_ms
                  << std::setw(10) << s.stddev_ms
                  << std::setw(8)  << s.min_ms
                  << std::setw(8)  << s.max_ms
                  << "\n";
    };
    pr("Fisher-Yates", fy);
    pr("Sattolo",      sa);

    std::cout << "\n── 100 sekans toplam süresi ────────────────────────────────────\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Fisher-Yates  |  Sadece shuffle       : " << fy_shuffle_sum_ms << " ms\n";
    std::cout << "  Fisher-Yates  |  Tüm işlemler (IO+dön): " << fy_total_ms       << " ms\n";
    std::cout << "  Sattolo       |  Sadece shuffle       : " << sa_shuffle_sum_ms << " ms\n";
    std::cout << "  Sattolo       |  Tüm işlemler (IO+dön): " << sa_total_ms       << " ms\n";

    double delta_shuffle = (sa_shuffle_sum_ms - fy_shuffle_sum_ms) / fy_shuffle_sum_ms * 100.0;
    double delta_total   = (sa_total_ms - fy_total_ms) / fy_total_ms * 100.0;
    std::cout << "\n  Delta (shuffle toplam) : " << std::showpos << std::setprecision(2)
              << delta_shuffle << "%\n";
    std::cout << "  Delta (tüm işlem toplam): " << delta_total << "%\n";

    // CSV
    std::cout << std::noshowpos;
    std::cout << "\nCSV_START\n";
    std::cout << "seq,fisher_yates_ms,sattolo_ms\n";
    for (int i = 0; i < N_SEQ; ++i)
        std::cout << i << "," << std::fixed << std::setprecision(4)
                  << fy_ns[i]/1e6 << "," << sa_ns[i]/1e6 << "\n";
    std::cout << "CSV_END\n";

    return 0;
}
