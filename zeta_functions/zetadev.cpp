#include <cstdio>
#include <cstdint>
#include <fstream>
#include <vector>
#include "tableio.h"
#include <array>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <numeric>

// Exponent array for a (degree 3) 
static const int CUBIC_EXP[10][3] = {
    {3,0,0}, {2,1,0}, {2,0,1}, {1,2,0}, {1,1,1},
    {1,0,2}, {0,3,0}, {0,2,1}, {0,1,2}, {0,0,3}
};
// Exponent array for b (degree 6) 
static const int SEXTIC_EXP[28][3] = {
    {6,0,0}, {5,1,0}, {5,0,1}, {4,2,0}, {4,1,1}, {4,0,2}, {3,3,0},
    {3,2,1}, {3,1,2}, {3,0,3}, {2,4,0}, {2,3,1}, {2,2,2}, {2,1,3},
    {2,0,4}, {1,5,0}, {1,4,1}, {1,3,2}, {1,2,3}, {1,1,4}, {1,0,5},
    {0,6,0}, {0,5,1}, {0,4,2}, {0,3,3}, {0,2,4}, {0,1,5}, {0,0,6}
};

struct K3 {
    uint8_t a_coef[10];   // 0/1 coefficients of cubic
    uint8_t b_coef[28];   // 0/1 coefficients of sextic
};

// Loading in binary file
std::vector<K3> load_k3s(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) { fprintf(stderr, "Cannot open %s\n", path.c_str()); std::exit(1); }
    std::vector<K3> out;
    uint8_t buf[5];
    while (in.read(reinterpret_cast<char*>(buf), 5)) {
        uint8_t bits[40];
        for (int byte = 0; byte < 5; byte++)
            for (int j = 0; j < 8; j++)
                bits[byte*8 + j] = (buf[byte] >> (7 - j)) & 1;
        K3 s;
        for (int i = 0; i < 10; i++) s.a_coef[i] = bits[i];
        for (int i = 0; i < 28; i++) s.b_coef[i] = bits[10 + i];
        out.push_back(s);
    }
    return out;
}


unsigned** mult;
unsigned** quadratic_roots;

unsigned eval_at(const uint8_t* coef, const int exp_table[][3], int nmono,
                 unsigned x0, unsigned x1, unsigned x2) {
    // cache powers up to 6
    unsigned p0[7] = {1}, p1[7] = {1}, p2[7] = {1};
    for (int e = 1; e <= 6; e++) {
        p0[e] = mult[p0[e-1]][x0];
        p1[e] = mult[p1[e-1]][x1];
        p2[e] = mult[p2[e-1]][x2];
    }
    unsigned acc = 0;
    for (int i = 0; i < nmono; i++) {
        if (!coef[i]) continue;
        unsigned m = mult[ p0[exp_table[i][0]] ][
                     mult[ p1[exp_table[i][1]] ][ p2[exp_table[i][2]] ] ];
        acc ^= m;     // sum in characteristic 2
    }
    return acc;
}

int count_points(const K3& s, unsigned q) {
    auto fibre = [&](unsigned x, unsigned y, unsigned z) -> int {
        unsigned f3 = eval_at(s.a_coef, CUBIC_EXP, 10, x, y, z);
        unsigned f6 = eval_at(s.b_coef, SEXTIC_EXP, 28, x, y, z);
        return (int)quadratic_roots[f3][f6];
    };
    int total = 0;
    total += fibre(0, 0, 1);
    for (unsigned z = 0; z < q; z++) total += fibre(0, 1, z);
    for (unsigned y = 0; y < q; y++)
        for (unsigned z = 0; z < q; z++)
            total += fibre(1, y, z);
    return total;
}

// orbit optimization: ask
int count_points_fast(const K3& s, unsigned q,
    const unsigned* orbit_rep,
    const unsigned* orbit_size) {
        auto fibre = [&](unsigned x, unsigned y, unsigned z) -> int {
            unsigned f3 = eval_at(s.a_coef, CUBIC_EXP, 10, x, y, z);
            unsigned f6 = eval_at(s.b_coef, SEXTIC_EXP, 28, x, y, z);
            return (int)quadratic_roots[f3][f6];
        };
        int total = fibre(0, 0, 1);
        for (unsigned z = 0; z < q; z++)
            if (orbit_rep[z] == z)
                total += fibre(0, 1, z) * (int)orbit_size[z];
        for (unsigned y = 0; y < q; y++) {
            if (orbit_rep[y] != y) continue;
            for (unsigned z = 0; z < q; z++)
                total += fibre(1, y, z) * (int)orbit_size[y];
        }
        return total;
}


int main() {
    auto v = load_k3s("/Users/jordanmandel/Desktop/allSmoothList.bin");
    printf("loaded %zu surfaces\n", v.size());
    for (int s = 0; s < 3; s++) {
        printf("surface %d  a:", s);
        for (int i = 0; i < 10; i++) printf(" %d", v[s].a_coef[i]);
        printf("\n           b:");
        for (int i = 0; i < 28; i++) printf(" %d", v[s].b_coef[i]);
        printf("\n");
    }

    printf("surface %d  a:", 10512);
    for (int i = 0; i < 10; i++) printf(" %d", v[10512].a_coef[i]);
    printf("\n           b:");
    for (int i = 0; i < 28; i++) printf(" %d", v[10512].b_coef[i]);

    auto surfaces = load_k3s("/Users/jordanmandel/Desktop/allSmoothList.bin");
    // printf("loaded %zu surfaces\n", v.size());

    // // ----- Test over F_2 -----
    mult = read_table(2, 2, std::string("mult_2"));
    printf("\n=== F_2 evaluation of surface 0 ===\n");
    int pts[7][3] = {{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}};
    for (int k = 0; k < 7; k++) {
        unsigned f3 = eval_at(surfaces[0].a_coef, CUBIC_EXP,  10, pts[k][0], pts[k][1], pts[k][2]);
        unsigned f6 = eval_at(surfaces[0].b_coef, SEXTIC_EXP, 28, pts[k][0], pts[k][1], pts[k][2]);
        printf("  [%u:%u:%u]  f3=%u  f6=%u\n", pts[k][0], pts[k][1], pts[k][2], f3, f6);
    }

    const int NUM_SURFACES_TO_TEST = 10;
    const int Q_VALUES[] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096};   // start small; 16 is q=2^4
    const int NUM_QS = sizeof(Q_VALUES) / sizeof(Q_VALUES[0]);

    int total_tests = 0, mismatches = 0;
    double naive_time_total = 0.0, fast_time_total = 0.0;

    for (int qi = 0; qi < NUM_QS; qi++) {
        unsigned q = Q_VALUES[qi];
        std::string qq = std::to_string(q);
        printf("\n--- Testing q = %u ---\n", q);

        // Load tables for this q
        mult            = read_table(q, q, "mult_"            + qq);
        quadratic_roots = read_table(q, q, "quadratic_roots_" + qq);
        unsigned* orbit_rep  = read_table(q, "orbit_rep_"  + qq);
        unsigned* orbit_size = read_table(q, "orbit_size_" + qq);

        for (int s = 0; s < NUM_SURFACES_TO_TEST && s < (int)surfaces.size(); s++) {
            const K3& surf = surfaces[s];

            auto t1 = std::chrono::high_resolution_clock::now();
            int naive = count_points(surf, q);
            auto t2 = std::chrono::high_resolution_clock::now();
            int fast  = count_points_fast(surf, q, orbit_rep, orbit_size);
            auto t3 = std::chrono::high_resolution_clock::now();

            double naive_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
            double fast_ms  = std::chrono::duration<double, std::milli>(t3 - t2).count();
            naive_time_total += naive_ms;
            fast_time_total  += fast_ms;

            total_tests++;
            bool ok = (naive == fast);
            if (!ok) mismatches++;

            // Weil-bound sanity check: |N - (q^2+q+1)| <= 22 q
            long expected_mid = (long)q*q + q + 1;
            long weil_radius  = 22L * q;
            bool weil_ok = std::abs((long)naive - expected_mid) <= weil_radius;

            printf("  surface %2d: naive=%d  fast=%d  %s  weil_%s  (naive %.2fms, fast %.2fms)\n",
                   s, naive, fast,
                   ok ? "OK" : "*** MISMATCH ***",
                   weil_ok ? "ok" : "FAIL",
                   naive_ms, fast_ms);
        }

        const K3& surf = surfaces[10512];

        auto t1 = std::chrono::high_resolution_clock::now();
        int naive = count_points(surf, q);
        auto t2 = std::chrono::high_resolution_clock::now();
        int fast  = count_points_fast(surf, q, orbit_rep, orbit_size);
        auto t3 = std::chrono::high_resolution_clock::now();

        double naive_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
        double fast_ms  = std::chrono::duration<double, std::milli>(t3 - t2).count();
        naive_time_total += naive_ms;
        fast_time_total  += fast_ms;

            total_tests++;
            bool ok = (naive == fast);
            if (!ok) mismatches++;

            // Weil-bound sanity check: |N - (q^2+q+1)| <= 22 q
            long expected_mid = (long)q*q + q + 1;
            long weil_radius  = 22L * q;
            bool weil_ok = std::abs((long)naive - expected_mid) <= weil_radius;

            printf("  surface %2d: naive=%d  fast=%d  %s  weil_%s  (naive %.2fms, fast %.2fms)\n",
                10512, naive, fast,
                   ok ? "OK" : "*** MISMATCH ***",
                   weil_ok ? "ok" : "FAIL",
                   naive_ms, fast_ms);
        

        // Free tables for this q
        // (slight leak: we should also delete the row pointers; ignore for testing)
        delete[] orbit_rep;
        delete[] orbit_size;
    }

    printf("\n=== Summary ===\n");
    printf("Total tests: %d\n", total_tests);
    printf("Mismatches:  %d\n", mismatches);
    printf("Naive total time: %.1f ms\n", naive_time_total);
    printf("Fast total time:  %.1f ms\n", fast_time_total);
    if (fast_time_total > 0)
        printf("Speedup factor:   %.2fx\n", naive_time_total / fast_time_total);

    return mismatches == 0 ? 0 : 1;   // nonzero exit if any test failed

    

}