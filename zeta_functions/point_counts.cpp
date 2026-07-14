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

// simple helper to free a q x q table
static void free_table(unsigned** T, unsigned q) {
    if (!T) return;
    for (unsigned r = 0; r < q; r++) delete[] T[r];
    delete[] T;
}


int main() {
    const int P_CHAR = 2;
    auto surfaces = load_k3s("/Users/jordanmandel/Desktop/allSmoothList.bin");
    printf("loaded %zu surfaces\n", surfaces.size());

    // q = p^1 .. p^11
    const int Q_VALUES[] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048};
    const int NUM_QS = sizeof(Q_VALUES) / sizeof(Q_VALUES[0]);

    std::vector<std::vector<long long>> counts(
        surfaces.size(), std::vector<long long>(NUM_QS, 0));

    for (int qi = 0; qi < NUM_QS; qi++) {
        unsigned q = Q_VALUES[qi];
        std::string qq = std::to_string(q);
        printf("Computing q = %u (k = %d) ...\n", q, qi + 1);

        mult            = read_table(q, q, "mult_"            + qq);
        quadratic_roots = read_table(q, q, "quadratic_roots_" + qq);
        unsigned* orbit_rep  = read_table(q, "orbit_rep_"  + qq);
        unsigned* orbit_size = read_table(q, "orbit_size_" + qq);

        auto t0 = std::chrono::high_resolution_clock::now();
        for (size_t s = 0; s < surfaces.size(); s++) {
            counts[s][qi] = (long long)
                count_points_fast(surfaces[s], q, orbit_rep, orbit_size);
            if ((s % 1000) == 0 && qi >= 6) {
                printf("    surface %zu / %zu\n", s, surfaces.size());
            }
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        double sec = std::chrono::duration<double>(t1 - t0).count();
        printf("  done q = %u in %.1f s\n", q, sec);

        free_table(mult, q);            mult = nullptr;
        free_table(quadratic_roots, q); quadratic_roots = nullptr;
        delete[] orbit_rep;
        delete[] orbit_size;
    }

    // ---- Write Magma file ----
    std::ofstream out("/Users/jordanmandel/k3-research/zeta_functions/pointCountList.m");
    out << "// Auto-generated by count_k3.cpp\n";
    out << "// pointcountlist[j][k] = #X(F_{p^k}) for K3 with Magma index j.\n";
    out << "// Magma index j corresponds to C++ surface index (j-1).\n";
    out << "// Number of surfaces: " << surfaces.size() << "\n";
    out << "// k ranges over 1.." << NUM_QS
        << "  (so q ranges over p^1 .. p^" << NUM_QS << ")\n";
    out << "p := " << P_CHAR << ";\n";
    out << "pointcountlist := [\n";
    for (size_t s = 0; s < surfaces.size(); s++) {
        out << "  [ ";
        for (int qi = 0; qi < NUM_QS; qi++) {
            out << counts[s][qi];
            if (qi + 1 < NUM_QS) out << ", ";
        }
        out << " ]";
        if (s + 1 < surfaces.size()) out << ",";
        out << "\n";
    }
    out << "];\n";
    out.close();

    printf("Wrote %zu surfaces x %d point counts to pointCountList.m\n",
           surfaces.size(), NUM_QS);
    return 0;
}