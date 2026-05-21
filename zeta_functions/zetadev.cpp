#include <cstdio>
#include <cstdint>
#include <fstream>
#include <vector>
#include "tableio.h"

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


int main() {
    // auto v = load_k3s("/Users/jordanmandel/Desktop/allSmoothList.bin");
    // printf("loaded %zu surfaces\n", v.size());
    // for (int s = 0; s < 3; s++) {
    //     printf("surface %d  a:", s);
    //     for (int i = 0; i < 10; i++) printf(" %d", v[s].a_coef[i]);
    //     printf("\n           b:");
    //     for (int i = 0; i < 28; i++) printf(" %d", v[s].b_coef[i]);
    //     printf("\n");
    // }

    auto v = load_k3s("/Users/jordanmandel/Desktop/allSmoothList.bin");
    printf("loaded %zu surfaces\n", v.size());

    // ----- Test over F_2 -----
    mult = read_table(2, 2, std::string("mult_2"));
    printf("\n=== F_2 evaluation of surface 0 ===\n");
    int pts[7][3] = {{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}};
    for (int k = 0; k < 7; k++) {
        unsigned f3 = eval_at(v[0].a_coef, CUBIC_EXP,  10, pts[k][0], pts[k][1], pts[k][2]);
        unsigned f6 = eval_at(v[0].b_coef, SEXTIC_EXP, 28, pts[k][0], pts[k][1], pts[k][2]);
        printf("  [%u:%u:%u]  f3=%u  f6=%u\n", pts[k][0], pts[k][1], pts[k][2], f3, f6);
    }

}