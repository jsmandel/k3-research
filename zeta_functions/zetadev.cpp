#include <cstdio>
#include <cstdint>
#include <fstream>
#include <vector>

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
}