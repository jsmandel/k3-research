// k3_count.cpp
// Self-contained K3 point counting over F_{2^n} for n = 1..12.
// Surface: w^2 + L(x,y,z) w + M(x,y,z) = 0 in P(3,1,1,1).

#include <cstdio>
#include <cstdint>
#include <vector>
#include <cassert>

// =====================================================================
// Irreducible polynomials over F_2 of degree 1..12.
// Each is encoded as a bitmask of its coefficients (bit i = coeff of t^i).
// E.g., t^2 + t + 1 = 0b111 = 7.
// =====================================================================
static const uint64_t IRRED_POLY[13] = {
    0,                          // n=0 unused
    /*  1 */ 0b11,              // t + 1
    /*  2 */ 0b111,             // t^2 + t + 1
    /*  3 */ 0b1011,            // t^3 + t + 1
    /*  4 */ 0b10011,           // t^4 + t + 1
    /*  5 */ 0b100101,          // t^5 + t^2 + 1
    /*  6 */ 0b1000011,         // t^6 + t + 1
    /*  7 */ 0b10000011,        // t^7 + t + 1
    /*  8 */ 0b100011101,       // t^8 + t^4 + t^3 + t^2 + 1
    /*  9 */ 0b1000010001,      // t^9 + t^4 + 1
    /* 10 */ 0b10000001001,     // t^10 + t^3 + 1
    /* 11 */ 0b100000000101,    // t^11 + t^2 + 1
    /* 12 */ 0b1000001010011    // t^12 + t^6 + t^4 + t + 1
};

// =====================================================================
// F_{2^n} arithmetic.  Elements are uint64_t with the low n bits used.
// =====================================================================

static int N_DEG;          // current n (degree of extension)
static uint64_t MOD_POLY;  // irreducible polynomial p(t) for F_{2^n}
static uint64_t Q_SIZE;    // 2^n

static inline uint64_t f_add(uint64_t a, uint64_t b) { return a ^ b; }

// Carry-less multiply of two polynomials over F_2 (small-degree, so just naive).
static uint64_t f_mul_raw(uint64_t a, uint64_t b) {
    uint64_t r = 0;
    while (b) {
        if (b & 1) r ^= a;
        a <<= 1;
        b >>= 1;
    }
    return r;
}

// Reduce a polynomial of degree < 2n modulo MOD_POLY, returning element of F_{2^n}.
static uint64_t f_reduce(uint64_t a) {
    int deg_mod = N_DEG;  // degree of MOD_POLY
    // For each bit at position >= deg_mod, XOR in MOD_POLY shifted appropriately.
    for (int i = 2*N_DEG - 1; i >= deg_mod; i--) {
        if ((a >> i) & 1ULL) {
            a ^= (MOD_POLY << (i - deg_mod));
        }
    }
    return a;
}

static inline uint64_t f_mul(uint64_t a, uint64_t b) {
    return f_reduce(f_mul_raw(a, b));
}

static inline uint64_t f_sqr(uint64_t a) { return f_mul(a, a); }

// Fast exponentiation a^e in F_{2^n}.
static uint64_t f_pow(uint64_t a, uint64_t e) {
    uint64_t r = 1;
    while (e) {
        if (e & 1) r = f_mul(r, a);
        a = f_mul(a, a);
        e >>= 1;
    }
    return r;
}

// Multiplicative inverse via Fermat's little theorem: a^(2^n - 2).
static inline uint64_t f_inv(uint64_t a) {
    assert(a != 0);
    return f_pow(a, Q_SIZE - 2);
}

static inline uint64_t f_div(uint64_t a, uint64_t b) {
    return f_mul(a, f_inv(b));
}

// Trace from F_{2^n} to F_2: tr(c) = c + c^2 + c^4 + ... + c^(2^(n-1)).
// In char 2 this is in {0, 1}. Returns 0 or 1.
static uint64_t f_trace(uint64_t c) {
    uint64_t s = 0;
    uint64_t x = c;
    for (int i = 0; i < N_DEG; i++) {
        s ^= x;
        x = f_sqr(x);
    }
    // s should now lie in F_2 = {0, 1}.
    return s;  // 0 or 1
}

// =====================================================================
// Surface: K3 over F_2 defined by L (cubic) and M (sextic).
//
//   L = x^3 + x^2 y + x^2 z + x y z + y^3 + y^2 z + z^3
//   M = x^6 + x^5 z + x^4 y z + x^2 y^4 + y^6 + y^5 z + z^6
// =====================================================================

static const unsigned L_COEFFS[10] = {
    1, 1, 1, 0, 1, 0, 1, 1, 0, 1
};
static const unsigned L_EXP[10][3] = {
    {3,0,0},{2,1,0},{2,0,1},{1,2,0},{1,1,1},
    {1,0,2},{0,3,0},{0,2,1},{0,1,2},{0,0,3}
};

static const unsigned M_COEFFS[28] = {
    1,0,1,0,1,0,0,0,0,0, 1,0,0,0,0, 0,0,0,0,0,0, 1,1,0,0,0,0,1
};
static const unsigned M_EXP[28][3] = {
    {6,0,0},{5,1,0},{5,0,1},{4,2,0},{4,1,1},{4,0,2},
    {3,3,0},{3,2,1},{3,1,2},{3,0,3},
    {2,4,0},{2,3,1},{2,2,2},{2,1,3},{2,0,4},
    {1,5,0},{1,4,1},{1,3,2},{1,2,3},{1,1,4},{1,0,5},
    {0,6,0},{0,5,1},{0,4,2},{0,3,3},{0,2,4},{0,1,5},{0,0,6}
};

// Precompute powers a^0, a^1, ..., a^maxExp into out[].
static void precompute_powers(uint64_t a, int maxExp, uint64_t* out) {
    out[0] = 1;
    for (int e = 1; e <= maxExp; e++) out[e] = f_mul(out[e-1], a);
}

// Evaluate a homogeneous polynomial with F_2 coeffs at (y0,y1,y2),
// using precomputed power tables p0,p1,p2 up to degree 6.
static uint64_t eval_poly(const unsigned* coeffs, const unsigned exps[][3],
                          int num_terms,
                          const uint64_t* p0, const uint64_t* p1, const uint64_t* p2) {
    uint64_t acc = 0;
    for (int t = 0; t < num_terms; t++) {
        if (!coeffs[t]) continue;
        uint64_t term = f_mul(f_mul(p0[exps[t][0]], p1[exps[t][1]]), p2[exps[t][2]]);
        acc ^= term;
    }
    return acc;
}

// Number of w in F_{2^n} solving w^2 + L0*w + M0 = 0.
//   - L0 = 0:  exactly 1 root (w = sqrt(M0); Frobenius is bijective).
//   - L0 != 0: substitute w = L0*t to get t^2 + t + M0/L0^2 = 0.
//              has 2 roots iff Tr(M0/L0^2) = 0, else 0 roots.
static int contribution(uint64_t L0, uint64_t M0) {
    if (L0 == 0) return 1;
    uint64_t L0sq = f_sqr(L0);
    uint64_t c = f_div(M0, L0sq);
    return (f_trace(c) == 0) ? 2 : 0;
}

// =====================================================================
// Frobenius orbits over F_2: orbits of squaring on F_{2^n}.
// Returns rep[a] = orbit representative of a, size[r] = orbit size at rep r.
// =====================================================================
static void compute_frobenius_orbits(std::vector<uint64_t>& rep,
                                     std::vector<uint64_t>& size) {
    rep.assign(Q_SIZE, 0);
    size.assign(Q_SIZE, 0);
    std::vector<char> seen(Q_SIZE, 0);
    for (uint64_t a = 0; a < Q_SIZE; a++) {
        if (seen[a]) continue;
        // Trace orbit of a under squaring.
        uint64_t x = a;
        uint64_t cnt = 0;
        do {
            seen[x] = 1;
            rep[x] = a;
            cnt++;
            x = f_sqr(x);
        } while (x != a);
        size[a] = cnt;
    }
}

// =====================================================================
// Main loop.
// =====================================================================
int main() {
    for (int n = 1; n <= 12; n++) {
        N_DEG = n;
        MOD_POLY = IRRED_POLY[n];
        Q_SIZE = 1ULL << n;

        std::vector<uint64_t> rep, osize;
        compute_frobenius_orbits(rep, osize);

        long long count = 0;
        uint64_t p0[7], p1[7], p2[7];

        // (1) The point (0:0:1).
        precompute_powers(0, 6, p0);
        precompute_powers(0, 6, p1);
        precompute_powers(1, 6, p2);
        {
            uint64_t L0 = eval_poly(L_COEFFS, L_EXP, 10, p0, p1, p2);
            uint64_t M0 = eval_poly(M_COEFFS, M_EXP, 28, p0, p1, p2);
            count += contribution(L0, M0);
        }

        // (2) Hyperplane at infinity: (0:1:y_2), y_2 in K.
        precompute_powers(0, 6, p0);
        precompute_powers(1, 6, p1);
        for (uint64_t y2 = 0; y2 < Q_SIZE; y2++) {
            if (rep[y2] != y2) continue;
            precompute_powers(y2, 6, p2);
            uint64_t L0 = eval_poly(L_COEFFS, L_EXP, 10, p0, p1, p2);
            uint64_t M0 = eval_poly(M_COEFFS, M_EXP, 28, p0, p1, p2);
            count += (long long)contribution(L0, M0) * (long long)osize[y2];
        }

        // (3) Affine patch: (1:y_1:y_2).
        precompute_powers(1, 6, p0);
        for (uint64_t y1 = 0; y1 < Q_SIZE; y1++) {
            if (rep[y1] != y1) continue;
            precompute_powers(y1, 6, p1);
            for (uint64_t y2 = 0; y2 < Q_SIZE; y2++) {
                precompute_powers(y2, 6, p2);
                uint64_t L0 = eval_poly(L_COEFFS, L_EXP, 10, p0, p1, p2);
                uint64_t M0 = eval_poly(M_COEFFS, M_EXP, 28, p0, p1, p2);
                count += (long long)contribution(L0, M0) * (long long)osize[y1];
            }
        }

        printf("|X(F_%llu)| = %lld\n", (unsigned long long)Q_SIZE, count);
        fflush(stdout);
    }
    return 0;
}