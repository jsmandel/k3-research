p := 2;
P<x> := PolynomialRing(Rationals());
C<i> := ComplexField(30);
PC := PolynomialRing(C);

// Load the point counts from the C++ output.
// pointcountlist[j][k] = #X(F_{p^k}) for K3 index j (Magma is 1-indexed).
// The C++ file wrote entries with 0-indexed surface positions; index j here
// corresponds to C++ surface index (j-1).
load "/Users/jordanmandel/k3-research/zeta_functions/pointCountList.m";

TOL := 0.0001;

// --- helper: absolute values of the roots (one per conjugate pair) ---
RightModulus := function(poly)
    polyC := PC ! poly;
    rts   := Roots(polyC, C);
    R     := [];
    for r in rts do
        val := Abs(r[1]);
        for m in [1..r[2]] do
            Append(~R, val);
        end for;
    end for;
    return R;
end function;

// --- helper: given the first 11 CharPolyCoeffs and a sign eps in {+1,-1},
// build the full degree-22 char poly (of Frobenius on H^2) and return the
// "modified" polynomial p^(-22) * f(p*x) whose roots should have |.| = 1. ---
BuildModified := function(coeffs11, eps)
    coeffs := coeffs11;                       // length 11
    // coefficient indices 12..21 come from the functional equation
    // (Magma 1-indexed: coeffs[i] is coefficient of x^(22-i))
    for k in [12..21] do
        // c_k = eps * c_{22-k} * p^(2k-22)
        Append(~coeffs, eps * coeffs[22-k] * p^(2*k - 22));
    end for;
    Append(~coeffs, eps * p^22);              // c_22
    assert #coeffs eq 22;
    CharPoly         := x^22 + &+[ coeffs[i]*x^(22-i) : i in [1..22] ];
    CharPolyModified := p^(-22) * Evaluate(CharPoly, p*x);
    return CharPolyModified, coeffs;
end function;

// --- helper: are all moduli equal to 1 within TOL? ---
AllUnitModulus := function(mods)
    for n in mods do
        if Abs(n - 1) gt TOL then
            return false;
        end if;
    end for;
    return true;
end function;

// --- outputs ---
WeilPolys := [];                              // list of <K3index, WeilPoly, sign>
NeedMore  := [ [] : k in [1..11] ];           // NeedMore[k]: needs count up to F_{p^k} extra

for j in [1..#pointcountlist] do
    FpPoints := pointcountlist[j];
    N := #FpPoints;   // should be 11

    // Traces: T_k = #X(F_{p^k}) - p^{2k} - 1
    Traces := [ Rationals() ! (FpPoints[k] - p^(2*k) - 1) : k in [1..N] ];

    // Newton's identities to get the first N coefficients of char poly of Frob.
    CharPolyCoeffs := [ -Traces[1] ];
    for k in [2..N] do
        ck := -(1/k) * ( Traces[k]
              + &+[ CharPolyCoeffs[i]*Traces[k-i] : i in [1..k-1] ] );
        Append(~CharPolyCoeffs, ck);
    end for;
    // Now CharPolyCoeffs has length 11: entries c_1, ..., c_11.

    determined := false;
    resolved_poly := 0;
    resolved_sign := 0;

    if CharPolyCoeffs[11] ne 0 then
        // Functional-equation sign is FORCED to be +1 by c_11 ne 0
        // (since c_11 must satisfy c_11 = eps * c_11 * p^0 = eps * c_11)
        modPoly, _ := BuildModified(CharPolyCoeffs, 1);
        resolved_poly := modPoly;
        resolved_sign := 1;
        determined := true;
    else
        // c_11 = 0 so both signs are a priori consistent with the FE.
        // Try eps = +1 first.
        modPolyPlus, _ := BuildModified(CharPolyCoeffs, 1);
        modsPlus := RightModulus(modPolyPlus);
        plus_ok  := AllUnitModulus(modsPlus);

        modPolyMinus, _ := BuildModified(CharPolyCoeffs, -1);
        modsMinus := RightModulus(modPolyMinus);
        minus_ok  := AllUnitModulus(modsMinus);

        if plus_ok and not minus_ok then
            resolved_poly := modPolyPlus;
            resolved_sign := 1;
            determined := true;
        elif minus_ok and not plus_ok then
            resolved_poly := modPolyMinus;
            resolved_sign := -1;
            determined := true;
        elif plus_ok and minus_ok then
            // AMBIGUOUS: need more data.
            // Find the smallest index k > 11 that would first show a
            // difference between the +/- polynomials.
            // For k in [12..22], BuildModified gave
            //   c_k^(+) = c_{22-k} * p^{2k-22}
            //   c_k^(-) = -c_{22-k} * p^{2k-22}
            // These differ iff c_{22-k} != 0.  So the smallest disambiguating
            // k > 11 corresponds to the smallest i = 22-k < 11 with c_i != 0,
            // equivalently the LARGEST i < 11 (since k = 22-i, smaller k = larger i).
            //
            // Meaning: the disambiguating extra point count is over F_{p^k}
            // where k = 22 - (largest i < 11 with c_i != 0).
            //
            // Since we already have counts up to k = 11, we need EXTRA counts
            // for k = 12, 13, ..., 22. But we only need up through the first
            // disambiguating k. The "amount extra" is (k - 11).
            found := false;
            need_k := 0;
            for i in [10..1 by -1] do
                if CharPolyCoeffs[i] ne 0 then
                    need_k := 22 - i;      // this is in [12..21]
                    found := true;
                    break;
                end if;
            end for;
            if not found then
                // All of c_1..c_11 are zero: extreme edge case, would need c_22
                need_k := 22;
            end if;
            extra := need_k - 11;          // in [1..11]
            Append(~NeedMore[extra], j);
        else
            // Neither sign works: something is wrong (shouldn't happen
            // for a genuine smooth K3).
            printf "WARNING: surface %o has neither sign consistent with Weil bounds.\n", j;
        end if;
    end if;

    if determined then
        Append(~WeilPolys, <j, resolved_poly, resolved_sign>);
    end if;
end for;

// --- report ---
printf "\n===== SUMMARY =====\n";
printf "Total surfaces:        %o\n", #pointcountlist;
printf "Weil poly determined:  %o\n", #WeilPolys;
total_pending := 0;
for k in [1..11] do
    total_pending +:= #NeedMore[k];
end for;
printf "Ambiguous (need more): %o\n", total_pending;
for k in [1..11] do
    printf "  need counts up to F_{p^%o}: %o surfaces\n", 11+k, #NeedMore[k];
end for;

// --- save results ---
PrintFile("weil_results.m", "WeilPolys := " * Sprint(WeilPolys) * ";\n" :
          Overwrite := true);
PrintFile("weil_results.m", "NeedMore := " * Sprint(NeedMore) * ";\n");