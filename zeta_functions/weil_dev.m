// ============================================================
//  Weil polynomial computation from point counts
//  (test mode: first NUM_TO_RUN surfaces only)
// ============================================================

p := 2;
P<x> := PolynomialRing(Rationals());
C<i> := ComplexField(30);
PC := PolynomialRing(C);

load "/Users/jordanmandel/k3-research/zeta_functions/pointCountList.m";

TOL := 0.0001;
NUM_TO_RUN := 10;        // <-- change this to #pointcountlist for full run

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
// build the full degree-22 char poly of Frob on H^2 and return the
// "modified" polynomial p^(-22) * f(p*x) whose roots should have |.| = 1. ---
BuildModified := function(coeffs11, eps)
    coeffs := coeffs11;                       // length 11
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
WeilPolys := [];                     // list of <K3index, WeilPoly, sign>
NeedMore  := [ [] : k in [1..11] ];  // NeedMore[k]: surfaces needing k extra counts
                                     //  (i.e. point counts up through F_{p^(11+k)})

// -------- main loop (test mode: only first NUM_TO_RUN surfaces) --------

upper := Min(NUM_TO_RUN, #pointcountlist);
printf "Running on surfaces 1..%o (of %o total)\n", upper, #pointcountlist;

for j in [1..upper] do
    FpPoints := pointcountlist[j];
    N := #FpPoints;   // should be 11

    // Traces: T_k = #X(F_{p^k}) - p^{2k} - 1
    Traces := [ Rationals() ! (FpPoints[k] - p^(2*k) - 1) : k in [1..N] ];

    // Newton's identities: first N coefficients of char poly of Frob.
    CharPolyCoeffs := [ -Traces[1] ];
    for k in [2..N] do
        ck := -(1/k) * ( Traces[k]
              + &+[ CharPolyCoeffs[i]*Traces[k-i] : i in [1..k-1] ] );
        Append(~CharPolyCoeffs, ck);
    end for;

    determined := false;
    resolved_poly := 0;
    resolved_sign := 0;

    if CharPolyCoeffs[11] ne 0 then
        // c_11 != 0 forces sign = +1
        modPoly, _ := BuildModified(CharPolyCoeffs, 1);
        resolved_poly := modPoly;
        resolved_sign := 1;
        determined := true;
    else
        // c_11 = 0: try both signs, check moduli
        modPolyPlus,  _ := BuildModified(CharPolyCoeffs,  1);
        modPolyMinus, _ := BuildModified(CharPolyCoeffs, -1);
        plus_ok  := AllUnitModulus(RightModulus(modPolyPlus));
        minus_ok := AllUnitModulus(RightModulus(modPolyMinus));

        if plus_ok and not minus_ok then
            resolved_poly := modPolyPlus;
            resolved_sign := 1;
            determined := true;
        elif minus_ok and not plus_ok then
            resolved_poly := modPolyMinus;
            resolved_sign := -1;
            determined := true;
        elif plus_ok and minus_ok then
            // Ambiguous - find smallest disambiguating k > 11.
            // (Diff at k = 22-i is 2 * c_i * p^{2k-22}, nonzero iff c_i != 0.)
            found  := false;
            need_k := 0;
            for i in [10..1 by -1] do
                if CharPolyCoeffs[i] ne 0 then
                    need_k := 22 - i;
                    found := true;
                    break;
                end if;
            end for;
            if not found then
                need_k := 22;
            end if;
            extra := need_k - 11;         // in [1..11]
            Append(~NeedMore[extra], j);
        else
            printf "WARNING: surface %o has neither sign consistent.\n", j;
        end if;
    end if;

    if determined then
        Append(~WeilPolys, <j, resolved_poly, resolved_sign>);
    end if;
end for;

// -------- report --------
printf "\n===== SUMMARY =====\n";
printf "Surfaces processed:    %o\n", upper;
printf "Weil poly determined:  %o\n", #WeilPolys;
total_pending := 0;
for k in [1..11] do
    total_pending +:= #NeedMore[k];
end for;
printf "Ambiguous (need more): %o\n", total_pending;
for k in [1..11] do
    printf "  need counts up to F_{p^%o}: %o surfaces\n", 11+k, #NeedMore[k];
end for;

// -------- save results --------

// Named aliases for readability
Need12 := NeedMore[1];
Need13 := NeedMore[2];
Need14 := NeedMore[3];
Need15 := NeedMore[4];
Need16 := NeedMore[5];
Need17 := NeedMore[6];
Need18 := NeedMore[7];
Need19 := NeedMore[8];
Need20 := NeedMore[9];
Need21 := NeedMore[10];
Need22 := NeedMore[11];

// Write the WeilPolys and the aggregate NeedMore list
PrintFile("weil_results.m",
    "WeilPolys := " * Sprint(WeilPolys) * ";\n" : Overwrite := true);
PrintFile("weil_results.m",
    "NeedMore := " * Sprint(NeedMore) * ";\n");

// Write each individual NeedMore list as its own file
PrintFile("need_p12.m", "Need12 := " * Sprint(Need12) * ";\n" : Overwrite := true);
PrintFile("need_p13.m", "Need13 := " * Sprint(Need13) * ";\n" : Overwrite := true);
PrintFile("need_p14.m", "Need14 := " * Sprint(Need14) * ";\n" : Overwrite := true);
PrintFile("need_p15.m", "Need15 := " * Sprint(Need15) * ";\n" : Overwrite := true);
PrintFile("need_p16.m", "Need16 := " * Sprint(Need16) * ";\n" : Overwrite := true);
PrintFile("need_p17.m", "Need17 := " * Sprint(Need17) * ";\n" : Overwrite := true);
PrintFile("need_p18.m", "Need18 := " * Sprint(Need18) * ";\n" : Overwrite := true);
PrintFile("need_p19.m", "Need19 := " * Sprint(Need19) * ";\n" : Overwrite := true);
PrintFile("need_p20.m", "Need20 := " * Sprint(Need20) * ";\n" : Overwrite := true);
PrintFile("need_p21.m", "Need21 := " * Sprint(Need21) * ";\n" : Overwrite := true);
PrintFile("need_p22.m", "Need22 := " * Sprint(Need22) * ";\n" : Overwrite := true);

printf "\nDone. Files written:\n";
printf "  weil_results.m  (WeilPolys, NeedMore)\n";
printf "  need_p12.m ... need_p22.m  (individual bucket lists)\n";