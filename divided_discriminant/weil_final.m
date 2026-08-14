// weil_dev.m -- Weil polynomials via divided discriminant sign.
// Reads:
//   pointCountList.m  (defines p, pointcountlist)
//   discriminants20.m   (defines rvals :: rvals[j] = r mod 8 for surface j)

// CHANGE TO PATH WHERE RUNNING CODE
load "/var/lib/private/f006vkh/m2/discriminants.m";    // defines rvals (list or assoc array, r in {3,7})
load "/var/lib/private/f006vkh/m2/pointCountList.m";


P<x> := PolynomialRing(Rationals());

MAXK := 11;
NUMSURFACES := #pointcountlist;         // set to #pointcountlist for full run

// Build full degree-22 Weil poly from c[1..11] and sign eps.
BuildWeil := function(c11, eps)
    c := c11;
    for n in [12..21] do
        k := 22 - n;
        Append(~c, eps * p^(22 - 2*k) * c[k]);
    end for;
    Append(~c, eps * p^22);
    return x^22 + &+[ c[m]*x^(22-m) : m in [1..22] ];
end function;

// map r mod 8 to sign of functional equation
SignFromR := function(r)
    rr := r mod 8;
    if rr eq 7 then return 1; end if;
    if rr eq 3 then return -1; end if;
    return 0;   // error sentinel: should never happen
end function;

// ---------- main ----------
WeilPolys := [];                    //  <sign, CharPoly, CharPoly factorized>
Errors    := [];                    // surfaces where r not in {3,7} mod 8
SignMismatch := [];                 // c11 != 0 but disc-sign != +1 (consistency check)

nproc := Min(NUMSURFACES, #pointcountlist);

for j in [1..nproc] do
    FpPoints := pointcountlist[j];
    kmax := Min(MAXK, #FpPoints);

    Traces := [ Rationals()!(FpPoints[k] - p^(2*k) - 1) : k in [1..kmax] ];
    c := [ -Traces[1] ];
    for k in [2..kmax] do
        ck := -(1/k)*(Traces[k] + &+[ c[m]*Traces[k-m] : m in [1..k-1] ]);
        Append(~c, ck);
    end for;

    eps := SignFromR(rvals[j]);
    if eps eq 0 then
        printf "Surface #%o: ERROR r = %o mod 8 (expected 3 or 7)\n", j, rvals[j] mod 8;
        Append(~Errors, j);
        Append(~WeilPolys, <0, 0, 0>);   // placeholder keeps indices aligned
        continue;
    end if;

    // consistency check: if c11 != 0, sign MUST be +1
    if c[11] ne 0 and eps ne 1 then
        printf "Surface #%o: MISMATCH c11 = %o != 0 but disc gives sign %o\n",
               j, c[11], eps;
        Append(~SignMismatch, j);
    end if;

    CharPoly := BuildWeil(c, eps);
    CharPolyModified := p^(-22)*Evaluate(CharPoly,p*x);
    CharPolyFactorized := Factorization(CharPolyModified);
    Append(~WeilPolys, <eps, CharPolyModified, CharPolyFactorized>);

    // ---- progress report every 1000 surfaces ----
    if j mod 1000 eq 0 then
        printf "Computed %o / %o Weil polynomials\n", j, nproc;
    end if;
end for;

// ---------- report ----------
printf "\n========== SUMMARY ==========\n";
printf "Surfaces processed:  %o\n", nproc;
printf "Weil polys computed: %o\n", #WeilPolys;
printf "Sign +1: %o\n", #[ j : j in [1..#WeilPolys] | WeilPolys[j][1] eq  1 ];
printf "Sign -1: %o\n", #[ j : j in [1..#WeilPolys] | WeilPolys[j][1] eq -1 ];
printf "Errors (r not 3/7):  %o  %o\n", #Errors, Errors;
printf "Consistency mismatches: %o  %o\n", #SignMismatch, SignMismatch;

// ---------- write results (single sequence: sign, poly, factorization) ----------
out := "weilResults.m";
fh := Open(out, "w");
fprintf fh, "// Weil polynomials via divided discriminant.\n";
fprintf fh, "// weilData[j] = <sign, WeilPoly, factorization>.\n";
fprintf fh, "P<x> := PolynomialRing(Rationals());\n";
fprintf fh, "weilData := [\n";
for idx in [1..#WeilPolys] do
    eps   := WeilPolys[idx][1];
    f     := WeilPolys[idx][2];
    facts := WeilPolys[idx][3];
    fprintf fh, "<%o, %o, %o>", eps, f, facts;
    if idx lt #WeilPolys then fprintf fh, ",\n"; else fprintf fh, "\n"; end if;
end for;
fprintf fh, "];\n";
delete fh;
printf "\nWrote results to %o\n", out;