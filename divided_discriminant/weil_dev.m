// weil_dev.m -- Weil polynomials via divided discriminant sign.
// Reads:
//   pointCountList.m  (defines p, pointcountlist)
//   discriminants20.m   (defines rvals :: rvals[j] = r mod 8 for surface j)

load "../zeta_functions/pointCountList.m";
load "discriminants20.m";    // defines rvals (list or assoc array, r in {3,7})

P<x> := PolynomialRing(Rationals());

MAXK := 11;
NUMSURFACES := 20;         // set to #pointcountlist for full run

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
WeilPolys := AssociativeArray();    // j -> <sign, CharPoly>
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
        continue;
    end if;

    // consistency check: if c11 != 0, sign MUST be +1
    if c[11] ne 0 and eps ne 1 then
        printf "Surface #%o: MISMATCH c11 = %o != 0 but disc gives sign %o\n",
               j, c[11], eps;
        Append(~SignMismatch, j);
    end if;

    CharPoly := BuildWeil(c, eps);
    WeilPolys[j] := <eps, CharPoly>;
end for;

// ---------- report ----------
printf "\n========== SUMMARY ==========\n";
printf "Surfaces processed:  %o\n", nproc;
printf "Weil polys computed: %o\n", #Keys(WeilPolys);
printf "Sign +1: %o\n", #[ j : j in Keys(WeilPolys) | WeilPolys[j][1] eq  1 ];
printf "Sign -1: %o\n", #[ j : j in Keys(WeilPolys) | WeilPolys[j][1] eq -1 ];
printf "Errors (r not 3/7):  %o  %o\n", #Errors, Errors;
printf "Consistency mismatches: %o  %o\n", #SignMismatch, SignMismatch;

// ---------- write results ----------
out := "weilResults.m";
fh := Open(out, "w");
fprintf fh, "// Weil polynomials via divided discriminant.\n";
fprintf fh, "// weilPolys[j] = <sign, [c1..c22]>, P(T)=T^22+c1 T^21+...+c22.\n";
fprintf fh, "weilPolys := AssociativeArray();\n";
for j in Sort(Setseq(Keys(WeilPolys))) do
    eps := WeilPolys[j][1]; f := WeilPolys[j][2];
    coeffs := [ Coefficient(f, 22 - m) : m in [1..22] ];
    fprintf fh, "weilPolys[%o] := <%o, %o>;\n", j, eps, coeffs;
end for;
delete fh;
printf "\nWrote results to %o\n", out;

// ---------- write polynomial results ----------
outPoly := "weilResultsPolys.m";
fhPoly := Open(outPoly, "w");

fprintf fhPoly, "// Weil polynomials via divided discriminant.\n";
fprintf fhPoly, "P<x> := PolynomialRing(Rationals());\n";
fprintf fhPoly, "weilPolys := AssociativeArray();\n";

for j in Sort(Setseq(Keys(WeilPolys))) do
    eps := WeilPolys[j][1];
    f := WeilPolys[j][2];
    fprintf fhPoly, "weilPolys[%o] := <%o, %o>;\n", j, eps, f;
end for;

delete fhPoly;
printf "Wrote polynomial results to %o\n", outPoly;