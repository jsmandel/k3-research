k := FiniteField(2);

// ring of polynomials 
R<[x]> := PolynomialRing(k,3);

// matrix for group action
G := GL(3, k);
W,g:=GModule(G,R,6);

// mons := MonomialsOfDegree(R, 6);
deg6PolyList := [];

// for coeffs in CartesianPower(k, 28) do
//     poly := R!0;  // Initialize a polynomial to zero
//     for i in [1..28] do
//         poly +:= coeffs[i] * mons[i];  // Add the monomial corresponding to the coefficients
//     end for;
//     Append(~deg6PolyList, poly);  // Add the constructed polynomial to the list
// end for;


for w in W do
    x := w@@g;
    Append(~deg6PolyList, x);
end for;
#deg6PolyList;

f := Open("deg6Polys2.m", "w");
Puts(f, "deg6PolyList := " cat Sprint(deg6PolyList) cat ";");
Flush(f);