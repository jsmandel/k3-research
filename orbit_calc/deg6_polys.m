k := FiniteField(2);

// ring of polynomials 
R<[x]> := PolynomialRing(k,3);

// matrix for group action
G := GL(3, k);
W,g:=GModule(G,R,6);

deg6PolyList := [];
for w in W do
    Append(~deg6PolyList, w);
end for;
#deg6PolyList;
f := Open("deg6Polys.m", "w");
Puts(f, "deg6PolyList := " cat Sprint(deg6PolyList) cat ";");
Flush(f);