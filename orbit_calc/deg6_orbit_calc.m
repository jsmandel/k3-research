k := FiniteField(2);

// ring of polynomials 
R<[x]> := PolynomialRing(k,3);

// matrix for group action
G := GL(3, k);

V,f:=GModule(G,R,3);
W,g:=GModule(G,R,6);

load "deg6Polys.m";

for w in deg6PolyList o
    orbit := [];
    if w in deg6PolyList then
        Append(orbit, w);
        for a in deg6PolyList do
            for M in G do
                for c in deg3PolyList do
                    term := a*M;
                    termPoly := term@@f;
                    termPoly := termPoly*c + c*c;
                    term2 := w*M;
                    termPoly := term2@@g + termPoly;
                    orbitp := g(termPoly);
                end for;
            end for;
            i := Index(deg6PolyList, orbitp);
            if i ne 0 then
                Remove(~deg6PolyList, i);
            end if;
        end for;
        Append(~orbitsList, orbit);
    end if;
end for;