k := FiniteField(2);

// ring of polynomials 
R<[x]> := PolynomialRing(k,3);

// matrix for group action
G := GL(3, k);

W,g:=GModule(G,R,6);
IsIsomorphism(f);
IsIsomorphism(g);


for w in W do
    orbit := [];
    if w in uniqueOrbitPolys then
        Append(orbit, w);
        for a in uniqueOrbitPolys do
            for M in G do
                for c in polys do
                    term := a*M;
                    termPoly := term@@f;
                    termPoly := termPoly*c + c*c;
                    term2 := w*M;
                    termPoly := term2@@g + termPoly;
                    orbitp := g(termPoly);
                end for;
            end for;
            i := Index(uniqueOrbitPolys, orbitp);
            if i ne 0 then
                Remove(~uniqueOrbitPolys, i);
            end if;
        end for;
        Append(~orbitsList, orbit);
    end if;
end for;

#orbitsList;