k := FiniteField(2);

// ring of polynomials 
R<[x]> := PolynomialRing(k,3);

// matrix for group action
G := GL(3, k);
W,g:=GModule(G,R,6);

load "deg3Polys.m";
load "uniqueOrbitPolys3.m";
load "deg6Polys.m";

uniqueOrbitList6 := [];
orbitList := [];
orbit := [];
pair := [];

for a in uniqueOrbitList3 do
    Append(~pair, a);
    for w in W do
        x = g(w);
        string = "W: (" cat Sprint(x) cat "), ";
        if string in deg6PolyList then
            Append(~pair, w);
            for M in G do
            for c in deg3PolyList do
                poly = x*M;
                poly = poly@@g;
                aVec = f(a);
                term = aVec*M;
                term = term@@f;
                term = term * c;
                poly = poly + term;
                poly = poly + (c*c);
                string = "W: " cat Sprint(poly);
                i := Index(deg6PolyList, string);
                if i ne 0 then
                    Remove(~deg6PolyList, i);
                end if;
            end for;
        end for;
        end if;



        Append(~pair, w);
        
        Append(~orbit, pair);
    end for;
    Append(~orbitList, orbit)
end for;


// for w in deg6PolyList do
//     orbit := [];
//     if w in uniqueOrbitList3 then
//         Append(orbit, w);
//         for a in deg6PolyList do
//             for M in G do
//                 for c in deg3PolyList do
//                     term := a*M;
//                     termPoly := term@@g;
//                     termPoly := termPoly*c + c*c;
//                     term2 := w*M;
//                     termPoly := term2@@g + termPoly;
//                     orbitp := g(termPoly);
//                 end for;
//             end for;
//             i := Index(deg6PolyList, orbitp);
//             if i ne 0 then
//                 Remove(~deg6PolyList, i);
//             end if;
//         end for;
//         Append(~orbitsList, orbit);
//     end if;
// end for;