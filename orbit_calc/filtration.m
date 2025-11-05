k := FiniteField(2);
G := GL(3, k);
R<[x]> := PolynomialRing(k,3);
V, f := GModule(G, R, 3);
GV := ActionGroup(V);

// Computes the orbits of G on the degree 3 polynomials
GVorbits := Orbits(GV);

// Returns a list of orbit representatives, as vectors
GVorbreps := [GVorbits[i][1] : i in [1..#GVorbits]];

// Returns a list of orbit representatives, as polynomials
GVorbrepspoly := [(GVorbreps[i])@@f : i in [1..#GVorbreps]];

orbitsList := [**];
orbitRepList := [**];
num_orbits := 0;

for a in GVorbreps do
    ap := a@@f;
    
    stabA := Stabilizer(GV, a); // stabilizer of a
    gvMap := GModuleAction(V); // map from G to GV
    stabG := stabA @@ gvMap; // makes the stabilizer subgroup as a subgroup of G

    // stabilizer subgroup acting on W
    if not IsTrivial(stabG) then
        W, g := GModule(stabG, R, 6);
        GW := ActionGroup(W);

        Wa := sub<W | [g((a@@f)*(c@@f) + (c@@f)^2) : c in V]>;
        Q, q := quo<W | Wa>;
        GQ := ActionGroup(Q);

        // list of all orbits of group action
        orbit := Orbits(GQ);

        // adds one orbit representative from each orbit to list
        for rep in orbit do
            orbitRep := [*a*];

            Append(~orbitRep, rep[1]@@q);
            Append(~orbitRepList, orbitRep);
        end for;
        Append(~orbitsList, orbit);
        num_orbits := num_orbits + #orbit;
    else
        W, g := GModule(stabG, R, 6);
        // GW := ActionGroup(W);

        Wa := sub<W | [g((a@@f)*(c@@f) + (c@@f)^2) : c in V]>;
        Q, q := quo<W | Wa>;
        // GQ := ActionGroup(Q);

        for b in Q do
            orbitRep := [*a*];
            Append(~orbitRep, b@@q);
            Append(~orbitRepList, orbitRep);
        end for;
        
        num_orbits := num_orbits + #Q;  // for now
    end if;

end for;

// Orbit List
// f := Open("orbitList2.m", "w");
// Puts(f, "orbitsList := " cat Sprint(orbitsList, "Magma") cat ";");
// Flush(f);

// Orbit Rep List
// f := Open("orbitReps2.m", "w");
// Puts(f, "orbitRepList := " cat Sprint(orbitRepList, "Magma") cat ";");
// Flush(f);
// num_orbits;