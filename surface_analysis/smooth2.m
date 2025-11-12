k := FiniteField(2);
G := GL(3, k);
R<[x]> := PolynomialRing(k,3);
V, f := GModule(G, R, 3);
GV := ActionGroup(V);

Weights:=[1,1,1,3];
P<x,y,z,w>:=WeightedProjectiveSpace(k,Weights);
PRing := CoordinateRing(P);
ringHom := hom<R -> PRing | x,y,z>;

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
        orbitfora := [**];
        for rep in orbit do
            b := (rep[1])@@q;
            bp := b@@g;
            orbitRep := [*a*];
            Append(~orbitRep, bp);
            Append(~orbitRepList, orbitRep);
            Append(~orbitfora, bp);
        end for;
        Append(~orbitsList, orbitfora);
        num_orbits := num_orbits + #orbit;
    else
        W, g := GModule(stabG, R, 6);
        // GW := ActionGroup(W);

        Wa := sub<W | [g((a@@f)*(c@@f) + (c@@f)^2) : c in V]>;
        Q, q := quo<W | Wa>;
        // GQ := ActionGroup(Q);
        orbitfora := [**];
        for rep in Q do
            b := (rep)@@q;
            bp := b@@g;
            orbitRep := [*a*];
            Append(~orbitRep, bp);
            Append(~orbitRepList, orbitRep);
            Append(~orbitfora, bp);
        end for;
        Append(~orbitsList, orbitfora);
        num_orbits := num_orbits + #Q;  // for now
    end if;

end for;

a := (GVorbreps[12])@@f;
blist := orbitsList[12];

smooth_count := 0;
smoothList := [**];
t := Cputime();
for i in [1..#blist] do
    ap := ringHom(a);
    pair := [*ap*];
    bp := ringHom(blist[i]);
    surface_poly:=w^2+ap*w+bp;
    k3_surface := Scheme(P,surface_poly); // degree 2 K3 surface in P(1,1,1,3)
    if IsNonsingular(k3_surface) then
        Append(~pair, bp);
        Append(~smoothList, pair);
        smooth_count := smooth_count + 1;
    end if;
end for;
Cputime(t);