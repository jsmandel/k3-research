k := FiniteField(2);
G := GL(3, k);
R<[x]> := PolynomialRing(k,3);
V, f := GModule(G, R, 3);
GV := ActionGroup(V);

// makes weighted projective space
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

smoothList := [];
num_orbits := 0;
smooth_count := 0;

for a in GVorbreps do
    ap := ringHom(a@@f);
    
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
        for b in orbit do
            bp := ringHom(Random(b@@q)@@g);

            surface_poly:=w^2+ap*w+bp;
            k3_surface := Scheme(P,surface_poly); // degree 2 K3 surface in P(1,1,1,3)
            if IsNonsingular(k3_surface) then
                Append(~smoothList, k3_surface);
                smooth_count := smooth_count + 1;
            end if;
        end for;
        num_orbits := num_orbits + #orbit;
    else
        b := Random(W);
        bp := b@@q;
        surface_poly:=w^2+ap*w+bp;
        k3_surface := Scheme(P,surface_poly); // degree 2 K3 surface in P(1,1,1,3)
        if IsNonsingular(k3_surface) then
            Append(~smoothList, k3_surface);
            smooth_count := smooth_count +1;
        end if;
        num_orbits := num_orbits + 1;
    end if;
    smooth_count;
end for;

// Smooth List
file := Open("smoothList.m", "w");
Puts(file, "smoothList := " cat Sprint(smoothList, "Magma") cat ";");
Flush(file);
smooth_count;

// > k:=GF(2);
// > Weights:=[1,1,1,3];
// > P<x,y,z,w>:=WeightedProjectiveSpace(k,Weights);
// > a:=x^3-2*x*y^2+z^3;
// > b:=x^6-5*x^2*y*z^3+y^5*z-y^3*z^3;
// > f:=w^2+a*w+b;
// > S:=Scheme(P,f); // degree 2 K3 surface in P(1,1,1,3)
// > S;
// Scheme over GF(2) defined by
// x^6 + x^3*w + x^2*y*z^3 + y^5*z + y^3*z^3 + z^3*w + w^2
// > IsNonsingular(S);
// true
// > 