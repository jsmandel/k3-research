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

orbitsList := [**];    // full orbits for each a orbit representative (only includes b polynomials)

//////////////////////////////////////////////
// Creating the list of orbits, separated into the 22 a orbit representatives
//////////////////////////////////////////////
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

        // creates lists of orbit representatives for b, separated by the 22 a representatives
        orbitfora := [**];
        for rep in orbit do
            b := (rep[1])@@q;
            bp := b@@g;
            Append(~orbitfora, bp);
        end for;
        Append(~orbitsList, orbitfora);
    else
        W, g := GModule(stabG, R, 6);

        Wa := sub<W | [g((a@@f)*(c@@f) + (c@@f)^2) : c in V]>;
        Q, q := quo<W | Wa>;
        
        orbitfora := [**];
        for rep in Q do
            b := (rep)@@q;
            bp := b@@g;
            Append(~orbitfora, bp);
        end for;
        Append(~orbitsList, orbitfora);
    end if;

end for;

//////////////////////////////////////////////
// Creating a list of all smooth degree 2 K3s
//////////////////////////////////////////////

smooth_count := 0;
smoothList := [**];
for i in [1..#GVorbreps] do
    // choose a and corresponding b orbits
    a := (GVorbreps[i])@@f;
    blist := orbitsList[i];

    // loops through all b orbit representatives
    for j in [1..#blist] do
        ap := ringHom(a);
        pair := [*a*];
        bp := ringHom(blist[j]);
        surface_poly:=w^2+ap*w+bp;  // creating polynomial for K3
        k3_surface := Scheme(P,surface_poly); // degree 2 K3 surface in P(1,1,1,3)
        if IsNonsingular(k3_surface) then
            Append(~pair, blist[j]);
            Append(~smoothList, pair);
            smooth_count := smooth_count + 1;
        end if;
    end for;
end for;

//////////////////////////////////////////////
// Converting a pair (a,b) of polynomials as vectors into bytes.
// We need to store our data as multiples of 8. But a is a degree 3 polynomial, 
// corresponding to a vector of length 10, and b is a degree 6 polynomial, 
// corresponding to a vector of length 28, so total we have 38 0s and 1s. 
// Concatenate two mores 0s to get a list of 40 and then forget the last 2 0s when we convert back from bytes to the list
//////////////////////////////////////////////

// Given a and b as module elements, turn a, b into lists, concatenate them, 
// and then concatenate [0,0] to get a list of 40
polypairtolist := function(a,b)
  lista := [a[i] : i in [1..10]];
  listb := [b[i] : i in [1..28]];
  seqout :=[];
  seqout cat:=lista;
  seqout cat:=listb;
  seqout cat:=[0,0];
  return seqout;
end function;

// Given a list of 40 0s and 1s representing (a,b) as above,
// returns a customized serialization into bytes.
serialize := function(f)
    byteseq := [];
    for i in [0..4] do
        n := 1;
        bitstring := 0;
        for j in [1..8] do
            if f[i*8 + j] eq 1 then
                bitstring := BitwiseOr(bitstring, ShiftLeft(1, 8-n));
            end if;
            n := n + 1;
        end for;
        Append(~byteseq, bitstring);
    end for;

    return byteseq;
end function;

//////////////////////////////////////////////
// writing binary file
//////////////////////////////////////////////

B, h :=GModule(G, R, 6);
binFile := "allSmoothList.bin";
out := Open(binFile, "wb");

for pair in smoothList do
    a := f(pair[1]);
    b := h(pair[2]);
    bytes := serialize(polypairtolist(a,b));
    WriteBytes(out, bytes);
end for;
delete out;