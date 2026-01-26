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

//////////////////////////////////////////////
// Now, I'll make the list of smooth K3s coming from the 2nd a in the list above
// There is a small change: Before, I was saving the polynomials in x, y, z, but to 
// convert them to the vector of 0s and 1s, I need them as polynomials in x[1],
// x[2], x[3]. So this is slightly different than what we emailed about before.
//////////////////////////////////////////////

a := (GVorbreps[2])@@f;
blist := orbitsList[2];

smooth_count := 0;
smoothList := [**];
for i in [1..#blist] do
ap := ringHom(a);
pair := [*a*];
bp := ringHom(blist[i]);
    surface_poly:=w^2+ap*w+bp;
    k3_surface := Scheme(P,surface_poly); // degree 2 K3 surface in P(1,1,1,3)
    if IsNonsingular(k3_surface) then
        Append(~pair, blist[i]);
        Append(~smoothList, pair);
        smooth_count := smooth_count + 1;
    end if;
end for;

//////////////////////////////////////////////
// Here is code to convert a pair (a,b) of polynomials as vectors into bytes.
// A byte is some magical way of taking 8 bits and turning it into a number,
// so to convert to bytes, we need to store our data as multiples of 8. But a is
// a degree 3 polynomial, corresponding to a vector of length 10, and b is a 
// degree 6 polynomial, corresponding to a vector of length 28, so total we have
// 38 0s and 1s. I'm just concatenating two mores 0s to get a list of 40, and then when we convert back from bytes to the list, we just forget the last 2 0s.
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

// Takes a byte sequence and returns the corresponding pair (a,b) of degree 3,
// degree 6 polynomials as module elements
deserialize := function(byteseq)    seq := [];
    for byte in byteseq do
seq cat:=[BitwiseAnd(ShiftRight(byte, 8-i), 1) : i in [1..8]];;
    end for;
   
    alist := [seq[i] : i in [1..10]];
    blist := [seq[i] : i in [11..38]];
    return [alist, blist];
end function;

//////////////////////////////////////////////
// writing binary file
//////////////////////////////////////////////

B, h :=GModule(G, R, 6);
binFile := "smoothList.bin";
out := Open(binFile, "wb");

// for pair in smoothList do
//     a := pair[1];
//     b := pair[2];
//     bytes := serialize(polypairtolist(a,b));
//     Write(out, bytes);
// end for;

testa := f((smoothList[5])[1]);
testa;
testb := h((smoothList[5])[2]); 
testb;
bytes := serialize(polypairtolist(testa,testb));
WriteBytes(out, bytes);
delete out;

// delete out;
print "Wrote";

test := Open(binFile, "rb");
smoothBin := ReadBytes(test, 5);
delete test;

smoothResult := [];
// pos := 1;
block := smoothBin[1..5];
result := deserialize(block);
resulta := result[1]; resultb := result[2];
polya := &+[resulta[i]*MonomialsOfDegree(R,3)[i] : i in [1..10]];
polyb := &+[resultb[i]*MonomialsOfDegree(R,6)[i] : i in [1..28]];
polya; polyb;