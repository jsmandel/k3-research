
// First, load the orbits of K3s.
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
orbitRepList := [**];  // (a,b) pairs for orbit representatives
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
            orbitRep := [*ap*];
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
            orbitRep := [*ap*];
            Append(~orbitRep, bp);
            Append(~orbitRepList, orbitRep);
            Append(~orbitfora, bp);
        end for;
        Append(~orbitsList, orbitfora);
        num_orbits := num_orbits + #Q;  // for now
    end if;

end for;

// Generate bases for V (degree 3) and W (degree 6)
V_basis := MonomialsOfDegree(R, 3);
W_basis := MonomialsOfDegree(R, 6);
V_space := VectorSpace(k, #V_basis);
W_space := VectorSpace(k, #W_basis);


// Define the action of a matrix M in GL_3(F_2) on a polynomial f
ApplyMatrix := function(f, M)
    vars := MonomialsOfDegree(R,1);
    // Substitute x[i] with sum_j M[i,j]*x[j]
    newvars := [ &+[M[i,j]*vars[j] : j in [1..3]] : i in [1..3] ];
    return Evaluate(f, newvars);
end function;

// Convert a polynomial to a vector in k^n based on a given monomial basis
PolyToVector := function(f, basis, space)
    return space ! [MonomialCoefficient(f, m) : m in basis];
end function;

// Convert a vector in k^n back to a polynomial
VectorToPoly := function(v, basis)
    return &+[ v[i]*basis[i] : i in [1..#basis] ];
end function;

// Compute the stabilizer subgroup (as a set of ordered pairs)
ComputeStabilizer := function(a, b)
    // Find Stab_G(a), the stabilizer subgroup of a in G
    Ha := [ M : M in G | ApplyMatrix(a, M) eq a ];
    
    // Construct the matrix for phi_a : V -> W, phi_a(c) = c^2 + ac
    phia_matrixrows := [];
    for mon in V_basis do
        imgpoly := mon^2 + a*mon;
        Append(~phia_matrixrows, PolyToVector(imgpoly, W_basis, W_space));
    end for;
    phia_mat := Matrix(phia_matrixrows); 
    
    Stabilizer := [];
    
    // Run through Ha and find valid c in V
    for M in Ha do
        target := b + ApplyMatrix(b, M); 
        targetvec := PolyToVector(target, W_basis, W_space);
        
        // Solve cvec * phia_mat = targetvec
        // IsConsistent returns: boolean (is solvable), a specific solution, and the nullspace
        ok, c_vec, N := IsConsistent(phia_mat, targetvec);
        
        if ok then
            // If solvable, add (c, M) for all c in the solution space
            for n in N do
                c_sol_vec := c_vec + n;
                c_sol := VectorToPoly(c_sol_vec, V_basis);
                Append(~Stabilizer, <c_sol, M>);
            end for;
        end if;
    end for;
    
    return Stabilizer;
end function;


a := GVorbreps[22]@@f;
blist := orbitsList[22];
maxorblist := [];

i := 1;
for bpoly in blist do
    if (#ComputeStabilizer(a,bpoly) eq 42) then
        pair := [a,bpoly];
        Append(~maxorblist, pair);
        i;
    end if;
    i := i + 1;
end for;

for maxpair in maxorblist do
    maxpair;
    ap := ringHom(a);
    bp := ringHom(maxpair[2]);
    surface_poly:=w^2+ap*w+bp;
    k3_surface := Scheme(P,surface_poly);
    IsNonsingular(k3_surface);
end for;