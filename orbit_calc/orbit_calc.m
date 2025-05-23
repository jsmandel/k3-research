k := FiniteField(2);

// ring of polynomials 
R<[x]> := PolynomialRing(k,3);

// matrix for group action
G := GL(3, k);

// GModule to turn polynomial into proper vector
V,f:=GModule(G,R,3);

// get all monomials of degree 3
mons := MonomialsOfDegree(R, 3);
// find all possible polynomials
polys := [];

// CartesianPower generates all possible n-tuples over a field
// so 10 variables over F2 for us
for coeffs in CartesianPower(k, 10) do
    poly := R!0;  // Initialize a polynomial to zero
    for i in [1..10] do
        poly +:= coeffs[i] * mons[i];  // Add the monomial corresponding to the coefficients
    end for;
    Append(~polys, poly);  // Add the constructed polynomial to the list
end for;
#polys;

// store degree3 polynomials
file := Open("deg3Polys.m", "w");
Puts(file, "deg3PolyList := " cat Sprint(polys) cat ";");
Flush(file);

// create list of lists for orbits
orbitsList := [];

// store polynomials we haven't seen
uniqueOrbitPolys := [];
for p in polys do
    v := f(p);
    Append(~uniqueOrbitPolys, v);
end for;

// calculate orbits
for p in polys do
    orbit := [];
    x := f(p);
    if x in uniqueOrbitPolys then
        Append(orbit, x);
        for M in G do
            orbitp := x*M;
            Append(~orbit, orbitp);
            i := Index(uniqueOrbitPolys, orbitp);
            if i ne 0 then
                Remove(~uniqueOrbitPolys, i);
            end if;
        end for;
        Append(~orbitsList, orbit);
    end if;
end for;

#orbitsList;