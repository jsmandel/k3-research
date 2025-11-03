k := GF(2);
G := GL(3, k);
R<[x]> := PolynomialRing(k,3);
V, f := GModule(G, R, 3);
GV := ActionGroup(V);
AttachSpec("../orbit_calc/k3.spec");
load "../surface_analysis/smooth.m";

/* The following functions are for point counting. The first function calculates the orbits
of the Frobenius map in a field K. It returns a representative for each orbit, as well
as the sizes of the orbits (in two different arrays that have compatible indexing).) 
From file */
CalculateFrobeniusOrbit := function(K)
	orbits := [];
	for a in K do
		orbit := [a];
		for r in [1..#K] do
			if Frobenius(a,r) ne a then
				Append(~orbit,Frobenius(a,r));
			else
				break;
			end if;
		end for;
	Append(~orbits,Seqset(orbit));
	end for;
	orbits := Setseq(Seqset(orbits));
	orbits := [ Setseq(orbits[i]) : i in [1..#orbits] ];
	FrobeniusReps := [ orbits[i][1] : i in [1..#orbits]];
	OrbitSizes := [ #orbits[i] : i in [1..#orbits]];
	return FrobeniusReps, OrbitSizes;
end function;

/* The following function counts points on a K3 surface of the form 

w^2 + L*w + M = 0    (*)

over a finite field of characteristic 2.  
The idea is as follows:
(1) Set x = 1. Look at y and z in K and evaluate L(1,y,z). If L is zero then you
have two points if M(1,y,z) is a square in K and zero otherwise.
(2) If L(1,y,z) is not zero then (*) has a solution iff

t^2 + t + M/L^2

has a solution.  This is an Artin-Schreier polynomial.  It has 2 solutions if
M/L^2 is in S := {x^2 + x : x in K}
and no solutions otherwise. These first two steps are sped up as follows:  
(2')    Since (*) is defined over the prime subfield, we
	may evaluate the sextic at y = ``representative of a Frobenius orbit'' and
	z arbitrary and add the orbit size of y every time we hit a solution.
(3) Set x = 0. Take y = 1. Run steps (1) and (2) as z ranges over K.
(4) Set x = 0 and y = 0.  Then we must have z = 1. Run step (1) and (2) for this
single triplet (x,y,z). 
From file */


PointCount := function(FieldCardinality,L,M)
        K := GF(FieldCardinality);
        Embed(k,K);
	S := {x^2 + x : x in K};
        FrobeniusReps, OrbitSizes := CalculateFrobeniusOrbit(K);
        // BaseChangeSextic := ChangeRing(PF,K);
	r1 := PolynomialRing(K);
	// first look at patch where x = 1
	s1 := 0;
	for i in [1..#FrobeniusReps] do
		y := FrobeniusReps[i];
		Osize := OrbitSizes[i];
		NewPolyL := Evaluate(L,[1,y,r1.1]);
		NewPolyM := Evaluate(M,[1,y,r1.1]);
		for z in K do
			if Evaluate(NewPolyL,z) eq 0 then
				if IsSquare(Evaluate(NewPolyM,z)) then
					s1 := s1 + Osize;
				end if;
			/* The equation w^2 + w*L(1,y,z) + M(1,y,z) has
			two solutions if M/L^2 is in S and zero solutions otherwise */
			elif  Evaluate(NewPolyM,z)/Evaluate(NewPolyL^2,z) in S then
				s1 := s1 + 2*Osize;
			end if;
		end for;
	end for;
	// Now look at patch where x = 0.
	// First take y = 1.
	s2 := 0;
	NewPolyL := Evaluate(L,[0,1,r1.1]);
	NewPolyM := Evaluate(M,[0,1,r1.1]);
	for z in K do
		if Evaluate(NewPolyL,z) eq 0 then
			if IsSquare(Evaluate(NewPolyM,z)) then		
				s2 := s2 + 1;
			end if;
		elif  Evaluate(NewPolyM,z)/Evaluate(NewPolyL^2,z) in S then
			s2 := s2 + 2;
		end if;
	end for;
	// Finally look at x = 0, y = 0, z = 1 (note that x = y = z = 0 can't be a point on the surface)
	if Evaluate(L,[0,0,1]) eq 0 then
		if IsSquare(Evaluate(M,[0,0,1])) then
			s3 := 1;
		end if;
	elif  Evaluate(M,[0,0,1])/Evaluate(L^2,[0,0,1]) in S then
		s3 := 2;
	else
		s3 := 0;
	end if;
        return (s1 + s2 + s3);
end function;

for rep in orbitRepList do
	PointCount(2, rep[1], rep[2]);
end for;