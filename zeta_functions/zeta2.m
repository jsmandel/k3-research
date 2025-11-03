load "/Users/jordanmandel/Downloads/smoothList.m";
k := GF(2);
G := GL(3, k);
R<[x]> := PolynomialRing(k,3);
V, f := GModule(G, R, 3);
GV := ActionGroup(V);
PP<x,y,z> := PolynomialRing(k,3);
PF<w,x,y,z> := PolynomialRing(k,[3,1,1,1]);
PPtoPF := hom<PP -> PF | x, y, z>;
Q<u,v,w,a,b,c,x,y,z> := PolynomialRing(k,9);
f := hom< PP -> Q | x, y, z >;
NewQ<U,V,W,AA,BB,CC> := PolynomialRing(k,6); 
g := hom< Q -> NewQ | U,V,W,AA,BB,CC,0,0,0 >;
FQ<u,v,w,a,b,c,x,y,z> := FieldOfFractions(Q);


/* The following functions are for point counting. The first function calculates the orbits
of the Frobenius map in a field K. It returns a representative for each orbit, as well
as the sizes of the orbits (in two different arrays that have compatible indexing).) */

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
single triplet (x,y,z). */


PointCount := function(FieldCardinality,L,M)
        K := GF(FieldCardinality);
        Embed(k,K);
	S := {x^2 + x : x in K};
        FrobeniusReps, OrbitSizes := CalculateFrobeniusOrbit(K);
        BaseChangeSextic := ChangeRing(PF,K);
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



/* Computing Characteristic polynomials in characteristic 2. The input is a list of 12 numbers:
X(F_2^n) for n = 1,...,12 inside ``Uberlist'' below.  */

p := 2;

// Example curve
Uberlist := [ [ 7, 25, 73, 249, 1137, 4273, 16737, 65313, 264385, 1047745, 4203393, 16767105] ];

P<x> := PolynomialRing(Rationals());
C<i> := ComplexField(30);

for j in [1..#Uberlist] do
	FpPoints := Uberlist[j];
	Traces := [ Rationals()!(FpPoints[i] - p^(2*i) - 1) : i in [1..#FpPoints] ];

	CharPolyCoeffs := [ -Traces[1] ];
	for k in [2..#FpPoints] do
		ck := -(1/k)*(Traces[k] + &+[ CharPolyCoeffs[i]*Traces[k-i] : i in [1..k-1] ]);
		Append(~CharPolyCoeffs,ck);
	end for;

	ResetCharPolyCoeffs := CharPolyCoeffs;

	// check root of a poly are root of unity
	RightModulus := function(ModifiedPolynomial)
	        ChangeRing(P,C);
        	rts := Roots(ModifiedPolynomial,C);
        	time R := [ Abs(i) : i in rts[j], j in [1..#rts] ];
        	S := [ R[2*i - 1] : i in [1..#rts] ];
        	return S;
	end function;

	// decide the sign of the functional equation:

	if CharPolyCoeffs[12] eq p^2*CharPolyCoeffs[10] then
		// the sign is positive
		for k in [13..21] do
			ck := CharPolyCoeffs[22-k]*p^(2*k - 22);
			Append(~CharPolyCoeffs,ck);
		end for;
		Append(~CharPolyCoeffs,p^22);
		CharPolyCoeffs;
		CharPoly := x^22 + &+[ CharPolyCoeffs[i]*x^(22-i) : i in [1..22] ];
		CharPolyModified := p^(-22)*Evaluate(CharPoly,p*x);
		if CharPolyCoeffs[12] ne 0 then
			print "The factorization of the characteristic polynomial of surface #", j, "is :";
			Factorization(CharPolyModified);
			print "The sign of the functional equation is positive; c11 and c12 are", CharPolyCoeffs[11],CharPolyCoeffs[12];
			print "The absolute values of the roots are", RightModulus(CharPolyModified);
			print "----------------------------";
		else
		        print "For surface #", j, "we have c12 = 0, so can't determine functional equation from data given";
        		print "Assume the sign is positive. Then we obtain:";
			Factorization(CharPolyModified);
			print "The absolute values of the roots are", RightModulus(CharPolyModified);
			print "----------------------------";
		end if;
	end if;

	CharPolyCoeffs := ResetCharPolyCoeffs;
	if CharPolyCoeffs[12] eq -p^2*CharPolyCoeffs[10] then
	// the sign is negative
		for k in [13..21] do
        		ck := -CharPolyCoeffs[22-k]*p^(2*k - 22);
        		Append(~CharPolyCoeffs,ck);
       		 end for;
        	Append(~CharPolyCoeffs,-p^22);
		CharPolyCoeffs;
        	CharPoly := x^22 + &+[ CharPolyCoeffs[i]*x^(22-i) : i in [1..22] ]; 
        	CharPolyModified := p^(-22)*Evaluate(CharPoly,p*x);
        	if CharPolyCoeffs[12] ne 0 then
                	print "The factorization of the characteristic polynomial of surface #", j, "is :";
        	        Factorization(CharPolyModified);
        	        print "The sign of the functional equation is negative; c11 and c12 are", CharPolyCoeffs[11], CharPolyCoeffs[12];
        	        print "The absolute values of the roots are", RightModulus(CharPolyModified);
        	        print "----------------------------";
		else
   	             print "For surface #", j, "we have c12 = 0, so can't determine functional equation from data given";
   	             print "Assume the sign is negative. Then we obtain:";
    	            Factorization(CharPolyModified);
    	            print "The absolute values of the roots are", RightModulus(CharPolyModified);
     	           print "----------------------------";
     	   end if;
	end if;
end for;

for rep in smoothList do
	PointCount(2, rep[1], rep[2]);
end for;