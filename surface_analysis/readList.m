k := FiniteField(2);
G := GL(3, k);
R<[x]> := PolynomialRing(k,3);

binFile := "smoothOrbitList.bin";
f := Open(binFile, "rb");
smoothBin := ReadBytes(f, 5);
delete f;

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

smoothResult := [];
pos := 1;
block := smoothBin[pos..pos+4];
result := deserialize(block);
resulta := result[1]; resultb := result[2];
polya := &+[resulta[i]*MonomialsOfDegree(R,3)[i] : i in [1..10]];
polyb := &+[resultb[i]*MonomialsOfDegree(R,6)[i] : i in [1..28]];
polya; polyb;