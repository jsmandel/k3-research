k := FiniteField(2);
G := GL(3, k);
R<[x]> := PolynomialRing(k,3);

binFile := "/Users/jordanmandel/Desktop/allSmoothList.bin";
f := Open(binFile, "rb");
smoothList := [];

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

// read binary file, 5 bytes at a time, and reconstruct smooth list
while true do
    block := ReadBytes(f, 5);  // try to read exactly 5 bytes at a time
    if #block eq 0 then        // EOF reached
        break;
    end if;

    // discard any final incomplete block (this shouldn't be possible if binary file correctly coded, but just in case)
    if #block lt 5 then
        break;
    end if;

    pair := [];
    result := deserialize(block);
    resulta := result[1]; resultb := result[2];
    polya := &+[resulta[i]*MonomialsOfDegree(R,3)[i] : i in [1..10]];
    polyb := &+[resultb[i]*MonomialsOfDegree(R,6)[i] : i in [1..28]];
    Append(~pair, polya); Append(~pair, polyb); 
    Append(~smoothList, pair);

end while;

delete f;