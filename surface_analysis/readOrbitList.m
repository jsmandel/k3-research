k := FiniteField(2);
G := GL(3, k);
R<[x]> := PolynomialRing(k,3);

binFile := "smoothOrbitList.bin";
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
    polya; polyb;
    Append(~pair, polya); Append(~pair, polyb); 
    Append(~smoothList, pair);

end while;

delete f;

// // WORKING CASE FOR 1 POLYNOMIAL PAIR
// smoothResult := [];
// pos := 1;
// block := smoothBin[pos..pos+4];
// result := deserialize(block);
// resulta := result[1]; resultb := result[2];
// polya := &+[resulta[i]*MonomialsOfDegree(R,3)[i] : i in [1..10]];
// polyb := &+[resultb[i]*MonomialsOfDegree(R,6)[i] : i in [1..28]];
// polya; polyb;


// WORKING CASE FOR SPECIFIC NUMBER OF BLOCKS
// smoothBin := ReadBytes(f, 10000);
// nBlocks := #smoothBin div 5;  // determine number of blocks

// for idx in [1..nBlocks] do
//     pair := [];
//     block := smoothBin[(5*idx-4)..(5*idx)]; // each 5-byte sequence
//     result := deserialize(block);
//     resulta := result[1]; resultb := result[2];
//     polya := &+[resulta[i]*MonomialsOfDegree(R,3)[i] : i in [1..10]];
//     polyb := &+[resultb[i]*MonomialsOfDegree(R,6)[i] : i in [1..28]];
//     polya; polyb;
//     Append(~pair, polya); Append(~pair, polyb); 
//     Append(~smoothList, pair);
// end for;