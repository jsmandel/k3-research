testList :=[];

for i in [1..10] do
    pair := [i, i];
    Append(~testList, pair);
end for;


f := Open("testList.m", "w");
Puts(f, "testList := " cat Sprint(testList, "Magma") cat ";");
Flush(f);