load "orbitList2.m";
// #orbitsList;
totalSets := 0;

for orbit in orbitsList do
    totalSets := totalSets + #orbit;
end for;