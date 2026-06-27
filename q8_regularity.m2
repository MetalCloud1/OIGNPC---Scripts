kk = ZZ/32749
R9 = kk[w_0..w_8]
fU = w_0*w_5*w_7 + w_1*w_3*w_8

rowAct = (i, j, f) -> sum(toList(0..2), k ->
    diff(w_(3*i+k), f) * w_(3*j+k))

allPairs9 = flatten apply(toList(0..8), i ->
    apply(toList(0..8), j -> (i,j)))
orbitGens = select(
    apply(allPairs9, ij -> diff(w_(ij#0), fU) * w_(ij#1)),
    g -> g != 0_R9)

rowPairs = {{0,1},{1,0},{0,2},{2,0},{1,2},{2,1}}

-- Work with the module in degree 3 as a vector space 
-- Base of Sym^3(R9) in grade 3
B = basis(3, R9)
-- Current space: columns of a matrix over kk
V = (map(kk^(numColumns B), kk^(#orbitGens), 
    matrix apply(toList(0..numColumns B - 1), r ->
        apply(orbitGens, g -> 
            coefficient(B_(0,r), g)))))

-- Iterate: apply rowAct to current generators and expand
currentGens = orbitGens
changed = true
iter = 0
while changed do (
    iter = iter + 1;
    newGens := flatten apply(currentGens, f ->
        apply(rowPairs, ij -> rowAct(ij#0, ij#1, f)));
    newGens = select(newGens, g -> g != 0_R9);
    Inew := ideal(currentGens | newGens);
    dimNew := numColumns basis(3, Inew);
    dimOld := numColumns basis(3, ideal currentGens);
    print("Iter " | toString iter | ": dim=" | toString dimNew);
    if dimNew > dimOld then (
        currentGens = currentGens | newGens;
        changed = true)
    else
        changed = false)

IU_GL3 = ideal currentGens
print("final dim in grade 3: " | toString(numColumns basis(3, IU_GL3)))
print("hilbertFunction(3, R9/IU_GL3) = " | 
      toString(hilbertFunction(3, R9^1/IU_GL3)))
print("Expected: 135 ideal, 30 quotient") 
print "Calculating regularity..."
print("reg(IU_GL3) = " | toString(regularity IU_GL3))