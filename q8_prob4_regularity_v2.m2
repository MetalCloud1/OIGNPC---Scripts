-- ============================================================
-- q8_prob4_regularity_v2.m2
-- Problem 4: reg(I(V_U^(n))) = O(n^c)?
--
-- KEY INSIGHT from q9_cylindrical_embedding.m2:
-- The cylindrical embedding is NOT fU^(n) as an nxn permanent.
-- It is fU^(3) embedded in larger GL_{m^2} by adding
-- inert auxiliary variables. fCyl_m = w_0*w_5*w_7 + w_1*w_3*w_8
-- living in R_{m^2} with w_9..w_{m^2-1} as auxiliary variables.
--
-- So I(V_U^{cyl,(m)}) is the GL_{m^2}-stable ideal of fU3
-- viewed in the ambient space C^{m^2}.
--
-- Problem 4 asks: does reg(I(V_U^{cyl,(m)})) = O(1) in m?
-- The paper conjectures reg = 3 for all m >= 3.
--
-- This script computes reg for m=3,4,5,6 using the same
-- linear algebra approach as q8_regularity.m2.
--
-- Optimization: rankInDeg uses coefficient matrix rank
-- (pure linear algebra, no Groebner) for the convergence loop.
-- The final regularity call uses the ideal directly.
-- ============================================================

kk = ZZ/32749

-- Fast rank check via coefficient matrix (no Groebner)
rankInDeg = (gens, Bd) -> (
    pClean := select(gens, g -> g != 0);
    if #pClean == 0 then return 0;
    C := (coefficients(matrix {pClean}, Monomials => Bd))#1;
    rank C)

-- Build I(V_U^{cyl,(m)}): GL_{m^2}-stable closure of fU3 in R_{m^2}
-- fU3 is fixed = w_0*w_5*w_7 + w_1*w_3*w_8, degree 3
-- Row action of GL_m on C^m otimes C^m: e_{ij} sends w_(m*i+c) -> w_(m*j+c)
buildCylIdeal = (m, Rm, fCyl) -> (
    mm := m*m;
    rowAct := (ii, jj, f) -> sum(toList(0..m-1), c ->
        diff(Rm_(m*ii+c), f) * Rm_(m*jj+c));
    rowPairs := flatten apply(toList(0..m-1), ii ->
        apply(select(toList(0..m-1), jj -> jj != ii), jj -> {ii,jj}));
    -- gl_{m^2} orbit of fCyl (all diff * var pairs)
    orbitGens := select(
        flatten apply(toList(0..mm-1), i ->
            apply(toList(0..mm-1), j -> diff(Rm_i, fCyl) * Rm_j)),
        g -> g != 0_Rm);
    print("  gl_" | toString mm | " orbit: " | toString(#orbitGens) | " gens");
    -- Close under GL_m row action in degree 3
    Bd := basis(3, Rm);
    currentGens := orbitGens;
    changed := true; iter := 0;
    while changed do (
        iter = iter + 1;
        newG := flatten apply(currentGens, f ->
            apply(rowPairs, ij -> rowAct(ij#0, ij#1, f)));
        newG = select(newG, g -> g != 0_Rm);
        dimOld := rankInDeg(currentGens, Bd);
        dimNew := rankInDeg(currentGens | newG, Bd);
        print("  iter " | toString iter | ": dim=" | toString dimNew);
        if dimNew > dimOld
            then (currentGens = currentGens | newG; changed = true)
            else changed = false;);
    ideal currentGens)

-- ============================================================
-- m=3 (baseline: reproduces q8_regularity.m2 exactly)
-- ============================================================
print "=== m=3 (baseline) ==="
R9  = kk[w_0..w_8]
fU3 = w_0*w_5*w_7 + w_1*w_3*w_8
IU3 = buildCylIdeal(3, R9, fU3)
print("hilbertFunction(3, R9/IU3) = " | toString(hilbertFunction(3, R9^1/IU3)) | "  (expected 30)")
print "Computing reg..."
reg3 = regularity IU3
print("reg(m=3) = " | toString reg3 | "  (expected 3)")
print ""

-- ============================================================
-- m=4: fU3 embedded in GL_16
-- fCyl4 = w_0*w_5*w_7 + w_1*w_3*w_8 in R16
-- dim(Stab) = 31 + (16-9)^2 = 80 (verified in q9_cylindrical_embedding.m2)
-- ============================================================
print "=== m=4 (fU3 in GL_16) ==="
R16   = kk[w_0..w_15]
fCyl4 = w_0*w_5*w_7 + w_1*w_3*w_8
IU4   = buildCylIdeal(4, R16, fCyl4)
print("hilbertFunction(3, R16/IU4) = " | toString(hilbertFunction(3, R16^1/IU4)))
print "Computing reg..."
reg4 = regularity IU4
print("reg(m=4) = " | toString reg4 | "  (expected 3)")
print ""

-- ============================================================
-- m=5: fU3 embedded in GL_25
-- dim(Stab) = 31 + (25-9)^2 = 287
-- ============================================================
print "=== m=5 (fU3 in GL_25) ==="
R25   = kk[w_0..w_24]
fCyl5 = w_0*w_5*w_7 + w_1*w_3*w_8
IU5   = buildCylIdeal(5, R25, fCyl5)
print("hilbertFunction(3, R25/IU5) = " | toString(hilbertFunction(3, R25^1/IU5)))
print "Computing reg..."
reg5 = regularity IU5
print("reg(m=5) = " | toString reg5 | "  (expected 3)")
print ""

-- ============================================================
-- m=6: fU3 embedded in GL_36
-- dim(Stab) = 31 + (36-9)^2 = 760
-- ============================================================
print "=== m=6 (fU3 in GL_36) ==="
R36   = kk[w_0..w_35]
fCyl6 = w_0*w_5*w_7 + w_1*w_3*w_8
IU6   = buildCylIdeal(6, R36, fCyl6)
print("hilbertFunction(3, R36/IU6) = " | toString(hilbertFunction(3, R36^1/IU6)))
print "Computing reg..."
reg6 = regularity IU6
print("reg(m=6) = " | toString reg6 | "  (expected 3)")
print ""

-- ============================================================
-- SUMMARY
-- ============================================================
print "============================================================"
print "Problem 4: reg(I(V_U^{cyl,(m)})) = O(m^c)?"
print "============================================================"
print("m=3: reg=" | toString reg3)
print("m=4: reg=" | toString reg4)
print("m=5: reg=" | toString reg5)
print("m=6: reg=" | toString reg6)
print ""
if reg3 == 3 and reg4 == 3 and reg5 == 3 and reg6 == 3 then (
    print "RESULT: reg = 3 for m=3,4,5,6.";
    print "Strong evidence: reg(I(V_U^{cyl,(m)})) = 3 = O(1) for all m.";
    print "This is stronger than O(m^c) -- confirms Problem 4 conjecture.")
else print("RESULT: regs = " | toString reg3 | "," | toString reg4 | "," | toString reg5 | "," | toString reg6)
print "============================================================"
print "q8_prob4_regularity_v2.m2 COMPLETE"
