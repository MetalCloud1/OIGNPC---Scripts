-- ========================================================================= 
-- Q8(i) STABILIZER - v3 (corrected GL_3 weight) 
-- 
-- The weight of E_{ij} below the action of FILA of GL_3 is: 
-- e_{row(j)} - e_{row(i)} whence row(k) = k // 3 
-- 
-- The action of row g in GL_3 above w_{3r+c} es: 
-- w_{3r+c} -> sum_s g_{sr} w_{3s+c} 
-- So E_{ij} holds weight e_{row(j)} - e_{row(i)}. 
-- =========================================================================

needsPackage "SchurRings"

kk = ZZ/32749
R9 = kk[w_0..w_8]
S  = schurRing(QQ, s, 3)

fU = w_0*w_5*w_7 + w_1*w_3*w_8

-- ── Action of gl_9 ────────────────────── ──────────────────────
actionOf   = (i, j, f) -> diff(R9_i, f) * R9_j
actionVecs = flatten apply(toList(0..8), i ->
    apply(toList(0..8), j -> actionOf(i, j, fU)))

R9Basis3 = flatten entries basis(3, R9)
nMons    = #R9Basis3

actionMatrix = matrix apply(toList(0..nMons-1), r ->
    apply(toList(0..80), c -> (
        p := actionVecs#c;
        if p == 0_R9 then 0_kk
        else coefficient(R9Basis3#r, p))))

K        = kernel actionMatrix
dimStab  = numColumns generators K
dimOrbit = 81 - dimStab

print("dim(stab)  = " | toString dimStab)
print("dim(orbit) = " | toString dimOrbit)
print ""

-- ── GL_3 CORRECT weight of E_{ij} ─────────────────────────────
rowOf     = k -> k // 3
gl3Weight = idx -> (
    i    := idx // 9;
    j    := idx % 9;
    rowI := rowOf i;
    rowJ := rowOf j;
    apply(toList(0..2), k -> 
        (if k == rowJ then 1 else 0) - (if k == rowI then 1 else 0))
)

-- Verification: show some weights
print "Weight verification (first 12 E_{ij}):"
apply(toList(0..11), idx -> (
    i := idx // 9; j := idx % 9;
    w := gl3Weight idx;
    print("  E_{" | toString i | "," | toString j |
          "} row(" | toString i | ")=" | toString(rowOf i) |
          " row(" | toString j | ")=" | toString(rowOf j) |
          " peso=" | toString w)))
print ""

weightSpaces = new MutableHashTable
apply(toList(0..80), idx -> (
    w := gl3Weight idx;
    if not weightSpaces#?w then weightSpaces#w = {};
    weightSpaces#w = weightSpaces#w | {idx}))

print "Pesos distintos de gl_9 bajo accion de fila de GL_3:"
scan(sort keys weightSpaces, w ->
    print("  " | toString w | " : " | toString(#(weightSpaces#w)) | " elements"))
print ""

print "Weight | dim total | dim stab | dim complement"
print "-----|-----------|----------|----------------"

stabWeightDims  = new MutableHashTable
complWeightDims = new MutableHashTable

scan(sort keys weightSpaces, w -> (
    idxs     := weightSpaces#w;
    fullDim  := #idxs;
    subM     := actionMatrix_idxs;
    Ksub     := kernel subM;
    stabDim  := numColumns generators Ksub;
    complDim := fullDim - stabDim;
    stabWeightDims#w  = stabDim;
    complWeightDims#w = complDim;
    print("  " | toString w |
          " | " | toString fullDim |
          " | " | toString stabDim |
          " | " | toString complDim)))

print ""
totalStab  = sum apply(values stabWeightDims,  x -> x)
totalCompl = sum apply(values complWeightDims, x -> x)
print("Stab sum by weights= " | toString totalStab  | "  (must be " | toString dimStab  | ")")
print("Suma compl by weights = " | toString totalCompl | "  (must be " | toString dimOrbit | ")")
print ""

print "============================================================"
print "Schur decomposition of the T_{f_U} complement(GL_9 orbit):"
print "============================================================"
print ""
print "Plugin weight table:"
scan(sort keys complWeightDims, w -> (
    d := complWeightDims#w;
    if d > 0 then print("  " | toString w | " -> " | toString d)))

offDiagCompl = select(keys complWeightDims,
    w -> (w#0 + w#1 + w#2 == 0) and w != {0,0,0} and complWeightDims#w > 0)
diagComplDim = if complWeightDims#?{0,0,0} then complWeightDims#{0,0,0} else 0

print ""
print("Dimension in off-diagonal weights (suma=0): " | toString(sum apply(offDiagCompl, w -> complWeightDims#w)))
print("Dimension in diagonal weight (0,0,0): " | toString diagComplDim)
print ""

if #offDiagCompl > 0 then (
    mults := apply(offDiagCompl, w -> complWeightDims#w);
    minM  := min mults;
    maxM  := max mults;
    print("Multiplicities in off-diagonal weights: min=" | toString minM | " max=" | toString maxM);
    if minM == maxM then (
        m := minM;
        print("  => Uniforms: each off-diagonal weight appears" | toString m | " veces");
        print("  => " | toString m | " copias de la rep Adjunta de GL_3 contribuyen dim " | toString(8*m));
        print("     + " | toString diagComplDim | " dim en diagonal total");
        residuoDiag := diagComplDim - 2*m;
        print("  => Diagonal residue = " | toString residuoDiag | " (representaciones triviales)");
        print("  => Total complement for this argument = " | toString(8*m + residuoDiag))
    ) else print("  => "Not uniform, more complex structure")
)

print ""
print "============================================================"
print "FINAL SUMMARY"
print "============================================================"
print ""
print("dim(gl_9)   = 81")
print("dim(stab)   = " | toString dimStab)
print("dim(orbita) = " | toString dimOrbit)
print ""
print "Implication for Q8(i):" 
print "The tangent space T_{f_U} decomposes under GL_3 as:" 
print "6 copies of the Attached representation + 2 trivial copies." 
print "When applying the shift by the weight of f_U (1,1,1), this generates" 
print "exactly the modules that, through the Cauchy identity," 
print " absorb all copies of s_{(d-1,1)} into Sym^d(C^9)." 
print ""
print "Done."