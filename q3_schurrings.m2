-- ============================================================
-- Q3 CLEAN: GL_3 Schur decomposition via SchurRings package
--
-- SchurRings computes multiplicities of GL_n irreps in graded
-- modules correctly, without needing a GL_3-stable ideal.
-- This is the mathematically correct approach for Q3.
-- ============================================================

needsPackage "SchurRings"

print "============================================================"
print "Q3: GL_3 Schur decomposition via SchurRings"
print "n=3 | Corrected lambda via Excess(V_P over V_U)"
print "============================================================"
print ""

-- ── 0. Rings ──────────────────────────────────────────────────
kk = ZZ/32749
R9 = kk[w_0..w_8]

-- SchurRing for GL_3 (rank 3)
S = schurRing(QQ, s, 3)

-- ── 1. Polynomials ────────────────────────────────────────────
fU     = w_0*w_5*w_7 + w_1*w_3*w_8
fISr   = w_0*w_3*w_5 + w_1*w_2*w_4
fIS    = fISr^2
fClique = w_0^2*w_5^2*w_6^2
        + 4*w_0*w_1*w_3*w_5*w_7*w_8
        + 2*w_0*w_1*w_4*w_5*w_6*w_7
        + 2*w_0*w_2*w_3*w_5*w_6*w_8
        + 2*w_0*w_2*w_4*w_5*w_6^2
        + w_1^2*w_4^2*w_7^2
        + 2*w_1*w_2*w_3*w_4*w_7*w_8
        + 2*w_1*w_2*w_4^2*w_6*w_7
        + w_2^2*w_3^2*w_8^2
        + 2*w_2^2*w_3*w_4*w_6*w_8
        + w_2^2*w_4^2*w_6^2
-- ── 2. Quotient rings (Raw Ideals) ───────────────────────────
allPairs9 = flatten apply(toList(0..8), i -> apply(toList(0..8), j -> (i,j)))
buildLie = f -> select(apply(allPairs9, ij -> diff(R9_(ij#0), f) * R9_(ij#1)), g -> g != 0_R9)

-- Simplemente creamos el ideal raw original.
mkI_raw = f -> ideal({f} | buildLie f)

print "Building raw ideals..."
IU_raw = mkI_raw fU
IIS_raw = mkI_raw fIS
IClique_raw = mkI_raw fClique
print "  Done."
print ""

-- ── 3. Helpers & Linear Algebra Vector Spaces (Graded Gröbner) ──
gl3Pairs = flatten apply(toList(0..2), i ->
    apply(toList(0..2), j -> if i != j then (i,j) else null))
gl3Pairs = select(gl3Pairs, x -> x =!= null)

applyGL3gen = (eij, f) -> (
    i := eij#0; j := eij#1;
    sum apply(toList(0..2), k -> diff(R9_(3*j+k), f) * R9_(3*i+k)))

monoWeight = m -> (e := (exponents m)#0; {e#0+e#1+e#2, e#3+e#4+e#5, e#6+e#7+e#8})
pad3 = p -> (if #p == 0 then {0,0,0} else if #p == 1 then {p#0, 0, 0} else if #p == 2 then {p#0, p#1, 0} else {p#0, p#1, p#2})
polyDegree = p -> sum(first exponents leadMonomial p)

reduceGlobal = (polys) -> (
    pClean := select(polys, p -> p != 0);
    if #pClean == 0 then return {};
    (M, C) := coefficients(matrix {pClean});
    basisC := mingens image C;
    flatten entries (M * basisC)
)

-- Variable global para heredar grados. Se vaciará al cambiar de variedad.
globalCache = new MutableHashTable

getGradedStableComponent = (Iraw, d) -> (
    if globalCache#?d then return globalCache#d;
    
    inherited := {};
    if d > 3 then (
        -- Magia de Gröbner Truncado: Heredar relaciones del grado d-1
        prevBasis := getGradedStableComponent(Iraw, d-1);
        varsR1 := flatten entries basis(1, R9);
        inherited = flatten apply(prevBasis, p -> apply(varsR1, v -> v * p));
    );
    
    rawGens := select(flatten entries gens Iraw, g -> polyDegree g == d);
    IdBasis := reduceGlobal(inherited | rawGens);
    
    changed := true;
    while changed do (
        changed = false;
        newGens := flatten apply(gl3Pairs, eij ->
            apply(IdBasis, p -> applyGL3gen(eij, p)));
            
        oldDim := length IdBasis;
        IdBasis = reduceGlobal(IdBasis | newGens);
        
        if length IdBasis > oldDim then changed = true;
    );
    
    globalCache#d = IdBasis;
    return IdBasis;
)

-- ── 4. Schur Processing Helpers ────────────────────────────────
lexGT = (p, q) -> (
    n := min(#p, #q); found := false; result := false; i := 0;
    while i < n and not found do (
        if p#i > q#i then (found = true; result = true)
        else if p#i < q#i then (found = true; result = false);
        i = i+1);
    result)

sortDesc = lst -> (
    arr := new MutableList from lst;
    apply(toList(1..#arr-1), i -> (
        key := arr#i; j := i-1;
        while j >= 0 and lexGT(key, arr#j) do (arr#(j+1) = arr#j; j = j-1);
        arr#(j+1) = key));
    toList arr)

Rxyz = QQ[x,y,z]
schurPoly3 = mu -> (
    mu1 := if #mu>=1 then mu#0 else 0;
    mu2 := if #mu>=2 then mu#1 else 0;
    mu3 := if #mu>=3 then mu#2 else 0;
    l1 := mu1+2; l2 := mu2+1; l3 := mu3;
    numer := x^l1*y^l2*z^l3 - x^l1*y^l3*z^l2 - x^l2*y^l1*z^l3 + x^l2*y^l3*z^l1 + x^l3*y^l1*z^l2 - x^l3*y^l2*z^l1;
    denom := (x-y)*(x-z)*(y-z);
    numer // denom)

weightMult = (mu, a, b, c) -> (
    sp := schurPoly3 mu;
    cf := coefficient(x^a*y^b*z^c, sp);
    if cf === null then 0 else lift(cf, ZZ))

-- ── 5. Decompose by subtracting characters ──────────────────
decomposeDegree = (Iraw, d) -> (
    IdBasis := getGradedStableComponent(Iraw, d);
    R9BasisD := flatten entries basis(d, R9);
    
    wtTbl := new MutableHashTable;
    
    -- Sumamos la contribución del anillo entero
    scan(R9BasisD, m -> (
        w := pad3(monoWeight m);
        wtTbl#w = if wtTbl#?w then wtTbl#w+1 else 1;
    ));
    
    -- Restamos la contribución de los vectores en el ideal
    scan(IdBasis, p -> (
        w := pad3(monoWeight(leadMonomial p));
        wtTbl#w = wtTbl#w - 1;
    ));
        
    parts := sortDesc apply(select(partitions d, p -> #p <= 3), p -> toList p);
    charElt := 0_S;
    
    scan(parts, p -> (
        hw := pad3 p;
        k  := if wtTbl#?hw then wtTbl#hw else 0;
        if k > 0 then (
            charElt = charElt + k * s_p;
            scan(toList(0..d), a ->
                scan(toList(0..d-a), b -> (
                    c  := d-a-b;
                    wt := {a,b,c};
                    km := weightMult(p, a, b, c);
                    if km > 0 then wtTbl#wt = (if wtTbl#?wt then wtTbl#wt else 0) - k*km
                )))
        )
    ));
    
    dimQ := (#R9BasisD) - (#IdBasis);
    {charElt, wtTbl, dimQ}
)

-- ── 6. Full decomposition ─────────────────────────────────────
DCOMP = 15

print "============================================================"
print ("Full GL_3 decomposition d=3.." | toString DCOMP)
print "============================================================"
print ""

joinStr = (sep, lst) -> (
    if #lst == 0 then ""
    else if #lst == 1 then lst#0
    else (s := lst#0; scan(toList(1..#lst-1), i -> s = s|sep|lst#i); s))

weylDim = mu -> (
    m1 := if #mu>=1 then mu#0 else 0;
    m2 := if #mu>=2 then mu#1 else 0;
    m3 := if #mu>=3 then mu#2 else 0;
    (m1-m2+1)*(m1-m3+2)*(m2-m3+1)//2)

-- Extract k_mu from SchurRing element
extractMults = chiElt -> (
    -- chiElt is an element of S = schurRing(QQ, s, 3)
    -- Its monomials are s_{mu} with coefficients k_mu
    result := new HashTable from apply(listForm chiElt, t -> (
        mu := toList(t#0);  -- partition as list
        k  := t#1;          -- multiplicity
        (mu, lift(k, ZZ))));
    result)

processVariety = (Iraw, label) -> (
    print("--- " | label | " ---");
    globalCache = new MutableHashTable;
    rows := apply(toList(3..DCOMP), d -> (
        res := decomposeDegree(Iraw, d);
        chi := res#0;
        ht  := extractMults chi;
        
        act := res#2; 
        
        cmp := sum apply(keys ht, mu -> (ht#mu) * weylDim(toList mu));
        ok  := cmp == act;
        ks  := sort keys ht;
        es  := apply(ks, k -> toString(toList k) | "->" | toString(ht#k));
        print("  d=" | toString d |
              "  dim=" | toString act |
              "  computed=" | toString cmp |
              "  " | if ok then "OK" else "FAIL" |
              "  {" | joinStr(", ", es) | "}");
        {d, ht}));
    rows)

resU = processVariety(IU_raw,      "V_U (SAT)")
print ""
resIS = processVariety(IIS_raw,     "V_IS (Ind.Set)")
print ""
resC  = processVariety(IClique_raw, "V_Clique")

-- ── 7. Excess and lambda ──────────────────────────────────────
print ""
print "============================================================"
print "Excess analysis: Excess(V_P over V_U)"
print "rho(P,d) = ExcessMult / TotalMult"
print "lambda(P) = limsup rho(P,d)"
print "============================================================"
print ""

computeExcess = (resP, resU) -> (
    apply(#resP, i -> (
        d   := (resP#i)#0;
        htP := (resP#i)#1;
        htU := (resU#i)#1;
        allMu  := unique(keys htP | keys htU);
        excess := select(allMu, mu ->
            (if htP#?mu then htP#mu else 0) > 0 and
            (if htU#?mu then htU#mu else 0) == 0);
        exMult  := sum(apply(excess, mu -> if htP#?mu then htP#mu else 0), x->x);
        totMult := sum apply(keys htP, mu -> htP#mu);
        rho     := if totMult > 0 then exMult/totMult else 0;
        {d, excess, exMult, totMult, rho})))

printExcess = (label, lst) -> (
    print("--- " | label | " ---");
    scan(lst, row -> (
        d   := row#0; exc := row#1; eM := row#2; tM := row#3; rho := row#4;
        excStr := if #exc==0 then "{}" else "{"|joinStr(", ",apply(exc,mu->toString mu))|"}";
        print("  d=" | toString d |
              "  excess=" | toString(#exc) | excStr |
              "  exMult=" | toString eM |
              "  total=" | toString tM |
              "  rho=" | toString rho))))

exIS = computeExcess(resIS, resU)
exC  = computeExcess(resC,  resU)

printExcess("V_IS over V_U", exIS)
print ""
printExcess("V_Clique over V_U", exC)

print ""
print "============================================================"
print "Lambda estimates"
print "============================================================"
print ""

totExIS = sum apply(exIS, r -> r#2)
totExC  = sum apply(exC,  r -> r#2)
totMIS  = sum apply(exIS, r -> r#3)
totMC   = sum apply(exC,  r -> r#3)
rhoIS   = if totMIS > 0 then totExIS/totMIS else 0
rhoC    = if totMC  > 0 then totExC/totMC   else 0

print ("lambda(SAT)    = 0  [reference]")
print ("rho(IS)        = " | toString rhoIS | "  =>  lambda(IS) ~ " | toString(rhoIS))
print ("rho(Clique)    = " | toString rhoC  | "  =>  lambda(Clique) ~ " | toString(rhoC))
print ""
print "d  | rho(IS,d) | rho(C,d)"
print "---|-----------|----------"
apply(#exIS, i -> (
    d  := (exIS#i)#0;
    rI := (exIS#i)#4;
    rC := (exC#i)#4;
    print(toString d | "  | " | toString rI | "  | " | toString rC)))

print ""
print "============================================================"
print "Done."
print "============================================================"
