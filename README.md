# OIGNPC: Scripts
# Computational Evidence for the Trivial-Multiplicity Obstruction

## Core Scripts (referenced in the paper)

q8_regularity.m2
  Verifies reg(I_U^{GL3}) = 3 over Z/32749Z.
  Expected output: regularity = 3
  Referenced in: Theorem 4.1 (Rigidity Theorem), Step 2

q8_stabilizer_v3.m2
  Computes dim(Stab_{GL9}(f_U)) = 31 and the tangent space
  decomposition s_{(2,1)}^6 ⊕ s_{(1,1,1)}^2.
  Expected output: stabilizer dimension 31, multiplicities (6,2)
  Referenced in: Theorem 3.2 (Exact Stabilizer Formula)

q3_schurrings.m2
  Verifies H(d) = 3·C(d+2,2) for d = 3,...,10.
  Expected output: Hilbert function values 30, 45, 63, 84, 108, 135, 165, 198
  Referenced in: Theorem 4.1 (Rigidity Theorem), Step 3b

q8_prob4_regularity_v2.m2
  Verifies reg(I_m) = 3 for m = 3, 4, 5.
  Expected output: regularity 3 in each case
  Referenced in: Theorem 5.2 (Regularity Stabilization)

## Environment
Macaulay2 version: 1.22
All computations run over Z/32749Z (prime chosen to avoid
characteristic-dependent behavior; coefficients of all generators
lie in {0,1}, so results are exact over this field).

## Auxiliary Scripts (exploratory, grouped by research question)

### Q1 — VARIETY NOVELTY & CANONICAL CONSTRUCTION
q1_variety_novelty.m2

### Q2 — ALGORITHM INDEPENDENCE (CIAP+ Invariance)
q2_algorithm_independence.m2

### Q3 — Schur ring structure and Hilbert functions
q3_kronecker_lambda_v8.m2
q3_schurrings.m2

### Q7 — Plethystic Filter
q7_plethystic_filter.m2

### Q8 — Regularity and resolution
q8_betti_degrees.m2
q8_capacity_check.m2
q8_clique_stabilizer.m2
q8_free_resolution.m2
q8_ideal_character.m2
q8_LR_verification.m2
q8_regularity_crosscheck.m2

### Q9 — Cylindrical embedding and padding
q9_cylindrical_embedding.m2
