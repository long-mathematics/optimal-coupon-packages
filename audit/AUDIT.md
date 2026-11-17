# Mathematical audit

Audit date: 17 November 2026. Author: Christopher D. Long.

## Outcome

No gap was identified in the manuscript's analytic and combinatorial proof
chains. The principal theorems were checked from their hypotheses through
the equality arguments. This is an assisted mathematical audit, not
independent refereeing or proof-assistant certification.

## Proof checks

1. **Finner and the partition bound (Sections 1–2).** The independent
   coordinates are the finitely many package-type Poisson counts. A count
   occurs in at most `s` coverage indicators, making weights `1/s` valid.
   The included induction proof of Finner was checked, including unused
   Holder weight and monotone truncation. Strict concavity and the size
   budget give the bound for real `q=n/s`. Independent exponential waiting
   times preserve the expectation, including infinity.
2. **Equality (Section 2).** Equality of the integrated bound forces equality
   of the continuous CDFs, balanced marginals, and exact-size packages.
   Distinct incidence classes have union-arrival rate strictly above the
   single-class rate. Bonferroni makes the leading failure-tail coefficient
   equal to the integer number of classes. Matching it with `q` forces an
   integral `q`, exactly `q` classes of size `s`, and uniform class masses.
   This also handles `q=1` and proves nonattainment for nonintegral `q`.
3. **Discrete complementary optimum (Section 3).** The missing count is
   bounded by `n-s` after the first exact-size draw. Its expectation and
   Jensen give the geometric tail. Equality at two draws makes distinct
   complements disjoint; balanced omissions make them a uniform partition.
   Deterministic enlargement is a valid iid coupling. A smaller original
   package creates a positive-probability strict completion-time gap, so it
   cannot attain equality. The assumption `s<n` is essential and present.
4. **Repeated collection (Section 4).** The derivative calculation proves
   strict log-concavity of the Poisson upper tail for every positive integer
   quota. Simultaneous failure of distinct classes requires at most `2r-2`
   union arrivals and hence decays at a faster exponential rate. The same
   tail-coefficient rigidity proves uniqueness. No unequal-quota claim is
   inferred.
5. **Paired perturbations (Sections 5–6).** Pair swaps prove normalization
   and preservation of all inclusion marginals below degree `d`. Only
   omitted sets selecting one label from each pair contribute to the
   change. Sign pairing gives the stated even-power expansion. Each
   coefficient is an actual covering-family count. The capacity argument
   proves strict positivity of the quadratic coefficient at every feasible
   count, including `m=2`, `n=2d`, and endpoint amplitudes. Odd degree can
   create infinite completion at a degree-one endpoint; the text permits
   this. The Schilling specialization has preserved positive marginals and
   finite mean even at amplitude one.
6. **Hessians (Section 7).** Symmetry makes the uniform law stationary and
   the Hessians scalar on Johnson components. Pair cancellation places
   `chi_d` in `V_d`, with squared Euclidean norm `2^d binom(n-2d,s-d)`.
   Dividing directional curvature by the squared tangent norm gives the
   displayed eigenvalues and signs. The mean series converges absolutely
   by the geometric counting bound. Fixed-horizon local improvement is
   not promoted to one uniform perturbation amplitude for all horizons.
7. **Scope and counterexample (Section 8 and Appendix A).** Positive
   Poisson mixing is used only in the valid direction. Only the first
   feasible coefficient is extracted from the small-time bound. The
   three-pair law on three coupons has mean `5/2`, whereas a uniform
   minimum two-package cover has mean `3`. General fixed-draw partition
   order, arbitrary radial monotonicity, and variance minimization are not
   claimed.

## Changes in this revision

The core theorem statements and proof chains are retained. The manuscript
date, repository paths, and reproducibility appendix are updated. A new
verification program covers the previously untested Poisson chain and
common-quota means. The minimum-cover appendix is now self-contained:
its criticism of a specific external lemma was removed because the source
passage could not be freshly checked. See `CITATION_AUDIT.md`.

The later fixed-draw working notes are outside this release; their claims
have not been imported without a separate audit.

## Exact computation

All three suites pass with standard-library integer/rational arithmetic.

| Checks | Count |
| --- | ---: |
| Perturbed-law CDF identities, independent union-state DP | 1,890 |
| Preserved lower-order inclusion marginals | 939 |
| Covering coefficient nonnegativity/strictness, through `n=35` | 40,835 |
| Paired-ray comparisons | 1,350 |
| Geometric tails | 896 |
| Complementary-partition identities | 108 |
| Earliest-draw inequalities and attainment | 60 |
| Mean-Hessian coordinate identities | 1,262 |
| CDF-Hessian coordinate identities | 8,182 |
| Tangent trace/multiplicity identities | 28 |
| Exact Poisson/Finner/Jensen chains | 168 |
| Partition Poisson CDF attainment | 32 |
| Common-quota comparisons and partition identities | 33 |

Explicit rational examples and the four `(8,4)` Hessian eigenvalues also
match. The logs give the exact parameter ranges. These checks cannot
replace the universal proofs.

## Build

The 17-page PDF compiles with pdfLaTeX and latexmk. The final LaTeX log has
no unresolved citations/references, LaTeX/package warnings, or overfull
boxes. The environment initially lacked TeX's generated format/database
and font map; local runtime initialization fixed that environment issue
without changing the mathematical source or disabling microtype. The
repository build targets a normally installed TeX Live environment.

All 17 pages were rendered and visually inspected, including the formula-dense
Hessian sections and bibliography; no clipping or overlap was observed.
