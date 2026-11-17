# Optimal coupon packages

**Christopher D. Long**

Research manuscript on optimal iid coupon packages, the Fang–Chang
(Chang–Fang) conjecture, and a marginal-preserving stochastic strengthening
of Schilling's conjecture.

[Read the paper](paper/optimal_coupon_packages.pdf) ·
[Download PDF](https://github.com/long-mathematics/optimal-coupon-packages/raw/refs/heads/main/paper/optimal_coupon_packages.pdf) ·
[LaTeX source](paper/optimal_coupon_packages.tex) ·
[Mathematical audit](audit/AUDIT.md) · [Reference audit](audit/CITATION_AUDIT.md)

## Principal results

Let packages be independent draws from a fixed law on subsets of an
`n`-coupon universe, each of size at most `s`. Put `q = n/s`.

| Result | Scope and conclusion |
| --- | --- |
| Fang–Chang expected-time optimum | If `s` divides `n`, the uniform law on a partition into `s`-sets uniquely minimizes the mean, with value `q H_q`, up to relabeling. |
| Poissonized stochastic optimum | Under rate-one package arrivals, the completion CDF is at most `(1-exp(-t/q))^q`. The bound holds also for nonintegral `q`; equality in its mean bound requires an integral `q` and a partition law. |
| Repeated collection | For every common positive integer quota, the partition law is again the unique mean optimizer and is optimal in Poissonized stochastic order when `s` divides `n`. |
| Schilling strengthening | For `2 <= s <= n-2`, an explicit four-coupon perturbation preserves every singleton marginal and strictly improves the uniform law's completion CDF at every feasible draw count. |
| Complementary partitions | If `n-s` divides `n` and `s<n`, the uniform complementary-partition law is optimal in discrete stochastic order, with mean `1+n/s`. |
| Hessian signature | At the uniform law, the mean Hessian has sign `(-1)^(d+1)` on Johnson component `V_d`; the feasible fixed-draw CDF Hessian has the opposite sign. |

The partition proof combines independent Poisson package counts, Finner's
inequality, strict concavity, and an incidence-class tail argument for
uniqueness. The Schilling proof uses signed paired perturbations whose
inclusion–exclusion coefficients count covering families.

## Attribution and status

This is an audited research draft, not a peer-reviewed publication or a
proof-assistant formalization. The finite checks supplement the proofs.

Fang and Chang already proved fixed-draw stochastic optimality for `n=2s`
in [their 2014 paper](https://doi.org/10.1016/j.ipl.2014.06.012).
Berend and Sher already proved Schilling's expectation conjecture in
[their paper](https://doi.org/10.1080/15326349.2026.2671340)
([arXiv version](https://arxiv.org/abs/2509.17201)). The present manuscript
gives an alternative proof and a strict stochastic strengthening.

**The general fixed-draw partition comparison remains open.** The
Poissonized comparison does not imply a coefficientwise discrete comparison.
The manuscript proves the earliest feasible draw count for every integral
`q` and recovers the known `q=2` theorem. Separate later working notes on
additional fixed-draw cases are outside this audited release.

Adaptive sampling, heterogeneous quotas, and a full classification of
nondivisible optimizers are also outside the proved scope.

## Reproduce

The verification scripts require Python 3.10 or newer and only its standard
library. PDF compilation requires `latexmk`, pdfLaTeX, and standard TeX Live
packages (`texlive-latex-extra` and `lmodern` on Debian/Ubuntu).

```sh
make verify
make pdf
```

The three scripts independently check the perturbation identities,
preserved marginals, geometric tails, Hessian actions, Poissonized bounds,
and repeated-collection expectations. Recorded outputs are in `audit/`.
GitHub Actions repeats verification and PDF compilation on pushes and pull
requests. The compiled PDF is committed for direct download.

## Citation and license

Use [CITATION.cff](CITATION.cff) for the manuscript's citation metadata.
The manuscript and verification code are released under the [MIT license](LICENSE).
