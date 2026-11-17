#!/usr/bin/env python3
"""Exact-arithmetic checks for the coupon-package draft (Python 3.10+).

No third-party modules are required. The proofs do not depend on these checks.
The default run compares an independent union-state dynamic program with the
paired-perturbation identity, checks preserved low-order inclusion probabilities,
and tests nonnegative covering coefficients. Run with --max-n to change the
small-state exhaustive range.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from itertools import combinations
from math import comb
from time import perf_counter


def choose(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0


def coefficient(n: int, s: int, d: int, m: int, ell: int) -> int:
    """The nonnegative integer K^{(d)}_{m,ell} in the manuscript."""
    v = n - 2 * d
    return sum(
        (-1) ** r * choose(v, r)
        * choose(v - r, s - d) ** (2 * ell)
        * choose(v - r + d, s) ** (m - 2 * ell)
        for r in range(v + 1)
    )


def packages(n: int, s: int, d: int) -> list[tuple[int, int]]:
    out = []
    for a in combinations(range(n), s):
        mask = sum(1 << i for i in a)
        chi = 1
        for j in range(d):
            chi *= ((mask >> (2*j)) & 1) - ((mask >> (2*j+1)) & 1)
        out.append((mask, chi))
    return out


def union_dp(n: int, packs: list[tuple[int, int]], a: int, b: int,
             max_m: int) -> list[int]:
    """CDF numerators for theta=a/b, denominator (b*binom(n,s))**m."""
    full = (1 << n) - 1
    law = {0: 1}
    values = [0]
    weighted = [(mask, b + a*chi) for mask, chi in packs if b + a*chi]
    for _ in range(max_m):
        nxt: dict[int, int] = {}
        for old, mass in law.items():
            for mask, weight in weighted:
                new = old | mask
                nxt[new] = nxt.get(new, 0) + mass*weight
        law = nxt
        values.append(law.get(full, 0))
    return values


def exact_expectation(n: int, s: int, d: int, theta: Fraction) -> Fraction:
    packs = packages(n, s, d)
    total = Fraction(0)
    N = choose(n, s)
    for missing in range(1, 1 << n):
        avoid = sum(((1 + theta*chi) for mask, chi in packs
                     if not (mask & missing)), Fraction(0)) / N
        if avoid == 1:
            raise ValueError('Completion expectation is infinite.')
        total += (-1) ** (missing.bit_count() + 1) / (1 - avoid)
    return total


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-n', type=int, default=9)
    parser.add_argument('--max-m', type=int, default=8)
    parser.add_argument('--coefficient-max-n', type=int, default=35)
    args = parser.parse_args()
    if args.max_n < 2 or args.max_m < 2 or args.coefficient_max_n < 2:
        parser.error('All range limits must be at least 2.')
    start = perf_counter()
    laws = cdf_checks = marginal_checks = coefficient_checks = 0
    for n in range(2, args.max_n + 1):
        for s in range(1, n):
            for d in range(1, min(s, n-s) + 1):
                packs = packages(n, s, d)
                assert sum(chi for _, chi in packs) == 0
                for j in range(d):
                    for subset in combinations(range(n), j):
                        mask = sum(1 << i for i in subset)
                        assert sum(chi for p, chi in packs if p & mask == mask) == 0
                        marginal_checks += 1
                base = union_dp(n, packs, 0, 1, args.max_m)
                for a, b in [(1, 2), (-1, 2), (1, 1)]:
                    perturbed = union_dp(n, packs, a, b, args.max_m)
                    laws += 1
                    for m in range(args.max_m + 1):
                        diff = perturbed[m] - b**m * base[m]
                        formula = (-1)**d * 2**d * sum(
                            choose(m, 2*ell) * a**(2*ell) * b**(m-2*ell)
                            * coefficient(n, s, d, m, ell)
                            for ell in range(1, m//2 + 1)
                        )
                        assert diff == formula, (n, s, d, m, a, b)
                        assert (-1)**d * diff >= 0
                        assert (diff != 0) == (m*s >= n), (n,s,d,m)
                        cdf_checks += 1
    for n in range(2, args.coefficient_max_n + 1):
        for s in range(1, n):
            for d in range(1, min(s, n-s) + 1):
                for m in range(2, max(7, (n+s-1)//s + 2)):
                    for ell in range(1, m//2 + 1):
                        k = coefficient(n,s,d,m,ell)
                        assert k >= 0, (n,s,d,m,ell,k)
                        if ell == 1:
                            assert (k > 0) == (m*s >= n)
                        coefficient_checks += 1
    assert exact_expectation(4,2,2,Fraction(0)) == Fraction(19,5)
    assert exact_expectation(4,2,2,Fraction(1)) == Fraction(18,5)
    assert exact_expectation(3,2,1,Fraction(0)) == Fraction(5,2)
    assert exact_expectation(4,2,2,Fraction(1,2)) == Fraction(619,165)
    for theta in [Fraction(0), Fraction(1,3), Fraction(1,2), Fraction(1)]:
        assert exact_expectation(4,2,2,theta) == (
            Fraction(19,5) - 24*theta**2 / (5*(25-theta**2)))
    print('PASS: exact paired-switch dynamic-programming identities')
    print(f'  n=2..{args.max_n}; m=0..{args.max_m}; all s and admissible d')
    print(f'  {laws} perturbed laws; {cdf_checks} exact CDF comparisons')
    print(f'  {marginal_checks} low-order inclusion-preservation checks')
    print(f'  {coefficient_checks} coefficient nonnegativity / strictness checks')
    print(f'  coefficient grid: n=2..{args.coefficient_max_n}')
    print('PASS: exact expectations 19/5, 18/5, 619/165, and 5/2; rational example formula')
    print(f'Elapsed: {perf_counter()-start:.3f} seconds')


if __name__ == '__main__':
    main()
