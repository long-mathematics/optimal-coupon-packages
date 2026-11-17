#!/usr/bin/env python3
"""Exact checks for the revised coupon-package manuscript (Python 3.10+).

Uses only the standard library. All comparisons are exact rational/integer
comparisons; these finite tests supplement, rather than establish, the proofs.
Run from this directory with: python3 verify_strengthenings.py
"""
from __future__ import annotations

from fractions import Fraction as Q
from itertools import combinations
from math import factorial
from random import Random
from time import perf_counter

from verify_results import choose, coefficient, packages, union_dp


def weighted_cdfs(n: int, law: list[tuple[int, int]], max_m: int) -> list[Q]:
    """Union-state DP for arbitrary nonnegative integer package weights."""
    if not law or any(w < 0 for _, w in law):
        raise ValueError("A nonempty, nonnegative weighted law is required.")
    total = sum(w for _, w in law)
    if total <= 0:
        raise ValueError("Total weight must be positive.")
    full = (1 << n) - 1
    state = {0: 1}
    result = [Q(0)]
    for m in range(1, max_m + 1):
        nxt: dict[int, int] = {}
        for seen, mass in state.items():
            for package, weight in law:
                if weight:
                    new = seen | package
                    nxt[new] = nxt.get(new, 0) + mass * weight
        state = nxt
        result.append(Q(state.get(full, 0), total ** m))
    return result


def exact_mean(n: int, law: list[tuple[int, int]]) -> Q:
    """Expectation from an independent finite inclusion--exclusion formula."""
    total = sum(w for _, w in law)
    answer = Q(0)
    for missing in range(1, 1 << n):
        hit = sum(w for A, w in law if A & missing)
        if hit == 0:
            raise ValueError("Infinite completion mean.")
        answer += (-1) ** (missing.bit_count() + 1) * Q(total, hit)
    return answer


def mean_eigenvalue(n: int, s: int, d: int) -> Q:
    N, M = choose(n, s), n - 2*d
    L = choose(M, s-d)
    return (-1)**(d+1) * Q(2*N**3, L) * sum(
        ((-1)**j * choose(M, j)
         * Q(choose(M-j, s-d)**2, (N-choose(M-j+d, s))**3)
         for j in range(M+1)), Q(0))


def hessian_kernels(n: int, s: int, m: int | None = None) -> list[Q]:
    """H[A,B] depends on k = n - |A union B|.

    Derived directly by differentiating the inclusion--exclusion formula,
    independently of the paired-perturbation coefficient formula.
    """
    N = choose(n, s)
    out = []
    for k in range(n-s+1):
        if m is None:
            value = sum((2 * (-1)**(j+1) * choose(k, j)
                         * Q(N**3, (N-choose(n-j, s))**3)
                         for j in range(1, k+1)), Q(0))
        else:
            value = m*(m-1) * sum(
                ((-1)**j * choose(k, j) * Q(choose(n-j, s), N)**(m-2)
                 for j in range(k+1)), Q(0))
        out.append(value)
    return out


def check_eigenvector(n: int, packs: list[tuple[int, int]],
                      kernel: list[Q], eigenvalue: Q) -> int:
    count = 0
    for A, chi_A in packs:
        actual = sum((kernel[n-(A | B).bit_count()] * chi_B
                      for B, chi_B in packs), Q(0))
        assert actual == eigenvalue * chi_A, (n, A, actual, eigenvalue*chi_A)
        count += 1
    return count


def main() -> None:
    started = perf_counter()
    rng = Random(20261112)
    ray_checks = tail_checks = first_checks = complement_checks = 0
    hessian_coordinates = cdf_hessian_coordinates = trace_checks = 0
    hessian_blocks = 0
    max_n, max_m = 8, 8

    # Monotonicity in absolute amplitude. CDFs are computed by union-state DP.
    for n in range(2, max_n+1):
        for s in range(1, n):
            N = choose(n, s)
            for d in range(1, min(s, n-s)+1):
                packs = packages(n, s, d)
                nums = [union_dp(n, packs, a, 3, max_m) for a in range(4)]
                for m in range(max_m+1):
                    for a in range(3):
                        diff = (-1)**d * (nums[a+1][m] - nums[a][m])
                        assert diff >= 0
                        assert (diff > 0) == (m*s >= n)
                        ray_checks += 1

            # Random but reproducible laws on all sizes <=s, including empty.
            allowed = [A for A in range(1 << n) if A.bit_count() <= s]
            for _ in range(4):
                law = [(A, rng.randrange(5)) for A in allowed]
                if not any(w for _, w in law):
                    law[0] = (law[0][0], 1)
                values = weighted_cdfs(n, law, max(max_m, n))
                for m in range(1, max_m+1):
                    assert 1-values[m] >= Q(n-s, n)**(m-1), (n,s,m)
                    tail_checks += 1
                if n % s == 0:
                    q = n//s
                    assert values[q] <= Q(factorial(q), q**q), (n,s,q)
                    first_checks += 1

            # Exact complementary-partition attainment, including endpoints.
            h = n-s
            if n % h == 0:
                k = n//h
                full = (1 << n)-1
                blocks = [sum(1 << i for i in range(j*h,(j+1)*h))
                          for j in range(k)]
                law = [(full ^ B, 1) for B in blocks]
                values = weighted_cdfs(n, law, max_m)
                for m in range(1, max_m+1):
                    assert 1-values[m] == Q(1, k)**(m-1)
                    complement_checks += 1
                assert exact_mean(n, law) == 1+Q(n,s)
                complement_checks += 1

            # Partition attainment at the earliest feasible draw count.
            if n % s == 0:
                q = n//s
                law = [(sum(1 << i for i in range(j*s,(j+1)*s)),1)
                       for j in range(q)]
                values = weighted_cdfs(n,law,q)
                assert values[q] == Q(factorial(q),q**q)
                first_checks += 1

            # Mean Hessian and each feasible fixed-draw Hessian through m=8.
            packs0 = packages(n,s,1)
            mean_kernel = hessian_kernels(n,s)
            eigenvalues = []
            for d in range(1,min(s,n-s)+1):
                packs = packages(n,s,d)
                eigenvalue = mean_eigenvalue(n,s,d)
                assert (-1)**(d+1)*eigenvalue > 0
                assert sum(c*c for _,c in packs) == 2**d*choose(n-2*d,s-d)
                hessian_coordinates += check_eigenvector(n,packs,mean_kernel,eigenvalue)
                eigenvalues.append((d,eigenvalue))
                hessian_blocks += 1
                for m in range(max(2,(n+s-1)//s),max_m+1):
                    gamma = (-1)**d * Q(m*(m-1)*coefficient(n,s,d,m,1),
                                        N**(m-2)*choose(n-2*d,s-d))
                    assert (-1)**d*gamma > 0
                    cdf_hessian_coordinates += check_eigenvector(
                        n,packs,hessian_kernels(n,s,m),gamma)
            # Trace minus the constant-space eigenvalue equals all tangent blocks.
            A = packs0[0][0]
            constant_eigenvalue = sum((mean_kernel[n-(A|B).bit_count()]
                                      for B,_ in packs0),Q(0))
            tangent_trace = N*mean_kernel[n-s]-constant_eigenvalue
            predicted = sum(((choose(n,d)-choose(n,d-1))*e
                             for d,e in eigenvalues),Q(0))
            assert tangent_trace == predicted
            trace_checks += 1

    print('PASS: all strengthening checks use exact rational or integer arithmetic')
    print(f'  n=2..{max_n}, all 1<=s<n and admissible d; draw counts through {max_m}')
    print(f'  paired-ray adjacent-amplitude CDF comparisons: {ray_checks}')
    print('  ray amplitudes: 0, 1/3, 2/3, 1')
    print(f'  geometric-tail checks for seeded at-most-size laws: {tail_checks}')
    print('  4 seeded laws per (n,s); empty and smaller packages permitted')
    print(f'  complementary-partition tail/mean equalities: {complement_checks}')
    print(f'  first-feasible-draw comparisons/partition equalities: {first_checks}')
    print(f'  mean-Hessian blocks: {hessian_blocks}; exact coordinate checks: {hessian_coordinates}')
    print(f'  CDF-Hessian exact coordinate checks: {cdf_hessian_coordinates}')
    print(f'  exact tangent-trace/multiplicity identities: {trace_checks}')
    print('  (8,4) mean-Hessian blocks (d, eigenvalue, multiplicity):')
    for d in range(1,5):
        e=mean_eigenvalue(8,4,d)
        print(f'    {d}: {float(e):.9f}; multiplicity {choose(8,d)-choose(8,d-1)}')
    print(f'Elapsed: {perf_counter()-started:.3f} seconds')


if __name__ == '__main__':
    main()
