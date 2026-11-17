#!/usr/bin/env python3
"""Independent exact checks of the Poisson bound and common-quota means.

Python 3.10+, standard library only. Finite checks do not prove the theorems.
At time t=n*W*log(b), all absence probabilities for integer package weights
with total W are rational. Raising the bound to power s eliminates n/s.
The quota expectation is calculated independently by a capped-count Markov
recursion, eliminating self-loops, and compared with the partition process.
"""
from fractions import Fraction as Q
from functools import lru_cache
from math import prod
from random import Random


def check_poisson(n, s, law, b):
    W = sum(w for _, w in law)
    D = b ** (n*W)
    numerator = sum((-1)**R.bit_count() * b**(n*sum(w for A,w in law if not A&R))
                    for R in range(1 << n))
    marginals = [sum(w for A,w in law if A & (1 << i)) for i in range(n)]
    assert 0 <= numerator <= D
    # Finner: F(t)^s <= product_i (1-exp(-t*p_i)).
    assert numerator**s * D**(n-s) <= prod(D-b**(n*(W-v)) for v in marginals)
    # Jensen: F(t)^s <= (1-exp(-t*s/n))^n.
    assert numerator**s <= (b**(s*W)-1)**n
    return numerator**s == (b**(s*W)-1)**n


def quota_mean(n, law, r):
    W = sum(w for _, w in law)
    assert all(any(w and A & (1 << i) for A,w in law) for i in range(n))
    increments = [(tuple((A >> i) & 1 for i in range(n)),w) for A,w in law if w]

    @lru_cache(None)
    def remaining(state):
        if all(c == r for c in state):
            return Q(0)
        self_weight = 0
        future = Q(W)
        for inc,w in increments:
            nxt = tuple(min(r,c+a) for c,a in zip(state,inc))
            if nxt == state:
                self_weight += w
            else:
                future += w*remaining(nxt)
        return future / (W-self_weight)

    return remaining((0,)*n)


def partition(n,s):
    return [(sum(1 << i for i in range(j,j+s)),1) for j in range(0,n,s)]


def main():
    rng = Random(20261117)
    poisson_checks = attainment_checks = quota_checks = 0
    for n in range(1,8):
        for s in range(1,n+1):
            allowed = [A for A in range(1 << n) if A.bit_count() <= s]
            # Include impossible completion and empty packages.
            laws = [[(0,1)], [(A,rng.randrange(4)) for A in allowed],
                    [(A,1) for A in allowed if A.bit_count() == s]]
            for law in laws:
                if not sum(w for _,w in law):
                    law = [(0,1)]
                for b in (2,3):
                    check_poisson(n,s,law,b)
                    poisson_checks += 1
            if n % s == 0:
                for b in (2,3):
                    assert check_poisson(n,s,partition(n,s),b)
                    attainment_checks += 1
    for n,s in [(1,1),(2,1),(2,2),(3,1),(3,3),(4,1),(4,2),(4,4),(6,2),(6,3),(6,6)]:
        q = n//s
        allowed = [A for A in range(1 << n) if A.bit_count() <= s]
        law = [(A,rng.randrange(1,4)) for A in allowed]
        for r in (1,2,3):
            optimum = quota_mean(q,[(1 << i,1) for i in range(q)],r)
            assert quota_mean(n,partition(n,s),r) == optimum
            assert quota_mean(n,law,r) > optimum
            if r == 1:
                assert optimum == q*sum((Q(1,j) for j in range(1,q+1)),Q(0))
            quota_checks += 1
    print('PASS: independent exact Poisson and common-quota checks')
    print(f'  Poisson CDF/Finner/Jensen chains: {poisson_checks}; n=1..7, all 1<=s<=n')
    print(f'  partition CDF attainment identities: {attainment_checks}')
    print(f'  capped-count quota comparisons and partition identities: {quota_checks}; r=1,2,3')


if __name__ == '__main__':
    main()
