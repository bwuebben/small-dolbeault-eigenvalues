#!/usr/bin/env python3
"""Exact checks of the finite triangular gauge induction.

The coefficient ring is Q[z,z^-1], with derivation d(z^k)=k*z^k.
Constant terms model the harmonic projection; nonconstant terms have a
unique primitive with zero constant term. This checks the nonlinear gauge
algebra, not Dolbeault Hodge theory or the spectral theorem on a curve.
"""

from fractions import Fraction as F
from random import Random


def add(a, b):
    c = dict(a)
    for k, v in b.items():
        c[k] = c.get(k, F(0)) + v
        if not c[k]:
            del c[k]
    return c


def scale(a, q):
    return {k: q * v for k, v in a.items() if q * v}


def mul(a, b):
    c = {}
    for i, x in a.items():
        for j, y in b.items():
            c[i + j] = c.get(i + j, F(0)) + x * y
    return {k: v for k, v in c.items() if v}


def derivative(a):
    return {k: k * v for k, v in a.items() if k}


def primitive(a):
    return {k: v / k for k, v in a.items() if k}


def zero(n):
    return [[{} for _ in range(n)] for _ in range(n)]


def identity(n):
    a = zero(n)
    for i in range(n):
        a[i][i] = {0: F(1)}
    return a


def madd(a, b):
    n = len(a)
    return [[add(a[i][j], b[i][j]) for j in range(n)] for i in range(n)]


def mscale(a, q):
    return [[scale(x, q) for x in row] for row in a]


def mmul(a, b):
    n = len(a)
    c = zero(n)
    for i in range(n):
        for k in range(n):
            if a[i][k]:
                for j in range(n):
                    c[i][j] = add(c[i][j], mul(a[i][k], b[k][j]))
    return c


def md(a):
    return [[derivative(x) for x in row] for row in a]


def unipotent_inverse(u):
    n = len(u)
    inverse = identity(n)
    power = identity(n)
    for k in range(1, n):
        power = mmul(power, u)
        inverse = madd(inverse, mscale(power, F((-1) ** k)))
    assert mmul(power, u) == zero(n)
    return inverse


def normal_form(a):
    n = len(a)
    b, total = a, identity(n)
    for distance in range(1, n):
        correction = zero(n)
        for i in range(n - distance):
            j = i + distance
            correction[i][j] = scale(primitive(b[i][j]), F(-1))
        g = madd(identity(n), correction)
        inverse = unipotent_inverse(correction)
        assert mmul(inverse, g) == identity(n)
        assert mmul(g, inverse) == identity(n)
        after = mmul(inverse, madd(mmul(b, g), md(g)))
        for i in range(n):
            for j in range(n):
                if j - i < distance:
                    assert after[i][j] == b[i][j]
                if j - i == distance:
                    expected = {0: b[i][j][0]} if 0 in b[i][j] else {}
                    assert after[i][j] == expected
        total = mmul(total, g)
        b = after
        # Check the full accumulated gauge relation at every stage.
        assert madd(md(total), mmul(a, total)) == mmul(total, b)
    assert all(not derivative(x) for row in b for x in row)
    assert all(not b[i][j] for i in range(n) for j in range(i + 1))
    return b, total


def conjugate_weights(a, weights, t):
    return [[scale(x, t ** (weights[i] - weights[j]))
             for j, x in enumerate(row)] for i, row in enumerate(a)]


def main():
    # Lower entries can generate a harmonic higher entry even when all
    # original entries have zero harmonic projection.
    a = zero(3)
    a[0][1], a[1][2] = {1: F(1)}, {-1: F(1)}
    b, _ = normal_form(a)
    expected = zero(3)
    expected[0][2] = {0: F(1)}
    assert b == expected
    assert all(0 not in x for row in a for x in row)
    print('PASS: three-block interaction creates beta_13 = 1 from A_12=z, A_23=z^-1.')

    rng = Random(20260920)
    cases = 0
    for n in range(2, 8):
        for _ in range(6):
            a = zero(n)
            for i in range(n):
                for j in range(i + 1, n):
                    a[i][j] = {k: F(rng.randint(-3, 3), rng.randint(1, 4))
                               for k in (-2, -1, 0, 1, 2)}
                    a[i][j] = {k: v for k, v in a[i][j].items() if v}
            b, u = normal_form(a)
            weights = [sum(range(n - i)) for i in range(n)]
            t = F(2, 5)
            at, bt, ut = [conjugate_weights(x, weights, t) for x in (a, b, u)]
            assert madd(md(ut), mmul(at, ut)) == mmul(ut, bt)
            cases += 1
    print(f'PASS: {cases} exact rational Laurent-polynomial cases, 2-7 blocks.')
    print('PASS: each earlier diagonal is preserved; full gauge identities and weighted scaling hold.')
    print('Scope: finite gauge algebra only; Hodge theory and spectral recovery require the written proof.')


if __name__ == '__main__':
    main()
