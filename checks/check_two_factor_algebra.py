#!/usr/bin/env python3
"""Check block signs and perturbation coefficients; not a geometric proof.

Uses only the Python standard library. Matrix identities are checked exactly
over the rationals. The final check is a finite-dimensional operator model,
not a numerical discretization of a bundle on a Kahler manifold.
"""

from decimal import Decimal, localcontext
from fractions import Fraction as F


def zero(rows, cols):
    return [[F(0) for _ in range(cols)] for _ in range(rows)]


def transpose(a):
    return [list(row) for row in zip(*a)]


def mul(a, b):
    return [[sum(x * y for x, y in zip(row, col))
             for col in zip(*b)] for row in a]


def add(a, b, sign=1):
    return [[x + sign * y for x, y in zip(ar, br)]
            for ar, br in zip(a, b)]


def comm(a, b):
    return add(mul(a, b), mul(b, a), -1)


def inner(a, b):
    return sum(x * y for ar, br in zip(a, b) for x, y in zip(ar, br))


def check_blocks(a, b, n):
    r = a + b
    q = zero(r, r)
    for i in range(r):
        q[i][i] = F(b, r) if i < a else -F(a, r)
    assert inner(q, q) == F(a * b, r)
    curvature = zero(r, r)
    a2q = zero(r, r)
    energy = F(0)
    for nu in range(n):
        beta = [[F((i + 1) * (j + 2) + 2 * nu - 3, nu + i + j + 1)
                 for j in range(b)] for i in range(a)]
        upper = zero(r, r)
        for i in range(a):
            for j in range(b):
                upper[i][a + j] = beta[i][j]
        bq = comm(upper, q)
        assert bq == [[-x for x in row] for row in upper]
        a2q = add(a2q, comm(transpose(upper), bq))
        ff, gg = mul(beta, transpose(beta)), mul(transpose(beta), beta)
        for i in range(a):
            for j in range(a):
                curvature[i][j] += ff[i][j]
        for i in range(b):
            for j in range(b):
                curvature[a + i][a + j] -= gg[i][j]
        energy += inner(beta, beta)
        v = [[F(2 * i - j + nu, 1 + i + j) for j in range(r)]
             for i in range(r)]
        for i in range(r):
            for j in range(r):
                elementary = zero(r, r)
                elementary[i][j] = F(1)
                assert inner(comm(upper, elementary), v) == inner(
                    elementary, comm(transpose(upper), v))
    assert a2q == curvature
    assert sum(curvature[i][i] for i in range(r)) == 0
    assert inner(q, curvature) == energy
    c = energy / inner(q, q)
    z = add(curvature, [[c * x for x in row] for row in q], -1)
    assert inner(q, z) == 0


def check_operator_model():
    # D0 = diag(0,rho), B = [[alpha,beta],[0,0]], q = (1,0).
    # A1 q = 0 and the spectrum is even, as in the geometric proof.
    alpha, beta, rho = F(2), F(3), F(5)
    c, d = alpha**2, alpha**2 * beta**2 / rho**2
    # Exact t^4 coefficient in det(Delta_t - (c*t^2-d*t^4)I).
    assert c**2 + rho**2 * d - c * (alpha**2 + beta**2) == 0
    # Here the t^6 eigenvalue coefficient is also known explicitly.
    expected_sixth = d * (beta**2 - alpha**2) / rho**2
    with localcontext() as ctx:
        ctx.prec = 90
        dc = Decimal(c.numerator) / Decimal(c.denominator)
        dd = Decimal(d.numerator) / Decimal(d.denominator)
        expected = Decimal(expected_sixth.numerator) / Decimal(expected_sixth.denominator)
        errors = []
        for exponent in (2, 3, 4):
            t = Decimal(10) ** -exponent
            trace = Decimal(25) + Decimal(13) * t**2
            determinant = Decimal(100) * t**2
            large = (trace + (trace**2 - 4 * determinant).sqrt()) / 2
            small = determinant / large
            assert 0 < small <= dc * t**2
            residual = (small - dc * t**2 + dd * t**4) / t**6
            errors.append(abs(residual - expected))
        assert errors[2] < errors[1] < errors[0]
        assert errors[-1] < Decimal('1e-8')
        print(f'Operator model: c={c}, d={d}, t^6 coefficient={expected_sixth}')


def main():
    count = 0
    for a in range(1, 5):
        for b in range(1, 5):
            for n in (1, 2, 3):
                check_blocks(a, b, n)
                count += 1
    print(f'PASS: {count} rational block cases, ranks 2 through 8, 1–3 form components')
    check_operator_model()
    print('PASS: fourth-order sign and sixth-order remainder in the operator model')
    print('Scope: algebra checks only; elliptic analysis and geometric sharpness use the written proof.')


if __name__ == '__main__':
    main()
