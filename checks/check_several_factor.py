#!/usr/bin/env python3
"""Finite-dimensional checks for the several-factor spectral comparison.

Standard library only. The exact checks use rational quadratic forms; the
spectral checks use 110-digit Decimal arithmetic. These are finite-dimensional
models, not geometric discretizations or a proof of the elliptic statements.
"""

from decimal import Decimal, localcontext
from fractions import Fraction as F

from check_two_factor_algebra import add, mul, transpose, zero


def eye(n):
    a = zero(n, n)
    for i in range(n):
        a[i][i] = F(1)
    return a


def scaled(a, s):
    return [[s * x for x in row] for row in a]


def inverse(a):
    n = len(a)
    rows = [row[:] + ident for row, ident in zip(a, eye(n))]
    for i in range(n):
        pivot = next(j for j in range(i, n) if rows[j][i])
        rows[i], rows[pivot] = rows[pivot], rows[i]
        d = rows[i][i]
        rows[i] = [x / d for x in rows[i]]
        for j in range(n):
            if j != i:
                d = rows[j][i]
                rows[j] = [x - d * y for x, y in zip(rows[j], rows[i])]
    return [row[n:] for row in rows]


def is_psd(a):
    """Exact symmetric elimination, including a possible common kernel."""
    a = [row[:] for row in a]
    assert a == transpose(a)
    for i in range(len(a)):
        if a[i][i] < 0:
            return False
        if a[i][i] == 0:
            if any(a[j][i] != 0 for j in range(i + 1, len(a))):
                return False
            continue
        for j in range(i + 1, len(a)):
            for k in range(j, len(a)):
                value = a[j][k] - a[j][i] * a[i][k] / a[i][i]
                a[j][k] = a[k][j] = value
    return True


def dec(x):
    if isinstance(x, F):
        return Decimal(x.numerator) / Decimal(x.denominator)
    return Decimal(x)


def eigenvalues(a):
    """Jacobi diagonalization for the small symmetric matrices used here."""
    a = [[dec(x) for x in row] for row in a]
    n = len(a)
    for _ in range(20000):
        p, q = max(((i, j) for i in range(n) for j in range(i + 1, n)),
                   key=lambda ij: abs(a[ij[0]][ij[1]]))
        off = a[p][q]
        if abs(off) < Decimal('1e-90'):
            return sorted(a[i][i] for i in range(n))
        theta = (a[q][q] - a[p][p]) / (2 * off)
        sign = Decimal(1) if theta >= 0 else Decimal(-1)
        tangent = sign / (abs(theta) + (1 + theta * theta).sqrt())
        cosine = 1 / (1 + tangent * tangent).sqrt()
        sine = tangent * cosine
        a[p][p] -= tangent * off
        a[q][q] += tangent * off
        a[p][q] = a[q][p] = Decimal(0)
        for k in range(n):
            if k not in (p, q):
                kp, kq = a[k][p], a[k][q]
                a[k][p] = a[p][k] = cosine * kp - sine * kq
                a[k][q] = a[q][k] = sine * kp + cosine * kq
    raise AssertionError('Jacobi iteration did not converge')


def graph_matrix(mass_roots, weights, edges, t):
    """Mass-normalized incidence matrix W and Laplacian W*W.

    edges are (i,j,amplitude), giving conductance amplitude^2*t^(2p).
    masses are squares of mass_roots, so the matrix is rational.
    """
    w = zero(len(edges), len(weights))
    for row, (i, j, amplitude) in enumerate(edges):
        coefficient = F(amplitude) * t ** (weights[i] - weights[j])
        w[row][i] = coefficient / mass_roots[i]
        w[row][j] = -coefficient / mass_roots[j]
    return w, mul(transpose(w), w)


def coupled_model(mass_roots, weights, edges, t):
    """D0 vanishes on the vertex space; D0=2I on its complement.

    In target harmonic/range coordinates,
      B_t = [[W_t, Z_t], [0, Y_t]],  D0 = [[0,0],[0,2I]].
    Thus D0* B_t kills the limiting kernel exactly; Z_t gives nonzero coupling.
    """
    w, graph = graph_matrix(mass_roots, weights, edges, t)
    m, h, q = len(weights), len(edges), 2
    pstar = min(weights[i] - weights[j] for i, j, _ in edges)
    b = zero(h + q, m + q)
    d0 = zero(h + q, m + q)
    for i in range(h):
        for j in range(m):
            b[i][j] = w[i][j]
        for j in range(q):
            b[i][m + j] = t**pstar * F((i + 1) * (j + 2) - 1, 10 * (m + q))
    for i in range(q):
        d0[h + i][m + i] = F(2)
        for j in range(q):
            b[h + i][m + j] = t**pstar * F(i - 2 * j - 1, 20 * (m + q))
    d = add(d0, b)
    delta = mul(transpose(d), d)
    eta_bound = sum(abs(x) for row in b for x in row)
    assert eta_bound <= F(1, 2)
    g0 = F(9, 4)  # (sqrt(gamma)-eta)^2 >= (2-1/2)^2, gamma=4.
    epsilon = 2 * eta_bound**2 / g0
    comparison = zero(m + q, m + q)
    for i in range(m):
        for j in range(m):
            comparison[i][j] = (1 - epsilon) * graph[i][j]
    for i in range(q):
        comparison[m + i][m + i] = g0 / 2
    assert is_psd(add(delta, comparison, -1))

    # Exact Schur correction at lambda=0, with a nonzero cross block.
    cross = [row[:m] for row in delta[m:]]
    complement = [row[m:] for row in delta[m:]]
    correction = mul(mul(transpose(cross), inverse(complement)), cross)
    assert is_psd(correction)
    assert is_psd(add(scaled(graph, eta_bound**2 / g0), correction, -1))
    assert any(x for row in cross for x in row)

    graph_spectrum = eigenvalues(graph)
    full_spectrum = eigenvalues(delta)
    tolerance = Decimal('1e-75')
    for mu, lam in zip(graph_spectrum, full_spectrum):
        if abs(mu) < tolerance:
            assert abs(lam) < tolerance
        else:
            assert dec(1 - epsilon) * mu - tolerance <= lam <= mu + tolerance
    assert full_spectrum[m] >= dec(g0 / 2) - tolerance
    return full_spectrum


def check_three_factor():
    """Compare two independent formulas: graph diagonalization and minors."""
    count = 0
    for roots in ([1, 1, 1], [1, 2, 3], [2, 3, 5]):
        masses = [F(x * x) for x in roots]
        m1, m2, m3 = masses
        for p, q in ((1, 2), (1, 5), (2, 5)):
            for chord in (0, 3):
                edges = [(0, 1, 1), (1, 2, 2)]
                if chord:
                    edges.append((0, 2, chord))
                for denominator in (100, 1000):
                    t = F(1, denominator)
                    _, graph = graph_matrix(roots, [p + q, q, 0], edges, t)
                    spectrum = eigenvalues(graph)
                    k12, k23 = t**(2*p), 4*t**(2*q)
                    k13 = chord**2*t**(2*(p+q))
                    tr = k12*(1/m1+1/m2)+k23*(1/m2+1/m3)+k13*(1/m1+1/m3)
                    prod = (m1+m2+m3)/(m1*m2*m3)*(k12*k23+k12*k13+k23*k13)
                    larger = (dec(tr)+(dec(tr)**2-4*dec(prod)).sqrt())/2
                    smaller = dec(prod)/larger
                    assert abs(spectrum[1]/smaller-1) < Decimal('1e-50')
                    assert abs(spectrum[2]/larger-1) < Decimal('1e-50')
                    fast = (1/m1+1/m2)*t**(2*p)
                    slow = 4*(1/(m1+m2)+1/m3)*t**(2*q)
                    assert abs(larger/dec(fast)-1) < Decimal('0.02')
                    assert abs(smaller/dec(slow)-1) < Decimal('0.02')
                    count += 1
    print(f'PASS: {count} three-factor spectra, unequal masses and exponents, optional third edge')


def check_contractions():
    # Equal-mass path with exponents 1,3,2,1. A p=7 chord changes no leading term.
    weights = [7, 6, 3, 1, 0]
    edges = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (0, 4, 2)]
    predictions = [(3, F(5, 6)), (2, F(3, 2)), (1, F(2)), (1, F(2))]
    errors = []
    for denominator in (100, 1000):
        t = F(1, denominator)
        _, graph = graph_matrix([1]*5, weights, edges, t)
        spectrum = eigenvalues(graph)[1:]
        error = max(abs(value/dec(coefficient*t**(2*power))-1)
                    for value, (power, coefficient) in zip(spectrum, predictions))
        errors.append(error)
    assert errors[1] < errors[0] < Decimal('0.001')
    print('PASS: five-factor graph coefficients 5/6, 3/2, 2, 2 at exponents 6, 4, 2, 2')


def check_missing_harmonicity():
    # D0=[0,1], B=[t,0]. Restricted energy on ker(D0) is t^2, but
    # D_t=[t,1] has an exact zero mode (-1,t). D0*B ker(D0) != 0.
    t = F(1, 10)
    delta = [[t*t, t], [t, F(1)]]
    assert mul(delta, [[F(-1)], [t]]) == [[0], [0]]
    assert delta[0][0] > 0
    print('PASS: removing the harmonic cancellation can invalidate relative comparison')


def main():
    with localcontext() as ctx:
        ctx.prec = 110
        cases = [
            ([3, 2, 0], [(0, 1, 1), (1, 2, 1), (0, 2, 2)]),
            ([6, 5, 0], [(0, 1, 1), (1, 2, 2)]),
            ([3, 2, 1, 0], [(0, 1, 1), (2, 3, 1)]),
            ([7, 6, 3, 1, 0], [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1)]),
        ]
        count = 0
        for weights, edges in cases:
            for roots in ([1]*len(weights), list(range(1, len(weights)+1))):
                for denominator in (50, 200, 1000):
                    coupled_model(roots, weights, edges, F(1, denominator))
                    count += 1
        print(f'PASS: {count} exact rational form and Schur bounds, plus coupled spectral comparisons')
        check_three_factor()
        check_contractions()
        check_missing_harmonicity()
        print('Illustration: coupled model with unit masses, p=1, q=2, edge norms=1')
        for denominator in (50, 200, 1000):
            t = F(1, denominator)
            spectrum = coupled_model([1]*3, [3, 2, 0], [(0, 1, 1), (1, 2, 1)], t)
            print(f'  t=1/{denominator}: lambda_slow/t^4={spectrum[1]/dec(t**4):.12f}, '
                  f'lambda_fast/t^2={spectrum[2]/dec(t**2):.12f}')
    print('Scope: finite-dimensional checks only; the analytic arguments are proved in the paper.')


if __name__ == '__main__':
    main()
