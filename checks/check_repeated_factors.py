#!/usr/bin/env python3
"""Exact multiplicity-matrix checks and finite spectral models.

Formal orthonormal harmonic channels retain contributions landing in the
same matrix block. All form and rank checks use rational arithmetic;
eigenvalues use 110-digit Decimal arithmetic. No elliptic or geometric
discretization is involved.
"""

from decimal import Decimal, localcontext
from fractions import Fraction as F

from check_two_factor_algebra import add, mul, transpose, zero
from check_several_factor import dec, eigenvalues, inverse, is_psd, scaled


def unit(n, i, j):
    a = zero(n, n)
    a[i][j] = F(1)
    return a


def inner(a, b):
    return sum(x * y for ar, br in zip(a, b) for x, y in zip(ar, br))


def basis(types):
    """Orthogonal rational basis for trace-free isotypical matrices.

    All representative ranks and the volume are one in these models.
    The diagonal Helmert vectors have squared lengths k*(k+1).
    """
    n = len(types)
    answer = []
    for k in range(1, n):
        a = zero(n, n)
        for i in range(k):
            a[i][i] = F(1)
        a[k][k] = F(-k)
        answer.append(a)
    answer.extend(unit(n, i, j) for i in range(n) for j in range(n)
                  if i != j and types[i] == types[j])
    norms = [inner(x, x) for x in answer]
    assert all(inner(x, y) == (norms[i] if i == j else 0)
               for i, x in enumerate(answer) for j, y in enumerate(answer))
    return answer, norms


def commutator_map(types, weights, channels, t):
    """A channel lists (row,column,amplitude) sharing one harmonic form.

    Channels are orthonormal. Entries within one channel must have the
    same source and target stable types. They are summed before computing
    the squared norm, so matrix interference is preserved.
    """
    n = len(types)
    base, norms = basis(types)
    matrices = []
    for channel in channels:
        assert len({(types[i], types[j]) for i, j, _ in channel}) == 1
        a = zero(n, n)
        for i, j, coefficient in channel:
            assert i < j and weights[i] > weights[j]
            a[i][j] += F(coefficient) * t ** (weights[i] - weights[j])
        matrices.append(a)
    columns = []
    for x in base:
        column = []
        for a in matrices:
            comm = add(mul(a, x), mul(x, a), -1)
            column.extend(v for row in comm for v in row)
        columns.append(column)
    w = transpose(columns)
    return w, mul(transpose(w), w), norms


def rank(a):
    a = [row[:] for row in a]
    r = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        value = a[r][c]
        a[r] = [x / value for x in a[r]]
        for i in range(r + 1, len(a)):
            value = a[i][c]
            a[i] = [x - value*y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def spectrum(form, norms):
    roots = [dec(x).sqrt() for x in norms]
    return eigenvalues([[dec(x) / (roots[i]*roots[j])
                         for j, x in enumerate(row)]
                        for i, row in enumerate(form)])


def check_four_factor():
    types, weights = [0, 1, 0, 2], [4, 3, 1, 0]
    channels = [[(0, 1, 1)], [(1, 2, 1)], [(2, 3, 1)]]
    base, _ = basis(types)
    assert len(base) == 5  # 2^2+1^2+1^2-1, not four minus one.
    for denominator in (10, 50, 200):
        t = F(1, denominator)
        _, form, norms = commutator_map(types, weights, channels, t)

        def energy(x):
            diagonal = [x[i][i] for i in range(4)]
            x1, x2, x3, x4 = diagonal
            u, v = x[0][2], x[2][0]
            return (t*t*((x2-x1)**2+(x4-x3)**2+u*u+v*v)
                    + t**4*((x3-x2)**2+v*v))

        predicted = [[(energy(add(x, y))-energy(x)-energy(y))/2
                      for y in base] for x in base]
        assert form == predicted
        assert rank(form) == 5
        values = spectrum(form, norms)
        a, b = dec(t*t), dec(t**4)
        fast = a+b+(a*a+b*b).sqrt()
        slow = 2*a*b/fast  # Avoid cancellation in a+b-sqrt(a^2+b^2).
        expected = sorted([slow, a, a+b, 2*a, fast])
        assert max(abs(x/y-1) for x, y in zip(values, expected)) < Decimal('1e-70')
        leading = [b, a, a, 2*a, 2*a]
        assert max(abs(x/y-1) for x, y in zip(values, leading)) < 2*dec(t*t)
    print('PASS: repeated A,B,A,C example; exact form, five positive modes, coefficients 1,1,1,2,2')


def check_connected_self_extension():
    _, form, norms = commutator_map([0, 0], [1, 0], [[(0, 1, 1)]], F(1, 10))
    # Basis: diag(1,-1), E12, E21. E12 is the extra commuting nilpotent.
    assert len(norms) == 3 and rank(form) == 2
    assert all(row[1] == 0 for row in form)
    values = spectrum(form, norms)
    assert values[0] == 0
    assert all(abs(x-Decimal('0.02')) < Decimal('1e-80') for x in values[1:])
    print('PASS: connected self-extension is nonsimple; nilpotent zero mode retained')


def check_interference():
    types, weights, t = [0, 0, 0], [3, 2, 0], F(1, 7)
    _, combined, _ = commutator_map(types, weights, [[(0, 1, 1), (1, 2, 1)]], t)
    _, separated, _ = commutator_map(types, weights, [[(0, 1, 1)], [(1, 2, 1)]], t)
    assert combined != separated
    # The same geometric form cannot be replaced by orthogonal edge labels.
    n = add(scaled(unit(3, 0, 1), t), scaled(unit(3, 1, 2), t*t))
    base, norms = basis(types)
    vector = [[inner(n, x)/norm] for x, norm in zip(base, norms)]
    combined_energy = mul(mul(transpose(vector), combined), vector)[0][0]
    separated_energy = mul(mul(transpose(vector), separated), vector)[0][0]
    assert combined_energy == 0 and separated_energy == 2*t**6
    print('PASS: shared harmonic entries have nonzero cross terms; separate-edge energy fails')


def coupled_model(types, weights, channels, t):
    w, form, norms = commutator_map(types, weights, channels, t)
    h, d, q = len(w), len(norms), 2
    pstar = min(weights[i]-weights[j] for channel in channels for i, j, _ in channel)
    b, d0 = zero(h+q, d+q), zero(h+q, d+q)
    for i in range(h):
        b[i][:d] = w[i]
        for j in range(q):
            b[i][d+j] = t**pstar * F((i+2*j) % 5-2, 100*(h+d))
    for i in range(q):
        d0[h+i][d+i] = F(2)
        for j in range(q):
            b[h+i][d+j] = t**pstar * F(i-j+1, 100*(h+d))
    delta = mul(transpose(add(d0, b)), add(d0, b))
    # Basis vectors have norms >=1; entrywise l1 norm bounds the operator norm.
    eta_bound = sum(abs(x) for row in b for x in row)
    assert eta_bound <= F(1, 2)
    g0 = F(9, 4)  # gamma=4 and eta<=1/2.
    epsilon = 2*eta_bound**2/g0
    lower = zero(d+q, d+q)
    for i in range(d):
        for j in range(d):
            lower[i][j] = (1-epsilon)*form[i][j]
    for i in range(q):
        lower[d+i][d+i] = g0/2
    assert is_psd(add(delta, lower, -1))
    cross = [row[:d] for row in delta[d:]]
    complement = [row[d:] for row in delta[d:]]
    correction = mul(mul(transpose(cross), inverse(complement)), cross)
    assert is_psd(correction)
    assert is_psd(add(scaled(form, eta_bound**2/g0), correction, -1))
    assert any(x for row in cross for x in row)
    assert d+q-rank(delta) == d-rank(form)
    finite, full = spectrum(form, norms), spectrum(delta, norms+[F(1)]*q)
    tolerance = Decimal('1e-70')
    for nu, lam in zip(finite, full):
        assert dec(1-epsilon)*nu-tolerance <= lam <= nu+tolerance
    assert full[d] >= dec(g0/2)-tolerance


def main():
    with localcontext() as ctx:
        ctx.prec = 110
        check_four_factor()
        check_connected_self_extension()
        check_interference()
        cases = [
            ([0, 1, 0, 2], [4, 3, 1, 0], [[(0, 1, 1)], [(1, 2, 1)], [(2, 3, 1)]]),
            ([0, 0, 0], [3, 2, 0], [[(0, 1, 1), (1, 2, 1)]]),
            ([0, 0, 0], [4, 1, 0], [[(0, 1, 2), (1, 2, 1)], [(0, 2, 1)]]),
            ([0, 0], [2, 0], [[(0, 1, 2)]]),
        ]
        count = 0
        for types, weights, channels in cases:
            for denominator in (100, 300, 1000):
                coupled_model(types, weights, channels, F(1, denominator))
                count += 1
        print(f'PASS: {count} exact rational relative-form, Schur and kernel checks; coupled spectra')
    print('Scope: finite matrix algebra only; the geometric and elliptic arguments are proved in the paper.')


if __name__ == '__main__':
    main()
