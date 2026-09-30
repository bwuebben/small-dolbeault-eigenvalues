#!/usr/bin/env python3
"""Finite checks for the three- and four-factor flow calculation.

Exact rational checks cover the moment map, time normalization, energy
dissipation, and changes of variables. High-precision graph spectra and
the closed scalar solution check the separated asymptotic coefficients.
These checks do not verify the geometric PDE, metric matching, or elliptic
and parabolic estimates in the manuscript.
"""

from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from itertools import product
from math import exp, log, sqrt

from check_two_factor_algebra import add, mul, transpose, zero
from check_several_factor import dec, eigenvalues


def edges(m):
    return [(0, 1), (2, 1)] + ([(2, 3)] if m == 4 else [])


def moments(m, weights):
    p = [F(0)] * m
    for (i, j), weight in zip(edges(m), weights):
        p[i] += weight
        p[j] -= weight
    return p


def vector_field(m, weights):
    p = moments(m, weights)
    return [-2 * weight * (p[i] - p[j])
            for (i, j), weight in zip(edges(m), weights)]


def graph(m, weights):
    result = zero(m, m)
    for (i, j), weight in zip(edges(m), weights):
        result[i][i] += weight
        result[j][j] += weight
        result[i][j] -= weight
        result[j][i] -= weight
    return result


def check_moments_and_energy():
    count, noncommuting = 0, 0
    for m in (3, 4):
        for amplitudes in product((F(1, 3), F(1), F(3, 2)), repeat=m-1):
            n = zero(m, m)
            for (i, j), amplitude in zip(edges(m), amplitudes):
                n[i][j] = amplitude
            curvature = add(mul(n, transpose(n)), mul(transpose(n), n), -1)
            weights = [value**2 for value in amplitudes]
            p = moments(m, weights)
            assert [curvature[i][i] for i in range(m)] == p
            assert sum(p) == 0
            assert any(curvature[i][j] != 0 for i in range(m)
                       for j in range(m) if i != j)
            if any((p[i]-p[j])*curvature[i][j] != 0
                   for i in range(m) for j in range(m)):
                noncommuting += 1
            prime = vector_field(m, weights)
            energy = sum(weights)
            assert sum(prime) == -2*sum(value**2 for value in p)
            assert sum(prime) <= -F(8, m)*energy**2
            a, b = weights[:2]
            assert prime[0] == -2*a*(2*a+b)
            if m == 3:
                assert prime[1] == -2*b*(a+2*b)
            else:
                c = weights[2]
                assert prime[1] == -2*b*(a+2*b+c)
                assert prime[2] == -2*c*(2*c+b)
            count += 1
    assert noncommuting > 0
    print(f'PASS: {count} exact curvature, projected-flow, and energy cases; '
          f'{noncommuting} have noncommuting diagonal and off-diagonal curvature')


def check_changes_of_variables():
    count3, count4 = 0, 0
    values = (F(1, 5), F(2, 3), F(1), F(7, 2))
    for a, b in product(values, repeat=2):
        ap, bp = vector_field(3, [a, b])
        e, ep = a+b, ap+bp
        r = (a-b)/e
        rp = (ap-bp)/e-r*ep/e
        assert ep == -(3+r*r)*e*e
        assert rp == -e*r*(1-r*r)
        count3 += 1
    for u, r, z in product(values, repeat=3):
        a, b, c = u*z, u*r, u/z
        ap, bp, cp = vector_field(4, [a, b, c])
        up = u*(ap/a+cp/c)/2
        rp = bp/u-r*up/u
        zp = z*(ap/a-cp/c)/2
        w, wp = (z-1)/(z+1), 2*zp/(z+1)**2
        # Divide by ds/dtau=2u to obtain all three printed identities.
        assert rp/(2*u) == -r*r
        assert wp/(2*u) == -2*w
        assert up/(2*u*u) == -(2+4*w*w/(1-w*w)+r)
        count4 += 1
    print(f'PASS: {count3} three-factor and {count4} four-factor exact changes of variables')


def check_graph_spectra():
    matrix = graph(3, [F(1, 6), F(1, 6)])
    for vector, eigenvalue in (([1, 0, -1], F(1, 6)),
                              ([1, -2, 1], F(1, 2))):
        assert [sum(matrix[i][j]*vector[j] for j in range(3))
                for i in range(3)] == [eigenvalue*x for x in vector]
    cases = 0
    for scale in (F(1), F(1, 1000)):
        for size in (10, 100, 1000, 10000):
            # The strong edges are also unequal; this checks cluster
            # asymptotics without imposing the symmetry a=c.
            for sign in (-1, 0, 1):
                a, b = scale, scale/F(size)
                c = scale*(1+F(sign, size))
                spectrum = eigenvalues(graph(4, [a, b, c]))
                assert abs(spectrum[0]) < D('1e-85')
                ratios = [spectrum[1]/dec(b), spectrum[2]/dec(2*a),
                          spectrum[3]/dec(2*a)]
                assert all(abs(ratio-1) < D(4)/D(size) for ratio in ratios)
                cases += 1
    print(f'PASS: exact three-factor eigenvectors and {cases} separated four-factor graph spectra')


def scalar_solution(a0, b0, c0, s):
    """Exact positive scalar solution, parametrized by s; tau(0)=1.

    Integrating d tau/ds=1/(2u) is elementary. This is not a PDE solve.
    """
    u0, z0 = sqrt(a0*c0), sqrt(a0/c0)
    offset, w0 = u0/b0, (z0-1)/(z0+1)
    constant = u0*offset*(1-w0*w0)
    w = w0*exp(-2*s)
    z = (1+w)/(1-w)
    u = constant*exp(-2*s)/((s+offset)*(1-w*w))
    plus = (exp(2*s)*(s+offset-0.5)-(offset-0.5))/2
    minus = ((offset+0.5)-exp(-2*s)*(s+offset+0.5))/2
    tau = 1+(plus-w0*w0*minus)/(2*constant)
    return tau, u*z, u/(s+offset), u/z


def check_scalar_asymptotics():
    cases = ((1., 1., 1.), (0.1, 0.3, 2.), (2., 0.4, 0.2), (0.3, 2., 0.7))
    worst_error = 0.
    for initial in cases:
        tau, a, b, c = scalar_solution(*initial, 0.)
        assert abs(tau-1) < 1e-12
        assert max(abs(x-y) for x, y in zip((a, b, c), initial)) < 1e-12
        # Finite-difference verification of the closed solution in the
        # original time variable, at a moderate parameter value.
        step, s = 1e-5, 2.
        before = scalar_solution(*initial, s-step)
        after = scalar_solution(*initial, s+step)
        now = scalar_solution(*initial, s)
        numerical = [(after[i]-before[i])/(after[0]-before[0]) for i in (1, 2, 3)]
        expected = vector_field(4, list(now[1:]))
        assert max(abs(x/y-1) for x, y in zip(numerical, expected)) < 1e-7
        tau, a, b, c = scalar_solution(*initial, 300.)
        ratios = (4*tau*a, 2*tau*log(tau)*b, 4*tau*c)
        error = max(abs(ratio-1) for ratio in ratios)
        assert error < 0.015
        worst_error = max(worst_error, error)
    print(f'PASS: {len(cases)} unequal-initial-data closed scalar solutions; '
          f'largest relative coefficient error at s=300 is {worst_error:.6f}')


def main():
    check_moments_and_energy()
    check_changes_of_variables()
    with localcontext() as context:
        context.prec = 110
        check_graph_spectra()
    check_scalar_asymptotics()
    print('Scope: finite algebra and scalar/graph diagnostics only; '
          'the PDE estimates, matching lemma, and spectral-subspace convergence are proved in the paper.')


if __name__ == '__main__':
    main()
