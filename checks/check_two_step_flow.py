#!/usr/bin/env python3
"""Finite checks for two_step_flow.tex; standard library only.

Exact rational identities test ranks, volume, dissipation, and variable
changes. Independent matrix scaling and numerical ODE/graph calculations
test the model coefficients. None is a verification of the bundle PDE,
late-time matching, or uniform elliptic/parabolic comparison.
"""

from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from itertools import product
from math import exp, log, sqrt

from check_two_factor_algebra import add, mul, transpose, zero
from check_several_factor import dec, eigenvalues


def moments(masses, edges, weights):
    result = [F(0)] * len(masses)
    for (i, j), weight in zip(edges, weights):
        result[i] += weight / masses[i]
        result[j] -= weight / masses[j]
    return result


def field(masses, edges, weights):
    p = moments(masses, edges, weights)
    return [-2 * k * (p[i] - p[j]) for (i, j), k in zip(edges, weights)]


def laplacian(n, edges, weights):
    result = zero(n, n)
    for (i, j), k in zip(edges, weights):
        result[i][i] += k
        result[j][j] += k
        result[i][j] -= k
        result[j][i] -= k
    return result


def check_block_masses():
    cases = 0
    for ranks in ((2, 1, 1, 2), (1, 3, 2, 4), (2, 2, 3, 1)):
        offsets = [sum(ranks[:i]) for i in range(5)]
        for edges in ([(0, 1), (2, 1), (2, 3)],
                      [(0, 1), (0, 3), (2, 1), (2, 3)]):
            for volume, weak in product((F(1, 3), F(1), F(5)),
                                        (F(1), F(1, 10**8))):
                n = zero(sum(ranks), sum(ranks))
                for index, (i, j) in enumerate(edges):
                    for a in range(ranks[i]):
                        for b in range(ranks[j]):
                            n[offsets[i]+a][offsets[j]+b] = (
                                F(1+a+2*b, 3+index) * (weak if index == 1 else 1))
                curvature = add(mul(n, transpose(n)), mul(transpose(n), n), -1)
                weights = [volume * sum(n[a][b]**2
                           for a in range(offsets[i], offsets[i+1])
                           for b in range(offsets[j], offsets[j+1]))
                           for i, j in edges]
                masses = [volume*r for r in ranks]
                p = moments(masses, edges, weights)
                assert p == [sum(curvature[a][a]
                             for a in range(offsets[i], offsets[i+1])) / ranks[i]
                             for i in range(4)]
                assert sum(m*x for m, x in zip(masses, p)) == 0
                assert sum(m*abs(x) for m, x in zip(masses, p)) == 2*sum(weights)
                ep = sum(field(masses, edges, weights))
                assert ep == -2*sum(m*x*x for m, x in zip(masses, p))
                assert ep <= -8*sum(weights)**2/sum(masses)
                cases += 1
    print(f'PASS: {cases} rectangular-block rank/volume and energy identities, '
          'including edge-amplitude ratio 1e-8')


def check_unequal_rank_changes():
    edges, masses = [(0, 1), (2, 1), (2, 3)], [2, 1, 1, 2]
    cases = 0
    for u, r, z in product((F(1, 5), F(2, 3), F(1), F(7, 2)), repeat=3):
        a, b, c = u*z, u*r, u/z
        ap, bp, cp = field(masses, edges, [a, b, c])
        assert (ap, bp, cp) == (-3*a*a-2*a*b, -2*(a+c)*b-4*b*b, -3*c*c-2*b*c)
        up = u*(ap/a+cp/c)/2
        rp = bp/u-r*up/u
        zp = z*(ap/a-cp/c)/2
        w, wp = (z-1)/(z+1), 2*zp/(z+1)**2
        cosh = (z+1/z)/2
        assert wp/u == -3*w
        assert rp/u == -r*(cosh+2*r)
        assert up/u**2 == -3*cosh-2*r
        assert cosh == 1+2*w*w/(1-w*w)
        cases += 1
    print(f'PASS: {cases} exact unequal-rank ODE and time-change identities')


def rk4(function, state, step):
    k1 = function(state)
    k2 = function([x+step*v/2 for x, v in zip(state, k1)])
    k3 = function([x+step*v/2 for x, v in zip(state, k2)])
    k4 = function([x+step*v for x, v in zip(state, k3)])
    return [x+step*(a+2*b+2*c+d)/6
            for x, a, b, c, d in zip(state, k1, k2, k3, k4)]


def check_balance():
    cases = [([1, 1, 1, 1], [0, 2], [1, 3], [(0, 1), (0, 3), (2, 1), (2, 3)]),
             ([1, 2, 3, 2, 1], [0, 1], [2, 3, 4],
              [(i, j) for i in (0, 1) for j in (2, 3, 4)]),
             ([1]*6, [0, 1, 2], [3, 4, 5],
              [(0, 3), (0, 4), (1, 4), (1, 5), (2, 5), (2, 3)])]
    for masses, targets, sources, edges in cases:
        total = sum(masses)
        mi, mj = sum(masses[i] for i in targets), sum(masses[j] for j in sources)
        q = [F(mj, 2*total) if i in targets else -F(mi, 2*total)
             for i in range(len(masses))]
        assert sum(m*x for m, x in zip(masses, q)) == 0
        assert all(q[i]-q[j] == F(1, 2) for i, j in edges)
        base = [0.07*(1+index) for index in range(len(edges))]
        # Iterative proportional fitting is separate from the ODE solver.
        balanced = base[:]
        for _ in range(4000):
            for vertices, endpoint in ((targets, 0), (sources, 1)):
                for i in vertices:
                    indices = [k for k, edge in enumerate(edges) if edge[endpoint] == i]
                    ratio = abs(masses[i]*float(q[i]))/sum(balanced[k] for k in indices)
                    for k in indices:
                        balanced[k] *= ratio
        def gradient_flow(y):
            weights = [k*exp(2*(y[i]-y[j])) for (i, j), k in zip(edges, base)]
            p = moments(masses, edges, weights)
            return [float(a)-b for a, b in zip(q, p)]
        for initial in (0., 0.12):
            y = [initial*(i-2) for i in range(len(masses))]
            mean = sum(m*x for m, x in zip(masses, y))/total
            y = [x-mean for x in y]
            for _ in range(12000):
                y = rk4(gradient_flow, y, 0.025)
            actual = [k*exp(2*(y[i]-y[j])) for (i, j), k in zip(edges, base)]
            assert max(abs(a/b-1) for a, b in zip(actual, balanced)) < 1e-9
            assert abs(sum(m*x for m, x in zip(masses, y))) < 1e-10
    # Strict positivity fails at the equal-rank zigzag boundary, and even
    # nonnegative feasibility fails for ranks (2,1,1,2).
    for masses, expected in (([1, 1, 1, 1], F(0)), ([2, 1, 1, 2], -F(1, 4))):
        # At equal side masses the required incident sums are M_i/4.
        forced_a = F(masses[0], 4)
        forced_middle = F(masses[1], 4)-forced_a
        assert forced_middle == expected
    print('PASS: independent matrix scaling versus log-time flow on 3 graphs, '
          '2 initial conditions each; both zigzag balance obstructions')


def check_graph_coefficients():
    edges = [(0, 1), (0, 3), (2, 1), (2, 3)]
    for ratio in (F(1, 7), F(1), F(3), F(10)):
        x, y = ratio/(4*(1+ratio)), 1/(4*(1+ratio))
        matrix = laplacian(4, edges, [x, y, y, x])
        for vector, value in (([1, 1, -1, -1], 2*y),
                              ([1, -1, -1, 1], 2*x),
                              ([1, -1, 1, -1], F(1, 2))):
            assert [sum(matrix[i][j]*vector[j] for j in range(4))
                    for i in range(4)] == [value*v for v in vector]
        assert x*x/(y*y) == ratio**2
    edges, masses = [(0, 1), (2, 1), (2, 3)], [2, 1, 1, 2]
    with localcontext() as context:
        context.prec = 110
        for size, sign in product((100, 1000, 10000), (-1, 0, 1)):
            a, b, c = F(1, 3), F(1, size), F(1, 3)*(1+F(sign, size))
            matrix = laplacian(4, edges, [a, b, c])
            normalized = [[dec(matrix[i][j])/(D(masses[i])*D(masses[j])).sqrt()
                           for j in range(4)] for i in range(4)]
            values = eigenvalues(normalized)
            assert abs(values[0]) < D('1e-85')
            ratios = (values[1]/dec(F(2, 3)*b), 2*values[2], 2*values[3])
            assert max(abs(v-1) for v in ratios) < D(12)/size
    print('PASS: 4 exact cycle spectra and 9 high-precision unequal-mass separated spectra')


def check_scalar_asymptotics():
    for initial in ((1., 1., 1.), (0.1, 0.3, 2.), (2., 0.4, 0.2), (0.3, 2., 0.7)):
        a, b, c = initial
        u, z = sqrt(a*c), sqrt(a/c)
        # Integrate log u, log r, w, and t=tau*u in rescaled time.
        state = [log(u), log(b/u), (z-1)/(z+1), 0.]
        def reduced(state):
            lu, lr, w, t = state
            r, cosh = exp(lr), 1+2*w*w/(1-w*w)
            return [-3*cosh-2*r, -cosh-2*r, -3*w, 1-t*(3*cosh+2*r)]
        observed = []
        for index in range(12000):
            state = rk4(reduced, state, 0.002)
            if index in (5999, 11999):
                lu, lr, w, t = state
                z = (1+w)/(1-w)
                observed.append(exp(lr+(4*log(t)-lu)/3))
                assert max(abs(3*t*z-1), abs(3*t/z-1)) < 2e-5
        assert abs(observed[1]/observed[0]-1) < 2e-5
        if initial == (1., 1., 1.):
            assert abs(observed[1]/(1/9)-1) < 1e-8
    # The explicit symmetric solution verifies the memory coefficient
    # separately at several initial scales, without a bundle PDE solve.
    for rho in (1., 0.1, 0.01, 0.001):
        u0, r0, s = rho*rho, 1., 25.
        denom = 1+2*r0-2*r0*exp(-s)
        b = u0*r0*exp(-4*s)/denom**2
        tau = ((1+2*r0)*(exp(3*s)-1)/3-r0*(exp(2*s)-1))/u0
        coefficient = F(2, 3)*b*tau**(4/3)
        expected = 2/(27*rho**(2/3))
        assert abs(coefficient/expected-1) < 1e-9
    print('PASS: 4 unequal-initial-data rescaled ODE solutions; '
          'explicit coefficient 2/(27 rho^(2/3)) at 4 initial scales')


def main():
    check_block_masses()
    check_unequal_rank_changes()
    check_balance()
    check_graph_coefficients()
    check_scalar_asymptotics()
    print('Scope: finite algebra and model diagnostics only; '
          'the geometric flow and transfer arguments are proved in the paper.')


if __name__ == '__main__':
    main()
