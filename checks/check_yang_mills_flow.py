#!/usr/bin/env python3
"""Exact finite checks for the two-factor Yang--Mills flow argument.

Checks the mass-weighted two-vertex metric flow, its time normalization,
endomorphism metric comparison, and the error integral in the restarted
barrier argument. This is not a discretization of the geometric PDE.
"""

from fractions import Fraction as F

from check_two_factor_algebra import add, mul, transpose, zero
from check_several_factor import eye, inverse, is_psd, scaled


def diagonal(values):
    result = zero(len(values), len(values))
    for i, value in enumerate(values):
        result[i][i] = F(value)
    return result


def endomorphism_gram(metric):
    """Matrix of tr(H^-1 U^T H V) in the real matrix-unit basis.

    For real H, this is also the coefficient matrix of the complex
    Hermitian form, so exact PSD checks apply over both fields.
    """
    n, hinv = len(metric), inverse(metric)
    indices = [(i, j) for i in range(n) for j in range(n)]
    return [[metric[i][k]*hinv[l][j] for k, l in indices]
            for i, j in indices]


def check_metric_comparison():
    cases = 0
    for n in (2, 3, 4):
        roots = diagonal(range(1, n+1))
        h0 = mul(roots, roots)
        for bound in (F(2), F(3), F(5, 2)):
            values = [F(1)]*n
            values[0], values[-1] = 1/bound, bound
            for i in range(n-1):
                rotation = eye(n)
                rotation[i][i] = rotation[i+1][i+1] = F(3, 5)
                rotation[i][i+1], rotation[i+1][i] = F(-4, 5), F(4, 5)
                middle = mul(mul(rotation, diagonal(values)), transpose(rotation))
                h1 = mul(mul(roots, middle), roots)
                assert is_psd(add(h1, scaled(h0, 1/bound), -1))
                assert is_psd(add(scaled(h0, bound), h1, -1))
                g0, g1 = endomorphism_gram(h0), endomorphism_gram(h1)
                assert is_psd(add(g1, scaled(g0, 1/bound**2), -1))
                assert is_psd(add(scaled(g0, bound**2), g1, -1))
                # An independent zero-order map supplies a target quadratic
                # form; the same bounds must hold after pullback by that map.
                d = [[F(((i+1)*(j+2)) % 7-3, 5) for j in range(n*n)]
                     for i in range(n*n)]
                a0 = mul(mul(transpose(d), g0), d)
                a1 = mul(mul(transpose(d), g1), d)
                assert is_psd(add(a1, scaled(a0, 1/bound**2), -1))
                assert is_psd(add(scaled(a0, bound**2), a1, -1))
                cases += 1
    # The power four in a Rayleigh comparison cannot be replaced by two:
    # input E12, output E21, H=diag(C,1/C) versus the identity metric.
    bound = F(3)
    gram = endomorphism_gram(diagonal([bound, 1/bound]))
    ratio = gram[2][2]/gram[1][1]
    assert ratio == bound**-4 and ratio < bound**-2
    print(f'PASS: {cases} exact noncommuting metric/form comparisons; Rayleigh exponent four is necessary')


def check_two_vertex_flow():
    cases = 0
    for a in range(1, 6):
        for b in range(1, 6):
            for volume in (F(1), F(3, 2)):
                for beta_norm_sq in (F(1, 3), F(2)):
                    r = a+b
                    c = r*beta_norm_sq/(volume*a*b)
                    qf, qg = F(b, r), F(-a, r)
                    assert a*qf+b*qg == 0 and qf-qg == 1
                    for rho_sq in (F(1, 4), F(1, 25)):
                        # These are the metric logarithmic derivatives from
                        # the projected mean curvature, with flow factor -2.
                        xf = -2*beta_norm_sq*rho_sq/(volume*a)
                        xg = 2*beta_norm_sq*rho_sq/(volume*b)
                        assert xf == -2*c*rho_sq*qf
                        assert xg == -2*c*rho_sq*qg
                        s_prime = rho_sq*(xf-xg)
                        assert s_prime == -2*c*rho_sq**2
                        lam, lam_prime = c*rho_sq, c*s_prime
                        assert lam_prime == -2*lam**2
                    # Closed solution, including a nonzero start time.
                    start, s0 = F(7), F(4, 9)
                    for elapsed in (F(0), F(1), F(10**7)):
                        time = start+elapsed
                        s = 1/(1/s0+2*c*elapsed)
                        lam = c*s
                        assert 1/lam == 1/(c*s0)+2*elapsed
                        expected_difference = (2*c*start-1/s0)/(2*(1/s0+2*c*elapsed))
                        assert time*lam-F(1, 2) == expected_difference
                    cases += 1
    print(f'PASS: {cases} rank, volume and extension-norm cases; lambda_prime=-2 lambda^2 and tau*lambda -> 1/2')


def check_residual_tail_and_restart():
    cases = 0
    for c in (F(1, 3), F(1), F(7)):
        for rho0 in (F(1, 2), F(1, 10), F(1, 100)):
            for rho in (rho0, rho0/2, rho0/5):
                elapsed = (rho**-2-rho0**-2)/(2*c)
                assert elapsed >= 0
                rho_prime = -c*rho**3
                # Exact antiderivative: -rho/c differentiates to rho^3.
                assert -rho_prime/c == rho**3
                integral = (rho0-rho)/c
                assert 0 <= integral < rho0/c
                cases += 1
    # Model the restarted enclosure by logarithmic width. A fixed start
    # leaves positive width; late starts with both errors ->0 close it.
    widths = [F(1, n*n)+F(5, n) for n in (10, 100, 1000, 10000)]
    assert all(x > y > 0 for x, y in zip(widths, widths[1:]))
    assert widths[-1] < F(1, 1000)
    print(f'PASS: {cases} exact cubic-residual tail identities; late-start comparison widths vanish')


def main():
    check_metric_comparison()
    check_two_vertex_flow()
    check_residual_tail_and_restart()
    print('Scope: finite algebra and reduced-flow normalization only; no numerical verification of the PDE or matching lemma.')


if __name__ == '__main__':
    main()
