# Exact matrix computations

Run `python3 check_rank_matrix.py`. Python's standard library is sufficient.

The program uses rational arithmetic to check the reduced mass inverse,
every principal-minor coefficient, and equality of the full and reduced
nonzero characteristic polynomials for all 360 rank tuples in `{1,2,3}`
with one to four quotient factors. It also checks the determinant, trace,
two-leaf discriminant, and the separated toric constants `6/5, 5/4, 4/3`.

The general identities and the analytic estimates are proved in the
[manuscript](../paper/main.pdf). These finite computations do not solve
the Hermitian–Einstein equation or compute the Green operator.
