#!/usr/bin/env python3
"""Exact finite checks for the star-graph matrices in the manuscript.

Standard library only. These check algebraic examples, not the HE analysis.
For k=1,...,4 and every factor rank in {1,2,3}, compare the full mass
Laplacian with the reduced matrix, including every characteristic coefficient.
"""

from fractions import Fraction as F
from itertools import combinations, product
from math import prod


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def multiply(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)]
            for row in a]


def determinant(a):
    a = [[F(x) for x in row] for row in a]
    answer = F(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            answer = -answer
        p = a[j][j]
        answer *= p
        for i in range(j + 1, len(a)):
            ratio = a[i][j] / p
            for h in range(j + 1, len(a)):
                a[i][h] -= ratio * a[j][h]
    return answer


def elementary_coefficients(a):
    """Characteristic coefficients via principal minors, no eigenvalue solver."""
    n = len(a)
    return [F(1)] + [sum(determinant([[a[i][j] for j in ids] for i in ids])
                           for ids in combinations(range(n), m))
                        for m in range(1, n + 1)]


def check_ranks(ranks):
    r0, *leaves = ranks
    k, total = len(leaves), sum(ranks)
    mass = [[F(leaves[i] if i == j else 0) - F(leaves[i] * leaves[j], total)
             for j in range(k)] for i in range(k)]
    inverse = [[F(1, leaves[i]) * (i == j) + F(1, r0)
                for j in range(k)] for i in range(k)]
    require(multiply(mass, inverse) == [[F(i == j) for j in range(k)]
                                      for i in range(k)], (ranks, "inverse"))
    # Similar to S(epsilon), but rational and with epsilon factors on columns.
    coefficient = [[inverse[i][j] * leaves[j] for j in range(k)]
                   for i in range(k)]
    for size in range(1, k + 1):
        for ids in combinations(range(k), size):
            actual = determinant([[coefficient[i][j] for j in ids] for i in ids])
            expected = F(r0 + sum(leaves[i] for i in ids), r0)
            require(actual == expected, (ranks, ids, "principal coefficient"))
    epsilon = [F(1, 2 ** (j + 1)) for j in range(k)]
    reduced = [[coefficient[i][j] * epsilon[j] for j in range(k)] for i in range(k)]
    full = [[F(0) for _ in ranks] for _ in ranks]
    for j in range(1, k + 1):
        weight = leaves[j - 1] * epsilon[j - 1]
        full[0][0] += weight / r0
        full[0][j] -= weight / r0
        full[j][0] -= weight / ranks[j]
        full[j][j] += weight / ranks[j]
    coeff = elementary_coefficients(reduced)
    require(elementary_coefficients(full) == coeff + [F(0)],
            (ranks, "full and reduced characteristic polynomials"))
    require(coeff[-1] == F(total, r0) * prod(epsilon), (ranks, "determinant"))
    require(coeff[1] == sum((1 + F(rj, r0)) * e for rj, e in zip(leaves, epsilon)),
            (ranks, "trace"))


def main():
    count = 0
    for k in range(1, 5):
        for ranks in product((1, 2, 3), repeat=k + 1):
            check_ranks(ranks)
            count += 1
    require(count == 360, "domain size")
    # X_2: trace 3/2 (a+b), determinant 2ab. The squared root term
    # in the displayed eigenvalues is trace^2/4 - determinant.
    for a, b in product((F(1, 7), F(2), F(11)), repeat=2):
        discriminant = (F(3, 2) * (a + b)) ** 2 / 4 - 2 * a * b
        require(discriminant == (9 * (a - b) ** 2 + 4 * a * b) / 16,
                "X_2 discriminant")
    # X_3: the coefficient for any m-element subset is (3+m)/3.
    # For epsilon=(s^6,s^4,s^2), the leading elementary terms are
    # e_1 ~ (4/3)s^2, e_2 ~ (5/3)s^6, e_3 ~ 2s^12.
    orders = (6, 4, 2)
    leading = []
    for m in (1, 2, 3):
        powers = [sum(orders[i] for i in ids) for ids in combinations(range(3), m)]
        degree = min(powers)
        leading.append((degree, F(3 + m, 3) * powers.count(degree)))
    require(leading == [(2, F(4, 3)), (6, F(5, 3)), (12, F(2))], "X_3 polynomial")
    constants = [leading[2][1] / leading[1][1],
                 leading[1][1] / leading[0][1], leading[0][1]]
    require(constants == [F(6, 5), F(5, 4), F(4, 3)], "X_3 separated constants")
    print(f"PASS: {count} rank tuples; all principal coefficients, inverse, full/reduced spectra.")
    print("PASS: X_2 discriminant and X_3 separated constants 6/5, 5/4, 4/3.")
    print("Scope: finite exact matrix checks; no PDE or priority certification.")


if __name__ == "__main__":
    main()
