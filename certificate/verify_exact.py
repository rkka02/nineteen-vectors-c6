#!/usr/bin/env python3
"""Exact verifier for the two-line d=6, m=19 phase-retrieval certificate.

Mathematical convention:
    A = [a_1 ... a_m] in C^{d x m}; measurement vectors are columns;
    Phi_A(x)_j = |a_j^dagger x|^2.

The script reconstructs the five quadratic Grassmann-chart equations exactly,
computes a lexicographic Groebner basis, extracts the degree-14 eliminant in
beta_5, and proves that the eliminant has no real root by an exact Sturm chain.
It also constructs one explicit Gaussian-integer frame (up to harmless
nonzero column scalings) and checks rank_R L_A = 19 exactly.

Requirements: Python 3.10+ and SymPy.
"""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import json
import math
from typing import Iterable, Sequence

import sympy as sp


D = 6
M = 19
R_POINTS = tuple(range(-5, 6))
S_POINTS = (-5, -4, -3, -2, -1, 1, 3, 4)
COS_PHI = sp.Rational(-7, 25)
SIN_PHI = sp.Rational(24, 25)
ZETA = COS_PHI + sp.I * SIN_PHI


def primitive_part(poly: sp.Poly) -> sp.Poly:
    """Return the primitive part over QQ with the same sign."""
    return poly.primitive()[1]


def sign_variations(signs: Iterable[int]) -> int:
    nonzero = [s for s in signs if s]
    return sum(a * b < 0 for a, b in zip(nonzero, nonzero[1:]))


def trim(coeffs: Sequence[Fraction]) -> list[Fraction]:
    out = list(coeffs)
    while out and out[0] == 0:
        out.pop(0)
    return out or [Fraction(0)]


def derivative(coeffs: Sequence[Fraction]) -> list[Fraction]:
    coeffs = trim(coeffs)
    degree = len(coeffs) - 1
    if degree <= 0:
        return [Fraction(0)]
    return [coeffs[i] * (degree - i) for i in range(degree)]


def polynomial_remainder(
    dividend: Sequence[Fraction], divisor: Sequence[Fraction]
) -> list[Fraction]:
    """Exact remainder for descending-coefficient polynomials over Q."""
    rem = trim(dividend)
    div = trim(divisor)
    if div == [Fraction(0)]:
        raise ZeroDivisionError("Polynomial division by zero.")
    while rem != [Fraction(0)] and len(rem) >= len(div):
        factor = rem[0] / div[0]
        shift = len(rem) - len(div)
        subtractor = [factor * c for c in div] + [Fraction(0)] * shift
        rem = trim([a - b for a, b in zip(rem, subtractor)])
    return rem


def sturm_variations_fraction(integer_coeffs: Sequence[int]) -> tuple[int, int]:
    """Compute V(-infinity), V(+infinity) using only Fraction arithmetic."""
    p0 = [Fraction(c) for c in integer_coeffs]
    p1 = derivative(p0)
    chain = [trim(p0), trim(p1)]
    while chain[-1] != [Fraction(0)]:
        rem = polynomial_remainder(chain[-2], chain[-1])
        if rem == [Fraction(0)]:
            break
        chain.append([-c for c in rem])

    plus_signs: list[int] = []
    minus_signs: list[int] = []
    for poly in chain:
        lead = poly[0]
        sign = 1 if lead > 0 else -1
        degree = len(poly) - 1
        plus_signs.append(sign)
        minus_signs.append(sign if degree % 2 == 0 else -sign)
    return sign_variations(minus_signs), sign_variations(plus_signs)


def sin_multiples() -> tuple[sp.Rational, ...]:
    """Return sin(k phi), k=1,...,5, exactly."""
    out: list[sp.Rational] = []
    power = sp.Integer(1)
    for _ in range(5):
        power = sp.expand(power * ZETA)
        out.append(sp.Rational(sp.im(power)))
    return tuple(out)


def build_reduced_system() -> tuple[
    tuple[sp.Symbol, ...], list[sp.Expr], list[sp.Rational], tuple[sp.Rational, ...]
]:
    """Build the five exact quadratic equations in the p_12=1 chart."""
    beta4, beta5, alpha5, alpha6, beta6 = sp.symbols(
        "beta4 beta5 alpha5 alpha6 beta6"
    )
    variables = (beta4, beta5, alpha5, alpha6, beta6)
    t = sp.symbols("t")

    target = sp.Poly(
        sp.prod(t - sp.Integer(s) for s in S_POINTS), t, domain=sp.QQ
    )
    g = [sp.Rational(target.nth(k)) for k in range(9)]
    if g[0] == 0:
        raise AssertionError("The target polynomial must have nonzero constant term.")

    s1, s2, s3, s4, s5 = sin_multiples()
    if s1 == 0 or s2 == 0:
        raise AssertionError("This chart reduction requires sin(phi) sin(2phi) != 0.")

    ratios = [sp.cancel(c / g[0]) for c in g]
    beta3 = sp.cancel((s1 / s2) * ratios[1])
    alpha3 = sp.cancel((s3 / s1) * beta4 - ratios[2])
    alpha4 = sp.cancel((s4 * beta5 - s1 * ratios[3]) / s2)

    # Correct coefficient formulas for d_5,...,d_9.
    d5 = s1 * (alpha3 * beta4 - alpha4 * beta3) - s3 * alpha5 + s5 * beta6
    d6 = s2 * (alpha3 * beta5 - alpha5 * beta3) - s4 * alpha6
    d7 = s3 * (alpha3 * beta6 - alpha6 * beta3) + s1 * (
        alpha4 * beta5 - alpha5 * beta4
    )
    d8 = s2 * (alpha4 * beta6 - alpha6 * beta4)
    d9 = s1 * (alpha5 * beta6 - alpha6 * beta5)

    equations: list[sp.Expr] = []
    for n, dn in zip(range(5, 10), (d5, d6, d7, d8, d9)):
        numerator = sp.together(dn - s1 * ratios[n - 1]).as_numer_denom()[0]
        poly = sp.Poly(numerator, *variables, domain=sp.QQ)
        equations.append(primitive_part(poly).as_expr())

    return variables, equations, g, (s1, s2, s3, s4, s5)


def extract_eliminant(
    variables: tuple[sp.Symbol, ...], equations: Sequence[sp.Expr]
) -> tuple[sp.Poly, sp.GroebnerBasis]:
    """Compute the exact beta_5 eliminant via grevlex + FGLM."""
    beta4, beta5, alpha5, alpha6, beta6 = variables
    order = (alpha5, alpha6, beta6, beta4, beta5)
    grevlex = sp.groebner(
        list(equations), *order, order="grevlex", method="f5b"
    )
    if not grevlex.is_zero_dimensional:
        raise AssertionError("Expected a zero-dimensional complex fiber.")
    lex_basis = grevlex.fglm(order="lex")

    candidates = [
        sp.Poly(p.as_expr(), beta5, domain=sp.QQ)
        for p in lex_basis.polys
        if p.as_expr().free_symbols == {beta5}
    ]
    if not candidates:
        raise AssertionError("No univariate beta_5 eliminant was found.")
    eliminant = primitive_part(max(candidates, key=lambda p: p.degree()))
    return eliminant, lex_basis


def hcoord_row(a: sp.Matrix) -> list[sp.Expr]:
    """Exact hcoord row representing Q -> a^dagger Q a."""
    pmat = a * a.conjugate().T
    row: list[sp.Expr] = [sp.simplify(pmat[p, p]) for p in range(D)]
    for p in range(D):
        for q in range(p + 1, D):
            row.extend(
                [
                    sp.simplify(2 * sp.re(pmat[p, q])),
                    sp.simplify(2 * sp.im(pmat[p, q])),
                ]
            )
    return row


def build_explicit_frame() -> sp.Matrix:
    """Build the unscaled rational two-line frame A in C^{6 x 19}."""
    points = [sp.Integer(r) for r in R_POINTS] + [ZETA * s for s in S_POINTS]
    columns = [
        sp.Matrix([sp.conjugate(z) ** k for k in range(D)]) for z in points
    ]
    return sp.Matrix.hstack(*columns)


def main() -> None:
    certificate_path = Path(__file__).with_name("certificate.json")
    certificate = json.loads(certificate_path.read_text(encoding="utf-8"))

    variables, equations, g, sins = build_reduced_system()
    eliminant, lex_basis = extract_eliminant(variables, equations)
    beta5 = variables[1]

    expected_coeffs = [
        sp.Integer(c) for c in certificate["eliminant_coeffs_descending"]
    ]
    if eliminant.all_coeffs() != expected_coeffs:
        raise AssertionError("The reconstructed eliminant does not match certificate.json.")

    # Exact ideal-membership check: the eliminant reduces to zero modulo the
    # Groebner basis generated from the five quadratics.
    _, remainder = lex_basis.reduce(eliminant.as_expr())
    if sp.expand(remainder) != 0:
        raise AssertionError("Eliminant ideal-membership verification failed.")

    if eliminant.degree() != 14:
        raise AssertionError("Expected a degree-14 eliminant.")
    if sp.gcd(eliminant, eliminant.diff()).degree() != 0:
        raise AssertionError("Expected a squarefree eliminant.")

    # SymPy's exact Sturm count.
    sympy_real_count = int(eliminant.count_roots(-sp.oo, sp.oo))

    # Independent pure-Python/Fraction Sturm calculation.
    integer_coeffs = [int(c) for c in eliminant.all_coeffs()]
    variations_minus, variations_plus = sturm_variations_fraction(integer_coeffs)
    fraction_real_count = variations_minus - variations_plus

    if (sympy_real_count, fraction_real_count) != (0, 0):
        raise AssertionError(
            f"Real-root exclusion failed: SymPy={sympy_real_count}, "
            f"Fraction-Sturm={fraction_real_count}."
        )

    # Exact lifted rank for one explicit frame.
    frame = build_explicit_frame()
    lifted = sp.Matrix([hcoord_row(frame[:, j]) for j in range(M)])
    lifted_rank = int(lifted.rank())
    if lifted_rank != M:
        raise AssertionError(f"Expected rank_R L_A = 19, got {lifted_rank}.")

    print("EXACT-CERTIFIED")
    print(f"d={D}, m={M}")
    print(f"zeta={ZETA}")
    print(f"r={R_POINTS}")
    print(f"s={S_POINTS}")
    print(f"g(t) coefficients, ascending={g}")
    print(f"sin(k phi), k=1,...,5={sins}")
    print(f"quadratic equations={len(equations)}")
    print(f"eliminant degree={eliminant.degree()}")
    print("eliminant squarefree=True")
    print(f"Sturm variations: V(-infinity)={variations_minus}, "
          f"V(+infinity)={variations_plus}")
    print("real roots=0")
    print(f"rank_R L_A={lifted_rank}")


if __name__ == "__main__":
    main()
