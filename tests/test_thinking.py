"""Symbolic verification of derivations in thinking.qmd."""

from sympy import Rational, cos, pi, simplify, sin, solve, sqrt, symbols


# ── Example Problem 1: Spring compression ──────────────────────────


def test_spring_compression_energy_conservation():
    """Verify: mg(H-h) = ½kx² → x = √(2mg(H-h)/k)."""
    m, g, H, h, k, x = symbols("m g H h k x", positive=True)

    # Conservation of energy: GPE lost = elastic PE gained
    energy_eq = m * g * (H - h) - Rational(1, 2) * k * x**2

    solutions = solve(energy_eq, x)
    expected = sqrt(2 * m * g * (H - h) / k)

    assert any(simplify(sol - expected) == 0 for sol in solutions)


def test_spring_compression_limiting_cases():
    """Verify limiting cases from the book."""
    m, g, H, h, k = symbols("m g H h k", positive=True)
    x = sqrt(2 * m * g * (H - h) / k)

    # H → h: x → 0
    assert x.subs(H, h) == 0

    # m → 0: x → 0
    assert x.subs(m, 0) == 0


# ── Example Problem 2: Piano on incline ────────────────────────────
#
# NOTE: The book's boxed result has an INTENTIONAL sign error, left in
# as a pedagogical example of how limiting-case checks catch mistakes.
# (See thinking.qmd line 177: "This working includes an error caught
# during the limiting cases check, and left in to show how this works.")
#
# The book gives:  a = g(Msinθ + μ_k Mcosθ - m) / (M + m)
# Correct answer:  a = g(m - Msinθ - μ_k Mcosθ) / (M + m)
#
# Convention: positive a = man descends, piano ascends.
# Man:   mg - T = ma  →  T = m(g - a)
# Piano: T - Mgsinθ - μ_k Mgcosθ = Ma
# ────────────────────────────────────────────────────────────────────


def test_piano_incline_acceleration():
    """Verify correct derivation for the piano-on-incline problem."""
    M, m_man, g, theta, mu_k, a = symbols("M m g theta mu_k a", positive=True)

    # Man (descending): mg - T = ma → T = m(g - a)
    T_expr = m_man * (g - a)

    # Piano (ascending along incline): T - Mgsinθ - μ_k Mgcosθ = Ma
    piano_eq = T_expr - M * g * sin(theta) - mu_k * M * g * cos(theta) - M * a

    a_solutions = solve(piano_eq, a)
    assert len(a_solutions) == 1

    a_solved = a_solutions[0]
    expected = g * (m_man - M * sin(theta) - mu_k * M * cos(theta)) / (M + m_man)

    assert simplify(a_solved - expected) == 0


def test_piano_incline_limiting_atwood():
    """At θ=90°, μ_k=0 this reduces to Atwood machine: a = g(m-M)/(M+m)."""
    M, m_man, g = symbols("M m g", positive=True)

    a = g * (m_man - M * sin(pi / 2) - 0 * M * cos(pi / 2)) / (M + m_man)
    expected_atwood = g * (m_man - M) / (M + m_man)

    assert simplify(a - expected_atwood) == 0


def test_piano_incline_flat_no_friction():
    """At θ=0, μ_k=0: a = gm/(M+m) — man's weight accelerates system."""
    M, m_man, g = symbols("M m g", positive=True)

    a = g * (m_man - M * sin(0) - 0 * M * cos(0)) / (M + m_man)
    expected = g * m_man / (M + m_man)

    assert simplify(a - expected) == 0


def test_book_sign_error_is_negative_of_correct():
    """Confirm the book's boxed result is exactly the negative of the correct one.

    This documents the intentional pedagogical error in thinking.qmd.
    """
    M, m_man, g, theta, mu_k = symbols("M m g theta mu_k", positive=True)

    correct = g * (m_man - M * sin(theta) - mu_k * M * cos(theta)) / (M + m_man)
    book = g * (M * sin(theta) + mu_k * M * cos(theta) - m_man) / (M + m_man)

    assert simplify(correct + book) == 0
