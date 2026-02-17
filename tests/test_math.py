"""Symbolic verification of derivations in math.qmd."""

from sympy import (
    Rational,
    Symbol,
    integrate,
    simplify,
    symbols,
)


# ── Kinematic equations ─────────────────────────────────────────────


def test_velocity_from_acceleration():
    """Verify: v = u + at (integrate a dt from 0 to t)."""
    a, u, t = symbols("a u t")

    v = u + integrate(a, (t, 0, t))
    assert simplify(v - (u + a * t)) == 0


def test_position_from_velocity():
    """Verify: x = x₀ + ut + ½at² (integrate v dt from 0 to t)."""
    a, u, t, x_0 = symbols("a u t x_0")

    v = u + a * t
    x = x_0 + integrate(v, (t, 0, t))
    expected = x_0 + u * t + Rational(1, 2) * a * t**2

    assert simplify(x - expected) == 0


def test_velocity_squared():
    """Verify: v² = u² + 2ax.

    From v = u + at and x = ut + ½at², eliminate t.
    """
    a, u, v, t = symbols("a u v t")

    # v = u + at → t = (v-u)/a
    t_expr = (v - u) / a

    # x = ut + ½at²
    x_expr = u * t_expr + Rational(1, 2) * a * t_expr**2
    x_simplified = simplify(x_expr)

    # Should equal (v²-u²)/(2a)
    expected_x = (v**2 - u**2) / (2 * a)

    assert simplify(x_simplified - expected_x) == 0


# ── Differential equations ──────────────────────────────────────────


def test_separable_de_exponential():
    """Verify: dx/x = kt dt → x = x₀ e^{kt²/2}."""
    from sympy import exp

    k, t, x_0 = symbols("k t x_0", positive=True)
    x = Symbol("x", positive=True)

    # ∫dx/x from x₀ to x = ∫kt dt from 0 to t
    # ln(x/x₀) = kt²/2
    # x = x₀ e^(kt²/2)
    result = x_0 * exp(k * t**2 / 2)

    # Verify by differentiating: dx/dt should equal k*t*x
    from sympy import diff

    dx_dt = diff(result, t)
    assert simplify(dx_dt - k * t * result) == 0


# ── Spring potential energy (used across chapters) ──────────────────


def test_spring_pe_from_force():
    """Verify: F = -kx → U = ½kx² (integrate from 0 to x_f)."""
    k, x = symbols("k x", positive=True)

    F = -k * x
    U = -integrate(F, (x, 0, x))
    expected = Rational(1, 2) * k * x**2

    assert simplify(U - expected) == 0
