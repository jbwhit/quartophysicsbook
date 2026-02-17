"""Symbolic verification of derivations in electrostatics.qmd."""

from sympy import (
    Rational,
    cos,
    integrate,
    limit,
    oo,
    pi,
    simplify,
    sin,
    solve,
    sqrt,
    symbols,
)


# ── Coulomb's law ───────────────────────────────────────────


def test_coulomb_inverse_square():
    """Verify: doubling distance quarters the force."""
    q1, q2, r, eps0 = symbols("q1 q2 r epsilon_0", positive=True)

    F = q1 * q2 / (4 * pi * eps0 * r**2)
    F_double = F.subs(r, 2 * r)

    assert simplify(F_double / F) == Rational(1, 4)


def test_coulomb_superposition():
    """Verify: force from two equal charges at ±d on axis is zero at origin."""
    q, Q, d, eps0 = symbols("q Q d epsilon_0", positive=True)

    # Charge Q at +d and +Q at -d, force on q at origin
    F_right = q * Q / (4 * pi * eps0 * d**2)  # points left (toward -d)
    F_left = q * Q / (4 * pi * eps0 * d**2)  # points right (toward +d)

    # Equal and opposite → net force = 0
    assert simplify(F_right - F_left) == 0


# ── Electric field ──────────────────────────────────────────


def test_electric_field_point_charge():
    """Verify: E = Q/(4πε₀r²) from F = qE."""
    q, Q, r, eps0 = symbols("q Q r epsilon_0", positive=True)

    F = q * Q / (4 * pi * eps0 * r**2)
    E = F / q
    expected = Q / (4 * pi * eps0 * r**2)

    assert simplify(E - expected) == 0


# ── Gauss's law applications ────────────────────────────────


def test_gauss_point_charge():
    """Verify: Gauss's law recovers Coulomb's law for a point charge."""
    q, r, eps0 = symbols("q r epsilon_0", positive=True)

    # E * 4πr² = q/ε₀ → E = q/(4πε₀r²)
    E = solve(symbols("E") * 4 * pi * r**2 - q / eps0, symbols("E"))
    expected = q / (4 * pi * eps0 * r**2)

    assert any(simplify(sol - expected) == 0 for sol in E)


def test_gauss_line_charge():
    """Verify: E = λ/(2πε₀r) for infinite line charge."""
    lam, r, L, eps0 = symbols("lambda r L epsilon_0", positive=True)

    # E * 2πrL = λL/ε₀
    E_val = solve(symbols("E") * 2 * pi * r * L - lam * L / eps0, symbols("E"))
    expected = lam / (2 * pi * eps0 * r)

    assert any(simplify(sol - expected) == 0 for sol in E_val)


def test_gauss_line_charge_alt_form():
    """Verify: λ/(2πε₀r) = (1/4πε₀)(2λ/r)."""
    lam, r, eps0 = symbols("lambda r epsilon_0", positive=True)

    form1 = lam / (2 * pi * eps0 * r)
    form2 = 2 * lam / (4 * pi * eps0 * r)

    assert simplify(form1 - form2) == 0


def test_gauss_conducting_sphere_inside():
    """Inside a conducting sphere: E = 0 (no enclosed charge)."""
    # q_enc = 0 inside conductor → E = 0
    # This is a conceptual test — verify q_enc = 0 → E = 0
    eps0 = symbols("epsilon_0", positive=True)
    r = symbols("r", positive=True)

    q_enc = 0
    E = q_enc / (4 * pi * eps0 * r**2)
    assert E == 0


# ── Electric potential ──────────────────────────────────────


def test_potential_point_charge():
    """Verify: V = Q/(4πε₀r) by integrating E."""
    Q, r, eps0 = symbols("Q r epsilon_0", positive=True)
    r_var = symbols("r_var", positive=True)

    E = Q / (4 * pi * eps0 * r_var**2)

    # V = -∫_∞^r E dr' = ∫_r^∞ E dr'
    V = integrate(E, (r_var, r, oo))
    expected = Q / (4 * pi * eps0 * r)

    assert simplify(V - expected) == 0


def test_potential_superposition():
    """Verify: potential from two charges is sum of individual potentials."""
    q1, q2, r1, r2, eps0 = symbols("q1 q2 r1 r2 epsilon_0", positive=True)

    V1 = q1 / (4 * pi * eps0 * r1)
    V2 = q2 / (4 * pi * eps0 * r2)
    V_total = V1 + V2
    expected = (q1 * r2 + q2 * r1) / (4 * pi * eps0 * r1 * r2)

    assert simplify(V_total - expected) == 0


def test_field_from_potential_1d():
    """Verify: E = -dV/dr recovers Coulomb field from point charge potential."""
    from sympy import diff

    Q, r, eps0 = symbols("Q r epsilon_0", positive=True)

    V = Q / (4 * pi * eps0 * r)
    E = -diff(V, r)
    expected = Q / (4 * pi * eps0 * r**2)

    assert simplify(E - expected) == 0


# ── Capacitors ──────────────────────────────────────────────


def test_parallel_plate_capacitance():
    """Verify: C = ε₀A/d for parallel plates."""
    eps0, A, d = symbols("epsilon_0 A d", positive=True)
    sigma = symbols("sigma", positive=True)

    E = sigma / eps0
    V = E * d
    q = sigma * A
    C = q / V
    expected = eps0 * A / d

    assert simplify(C - expected) == 0


def test_capacitor_energy():
    """Verify: U = ½Q²/C = ½CV²."""
    Q, C, V = symbols("Q C V", positive=True)

    U1 = Rational(1, 2) * Q**2 / C
    U2 = Rational(1, 2) * C * V**2

    # With Q = CV: U1 should equal U2
    assert simplify(U1.subs(Q, C * V) - U2) == 0


def test_capacitor_energy_from_integration():
    """Verify: U = ½Q²/C by integrating q/C dq from 0 to Q."""
    q, Q, C = symbols("q Q C", positive=True)

    U = integrate(q / C, (q, 0, Q))
    expected = Rational(1, 2) * Q**2 / C

    assert simplify(U - expected) == 0


# ── Dielectrics ─────────────────────────────────────────────


def test_dielectric_capacitance():
    """Verify: C = κC₀ with dielectric."""
    eps0, A, d, kappa = symbols("epsilon_0 A d kappa", positive=True)

    C0 = eps0 * A / d
    C_dielectric = kappa * eps0 * A / d

    assert simplify(C_dielectric / C0 - kappa) == 0


def test_dielectric_field_reduction():
    """Verify: E_total = E₀/κ inside a linear dielectric."""
    E0, kappa = symbols("E_0 kappa", positive=True)

    E_total = E0 / kappa

    # For κ = 1 (vacuum): E_total = E₀
    assert simplify(E_total.subs(kappa, 1) - E0) == 0

    # For κ = 2: E_total = E₀/2
    assert simplify(E_total.subs(kappa, 2) - E0 / 2) == 0


# ── Dipoles ─────────────────────────────────────────────────


def test_dipole_torque_maximum():
    """Verify: torque is maximum when p ⊥ E (θ = π/2)."""
    p, E = symbols("p E", positive=True)
    theta = symbols("theta")

    tau = p * E * sin(theta)

    # Maximum at θ = π/2
    assert simplify(tau.subs(theta, pi / 2) - p * E) == 0

    # Zero at θ = 0 (aligned)
    assert tau.subs(theta, 0) == 0


def test_dipole_energy_extrema():
    """Verify: U = -pE cosθ has min at θ=0 and max at θ=π."""
    p, E = symbols("p E", positive=True)
    theta = symbols("theta")

    U = -p * E * cos(theta)

    # Minimum at θ = 0: U = -pE
    assert simplify(U.subs(theta, 0) - (-p * E)) == 0

    # Maximum at θ = π: U = +pE
    assert simplify(U.subs(theta, pi) - p * E) == 0


def test_dipole_potential():
    """Verify: V = p cosθ / (4πε₀r²) falls off as 1/r²."""
    p, r, eps0 = symbols("p r epsilon_0", positive=True)
    theta = symbols("theta")

    V = p * cos(theta) / (4 * pi * eps0 * r**2)

    # Doubling r → V/4
    V_double = V.subs(r, 2 * r)
    assert simplify(V_double / V) == Rational(1, 4)


def test_dipole_field_on_axis():
    """Verify: on-axis (θ=0) dipole field is E = 2p/(4πε₀r³)."""
    p, r, eps0 = symbols("p r epsilon_0", positive=True)

    # E = p/(4πε₀r³)(2cosθ r̂ + sinθ θ̂), at θ = 0:
    E_r = p * 2 * cos(0) / (4 * pi * eps0 * r**3)
    E_theta = p * sin(0) / (4 * pi * eps0 * r**3)

    expected_r = 2 * p / (4 * pi * eps0 * r**3)
    assert simplify(E_r - expected_r) == 0
    assert E_theta == 0


def test_dipole_field_equatorial():
    """Verify: equatorial (θ=π/2) dipole field is E = p/(4πε₀r³)."""
    p, r, eps0 = symbols("p r epsilon_0", positive=True)

    E_r = p * 2 * cos(pi / 2) / (4 * pi * eps0 * r**3)
    E_theta = p * sin(pi / 2) / (4 * pi * eps0 * r**3)

    # Radial component vanishes
    assert E_r == 0

    # θ component
    expected = p / (4 * pi * eps0 * r**3)
    assert simplify(E_theta - expected) == 0


def test_dipole_field_from_potential():
    """Verify: E = -∇V gives the correct dipole field components."""
    from sympy import diff

    p, r, eps0 = symbols("p r epsilon_0", positive=True)
    theta = symbols("theta")

    V = p * cos(theta) / (4 * pi * eps0 * r**2)

    # In spherical: E_r = -∂V/∂r, E_θ = -(1/r)∂V/∂θ
    E_r = -diff(V, r)
    E_theta = -diff(V, theta) / r

    expected_r = 2 * p * cos(theta) / (4 * pi * eps0 * r**3)
    expected_theta = p * sin(theta) / (4 * pi * eps0 * r**3)

    assert simplify(E_r - expected_r) == 0
    assert simplify(E_theta - expected_theta) == 0


# ── Limiting cases ──────────────────────────────────────────


def test_coulomb_force_at_infinity():
    """As r → ∞, F → 0."""
    q1, q2, r, eps0 = symbols("q1 q2 r epsilon_0", positive=True)

    F = q1 * q2 / (4 * pi * eps0 * r**2)
    assert limit(F, r, oo) == 0


def test_potential_at_infinity():
    """As r → ∞, V → 0."""
    Q, r, eps0 = symbols("Q r epsilon_0", positive=True)

    V = Q / (4 * pi * eps0 * r)
    assert limit(V, r, oo) == 0


def test_capacitor_energy_zero_charge():
    """No charge → no energy."""
    C = symbols("C", positive=True)

    U = Rational(1, 2) * 0**2 / C
    assert U == 0


def test_dipole_potential_at_90_degrees():
    """Dipole potential vanishes in equatorial plane (θ = π/2)."""
    p, r, eps0 = symbols("p r epsilon_0", positive=True)

    V = p * cos(pi / 2) / (4 * pi * eps0 * r**2)
    assert V == 0
