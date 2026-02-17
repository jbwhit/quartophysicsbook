"""Symbolic verification of derivations in mechanics.qmd."""

from sympy import (
    Rational,
    cos,
    pi,
    simplify,
    sin,
    solve,
    sqrt,
    symbols,
)


# ── Gravitational PE ──────────────────────────────────────────────


def test_gravitational_pe_from_force():
    """Verify: F = -Gm1m2/r² → U = -Gm1m2/r (integrating from ∞)."""
    from sympy import integrate, oo

    G, m1, m2, r = symbols("G m1 m2 r", positive=True)

    F = -G * m1 * m2 / r**2
    # U = -∫_∞^r F dr  (work done by gravity bringing mass from ∞)
    U = -integrate(F, (r, oo, r))
    expected = -G * m1 * m2 / r

    assert simplify(U - expected) == 0


# ── Spring / Hooke's law ──────────────────────────────────────────


def test_spring_work():
    """Verify: F = -kx → W = ½kx² (integrating from 0 to x)."""
    from sympy import integrate

    k, x = symbols("k x", positive=True)

    F = k * x  # magnitude of applied force
    W = integrate(F, (x, 0, x))
    expected = Rational(1, 2) * k * x**2

    assert simplify(W - expected) == 0


# ── Elastic collisions ───────────────────────────────────────────


def test_elastic_collision_final_velocities():
    """Verify elastic collision formulas conserve both momentum and KE."""
    mA, mB, vAi, vBi = symbols("m_A m_B v_Ai v_Bi")

    vAf = (mA - mB) / (mA + mB) * vAi + 2 * mB / (mA + mB) * vBi
    vBf = 2 * mA / (mA + mB) * vAi + (mB - mA) / (mA + mB) * vBi

    # Check momentum conservation
    p_initial = mA * vAi + mB * vBi
    p_final = mA * vAf + mB * vBf
    assert simplify(p_initial - p_final) == 0

    # Check KE conservation
    ke_initial = Rational(1, 2) * mA * vAi**2 + Rational(1, 2) * mB * vBi**2
    ke_final = Rational(1, 2) * mA * vAf**2 + Rational(1, 2) * mB * vBf**2
    assert simplify(ke_initial - ke_final) == 0


def test_elastic_collision_equal_masses_swap():
    """If mA = mB, they swap velocities."""
    m, vAi, vBi = symbols("m v_Ai v_Bi")

    vAf = (m - m) / (m + m) * vAi + 2 * m / (m + m) * vBi
    vBf = 2 * m / (m + m) * vAi + (m - m) / (m + m) * vBi

    assert simplify(vAf - vBi) == 0
    assert simplify(vBf - vAi) == 0


# ── Inelastic collisions ─────────────────────────────────────────


def test_perfectly_inelastic_collision():
    """Verify: vf = (mA*vAi + mB*vBi)/(mA+mB) conserves momentum."""
    mA, mB, vAi, vBi = symbols("m_A m_B v_Ai v_Bi")

    vf = (mA * vAi + mB * vBi) / (mA + mB)

    p_initial = mA * vAi + mB * vBi
    p_final = (mA + mB) * vf

    assert simplify(p_initial - p_final) == 0


# ── Moment of inertia of sphere ───────────────────────────────────


def test_sphere_moment_of_inertia():
    """Verify I = 2/5 MR² by direct integration."""
    from sympy import integrate

    R_val, rho = symbols("R rho", positive=True)

    # Integrate r_perp² dm over sphere using cylindrical shells
    # For a sphere of radius R, at height z, the cross-sectional radius is
    # r_perp = sqrt(R² - z²). Using disk method:
    # dI = ½ dm * r_perp² where dm = ρπr_perp²dz (disk of radius r_perp)
    # Actually for a disk: I_disk = ½ m r² = ½ (ρπr²dz) r²
    from sympy import Symbol

    z = Symbol("z")
    r_perp_sq = R_val**2 - z**2
    # Moment of inertia of thin disk about its axis: ½mr² where m = ρπr²dz
    dI = Rational(1, 2) * rho * pi * r_perp_sq**2

    I_total = integrate(dI, (z, -R_val, R_val))

    M = rho * Rational(4, 3) * pi * R_val**3
    expected = Rational(2, 5) * M * R_val**2

    assert simplify(I_total - expected) == 0


# ── Rolling without slipping ─────────────────────────────────────


def test_rolling_sphere_acceleration():
    """Verify: a = 5/7 g sinθ for a solid sphere rolling without slipping."""
    M, g, theta, a = symbols("M g theta a", positive=True)
    R = symbols("R", positive=True)

    I = Rational(2, 5) * M * R**2

    # From Mg sinθ - f = Ma and fR = Ia/R → f = Ia/R²
    # Mg sinθ - Ia/R² = Ma
    # a(M + I/R²) = Mg sinθ
    eq = M * g * sin(theta) - (M + I / R**2) * a
    a_solved = solve(eq, a)[0]

    expected = Rational(5, 7) * g * sin(theta)

    assert simplify(a_solved - expected) == 0


# ── Precession ────────────────────────────────────────────────────


def test_precession_rate():
    """Verify: ωp = rmg/(Iω) = rmg/L."""
    r, m, g, I, omega = symbols("r m g I omega", positive=True)

    L = I * omega
    omega_p = r * m * g / L
    expected = r * m * g / (I * omega)

    assert simplify(omega_p - expected) == 0


# ── Escape velocity ───────────────────────────────────────────────


def test_escape_velocity():
    """Verify: ½mv² - GMm/r = 0 → v = √(2GM/r)."""
    G, M, m, r, v = symbols("G M m r v", positive=True)

    energy_eq = Rational(1, 2) * m * v**2 - G * M * m / r
    v_solutions = solve(energy_eq, v)
    expected = sqrt(2 * G * M / r)

    assert any(simplify(sol - expected) == 0 for sol in v_solutions)


# ── Kepler's third law ────────────────────────────────────────────


def test_kepler_third_law():
    """Verify: GMm/r² = mω²r → T²/r³ = 4π²/(GM)."""
    G, M, m, r, omega = symbols("G M m r omega", positive=True)

    # Gravitational = centripetal
    eq = G * M * m / r**2 - m * r * omega**2
    omega_solved = solve(eq, omega)

    # T = 2π/ω → T² = 4π²/ω²
    # Take positive solution
    omega_pos = [s for s in omega_solved if s.is_positive][0]
    T_sq = (2 * pi / omega_pos) ** 2

    ratio = simplify(T_sq / r**3)
    expected = 4 * pi**2 / (G * M)

    assert simplify(ratio - expected) == 0


# ── Hydrostatic pressure ─────────────────────────────────────────


def test_hydrostatic_pressure():
    """Verify: P = ρgh + P₀ from F = ρVg."""
    rho, g, h, A, P0 = symbols("rho g h A P_0", positive=True)

    F = rho * A * h * g  # weight of fluid column
    P = F / A + P0

    expected = rho * g * h + P0
    assert simplify(P - expected) == 0


# ── Continuity equation ──────────────────────────────────────────


def test_continuity_equation():
    """Verify: A₁v₁ = A₂v₂ (mass conservation for incompressible flow)."""
    A1, v1, A2, v2, rho = symbols("A1 v1 A2 v2 rho", positive=True)

    # Mass flow rate: ρA₁v₁ = ρA₂v₂ → A₁v₁ = A₂v₂
    assert simplify(rho * A1 * v1 - rho * A2 * v2) == simplify(
        rho * (A1 * v1 - A2 * v2)
    )
    # Given A1*v1 = A2*v2, solve for v2:
    v2_solved = solve(A1 * v1 - A2 * v2, v2)[0]
    assert simplify(v2_solved - A1 * v1 / A2) == 0


# ── Poiseuille's law ─────────────────────────────────────────────


def test_poiseuille_r4_dependence():
    """Verify Poiseuille: doubling radius → 16× flow rate."""
    r, eta, dP, dL = symbols("r eta dP dL", positive=True)

    Q = pi * r**4 * dP / (8 * eta * dL)
    Q_doubled = Q.subs(r, 2 * r)

    ratio = simplify(Q_doubled / Q)
    assert ratio == 16


# ── Limiting cases ────────────────────────────────────────────────


def test_escape_velocity_limiting_cases():
    """Limiting cases for escape velocity."""
    G, M, r = symbols("G M r", positive=True)

    v_esc = sqrt(2 * G * M / r)

    # As M → 0, v_esc → 0
    assert v_esc.subs(M, 0) == 0

    # As r → ∞, v_esc → 0
    from sympy import limit, oo

    assert limit(v_esc, r, oo) == 0


def test_elastic_collision_heavy_target_at_rest():
    """If M >> m and target at rest: projectile bounces back, target barely moves."""
    M, m, vAi = symbols("M_big m v_Ai", positive=True)

    # B at rest (vBi = 0), A = m (light), B = M (heavy)
    vAf = (m - M) / (m + M) * vAi
    vBf = 2 * m / (m + M) * vAi

    from sympy import limit, oo

    # As M → ∞: vAf → -vAi (bounces back)
    assert limit(vAf, M, oo) == -vAi

    # As M → ∞: vBf → 0 (heavy target doesn't move)
    assert limit(vBf, M, oo) == 0
