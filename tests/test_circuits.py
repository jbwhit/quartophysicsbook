"""Symbolic verification of derivations in circuits.qmd."""

from sympy import (
    I as j,
    Rational,
    cos,
    exp,
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


# ── Ohm's law & power ─────────────────────────────────────────


def test_ohms_law():
    """Verify: V = IR and the three power forms are equivalent."""
    V, I_sym, R = symbols("V I R", positive=True)

    # P = IV = I²R = V²/R (all equivalent via V = IR)
    P1 = I_sym * V
    P2 = I_sym**2 * R
    P3 = V**2 / R

    # Substitute V = IR into P1
    assert simplify(P1.subs(V, I_sym * R) - P2) == 0
    # Substitute I = V/R into P1
    assert simplify(P1.subs(I_sym, V / R) - P3) == 0


def test_power_dissipation():
    """Verify: doubling voltage quadruples power (at constant R)."""
    V, R = symbols("V R", positive=True)

    P = V**2 / R
    P_double = P.subs(V, 2 * V)

    assert simplify(P_double / P) == 4


# ── Resistor combinations ─────────────────────────────────────


def test_resistors_in_series():
    """Verify: series resistance R = R₁ + R₂."""
    R1, R2, V = symbols("R1 R2 V", positive=True)

    I_total = V / (R1 + R2)
    V1 = I_total * R1
    V2 = I_total * R2

    # Voltages add to total
    assert simplify(V1 + V2 - V) == 0


def test_resistors_in_parallel():
    """Verify: 1/R_par = 1/R₁ + 1/R₂ → R_par = R₁R₂/(R₁+R₂)."""
    R1, R2 = symbols("R1 R2", positive=True)

    R_par = 1 / (1 / R1 + 1 / R2)
    expected = R1 * R2 / (R1 + R2)

    assert simplify(R_par - expected) == 0


def test_equal_parallel_resistors():
    """Two equal resistors in parallel: R_par = R/2."""
    R = symbols("R", positive=True)

    R_par = 1 / (1 / R + 1 / R)
    assert simplify(R_par - R / 2) == 0


# ── Capacitor & inductor combinations ─────────────────────────


def test_capacitors_in_series():
    """Verify: 1/C_s = 1/C₁ + 1/C₂ (opposite of resistors)."""
    C1, C2 = symbols("C1 C2", positive=True)

    C_s = 1 / (1 / C1 + 1 / C2)
    expected = C1 * C2 / (C1 + C2)

    assert simplify(C_s - expected) == 0


def test_capacitors_in_parallel():
    """Verify: C_par = C₁ + C₂."""
    C1, C2, V = symbols("C1 C2 V", positive=True)

    q1 = C1 * V
    q2 = C2 * V
    q_total = q1 + q2
    C_par = q_total / V

    assert simplify(C_par - (C1 + C2)) == 0


def test_inductors_in_series():
    """Verify: L_s = L₁ + L₂ (like resistors)."""
    L1, L2 = symbols("L1 L2", positive=True)

    L_s = L1 + L2

    # Doubling identical inductors doubles inductance
    assert simplify(L_s.subs(L2, L1) / L1) == 2


def test_inductors_in_parallel():
    """Verify: 1/L_par = 1/L₁ + 1/L₂ (like resistors)."""
    L1, L2 = symbols("L1 L2", positive=True)

    L_par = 1 / (1 / L1 + 1 / L2)
    expected = L1 * L2 / (L1 + L2)

    assert simplify(L_par - expected) == 0


# ── RC circuit ─────────────────────────────────────────────────


def test_rc_charging():
    """Verify: q(t) = CV(1 - e^{-t/RC}) satisfies the RC ODE."""
    from sympy import diff

    R, C, V_s, t = symbols("R C V t", positive=True)

    q = C * V_s * (1 - exp(-t / (R * C)))
    I_sym = diff(q, t)

    # ODE: V = IR + q/C → V - IR - q/C = 0
    residual = V_s - I_sym * R - q / C
    assert simplify(residual) == 0


def test_rc_time_constant():
    """At t = RC, charge is (1 - 1/e) ≈ 63% of final value."""
    R, C, V_s = symbols("R C V", positive=True)

    q_final = C * V_s
    q_at_tau = C * V_s * (1 - exp(-1))
    ratio = q_at_tau / q_final

    # 1 - 1/e ≈ 0.632
    assert simplify(ratio - (1 - exp(-1))) == 0


def test_rc_initial_and_final():
    """RC charging: q(0) = 0 and q(∞) = CV."""
    R, C, V_s, t = symbols("R C V t", positive=True)

    q = C * V_s * (1 - exp(-t / (R * C)))

    assert q.subs(t, 0) == 0
    assert limit(q, t, oo) == C * V_s


# ── RL circuit ─────────────────────────────────────────────────


def test_rl_charging():
    """Verify: I(t) = (V/R)(1 - e^{-Rt/L}) satisfies the RL ODE."""
    from sympy import diff

    R, L, V_s, t = symbols("R L V t", positive=True)

    I_sym = (V_s / R) * (1 - exp(-R * t / L))
    dI_dt = diff(I_sym, t)

    # ODE: V = IR + L dI/dt
    residual = V_s - I_sym * R - L * dI_dt
    assert simplify(residual) == 0


def test_rl_initial_and_final():
    """RL charging: I(0) = 0 and I(∞) = V/R."""
    R, L, V_s, t = symbols("R L V t", positive=True)

    I_sym = (V_s / R) * (1 - exp(-R * t / L))

    assert I_sym.subs(t, 0) == 0
    assert simplify(limit(I_sym, t, oo) - V_s / R) == 0


# ── LC circuit ─────────────────────────────────────────────────


def test_lc_frequency():
    """Verify: ω = 1/√(LC) from d²q/dt² + q/(LC) = 0."""
    L, C, omega = symbols("L C omega", positive=True)

    # Characteristic equation: omega² = 1/(LC)
    solutions = solve(omega**2 - 1 / (L * C), omega)
    expected = 1 / sqrt(L * C)

    assert any(simplify(sol - expected) == 0 for sol in solutions)


def test_lc_energy_conservation():
    """Verify: total energy in LC circuit is constant."""
    L, C, q0, omega, t = symbols("L C q_0 omega t", real=True)

    q = q0 * cos(omega * t)
    I_sym = -q0 * omega * sin(omega * t)

    U_C = Rational(1, 2) * q**2 / C
    U_L = Rational(1, 2) * L * I_sym**2
    E_total = U_C + U_L

    # Substitute omega² = 1/(LC)
    E_sub = E_total.subs(omega**2, 1 / (L * C))
    E_simplified = simplify(E_sub)

    expected = Rational(1, 2) * q0**2 / C
    assert simplify(E_simplified - expected) == 0


# ── Biot–Savart & Ampère ──────────────────────────────────────


def test_wire_field_inverse_distance():
    """Verify: B = μ₀I/(2πd) — doubling distance halves field."""
    mu0, I_sym, d = symbols("mu_0 I d", positive=True)

    B = mu0 * I_sym / (2 * pi * d)
    B_double = B.subs(d, 2 * d)

    assert simplify(B_double / B) == Rational(1, 2)


def test_parallel_conductors_force():
    """Verify: F/L = μ₀I₁I₂/(2πd)."""
    mu0, I1, I2, d = symbols("mu_0 I1 I2 d", positive=True)

    F_per_L = mu0 * I1 * I2 / (2 * pi * d)

    # Doubling separation halves force
    F_double_d = F_per_L.subs(d, 2 * d)
    assert simplify(F_double_d / F_per_L) == Rational(1, 2)


def test_solenoid_field():
    """Verify: B = μ₀nI doesn't depend on radius."""
    mu0, n, I_sym, R = symbols("mu_0 n I R", positive=True)

    B = mu0 * n * I_sym

    # No R dependence
    from sympy import diff

    assert diff(B, R) == 0


# ── Solenoid inductance ────────────────────────────────────────


def test_solenoid_inductance():
    """Verify: L = μ₀n²πR²ℓ for a solenoid."""
    mu0, n, R, ell = symbols("mu_0 n R ell", positive=True)

    B = mu0 * n  # B per unit current (B = μ₀nI, factor out I)
    flux_per_turn = B * pi * R**2
    N = n * ell
    L = N * flux_per_turn
    expected = mu0 * n**2 * pi * R**2 * ell

    assert simplify(L - expected) == 0


# ── Magnetic energy density ────────────────────────────────────


def test_magnetic_energy_density():
    """Verify: u_B = B²/(2μ₀) from U = ½LI² for a solenoid."""
    mu0, n, I_sym, R, ell = symbols("mu_0 n I R ell", positive=True)

    L = mu0 * n**2 * pi * R**2 * ell
    U = Rational(1, 2) * L * I_sym**2
    volume = pi * R**2 * ell
    u = U / volume

    # B = μ₀nI
    B = mu0 * n * I_sym
    expected = B**2 / (2 * mu0)

    assert simplify(u - expected) == 0


def test_electric_energy_density():
    """Verify: u_E = ½ε₀E² from U = ½CV² for parallel plates."""
    eps0, A, d, E = symbols("epsilon_0 A d E", positive=True)

    C = eps0 * A / d
    V = E * d
    U = Rational(1, 2) * C * V**2
    volume = A * d
    u = U / volume

    expected = Rational(1, 2) * eps0 * E**2
    assert simplify(u - expected) == 0


# ── Faraday's law — motional EMF ──────────────────────────────


def test_motional_emf_force():
    """Verify: F = B²h²v₀/R for a wire moving in a field."""
    B, h, v0, R = symbols("B h v_0 R", positive=True)

    emf = B * h * v0
    I_sym = emf / R
    F = I_sym * B * h  # Force on current-carrying wire in B field
    expected = B**2 * h**2 * v0 / R

    assert simplify(F - expected) == 0


# ── AC circuits ────────────────────────────────────────────────


def test_resonance_frequency():
    """Verify: at resonance ωL = 1/(ωC) → ω = 1/√(LC)."""
    L, C, omega = symbols("L C omega", positive=True)

    # Resonance: ωL = 1/(ωC)
    solutions = solve(omega * L - 1 / (omega * C), omega)
    expected = 1 / sqrt(L * C)

    assert any(simplify(sol - expected) == 0 for sol in solutions)


def test_rlc_impedance_at_resonance():
    """At resonance, |Z| = R (purely resistive)."""
    R, L, C = symbols("R L C", positive=True)

    omega_res = 1 / sqrt(L * C)
    Z_imag = omega_res * L - 1 / (omega_res * C)

    assert simplify(Z_imag) == 0


def test_rms_voltage():
    """Verify: V_rms = V₀/√2 for sinusoidal voltage."""
    V0, T, t = symbols("V_0 T t", positive=True)
    omega = symbols("omega", positive=True)

    # V²_rms = (1/T) ∫₀ᵀ V₀²cos²(ωt) dt, with T = 2π/ω
    T_val = 2 * pi / omega
    V_sq_avg = integrate(V0**2 * cos(omega * t) ** 2, (t, 0, T_val)) / T_val

    expected = V0**2 / 2
    assert simplify(V_sq_avg - expected) == 0


def test_average_power():
    """Verify: P_avg = V_rms I_rms cos(φ) = V₀I₀cos(φ)/2."""
    V0, I0, phi = symbols("V_0 I_0 phi", positive=True)

    V_rms = V0 / sqrt(2)
    I_rms = I0 / sqrt(2)
    P_avg = V_rms * I_rms * cos(phi)
    expected = Rational(1, 2) * V0 * I0 * cos(phi)

    assert simplify(P_avg - expected) == 0


# ── Maxwell & EM waves ─────────────────────────────────────────


def test_speed_of_light():
    """Verify: c = 1/√(μ₀ε₀)."""
    mu0, eps0, c = symbols("mu_0 epsilon_0 c", positive=True)

    c_derived = 1 / sqrt(mu0 * eps0)

    # c² = 1/(μ₀ε₀)
    assert simplify(c_derived**2 - 1 / (mu0 * eps0)) == 0


def test_em_wave_amplitude_relation():
    """Verify: E₀ = cB₀ from Faraday's law applied to EM wave."""
    E0, B0, c, k, omega = symbols("E_0 B_0 c k omega", positive=True)

    # From ∂E/∂x = -∂B/∂t: kE₀ = ωB₀ → E₀ = (ω/k)B₀ = cB₀
    E0_derived = (omega / k) * B0

    # With c = ω/k
    assert simplify(E0_derived.subs(omega, c * k) - c * B0) == 0


def test_poynting_intensity():
    """Verify: I = cε₀E²_rms = E²_rms/(cμ₀)."""
    eps0, mu0, c, E_rms = symbols("epsilon_0 mu_0 c E_rms", positive=True)

    I1 = c * eps0 * E_rms**2
    I2 = E_rms**2 / (c * mu0)

    # These are equal when c² = 1/(μ₀ε₀) → cε₀ = 1/(cμ₀)
    ratio = simplify(I1 / I2)
    # ratio = c²μ₀ε₀ = 1
    assert simplify(ratio.subs(c**2, 1 / (mu0 * eps0))) == 1


# ── Limiting cases ─────────────────────────────────────────────


def test_capacitor_impedance_dc():
    """At DC (ω → 0), capacitor impedance → ∞ (blocks DC)."""
    C, omega = symbols("C omega", positive=True)

    Z_C = 1 / (omega * C)
    assert limit(Z_C, omega, 0, "+") is oo


def test_inductor_impedance_dc():
    """At DC (ω → 0), inductor impedance → 0 (passes DC)."""
    L, omega = symbols("L omega", positive=True)

    Z_L = omega * L
    assert Z_L.subs(omega, 0) == 0


def test_capacitor_impedance_high_freq():
    """At high frequency (ω → ∞), capacitor impedance → 0."""
    C, omega = symbols("C omega", positive=True)

    Z_C = 1 / (omega * C)
    assert limit(Z_C, omega, oo) == 0


def test_inductor_impedance_high_freq():
    """At high frequency (ω → ∞), inductor impedance → ∞."""
    L, omega = symbols("L omega", positive=True)

    Z_L = omega * L
    assert limit(Z_L, omega, oo) is oo


def test_internal_resistance_open_circuit():
    """With no current, terminal voltage equals EMF."""
    emf, r, I_sym = symbols("emf r I", positive=True)

    V = emf - I_sym * r
    assert V.subs(I_sym, 0) == emf


def test_rc_discharging():
    """RC discharging: q → 0 as t → ∞."""
    R, C, V_s, t = symbols("R C V t", positive=True)

    q = C * V_s * exp(-t / (R * C))

    assert limit(q, t, oo) == 0
    assert q.subs(t, 0) == C * V_s
