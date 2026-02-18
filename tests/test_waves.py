"""Symbolic verification of derivations in waves.qmd."""

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


# ── SHM angular frequency ─────────────────────────────────────


def test_shm_angular_frequency():
    """Verify: ω = √(k/m) from m*x'' + k*x = 0."""
    k, m, omega = symbols("k m omega", positive=True)

    # Characteristic equation: m*omega^2 = k
    solutions = solve(m * omega**2 - k, omega)
    expected = sqrt(k / m)

    assert any(simplify(sol - expected) == 0 for sol in solutions)


# ── Spring PE ──────────────────────────────────────────────────


def test_spring_potential_energy():
    """Verify: U = ½kx² from integrating F = -kx."""
    from sympy import integrate

    k, x = symbols("k x", positive=True)

    F_magnitude = k * x
    U = integrate(F_magnitude, (x, 0, x))
    expected = Rational(1, 2) * k * x**2

    assert simplify(U - expected) == 0


# ── Energy conservation in SHM ─────────────────────────────────


def test_shm_energy_conservation():
    """Verify: E_T = ½kx₀² is constant (using m*ω² = k)."""
    k, m, x0, omega, t, phi = symbols("k m x_0 omega t phi", real=True)

    x = x0 * cos(omega * t + phi)
    v = -x0 * omega * sin(omega * t + phi)

    U = Rational(1, 2) * k * x**2
    K = Rational(1, 2) * m * v**2
    E_total = U + K

    # Substitute k = m*omega^2
    E_total_sub = E_total.subs(k, m * omega**2)
    E_total_simplified = simplify(E_total_sub)

    expected = Rational(1, 2) * m * omega**2 * x0**2

    assert simplify(E_total_simplified - expected) == 0


# ── Simple pendulum ────────────────────────────────────────────


def test_simple_pendulum_period():
    """Verify: T = 2π√(L/g) from k_eff = mg/L."""
    m, g, L = symbols("m g L", positive=True)

    k_eff = m * g / L
    omega = sqrt(k_eff / m)
    T = 2 * pi / omega
    expected = 2 * pi * sqrt(L / g)

    assert simplify(T - expected) == 0


# ── Wave speed ─────────────────────────────────────────────────


def test_wave_speed():
    """Verify: v = ω/k = λf."""
    omega, k, lam, f = symbols("omega k lambda f", positive=True)

    # v = ω/k, and ω = 2πf, k = 2π/λ
    v_from_omega_k = omega / k
    v_from_lambda_f = lam * f

    v_sub = v_from_omega_k.subs([(omega, 2 * pi * f), (k, 2 * pi / lam)])

    assert simplify(v_sub - v_from_lambda_f) == 0


# ── String wave speed ─────────────────────────────────────────


def test_string_wave_speed():
    """Verify: v = √(F_T/μ) — doubling tension increases speed by √2."""
    F_T, mu = symbols("F_T mu", positive=True)

    v = sqrt(F_T / mu)
    v_doubled = sqrt(2 * F_T / mu)

    ratio = simplify(v_doubled / v)
    assert simplify(ratio - sqrt(2)) == 0


# ── Standing wave from superposition ──────────────────────────


def test_standing_wave_superposition():
    """Verify: y₁ + y₂ = 2y₀sin(kx)cos(ωt)."""
    y0, k, x, omega, t = symbols("y_0 k x omega t", real=True)

    y1 = y0 * sin(k * x - omega * t)
    y2 = y0 * sin(k * x + omega * t)

    total = simplify(y1 + y2)
    expected = 2 * y0 * sin(k * x) * cos(omega * t)

    assert simplify(total - expected) == 0


# ── Resonance: two fixed ends ─────────────────────────────────


def test_resonance_two_fixed_ends():
    """Verify: f_n = nv/(2L) from L = nλ/2."""
    n, v, L = symbols("n v L", positive=True)

    # L = n*λ/2 → λ = 2L/n
    lam = 2 * L / n
    f_n = v / lam
    expected = n * v / (2 * L)

    assert simplify(f_n - expected) == 0


# ── Resonance: one closed end ─────────────────────────────────


def test_resonance_one_closed_end():
    """Verify: f_n = (2n-1)v/(4L) from L = (2n-1)λ/4."""
    n, v, L = symbols("n v L", positive=True)

    lam = 4 * L / (2 * n - 1)
    f_n = v / lam
    expected = (2 * n - 1) * v / (4 * L)

    assert simplify(f_n - expected) == 0


# ── Beats ─────────────────────────────────────────────────────


def test_beat_superposition():
    """Verify: cos(ω₁t) + cos(ω₂t) = 2cos((ω₁-ω₂)/2 t)cos((ω₁+ω₂)/2 t)."""
    omega1, omega2, t = symbols("omega_1 omega_2 t", real=True)

    lhs = cos(omega1 * t) + cos(omega2 * t)
    rhs = 2 * cos((omega1 - omega2) / 2 * t) * cos((omega1 + omega2) / 2 * t)

    assert simplify(lhs - rhs) == 0


# ── Doppler effect ────────────────────────────────────────────


def test_doppler_moving_source():
    """Verify: f' = f*V/(V - V_s) for moving source, stationary observer."""
    f, V, V_s = symbols("f V V_s", positive=True)

    # Wavelength compressed: λ' = λ(1 - V_s/V) = (V - V_s)/f
    lam = V / f
    lam_prime = lam * (1 - V_s / V)
    f_prime = V / lam_prime
    expected = f * V / (V - V_s)

    assert simplify(f_prime - expected) == 0


def test_doppler_moving_observer():
    """Verify: f' = f*(V + V_o)/V for moving observer, stationary source."""
    f, V, V_o = symbols("f V V_o", positive=True)

    lam = V / f
    # Observer approaching: encounters crests at rate (V + V_o)/λ
    f_prime = (V + V_o) / lam
    expected = f * (V + V_o) / V

    assert simplify(f_prime - expected) == 0


def test_doppler_general_formula():
    """Verify: general Doppler reduces to special cases."""
    f, V, V_s, V_o = symbols("f V V_s V_o", positive=True)

    f_general = f * (V - V_o) / (V - V_s)

    # Case 1: stationary observer (V_o = 0)
    assert simplify(f_general.subs(V_o, 0) - f * V / (V - V_s)) == 0

    # Case 2: stationary source (V_s = 0)
    assert simplify(f_general.subs(V_s, 0) - f * (V - V_o) / V) == 0


# ── Limiting cases ────────────────────────────────────────────


def test_pendulum_limiting_cases():
    """Limiting cases for T = 2π√(L/g)."""
    from sympy import limit, oo

    L, g = symbols("L g", positive=True)

    T = 2 * pi * sqrt(L / g)

    # As L → 0: T → 0 (short pendulum oscillates fast)
    assert T.subs(L, 0) == 0

    # As g → ∞: T → 0 (strong gravity → fast oscillation)
    assert limit(T, g, oo) == 0

    # As L → ∞: T → ∞ (long pendulum is slow)
    assert limit(T, L, oo) is oo


def test_doppler_sonic_boom():
    """As V_s → V, f' → ∞ (sonic boom)."""
    from sympy import limit, oo

    f, V, V_s = symbols("f V V_s", positive=True)

    f_prime = f * V / (V - V_s)

    assert limit(f_prime, V_s, V, "+") is oo or limit(f_prime, V_s, V, "-") is oo


def test_shm_amplitude_zero_energy():
    """If amplitude is zero, total energy is zero."""
    k, x0 = symbols("k x_0", positive=True)

    E = Rational(1, 2) * k * x0**2

    assert E.subs(x0, 0) == 0


def test_wave_power_proportional_to_amplitude_squared():
    """Verify P ∝ y₀² and P ∝ ω²."""
    mu, omega, y0, v = symbols("mu omega y_0 v", positive=True)

    P = Rational(1, 2) * mu * omega**2 * y0**2 * v

    # Doubling amplitude → 4× power
    P_doubled_amp = P.subs(y0, 2 * y0)
    assert simplify(P_doubled_amp / P) == 4

    # Doubling frequency → 4× power
    P_doubled_freq = P.subs(omega, 2 * omega)
    assert simplify(P_doubled_freq / P) == 4


# ── Guitar string fundamental frequency ─────────────────────


def test_guitar_string_fundamental():
    """Verify: f₁ = (1/2L)√(F_T/μ) from wave speed and resonance."""
    L, F_T, mu = symbols("L F_T mu", positive=True)

    v = sqrt(F_T / mu)
    f1 = v / (2 * L)
    expected = 1 / (2 * L) * sqrt(F_T / mu)

    assert simplify(f1 - expected) == 0


def test_guitar_string_limiting_cases():
    """Limiting cases for guitar string frequency."""
    from sympy import limit, oo

    L, F_T, mu = symbols("L F_T mu", positive=True)

    f1 = sqrt(F_T / mu) / (2 * L)

    # Doubling tension → frequency increases by √2
    f1_double_tension = f1.subs(F_T, 2 * F_T)
    assert simplify(f1_double_tension / f1 - sqrt(2)) == 0

    # Doubling length → frequency halves
    f1_double_length = f1.subs(L, 2 * L)
    assert simplify(f1_double_length / f1 - Rational(1, 2)) == 0

    # As tension → 0: f → 0 (no restoring force, no vibration)
    assert f1.subs(F_T, 0) == 0


# ── Beat frequency and period ────────────────────────────────


def test_beat_frequency():
    """Verify: f_beat = |f1 - f2|."""
    from sympy import Abs

    f1, f2 = symbols("f1 f2", positive=True)

    f_beat = Abs(f1 - f2)
    T_beat = 1 / f_beat

    # Numerical check: 440 and 444 Hz
    assert f_beat.subs([(f1, 440), (f2, 444)]) == 4
    assert T_beat.subs([(f1, 440), (f2, 444)]) == Rational(1, 4)


def test_beat_frequency_equal_frequencies():
    """If f1 = f2, no beats (f_beat = 0)."""
    f = symbols("f", positive=True)

    f_beat = f - f
    assert f_beat == 0
