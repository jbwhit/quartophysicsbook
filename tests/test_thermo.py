"""Symbolic verification of derivations in thermo.qmd."""

from sympy import (
    Rational,
    log,
    simplify,
    solve,
    sqrt,
    symbols,
)


# ── Ideal gas law ───────────────────────────────────────────────


def test_ideal_gas_nk_equals_nr():
    """Verify: PV = NkT = nRT with R = N_A * k."""
    N, k, T, n, R, N_A = symbols("N k T n R N_A", positive=True)

    # N = n * N_A, R = N_A * k
    pv_molecular = N * k * T
    pv_molar = n * R * T

    # Substituting N = n*N_A and R = N_A*k:
    assert simplify(pv_molecular.subs(N, n * N_A) - pv_molar.subs(R, N_A * k)) == 0


# ── Internal energy ─────────────────────────────────────────────


def test_ideal_gas_internal_energy():
    """Verify: U = f/2 * NkT from equipartition."""
    f, N, k, T = symbols("f N k T", positive=True)

    # Each DOF contributes 1/2 kT per molecule, N molecules
    U = f * Rational(1, 2) * N * k * T

    # For monatomic (f=3): U = 3/2 NkT
    assert simplify(U.subs(f, 3) - Rational(3, 2) * N * k * T) == 0

    # For diatomic (f=5): U = 5/2 NkT
    assert simplify(U.subs(f, 5) - Rational(5, 2) * N * k * T) == 0


# ── Kinetic theory ──────────────────────────────────────────────


def test_kinetic_theory_pressure():
    """Verify: P = 1/3 * N*m*v_rms^2 / V and PV = NkT → v_rms = √(3kT/m)."""
    N, m, v_rms, V, k, T = symbols("N m v_rms V k T", positive=True)

    # Pressure from kinetic theory
    P_kinetic = Rational(1, 3) * N * m * v_rms**2 / V

    # Ideal gas: P = NkT/V
    P_ideal = N * k * T / V

    # Equating: 1/3 m v_rms^2 = kT → v_rms = √(3kT/m)
    eq = P_kinetic - P_ideal
    v_solutions = solve(eq, v_rms)
    expected = sqrt(3 * k * T / m)

    assert any(simplify(sol - expected) == 0 for sol in v_solutions)


def test_kinetic_theory_rho_form():
    """Verify: P = 1/3 ρ v_rms^2."""
    N, m, v_rms, V = symbols("N m v_rms V", positive=True)

    rho = N * m / V
    P_kinetic = Rational(1, 3) * N * m * v_rms**2 / V
    P_rho = Rational(1, 3) * rho * v_rms**2

    assert simplify(P_kinetic - P_rho) == 0


# ── Specific heat capacities ────────────────────────────────────


def test_cv():
    """Verify: Cv = f/2 * R."""
    f, R = symbols("f R", positive=True)

    Cv = Rational(1, 2) * f * R

    # Monatomic (f=3): Cv = 3/2 R
    assert simplify(Cv.subs(f, 3) - Rational(3, 2) * R) == 0

    # Diatomic (f=5): Cv = 5/2 R
    assert simplify(Cv.subs(f, 5) - Rational(5, 2) * R) == 0


def test_cp():
    """Verify: Cp = (f/2 + 1) R = Cv + R."""
    f, R = symbols("f R", positive=True)

    Cv = Rational(1, 2) * f * R
    Cp = (Rational(1, 2) * f + 1) * R

    # Cp = Cv + R (Mayer relation)
    assert simplify(Cp - Cv - R) == 0


def test_gamma():
    """Verify: γ = Cp/Cv = (f+2)/f."""
    f, R = symbols("f R", positive=True)

    Cv = Rational(1, 2) * f * R
    Cp = (Rational(1, 2) * f + 1) * R

    gamma = Cp / Cv
    expected = (f + 2) / f

    assert simplify(gamma - expected) == 0

    # Monatomic (f=3): γ = 5/3
    assert simplify(expected.subs(f, 3) - Rational(5, 3)) == 0

    # Diatomic (f=5): γ = 7/5
    assert simplify(expected.subs(f, 5) - Rational(7, 5)) == 0


# ── Isothermal work ─────────────────────────────────────────────


def test_isothermal_work():
    """Verify: W_isotherm = -NkT ln(V2/V1) from integrating P = NkT/V."""
    from sympy import integrate

    N, k, T, V, V1, V2 = symbols("N k T V V1 V2", positive=True)

    P = N * k * T / V
    W = -integrate(P, (V, V1, V2))
    expected = -N * k * T * log(V2 / V1)

    assert simplify(W - expected) == 0


# ── Adiabatic process ───────────────────────────────────────────


def test_adiabatic_pv_gamma():
    """Verify: PV^γ = const → doubling V reduces P by factor 2^γ."""
    P1, V1, gamma = symbols("P1 V1 gamma", positive=True)

    # P1 * V1^gamma = P2 * (2*V1)^gamma
    P2 = P1 * V1**gamma / (2 * V1) ** gamma
    ratio = simplify(P1 / P2)

    assert simplify(ratio - 2**gamma) == 0


# ── Carnot efficiency ───────────────────────────────────────────


def test_carnot_efficiency():
    """Verify: e_Carnot = 1 - T_C/T_H."""
    T_H, T_C = symbols("T_H T_C", positive=True)

    e = 1 - T_C / T_H

    # T_C = 0 → e = 1 (perfect efficiency)
    assert e.subs(T_C, 0) == 1

    # T_C = T_H → e = 0 (no temperature difference → no work)
    assert e.subs(T_C, T_H) == 0

    # e is always < 1 for T_C > 0 (can verify symbolically)
    # 1 - T_C/T_H < 1 ⟺ T_C/T_H > 0, which is true for positive temps


def test_carnot_efficiency_from_heats():
    """Verify: e = 1 - |Q_C|/|Q_H| = 1 - T_C/T_H when |Q_C|/|Q_H| = T_C/T_H."""
    T_H, T_C, Q_H, Q_C = symbols("T_H T_C Q_H Q_C", positive=True)

    # For Carnot: |Q_C|/|Q_H| = T_C/T_H
    e_from_heats = 1 - Q_C / Q_H
    e_carnot = 1 - T_C / T_H

    # Substituting the Carnot relation:
    assert simplify(e_from_heats.subs(Q_C, Q_H * T_C / T_H) - e_carnot) == 0


# ── Carnot refrigerator COP ────────────────────────────────────


def test_carnot_cop():
    """Verify: K_Carnot = T_C / (T_H - T_C)."""
    T_H, T_C = symbols("T_H T_C", positive=True)

    K = T_C / (T_H - T_C)

    # As T_C → 0: K → 0 (can't pump heat from absolute zero)
    assert K.subs(T_C, 0) == 0

    # Verify K = |Q_C|/|W| = |Q_C|/(|Q_H| - |Q_C|) with Carnot relation
    Q_H, Q_C = symbols("Q_H Q_C", positive=True)
    K_from_heats = Q_C / (Q_H - Q_C)
    K_carnot = T_C / (T_H - T_C)

    assert simplify(K_from_heats.subs(Q_C, Q_H * T_C / T_H) - K_carnot) == 0


# ── Limiting cases ──────────────────────────────────────────────


def test_vrms_limiting_cases():
    """Limiting cases for v_rms = √(3kT/m)."""
    from sympy import limit, oo

    k, T, m = symbols("k T m", positive=True)

    v_rms = sqrt(3 * k * T / m)

    # As T → 0: v_rms → 0
    assert v_rms.subs(T, 0) == 0

    # As m → ∞: v_rms → 0 (heavy molecules move slowly)
    assert limit(v_rms, m, oo) == 0

    # As T → ∞: v_rms → ∞
    assert limit(v_rms, T, oo) is oo


def test_carnot_engine_plus_refrigerator():
    """Verify: Carnot engine driving Carnot refrigerator between same temps does no net work."""
    T_H, T_C, Q_H = symbols("T_H T_C Q_H", positive=True)

    e = 1 - T_C / T_H
    W = e * Q_H  # Work produced by engine

    # Refrigerator COP
    K = T_C / (T_H - T_C)
    Q_C_ref = K * W  # Heat pumped by refrigerator

    # Heat rejected by engine to cold reservoir
    Q_C_engine = Q_H * T_C / T_H

    # Refrigerator pumps same heat back: Q_C_ref should equal Q_C_engine
    assert simplify(Q_C_ref - Q_C_engine) == 0
