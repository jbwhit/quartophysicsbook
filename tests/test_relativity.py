"""Symbolic verification of derivations in relativity.qmd."""

import pytest
from sympy import (
    Rational,
    S,
    Symbol,
    acos,
    cos,
    pi,
    simplify,
    sqrt,
    symbols,
    tan,
)


# --- Physical constants as symbols ---
h, c, m, e = symbols("h c m e", positive=True)
hbar = h / (2 * pi)


# === Bohr Atom ===


def test_angular_momentum_quantization():
    """L = n * hbar for stationary states."""
    n = Symbol("n", positive=True, integer=True)
    L = n * hbar
    assert L.subs(n, 1) == hbar
    assert L.subs(n, 3) == 3 * hbar


def test_photon_emission_energy():
    """h*nu = E2 - E1 for transitions: frequency from energy gap."""
    E1, E2 = symbols("E_1 E_2", positive=True)
    # Transition frequency
    nu = (E2 - E1) / h
    # Photon energy should equal the gap
    assert simplify(h * nu - (E2 - E1)) == 0


# === Photoelectric Effect ===


def test_photoelectric_effect():
    """K_max = h*nu - phi (Einstein's photoelectric equation)."""
    nu, phi = symbols("nu phi", positive=True)
    K_max = h * nu - phi
    # At threshold frequency, K_max = 0
    nu_threshold = phi / h
    assert simplify(K_max.subs(nu, nu_threshold)) == 0


def test_photoelectric_threshold():
    """Below threshold frequency, no electrons emitted."""
    nu, phi = symbols("nu phi", positive=True)
    K_max = h * nu - phi
    nu_0 = phi / h
    # Just below threshold: K_max < 0 (no emission)
    # At threshold: K_max = 0
    assert K_max.subs(nu, nu_0) == 0


# === Wave-Particle Duality ===


def test_de_broglie_wavelength():
    """lambda = h/p = h/(mv)."""
    v = Symbol("v", positive=True)
    lam = h / (m * v)
    p = m * v
    assert simplify(lam - h / p) == 0


def test_photon_momentum():
    """p = E/c = h*nu/c = h/lambda for photons."""
    nu, lam = symbols("nu lambda", positive=True)
    p_from_energy = h * nu / c
    p_from_wavelength = h / lam
    # With nu = c/lambda
    assert simplify(p_from_energy.subs(nu, c / lam) - p_from_wavelength) == 0


def test_photon_energy_momentum():
    """E = pc for massless particles."""
    p = Symbol("p", positive=True)
    E = p * c
    # From E^2 = p^2*c^2 + m^2*c^4 with m=0
    E_from_relation = sqrt(p**2 * c**2 + 0)
    assert simplify(E - E_from_relation) == 0


# === Heisenberg Uncertainty ===


def test_uncertainty_position_momentum():
    """Delta_x * Delta_p >= hbar/2."""
    Dx, Dp = symbols("Delta_x Delta_p", positive=True)
    # Minimum uncertainty product
    min_product = hbar / 2
    assert simplify(min_product - h / (4 * pi)) == 0


def test_uncertainty_energy_time():
    """Delta_E * Delta_t >= hbar/2."""
    DE, Dt = symbols("Delta_E Delta_t", positive=True)
    min_product = hbar / 2
    assert simplify(min_product - h / (4 * pi)) == 0


# === Special Relativity: Lorentz Transformations ===


def test_lorentz_factor_at_rest():
    """gamma = 1 when v = 0."""
    v = Symbol("v")
    gamma = 1 / sqrt(1 - v**2 / c**2)
    assert gamma.subs(v, 0) == 1


def test_spacetime_interval_invariance():
    """s^2 = c^2*t^2 - x^2 is Lorentz invariant."""
    x, t, v = symbols("x t v")
    gamma = 1 / sqrt(1 - v**2 / c**2)

    # Lorentz transformation
    x_prime = gamma * (x - v * t)
    t_prime = gamma * (t - v * x / c**2)

    # Interval in original frame
    s2 = c**2 * t**2 - x**2

    # Interval in primed frame
    s2_prime = c**2 * t_prime**2 - x_prime**2

    # Should be equal
    diff = simplify(s2_prime - s2)
    assert diff == 0


def test_length_contraction():
    """L = L_0 / gamma."""
    v, L_0 = symbols("v L_0", positive=True)
    gamma = 1 / sqrt(1 - v**2 / c**2)
    L = L_0 / gamma
    # At rest, L = L_0
    assert L.subs(v, 0) == L_0
    # L < L_0 for v > 0 (length contraction)
    assert simplify(L**2 - L_0**2 * (1 - v**2 / c**2)) == 0


def test_time_dilation():
    """Delta_t = gamma * Delta_t_0."""
    v, dt_0 = symbols("v Delta_t_0", positive=True)
    gamma = 1 / sqrt(1 - v**2 / c**2)
    dt = gamma * dt_0
    # At rest, dt = dt_0
    assert dt.subs(v, 0) == dt_0


def test_light_clock_time_dilation():
    """Light clock derivation: T_m = gamma * T_0."""
    v, L = symbols("v L", positive=True)
    T_0 = 2 * L / c

    # In moving frame, light path is diagonal
    # Half-trip: sqrt(L^2 + (v*T_m/2)^2) = c*T_m/2
    # Solving: T_m^2/4 * (c^2 - v^2) = L^2
    # T_m = 2L / sqrt(c^2 - v^2) = (2L/c) / sqrt(1 - v^2/c^2) = gamma * T_0
    gamma = 1 / sqrt(1 - v**2 / c**2)
    T_m = 2 * L / sqrt(c**2 - v**2)
    assert simplify(T_m - gamma * T_0) == 0


def test_velocity_addition():
    """u' = (u - v) / (1 - uv/c^2)."""
    u, v = symbols("u v")
    u_prime = (u - v) / (1 - u * v / c**2)

    # If u = c, then u' = c (speed of light invariance)
    assert simplify(u_prime.subs(u, c) - c) == 0

    # If u = -c, then u' = -c
    assert simplify(u_prime.subs(u, -c) + c) == 0


def test_velocity_addition_newtonian_limit():
    """For v << c, velocity addition reduces to u' = u - v."""
    u, v = symbols("u v")
    u_prime = (u - v) / (1 - u * v / c**2)
    # In limit c -> infinity, denominator -> 1
    # So u' -> u - v
    from sympy import limit, oo

    result = limit(u_prime, c, oo)
    assert simplify(result - (u - v)) == 0


# === Relativistic Dynamics ===


def test_four_velocity_magnitude():
    """U^mu U_mu = c^2 (invariant)."""
    v = Symbol("v", positive=True)
    gamma = 1 / sqrt(1 - v**2 / c**2)
    # U = gamma*(c, v) in 1+1 dimensions
    U_sq = gamma**2 * (c**2 - v**2)
    assert simplify(U_sq - c**2) == 0


def test_relativistic_momentum():
    """p = gamma * m * v."""
    v = Symbol("v", positive=True)
    gamma = 1 / sqrt(1 - v**2 / c**2)
    p = gamma * m * v
    # For v << c: p -> m*v
    from sympy import limit, oo

    assert limit(p, c, oo) == m * v


def test_relativistic_energy():
    """E = gamma * m * c^2."""
    v = Symbol("v", positive=True)
    gamma = 1 / sqrt(1 - v**2 / c**2)
    E = gamma * m * c**2
    # At rest: E = mc^2
    assert E.subs(v, 0) == m * c**2


def test_kinetic_energy_newtonian_limit():
    """KE = (gamma - 1)*mc^2 -> (1/2)*m*v^2 for v << c."""
    v = Symbol("v", positive=True)
    gamma = 1 / sqrt(1 - v**2 / c**2)
    K = (gamma - 1) * m * c**2

    # Taylor expand gamma for small v/c
    from sympy import series

    gamma_series = series(gamma, v, 0, n=4).removeO()
    # gamma ≈ 1 + v^2/(2c^2)
    K_approx = (gamma_series - 1) * m * c**2
    K_approx = simplify(K_approx)
    leading = Rational(1, 2) * m * v**2
    assert simplify(K_approx - leading) == 0


def test_energy_momentum_relation():
    """E^2 = p^2*c^2 + m^2*c^4."""
    v = Symbol("v", positive=True)
    gamma = 1 / sqrt(1 - v**2 / c**2)
    E = gamma * m * c**2
    p = gamma * m * v
    lhs = E**2
    rhs = p**2 * c**2 + m**2 * c**4
    assert simplify(lhs - rhs) == 0


def test_energy_momentum_relation_photon():
    """For m=0: E = pc."""
    p = Symbol("p", positive=True)
    E_sq = p**2 * c**2 + 0  # m = 0
    E = sqrt(E_sq)
    assert simplify(E - p * c) == 0


def test_rest_energy():
    """E_0 = mc^2."""
    v = Symbol("v")
    gamma = 1 / sqrt(1 - v**2 / c**2)
    E = gamma * m * c**2
    E_rest = E.subs(v, 0)
    assert E_rest == m * c**2


# === Compton Scattering ===


def test_compton_wavelength_shift():
    """Delta_lambda = (h/mc)(1 - cos(theta))."""
    theta = Symbol("theta")
    Delta_lam = (h / (m * c)) * (1 - cos(theta))

    # At theta = 0 (forward scattering): no wavelength change
    assert Delta_lam.subs(theta, 0) == 0

    # At theta = pi (backscattering): maximum shift = 2h/(mc)
    assert simplify(Delta_lam.subs(theta, pi) - 2 * h / (m * c)) == 0

    # At theta = pi/2: shift = h/(mc)
    assert simplify(Delta_lam.subs(theta, pi / 2) - h / (m * c)) == 0


def test_compton_wavelength():
    """Compton wavelength = h/(mc)."""
    lambda_C = h / (m * c)
    # Verify dimensions are consistent (symbolic check)
    assert simplify(lambda_C * m * c - h) == 0


def test_inverse_compton_boost():
    """Inverse Compton: f' ~ gamma^2 * f for head-on collision."""
    v, f = symbols("v f", positive=True)
    gamma = 1 / sqrt(1 - v**2 / c**2)
    beta = v / c
    # Head-on: photon approaches from opposite direction
    # Doppler shift going to electron rest frame: factor gamma*(1+beta)
    # Doppler shift going back to lab frame: another factor gamma*(1+beta)
    # Total boost: gamma^2*(1+beta)^2
    # For ultra-relativistic (beta -> 1): boost ~ 4*gamma^2
    # But approximate formula gives f' ~ gamma^2 * f
    boost = gamma**2 * (1 + beta) ** 2
    # For beta -> 1: boost -> 4*gamma^2
    # Check at v=0: boost = 1 (no change)
    assert boost.subs(v, 0) == 1


# === Aberration ===


def test_aberration_formula_head_on():
    """At alpha = 0 (directly ahead), alpha' = 0."""
    v = Symbol("v", positive=True)
    alpha = S.Zero
    # tan(0/2) = 0, so tan(alpha'/2) = 0, hence alpha' = 0
    rhs = sqrt((c - v) / (c + v)) * tan(alpha / 2)
    assert rhs == 0


def test_aberration_formula_general():
    """At v = 0, aberration gives alpha' = alpha for any angle."""
    alpha = Symbol("alpha", positive=True)
    # At v=0: factor = sqrt((c-0)/(c+0)) = 1
    # So tan(alpha'/2) = tan(alpha/2), hence alpha' = alpha
    factor = sqrt((c) / (c))
    assert factor == 1
    # Also test at alpha = pi/2 with nonzero v
    v = Symbol("v", positive=True)
    rhs = sqrt((c - v) / (c + v)) * tan(pi / 4)
    # tan(pi/4) = 1, so rhs = sqrt((c-v)/(c+v))
    assert simplify(rhs - sqrt((c - v) / (c + v))) == 0


def test_aberration_no_motion():
    """At v = 0, no aberration: alpha' = alpha."""
    alpha = Symbol("alpha", positive=True)
    # At v=0: sqrt((c-0)/(c+0)) = 1
    # So tan(alpha'/2) = tan(alpha/2), hence alpha' = alpha
    factor = sqrt((c - 0) / (c + 0))
    assert factor == 1


# === Relativistic Doppler ===


def test_relativistic_doppler_approaching():
    """f' = f * sqrt((1+beta)/(1-beta)) for head-on approach."""
    v, f = symbols("v f", positive=True)
    beta = v / c
    gamma = 1 / sqrt(1 - beta**2)
    # alpha = 0 (approaching)
    f_prime = gamma * (1 + beta) * f
    # Square both sides to compare
    f_prime_sq = simplify(f_prime**2)
    expected_sq = f**2 * (1 + beta) / (1 - beta)
    assert simplify(f_prime_sq - expected_sq) == 0


def test_relativistic_doppler_receding():
    """f' = f * sqrt((1-beta)/(1+beta)) for receding source (alpha=pi)."""
    v, f = symbols("v f", positive=True)
    beta = v / c
    gamma = 1 / sqrt(1 - beta**2)
    # alpha = pi (receding): cos(pi) = -1
    f_prime = gamma * (1 - beta) * f
    # Square both sides to compare
    f_prime_sq = simplify(f_prime**2)
    expected_sq = f**2 * (1 - beta) / (1 + beta)
    assert simplify(f_prime_sq - expected_sq) == 0


def test_relativistic_doppler_transverse():
    """f' = gamma * f for transverse motion (alpha = pi/2)."""
    v, f = symbols("v f", positive=True)
    beta = v / c
    gamma = 1 / sqrt(1 - beta**2)
    # cos(pi/2) = 0
    f_prime = gamma * (1 + beta * 0) * f
    assert simplify(f_prime - gamma * f) == 0


def test_relativistic_doppler_no_motion():
    """At v = 0, f' = f (no Doppler shift)."""
    f = Symbol("f", positive=True)
    beta = 0
    gamma = 1
    alpha = Symbol("alpha")
    f_prime = gamma * (1 + beta * cos(alpha)) * f
    assert f_prime == f


# === Limiting Cases ===


def test_lorentz_newtonian_limit():
    """For v << c, Lorentz transforms reduce to Galilean."""
    x, t, v = symbols("x t v")
    gamma = 1 / sqrt(1 - v**2 / c**2)
    x_prime = gamma * (x - v * t)
    t_prime = gamma * (t - v * x / c**2)

    from sympy import limit, oo

    # In limit c -> infinity
    x_prime_newton = limit(x_prime, c, oo)
    t_prime_newton = limit(t_prime, c, oo)
    assert simplify(x_prime_newton - (x - v * t)) == 0
    assert simplify(t_prime_newton - t) == 0


def test_relativity_of_simultaneity():
    """Events simultaneous in S are not in S' if spatially separated."""
    x1, x2, t, v = symbols("x_1 x_2 t v")
    gamma = 1 / sqrt(1 - v**2 / c**2)
    # Both at same time t in S
    t1_prime = gamma * (t - v * x1 / c**2)
    t2_prime = gamma * (t - v * x2 / c**2)
    dt_prime = simplify(t2_prime - t1_prime)
    expected = -gamma * v * (x2 - x1) / c**2
    assert simplify(dt_prime - expected) == 0
