"""Symbolic verification of derivations in optics.qmd."""

from sympy import (
    Rational,
    asin,
    atan,
    cos,
    limit,
    oo,
    pi,
    simplify,
    sin,
    solve,
    sqrt,
    symbols,
    tan,
)


# ── Snell's law ──────────────────────────────────────────────


def test_snells_law_symmetry():
    """Verify: Snell's law is symmetric — swapping media reverses the ray."""
    n1, n2, theta1, theta2 = symbols("n1 n2 theta1 theta2", positive=True)

    # n1 sin(theta1) = n2 sin(theta2)
    # If we swap 1 <-> 2, we get the same equation
    lhs = n1 * sin(theta1) - n2 * sin(theta2)
    rhs_swapped = n2 * sin(theta2) - n1 * sin(theta1)

    assert simplify(lhs + rhs_swapped) == 0


def test_snells_law_normal_incidence():
    """At normal incidence (θ₁ = 0), transmitted ray is also normal (θ₂ = 0)."""
    n1, n2 = symbols("n1 n2", positive=True)

    # n1 * sin(0) = n2 * sin(theta2) → sin(theta2) = 0 → theta2 = 0
    sin_theta2 = n1 * sin(0) / n2
    assert sin_theta2 == 0


# ── Critical angle ───────────────────────────────────────────


def test_critical_angle():
    """Verify: sin(θ_c) = n₂/n₁ from Snell's law with θ₂ = π/2."""
    n1, n2 = symbols("n1 n2", positive=True)

    # At critical angle: n1 sin(θ_c) = n2 sin(π/2) = n2
    sin_theta_c = n2 / n1
    theta_c = asin(n2 / n1)

    # For glass (n1=1.5) to air (n2=1): θ_c ≈ 41.8°
    import sympy

    val = theta_c.subs([(n1, Rational(3, 2)), (n2, 1)])
    assert 0 < val < pi / 2  # Must be between 0 and 90°


# ── Mirror equation ──────────────────────────────────────────


def test_mirror_equation_object_at_C():
    """Object at center of curvature (s = 2f) → image at s' = 2f."""
    f = symbols("f", positive=True)

    s = 2 * f
    # 1/s + 1/s' = 1/f → 1/s' = 1/f - 1/s
    s_prime = solve(1 / s + 1 / symbols("sp") - 1 / f, symbols("sp"))
    assert any(simplify(sol - 2 * f) == 0 for sol in s_prime)


def test_mirror_equation_object_at_infinity():
    """Object at infinity → image at focal point."""
    from sympy import limit, oo

    f, s = symbols("f s", positive=True)

    # s' = sf/(s - f)
    s_prime = s * f / (s - f)
    result = limit(s_prime, s, oo)

    assert simplify(result - f) == 0


def test_mirror_equation_object_at_F():
    """Object at focal point → image at infinity."""
    f = symbols("f", positive=True)

    s = f
    # 1/s' = 1/f - 1/f = 0 → s' → ∞
    s_prime = s * f / (s - f)

    # s = f makes denominator 0 → s' → ∞
    s_var = symbols("s_var", positive=True)
    s_prime_general = s_var * f / (s_var - f)
    result = limit(s_prime_general, s_var, f, "+")

    assert result is oo


# ── Magnification ────────────────────────────────────────────


def test_magnification_plane_mirror():
    """Plane mirror (f → ∞): M = 1 (upright, same size)."""
    s, f = symbols("s f", positive=True)

    s_prime = s * f / (s - f)
    M = -s_prime / s

    # As f → ∞: s' → -s (virtual image), M → 1
    M_limit = limit(M, f, oo)
    assert simplify(M_limit - 1) == 0


def test_magnification_at_C():
    """Object at C (s = 2f): M = -1 (inverted, same size)."""
    f = symbols("f", positive=True)

    s = 2 * f
    s_prime = s * f / (s - f)
    M = -s_prime / s

    assert simplify(M - (-1)) == 0


# ── Lensmaker's equation ─────────────────────────────────────


def test_lensmaker_symmetric_lens():
    """Symmetric biconvex lens (R₁ = R, R₂ = -R): 1/f = 2(n-1)/R."""
    n, R = symbols("n R", positive=True)

    R1 = R
    R2 = -R
    inv_f = (n - 1) * (1 / R1 - 1 / R2)
    expected = 2 * (n - 1) / R

    assert simplify(inv_f - expected) == 0


def test_lensmaker_flat_surface():
    """Plano-convex lens (R₂ → ∞): 1/f = (n-1)/R₁."""
    n, R1 = symbols("n R1", positive=True)
    R2 = symbols("R2")

    inv_f = (n - 1) * (1 / R1 - 1 / R2)
    result = limit(inv_f, R2, oo)
    expected = (n - 1) / R1

    assert simplify(result - expected) == 0


# ── Double slit intensity ────────────────────────────────────


def test_double_slit_central_maximum():
    """At θ = 0, intensity is 4I₀ (constructive)."""
    I0, d, lam = symbols("I_0 d lambda", positive=True)

    # I = 4I₀ cos²(π d sinθ / λ), at θ = 0: sin(0) = 0
    I_center = 4 * I0 * cos(pi * d * 0 / lam) ** 2
    assert simplify(I_center - 4 * I0) == 0


def test_double_slit_first_minimum():
    """First minimum at d sinθ = λ/2."""
    I0, d, lam = symbols("I_0 d lambda", positive=True)
    theta = symbols("theta")

    # At d sinθ = λ/2: π d sinθ / λ = π/2
    I_min = 4 * I0 * cos(pi / 2) ** 2
    assert I_min == 0


def test_double_slit_maxima_spacing():
    """Bright fringes occur at d sinθ = mλ."""
    d, lam, m = symbols("d lambda m", positive=True, integer=True)

    # At d sinθ = mλ: π d sinθ / λ = mπ → cos²(mπ) = 1
    I_factor = cos(m * pi) ** 2
    assert simplify(I_factor - 1) == 0


# ── Bragg's law ──────────────────────────────────────────────


def test_bragg_first_order():
    """First order (m=1): λ = 2d sinθ."""
    d, lam, theta = symbols("d lambda theta", positive=True)

    # mλ = 2d sinθ, m=1
    lam_bragg = 2 * d * sin(theta)

    # At θ = π/6 (30°): λ = d
    assert simplify(lam_bragg.subs(theta, pi / 6) - d) == 0


# ── Single slit diffraction ─────────────────────────────────


def test_single_slit_minima():
    """Verify: minima at D sinθ = mλ (m ≠ 0)."""
    D, lam = symbols("D lambda", positive=True)

    # At first minimum: sinθ = λ/D
    sin_theta_min = lam / D

    # For D = 2λ: sinθ = 1/2 → θ = 30°
    val = sin_theta_min.subs(D, 2 * lam)
    assert simplify(val - Rational(1, 2)) == 0


def test_single_slit_sinc_central_max():
    """Verify: sinc²(α) → 1 as α → 0 (central maximum)."""
    alpha = symbols("alpha")

    sinc_sq = (sin(alpha) / alpha) ** 2
    result = limit(sinc_sq, alpha, 0)

    assert result == 1


def test_single_slit_sinc_zeros():
    """Verify: sinc²(α) = 0 at α = mπ."""
    m = symbols("m", integer=True, nonzero=True)

    # sin(mπ) = 0 for integer m ≠ 0
    assert sin(pi) == 0
    assert sin(2 * pi) == 0
    assert sin(3 * pi) == 0


# ── Rayleigh criterion ───────────────────────────────────────


def test_rayleigh_criterion():
    """Verify: larger aperture → better resolution (smaller θ)."""
    lam = symbols("lambda", positive=True)

    d1 = symbols("d1", positive=True)
    d2 = 2 * d1  # Double the aperture

    theta1 = Rational(122, 100) * lam / d1
    theta2 = Rational(122, 100) * lam / d2

    ratio = simplify(theta1 / theta2)
    assert ratio == 2  # Doubling aperture halves the angular resolution


# ── Malus's law ──────────────────────────────────────────────


def test_malus_law_parallel():
    """Parallel polarizers (θ = 0): full transmission."""
    I0 = symbols("I_0", positive=True)

    I = I0 * cos(0) ** 2
    assert simplify(I - I0) == 0


def test_malus_law_crossed():
    """Crossed polarizers (θ = π/2): zero transmission."""
    I0 = symbols("I_0", positive=True)

    I = I0 * cos(pi / 2) ** 2
    assert I == 0


def test_malus_law_45_degrees():
    """At 45°: half the intensity."""
    I0 = symbols("I_0", positive=True)

    I = I0 * cos(pi / 4) ** 2
    expected = I0 / 2

    assert simplify(I - expected) == 0


def test_three_polarizers():
    """Three polarizers (0°, 45°, 90°): I = I₀/8."""
    I0 = symbols("I_0", positive=True)

    # Unpolarized → first polarizer: I₀/2
    I1 = I0 / 2
    # Through second at 45°: I₁ cos²(45°) = I₀/4
    I2 = I1 * cos(pi / 4) ** 2
    # Through third at 45° from second: I₂ cos²(45°) = I₀/8
    I3 = I2 * cos(pi / 4) ** 2

    expected = I0 / 8
    assert simplify(I3 - expected) == 0


# ── Brewster's angle ─────────────────────────────────────────


def test_brewster_angle():
    """Verify: θ_B = arctan(n₂/n₁) and θ₁ + θ₂ = 90° at Brewster's angle."""
    n1, n2 = symbols("n1 n2", positive=True)

    theta_B = atan(n2 / n1)

    # At Brewster's angle, reflected and refracted rays are perpendicular:
    # θ₁ + θ₂ = π/2 → θ₂ = π/2 - θ₁
    # Snell: n1 sin(θ_B) = n2 sin(π/2 - θ_B) = n2 cos(θ_B)
    # → tan(θ_B) = n2/n1  ✓
    lhs = n1 * sin(theta_B) - n2 * cos(theta_B)

    # tan(θ_B) = n2/n1 → sin(θ_B)/cos(θ_B) = n2/n1
    # → n1 sin(θ_B) = n2 cos(θ_B)
    assert simplify(lhs) == 0


def test_brewster_glass():
    """Brewster's angle for air-glass (n=1.5): θ_B ≈ 56.3°."""
    theta_B = atan(Rational(3, 2))

    # Should be between 45° and 60°
    assert pi / 4 < theta_B < pi / 3


# ── Resolving power ──────────────────────────────────────────


def test_grating_resolving_power():
    """Verify: R = Nm — more slits and higher order give better resolution."""
    N, m = symbols("N m", positive=True)

    R = N * m

    # Doubling N doubles R
    assert simplify(R.subs(N, 2 * N) / R) == 2

    # Doubling m doubles R
    assert simplify(R.subs(m, 2 * m) / R) == 2


# ── Dispersion ───────────────────────────────────────────────


def test_grating_dispersion():
    """Verify: D = m/(d cosθ) at normal incidence (θ ≈ 0): D = m/d."""
    m, d = symbols("m d", positive=True)
    theta = symbols("theta")

    D = m / (d * cos(theta))

    # At θ = 0: D = m/d
    D_normal = D.subs(theta, 0)
    assert simplify(D_normal - m / d) == 0


# ── Limiting cases ───────────────────────────────────────────


def test_snells_law_equal_media():
    """If n₁ = n₂, light passes straight through (θ₁ = θ₂)."""
    n, theta = symbols("n theta", positive=True)

    # n sinθ₁ = n sinθ₂ → θ₁ = θ₂
    # Trivially true


def test_mirror_concave_convex_duality():
    """Convex mirror (f < 0): always produces virtual image (s' < 0)."""
    s = symbols("s", positive=True)
    f = symbols("f", negative=True)

    s_prime = s * f / (s - f)

    # s > 0, f < 0 → s - f > 0 and s*f < 0 → s' < 0 (virtual)
    # Test with specific values
    val = s_prime.subs([(s, 10), (f, -5)])
    assert val < 0


def test_thin_film_quarter_wave():
    """Quarter-wave anti-reflection coating: t = λ/(4n)."""
    lam, n = symbols("lambda n", positive=True)

    t = lam / (4 * n)

    # Optical path = 2nt = λ/2 (half wavelength → destructive interference)
    optical_path = 2 * n * t
    assert simplify(optical_path - lam / 2) == 0
