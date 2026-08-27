//! MAGSAC++ noise-marginalized weighting and the special functions it relies on.
//!
//! The MAGSAC++ weight of a residual `r`, given an upper noise bound `σ_max` and a residual
//! distribution with `ν` degrees of freedom, is (Barath et al., 2020):
//!
//! ```text
//! w(r) = C(ν) · 2^((ν-1)/2) / σ_max · [ Γ((ν-1)/2, r²/(2σ_max²)) − Γ((ν-1)/2, k²/2) ]
//! ```
//!
//! for `0 ≤ r ≤ k·σ_max` and `0` otherwise, where `Γ(a, x)` is the upper incomplete gamma
//! function, `C(ν) = 1 / (2^(ν/2−1) · Γ(ν/2))` is the χ-distribution normalization constant, and
//! `k` is the `0.99` quantile of the χ-distribution with `ν` degrees of freedom. Subtracting the
//! value at the cutoff makes the weight fall continuously to zero at `r = k·σ_max`. On `(0, k·σ_max)`
//! the weight is positive and strictly decreasing.
//!
//! For every geometric primitive in this crate the residual is a Euclidean point-to-model distance
//! embedded in `D`-dimensional space, so `ν = D ≥ 2` and the exponent `a = (ν−1)/2 ≥ 0.5`, which
//! keeps `Γ(a, 0) = Γ(a)` finite and the weight bounded as `r → 0`.

use crate::common::vec_f64::median;

const LANCZOS_G: f64 = 7.0;

/// Lanczos approximation coefficients (g = 7, n = 9). These are the standard published constants;
/// the extra digits are harmless and simply round to the nearest `f64`.
#[allow(clippy::excessive_precision)]
const LANCZOS_COEFFS: [f64; 9] = [
    0.999_999_999_999_809_93,
    676.520_368_121_885_1,
    -1_259.139_216_722_402_8,
    771.323_428_777_653_1,
    -176.615_029_162_140_6,
    12.507_343_278_686_905,
    -0.138_571_095_265_720_12,
    9.984_369_578_019_572e-6,
    1.505_632_735_149_311_6e-7,
];

/// Natural logarithm of the gamma function via the Lanczos approximation. Valid for `x > 0`.
pub(crate) fn ln_gamma(x: f64) -> f64 {
    // Reflection is not needed: every argument used here is positive (a ≥ 0.5).
    let x = x - 1.0;
    let mut a = LANCZOS_COEFFS[0];
    let t = x + LANCZOS_G + 0.5;
    for (i, &c) in LANCZOS_COEFFS.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// The regularized lower incomplete gamma function `P(a, x) = γ(a, x) / Γ(a)` for `a > 0`, `x ≥ 0`.
///
/// Uses the series expansion for `x < a + 1` and the complement of the continued fraction
/// otherwise (Numerical Recipes, `gammp`/`gammq`).
pub(crate) fn reg_lower_gamma(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x < a + 1.0 {
        gamma_series(a, x)
    } else {
        1.0 - gamma_cont_frac(a, x)
    }
}

/// The regularized upper incomplete gamma function `Q(a, x) = Γ(a, x) / Γ(a)`.
pub(crate) fn reg_upper_gamma(a: f64, x: f64) -> f64 {
    1.0 - reg_lower_gamma(a, x)
}

/// The (unregularized) upper incomplete gamma function `Γ(a, x)`.
pub(crate) fn upper_gamma(a: f64, x: f64) -> f64 {
    reg_upper_gamma(a, x) * ln_gamma(a).exp()
}

/// Series representation of `P(a, x)`, converging quickly for `x < a + 1`.
fn gamma_series(a: f64, x: f64) -> f64 {
    const MAX_ITER: usize = 200;
    const EPS: f64 = 1e-15;

    let mut ap = a;
    let mut sum = 1.0 / a;
    let mut del = sum;
    for _ in 0..MAX_ITER {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * EPS {
            break;
        }
    }
    sum * (-x + a * x.ln() - ln_gamma(a)).exp()
}

/// Continued-fraction representation of `Q(a, x)`, converging quickly for `x ≥ a + 1`.
fn gamma_cont_frac(a: f64, x: f64) -> f64 {
    const MAX_ITER: usize = 200;
    const EPS: f64 = 1e-15;
    const TINY: f64 = 1e-300;

    let mut b = x + 1.0 - a;
    let mut c = 1.0 / TINY;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..MAX_ITER {
        let an = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < TINY {
            d = TINY;
        }
        c = b + an / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < EPS {
            break;
        }
    }
    (-x + a * x.ln() - ln_gamma(a)).exp() * h
}

/// The `p`-quantile of the χ²-distribution with `dof` degrees of freedom, found by bisection on the
/// regularized lower incomplete gamma CDF `P(dof/2, x/2)`.
fn chi2_quantile(p: f64, dof: usize) -> f64 {
    let a = dof as f64 / 2.0;
    let target = p;
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    // Expand the upper bound until the CDF exceeds the target.
    while reg_lower_gamma(a, hi / 2.0) < target && hi < 1e6 {
        hi *= 2.0;
    }
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if reg_lower_gamma(a, mid / 2.0) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// A precomputed MAGSAC++ weighting function for a fixed `σ_max` and residual degrees of freedom.
///
/// All constants that depend only on `σ_max` and `dof` are computed once so that per-point weight
/// evaluation reduces to a single incomplete-gamma call.
pub(crate) struct MagsacWeight {
    /// Squared denominator `2 · σ_max²` used to form `x = r² / (2 σ_max²)`.
    two_sigma_sq: f64,
    /// The exponent `a = (ν − 1) / 2` of the incomplete gamma function.
    a: f64,
    /// The overall constant `C(ν) · 2^((ν−1)/2) / σ_max`.
    coeff: f64,
    /// The residual cutoff `k · σ_max`; residuals at or beyond it have zero weight.
    cutoff: f64,
    /// The incomplete gamma value at the cutoff, subtracted so the weight reaches zero there.
    gamma_at_cutoff: f64,
}

impl MagsacWeight {
    pub(crate) fn new(sigma_max: f64, dof: usize) -> Self {
        let nu = dof as f64;
        let a = (nu - 1.0) / 2.0;
        let k = chi2_quantile(0.99, dof).sqrt();
        let two_sigma_sq = 2.0 * sigma_max * sigma_max;

        let c_nu = 1.0 / (2.0_f64.powf(nu / 2.0 - 1.0) * ln_gamma(nu / 2.0).exp());
        let coeff = c_nu * 2.0_f64.powf((nu - 1.0) / 2.0) / sigma_max;

        let x_cutoff = k * k / 2.0;
        let gamma_at_cutoff = upper_gamma(a, x_cutoff);

        Self {
            two_sigma_sq,
            a,
            coeff,
            cutoff: k * sigma_max,
            gamma_at_cutoff,
        }
    }

    /// The residual magnitude at or beyond which a point receives zero weight (`k · σ_max`).
    pub(crate) fn cutoff(&self) -> f64 {
        self.cutoff
    }

    /// The MAGSAC++ weight of a residual magnitude `r` (`r ≥ 0`).
    pub(crate) fn weight(&self, r: f64) -> f64 {
        if r >= self.cutoff {
            return 0.0;
        }
        let x = r * r / self.two_sigma_sq;
        (self.coeff * (upper_gamma(self.a, x) - self.gamma_at_cutoff)).max(0.0)
    }
}

/// The scale factor that turns a median absolute deviation into a consistent estimate of the
/// standard deviation of normally distributed data.
pub(crate) const MAD_TO_SIGMA: f64 = 1.4826;

/// Estimates a MAGSAC++ `sigma_max` from a set of residuals via the median absolute deviation.
///
/// MAD is used rather than the standard deviation because it is insensitive to the gross outliers
/// the robust weighting exists to suppress: contaminating up to half the data cannot move it
/// arbitrarily, whereas a single distant point can dominate a standard deviation.
///
/// The residuals are expected to be *signed*, so that a well-fitted set is centered near zero.
/// Feeding unsigned distances in works, but measures the spread about the typical distance rather
/// than about zero, which biases the estimate low.
///
/// Returns `None` when the spread is zero or non-finite, which happens when the fit is already
/// essentially exact and there is nothing to reweight.
pub(crate) fn estimate_sigma_max(residuals: &[f64]) -> Option<f64> {
    let center = median(residuals)?;
    let deviations: Vec<f64> = residuals.iter().map(|r| (r - center).abs()).collect();
    let sigma = MAD_TO_SIGMA * median(&deviations)?;

    (sigma.is_finite() && sigma > 0.0).then_some(sigma)
}

#[cfg(test)]
mod tests {
    #[test]
    fn sigma_estimate_is_robust_to_outliers() {
        // Eleven values with unit spread, plus one enormous outlier. A standard deviation would
        // be dominated by the outlier; the MAD-based estimate should barely notice it.
        let mut values: Vec<f64> = (-5..=5).map(|i| i as f64).collect();
        let clean = estimate_sigma_max(&values).unwrap();

        values.push(10_000.0);
        let contaminated = estimate_sigma_max(&values).unwrap();

        assert!(
            (contaminated - clean).abs() < 0.5 * clean,
            "outlier moved the estimate too far: {clean} -> {contaminated}"
        );
    }

    #[test]
    fn sigma_estimate_rejects_degenerate_spread() {
        // Every residual identical means no spread to estimate from.
        assert_eq!(estimate_sigma_max(&[2.0; 10]), None);
        assert_eq!(estimate_sigma_max(&[]), None);
    }
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn ln_gamma_known_values() {
        // Γ(1) = 1, Γ(2) = 1, Γ(1/2) = √π, Γ(5) = 24
        assert_relative_eq!(ln_gamma(1.0).exp(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(ln_gamma(2.0).exp(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(
            ln_gamma(0.5).exp(),
            std::f64::consts::PI.sqrt(),
            epsilon = 1e-10
        );
        assert_relative_eq!(ln_gamma(5.0).exp(), 24.0, epsilon = 1e-9);
    }

    #[test]
    fn upper_gamma_known_values() {
        // Γ(1, x) = e^{-x}
        for &x in &[0.25_f64, 1.0, 3.0, 7.0] {
            assert_relative_eq!(upper_gamma(1.0, x), (-x).exp(), epsilon = 1e-9);
        }
        // Γ(a, 0) = Γ(a)
        assert_relative_eq!(
            upper_gamma(0.5, 0.0),
            std::f64::consts::PI.sqrt(),
            epsilon = 1e-9
        );
        // Regularized values sum to one
        assert_relative_eq!(
            reg_lower_gamma(2.5, 1.3) + reg_upper_gamma(2.5, 1.3),
            1.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn chi2_quantile_known_values() {
        // For dof = 2 the χ² CDF is 1 − e^{−x/2}, so the 0.99 quantile is −2 ln(0.01).
        let expected = -2.0 * 0.01_f64.ln();
        assert_relative_eq!(chi2_quantile(0.99, 2), expected, epsilon = 1e-6);
    }

    #[test]
    fn weight_is_positive_and_decreasing_then_zero() {
        let w = MagsacWeight::new(0.1, 2);
        let cutoff = w.cutoff();

        // Strictly positive and strictly decreasing over the interior of the support. We stop just
        // short of the cutoff, where floating-point cancellation between the two incomplete-gamma
        // terms legitimately clamps the weight to zero.
        let mut prev = f64::INFINITY;
        let mut r = 0.0;
        while r < 0.95 * cutoff {
            let val = w.weight(r);
            assert!(val > 0.0, "weight should be positive inside the support");
            assert!(val < prev, "weight should be strictly decreasing");
            prev = val;
            r += cutoff / 50.0;
        }

        // The weight decays to zero at the cutoff and stays there beyond it.
        assert_eq!(w.weight(cutoff), 0.0);
        assert_eq!(w.weight(cutoff * 2.0), 0.0);
    }
}
