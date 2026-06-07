//! Discrete Gaussian sampling for lattice-based cryptography.
//!
//! Implements the Box-Muller transform for approximate Gaussian sampling
//! and a rejection sampling method for precise discrete Gaussian sampling.

use crate::util::XorShift64;

/// A discrete Gaussian sampler parameterized by standard deviation.
#[derive(Debug, Clone)]
pub struct DiscreteGaussianSampler {
    /// Standard deviation (sigma) of the Gaussian distribution.
    pub sigma: f64,
    /// Center of the distribution.
    pub center: f64,
    rng: XorShift64,
}

impl DiscreteGaussianSampler {
    /// Create a new sampler with the given standard deviation and center.
    ///
    /// # Panics
    /// Panics if `sigma` is not positive.
    pub fn new(sigma: f64, center: f64, seed: u64) -> Self {
        assert!(sigma > 0.0, "sigma must be positive");
        Self {
            sigma,
            center,
            rng: XorShift64::new(seed),
        }
    }

    /// Generate a single discrete Gaussian sample using the Box-Muller transform.
    ///
    /// Returns an integer sampled approximately from the discrete Gaussian distribution
    /// D_{Z,center,sigma}.
    pub fn sample_box_muller(&mut self) -> i64 {
        let u1 = self.uniform_open();
        let u2 = self.uniform_open();
        let mag = self.sigma * (-2.0 * u1.ln()).sqrt();
        let z = mag * (2.0 * std::f64::consts::PI * u2).cos();
        (z + self.center).round() as i64
    }

    /// Generate a single discrete Gaussian sample using rejection sampling.
    ///
    /// This is slower but produces samples from the exact discrete Gaussian distribution.
    /// Samples are drawn from the interval `[center - tail_cutoff * sigma, center + tail_cutoff * sigma]`.
    pub fn sample_rejection(&mut self, tail_cutoff: f64) -> i64 {
        let center = self.center;
        let sigma = self.sigma;
        let rho = |x: i64| -> f64 {
            let d = x as f64 - center;
            (-std::f64::consts::PI * d * d / (sigma * sigma)).exp()
        };
        let max_rho = rho(center.round() as i64);
        let lo = (center - tail_cutoff * sigma).floor() as i64;
        let hi = (center + tail_cutoff * sigma).ceil() as i64;
        loop {
            let x = lo + (self.rng.next_u64() % ((hi - lo + 1) as u64)) as i64;
            let u = self.uniform_open();
            if u * max_rho <= rho(x) {
                return x;
            }
        }
    }

    /// Generate `n` samples using the Box-Muller method.
    pub fn sample_n(&mut self, n: usize) -> Vec<i64> {
        (0..n).map(|_| self.sample_box_muller()).collect()
    }

    /// Compute the probability mass function at integer `x`.
    pub fn pmf(&self, x: i64) -> f64 {
        let d = x as f64 - self.center;
        (-std::f64::consts::PI * d * d / (self.sigma * self.sigma)).exp()
    }

    /// Uniform sample in (0, 1) for internal use.
    fn uniform_open(&mut self) -> f64 {
        // Map to (0, 1) excluding endpoints
        let raw = self.rng.next_u64();
        ((raw >> 11) as f64) / (1u64 << 53) as f64
    }
}

/// Compute the statistical distance between two discrete distributions.
///
/// Both slices should cover the same range of integers.
pub fn statistical_distance(p: &[f64], q: &[f64]) -> f64 {
    assert_eq!(p.len(), q.len(), "distributions must have same length");
    let total: f64 = p.iter().zip(q.iter()).map(|(a, b)| (a - b).abs()).sum();
    total / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sampler_creation() {
        let sampler = DiscreteGaussianSampler::new(3.0, 0.0, 42);
        assert!((sampler.sigma - 3.0).abs() < 1e-10);
        assert!((sampler.center - 0.0).abs() < 1e-10);
    }

    #[test]
    #[should_panic(expected = "sigma must be positive")]
    fn test_sampler_zero_sigma() {
        DiscreteGaussianSampler::new(0.0, 0.0, 42);
    }

    #[test]
    fn test_box_muller_samples() {
        let mut sampler = DiscreteGaussianSampler::new(3.0, 0.0, 42);
        let samples: Vec<i64> = (0..100).map(|_| sampler.sample_box_muller()).collect();
        // Most samples should be within 4*sigma of center
        let within = samples.iter().filter(|&&x| x.abs() <= 12).count();
        assert!(within > 90, "Expected most samples within 4*sigma, got {}/100", within);
    }

    #[test]
    fn test_rejection_samples() {
        let mut sampler = DiscreteGaussianSampler::new(2.0, 0.0, 123);
        let samples: Vec<i64> = (0..50).map(|_| sampler.sample_rejection(6.0)).collect();
        let within = samples.iter().filter(|&&x| x.abs() <= 12).count();
        assert!(within == 50, "All rejection samples should be within tail_cutoff");
    }

    #[test]
    fn test_sample_n() {
        let mut sampler = DiscreteGaussianSampler::new(1.0, 5.0, 99);
        let samples = sampler.sample_n(20);
        assert_eq!(samples.len(), 20);
    }

    #[test]
    fn test_pmf_center() {
        let sampler = DiscreteGaussianSampler::new(1.0, 0.0, 42);
        let pmf_0 = sampler.pmf(0);
        let pmf_5 = sampler.pmf(5);
        assert!(pmf_0 > pmf_5, "PMF at center should be higher than far away");
        assert!((pmf_0 - 1.0).abs() < 1e-10, "PMF at center should be ~1.0");
    }

    #[test]
    fn test_statistical_distance_identical() {
        let p = [0.25, 0.25, 0.25, 0.25];
        assert!((statistical_distance(&p, &p) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_statistical_distance_different() {
        let p = [1.0, 0.0, 0.0, 0.0];
        let q = [0.0, 0.0, 0.0, 1.0];
        assert!((statistical_distance(&p, &q) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_centered_distribution() {
        let mut sampler = DiscreteGaussianSampler::new(2.0, 10.0, 42);
        let samples: Vec<i64> = (0..1000).map(|_| sampler.sample_box_muller()).collect();
        let mean = samples.iter().sum::<i64>() as f64 / samples.len() as f64;
        assert!((mean - 10.0).abs() < 1.0, "Mean should be close to center 10.0, got {}", mean);
    }
}
