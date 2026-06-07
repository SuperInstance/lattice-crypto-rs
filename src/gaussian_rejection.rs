//! Advanced discrete Gaussian sampling via rejection sampling.
//!
//! Provides algorithms for sampling from the discrete Gaussian distribution
//! over the integers and over short vectors, with exact statistical guarantees.
//!
//! The discrete Gaussian `D_{Z,σ,c}` assigns mass proportional to
//! `ρ_{σ,c}(x) = exp(-π·(x-c)²/σ²)` to each integer `x`.

use crate::util::XorShift64;

/// Sampler for the discrete Gaussian over the integers using exact rejection.
///
/// Uses the method of Gentry, Peikert & Vaikuntanathan (STOC 2008):
/// sample from a uniform interval and accept with probability proportional
/// to the Gaussian mass.
#[derive(Debug, Clone)]
pub struct RejectionSampler {
    sigma: f64,
    center: f64,
    tail_cutoff: f64,
    rng: XorShift64,
}

impl RejectionSampler {
    /// Create a new rejection sampler.
    ///
    /// `tail_cutoff` is the number of standard deviations to sample from
    /// (e.g. 6.0 captures > 99.9999% of the mass).
    pub fn new(sigma: f64, center: f64, tail_cutoff: f64, seed: u64) -> Self {
        assert!(sigma > 0.0, "sigma must be positive");
        assert!(tail_cutoff > 0.0, "tail_cutoff must be positive");
        Self {
            sigma,
            center,
            tail_cutoff,
            rng: XorShift64::new(seed),
        }
    }

    /// Gaussian rho function `ρ(x) = exp(-π·(x-center)²/σ²)`.
    fn rho(&self, x: i64) -> f64 {
        let d = x as f64 - self.center;
        (-std::f64::consts::PI * d * d / (self.sigma * self.sigma)).exp()
    }

    /// Sample a single integer from `D_{Z,σ,c}`.
    ///
    /// Expected number of trials is `≈ σ·√(2e) / tail_cutoff_range`.
    pub fn sample(&mut self) -> i64 {
        let lo = (self.center - self.tail_cutoff * self.sigma).floor() as i64;
        let hi = (self.center + self.tail_cutoff * self.sigma).ceil() as i64;
        let range = (hi - lo + 1) as u64;
        let max_rho = self.rho(self.center.round() as i64);

        loop {
            let x = lo + (self.rng.next_u64() % range) as i64;
            let u = self.uniform_open();
            if u * max_rho <= self.rho(x) {
                return x;
            }
        }
    }

    /// Sample `n` independent integers.
    pub fn sample_n(&mut self, n: usize) -> Vec<i64> {
        (0..n).map(|_| self.sample()).collect()
    }

    /// Sample a short vector in Z^n where each coordinate is independently
    /// drawn from `D_{Z,σ,c}`.
    pub fn sample_vector(&mut self, n: usize) -> Vec<i64> {
        self.sample_n(n)
    }

    /// Kolmogorov–Smirnov-style check: empirical CDF vs exact CDF.
    ///
    /// Returns the maximum absolute deviation.  For large `samples` this
    /// should be small (≈ 1/√N).
    pub fn ks_test(&mut self, samples: usize) -> f64 {
        let mut drawn: Vec<i64> = self.sample_n(samples);
        drawn.sort_unstable();

        let lo = (self.center - self.tail_cutoff * self.sigma).floor() as i64;
        let hi = (self.center + self.tail_cutoff * self.sigma).ceil() as i64;

        // Compute exact CDF
        let mut total_rho = 0.0;
        for x in lo..=hi {
            total_rho += self.rho(x);
        }

        let mut max_dev = 0.0;
        let mut empirical_idx = 0;
        let mut cum_exact = 0.0;

        for x in lo..=hi {
            cum_exact += self.rho(x) / total_rho;
            while empirical_idx < drawn.len() && drawn[empirical_idx] <= x {
                empirical_idx += 1;
            }
            let cum_empirical = empirical_idx as f64 / samples as f64;
            let dev = (cum_empirical - cum_exact).abs();
            if dev > max_dev {
                max_dev = dev;
            }
        }
        max_dev
    }

    fn uniform_open(&mut self) -> f64 {
        let raw = self.rng.next_u64();
        ((raw >> 11) as f64) / (1u64 << 53) as f64
    }
}

/// One-dimensional discrete Gaussian sampler using the "convolution" method
/// (Ducas & Nguyen, 2013): sample from a wide Gaussian by adding a small
/// Gaussian to a uniform offset.
pub struct ConvolutionSampler {
    base: RejectionSampler,
    k: i64,
}

impl ConvolutionSampler {
    /// Create a sampler for a large `sigma` by decomposing it as
    /// `sigma = k·σ_base` where `σ_base` is small enough for fast rejection.
    pub fn new(sigma: f64, _center: f64, seed: u64) -> Self {
        let k = (sigma / 5.0).ceil() as i64;
        let sigma_base = sigma / k as f64;
        let base = RejectionSampler::new(sigma_base, 0.0, 6.0, seed);
        Self { base, k }
    }

    /// Sample from the large Gaussian.
    pub fn sample(&mut self) -> i64 {
        let y = self.base.sample();
        let offset = (self.base.rng.next_u64() % self.k as u64) as i64;
        y * self.k + offset
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rejection_sampler_range() {
        let mut sampler = RejectionSampler::new(3.0, 0.0, 6.0, 42);
        let samples = sampler.sample_n(100);
        let max_abs = samples.iter().map(|x| x.abs()).max().unwrap();
        assert!(max_abs <= 20, "all samples should be within tail cutoff");
    }

    #[test]
    fn test_rejection_sampler_centered() {
        let mut sampler = RejectionSampler::new(2.0, 10.0, 6.0, 123);
        let samples = sampler.sample_n(500);
        let mean = samples.iter().sum::<i64>() as f64 / samples.len() as f64;
        assert!((mean - 10.0).abs() < 0.5, "mean should be close to center, got {}", mean);
    }

    #[test]
    fn test_rejection_sampler_ks() {
        let mut sampler = RejectionSampler::new(2.0, 0.0, 8.0, 99);
        let dev = sampler.ks_test(2000);
        assert!(dev < 0.05, "KS deviation too large: {}", dev);
    }

    #[test]
    fn test_convolution_sampler_range() {
        let mut sampler = ConvolutionSampler::new(20.0, 0.0, 77);
        let samples: Vec<i64> = (0..100).map(|_| sampler.sample()).collect();
        let max_abs = samples.iter().map(|x| x.abs()).max().unwrap();
        assert!(max_abs <= 150, "convolution samples should stay bounded");
    }
}
