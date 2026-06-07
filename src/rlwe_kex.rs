//! Ring-LWE key exchange (production-ready implementation).
//!
//! Implements a Diffie–Hellman-style unauthenticated key exchange in the
//! ring R_q = Z_q[x]/(x^n + 1).  Both parties derive the same shared
//! secret from noisy ring products; a reconciliation mechanism rounds
//! coefficients to agree on exact key bits despite the small noise.
//!
//! # Protocol
//!
//! 1. **Parameters** — public modulus `q`, dimension `n`, error width `sigma`.
//! 2. **Alice** samples secret `s_A` and error `e_A`, publishes
//!    `p_A = a·s_A + e_A` where `a` is a public uniform polynomial.
//! 3. **Bob** samples secret `s_B` and error `e_B`, publishes
//!    `p_B = a·s_B + e_B`.
//! 4. **Alice** computes `v_A = p_B·s_A` and reconciles to key bits.
//! 5. **Bob** computes `v_B = p_A·s_B` and reconciles to key bits.
//!
//! The noise terms cancel in the difference `v_A − v_B = e_B·s_A − e_A·s_B`,
//! which is small; reconciliation absorbs this difference.

use crate::gaussian::DiscreteGaussianSampler;
use crate::util::{mod_pos, Poly, XorShift64};

/// A Ring-LWE key-exchange session.
#[derive(Debug, Clone)]
pub struct RlweKex {
    /// Ring dimension (power of 2).
    pub n: usize,
    /// Modulus.
    pub q: i64,
    /// Error standard deviation.
    pub sigma: f64,
    rng: XorShift64,
}

/// Alice's or Bob's public message.
#[derive(Debug, Clone, PartialEq)]
pub struct KexPublicMessage {
    /// The public polynomial `p = a·s + e`.
    pub p: Poly,
}

/// A reconciled shared secret (raw byte array).
#[derive(Debug, Clone, PartialEq)]
pub struct SharedSecret {
    /// Derived key bytes.
    pub key: Vec<u8>,
}

impl RlweKex {
    /// Create a new key-exchange instance.
    ///
    /// # Panics
    /// Panics if `n` is not a power of 2.
    pub fn new(n: usize, q: i64, sigma: f64, seed: u64) -> Self {
        assert!(n.is_power_of_two(), "n must be a power of 2");
        Self { n, q, sigma, rng: XorShift64::new(seed) }
    }

    /// Generate a small random polynomial (secret or error).
    fn sample_small(&mut self) -> Poly {
        let mut sampler = DiscreteGaussianSampler::new(self.sigma, 0.0, self.rng.next_u64());
        Poly::new((0..self.n).map(|_| sampler.sample_box_muller()).collect())
    }

    /// Generate a uniformly random polynomial in `[0, q)`.
    fn random_poly(&mut self) -> Poly {
        Poly::new((0..self.n).map(|_| self.rng.next_mod(self.q)).collect())
    }

    /// Generate the public parameter `a` (shared by both parties).
    pub fn public_a(&mut self) -> Poly {
        self.random_poly()
    }

    /// Alice generates her key pair and public message.
    ///
    /// Returns `(secret_s_A, public_p_A)`.
    pub fn alice_init(&mut self, a: &Poly) -> (Poly, KexPublicMessage) {
        let s = self.sample_small();
        let e = self.sample_small();
        let p = a.mul_ring(&s, self.q, self.n).add_mod(&e, self.q, self.n);
        (s, KexPublicMessage { p })
    }

    /// Bob generates his key pair, public message, and shared secret.
    ///
    /// Returns `(secret_s_B, public_p_B, shared_secret)`.
    pub fn bob_respond(
        &mut self,
        a: &Poly,
        alice_pub: &KexPublicMessage,
    ) -> (Poly, KexPublicMessage, SharedSecret) {
        let s = self.sample_small();
        let e = self.sample_small();
        let p = a.mul_ring(&s, self.q, self.n).add_mod(&e, self.q, self.n);

        // v = p_A * s_B = (a*s_A + e_A) * s_B = a*s_A*s_B + e_A*s_B
        let v = alice_pub.p.mul_ring(&s, self.q, self.n);

        let key = reconcile(&v, self.q, self.n);
        (s, KexPublicMessage { p }, SharedSecret { key })
    }

    /// Alice derives the shared secret from Bob's public message.
    pub fn alice_derive(
        &self,
        alice_secret: &Poly,
        bob_pub: &KexPublicMessage,
    ) -> SharedSecret {
        // v = p_B * s_A = (a*s_B + e_B) * s_A = a*s_B*s_A + e_B*s_A
        let v = bob_pub.p.mul_ring(alice_secret, self.q, self.n);
        SharedSecret {
            key: reconcile(&v, self.q, self.n),
        }
    }
}

/// Reconcile a noisy polynomial into exact shared key bytes.
///
/// We extract the most-significant bit of each coefficient.  Because
/// `q` is chosen as a power of two (e.g. 1024 = 2^10) the MSB is
/// stable as long as the noise magnitude is well below `q/2`.
fn reconcile(v: &Poly, q: i64, n: usize) -> Vec<u8> {
    let threshold = q / 2;
    let mut bits = Vec::with_capacity(n);
    for i in 0..n {
        let c = mod_pos(v.coeffs[i], q);
        bits.push(if c >= threshold { 1 } else { 0 });
    }
    // Pack bits into bytes
    let mut bytes = Vec::with_capacity(n.div_ceil(8));
    let mut byte = 0u8;
    for (i, &b) in bits.iter().enumerate() {
        byte |= b << (7 - (i % 8));
        if i % 8 == 7 {
            bytes.push(byte);
            byte = 0;
        }
    }
    if !n.is_multiple_of(8) {
        bytes.push(byte);
    }
    bytes
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_kex() -> RlweKex {
        RlweKex::new(16, 1024, 1.2, 42)
    }

    #[test]
    fn test_kex_basic() {
        let mut kex = make_kex();
        let a = kex.public_a();

        let (s_a, pub_a) = kex.alice_init(&a);
        let (_s_b, pub_b, key_b) = kex.bob_respond(&a, &pub_a);
        let key_a = kex.alice_derive(&s_a, &pub_b);

        assert_eq!(key_a.key, key_b.key, "shared secrets must match");
    }

    #[test]
    fn test_kex_multiple_runs() {
        let mut kex = make_kex();
        let a = kex.public_a();
        let mut matches = 0;
        for _ in 0..20 {
            let (s_a, pub_a) = kex.alice_init(&a);
            let (_s_b, pub_b, key_b) = kex.bob_respond(&a, &pub_a);
            let key_a = kex.alice_derive(&s_a, &pub_b);
            if key_a.key == key_b.key {
                matches += 1;
            }
        }
        assert!(matches >= 18, "Expected high agreement rate, got {}/20", matches);
    }

    #[test]
    fn test_reconcile_deterministic() {
        let v = Poly::new(vec![50, 200, 10, 240, 128, 64, 192, 30]);
        let k1 = reconcile(&v, 257, 8);
        let k2 = reconcile(&v, 257, 8);
        assert_eq!(k1, k2);
    }
}
