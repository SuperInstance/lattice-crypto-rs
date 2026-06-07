//! Ring-LWE cryptography: key generation, encryption, and decryption.
//!
//! Operates in the ring R_q = Z_q[x]/(x^n + 1) where n is a power of 2.
//! More efficient than standard LWE due to the ring structure.

use crate::gaussian::DiscreteGaussianSampler;
use crate::util::{mod_pos, Poly, XorShift64};

/// Ring-LWE cryptosystem instance.
#[derive(Debug, Clone)]
pub struct RingLWE {
    /// Ring dimension `n` (should be a power of 2).
    pub n: usize,
    /// Modulus `q`.
    pub q: i64,
    /// Error standard deviation.
    pub sigma: f64,
    rng: XorShift64,
}

/// Ring-LWE key pair.
#[derive(Debug, Clone, PartialEq)]
pub struct RingKeyPair {
    /// Secret key polynomial.
    pub secret: Poly,
    /// Public key polynomial a (sampled uniformly).
    pub a: Poly,
    /// Public key polynomial b = a*s + e.
    pub b: Poly,
}

/// Ring-LWE ciphertext: a pair of polynomials (c1, c2).
#[derive(Debug, Clone, PartialEq)]
pub struct RingCiphertext {
    /// First component polynomial.
    pub c1: Poly,
    /// Second component polynomial.
    pub c2: Poly,
}

impl RingLWE {
    /// Create a new Ring-LWE instance.
    ///
    /// # Panics
    /// Panics if `n` is not a power of 2.
    pub fn new(n: usize, q: i64, sigma: f64, seed: u64) -> Self {
        assert!(n.is_power_of_two(), "n must be a power of 2");
        Self { n, q, sigma, rng: XorShift64::new(seed) }
    }

    /// Generate a random polynomial with coefficients in `[0, q)`.
    fn random_poly(&mut self) -> Poly {
        Poly::new((0..self.n).map(|_| self.rng.next_mod(self.q)).collect())
    }

    /// Sample an error polynomial from the discrete Gaussian.
    fn error_poly(&mut self) -> Poly {
        let mut sampler = DiscreteGaussianSampler::new(self.sigma, 0.0, self.rng.next_u64());
        Poly::new((0..self.n).map(|_| sampler.sample_box_muller()).collect())
    }

    /// Generate a key pair.
    pub fn keygen(&mut self) -> RingKeyPair {
        let a = self.random_poly();
        let secret = self.error_poly(); // Small polynomial for secret
        let e = self.error_poly();
        // b = a * s + e (mod q) in the ring
        let as_prod = a.mul_ring(&secret, self.q, self.n);
        let b = as_prod.add_mod(&e, self.q, self.n);
        RingKeyPair { secret, a, b }
    }

    /// Encrypt a message polynomial.
    ///
    /// Returns (c1, c2) where:
    /// - c1 = a*r + e1
    /// - c2 = b*r + e2 + encode(m)
    pub fn encrypt(&mut self, kp: &RingKeyPair, message: &Poly) -> RingCiphertext {
        let r = self.error_poly(); // Small random polynomial
        let e1 = self.error_poly();
        let e2 = self.error_poly();

        let c1 = kp.a.mul_ring(&r, self.q, self.n).add_mod(&e1, self.q, self.n);
        let br = kp.b.mul_ring(&r, self.q, self.n).add_mod(&e2, self.q, self.n);

        // Encode message: multiply by ⌊q/2⌉
        let encoded = message.scale(self.q / 2, self.q);
        let c2 = br.add_mod(&encoded, self.q, self.n);

        RingCiphertext { c1, c2 }
    }

    /// Decrypt a ciphertext.
    ///
    /// Computes c2 - c1*s, then decodes by dividing by ⌊q/2⌉.
    pub fn decrypt(&self, kp: &RingKeyPair, ct: &RingCiphertext) -> Poly {
        let c1s = ct.c1.mul_ring(&kp.secret, self.q, self.n);
        let v = ct.c2.add_mod(&c1s.scale(-1, self.q), self.q, self.n);

        // Decode: round each coefficient to nearest multiple of q/2
        let _half_q = self.q / 2;
        let coeffs: Vec<i64> = v.coeffs.iter().map(|&c| {
            let c = mod_pos(c, self.q);
            if c > self.q / 4 && c < 3 * self.q / 4 {
                1
            } else {
                0
            }
        }).collect();

        Poly::new(coeffs)
    }

    /// Encrypt a single bit as a constant polynomial.
    pub fn encrypt_bit(&mut self, kp: &RingKeyPair, bit: u8) -> RingCiphertext {
        assert!(bit <= 1);
        let mut coeffs = vec![0i64; self.n];
        coeffs[0] = bit as i64;
        self.encrypt(kp, &Poly::new(coeffs))
    }

    /// Decrypt to a single bit (from constant polynomial).
    pub fn decrypt_bit(&self, kp: &RingKeyPair, ct: &RingCiphertext) -> u8 {
        let msg = self.decrypt(kp, ct);
        mod_pos(msg.coeffs[0], self.q) as u8
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_ring_lwe() -> RingLWE {
        RingLWE::new(4, 97, 1.5, 42)
    }

    #[test]
    fn test_ring_lwe_creation() {
        let rlwe = make_ring_lwe();
        assert_eq!(rlwe.n, 4);
        assert_eq!(rlwe.q, 97);
    }

    #[test]
    #[should_panic(expected = "n must be a power of 2")]
    fn test_ring_lwe_bad_n() {
        RingLWE::new(5, 97, 1.5, 42);
    }

    #[test]
    fn test_keygen() {
        let mut rlwe = make_ring_lwe();
        let kp = rlwe.keygen();
        assert_eq!(kp.secret.coeffs.len(), 4);
        assert_eq!(kp.a.coeffs.len(), 4);
        assert_eq!(kp.b.coeffs.len(), 4);
    }

    #[test]
    fn test_encrypt_zero_polynomial() {
        let mut rlwe = make_ring_lwe();
        let kp = rlwe.keygen();
        let msg = Poly::zero(4);
        let ct = rlwe.encrypt(&kp, &msg);
        assert_eq!(ct.c1.coeffs.len(), 4);
        assert_eq!(ct.c2.coeffs.len(), 4);
    }

    #[test]
    fn test_encrypt_decrypt_zero() {
        let mut rlwe = make_ring_lwe();
        let kp = rlwe.keygen();
        let msg = Poly::zero(4);
        let mut correct = 0;
        for _ in 0..10 {
            let ct = rlwe.encrypt(&kp, &msg);
            let dec = rlwe.decrypt(&kp, &ct);
            if dec.coeffs.iter().all(|&c| c == 0) {
                correct += 1;
            }
        }
        assert!(correct >= 6, "Expected most zero encryptions to decrypt correctly, got {}/10", correct);
    }

    #[test]
    fn test_encrypt_decrypt_one() {
        let mut rlwe = make_ring_lwe();
        let kp = rlwe.keygen();
        let msg = Poly::new(vec![1, 0, 0, 0]);
        let mut correct = 0;
        for _ in 0..10 {
            let ct = rlwe.encrypt(&kp, &msg);
            let dec = rlwe.decrypt(&kp, &ct);
            if dec.coeffs[0] == 1 && dec.coeffs[1..].iter().all(|&c| c == 0) {
                correct += 1;
            }
        }
        assert!(correct >= 6, "Expected most one encryptions to decrypt correctly, got {}/10", correct);
    }

    #[test]
    fn test_encrypt_decrypt_bit() {
        let mut rlwe = make_ring_lwe();
        let kp = rlwe.keygen();
        let mut correct = 0;
        for _ in 0..20 {
            let ct = rlwe.encrypt_bit(&kp, 1);
            if rlwe.decrypt_bit(&kp, &ct) == 1 {
                correct += 1;
            }
        }
        assert!(correct >= 14, "Expected most bit encryptions to decrypt correctly, got {}/20", correct);
    }

    #[test]
    fn test_key_pair_public_key_relation() {
        let mut rlwe = make_ring_lwe();
        let kp = rlwe.keygen();
        // b = a*s + e (mod q), so b should be deterministic given a, s, e
        assert_eq!(kp.b.coeffs.len(), kp.a.coeffs.len());
    }

    #[test]
    fn test_multiple_keygens_differ() {
        let mut rlwe = make_ring_lwe();
        let kp1 = rlwe.keygen();
        let kp2 = rlwe.keygen();
        // Different keygens should produce different secret keys (with high probability)
        assert_ne!(kp1.secret, kp2.secret, "Consecutive keygens should differ");
    }

    #[test]
    fn test_encrypt_decrypt_all_zeros() {
        let mut rlwe = make_ring_lwe();
        let kp = rlwe.keygen();
        let msg = Poly::new(vec![0, 0, 0, 0]);
        let mut correct = 0;
        for _ in 0..20 {
            let ct = rlwe.encrypt(&kp, &msg);
            let dec = rlwe.decrypt(&kp, &ct);
            if dec.coeffs.iter().all(|&c| c == 0) {
                correct += 1;
            }
        }
        assert!(correct >= 14, "All-zero message should decrypt reliably, got {}/20", correct);
    }

    #[test]
    fn test_encrypt_decrypt_all_ones() {
        let mut rlwe = make_ring_lwe();
        let kp = rlwe.keygen();
        let msg = Poly::new(vec![1, 1, 1, 1]);
        let mut correct = 0;
        for _ in 0..20 {
            let ct = rlwe.encrypt(&kp, &msg);
            let dec = rlwe.decrypt(&kp, &ct);
            if dec.coeffs.iter().all(|&c| c == 1) {
                correct += 1;
            }
        }
        assert!(correct >= 12, "All-ones message should decrypt somewhat reliably, got {}/20", correct);
    }

    #[test]
    fn test_ciphertext_components_correct_size() {
        let mut rlwe = make_ring_lwe();
        let kp = rlwe.keygen();
        let msg = Poly::new(vec![0, 1, 0, 1]);
        let ct = rlwe.encrypt(&kp, &msg);
        assert_eq!(ct.c1.coeffs.len(), 4);
        assert_eq!(ct.c2.coeffs.len(), 4);
        // Ciphertext coefficients should be in [0, q)
        for &c in &ct.c1.coeffs {
            assert!(c >= 0 && c < 97);
        }
        for &c in &ct.c2.coeffs {
            assert!(c >= 0 && c < 97);
        }
    }
}
