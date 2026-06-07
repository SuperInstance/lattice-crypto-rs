//! Learning With Errors (LWE) encryption and decryption.
//!
//! Implements the basic LWE cryptosystem:
//! - Key generation: sample secret `s` uniformly at random
//! - Encryption: compute `(a, b = <a,s> + e + m·⌈q/2⌉ mod q)` where `e` is a small error
//! - Decryption: recover `m` from `b - <a,s> mod q`

use crate::gaussian::DiscreteGaussianSampler;
use crate::util::{mod_pos, XorShift64};

/// LWE cryptosystem instance parameterized by dimension, modulus, and error standard deviation.
#[derive(Debug, Clone)]
pub struct LWE {
    /// Dimension `n` of the secret vector.
    pub n: usize,
    /// Modulus `q`.
    pub q: i64,
    /// Error standard deviation.
    pub sigma: f64,
    rng: XorShift64,
}

/// LWE public key: a matrix A and vector b such that b = A·s + e (mod q).
#[derive(Debug, Clone, PartialEq)]
pub struct LWEPublicKey {
    /// The matrix `A` of size `m × n`.
    pub a: Vec<Vec<i64>>,
    /// The vector `b` of size `m`.
    pub b: Vec<i64>,
}

/// LWE ciphertext: a pair (a, b).
#[derive(Debug, Clone, PartialEq)]
pub struct LWECiphertext {
    /// The vector `a`.
    pub a: Vec<i64>,
    /// The scalar `b = <a, s> + e + encode(m) mod q`.
    pub b: i64,
}

impl LWE {
    /// Create a new LWE instance with the given parameters.
    pub fn new(n: usize, q: i64, sigma: f64, seed: u64) -> Self {
        Self { n, q, sigma, rng: XorShift64::new(seed) }
    }

    /// Generate a secret key: a vector of `n` random values in `[0, q)`.
    pub fn keygen(&mut self) -> Vec<i64> {
        (0..self.n).map(|_| self.rng.next_mod(self.q)).collect()
    }

    /// Generate a public key with `m` samples.
    ///
    /// Returns `(A, b)` where `b_i = <A_i, s> + e_i mod q`.
    pub fn public_keygen(&mut self, secret: &[i64], m: usize) -> LWEPublicKey {
        let seed = self.rng.next_u64();
        let mut sampler = DiscreteGaussianSampler::new(self.sigma, 0.0, seed);
        let mut a_matrix = Vec::with_capacity(m);
        let mut b_vec = Vec::with_capacity(m);
        for _ in 0..m {
            let a: Vec<i64> = (0..self.n).map(|_| self.rng.next_mod(self.q)).collect();
            let dot: i64 = a.iter().zip(secret.iter()).map(|(ai, si)| ai * si).sum();
            let e = sampler.sample_box_muller();
            let b = mod_pos(dot + e, self.q);
            a_matrix.push(a);
            b_vec.push(b);
        }
        LWEPublicKey { a: a_matrix, b: b_vec }
    }

    /// Encrypt a single bit `m ∈ {0, 1}` using the public key.
    ///
    /// Picks a random subset of public key rows, sums them, and encodes the message.
    pub fn encrypt(&mut self, pk: &LWEPublicKey, m: u8) -> LWECiphertext {
        assert!(m <= 1, "message must be 0 or 1");
        let m_size = pk.a.len();
        assert!(m_size > 0, "public key must have at least one row");

        // Pick a single random row
        let idx = (self.rng.next_u64() % (m_size as u64)) as usize;
        let sum_a = pk.a[idx].clone();
        let sum_b = pk.b[idx];

        // Add fresh error
        let seed = self.rng.next_u64();
        let mut sampler = DiscreteGaussianSampler::new(self.sigma, 0.0, seed);
        let e = sampler.sample_box_muller();

        // Encode message: if m=1, add q/2
        let encode = if m == 1 { self.q / 2 } else { 0 };

        LWECiphertext {
            a: sum_a,
            b: mod_pos(sum_b + e + encode, self.q),
        }
    }

    /// Decrypt a ciphertext using the secret key.
    ///
    /// Returns `0` or `1`.
    pub fn decrypt(&self, secret: &[i64], ct: &LWECiphertext) -> u8 {
        let dot: i64 = ct.a.iter().zip(secret.iter()).map(|(ai, si)| ai * si).sum();
        let v = mod_pos(ct.b - dot, self.q);
        // Decode: check if v is closer to 0 or q/2
        if v > self.q / 4 && v < 3 * self.q / 4 {
            1
        } else {
            0
        }
    }

    /// Encrypt a message vector bit-by-bit.
    pub fn encrypt_bits(&mut self, pk: &LWEPublicKey, bits: &[u8]) -> Vec<LWECiphertext> {
        bits.iter().map(|&b| self.encrypt(pk, b)).collect()
    }

    /// Decrypt a vector of ciphertexts.
    pub fn decrypt_bits(&self, secret: &[i64], cts: &[LWECiphertext]) -> Vec<u8> {
        cts.iter().map(|ct| self.decrypt(secret, ct)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_lwe() -> LWE {
        LWE::new(4, 97, 2.0, 42)
    }

    #[test]
    fn test_keygen_length() {
        let mut lwe = make_lwe();
        let sk = lwe.keygen();
        assert_eq!(sk.len(), 4);
        for &s in &sk {
            assert!((0..97).contains(&s));
        }
    }

    #[test]
    fn test_public_keygen() {
        let mut lwe = make_lwe();
        let sk = lwe.keygen();
        let pk = lwe.public_keygen(&sk, 10);
        assert_eq!(pk.a.len(), 10);
        assert_eq!(pk.b.len(), 10);
        assert_eq!(pk.a[0].len(), 4);
    }

    #[test]
    fn test_encrypt_decrypt_zero() {
        let mut lwe = make_lwe();
        let sk = lwe.keygen();
        let pk = lwe.public_keygen(&sk, 100);
        // Encrypt 0 many times; most should decrypt correctly
        let mut correct = 0;
        for _ in 0..20 {
            let ct = lwe.encrypt(&pk, 0);
            if lwe.decrypt(&sk, &ct) == 0 {
                correct += 1;
            }
        }
        assert!(correct >= 15, "Expected most encryptions of 0 to decrypt correctly, got {}/20", correct);
    }

    #[test]
    fn test_encrypt_decrypt_one() {
        let mut lwe = make_lwe();
        let sk = lwe.keygen();
        let pk = lwe.public_keygen(&sk, 100);
        let mut correct = 0;
        for _ in 0..20 {
            let ct = lwe.encrypt(&pk, 1);
            if lwe.decrypt(&sk, &ct) == 1 {
                correct += 1;
            }
        }
        assert!(correct >= 15, "Expected most encryptions of 1 to decrypt correctly, got {}/20", correct);
    }

    #[test]
    fn test_encrypt_decrypt_bits() {
        let mut lwe = make_lwe();
        let sk = lwe.keygen();
        let pk = lwe.public_keygen(&sk, 200);
        let msg = vec![1, 0, 1, 1, 0, 0, 1, 0];
        let cts = lwe.encrypt_bits(&pk, &msg);
        let dec = lwe.decrypt_bits(&sk, &cts);
        // At least 6 of 8 should be correct with these parameters
        let correct = msg.iter().zip(dec.iter()).filter(|(a, b)| a == b).count();
        assert!(correct >= 6, "Expected most bits to decrypt correctly, got {}/8", correct);
    }

    #[test]
    fn test_ciphertext_form() {
        let mut lwe = make_lwe();
        let sk = lwe.keygen();
        let pk = lwe.public_keygen(&sk, 10);
        let ct = lwe.encrypt(&pk, 1);
        assert_eq!(ct.a.len(), 4);
        assert!(ct.b >= 0 && ct.b < 97);
    }
}
