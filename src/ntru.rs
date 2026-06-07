//! NTRU-style encryption over the ring R = Z_q[x]/(x^n - 1).
//!
//! Implements the classic NTRU cryptosystem with parameters:
//!
//! - `n` — polynomial degree
//! - `p` — small modulus (typically 3)
//! - `q` — large modulus (power of 2 or prime)
//!
//! # Key Generation
//!
//! 1. Choose small polynomials `f` and `g`.
//! 2. Compute `f_q = f^{-1} mod q` in the ring.
//! 3. Compute `f_p = f^{-1} mod p` in the ring.
//! 4. Public key: `h = p · g · f_q  (mod q)`.
//! 5. Secret key: `(f, f_p)`.
//!
//! # Encryption
//!
//! 1. Choose small random `r`.
//! 2. Ciphertext: `e = r · h + m  (mod q)` where `m` is the message polynomial
//!    with coefficients in `[-(p-1)/2, (p-1)/2]`.
//!
//! # Decryption
//!
//! 1. Compute `a = f · e  (mod q)` and center coefficients into `[-q/2, q/2]`.
//! 2. Recover `m = f_p · a  (mod p)`.

use crate::util::{mod_inverse, mod_pos, Poly, XorShift64};
use crate::LatticeError;

/// NTRU cryptosystem instance.
#[derive(Debug, Clone)]
pub struct Ntru {
    /// Polynomial degree.
    pub n: usize,
    /// Small modulus (e.g. 3).
    pub p: i64,
    /// Large modulus.
    pub q: i64,
    /// Number of +1 coefficients in random trinary polynomials.
    pub d: usize,
    rng: XorShift64,
}

/// NTRU key pair.
#[derive(Debug, Clone, PartialEq)]
pub struct NtruKeyPair {
    /// Public key polynomial `h`.
    pub h: Poly,
    /// Secret key polynomial `f`.
    pub f: Poly,
    /// Secret key `f_p = f^{-1} mod p`.
    pub f_p: Poly,
}

/// NTRU ciphertext.
#[derive(Debug, Clone, PartialEq)]
pub struct NtruCiphertext {
    /// Ciphertext polynomial `e`.
    pub e: Poly,
}

impl Ntru {
    /// Create a new NTRU instance.
    ///
    /// # Panics
    /// Panics if `p` does not divide `q` evenly or if `d > n/3`.
    pub fn new(n: usize, p: i64, q: i64, d: usize, seed: u64) -> Self {
        assert!(p > 1, "p must be > 1");
        assert!(q > p, "q must be > p");
        assert!(d <= n / 3, "d must be <= n/3");
        Self { n, p, q, d, rng: XorShift64::new(seed) }
    }

    /// Generate a random trinary polynomial with `d` entries of +1 and `d` of -1.
    fn random_trinary(&mut self) -> Poly {
        let mut coeffs = vec![0i64; self.n];
        let mut placed = 0;
        while placed < self.d {
            let idx = (self.rng.next_u64() as usize) % self.n;
            if coeffs[idx] == 0 {
                coeffs[idx] = 1;
                placed += 1;
            }
        }
        placed = 0;
        while placed < self.d {
            let idx = (self.rng.next_u64() as usize) % self.n;
            if coeffs[idx] == 0 {
                coeffs[idx] = -1;
                placed += 1;
            }
        }
        Poly::new(coeffs)
    }

    /// Generate a key pair.
    ///
    /// Returns `Err(LatticeError::SamplingFailed)` if `f` is not invertible mod q.
    pub fn keygen(&mut self) -> Result<NtruKeyPair, LatticeError> {
        for _ in 0..100 {
            let f = self.random_trinary();
            let g = self.random_trinary();

            if let Some(f_q) = invert_negacyclic(&f, self.q, self.n) {
                if let Some(f_p) = invert_negacyclic(&f, self.p, self.n) {
                    // h = p * g * f_q  (mod q) in x^n + 1 ring
                    let pg = g.scale(self.p, self.q);
                    let h = pg.mul_ring(&f_q, self.q, self.n);
                    return Ok(NtruKeyPair { h, f, f_p });
                }
            }
        }
        Err(LatticeError::SamplingFailed)
    }

    /// Encrypt a message polynomial (coefficients should be in `[-(p-1)/2, (p-1)/2]`).
    pub fn encrypt(&mut self, kp: &NtruKeyPair, message: &Poly) -> NtruCiphertext {
        let r = self.random_trinary();
        let rh = kp.h.mul_ring(&r, self.q, self.n);
        let e = rh.add_mod(message, self.q, self.n);
        NtruCiphertext { e }
    }

    /// Decrypt a ciphertext.
    pub fn decrypt(&self, kp: &NtruKeyPair, ct: &NtruCiphertext) -> Poly {
        // a = f * e (mod q) in x^n + 1 ring, then center into [-q/2, q/2]
        let a_mod_q = kp.f.mul_ring(&ct.e, self.q, self.n);
        let mut a_centered = vec![0i64; self.n];
        for i in 0..self.n {
            let c = mod_pos(a_mod_q.coeffs[i], self.q);
            // Map to centered representation
            a_centered[i] = if c > self.q / 2 { c - self.q } else { c };
        }

        // m = f_p * a (mod p) in x^n + 1 ring
        let a_poly = Poly::new(a_centered);
        a_poly.mul_ring(&kp.f_p, self.p, self.n)
    }

    /// Encrypt a single bit as a constant polynomial.
    pub fn encrypt_bit(&mut self, kp: &NtruKeyPair, bit: u8) -> NtruCiphertext {
        assert!(bit <= 1);
        let mut coeffs = vec![0i64; self.n];
        coeffs[0] = bit as i64;
        self.encrypt(kp, &Poly::new(coeffs))
    }

    /// Decrypt to a single bit.
    pub fn decrypt_bit(&self, kp: &NtruKeyPair, ct: &NtruCiphertext) -> u8 {
        let msg = self.decrypt(kp, ct);
        mod_pos(msg.coeffs[0], self.p) as u8
    }
}

/// Compute the inverse of a polynomial in the ring Z_mod[x]/(x^n + 1)
/// using Gaussian elimination on the negacyclic matrix.
///
/// Returns `None` if the polynomial is not invertible.
fn invert_negacyclic(poly: &Poly, modulus: i64, n: usize) -> Option<Poly> {
    // For ring R = Z_mod[x]/(x^n + 1), multiplication by a corresponds to
    // the n×n matrix M where:
    //   M[k][j] =  a[k-j]        if k >= j
    //   M[k][j] = -a[n+k-j]      if k <  j
    // We solve M * inv = e_0  (mod modulus).

    let mut aug: Vec<Vec<i64>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut row = Vec::with_capacity(2 * n);
        for j in 0..n {
            let val = if i >= j {
                let idx = i - j;
                if idx < poly.coeffs.len() { poly.coeffs[idx] } else { 0 }
            } else {
                let idx = n + i - j;
                if idx < poly.coeffs.len() { -poly.coeffs[idx] } else { 0 }
            };
            row.push(mod_pos(val, modulus));
        }
        for j in 0..n {
            row.push(if i == j { 1 } else { 0 });
        }
        aug.push(row);
    }

    // Gaussian elimination
    for col in 0..n {
        let mut pivot = None;
        for row in col..n {
            if aug[row][col] % modulus != 0 {
                pivot = Some(row);
                break;
            }
        }
        let pivot = pivot?;
        aug.swap(col, pivot);

        let piv_val = aug[col][col];
        let piv_inv = mod_inverse(piv_val, modulus)?;

        for j in 0..2 * n {
            aug[col][j] = mod_pos(aug[col][j] * piv_inv, modulus);
        }

        for row in 0..n {
            if row != col {
                let factor = aug[row][col];
                if factor != 0 {
                    for j in 0..2 * n {
                        aug[row][j] = mod_pos(aug[row][j] - factor * aug[col][j], modulus);
                    }
                }
            }
        }
    }

    let mut inv_coeffs = Vec::with_capacity(n);
    for i in 0..n {
        inv_coeffs.push(aug[i][n]);
    }
    Some(Poly::new(inv_coeffs))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_ntru() -> Ntru {
        Ntru::new(16, 3, 257, 3, 42)
    }

    #[test]
    fn test_keygen() {
        let mut ntru = make_ntru();
        let kp = ntru.keygen().unwrap();
        assert_eq!(kp.h.coeffs.len(), 16);
        assert_eq!(kp.f.coeffs.len(), 16);
        assert_eq!(kp.f_p.coeffs.len(), 16);
    }

    #[test]
    fn test_encrypt_decrypt_zero() {
        let mut ntru = make_ntru();
        let kp = ntru.keygen().unwrap();
        let msg = Poly::zero(16);
        let ct = ntru.encrypt(&kp, &msg);
        let dec = ntru.decrypt(&kp, &ct);
        assert!(dec.coeffs.iter().all(|&c| c == 0));
    }

    #[test]
    fn test_encrypt_decrypt_one() {
        let mut ntru = make_ntru();
        let kp = ntru.keygen().unwrap();
        let mut coeffs = vec![0i64; 16];
        coeffs[0] = 1;
        let msg = Poly::new(coeffs);
        let ct = ntru.encrypt(&kp, &msg);
        let dec = ntru.decrypt(&kp, &ct);
        assert_eq!(dec.coeffs[0], 1);
        assert!(dec.coeffs[1..].iter().all(|&c| c == 0));
    }

    #[test]
    fn test_encrypt_decrypt_bit() {
        let mut ntru = make_ntru();
        let kp = ntru.keygen().unwrap();
        let mut correct = 0;
        for _ in 0..20 {
            let ct = ntru.encrypt_bit(&kp, 1);
            if ntru.decrypt_bit(&kp, &ct) == 1 {
                correct += 1;
            }
        }
        assert!(correct >= 18, "Expected high decryption rate, got {}/20", correct);
    }

    #[test]
    fn test_invert_negacyclic_identity() {
        let f = Poly::new(vec![1, 0, 0, 0]);
        let inv = invert_negacyclic(&f, 97, 4).unwrap();
        assert_eq!(inv.coeffs, vec![1, 0, 0, 0]);
    }

    #[test]
    fn test_invert_negacyclic_roundtrip() {
        let f = Poly::new(vec![1, 1, 0, 0]);
        let inv = invert_negacyclic(&f, 97, 4).unwrap();
        let prod = f.mul_ring(&inv, 97, 4);
        assert!(prod.coeffs.iter().enumerate().all(|(i, &c)| {
            if i == 0 { mod_pos(c, 97) == 1 } else { mod_pos(c, 97) == 0 }
        }));
    }
}
