//! Modular arithmetic and polynomial utilities used throughout the library.

/// Compute `a mod m`, always returning a non-negative result.
///
/// # Examples
/// ```
/// use lattice_crypto_rs::util::mod_pos;
/// assert_eq!(mod_pos(-3, 7), 4);
/// assert_eq!(mod_pos(10, 7), 3);
/// ```
pub fn mod_pos(a: i64, m: i64) -> i64 {
    ((a % m) + m) % m
}

/// Modular exponentiation: computes `base^exp mod modulus` efficiently.
///
/// # Examples
/// ```
/// use lattice_crypto_rs::util::mod_pow;
/// assert_eq!(mod_pow(2, 10, 1000), 1024 % 1000);
/// assert_eq!(mod_pow(3, 0, 7), 1);
/// ```
pub fn mod_pow(base: i64, exp: u64, modulus: i64) -> i64 {
    if modulus == 1 {
        return 0;
    }
    let mut result = 1i64;
    let mut b = mod_pos(base, modulus);
    let mut e = exp;
    while e > 0 {
        if e % 2 == 1 {
            result = mod_pos(result * b, modulus);
        }
        e >>= 1;
        b = mod_pos(b * b, modulus);
    }
    result
}

/// Modular inverse of `a` modulo `m` using the extended Euclidean algorithm.
///
/// Returns `None` if the inverse does not exist (i.e., gcd(a, m) ≠ 1).
///
/// # Examples
/// ```
/// use lattice_crypto_rs::util::mod_inverse;
/// let inv = mod_inverse(3, 7).unwrap();
/// assert_eq!((inv * 3) % 7, 1);
/// ```
pub fn mod_inverse(a: i64, m: i64) -> Option<i64> {
    let (g, x, _) = extended_gcd(mod_pos(a, m), m);
    if g != 1 {
        None
    } else {
        Some(mod_pos(x, m))
    }
}

/// Extended Euclidean algorithm. Returns (gcd, x, y) such that a*x + b*y = gcd.
pub fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if a == 0 {
        (b, 0, 1)
    } else {
        let (g, x, y) = extended_gcd(b % a, a);
        (g, y - (b / a) * x, x)
    }
}

/// Dot product of two vectors modulo `q`.
pub fn dot_mod(a: &[i64], b: &[i64], q: i64) -> i64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum::<i64>() % q
}

/// Add two vectors element-wise, reducing modulo `q`.
pub fn vec_add_mod(a: &[i64], b: &[i64], q: i64) -> Vec<i64> {
    a.iter().zip(b.iter()).map(|(x, y)| mod_pos(x + y, q)).collect()
}

/// A simple polynomial represented as coefficient vector (index = degree).
#[derive(Debug, Clone, PartialEq)]
pub struct Poly {
    /// Coefficients, coeffs[i] is the coefficient of x^i.
    pub coeffs: Vec<i64>,
}

impl Poly {
    /// Create a new polynomial from coefficients.
    pub fn new(coeffs: Vec<i64>) -> Self {
        Self { coeffs }
    }

    /// The zero polynomial of given degree.
    pub fn zero(n: usize) -> Self {
        Self { coeffs: vec![0; n] }
    }

    /// Add two polynomials modulo `q`, trimming to `n` coefficients (cyclotomic ring).
    pub fn add_mod(&self, other: &Poly, q: i64, n: usize) -> Poly {
        let mut result = vec![0i64; n];
        for i in 0..n {
            let a = if i < self.coeffs.len() { self.coeffs[i] } else { 0 };
            let b = if i < other.coeffs.len() { other.coeffs[i] } else { 0 };
            result[i] = mod_pos(a + b, q);
        }
        Poly::new(result)
    }

    /// Multiply two polynomials in the ring Z_q[x]/(x^n + 1).
    pub fn mul_ring(&self, other: &Poly, q: i64, n: usize) -> Poly {
        let mut result = vec![0i64; n];
        for i in 0..n {
            if i >= self.coeffs.len() || self.coeffs[i] == 0 {
                continue;
            }
            for j in 0..n {
                if j >= other.coeffs.len() || other.coeffs[j] == 0 {
                    continue;
                }
                let k = i + j;
                if k < n {
                    result[k] = mod_pos(result[k] + self.coeffs[i] * other.coeffs[j], q);
                } else {
                    // x^(n) = -1 in the cyclotomic ring
                    let red_k = k - n;
                    result[red_k] = mod_pos(result[red_k] - self.coeffs[i] * other.coeffs[j], q);
                }
            }
        }
        Poly::new(result)
    }

    /// Scale all coefficients by a scalar modulo `q`.
    pub fn scale(&self, s: i64, q: i64) -> Poly {
        Poly::new(self.coeffs.iter().map(|c| mod_pos(c * s, q)).collect())
    }
}

/// Deterministic pseudo-random number generator (xorshift64) for reproducible tests.
#[derive(Debug, Clone)]
pub struct XorShift64 {
    state: u64,
}

impl XorShift64 {
    /// Create a new PRNG with the given seed (must not be zero).
    pub fn new(seed: u64) -> Self {
        Self { state: if seed == 0 { 1 } else { seed } }
    }

    /// Generate the next random `u64`.
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    /// Generate a random `i64` in `[0, modulus)`.
    pub fn next_mod(&mut self, modulus: i64) -> i64 {
        (self.next_u64() % (modulus as u64)) as i64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mod_pos() {
        assert_eq!(mod_pos(10, 7), 3);
        assert_eq!(mod_pos(-3, 7), 4);
        assert_eq!(mod_pos(0, 5), 0);
        assert_eq!(mod_pos(-7, 7), 0);
    }

    #[test]
    fn test_mod_pow() {
        assert_eq!(mod_pow(2, 10, 1000), 24);
        assert_eq!(mod_pow(3, 0, 7), 1);
        assert_eq!(mod_pow(5, 3, 13), 125 % 13);
    }

    #[test]
    fn test_mod_inverse() {
        let inv = mod_inverse(3, 7).unwrap();
        assert_eq!((inv * 3) % 7, 1);
        assert!(mod_inverse(2, 4).is_none()); // gcd(2,4) != 1
    }

    #[test]
    fn test_extended_gcd() {
        let (g, x, y) = extended_gcd(240, 46);
        assert_eq!(g, 2);
        assert_eq!(240 * x + 46 * y, 2);
    }

    #[test]
    fn test_dot_mod() {
        let a = vec![1, 2, 3];
        let b = vec![4, 5, 6];
        assert_eq!(dot_mod(&a, &b, 100), (4 + 10 + 18));
    }

    #[test]
    fn test_vec_add_mod() {
        let a = vec![5, 10, 15];
        let b = vec![3, 6, 9];
        let result = vec_add_mod(&a, &b, 20);
        assert_eq!(result, vec![8, 16, 4]);
    }

    #[test]
    fn test_poly_add() {
        let p1 = Poly::new(vec![1, 2, 3]);
        let p2 = Poly::new(vec![4, 5, 6]);
        let sum = p1.add_mod(&p2, 100, 3);
        assert_eq!(sum.coeffs, vec![5, 7, 9]);
    }

    #[test]
    fn test_poly_mul_ring() {
        // (1 + x) * (1 + x) = 1 + 2x + x^2 in ring of degree 4
        let p1 = Poly::new(vec![1, 1, 0, 0]);
        let p2 = Poly::new(vec![1, 1, 0, 0]);
        let prod = p1.mul_ring(&p2, 1000, 4);
        assert_eq!(prod.coeffs, vec![1, 2, 1, 0]);
    }

    #[test]
    fn test_xorshift() {
        let mut rng = XorShift64::new(42);
        let v1 = rng.next_u64();
        let v2 = rng.next_u64();
        assert_ne!(v1, v2);
        assert_ne!(v1, 0);
    }
}
