//! Lattice basis operations: Gram-Schmidt orthogonalization, LLL reduction,
//! and approximate closest vector problem (CVP) solving.


/// A lattice basis represented as a matrix of row vectors.
#[derive(Debug, Clone, PartialEq)]
pub struct LatticeBasis {
    /// Row vectors of the basis matrix.
    pub rows: Vec<Vec<f64>>,
}

impl LatticeBasis {
    /// Create a new lattice basis from row vectors.
    pub fn new(rows: Vec<Vec<f64>>) -> Self {
        Self { rows }
    }

    /// Number of basis vectors (rank).
    pub fn rank(&self) -> usize {
        self.rows.len()
    }

    /// Dimension of the space.
    pub fn dim(&self) -> usize {
        self.rows.first().map_or(0, |r| r.len())
    }

    /// Compute the Gram-Schmidt orthogonalization.
    ///
    /// Returns (orthogonalized vectors, mu coefficients) where
    /// `mu[i][j]` is the projection coefficient of `b_i` onto `b*_j`.
    pub fn gram_schmidt(&self) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
        let n = self.rank();
        let d = self.dim();
        let mut gs = vec![vec![0.0; d]; n];
        let mut mu = vec![vec![0.0; n]; n];

        for i in 0..n {
            gs[i] = self.rows[i].clone();
            for j in 0..i {
                mu[i][j] = dot(&gs[j], &self.rows[i]) / dot(&gs[j], &gs[j]);
                for k in 0..d {
                    gs[i][k] -= mu[i][j] * gs[j][k];
                }
            }
        }

        (gs, mu)
    }

    /// Compute the determinant (volume) of the lattice.
    ///
    /// This is the product of the Gram-Schmidt vector norms.
    pub fn determinant(&self) -> f64 {
        let (gs, _) = self.gram_schmidt();
        gs.iter().map(|v| dot(v, v).sqrt()).product()
    }

    /// Perform LLL basis reduction with the given delta parameter (typically 0.75).
    ///
    /// Returns a new, more orthogonal basis for the same lattice.
    pub fn lll_reduce(&self, delta: f64) -> LatticeBasis {
        let mut b = self.rows.clone();
        let n = b.len();
        let d = if n > 0 { b[0].len() } else { 0 };
        if n == 0 || d == 0 {
            return LatticeBasis::new(b);
        }

        let mut k = 1;
        while k < n {
            // Compute GS on the fly
            let (_gs, mu) = gram_schmidt_inline(&b);

            // Size-reduce b[k]
            for j in (0..k).rev() {
                if mu[k][j].abs() > 0.5 {
                    let round_mu = mu[k][j].round();
                    for i in 0..d {
                        b[k][i] -= round_mu * b[j][i];
                    }
                }
            }

            // Recompute after size reduction
            let (gs, mu) = gram_schmidt_inline(&b);

            // Lovász condition
            let gs_k_norm_sq = dot(&gs[k], &gs[k]);
            let gs_k1_norm_sq = if k > 0 { dot(&gs[k - 1], &gs[k - 1]) } else { 0.0 };

            if gs_k_norm_sq >= (delta - mu[k][k - 1].powi(2)) * gs_k1_norm_sq {
                k += 1;
            } else {
                // Swap b[k] and b[k-1]
                b.swap(k, k - 1);
                k = if k > 1 { k - 1 } else { 1 };
            }
        }

        LatticeBasis::new(b)
    }

    /// Approximate CVP: find the closest lattice vector to the target.
    ///
    /// Uses Babai's nearest plane algorithm on an LLL-reduced basis.
    pub fn closest_vector(&self, target: &[f64]) -> Vec<f64> {
        let reduced = self.lll_reduce(0.75);
        let (gs, _mu) = reduced.gram_schmidt();
        let n = reduced.rank();
        let d = reduced.dim();
        let mut b = target.to_vec();

        for i in (0..n).rev() {
            let gs_norm_sq = dot(&gs[i], &gs[i]);
            if gs_norm_sq < 1e-12 {
                continue;
            }
            let coeff = dot(&gs[i], &b) / gs_norm_sq;
            let rounded = coeff.round();
            for j in 0..d {
                b[j] -= rounded * reduced.rows[i][j];
            }
        }

        // Closest vector = target - residual
        target.iter().zip(b.iter()).map(|(t, r)| t - r).collect()
    }

    /// Check if the basis vectors are linearly independent (non-zero GS norms).
    pub fn is_independent(&self) -> bool {
        let (gs, _) = self.gram_schmidt();
        gs.iter().all(|v| dot(v, v) > 1e-12)
    }

    /// Compute the Hadamard ratio (measure of orthogonality, 1.0 = orthogonal).
    pub fn hadamard_ratio(&self) -> f64 {
        let det = self.determinant();
        if det.abs() < 1e-15 {
            return 0.0;
        }
        let product_norms: f64 = self.rows.iter().map(|v| dot(v, v).sqrt()).product();
        if product_norms.abs() < 1e-15 {
            return 0.0;
        }
        (det.abs() / product_norms).powf(1.0 / self.rank() as f64)
    }
}

/// Dot product of two f64 vectors.
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Inline Gram-Schmidt for a vector of vectors.
fn gram_schmidt_inline(b: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = b.len();
    let d = if n > 0 { b[0].len() } else { 0 };
    let mut gs = vec![vec![0.0; d]; n];
    let mut mu = vec![vec![0.0; n]; n];

    for i in 0..n {
        gs[i] = b[i].clone();
        for j in 0..i {
            let gs_j_norm = dot(&gs[j], &gs[j]);
            if gs_j_norm.abs() < 1e-15 {
                mu[i][j] = 0.0;
                continue;
            }
            mu[i][j] = dot(&gs[j], &b[i]) / gs_j_norm;
            for k in 0..d {
                gs[i][k] -= mu[i][j] * gs[j][k];
            }
        }
    }

    (gs, mu)
}

/// Compute the GCD of two integers using the Euclidean algorithm.
pub fn gcd(a: i64, b: i64) -> i64 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// Compute the LCM of two integers.
pub fn lcm(a: i64, b: i64) -> i64 {
    if a == 0 || b == 0 {
        0
    } else {
        (a.abs() / gcd(a, b)) * b.abs()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lattice_creation() {
        let basis = LatticeBasis::new(vec![vec![1.0, 0.0], vec![0.0, 1.0]]);
        assert_eq!(basis.rank(), 2);
        assert_eq!(basis.dim(), 2);
    }

    #[test]
    fn test_gram_schmidt_identity() {
        let basis = LatticeBasis::new(vec![vec![1.0, 0.0], vec![0.0, 1.0]]);
        let (gs, mu) = basis.gram_schmidt();
        assert!((gs[0][0] - 1.0).abs() < 1e-10);
        assert!((gs[1][1] - 1.0).abs() < 1e-10);
        assert!((mu[1][0]).abs() < 1e-10);
    }

    #[test]
    fn test_gram_schmidt_nontrivial() {
        let basis = LatticeBasis::new(vec![vec![1.0, 1.0], vec![1.0, 0.0]]);
        let (gs, mu) = basis.gram_schmidt();
        // gs[0] = [1,1], gs[1] = [1,0] - mu[1][0]*[1,1]
        // mu[1][0] = dot([1,1],[1,0])/dot([1,1],[1,1]) = 1/2
        assert!((mu[1][0] - 0.5).abs() < 1e-10);
        // gs[1] = [1-0.5, 0-0.5] = [0.5, -0.5]
        assert!((gs[1][0] - 0.5).abs() < 1e-10);
        assert!((gs[1][1] - (-0.5)).abs() < 1e-10);
    }

    #[test]
    fn test_determinant() {
        let basis = LatticeBasis::new(vec![vec![3.0, 0.0], vec![0.0, 4.0]]);
        let det = basis.determinant();
        assert!((det - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_lll_simple() {
        let basis = LatticeBasis::new(vec![vec![1.0, 1.0], vec![1.0, 0.0]]);
        let reduced = basis.lll_reduce(0.75);
        assert_eq!(reduced.rank(), 2);
        // First vector should be shorter after reduction
        let orig_norm0 = dot(&basis.rows[0], &basis.rows[0]).sqrt();
        let red_norm0 = dot(&reduced.rows[0], &reduced.rows[0]).sqrt();
        assert!(red_norm0 <= orig_norm0 + 1e-10);
    }

    #[test]
    fn test_lll_preserves_lattice() {
        let basis = LatticeBasis::new(vec![vec![4.0, 1.0], vec![2.0, 3.0]]);
        let reduced = basis.lll_reduce(0.75);
        // Determinant should be preserved
        assert!((basis.determinant().abs() - reduced.determinant().abs()).abs() < 1e-6);
    }

    #[test]
    fn test_is_independent() {
        let basis = LatticeBasis::new(vec![vec![1.0, 0.0], vec![0.0, 1.0]]);
        assert!(basis.is_independent());
        let dependent = LatticeBasis::new(vec![vec![1.0, 0.0], vec![2.0, 0.0]]);
        assert!(!dependent.is_independent());
    }

    #[test]
    fn test_hadamard_ratio_orthogonal() {
        let basis = LatticeBasis::new(vec![vec![1.0, 0.0], vec![0.0, 1.0]]);
        let ratio = basis.hadamard_ratio();
        assert!((ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_hadamard_ratio_non_orthogonal() {
        let basis = LatticeBasis::new(vec![vec![10.0, 0.0], vec![10.0, 1.0]]);
        let ratio = basis.hadamard_ratio();
        assert!(ratio < 1.0 && ratio > 0.0);
    }

    #[test]
    fn test_closest_vector() {
        let basis = LatticeBasis::new(vec![vec![1.0, 0.0], vec![0.0, 1.0]]);
        let target = vec![0.3, 0.7];
        let closest = basis.closest_vector(&target);
        assert!((closest[0] - 0.0).abs() < 1e-10);
        assert!((closest[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_gcd() {
        assert_eq!(gcd(12, 8), 4);
        assert_eq!(gcd(7, 13), 1);
        assert_eq!(gcd(0, 5), 5);
    }

    #[test]
    fn test_lcm() {
        assert_eq!(lcm(4, 6), 12);
        assert_eq!(lcm(7, 13), 91);
    }
}
