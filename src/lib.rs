#![allow(clippy::needless_range_loop)]

//! # lattice-crypto-rs
//!
//! A pure-Rust library implementing lattice-based cryptography primitives for research and education.
//!
//! ## Modules
//!
//! - [`lwe`] — Learning With Errors (LWE) encryption and decryption
//! - [`ring_lwe`] — Ring-LWE key generation, encryption, and decryption
//! - [`gaussian`] — Discrete Gaussian sampling (Box-Muller, rejection sampling)
//! - [`lattice`] — Lattice basis operations (Gram-Schmidt, LLL reduction, CVP approximation)
//! - [`util`] — Modular arithmetic and polynomial utilities

pub mod gaussian;
pub mod gaussian_rejection;
pub mod lattice;
pub mod lwe;
pub mod ntru;
pub mod ring_lwe;
pub mod rlwe_kex;
pub mod util;
pub mod identity;

pub use gaussian::DiscreteGaussianSampler;
pub use lattice::LatticeBasis;
pub use lwe::LWE;
pub use ring_lwe::RingLWE;
pub use identity::{
    AgentKeyPair, IdentityToken,
    generate_agent_keypair, sign_token, verify_token, derive_shared_digest,
};

/// Common error type for lattice-crypto operations.
#[derive(Debug, Clone, PartialEq)]
pub enum LatticeError {
    /// A dimension mismatch was detected.
    DimensionMismatch { expected: usize, found: usize },
    /// A modulus-related error (e.g., non-positive modulus).
    InvalidModulus(i64),
    /// Sampling failed after maximum retries.
    SamplingFailed,
    /// Matrix is not square where required.
    NotSquare,
}

impl std::fmt::Display for LatticeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LatticeError::DimensionMismatch { expected, found } => {
                write!(f, "dimension mismatch: expected {}, found {}", expected, found)
            }
            LatticeError::InvalidModulus(m) => write!(f, "invalid modulus: {}", m),
            LatticeError::SamplingFailed => write!(f, "sampling failed after max retries"),
            LatticeError::NotSquare => write!(f, "matrix is not square"),
        }
    }
}

impl std::error::Error for LatticeError {}
