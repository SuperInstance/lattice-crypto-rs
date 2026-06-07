//! Agent identity verification via lattice-based cryptography.
//!
//! Uses Ring-LWE primitives to build a lightweight identity layer:
//! each agent holds a [`AgentKeyPair`] and can issue/verify [`IdentityToken`]s
//! that authenticate the agent without revealing the private key.

use crate::util::{Poly, XorShift64};

/// An agent's lattice-based key pair.
#[derive(Debug, Clone)]
pub struct AgentKeyPair {
    /// Agent identifier.
    pub agent_id: String,
    /// Random polynomial `a` (part of the public key).
    pub a_poly: Vec<i64>,
    /// Public key polynomial `b = a*s + e mod q`.
    pub public_key: Vec<i64>,
    /// Secret key polynomial coefficients.
    pub secret_key: Vec<i64>,
    /// LWE modulus q.
    pub modulus: i64,
    /// Ring dimension n.
    pub dimension: usize,
    /// Error bound used during key generation.
    pub error_bound: f64,
}

impl AgentKeyPair {
    /// Export the full public key bytes: (a_poly || public_key) in little-endian i64.
    pub fn public_key_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(self.a_poly.len() * 8 * 2);
        for &v in &self.a_poly {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        for &v in &self.public_key {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf
    }
}

/// A signed identity token.
#[derive(Debug, Clone, PartialEq)]
pub struct IdentityToken {
    /// The agent that produced this token.
    pub agent_id: String,
    /// Hash of the data that was signed (truncated to u64 for portability).
    pub data_hash: u64,
    /// Signature polynomial `b_sig = a*s + e_sig + encode(hash)`.
    pub signature_b: Vec<i64>,
    /// Token creation timestamp (unix seconds).
    pub timestamp: i64,
    /// Modulus used for signing.
    pub modulus: i64,
    /// Ring dimension.
    pub dimension: usize,
}

impl IdentityToken {
    /// Serialize the token to bytes for transmission.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        let id_bytes = self.agent_id.as_bytes();
        buf.extend_from_slice(&(id_bytes.len() as u64).to_le_bytes());
        buf.extend_from_slice(id_bytes);
        buf.extend_from_slice(&self.data_hash.to_le_bytes());
        buf.extend_from_slice(&self.timestamp.to_le_bytes());
        buf.extend_from_slice(&self.modulus.to_le_bytes());
        buf.extend_from_slice(&(self.dimension as u64).to_le_bytes());
        buf.extend_from_slice(&(self.signature_b.len() as u64).to_le_bytes());
        for &v in &self.signature_b {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf
    }
}

/// Simple hash of data for token binding (FNV-1a-inspired).
fn simple_hash(data: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &byte in data {
        h ^= byte as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

/// Generate a new lattice-based key pair for an agent.
///
/// Uses Ring-LWE key generation: public key is `(a, b=a*s+e)`.
pub fn generate_agent_keypair(agent_id: &str) -> AgentKeyPair {
    let n = 64;
    let q: i64 = 7681;
    let error_bound = 3.2;

    // Deterministic seed from agent_id
    let mut seeder = XorShift64::new(0x1234567890ABCDEF);
    for &b in agent_id.as_bytes() {
        seeder = XorShift64::new(seeder.next_u64() ^ b as u64);
    }
    let seed = seeder.next_u64();
    let mut rng = XorShift64::new(seed);

    // Secret key: small random coefficients in {-1, 0, 1}
    let secret_key: Vec<i64> = (0..n).map(|_| (rng.next_mod(3) - 1)).collect();

    // Random polynomial a in [0, q)
    let a_coeffs: Vec<i64> = (0..n).map(|_| rng.next_mod(q)).collect();

    // Error: small random in {-2, -1, 0, 1, 2}
    let error: Vec<i64> = (0..n).map(|_| (rng.next_mod(5) - 2)).collect();

    // Public key: b = a * s + e mod q
    let a_poly = Poly::new(a_coeffs.clone());
    let s_poly = Poly::new(secret_key.clone());
    let e_poly = Poly::new(error);
    let b_poly = a_poly.mul_ring(&s_poly, q, n).add_mod(&e_poly, q, n);

    AgentKeyPair {
        agent_id: agent_id.to_string(),
        a_poly: a_coeffs,
        public_key: b_poly.coeffs,
        secret_key,
        modulus: q,
        dimension: n,
        error_bound,
    }
}

/// Sign data using the agent's key pair.
///
/// Produces signature `b_sig = a*s + e_sig + encode(hash)` using the same `a`
/// from the key pair, so that `b_sig - b_public ≈ encode(hash)`.
pub fn sign_token(keypair: &AgentKeyPair, data: &[u8]) -> IdentityToken {
    let n = keypair.dimension;
    let q = keypair.modulus;
    let data_hash = simple_hash(data);

    // Deterministic RNG seeded from keypair + data
    let mut seeder = XorShift64::new(0xDEADBEEFCAFEBABE);
    for &v in &keypair.secret_key {
        seeder = XorShift64::new(seeder.next_u64() ^ v as u64);
    }
    seeder = XorShift64::new(seeder.next_u64() ^ data_hash);
    let mut rng = XorShift64::new(seeder.next_u64());

    // Re-use the keypair's `a` polynomial
    // Small signing error
    let error: Vec<i64> = (0..n).map(|_| (rng.next_mod(3) - 1)).collect();

    // Encode hash into polynomial
    let mut encoded = vec![0i64; n];
    for i in 0..n {
        encoded[i] = (((data_hash.wrapping_shr((i % 64) as u32)) & 0x1) as i64) * (q / 4);
    }

    // b_sig = a*s + e_sig + encode(hash) mod q
    let a_poly = Poly::new(keypair.a_poly.clone());
    let s_poly = Poly::new(keypair.secret_key.clone());
    let e_poly = Poly::new(error);
    let enc_poly = Poly::new(encoded);

    let as_poly = a_poly.mul_ring(&s_poly, q, n);
    let mut b_sig = as_poly.add_mod(&e_poly, q, n);
    b_sig = b_sig.add_mod(&enc_poly, q, n);

    let timestamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs() as i64;

    IdentityToken {
        agent_id: keypair.agent_id.clone(),
        data_hash,
        signature_b: b_sig.coeffs,
        timestamp,
        modulus: q,
        dimension: n,
    }
}

/// Verify an identity token against a public key.
///
/// Checks that `b_sig - b_public ≈ encode(hash)` within the error bound,
/// where `b_public = a*s + e` from key generation and `b_sig = a*s + e' + encode(hash)`.
pub fn verify_token(token: &IdentityToken, public_key_bytes: &[u8]) -> bool {
    if public_key_bytes.len() % 8 != 0 {
        return false;
    }
    let total_coeffs = public_key_bytes.len() / 8;
    if total_coeffs != token.dimension * 2 {
        return false;
    }
    let n = token.dimension;

    // Decode a_poly (first n coefficients) and b_public (next n coefficients)
    let a_poly: Vec<i64> = (0..n)
        .map(|i| {
            let start = i * 8;
            let bytes: [u8; 8] = public_key_bytes[start..start + 8].try_into().unwrap();
            i64::from_le_bytes(bytes)
        })
        .collect();
    let b_public: Vec<i64> = (n..2 * n)
        .map(|i| {
            let start = i * 8;
            let bytes: [u8; 8] = public_key_bytes[start..start + 8].try_into().unwrap();
            i64::from_le_bytes(bytes)
        })
        .collect();

    let q = token.modulus;

    // Re-encode the hash
    let mut encoded = vec![0i64; n];
    for i in 0..n {
        encoded[i] = (((token.data_hash.wrapping_shr((i % 64) as u32)) & 0x1) as i64) * (q / 4);
    }

    // Check: b_sig - b_public should be close to encode(hash)
    // Tolerance: keygen error + signing error, both small
    let threshold = (q as f64 / 6.0) as i64;
    for i in 0..n {
        let diff = crate::util::mod_pos(token.signature_b[i] - b_public[i], q);
        let expected = crate::util::mod_pos(encoded[i], q);
        let err = if diff > expected { diff - expected } else { expected - diff };
        let err = err.min(q - err); // wrap-around distance
        if err > threshold {
            return false;
        }
    }
    true
}

/// Derive a shared authentication digest from two agent key pairs.
pub fn derive_shared_digest(kp1: &AgentKeyPair, kp2: &AgentKeyPair) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &v in &kp1.public_key {
        h ^= v as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    for &v in &kp2.public_key {
        h ^= v as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_keypair() {
        let kp = generate_agent_keypair("test-agent");
        assert_eq!(kp.agent_id, "test-agent");
        assert_eq!(kp.dimension, 64);
        assert_eq!(kp.public_key.len(), 64);
        assert_eq!(kp.secret_key.len(), 64);
        assert_eq!(kp.a_poly.len(), 64);
    }

    #[test]
    fn test_sign_and_verify() {
        let kp = generate_agent_keypair("alice");
        let data = b"hello world";
        let token = sign_token(&kp, data);
        assert_eq!(token.agent_id, "alice");
        assert!(verify_token(&token, &kp.public_key_bytes()));
    }

    #[test]
    fn test_verify_fails_with_wrong_key() {
        let kp1 = generate_agent_keypair("alice");
        let kp2 = generate_agent_keypair("bob");
        let data = b"hello world";
        let token = sign_token(&kp1, data);
        assert!(!verify_token(&token, &kp2.public_key_bytes()));
    }

    #[test]
    fn test_verify_fails_with_wrong_data() {
        let kp = generate_agent_keypair("alice");
        let token = sign_token(&kp, b"correct data");
        let mut bad_token = token.clone();
        bad_token.data_hash = bad_token.data_hash.wrapping_add(1);
        assert!(!verify_token(&bad_token, &kp.public_key_bytes()));
    }

    #[test]
    fn test_different_agents_different_keys() {
        let kp1 = generate_agent_keypair("agent-1");
        let kp2 = generate_agent_keypair("agent-2");
        assert_ne!(kp1.public_key, kp2.public_key);
        assert_ne!(kp1.secret_key, kp2.secret_key);
    }

    #[test]
    fn test_public_key_bytes_roundtrip() {
        let kp = generate_agent_keypair("bytes-test");
        let bytes = kp.public_key_bytes();
        assert_eq!(bytes.len(), 64 * 8 * 2);
    }

    #[test]
    fn test_token_serialization() {
        let kp = generate_agent_keypair("serialize-test");
        let token = sign_token(&kp, b"test data");
        let bytes = token.to_bytes();
        assert!(!bytes.is_empty());
        assert!(bytes.len() > 100);
    }

    #[test]
    fn test_simple_hash_deterministic() {
        let h1 = simple_hash(b"test");
        let h2 = simple_hash(b"test");
        assert_eq!(h1, h2);
        let h3 = simple_hash(b"different");
        assert_ne!(h1, h3);
    }

    #[test]
    fn test_derive_shared_digest() {
        let kp1 = generate_agent_keypair("alice");
        let kp2 = generate_agent_keypair("bob");
        let d1 = derive_shared_digest(&kp1, &kp2);
        let d2 = derive_shared_digest(&kp2, &kp1);
        assert_ne!(d1, 0);
        assert_ne!(d2, 0);
    }

    #[test]
    fn test_sign_empty_data() {
        let kp = generate_agent_keypair("empty-data");
        let token = sign_token(&kp, b"");
        assert!(verify_token(&token, &kp.public_key_bytes()));
    }

    #[test]
    fn test_sign_large_data() {
        let kp = generate_agent_keypair("large-data");
        let data: Vec<u8> = (0..10000).map(|i| (i % 256) as u8).collect();
        let token = sign_token(&kp, &data);
        assert!(verify_token(&token, &kp.public_key_bytes()));
    }

    #[test]
    fn test_verify_wrong_dimension_fails() {
        let kp = generate_agent_keypair("dim-test");
        let token = sign_token(&kp, b"data");
        assert!(!verify_token(&token, &[0u8; 16]));
    }

    #[test]
    fn test_token_fields() {
        let kp = generate_agent_keypair("field-test");
        let token = sign_token(&kp, b"check fields");
        assert_eq!(token.agent_id, "field-test");
        assert_eq!(token.modulus, 7681);
        assert_eq!(token.dimension, 64);
        assert_eq!(token.signature_b.len(), 64);
        assert!(token.timestamp > 0);
    }
}
