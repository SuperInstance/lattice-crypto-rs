# Integration Guide: lattice-crypto

## What This Crate Provides

- **`LWE`** — Learning With Errors encryption: keygen, encrypt, decrypt with dimension n, modulus q, error σ
- **`LWEPublicKey`** / **`LWECiphertext`** — Public key and ciphertext types for LWE
- **`RingLWE`** — Ring-LWE encryption: more efficient than LWE via polynomial ring structure (n must be power of 2)
- **`RingKeyPair`** / **`RingCiphertext`** — Key pair and ciphertext for Ring-LWE
- **`Ntru`** — NTRU-style encryption over Z_q[x]/(x^n - 1) with small/large modulus (p, q)
- **`DiscreteGaussianSampler`** — Discrete Gaussian sampling (Box-Muller, rejection sampling) for error generation
- **`LatticeBasis`** — Lattice basis operations: Gram-Schmidt orthogonalization, LLL reduction, approximate CVP
- **`LatticeError`** — Error type: DimensionMismatch, InvalidModulus, SamplingFailed, NotSquare

This crate provides pure-Rust lattice-based cryptography primitives for research and education: LWE, Ring-LWE, NTRU encryption, discrete Gaussian sampling, and lattice basis reduction (LLL).

## How to Add This Crate

```bash
cargo add lattice-crypto
```

```rust
use lattice_crypto::lwe::LWE;

let mut lwe = LWE::new(4, 97, 2.0, 42);
let secret = lwe.keygen();
let pk = lwe.generate_public_key(&secret, 8);
let ct = lwe.encrypt(&pk.a[0], &pk.b, 1, &secret);
let decrypted = lwe.decrypt(&ct.a, ct.b, &secret);
println!("Decrypted: {}", decrypted);
```

## Integration Points

### agent-identity

- **Why**: agent-identity needs cryptographic identity verification; lattice-crypto provides post-quantum encryption (LWE, Ring-LWE, NTRU) that resists quantum attacks. Agent identity keys should be lattice-based.
- **How**: Use `RingLWE::keygen()` to generate agent identity key pairs. Use the public key for agent identification and the secret key for authentication.

```rust
use lattice_crypto::ring_lwe::RingLWE;

// Generate agent identity keys (post-quantum secure)
let mut rlwe = RingLWE::new(256, 7681, 3.0, 42);
let keypair = rlwe.keygen();

// Public key (a, b) = agent identity
// Secret key = authentication credential
println!("Agent public key a: {:?}", keypair.a);
println!("Agent public key b: {:?}", keypair.b);
```

### agent-handshake

- **Why**: agent-handshake needs a key exchange protocol; lattice-crypto provides `rlwe_kex` (Ring-LWE key exchange) — a post-quantum Diffie-Hellman analogue.
- **How**: Use the RLWE-KEX module for mutual key exchange between agents. Each agent generates a key pair, shares the public component, and derives a shared secret.

```rust
use lattice_crypto::ring_lwe::RingLWE;
use lattice_crypto::lwe::LWE;

// Agent A generates keys
let mut alice = RingLWE::new(512, 12289, 3.19, 0xDEAD);
let kp_a = alice.keygen();

// Agent B generates keys
let mut bob = RingLWE::new(512, 12289, 3.19, 0xBEEF);
let kp_b = bob.keygen();

// Exchange public keys → derive shared secret via RLWE-KEX
// (Full KEX protocol in rlwe_kex module)
```

## For AI Agents

- **Context needed**: Security parameter (dimension n), modulus q, error standard deviation σ, seed for deterministic key generation
- **Key imports**: `lattice_crypto::{LWE, RingLWE, LatticeBasis, LatticeError}`
- **Integration pattern**: Instantiate crypto system → `keygen()` → distribute public key → `encrypt()` / `decrypt()` → verify with `LatticeBasis` operations
- **Error handling**: `LatticeError::DimensionMismatch` (wrong vector sizes), `LatticeError::InvalidModulus` (non-positive q), `LatticeError::SamplingFailed` (Gaussian sampler exhaustion). Always handle `SamplingFailed` — retry with different seed.

## For Humans

- **Prerequisites**: Basic lattice cryptography (LWE problem), polynomial rings, modular arithmetic
- **Learning path**: Start with `util.rs` (modular arithmetic helpers), then `gaussian.rs` (error sampling), then `lwe.rs` (simplest encryption), then `ring_lwe.rs` (efficient), then `ntru.rs` (classic NTRU)
- **Common pitfalls**:
  - Ring-LWE dimension n MUST be a power of 2 — the crate panics otherwise
  - Error σ too small → decryption failures (noise overwhelms message); too large → security weakens
  - NTRU requires f to be invertible mod q AND mod p — not all random polynomials work; the crate handles this internally
  - Seeds are for deterministic testing — use true randomness in production
  - This is research/educational code — not audited for production cryptographic use
