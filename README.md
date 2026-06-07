# lattice-crypto-rs

Lattice-based cryptography primitives: LWE encryption/decryption, Ring-LWE structures, Gaussian sampling, and lattice basis operations.

## Features

- **LWE**: Learning With Errors encryption and decryption
- **Ring-LWE**: Ring-LWE key generation, encryption, and decryption  
- **Gaussian Sampling**: Box-Muller and rejection sampling from discrete Gaussian
- **Lattice Basis**: Gram-Schmidt orthogonalization, LLL reduction, CVP approximation
- **Utilities**: Modular arithmetic, polynomial ring operations

Pure Rust, no external dependencies.

## Usage

```rust
use lattice_crypto_rs::LWE;

let mut lwe = LWE::new(4, 97, 2.0, 42);
let secret = lwe.keygen();
let pk = lwe.public_keygen(&secret, 100);
let ct = lwe.encrypt(&pk, 1);
assert_eq!(lwe.decrypt(&secret, &ct), 1);
```

License: MIT OR Apache-2.0
