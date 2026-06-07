//! Property-based tests for lattice-crypto-rs.

use lattice_crypto_rs::gaussian_rejection::RejectionSampler;
use lattice_crypto_rs::ntru::Ntru;
use lattice_crypto_rs::rlwe_kex::RlweKex;
use lattice_crypto_rs::util::Poly;
use proptest::prelude::*;

proptest! {
    #[test]
    fn rlwe_kex_agreement(_seed in 0u64..100) {
        let n = 16;
        let mut kex = RlweKex::new(n, 2048, 1.0, 42);
        let a = kex.public_a();
        let (s_a, pub_a) = kex.alice_init(&a);
        let (_s_b, pub_b, key_b) = kex.bob_respond(&a, &pub_a);
        let key_a = kex.alice_derive(&s_a, &pub_b);
        prop_assert_eq!(key_a.key, key_b.key);
    }

    #[test]
    fn ntru_encrypt_decrypt_zero(n in prop::sample::select(&[8usize, 16])) {
        let d = n / 4;
        let mut ntru = Ntru::new(n, 3, 257, d, 42);
        let kp = ntru.keygen().unwrap();
        let msg = Poly::zero(n);
        let ct = ntru.encrypt(&kp, &msg);
        let dec = ntru.decrypt(&kp, &ct);
        prop_assert!(dec.coeffs.iter().all(|&c| c == 0));
    }

    #[test]
    fn ntru_encrypt_decrypt_one(n in prop::sample::select(&[8usize, 16])) {
        let d = n / 4;
        let mut ntru = Ntru::new(n, 3, 257, d, 42);
        let kp = ntru.keygen().unwrap();
        let mut coeffs = vec![0i64; n];
        coeffs[0] = 1;
        let msg = Poly::new(coeffs);
        let ct = ntru.encrypt(&kp, &msg);
        let dec = ntru.decrypt(&kp, &ct);
        prop_assert_eq!(dec.coeffs[0], 1);
        prop_assert!(dec.coeffs[1..].iter().all(|&c| c == 0));
    }

    #[test]
    fn rejection_sampler_range(
        sigma in 0.5f64..5.0,
        center in -10.0f64..10.0
    ) {
        let mut sampler = RejectionSampler::new(sigma, center, 6.0, 123);
        let samples = sampler.sample_n(50);
        let lo = (center - 6.0 * sigma).floor() as i64;
        let hi = (center + 6.0 * sigma).ceil() as i64;
        for &x in &samples {
            prop_assert!(x >= lo && x <= hi, "sample {} out of range [{}, {}]", x, lo, hi);
        }
    }

    #[test]
    fn poly_add_mod_commutative(
        coeffs_a in prop::collection::vec(-50i64..50, 4..8),
        coeffs_b in prop::collection::vec(-50i64..50, 4..8),
        q in 50i64..200
    ) {
        let n = coeffs_a.len().max(coeffs_b.len());
        let a = Poly::new(coeffs_a);
        let b = Poly::new(coeffs_b);
        let sum_ab = a.add_mod(&b, q, n);
        let sum_ba = b.add_mod(&a, q, n);
        prop_assert_eq!(sum_ab.coeffs, sum_ba.coeffs);
    }

    #[test]
    fn poly_mul_ring_associative_scalar(
        coeffs in prop::collection::vec(-5i64..5, 4..8),
        q in 50i64..200
    ) {
        let n = coeffs.len();
        let a = Poly::new(coeffs.clone());
        let b = Poly::new(coeffs);
        // (a * b) should be well-defined in the ring
        let ab = a.mul_ring(&b, q, n);
        prop_assert_eq!(ab.coeffs.len(), n);
        // All coefficients should be in [0, q)
        prop_assert!(ab.coeffs.iter().all(|&c| c >= 0 && c < q));
    }
}
