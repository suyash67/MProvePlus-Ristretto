/*

Copyright (c) 2018 Chain, Inc.
Some parts of this file are from the bulletproofs package by dalek-cryptography.
Link: https://github.com/dalek-cryptography/bulletproofs
Description: Utilities for Scalar operations.

*/

#![allow(non_snake_case)]

use curve25519_dalek::scalar::Scalar;
use alloc::vec;
use alloc::vec::Vec;
use libc_print::{libc_println};

/// Provides an iterator over the powers of a `Scalar`.
///
/// This struct is created by the `exp_iter` function.
pub struct ScalarExp {
    x: Scalar,
    next_exp_x: Scalar,
}

impl Iterator for ScalarExp {
    type Item = Scalar;

    fn next(&mut self) -> Option<Scalar> {
        let exp_x = self.next_exp_x;
        self.next_exp_x *= self.x;
        Some(exp_x)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (usize::max_value(), None)
    }
}

/// Return an iterator of the powers of `x`.
pub fn exp_iter(x: Scalar) -> ScalarExp {
    let next_exp_x = Scalar::one();
    ScalarExp { x, next_exp_x }
}

pub fn add_vec(a: &[Scalar], b: &[Scalar]) -> Vec<Scalar> {
    if a.len() != b.len() {
        // throw some error
        //println!("lengths of vectors don't match for vector addition");
    }
    let mut out = vec![Scalar::zero(); b.len()];
    for i in 0..a.len() {
        out[i] = a[i] + b[i];
    }
    out
}

/// Raises `x` to the power `n` using binary exponentiation,
/// with (1 to 2)*lg(n) scalar multiplications.
/// TODO: a consttime version of this would be awfully similar to a Montgomery ladder.
pub fn scalar_exp_vartime(x: &Scalar, mut n: u64) -> Scalar {
    let mut result = Scalar::one();
    let mut aux = *x; // x, x^2, x^4, x^8, ...
    while n > 0 {
        let bit = n & 1;
        if bit == 1 {
            result = result * aux;
        }
        n = n >> 1;
        aux = aux * aux; // FIXME: one unnecessary mult at the last step here!
    }
    result
}

// takes the sum of all of the powers of x, up to n
fn sum_of_powers_slow(x: &Scalar, n: usize) -> Scalar {
    exp_iter(*x).take(n).sum()
}

/// Takes the sum of all the powers of `x`, up to `n`
/// If `n` is a power of 2, it uses the efficient algorithm with `2*lg n` multiplications and additions.
/// If `n` is not a power of 2, it uses the slow algorithm with `n` multiplications and additions.
/// In the Bulletproofs case, all calls to `sum_of_powers` should have `n` as a power of 2.
pub fn sum_of_powers(x: &Scalar, n: usize) -> Scalar {
    if !n.is_power_of_two() {
        return sum_of_powers_slow(x, n);
    }
    if n == 0 || n == 1 {
        return Scalar::from(n as u64);
    }
    let mut m = n;
    let mut result = Scalar::one() + x;
    let mut factor = *x;
    while m > 2 {
        factor = factor * factor;
        result = result + factor * result;
        m = m / 2;
    }
    result
}

/// Batch inversion using Montgomery's trick. 
/// Note that this would throw an error if any of the inputs is zero.
pub fn batch_invert(inputs: &mut Vec<Scalar>) {

    let n: usize = inputs.len();
    let mut temp_products: Vec<Scalar> = Vec::with_capacity(n);
    temp_products.push(inputs[0]);
    for i in 1..n {
        assert!(inputs[i] != Scalar::zero());
        temp_products.push(temp_products[i - 1] * inputs[i]);
    }

    let inverted: Scalar = temp_products[n - 1].invert();
    let mut remaining = Scalar::one();
    for i in (1..n).rev() {
        let temp = inputs[i];
        inputs[i] = temp_products[i - 1] * inverted * remaining;
        remaining *= temp;
    }
    inputs[0] = remaining * inverted;
}

/// Given `data` with `len >= 32`, return the first 32 bytes.
pub fn read32(data: &[u8]) -> [u8; 32] {
    let mut buf32 = [0u8; 32];
    buf32[..].copy_from_slice(&data[..32]);
    buf32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn batch_invert_test() {
        let n: usize = 512;
        let mut rng = rand::thread_rng();
        let mut inputs_copy: Vec<Scalar> = Vec::new();
        let mut inputs = (0..n).map(|_| {
                let input = Scalar::random(&mut rng);
                inputs_copy.push(input);
                input.clone()
            }).collect::<Vec<Scalar>>();
        
        batch_invert(&mut inputs);
        for i in 0..n {
            assert_eq!(inputs_copy[i].invert(), inputs[i]);
        }
    }

}