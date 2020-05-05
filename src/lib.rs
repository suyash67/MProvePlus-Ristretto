/*

Copyright 2020 by Suyash Bagad, Saravanan Vijayakumaran

This file is part of MProvePlus-Ristretto library
Link: https://github.com/suyash67/MProvePlus-Ristretto

*/

#![allow(non_snake_case)]
#![allow(unused_imports)]

#![cfg_attr(not(feature = "std"), no_std)]
extern crate alloc;

#[macro_use]
extern crate serde_derive;
extern crate serde;

extern crate itertools;
extern crate rand;
extern crate curve25519_dalek;
extern crate sha2;
extern crate sha3;
extern crate bulletproofs;
extern crate merlin;
extern crate digest;
extern crate byteorder;
extern crate libc_print;
extern crate time;


pub mod proofs;
pub mod util;
pub mod errors;
pub mod transcript;
pub mod generators;

#[derive(Copy, PartialEq, Eq, Clone, Debug)]
pub enum Errors {
    InnerProductError,
    MProvePlusError,
    ProofError,
}
