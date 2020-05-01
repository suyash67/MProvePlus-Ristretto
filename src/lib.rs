#![allow(non_snake_case)]
/*

Copyright 2018 by Suyash Bagad, Saravanan Vijayakumaran

This file is part of revelioPlus library
(<add a link to github>)

*/

// based on the paper: <link to paper>
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
