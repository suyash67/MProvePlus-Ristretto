//! Errors related to proving and verifying proofs.

extern crate alloc;

#[cfg(feature = "std")]
use thiserror::Error;

/// Represents an error in proof creation, verification, or parsing.
#[derive(Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "std", derive(Error))]
pub enum ProofError {
    /// This error occurs when a proof failed to verify.
    #[cfg_attr(feature = "std", error("Proof verification failed."))]
    VerificationError,
    /// This error occurs when the proof encoding is malformed.
    #[cfg_attr(feature = "std", error("Proof data could not be parsed."))]
    FormatError,
    /// This error occurs during proving if the number of blinding
    /// factors does not match the number of values.
    #[cfg_attr(feature = "std", error("Wrong number of blinding factors supplied."))]
    WrongNumBlindingFactors,
    /// This error occurs when attempting to create a proof with
    /// bitsize other than \\(8\\), \\(16\\), \\(32\\), or \\(64\\).
    #[cfg_attr(feature = "std", error("Invalid bitsize, must have n = 8,16,32,64."))]
    InvalidBitsize,
    /// This error occurs when attempting to create an aggregated
    /// proof with non-power-of-two aggregation size.
    #[cfg_attr(
        feature = "std",
        error("Invalid aggregation size, m must be a power of 2.")
    )]
    InvalidAggregation,
    /// This error occurs when there are insufficient generators for the proof.
    #[cfg_attr(
        feature = "std",
        error("Invalid generators size, too few generators for proof")
    )]
    InvalidGeneratorsLength,
}