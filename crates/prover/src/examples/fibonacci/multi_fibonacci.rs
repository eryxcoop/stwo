
use std::iter::zip;

use super::air::MultiFibonacciAir;
use super::Fibonacci;

use crate::core::backend::cpu::CpuCircleEvaluation;
use crate::core::channel::{Blake2sChannel, Channel};
use crate::core::fields::m31::BaseField;
use crate::core::fields::IntoSlice;
use crate::core::poly::BitReversedOrder;
use crate::core::prover::{prove, verify, ProvingError, StarkProof, VerificationError};
use crate::core::vcs::blake2_hash::Blake2sHasher;
use crate::core::vcs::hasher::Hasher;

pub struct MultiFibonacci {
    pub air: MultiFibonacciAir,
    log_sizes: Vec<u32>,
    claims: Vec<BaseField>,
}

impl MultiFibonacci {
    pub fn new(log_sizes: Vec<u32>, claims: Vec<BaseField>) -> Self {
        assert!(!log_sizes.is_empty());
        assert_eq!(log_sizes.len(), claims.len());
        let air = MultiFibonacciAir::new(&log_sizes, &claims);
        Self {
            air,
            log_sizes,
            claims,
        }
    }

    pub fn get_trace(&self) -> Vec<CpuCircleEvaluation<BaseField, BitReversedOrder>> {
        zip(&self.log_sizes, &self.claims)
            .map(|(log_size, claim)| {
                let fib = Fibonacci::new(*log_size, *claim);
                fib.get_trace()
            })
            .collect()
    }

    pub fn prove(&self) -> Result<StarkProof, ProvingError> {
        let channel =
            &mut Blake2sChannel::new(Blake2sHasher::hash(BaseField::into_slice(&self.claims)));
        prove(&self.air, channel, self.get_trace())
    }

    pub fn verify(&self, proof: StarkProof) -> Result<(), VerificationError> {
        let channel =
            &mut Blake2sChannel::new(Blake2sHasher::hash(BaseField::into_slice(&self.claims)));
        verify(proof, &self.air, channel)
    }
}

#[cfg(test)]
mod tests {
    use super::MultiFibonacci;
    use crate::m31;

    #[test]
    fn test_rectangular_multi_fibonacci() {
        let multi_fib = MultiFibonacci::new(vec![5; 16], vec![m31!(443693538); 16]);
        let proof = multi_fib.prove().unwrap();
        multi_fib.verify(proof).unwrap();
    }

    #[test]
    fn test_mixed_degree_multi_fibonacci() {
        let multi_fib = MultiFibonacci::new(
            // TODO(spapini): Change order of log_sizes.
            vec![3, 5, 7],
            vec![m31!(1056169651), m31!(443693538), m31!(722122436)],
        );
        let proof = multi_fib.prove().unwrap();
        multi_fib.verify(proof).unwrap();
    }
}
