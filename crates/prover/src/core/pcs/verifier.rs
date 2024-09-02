use std::iter::zip;

use itertools::Itertools;

use super::super::circle::CirclePoint;
use super::super::fields::qm31::SecureField;
use super::super::fri::{CirclePolyDegreeBound, FriVerifier};
use super::quotients::{fri_answers, PointSample};
use super::utils::TreeVec;
use super::{CommitmentSchemeProof, PcsConfig};
use crate::core::channel::{Channel, MerkleChannel};
use crate::core::prover::VerificationError;
use crate::core::vcs::ops::MerkleHasher;
use crate::core::vcs::verifier::MerkleVerifier;
use crate::core::ColumnVec;

/// The verifier side of a FRI polynomial commitment scheme. See [super].
#[derive(Default)]
pub struct CommitmentSchemeVerifier<MC: MerkleChannel> {
    pub trees: TreeVec<MerkleVerifier<MC::H>>,
    pub config: PcsConfig,
}

impl<MC: MerkleChannel> CommitmentSchemeVerifier<MC> {
    pub fn new(config: PcsConfig) -> Self {
        Self {
            trees: TreeVec::default(),
            config,
        }
    }

    /// A [TreeVec<ColumnVec>] of the log sizes of each column in each commitment tree.
    fn column_log_sizes(&self) -> TreeVec<ColumnVec<u32>> {
        self.trees
            .as_ref()
            .map(|tree| tree.column_log_sizes.clone())
    }

    /// Reads a commitment from the prover.
    pub fn commit(
        &mut self,
        commitment: <MC::H as MerkleHasher>::Hash,
        log_sizes: &[u32],
        channel: &mut MC::C,
    ) {
        MC::mix_root(channel, commitment);
        let extended_log_sizes = log_sizes
            .iter()
            .map(|&log_size| log_size + self.config.fri_config.log_blowup_factor)
            .collect();
        let verifier = MerkleVerifier::new(commitment, extended_log_sizes);
        self.trees.push(verifier);
    }

    pub fn verify_values(
        &self,
        sampled_points: TreeVec<ColumnVec<Vec<CirclePoint<SecureField>>>>,
        proof: CommitmentSchemeProof<MC::H>,
        channel: &mut MC::C,
    ) -> Result<(), VerificationError> {
        channel.mix_felts(&proof.sampled_values.clone().flatten_cols());
        let random_coeff = channel.draw_felt();

        let bounds = self
            .column_log_sizes()
            .zip_cols(&sampled_points)
            .map_cols(|(log_size, sampled_points)| {
                vec![
                    CirclePolyDegreeBound::new(log_size - self.config.fri_config.log_blowup_factor);
                    sampled_points.len()
                ]
            })
            .flatten_cols()
            .into_iter()
            .sorted()
            .rev()
            .dedup()
            .collect_vec();

        // FRI commitment phase on OODS quotients.
        let mut fri_verifier =
            FriVerifier::<MC>::commit(channel, self.config.fri_config, proof.fri_proof, bounds)?;

        // Verify proof of work.
        channel.mix_nonce(proof.proof_of_work);
        if channel.trailing_zeros() < self.config.pow_bits {
            return Err(VerificationError::ProofOfWork);
        }

        // Get FRI query domains.
        let fri_query_domains = fri_verifier.column_query_positions(channel);

        // Verify merkle decommitments.
        self.trees
            .as_ref()
            .zip_eq(proof.decommitments)
            .zip_eq(proof.queried_values.clone())
            .map(|((tree, decommitment), queried_values)| {
                let queries = fri_query_domains
                    .iter()
                    .map(|(&log_size, domain)| (log_size, domain.flatten()))
                    .collect();
                tree.verify(queries, queried_values, decommitment)
            })
            .0
            .into_iter()
            .collect::<Result<_, _>>()?;

        // Answer FRI queries.
        let samples = sampled_points
            .zip_cols(proof.sampled_values)
            .map_cols(|(sampled_points, sampled_values)| {
                zip(sampled_points, sampled_values)
                    .map(|(point, value)| PointSample { point, value })
                    .collect_vec()
            })
            .flatten();

        // TODO(spapini): Properly defined column log size and dinstinguish between poly and
        // commitment.
        let fri_answers = fri_answers(
            self.column_log_sizes().flatten().into_iter().collect(),
            &samples,
            random_coeff,
            fri_query_domains,
            &proof.queried_values.flatten(),
        )?;

        fri_verifier.decommit(fri_answers)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use starknet_ff::FieldElement;

    use crate::core::{backend::CpuBackend, channel::Poseidon252Channel, circle::CirclePoint, fields::{m31::M31, qm31::SecureField}, fri::FriConfig, pcs::{CommitmentSchemeProver, PcsConfig, TreeVec}, poly::circle::{CanonicCoset, CircleEvaluation, PolyOps}, vcs::poseidon252_merkle::Poseidon252MerkleChannel};
    use super::CommitmentSchemeVerifier;

    #[test]
    fn test_commitment_scheme_verifier_commit() {
        let config = PcsConfig {
            pow_bits: 10,
            fri_config: FriConfig::new(5, 4, 64),
        };
        
        let channel = &mut Poseidon252Channel::default();
        let commitment_scheme = &mut CommitmentSchemeVerifier::<Poseidon252MerkleChannel>::new(config);

        let commitment_1 = FieldElement::from_hex_be("0x01cafae415ba4e4f6648b9c8d0c44a664060485580ac65ff8c161fb756836bd5").unwrap();
        let sizes_1 = vec![10, 10, 10, 10];

        assert_eq!(commitment_scheme.trees.len(), 0);
        commitment_scheme.commit(commitment_1, &sizes_1, channel);
        assert_eq!(commitment_scheme.trees.len(), 1);
        assert_eq!(commitment_scheme.trees[0].root, FieldElement::from_hex_be("0x01cafae415ba4e4f6648b9c8d0c44a664060485580ac65ff8c161fb756836bd5").unwrap());
        assert_eq!(commitment_scheme.trees[0].column_log_sizes, vec![14, 14, 14, 14]);

        let commitment_2 = FieldElement::from_hex_be("0x0478dd9207927ad71f7c07e332b745a3d3aa08f593fcb033ef756d7438994595").unwrap();
        let sizes_2 = vec![10, 10, 10, 10, 10, 10, 10, 10];
        assert_eq!(commitment_scheme.trees.len(), 1);
        commitment_scheme.commit(commitment_2, &sizes_2, channel);
        assert_eq!(commitment_scheme.trees.len(), 2);
        assert_eq!(commitment_scheme.trees[1].root, FieldElement::from_hex_be("0x0478dd9207927ad71f7c07e332b745a3d3aa08f593fcb033ef756d7438994595").unwrap());
        assert_eq!(commitment_scheme.trees[1].column_log_sizes, vec![14, 14, 14, 14, 14, 14, 14, 14]);

        let commitment_3 = FieldElement::from_hex_be("0x032fb1cb4a18da436f91b455ef3a8153b55eab841ba8b3fee7fa33ec050356bc").unwrap();
        let sizes_3 = vec![10, 10, 10, 10];
        assert_eq!(commitment_scheme.trees.len(), 2);
        commitment_scheme.commit(commitment_3, &sizes_3, channel);
        assert_eq!(commitment_scheme.trees.len(), 3);
        assert_eq!(commitment_scheme.trees[2].root, FieldElement::from_hex_be("0x032fb1cb4a18da436f91b455ef3a8153b55eab841ba8b3fee7fa33ec050356bc").unwrap());
        assert_eq!(commitment_scheme.trees[2].column_log_sizes, vec![14, 14, 14, 14]);
    }

    #[test]
    fn test_commitment_scheme_verifier_verify_values() {
        let log_n_rows = 3;
        let log_sizes_cols = vec![log_n_rows];

        let config = PcsConfig {
            pow_bits: 10,
            fri_config: FriConfig::new(0, 1, 1),
        };

        let domain = CanonicCoset::new(log_n_rows).circle_domain();
        let extended_domain = CanonicCoset::new(log_n_rows + config.fri_config.log_blowup_factor + 1).circle_domain();

        let circle_evaluation = CircleEvaluation::new(
            domain,
            vec![M31::from(0), M31::from(1), M31::from(2), M31::from(3), M31::from(4), M31::from(5), M31::from(6), M31::from(7)]
        );

        let trace = vec![circle_evaluation];

        // Prove
        let channel = &mut Poseidon252Channel::default();
        let twiddles = CpuBackend::precompute_twiddles(extended_domain.half_coset);
        
        let commitment_scheme_prover = &mut CommitmentSchemeProver::<_, Poseidon252MerkleChannel>::new(config, &twiddles);

        let mut tree_builder = commitment_scheme_prover.tree_builder();
        tree_builder.extend_evals(trace);
        tree_builder.commit(channel);

        let oods_point = CirclePoint::<SecureField>::get_random_point(channel);
        let sample_points = TreeVec::new(vec![vec![vec![oods_point]]]);

        let commitment_scheme_prover_proof = commitment_scheme_prover.prove_values(sample_points.clone(), channel);

        // Verify
        let channel = &mut Poseidon252Channel::default();
        let commitment_scheme_verifier = &mut CommitmentSchemeVerifier::<Poseidon252MerkleChannel>::new(config);
        commitment_scheme_verifier.commit(commitment_scheme_prover.roots()[0], &log_sizes_cols, channel);

        println!("Sample points: {:?}", sample_points);
        println!("Proof: {:?}", commitment_scheme_prover_proof);
        assert!(commitment_scheme_verifier.verify_values(sample_points, commitment_scheme_prover_proof, channel).is_ok());
    }
}
