use std::cmp::Reverse;
use std::collections::BTreeMap;
use std::iter::zip;

use itertools::{izip, multiunzip, Itertools};
use tracing::{span, Level};

use crate::core::backend::cpu::quotients::{accumulate_row_quotients, quotient_constants};
use crate::core::circle::CirclePoint;
use crate::core::fields::m31::BaseField;
use crate::core::fields::qm31::SecureField;
use crate::core::fri::SparseCircleEvaluation;
use crate::core::poly::circle::{
    CanonicCoset, CircleDomain, CircleEvaluation, PolyOps, SecureEvaluation,
};
use crate::core::poly::BitReversedOrder;
use crate::core::prover::VerificationError;
use crate::core::queries::SparseSubCircleDomain;
use crate::core::utils::bit_reverse_index;

pub trait QuotientOps: PolyOps {
    /// Accumulates the quotients of the columns at the given domain.
    /// For a column f(x), and a point sample (p,v), the quotient is
    ///   (f(x) - V0(x))/V1(x)
    /// where V0(p)=v, V0(conj(p))=conj(v), and V1 is a vanishing polynomial for p,conj(p).
    /// This ensures that if f(p)=v, then the quotient is a polynomial.
    /// The result is a linear combination of the quotients using powers of random_coeff.
    fn accumulate_quotients(
        domain: CircleDomain,
        columns: &[&CircleEvaluation<Self, BaseField, BitReversedOrder>],
        random_coeff: SecureField,
        sample_batches: &[ColumnSampleBatch],
        log_blowup_factor: u32,
    ) -> SecureEvaluation<Self, BitReversedOrder>;
}

/// A batch of column samplings at a point.
pub struct ColumnSampleBatch {
    /// The point at which the columns are sampled.
    pub point: CirclePoint<SecureField>,
    /// The sampled column indices and their values at the point.
    pub columns_and_values: Vec<(usize, SecureField)>,
}

impl ColumnSampleBatch {
    /// Groups column samples by sampled point.
    /// # Arguments
    /// samples: For each column, a vector of samples.
    pub fn new_vec(samples: &[&Vec<PointSample>]) -> Vec<Self> {
        // Group samples by point, and create a ColumnSampleBatch for each point.
        // This should keep a stable ordering.
        let mut grouped_samples = BTreeMap::new();
        for (column_index, samples) in samples.iter().enumerate() {
            for sample in samples.iter() {
                grouped_samples
                    .entry(sample.point)
                    .or_insert_with(Vec::new)
                    .push((column_index, sample.value));
            }
        }
        grouped_samples
            .into_iter()
            .map(|(point, columns_and_values)| ColumnSampleBatch {
                point,
                columns_and_values,
            })
            .collect()
    }
}

#[derive(Debug)]
pub struct PointSample {
    pub point: CirclePoint<SecureField>,
    pub value: SecureField,
}

pub fn compute_fri_quotients<B: QuotientOps>(
    columns: &[&CircleEvaluation<B, BaseField, BitReversedOrder>],
    samples: &[Vec<PointSample>],
    random_coeff: SecureField,
    log_blowup_factor: u32,
) -> Vec<SecureEvaluation<B, BitReversedOrder>> {
    let _span = span!(Level::INFO, "Compute FRI quotients").entered();
    zip(columns, samples)
        .sorted_by_key(|(c, _)| Reverse(c.domain.log_size()))
        .group_by(|(c, _)| c.domain.log_size())
        .into_iter()
        .map(|(log_size, tuples)| {
            let (columns, samples): (Vec<_>, Vec<_>) = tuples.unzip();
            let domain = CanonicCoset::new(log_size).circle_domain();
            // TODO: slice.
            let sample_batches = ColumnSampleBatch::new_vec(&samples);
            B::accumulate_quotients(
                domain,
                &columns,
                random_coeff,
                &sample_batches,
                log_blowup_factor,
            )
        })
        .collect()
}

pub fn fri_answers(
    column_log_sizes: Vec<u32>,
    samples: &[Vec<PointSample>],
    random_coeff: SecureField,
    query_domain_per_log_size: BTreeMap<u32, SparseSubCircleDomain>,
    queried_values_per_column: &[Vec<BaseField>],
) -> Result<Vec<SparseCircleEvaluation>, VerificationError> {
    izip!(column_log_sizes, samples, queried_values_per_column)
        .sorted_by_key(|(log_size, ..)| Reverse(*log_size))
        .group_by(|(log_size, ..)| *log_size)
        .into_iter()
        .map(|(log_size, tuples)| {
            let (_, samples, queried_valued_per_column): (Vec<_>, Vec<_>, Vec<_>) =
                multiunzip(tuples);
            fri_answers_for_log_size(
                log_size,
                &samples,
                random_coeff,
                &query_domain_per_log_size[&log_size],
                &queried_valued_per_column,
            )
        })
        .collect()
}

pub fn fri_answers_for_log_size(
    log_size: u32,
    samples: &[&Vec<PointSample>],
    random_coeff: SecureField,
    query_domain: &SparseSubCircleDomain,
    queried_values_per_column: &[&Vec<BaseField>],
) -> Result<SparseCircleEvaluation, VerificationError> {
    let commitment_domain = CanonicCoset::new(log_size).circle_domain();
    let sample_batches = ColumnSampleBatch::new_vec(samples);
    for queried_values in queried_values_per_column {
        if queried_values.len() != query_domain.flatten().len() {
            return Err(VerificationError::InvalidStructure(
                "Insufficient number of queried values".to_string(),
            ));
        }
    }
    let mut queried_values_per_column = queried_values_per_column
        .iter()
        .map(|q| q.iter())
        .collect_vec();

    let mut evals = Vec::new();
    for subdomain in query_domain.iter() {
        let domain = subdomain.to_circle_domain(&commitment_domain);
        let quotient_constants = quotient_constants(&sample_batches, random_coeff, domain);
        let mut column_evals = Vec::new();
        for queried_values in queried_values_per_column.iter_mut() {
            let eval = CircleEvaluation::new(
                domain,
                queried_values.take(domain.size()).copied().collect_vec(),
            );
            column_evals.push(eval);
        }

        let mut values = Vec::new();
        for row in 0..domain.size() {
            let domain_point = domain.at(bit_reverse_index(row, log_size));
            let value = accumulate_row_quotients(
                &sample_batches,
                &column_evals.iter().collect_vec(),
                &quotient_constants,
                row,
                domain_point,
            );
            values.push(value);
        }
        let eval = CircleEvaluation::new(domain, values);
        evals.push(eval);
    }

    let res = SparseCircleEvaluation::new(evals);
    if !queried_values_per_column.iter().all(|x| x.is_empty()) {
        return Err(VerificationError::InvalidStructure(
            "Too many queried values".to_string(),
        ));
    }
    Ok(res)
}

#[cfg(test)]
mod tests {

    use crate::core::backend::cpu::quotients::{accumulate_row_quotients, QuotientConstants};
    use crate::core::backend::cpu::{CpuCircleEvaluation, CpuCirclePoly};
    use crate::core::circle::{CirclePoint, CirclePointIndex, Coset, SECURE_FIELD_CIRCLE_GEN};
    use crate::core::fields::cm31::CM31;
    use crate::core::fields::qm31::QM31;
    use crate::core::pcs::quotients::{compute_fri_quotients, ColumnSampleBatch, PointSample};
    use crate::core::poly::circle::{CanonicCoset, CircleDomain, CircleEvaluation};
    use crate::{m31, qm31};

    #[test]
    fn test_quotients_are_low_degree() {
        const LOG_SIZE: u32 = 7;
        const LOG_BLOWUP_FACTOR: u32 = 1;
        let polynomial = CpuCirclePoly::new((0..1 << LOG_SIZE).map(|i| m31!(i)).collect());
        let eval_domain = CanonicCoset::new(LOG_SIZE + 1).circle_domain();
        let eval = polynomial.evaluate(eval_domain);
        let point = SECURE_FIELD_CIRCLE_GEN;
        let value = polynomial.eval_at_point(point);
        let coeff = qm31!(1, 2, 3, 4);
        let quot_eval = compute_fri_quotients(
            &[&eval],
            &[vec![PointSample { point, value }]],
            coeff,
            LOG_BLOWUP_FACTOR,
        )
        .pop()
        .unwrap();
        let quot_poly_base_field =
            CpuCircleEvaluation::new(eval_domain, quot_eval.values.columns[0].clone())
                .interpolate();
        assert!(quot_poly_base_field.is_in_fri_space(LOG_SIZE));
    }

    #[test]
    fn test_accumulate_row_quotients() {
        let sample_batches = [ColumnSampleBatch {
            point: CirclePoint {
                x: QM31(
                    CM31::from_m31(m31!(1395048677), m31!(640591314)),
                    CM31::from_m31(m31!(342871101), m31!(1049385418)),
                ),
                y: QM31(
                    CM31::from_m31(m31!(474688795), m31!(2119282552)),
                    CM31::from_m31(m31!(160740005), m31!(798859953)),
                ),
            },
            columns_and_values: vec![(
                0,
                QM31(
                    CM31::from_m31(m31!(2082657879), m31!(1175528048)),
                    CM31::from_m31(m31!(1000432343), m31!(763013627)),
                ),
            )],
        }];
        let column_evals = vec![CircleEvaluation::new(
            CircleDomain {
                half_coset: Coset::new(CirclePointIndex(1946157056), 0),
            },
            vec![m31!(1323727772), m31!(1323695004)],
        )];
        let quotient_constants = QuotientConstants {
            line_coeffs: vec![vec![(
                QM31(
                    CM31::from_m31(m31!(1507438696), m31!(1485287082)),
                    CM31::from_m31(m31!(462219637), m31!(1981546355)),
                ),
                QM31(
                    CM31::from_m31(m31!(747975780), m31!(1185576571)),
                    CM31::from_m31(m31!(718995102), m31!(494594746)),
                ),
                QM31(
                    CM31::from_m31(m31!(1714794137), m31!(1220651711)),
                    CM31::from_m31(m31!(676173891), m31!(2009908523)),
                ),
            )]],
            batch_random_coeffs: vec![QM31(
                CM31::from_m31(m31!(1884891309), m31!(790325555)),
                CM31::from_m31(m31!(1984078504), m31!(1012128825)),
            )],
            denominator_inverses: vec![vec![
                CM31::from_m31(m31!(2087781265), m31!(1887176073)),
                CM31::from_m31(m31!(1973951885), m31!(1467411716)),
            ]],
        }; //, array![CM31 { a: M31 { inner: 2087781265 }, b: M31 { inner: 1887176073 } }, CM31 { a:
           //, M31 { inner: 1973951885 }, b: M31 { inner: 1467411716 } }]] };
        let row = 0;
        let domain_point = CirclePoint {
            x: m31!(34602070),
            y: m31!(1415090252),
        };
        let expected_value = QM31(
            CM31::from_m31(m31!(50599093), m31!(1497794168)),
            CM31::from_m31(m31!(254541753), m31!(369788671)),
        );
        let value = accumulate_row_quotients(
            &sample_batches,
            &[&column_evals[0]],
            &quotient_constants,
            row,
            domain_point,
        );
        println!("{:?}", value);
        assert_eq!(expected_value, value);
    }
}
