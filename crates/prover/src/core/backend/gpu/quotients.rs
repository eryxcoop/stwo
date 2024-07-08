use super::GpuBackend;
use crate::core::fields::m31::BaseField;
use crate::core::fields::qm31::SecureField;
use crate::core::fields::secure_column::SecureColumn;
use crate::core::pcs::quotients::{ColumnSampleBatch, QuotientOps};
use crate::core::poly::circle::{CircleDomain, CircleEvaluation, SecureEvaluation};
use crate::core::poly::BitReversedOrder;

impl QuotientOps for GpuBackend {
    fn accumulate_quotients(
        domain: CircleDomain,
        columns: &[&CircleEvaluation<Self, BaseField, BitReversedOrder>],
        random_coeff: SecureField,
        sample_batches: &[ColumnSampleBatch],
    ) -> SecureEvaluation<Self> {
        let mut values = SecureColumn::<Self>::zeros(domain.size());
        
    }
}
