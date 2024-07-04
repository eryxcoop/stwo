use num_traits::One;

pub mod air;
mod component;

use self::air::FibonacciAir;
use self::component::FibonacciComponent;
use crate::core::backend::gpu::GpuCircleEvaluation;
use crate::core::fields::m31::BaseField;
use crate::core::fields::FieldExpOps;
use crate::core::poly::circle::{CanonicCoset, CircleEvaluation};
use crate::core::poly::BitReversedOrder;


pub struct Fibonacci {
    pub air: FibonacciAir,
}

impl Fibonacci {
    pub fn new(log_size: u32, claim: BaseField) -> Self {
        let component = FibonacciComponent::new(log_size, claim);
        Self {
            air: FibonacciAir::new(component),
        }
    }

    pub fn get_trace(&self) -> GpuCircleEvaluation<BaseField, BitReversedOrder> {
        // Trace.
        let trace_domain = CanonicCoset::new(self.air.component.log_size);
        // TODO(AlonH): Consider using Vec::new instead of Vec::with_capacity throughout file.
        let mut trace = Vec::with_capacity(trace_domain.size());

        // Fill trace with fibonacci squared.
        let mut a = BaseField::one();
        let mut b = BaseField::one();
        for _ in 0..trace_domain.size() {
            trace.push(a);
            let tmp = a.square() + b.square();
            a = b;
            b = tmp;
        }

        // Returns as a CircleEvaluation.
        CircleEvaluation::new_canonical_ordered(trace_domain, trace)
    }
}

#[cfg(test)]
mod tests {
    use std::iter::zip;

    use itertools::Itertools;

    use super::Fibonacci;
    use crate::core::air::accumulation::PointEvaluationAccumulator;
    use crate::core::air::{AirExt, AirProverExt, Component, ComponentTrace};
    use crate::core::circle::CirclePoint;
    use crate::core::fields::qm31::SecureField;
    use crate::core::poly::circle::CanonicCoset;
    use crate::{m31, qm31};

    #[test]
    fn test_composition_polynomial_is_low_degree() {
        let fib = Fibonacci::new(5, m31!(443693538));
        let trace = fib.get_trace();
        let trace_poly = trace.interpolate();
        let trace_eval =
            trace_poly.evaluate(CanonicCoset::new(trace_poly.log_size() + 1).circle_domain());
        let trace = ComponentTrace::new(vec![&trace_poly], vec![&trace_eval]);

        let random_coeff = qm31!(2213980, 2213981, 2213982, 2213983);
        let component_traces = vec![trace];
        let composition_polynomial_poly = fib
            .air
            .compute_composition_polynomial(random_coeff, &component_traces);

        // Evaluate this polynomial at another point out of the evaluation domain and compare to
        // what we expect.
        let point = CirclePoint::<SecureField>::get_point(98989892);

        let points = fib.air.mask_points(point);
        let mask_values = zip(&component_traces[0].polys, &points[0])
            .map(|(poly, points)| {
                points
                    .iter()
                    .map(|point| poly.eval_at_point(*point))
                    .collect_vec()
            })
            .collect_vec();

        let mut evaluation_accumulator = PointEvaluationAccumulator::new(random_coeff);
        fib.air.component.evaluate_constraint_quotients_at_point(
            point,
            &mask_values,
            &mut evaluation_accumulator,
        );
        let oods_value = evaluation_accumulator.finalize();
        assert_eq!(oods_value, composition_polynomial_poly.eval_at_point(point));
    }
}
