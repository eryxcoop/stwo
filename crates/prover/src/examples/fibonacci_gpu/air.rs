use super::component::FibonacciComponent;
use crate::core::air::{Air, AirProver, Component, ComponentProver};
use crate::core::backend::CpuBackend;

pub struct FibonacciAir {
    pub component: FibonacciComponent,
}

impl FibonacciAir {
    pub fn new(component: FibonacciComponent) -> Self {
        Self { component }
    }
}
impl Air for FibonacciAir {
    fn components(&self) -> Vec<&dyn Component> {
        vec![&self.component]
    }
}
impl AirProver<CpuBackend> for FibonacciAir {
    fn prover_components(&self) -> Vec<&dyn ComponentProver<CpuBackend>> {
        vec![&self.component]
    }
}
