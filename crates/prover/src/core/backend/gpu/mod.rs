mod accumulation;
mod bit_reverse;
mod circle;
pub mod column;
pub mod error;
mod fri;
pub mod m31;
mod quotients;
// pub mod packedm31;

use std::ffi::c_void;
use std::fmt::Debug;
use std::sync::Arc;

use cudarc::driver::{CudaDevice, DeviceRepr, LaunchConfig, ValidAsZeroBits};
// use error::Error;
use once_cell::sync::Lazy;

use self::m31::LoadBaseField;
use super::Backend;
use crate::core::fields::m31::M31;
use crate::core::fields::qm31::QM31;

// TODO:: cleanup unwraps with error handling?
// (We can replace lazy statics with unsafe global references)
static DEVICE: Lazy<Arc<CudaDevice>> = Lazy::new(|| CudaDevice::new(0).unwrap().load());
// static MODULUS: Lazy<CudaSlice<u32>> =
//     Lazy::new(|| DEVICE.htod_copy([P; VECTOR_SIZE].to_vec()).unwrap());

type Device = Arc<CudaDevice>;

trait Load {
    fn load(self) -> Self;
}

impl Load for Device {
    fn load(self) -> Self {
        bit_reverse::load_bit_reverse_ptx(&self);
        LoadBaseField::load(&self);
        column::load_batch_inverse_ptx(&self);
        circle::load_circle(&self);
        self
    }
}

unsafe impl DeviceRepr for M31 {
    fn as_kernel_param(&self) -> *mut c_void {
        self as *const Self as *mut c_void
    }
}

unsafe impl DeviceRepr for QM31 {
    fn as_kernel_param(&self) -> *mut c_void {
        self as *const Self as *mut c_void
    }
}

unsafe impl DeviceRepr for &mut QM31 {
    fn as_kernel_param(&self) -> *mut std::ffi::c_void {
        self as *const Self as *mut _
    }
}

unsafe impl ValidAsZeroBits for M31 {}
unsafe impl ValidAsZeroBits for QM31 {}

#[derive(Copy, Clone, Debug)]
pub struct GpuBackend;

impl Backend for GpuBackend {}

impl GpuBackend {
    /// Creates a [LaunchConfig] with:
    /// - block_dim == `1024`
    /// - grid_dim == `(n + 1023) / 1024`
    /// - shared_mem_bytes == `0`
    pub fn launch_config_for_num_elems(n: u32, num_threads: u32, shared_mem_bytes: u32) -> LaunchConfig {
        let num_blocks = (n + num_threads - 1) / num_threads;
        LaunchConfig {
            grid_dim: (num_blocks, 1, 1),
            block_dim: (num_threads, 1, 1),
            shared_mem_bytes,
        }
    }
}
