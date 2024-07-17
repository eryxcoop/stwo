#include "qm31.cuh"

typedef struct {
    uint32_t qm31;
    uint32_t qm31;
} circle_point;

/// Point vanishing for the packed representation of the points. skips the division.
/// See [crate::core::constraints::point_vanishing_fraction] for more details.
fn packed_point_vanishing_fraction(
    excluded: CirclePoint<SecureField>,
    p: (PackedM31, PackedM31),
) -> (PackedQM31, PackedQM31) {
    let e_conjugate = excluded.conjugate();
    let h_x = PackedQM31::broadcast(e_conjugate.x) * p.0
        - PackedQM31::broadcast(e_conjugate.y) * p.1;
    let h_y = PackedQM31::broadcast(e_conjugate.y) * p.0
        + PackedQM31::broadcast(e_conjugate.x) * p.1;
    (h_y, PackedQM31::one() + h_x)
}

__device__ 