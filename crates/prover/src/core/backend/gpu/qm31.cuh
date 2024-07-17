
// qm31.cuh
#ifndef QM31_H
#define QM31_H

#include "m31.cuh"

__device__ void mul_qm31(unsigned int *lhs, unsigned int *rhs, unsigned int *out);
__device__ void mul_cm31(unsigned int *lhs, unsigned int *rhs, unsigned int *out);

#endif 

struct cm31 {
    m31 a;
    m31 b;

    __device__ __host__ cm31() : a(0), b(0) {}
    __device__ __host__ cm31(m31 a, m31 b) : a(a), b(b) {}

    __device__ __inline__ cm31 operator*(const cm31& rhs) const {
        unsigned int ac = a.f * rhs.a.f;
        unsigned int bd = b.f * rhs.b.f;

        unsigned int ab_t_cd = (a.f + b.f) * (rhs.a.f + rhs.b.f);
        
        return cm31(ac - bd, ab_t_cd - (ac + bd));
    }

    __device__ __inline__ cm31 operator+(const cm31& rhs) const {
        return cm31(a + rhs.a, b + rhs.b);
    }

    __device__ __inline__ cm31 operator-(const cm31& rhs) const {
        return cm31(a - rhs.a, b - rhs.b);
    }
};

struct qm31 {
    cm31 a;
    cm31 b;

    __device__ __host__ qm31() : a(cm31()), b(cm31()) {}
    __device__ __host__ qm31(cm31 a, cm31 b) : a(a), b(b) {}

    __device__ __inline__ qm31 operator*(const qm31& rhs) const {
        cm31 ac, bd, bd_times_1_plus_i, ac_p_bd, ad_p_bc, l; 

        ac = a * rhs.b;
        bd = b * rhs.b; 

        bd_times_1_plus_i = cm31(bd.a - bd.b, bd.a + bd.b); 

        ac_p_bd = ac + bd; 
        ad_p_bc = ((a + b) * (rhs.a + rhs.b)) - ac_p_bd;

        l = cm31(ac_p_bd.a + bd_times_1_plus_i.a, ac_p_bd.b + bd_times_1_plus_i.b);
        return qm31(l, ad_p_bc);
    }
};


// TODO (Daniel): Eventually replace with struct implementation
__device__  void mul_cm31(unsigned int *lhs,  unsigned int *rhs,  unsigned int *out) {
    unsigned int ac = mul_m31(lhs[0], rhs[0]);
    unsigned int bd = mul_m31(lhs[1], rhs[1]);

    unsigned int ab_t_cd = mul_m31(add_m31(lhs[0], lhs[1]), add_m31(rhs[0], rhs[1])); 
    out[0] = sub_m31(ac, bd); 
    out[1] = sub_m31(ab_t_cd, add_m31(ac, bd)); 
}

__device__  void mul_qm31(unsigned int *lhs, unsigned int *rhs, unsigned int *out) {
    unsigned int ac[2];
    unsigned int bd[2];
    unsigned int bd_times_1_plus_i[2];
    unsigned int ac_p_bd[2];
    unsigned int ad_p_bc[2];
    unsigned int l[2];

    mul_cm31(lhs, rhs, ac);
    mul_cm31(lhs + 2, rhs + 2, bd);

    bd_times_1_plus_i[0] = sub_m31(bd[0], bd[1]);
    bd_times_1_plus_i[1] = add_m31(bd[0], bd[1]);

    ac_p_bd[0] = add_m31(ac[0], bd[0]);
    ac_p_bd[1] = add_m31(ac[1], bd[1]);

    unsigned int lhs_a_plus_b[2];
    unsigned int rhs_a_plus_b[2];
    unsigned int res[2]; 

    lhs_a_plus_b[0] = add_m31(lhs[0], lhs[2]);
    lhs_a_plus_b[1] = add_m31(lhs[1], lhs[3]);

    rhs_a_plus_b[0] = add_m31(rhs[0], rhs[2]);
    rhs_a_plus_b[1] = add_m31(rhs[1], rhs[3]);

    mul_cm31(lhs_a_plus_b, rhs_a_plus_b, res);

    ad_p_bc[0] = sub_m31(res[0], ac_p_bd[0]);
    ad_p_bc[1] = sub_m31(res[1], ac_p_bd[1]);

    l[0] = add_m31(ac_p_bd[0], bd_times_1_plus_i[0]);
    l[1] = add_m31(ac_p_bd[1], bd_times_1_plus_i[1]);
    
    out[0] = l[0];
    out[1] = l[1]; 
    out[2] = ad_p_bc[0];
    out[3] = ad_p_bc[1];
}