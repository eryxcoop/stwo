typedef unsigned int uint32_t;
typedef unsigned long long uint64_t;

typedef struct {
    uint32_t a;
    uint32_t b;
} cm31;

typedef struct {
    cm31 a;
    cm31 b;
} qm31;

const int LOG_MAX_NUM_CONCURRENT_THREADS = 14;
const uint32_t P = 2147483647;
const cm31 R = {2, 1};

extern "C"
__global__ void sort_values(uint32_t *from, uint32_t *dst, int size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < size) {
        if(idx < (size >> 1)) {
            dst[idx] = from[idx << 1];
        } else {
            int tmp = idx - (size >> 1);
            dst[idx] = from[size - (tmp << 1) - 1];
        }
    }
}

typedef struct {
    uint32_t x;
    uint32_t y;
} point;

__device__ uint32_t m31_mul(uint32_t a, uint32_t b) {
    // TODO: use mul from m31.cu
    uint64_t v = ((uint64_t) a * (uint64_t) b);
    uint64_t w = v + (v >> 31);
    uint64_t u = v + (w >> 31);
    return u & P;}

__device__ uint32_t m31_add(uint32_t a, uint32_t b) {
    // TODO: use add from m31.cu
    uint64_t sum = ((uint64_t) a + (uint64_t) b);
    return min(sum, sum - P);
}

__device__ uint32_t m31_sub(uint32_t a, uint32_t b) {
    // TODO: use sub from m31.cu
    return m31_add(a, P - b);
}

__device__ uint32_t m31_neg(uint32_t a) {
    // TODO: use neg from m31.cu
    return P - a;
}


/*##### CM1 ##### */

__device__ cm31 cm31_mul(cm31 x, cm31 y) {
    return {m31_sub(m31_mul(x.a, y.a), m31_mul(x.b, y.b)), m31_add(m31_mul(x.a, y.b), m31_mul(x.b, y.a))};
}

__device__ cm31 cm31_add(cm31 x, cm31 y) {
    return {m31_add(x.a, y.a), m31_add(x.b, y.b)};
}

__device__ cm31 cm31_sub(cm31 x, cm31 y) {
    return {m31_sub(x.a, y.a), m31_sub(x.b, y.b)};
}
/*##### Q31 ##### */

__device__ qm31 qm31_mul(qm31 x, qm31 y) {
    // Karatsuba multiplication
    cm31 v0 = cm31_mul(x.a, y.a);
    cm31 v1 = cm31_mul(x.b, y.b);
    cm31 v2 = cm31_mul(cm31_add(x.a, x.b), cm31_add(y.a, y.b));
    return {
        cm31_add(v0, cm31_mul(R, v1)),
        cm31_sub(v2, cm31_add(v0, v1))
    };
}

__device__ qm31 qm31_add(qm31 x, qm31 y) {
    return {cm31_add(x.a, y.a), cm31_add(x.b, y.b)};
}


/*##### Point ##### */

__device__ point point_mul(point &p1, point &p2) {
    return {
        m31_sub(m31_mul(p1.x, p2.x), m31_mul(p1.y, p2.y)),
        m31_add(m31_mul(p1.x, p2.y), m31_mul(p1.y, p2.x)),
    };
}

__device__ point point_square(point &p1) {
    return point_mul(p1, p1);
}

__device__ point one() {
    return {1, 0};
}

__device__ point pow_to_power_of_two(point p, int log_exponent) {
    int i = 0;
    while (i < log_exponent) {
        p = point_square(p);
        i++;
    }
    return p;
}

__device__ point point_pow(point p, int exponent) {
    point result = one();
    while (exponent > 0) {
        if (exponent & 1) {
            result = point_mul(p, result);
        }
        p = point_square(p);
        exponent >>= 1;
    }
    return result;
}

const point m31_circle_gen = {2, 1268011823};

__device__ unsigned int bit_reverse(unsigned int n, int bits) {
    unsigned int reversed_n = __brev(n);
    return reversed_n >> (32 - bits);
}

extern "C"
__global__ void put_one(uint32_t *dst, int offset) {
    dst[offset] = 1;
}

extern "C"
__global__ void precompute_twiddles(uint32_t *dst, point initial, point step, int offset, int size, int log_size) {
    // Computes one level of twiddles for a particular Coset.
    //      dst: twiddles array.
    //  initial: coset factor.
    //     step: generator of the group.
    //   offset: store values in dst[offset]
    //     size: coset size
    // log_size: log(size)

    // TODO: when size is larger than the max number of concurrent threads,
    //       consecutive numbers can me computed with a multiplication within the same thread,
    //       instead of using another pow.
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    size >>= 1;
    if (idx < size) {
        point pow = point_pow(step, bit_reverse(idx, log_size - 1));
        dst[offset + idx] = point_mul(initial, pow).x;
    }
}

__device__ int get_twiddle(uint32_t *twiddles, int index) {
    int k = index >> 2;
    if (index % 4 == 0) {
        return twiddles[2 * k + 1];
    } else if (index % 4 == 1) {
        return m31_neg(twiddles[2 * k + 1]);
    } else if (index % 4 == 2) {
        return m31_neg(twiddles[2 * k]);
    } else {
        return twiddles[2 * k];
    }
}

extern "C"
__global__ void ifft_circle_part(uint32_t *values, uint32_t *inverse_twiddles_tree, int values_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < (values_size >> 1)) {
        uint32_t val0 = values[2 * idx];
        uint32_t val1 = values[2 * idx + 1];
        uint32_t twiddle = get_twiddle(inverse_twiddles_tree, idx);

        values[2 * idx] = m31_add(val0, val1);
        values[2 * idx + 1] = m31_mul(m31_sub(val0, val1), twiddle);
    }
}


extern "C"
__global__ void ifft_line_part(uint32_t *values, uint32_t *inverse_twiddles_tree, int values_size, int inverse_twiddles_size, int layer_domain_offset, int layer) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < (values_size >> 1)) {
        int number_polynomials = 1 << layer;
        int h = idx >> layer;
        int l = idx & (number_polynomials - 1);
        int idx0 = (h << (layer + 1)) + l;
        int idx1 = idx0 + number_polynomials;

        uint32_t val0 = values[idx0];
        uint32_t val1 = values[idx1];
        uint32_t twiddle = inverse_twiddles_tree[layer_domain_offset + h];
        
        values[idx0] = m31_add(val0, val1);
        values[idx1] = m31_mul(m31_sub(val0, val1), twiddle);
    }
}

extern "C"
__global__ void rfft_circle_part(uint32_t *values, uint32_t *inverse_twiddles_tree, int values_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    
    if (idx < (values_size >> 1)) {
        uint32_t val0 = values[2 * idx];
        uint32_t val1 = values[2 * idx + 1];
        uint32_t twiddle = get_twiddle(inverse_twiddles_tree, idx);
        
        uint32_t temp = m31_mul(val1, twiddle);
        
        values[2 * idx] = m31_add(val0, temp);
        values[2 * idx + 1] = m31_sub(val0, temp);
    }
}


extern "C"
__global__ void rfft_line_part(uint32_t *values, uint32_t *inverse_twiddles_tree, int values_size, int inverse_twiddles_size, int layer_domain_offset, int layer) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < (values_size >> 1)) {
        int number_polynomials = 1 << layer;
        int h = idx / number_polynomials;
        int l = idx % number_polynomials;
        int idx0 = (h << (layer + 1)) + l;
        int idx1 = idx0 + number_polynomials;

        uint32_t val0 = values[idx0];
        uint32_t val1 = values[idx1];
        uint32_t twiddle = inverse_twiddles_tree[layer_domain_offset + h];
        
        uint32_t temp = m31_mul(val1, twiddle);
        
        values[idx0] = m31_add(val0, temp);
        values[idx1] = m31_sub(val0, temp);
    }
}

extern "C"
__global__ void rescale(uint32_t *values, int size, uint32_t factor) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx < size) {
        values[idx] = m31_mul(values[idx], factor);
    }
}

extern "C"
__global__ void eval_at_point_first_pass(uint32_t* g_coeffs, qm31 *temp, qm31 *factors, int coeffs_size, int factors_size, qm31 point, int output_offset) {
    int idx = threadIdx.x;

    qm31 *output = &temp[output_offset];

    // Thread syncing happens within a block. 
    // Split the problem to feed them to multiple blocks.
    if(coeffs_size >= 512) {
        coeffs_size = 512;
    }
    
    extern __shared__ uint32_t s_coeffs[];
    extern __shared__ qm31 s_level[];

    s_coeffs[idx] = g_coeffs[2 * blockIdx.x * blockDim.x + idx];
    s_coeffs[idx + blockDim.x] = g_coeffs[2 * blockIdx.x * blockDim.x + idx + blockDim.x];
    __syncthreads();
    
    int level_size = coeffs_size >> 1;
    int factor_idx = factors_size - 1;

    if(idx < level_size) {
        uint32_t alpha = s_coeffs[2 * idx];
        uint32_t beta = s_coeffs[2 * idx + 1];
        qm31 factor = factors[factor_idx];
        __syncthreads();
        qm31 result = { 
            {m31_add(m31_mul(beta, factor.a.a), alpha), m31_mul(factor.a.b, beta)}, 
            {m31_mul(beta,  factor.b.a), m31_mul(beta, factor.b.b)} 
        };
        s_level[idx] = result;
    }
    factor_idx -= 1;
    level_size >>= 1;

    while(level_size > 0) {
        if(idx < level_size) {
            __syncthreads();
            qm31 a = s_level[2 * idx];
            qm31 b = s_level[2 * idx + 1];
            __syncthreads();
            s_level[idx] = qm31_add(a, qm31_mul(b, factors[factor_idx]));
        }
        factor_idx -= 1;
        level_size >>= 1;
        
    }
    __syncthreads();

    if(idx == 0) {
        output[blockIdx.x] = s_level[0];
    }
}

extern "C"
__global__ void eval_at_point_second_pass(qm31* temp, qm31 *factors, int level_size, int factor_offset, qm31 point, int level_offset, int output_offset) {
    int idx = threadIdx.x;

    qm31 *level = &temp[level_offset];
    qm31 *output = &temp[output_offset];

    // Thread syncing happens within a block. 
    // Split the problem to feed them to multiple blocks.
    if(level_size >= 512) {
        level_size = 512;
    }
    
    extern __shared__ qm31 s_level[];

    s_level[idx] = level[2 * blockIdx.x * blockDim.x + idx];
    s_level[idx + blockDim.x] = level[2 * blockIdx.x * blockDim.x + idx + blockDim.x];

    level_size >>= 1;

    int factor_idx = factor_offset;

    while(level_size > 0) {
        if(idx < level_size) {
            __syncthreads();
            qm31 a = s_level[2 * idx];
            qm31 b = s_level[2 * idx + 1];
            __syncthreads();
            s_level[idx] = qm31_add(a, qm31_mul(b, factors[factor_idx]));
        }
        factor_idx -= 1;
        level_size >>= 1;
    }
    __syncthreads();

    if(idx == 0) {
        output[blockIdx.x] = s_level[0];
    }
}

extern "C"
__global__ void get_result_from_temp(qm31 *temp, qm31 *result) {
    if (threadIdx.x == 0) {
        result[0] = temp[0];
    }
}