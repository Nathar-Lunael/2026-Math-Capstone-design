// [제 1실험] secp256k1 + DaA ECC - ESP-IDF Framework 버전 (Performance Benchmark)
// (워치독 방지 최적화 반영 완료)

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>

// FreeRTOS 및 ESP32 커널 헤더 포함
#include "freertos/FreeRTOS.h"
#include "freertos/task.h"
#include "xtensa/core-macros.h"

#include "test_vectors.h"  // 파이썬으로 자동 생성한 100개의 무작위 스칼라 헤더 포함

// ====================================================================
// [0. 전역 상수 및 구조체 정의 (secp256k1)]
// ====================================================================

const uint64_t P[4] = {                                         
    0xFFFFFFFEFFFFFC2FULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL                                       
};
const uint64_t SECP_C = 0x1000003D1ULL;                         
const uint64_t B_VAL[4] = { 7, 0, 0, 0 };                       

const uint64_t ZERO[4] = { 0, 0, 0, 0 };
const uint64_t ONE[4] = { 1, 0, 0, 0 };
const uint64_t TWO[4] = { 2, 0, 0, 0 };
const uint64_t THREE[4] = { 3, 0, 0, 0 };
const uint64_t FOUR[4] = { 4, 0, 0, 0 };                        
const uint64_t EIGHT[4] = { 8, 0, 0, 0 };                       

typedef struct {                                         
    uint64_t X[4];
    uint64_t Y[4];
    int is_infty;                                        
} AffinePoint;

typedef struct {                                         
    uint64_t X[4];
    uint64_t Y[4];
    uint64_t Z[4];
    int is_infty;
} JacobianPoint;

// ====================================================================
// [1. 다중 정밀도 연산 함수 (Multi-Precision Basics)]
// ====================================================================

// ESP32 32-bit 컴파일러 호환용 128-bit 곱셈 함수 (MSVC _umul128 완벽 대체)
static inline uint64_t mock_umul128(uint64_t a, uint64_t b, uint64_t* high) {
    uint64_t a_lo = a & 0xFFFFFFFFULL;
    uint64_t a_hi = a >> 32;
    uint64_t b_lo = b & 0xFFFFFFFFULL;
    uint64_t b_hi = b >> 32;

    uint64_t p0 = a_lo * b_lo;
    uint64_t p1 = a_lo * b_hi;
    uint64_t p2 = a_hi * b_lo;
    uint64_t p3 = a_hi * b_hi;

    uint64_t mid = (p0 >> 32) + (p1 & 0xFFFFFFFFULL) + (p2 & 0xFFFFFFFFULL);
    *high = p3 + (p1 >> 32) + (p2 >> 32) + (mid >> 32);
    return (mid << 32) | (p0 & 0xFFFFFFFFULL);
}

// ESP32 호환용 addcarry / subborrow 대체 함수들
static inline uint8_t mock_addcarry_u64(uint8_t c_in, uint64_t a, uint64_t b, uint64_t *out) {
    uint64_t sum = a + c_in;
    uint8_t carry1 = (sum < a) ? 1 : 0;
    uint64_t total = sum + b;
    uint8_t carry2 = (total < b) ? 1 : 0;
    *out = total;
    return carry1 | carry2;
}

static inline uint8_t mock_subborrow_u64(uint8_t b_in, uint64_t a, uint64_t b, uint64_t *out) {
    uint64_t diff = a - b_in;
    uint8_t borrow1 = (a < b_in) ? 1 : 0;
    uint64_t total = diff - b;
    uint8_t borrow2 = (diff < b) ? 1 : 0;
    *out = total;
    return borrow1 | borrow2;
}

int geq_P(const uint64_t* a) {
    for (int i = 3; i >= 0; i--) {
        if (a[i] > P[i]) return 1; 
        if (a[i] < P[i]) return 0; 
    }
    return 1; 
}

int is_equal(const uint64_t* a, const uint64_t* b) {
    for (int i = 0; i < 4; i++) {
        if (a[i] != b[i]) return 0;
    }
    return 1;
}

int is_zero(const uint64_t* a) {
    for (int i = 0; i < 4; i++) {
        if (a[i] != 0) return 0;
    }
    return 1;
}

// ====================================================================
// [2. 모듈러 연산 함수]
// ====================================================================

void mod_add(uint64_t* res, const uint64_t* a, const uint64_t* b) {
    uint64_t temp[4];
    uint8_t carry = 0;

    for (int i = 0; i < 4; i++) {
        carry = mock_addcarry_u64(carry, a[i], b[i], &temp[i]);
    }

    while (geq_P(temp)) {
        uint8_t borrow = 0;
        for (int i = 0; i < 4; i++) {
            borrow = mock_subborrow_u64(borrow, temp[i], P[i], &temp[i]);
        }
    }
    for (int i = 0; i < 4; i++) res[i] = temp[i];
}

void mod_sub(uint64_t* res, const uint64_t* a, const uint64_t* b) {
    uint64_t temp[4];
    uint8_t borrow = 0;

    for (int i = 0; i < 4; i++) {
        borrow = mock_subborrow_u64(borrow, a[i], b[i], &temp[i]);
    }

    if (borrow) {
        uint8_t carry = 0;
        for (int i = 0; i < 4; i++) {
            carry = mock_addcarry_u64(carry, temp[i], P[i], &temp[i]);
        }
    }
    for (int i = 0; i < 4; i++) res[i] = temp[i];
}

void mod_mul(uint64_t* res, const uint64_t* a, const uint64_t* b) {
    uint64_t t[8] = { 0 };

    for (int i = 0; i < 4; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < 4; j++) {
            uint64_t hi, lo;
            lo = mock_umul128(a[i], b[j], &hi);

            uint8_t c1 = mock_addcarry_u64(0, t[i + j], lo, &t[i + j]);
            uint8_t c2 = mock_addcarry_u64(0, t[i + j], carry, &t[i + j]);
            carry = hi + c1 + c2;
        }
        t[i + 4] = carry;
    }

    uint64_t H[4] = { t[4], t[5], t[6], t[7] }; 
    uint64_t L[4] = { t[0], t[1], t[2], t[3] }; 

    uint64_t hiC[5] = { 0 };
    uint64_t carryC = 0;
    for (int i = 0; i < 4; i++) {
        uint64_t hi, lo;
        lo = mock_umul128(H[i], SECP_C, &hi);
        uint8_t c = mock_addcarry_u64(0, lo, carryC, &hiC[i]);
        carryC = hi + c;
    }
    hiC[4] = carryC;

    uint64_t L_res[4];
    uint8_t c_L = 0;
    for (int i = 0; i < 4; i++) {
        c_L = mock_addcarry_u64(c_L, L[i], hiC[i], &L_res[i]);
    }

    uint64_t over = hiC[4] + c_L;
    uint64_t temp[4];

    uint64_t c2_hi, c2_lo;
    c2_lo = mock_umul128(over, SECP_C, &c2_hi);

    uint8_t cc = mock_addcarry_u64(0, L_res[0], c2_lo, &temp[0]);
    cc = mock_addcarry_u64(cc, L_res[1], c2_hi, &temp[1]);
    cc = mock_addcarry_u64(cc, L_res[2], 0, &temp[2]);
    cc = mock_addcarry_u64(cc, L_res[3], 0, &temp[3]);

    if (cc) {
        cc = mock_addcarry_u64(0, temp[0], SECP_C, &temp[0]);
        for (int i = 1; i < 4; i++) cc = mock_addcarry_u64(cc, temp[i], 0, &temp[i]);
    }

    while (geq_P(temp)) {
        uint8_t br = 0;
        for (int i = 0; i < 4; i++) br = mock_subborrow_u64(br, temp[i], P[i], &temp[i]);
    }
    for (int i = 0; i < 4; i++) res[i] = temp[i];
}

void mod_inv(uint64_t* res, const uint64_t* a) {
    uint64_t r[4] = { 1, 0, 0, 0 };
    uint64_t base[4];
    for (int i = 0; i < 4; i++) base[i] = a[i];

    const uint64_t PM2[4] = {
        0xFFFFFFFEFFFFFC2DULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL
    };

    for (int i = 0; i < 256; i++) {
        uint64_t bit = (PM2[i / 64] >> (i % 64)) & 1;
        if (bit) mod_mul(r, r, base);
        mod_mul(base, base, base);
    }
    for (int i = 0; i < 4; i++) res[i] = r[i];
}

// ====================================================================
// [3. 타원곡선 점 연산 (Jacobian Coordinates)]
// ====================================================================

void point_double_jacobian(JacobianPoint* res, const JacobianPoint* pt) {
    if (pt->is_infty) {
        res->is_infty = 1;
        return;
    }

    JacobianPoint out; 
    out.is_infty = 0;

    uint64_t XX[4], YY[4], YYYY[4], S[4], M[4], T[4];
    uint64_t X3[4], Y3[4], Z3[4];

    mod_mul(XX, pt->X, pt->X);         
    mod_mul(YY, pt->Y, pt->Y);         
    mod_mul(YYYY, YY, YY);             

    mod_mul(S, pt->X, YY);
    mod_mul(S, S, FOUR);

    mod_mul(M, XX, THREE);

    mod_mul(X3, M, M);
    mod_mul(T, S, TWO);
    mod_sub(X3, X3, T);

    mod_sub(T, S, X3);
    mod_mul(Y3, M, T);
    mod_mul(T, YYYY, EIGHT);
    mod_sub(Y3, Y3, T);

    mod_mul(Z3, pt->Y, pt->Z);
    mod_mul(Z3, Z3, TWO);

    for (int i = 0; i < 4; i++) {
        out.X[i] = X3[i];
        out.Y[i] = Y3[i];
        out.Z[i] = Z3[i];
    }
    *res = out;
}

void point_add_jacobian(JacobianPoint* res, const JacobianPoint* p1, const JacobianPoint* p2) {
    if (p1->is_infty) { *res = *p2; return; }
    if (p2->is_infty) { *res = *p1; return; }

    JacobianPoint out; 
    out.is_infty = 0;

    uint64_t Z1Z1[4], Z2Z2[4], U1[4], U2[4], S1[4], S2[4], H[4], R[4];
    uint64_t H2[4], H3[4], U1H2[4], X3[4], Y3[4], Z3[4], T[4], Z1Z2[4];

    mod_mul(Z1Z1, p1->Z, p1->Z);       
    mod_mul(Z2Z2, p2->Z, p2->Z);       

    mod_mul(U1, p1->X, Z2Z2);
    mod_mul(U2, p2->X, Z1Z1);

    mod_mul(T, p2->Z, Z2Z2);           
    mod_mul(S1, p1->Y, T);             

    mod_mul(T, p1->Z, Z1Z1);           
    mod_mul(S2, p2->Y, T);             

    if (is_equal(U1, U2)) {
        if (is_equal(S1, S2)) {
            point_double_jacobian(res, p1);
            return;
        }
        else {
            res->is_infty = 1;
            return;
        }
    }

    mod_sub(H, U2, U1);                
    mod_sub(R, S2, S1);                

    mod_mul(H2, H, H);                 
    mod_mul(H3, H, H2);                
    mod_mul(U1H2, U1, H2);             

    mod_mul(X3, R, R);
    mod_sub(X3, X3, H3);
    mod_mul(T, U1H2, TWO);
    mod_sub(X3, X3, T);

    mod_sub(T, U1H2, X3);
    mod_mul(Y3, R, T);
    mod_mul(T, S1, H3);
    mod_sub(Y3, Y3, T);

    mod_mul(Z1Z2, p1->Z, p2->Z);
    mod_mul(Z3, Z1Z2, H);

    for (int i = 0; i < 4; i++) {
        out.X[i] = X3[i];
        out.Y[i] = Y3[i];
        out.Z[i] = Z3[i];
    }
    *res = out;
}

// ====================================================================
// [4. Double-and-Add 알고리즘]
// ====================================================================

void scalar_mult_jacobian(JacobianPoint* res, const uint64_t* scalar, const AffinePoint* base_point) {
    JacobianPoint P_base;
    for (int i = 0; i < 4; i++) {
        P_base.X[i] = base_point->X[i];
        P_base.Y[i] = base_point->Y[i];
        P_base.Z[i] = (i == 0) ? 1 : 0; 
    }
    P_base.is_infty = base_point->is_infty;

    JacobianPoint Q;
    Q.is_infty = 1;

    for (int i = 255; i >= 0; i--) {
        if (!Q.is_infty) {
            point_double_jacobian(&Q, &Q);
        }

        uint64_t bit = (scalar[i / 64] >> (i % 64)) & 1;
        if (bit) {
            point_add_jacobian(&Q, &Q, &P_base);
        }
    }
    *res = Q;
}

// ====================================================================
// [5. 최종 좌표 변환]
// ====================================================================

void jacobian_to_affine(AffinePoint* res, const JacobianPoint* pt) {
    if (pt->is_infty) {
        res->is_infty = 1;
        return;
    }

    uint64_t z_inv[4], z_inv2[4], z_inv3[4];

    mod_inv(z_inv, pt->Z);

    mod_mul(z_inv2, z_inv, z_inv);     
    mod_mul(z_inv3, z_inv2, z_inv);    

    mod_mul(res->X, pt->X, z_inv2);    
    mod_mul(res->Y, pt->Y, z_inv3);    

    res->is_infty = 0;
}

int is_on_curve_secp256k1(const AffinePoint* pt) {          
    if (pt->is_infty)
        return 1;

    uint64_t x2[4], x3[4], rhs[4], lhs[4];

    mod_mul(x2, pt->X, pt->X);    
    mod_mul(x3, x2, pt->X);       
    mod_add(rhs, x3, B_VAL);
    mod_mul(lhs, pt->Y, pt->Y);

    return is_equal(lhs, rhs);
}

void print_256(const char* label, const uint64_t* num) {    
    printf("%s:\n%016" PRIx64 " %016" PRIx64 " %016" PRIx64 " %016" PRIx64 "\n", label, num[3], num[2], num[1], num[0]);
}

// ====================================================================
// [FreeRTOS 벤치마킹 실행 Task]
// ====================================================================
void secp256k1_benchmark_task(void* pvParameters) {
    printf("=== secp256k1 Double-and-Add ESP-IDF Start ===\n\n");

    uint64_t G_X[4] = { 0x59F2815B16F81798ULL, 0x029BFCDB2DCE28D9ULL, 0x55A06295CE870B07ULL, 0x79BE667EF9DCBBACULL };
    uint64_t G_Y[4] = { 0x9C47D08FFB10D4B8ULL, 0xFD17B448A6855419ULL, 0x5DA4FBFC0E1108A8ULL, 0x483ADA7726A3C465ULL };

    AffinePoint Base_Point;
    for (int i = 0; i < 4; i++) {
        Base_Point.X[i] = G_X[i];
        Base_Point.Y[i] = G_Y[i];
    }
    Base_Point.is_infty = 0;

    if (!is_on_curve_secp256k1(&Base_Point)) {
        printf("[X] CRITICAL ERROR: Base point is NOT on the secp256k1 curve.\n");
        vTaskDelete(NULL);
    }

    uint64_t Private_Key_All_Ones[4] = { 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL };
    uint64_t Private_Key_All_Zeros[4] = { 0, 0, 0, 0 };

    JacobianPoint J_Res;
    AffinePoint A_Res;

    // Xtensa 하드웨어 타이머 변수 (32비트)
    uint32_t start_cycle, end_cycle;

    const int TEST1_SCALARS = 100;
    const int TEST23_ITERATIONS = 1000;

    print_256("Base_X", Base_Point.X);
    print_256("Base_Y", Base_Point.Y);
    printf("\n");

    // ==========================================================
    // [테스트 1] 무작위 256비트 스칼라 100개 평균 측정 (Average Case)
    // ==========================================================
    printf("--- Test 1 : 100 Random Scalars (Average Case) ---\n");
    uint64_t sum_cycles_test1 = 0;

    for (int iter = 0; iter < TEST1_SCALARS; iter++) {
        start_cycle = xthal_get_ccount(); // 사이클 카운터 시작
        scalar_mult_jacobian(&J_Res, Random_Scalars[iter], &Base_Point);
        end_cycle = xthal_get_ccount();   // 사이클 카운터 종료
        
        sum_cycles_test1 += (end_cycle - start_cycle);

        // 💡 핵심 최적화: 측정이 끝난 후 OS에게 잠시 CPU를 양보하여 워치독 리셋 방지
        vTaskDelay(1);
    }

    jacobian_to_affine(&A_Res, &J_Res);
    print_256("Result_1_X (100th run)", A_Res.X);
    printf("-> Average CPU Cycles (100 runs): %" PRIu64 "\n\n", sum_cycles_test1 / TEST1_SCALARS);

    // ==========================================================
    // [테스트 2] 모든 비트가 1인 키 스칼라 곱셈 (1,000번 반복 평균)
    // ==========================================================
    printf("--- Test 2 : All Bits are 1 (1,000 Iterations) ---\n");
    uint64_t sum_cycles_test2 = 0;

    for (int iter = 0; iter < TEST23_ITERATIONS; iter++) {
        start_cycle = xthal_get_ccount();
        scalar_mult_jacobian(&J_Res, Private_Key_All_Ones, &Base_Point);
        end_cycle = xthal_get_ccount();
        
        sum_cycles_test2 += (end_cycle - start_cycle);

        // 💡 워치독 방지
        vTaskDelay(1);
    }

    jacobian_to_affine(&A_Res, &J_Res);
    print_256("Result_2_X", A_Res.X);
    printf("-> Average CPU Cycles (1,000 runs): %" PRIu64 "\n\n", sum_cycles_test2 / TEST23_ITERATIONS);

    // ==========================================================
    // [테스트 3] 모든 비트가 0인 키 스칼라 곱셈 (1,000번 반복 평균)
    // ==========================================================
    printf("--- Test 3 : All Bits are 0 (1,000 Iterations) ---\n");
    uint64_t sum_cycles_test3 = 0;

    for (int iter = 0; iter < TEST23_ITERATIONS; iter++) {
        start_cycle = xthal_get_ccount();
        scalar_mult_jacobian(&J_Res, Private_Key_All_Zeros, &Base_Point);
        end_cycle = xthal_get_ccount();
        
        sum_cycles_test3 += (end_cycle - start_cycle);

        // 💡 워치독 방지
        vTaskDelay(1);
    }

    jacobian_to_affine(&A_Res, &J_Res);
    print_256("Result_3_X", A_Res.X);
    printf("-> Average CPU Cycles (1,000 runs): %" PRIu64 "\n\n", sum_cycles_test3 / TEST23_ITERATIONS);

    printf("-> ESP32 secp256k1 DaA Benchmark Complete!\n");

    // 태스크 자체 소멸 (안전한 종료)
    vTaskDelete(NULL);
}

// ====================================================================
// [ESP-IDF 메인 엔트리 포인트]
// ====================================================================
void app_main(void) {
    // 32KB 스택 할당 및 벤치마크 태스크 실행
    xTaskCreate(secp256k1_benchmark_task, "secp256k1_bench", 32768, NULL, 5, NULL);
}