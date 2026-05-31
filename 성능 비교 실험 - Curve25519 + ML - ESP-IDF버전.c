// [제 1실험] Curve25519 + ML ECC - ESP-IDF Framework 버전
// (워치독 방지 최적화 반영 완료)

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>

// FreeRTOS 및 ESP32 커널 헤더 포함
#include "freertos/FreeRTOS.h"
#include "freertos/task.h"
#include "xtensa/core-macros.h" // ESP32 하드웨어 사이클 카운터(xthal_get_ccount) 사용

#include "test_vectors.h"  // 파이썬으로 자동 생성한 100개의 무작위 스칼라 헤더 포함

// ====================================================================
// [Step0 : 전역 상수 정의]
// ====================================================================

const uint64_t P[4] = { 0xFFFFFFFFFFFFFFEDULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0x7FFFFFFFFFFFFFFFULL };
const uint64_t R2_MOD_P[4] = { 0x5A4, 0, 0, 0 };
const uint64_t ONE[4] = { 1, 0, 0, 0 };
const uint64_t ONE_MONT[4] = { 38, 0, 0, 0 };
const uint64_t P_MINUS_2[4] = { 0xFFFFFFFFFFFFFFEBULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0x7FFFFFFFFFFFFFFFULL };
const uint64_t A24_MONT[4] = { 4623308, 0, 0, 0 };
const uint64_t M0 = 9708812670373448219ULL;

typedef struct ProjPoint {                                                              
    uint64_t X[4];
    uint64_t Z[4];
} ProjPoint;

// ====================================================================
// [Step 1: 다중 정밀도 기본기 및 임베디드용 고속 연산]
// ====================================================================

// ESP32 아키텍처 호환을 위한 크로스 컴파일용 128비트 부호 없는 정수 곱셈 함수 (MSVC _umul128 대체)
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

uint64_t add_256(uint64_t* res, const uint64_t* a, const uint64_t* b) {
    uint64_t carry = 0;                                                                 
    for (int i = 0; i < 4; i++) {
        uint64_t sum = a[i] + carry;                                                    
        uint64_t carry1 = (sum < carry);                                                
        uint64_t total_sum = sum + b[i];
        uint64_t carry2 = (total_sum < b[i]);                                            
        res[i] = total_sum;
        carry = carry1 | carry2;                                                        
    }
    return carry;                                                                      
}

uint64_t sub_256(uint64_t* res, const uint64_t* a, const uint64_t* b) {
    uint64_t borrow = 0;
    for (int i = 0; i < 4; i++) {
        uint64_t diff = a[i] - borrow;
        uint64_t borrow1 = (a[i] < borrow);
        uint64_t total_diff = diff - b[i];
        uint64_t borrow2 = (diff < b[i]);
        res[i] = total_diff;
        borrow = borrow1 | borrow2;
    }
    return borrow;
}

void cswap_256(uint64_t swap_flag, uint64_t* a, uint64_t* b) {
    uint64_t mask = 0 - swap_flag;                                                  
    for (int i = 0; i < 4; i++) {
        uint64_t dummy = mask & (a[i] ^ b[i]);                                    
        a[i] ^= dummy;
        b[i] ^= dummy;
    }
}

void cswap_point(uint64_t swap_flag, ProjPoint* p1, ProjPoint* p2) {
    cswap_256(swap_flag, p1->X, p2->X);                                             
    cswap_256(swap_flag, p1->Z, p2->Z);                                             
}

void copy_point(ProjPoint* dest, const ProjPoint* src) {
    for (int i = 0; i < 4; i++) {
        dest->X[i] = src->X[i];
        dest->Z[i] = src->Z[i];
    }
}

// ====================================================================
// [Step 2: 고속 모듈러 연산 세트]
// ====================================================================

uint64_t mac_64(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t* carry) {
    uint64_t high, low;
    low = mock_umul128(a, b, &high);                                                         
    low += c;
    high += (low < c);
    low += d;                                                                            
    high += (low < d);
    *carry = high;
    return low;
}

void mod_add_256(uint64_t* res, const uint64_t* a, const uint64_t* b, const uint64_t* p) {
    uint64_t sum[4], diff[4];                                                                 
    uint64_t carry = add_256(sum, a, b);                                                       
    uint64_t borrow = sub_256(diff, sum, p);                                                   
    uint64_t mask = 0 - ((~carry & borrow) & 1);
    for (int i = 0; i < 4; i++) { res[i] = (sum[i] & mask) | (diff[i] & ~mask); }
}

void mod_sub_256(uint64_t* res, const uint64_t* a, const uint64_t* b, const uint64_t* p) {
    uint64_t diff[4], sum[4];
    uint64_t borrow = sub_256(diff, a, b);
    add_256(sum, diff, p);
    uint64_t mask = 0 - borrow;
    for (int i = 0; i < 4; i++) {
        res[i] = (sum[i] & mask) | (diff[i] & ~mask);                                                                        
    }
}

void mont_mul_256(uint64_t* res, const uint64_t* a, const uint64_t* b, const uint64_t* p) {
    uint64_t t[5] = { 0 };
    uint64_t carry = 0;
    for (int i = 0; i < 4; i++) {                                                                           
        carry = 0;
        for (int j = 0; j < 4; j++)
            t[j] = mac_64(a[i], b[j], t[j], carry, &carry);
        t[4] += carry;                                                                             

        uint64_t m = t[0] * M0;                                                                  
        carry = 0;
        mac_64(m, p[0], t[0], carry, &carry);

        for (int j = 1; j < 4; j++)
            t[j - 1] = mac_64(m, p[j], t[j], carry, &carry);                                     

        uint64_t sum = t[4] + carry; t[3] = sum; t[4] = (sum < carry);
    }
    uint64_t temp[4];                                                                            
    uint64_t borrow = sub_256(temp, t, p);
    uint64_t mask = 0 - borrow;
    for (int i = 0; i < 4; i++)
        res[i] = (t[i] & mask) | (temp[i] & ~mask);
}

void to_montgomery(uint64_t* res, const uint64_t* a, const uint64_t* p) {                               
    mont_mul_256(res, a, R2_MOD_P, p);                                                  
}

void from_montgomery(uint64_t* res, const uint64_t* a_mont, const uint64_t* p) {       
    mont_mul_256(res, a_mont, ONE, p);                                                 
}

void mod_inv_256(uint64_t* res, const uint64_t* a, const uint64_t* p) {                               
    uint64_t a_mont[4], r_mont[4] = { ONE_MONT[0], ONE_MONT[1], ONE_MONT[2], ONE_MONT[3] };
    to_montgomery(a_mont, a, p);
    for (int i = 254; i >= 0; i--) {
        mont_mul_256(r_mont, r_mont, r_mont, p);
        if (((P_MINUS_2[i / 64] >> (i % 64)) & 1) == 1) mont_mul_256(r_mont, r_mont, a_mont, p);
    }
    from_montgomery(res, r_mont, p);
}

// ====================================================================
// [Step 3: 타원곡선 점 연산]
// ====================================================================

void point_double(ProjPoint* res, const ProjPoint* pt, const uint64_t* p) {     
    uint64_t t0[4], t1[4], U[4], V[4], Diff[4], temp[4];
    mod_add_256(t0, pt->X, pt->Z, p); mod_sub_256(t1, pt->X, pt->Z, p);         
    mont_mul_256(U, t0, t0, p); mont_mul_256(V, t1, t1, p);                     
    mod_sub_256(Diff, U, V, p);                                                 
    mont_mul_256(res->X, U, V, p);                                              
    mont_mul_256(temp, Diff, A24_MONT, p); mod_add_256(temp, V, temp, p);       
    mont_mul_256(res->Z, Diff, temp, p);
}

void point_add(ProjPoint* res, const ProjPoint* P1, const ProjPoint* P2, const ProjPoint* P_diff, const uint64_t* p) {          
    uint64_t t0[4], t1[4], A[4], B[4], sum_sq[4], diff_sq[4];
    mod_add_256(t0, P1->X, P1->Z, p); mod_sub_256(t1, P2->X, P2->Z, p); mont_mul_256(A, t0, t1, p);
    mod_sub_256(t0, P1->X, P1->Z, p); mod_add_256(t1, P2->X, P2->Z, p); mont_mul_256(B, t0, t1, p);
    mod_add_256(t0, A, B, p); mont_mul_256(sum_sq, t0, t0, p);
    mod_sub_256(t1, A, B, p); mont_mul_256(diff_sq, t1, t1, p);
    mont_mul_256(res->X, P_diff->Z, sum_sq, p); mont_mul_256(res->Z, P_diff->X, diff_sq, p);
}

// ====================================================================
// [Step 4: 몽고메리 사다리 메인 루프] 
// ====================================================================
void curve25519_scalar_mult(ProjPoint* res, const uint64_t* scalar, const uint64_t* base_x) {
    ProjPoint R0, R1;
    uint64_t k[4];

    for (int i = 0; i < 4; i++)
        k[i] = scalar[i];
    k[0] &= 0xFFFFFFFFFFFFFFF8ULL; k[3] &= 0x7FFFFFFFFFFFFFFFULL; k[3] |= 0x4000000000000000ULL;     

    uint64_t zero[4] = { 0, 0, 0, 0 };
    to_montgomery(R0.X, ONE, P); to_montgomery(R0.Z, zero, P);
    to_montgomery(R1.X, base_x, P); to_montgomery(R1.Z, ONE, P);

    ProjPoint P_diff; copy_point(&P_diff, &R1);

    uint64_t prev_bit = 0;
    for (int i = 254; i >= 0; i--) {
        uint64_t bit = (k[i / 64] >> (i % 64)) & 1;
        uint64_t swap = bit ^ prev_bit;
        cswap_point(swap, &R0, &R1);

        ProjPoint next_R0, next_R1;
        point_add(&next_R1, &R0, &R1, &P_diff, P);
        point_double(&next_R0, &R0, P);

        copy_point(&R0, &next_R0); copy_point(&R1, &next_R1);
        prev_bit = bit;
    }
    cswap_point(prev_bit, &R0, &R1);
    copy_point(res, &R0);                                                                           
}

// ====================================================================
// [Step 5: 최종 좌표 변환 (Affine X)] 
// ====================================================================
void get_affine_x(uint64_t* affine_x, const ProjPoint* pt_mont, const uint64_t* p) {
    uint64_t z_norm[4], z_inv[4];
    from_montgomery(z_norm, pt_mont->Z, p);
    mod_inv_256(z_inv, z_norm, p);                                  
    mont_mul_256(affine_x, pt_mont->X, z_inv, p);                   
}

void print_256(const char* name, const uint64_t* num) {
    printf("%s:\n%016" PRIx64 " %016" PRIx64 " %016" PRIx64 " %016" PRIx64 "\n", name, num[3], num[2], num[1], num[0]);
}

// ====================================================================
// [FreeRTOS 벤치마킹 실행 Task]
// ====================================================================
void curve25519_benchmark_task(void* pvParameters) {
    printf("=== Curve25519 Montgomery Ladder ESP-IDF Start ===\n\n");

    uint64_t Base_X[4] = { 9, 0, 0, 0 };
    uint64_t Private_Key_All_Ones[4] = { 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL };
    uint64_t Private_Key_All_Zeros[4] = { 0, 0, 0, 0 };

    ProjPoint Projective_Result;
    uint64_t Final_Affine_X[4];

    // 하드웨어 타이머 레지스터값(32비트)을 담을 변수 선언
    uint32_t start_cycle, end_cycle;
    
    const int TEST1_SCALARS = 100;
    const int TEST23_ITERATIONS = 1000;

    print_256("Base_X", Base_X);
    printf("\n");

    // ==========================================================
    // [테스트 1] 무작위 256비트 스칼라 100개 평균 측정
    // ==========================================================
    printf("--- Test 1 : 100 Random Scalars (Average Case) ---\n");
    uint64_t sum_cycles_test1 = 0;

    for (int iter = 0; iter < TEST1_SCALARS; iter++) {
        start_cycle = xthal_get_ccount(); // Xtensa 하드웨어 사이클 카운터 시작
        curve25519_scalar_mult(&Projective_Result, Random_Scalars[iter], Base_X);
        end_cycle = xthal_get_ccount();   // 카운터 종료
        
        sum_cycles_test1 += (end_cycle - start_cycle);

        // 💡 핵심 최적화: OS에게 잠시 CPU를 양보하여 워치독 리셋 방지 (사이클 측정과는 무관하게 동작)
        vTaskDelay(1);
    }

    get_affine_x(Final_Affine_X, &Projective_Result, P);
    print_256("Result_1 (100th run)", Final_Affine_X);
    printf("-> Average CPU Cycles (100 runs): %" PRIu64 "\n\n", sum_cycles_test1 / TEST1_SCALARS);

    // ==========================================================
    // [테스트 2] 모든 비트가 1인 키 스칼라 곱셈 (1,000번 반복 평균)
    // ==========================================================
    printf("--- Test 2 : All Bits are 1 (1,000 Iterations) ---\n");
    uint64_t sum_cycles_test2 = 0;

    for (int iter = 0; iter < TEST23_ITERATIONS; iter++) {
        start_cycle = xthal_get_ccount();
        curve25519_scalar_mult(&Projective_Result, Private_Key_All_Ones, Base_X);
        end_cycle = xthal_get_ccount();
        
        sum_cycles_test2 += (end_cycle - start_cycle);

        // 💡 워치독 방지
        vTaskDelay(1);
    }

    get_affine_x(Final_Affine_X, &Projective_Result, P);
    print_256("Result_2", Final_Affine_X);
    printf("-> Average CPU Cycles (1,000 runs): %" PRIu64 "\n\n", sum_cycles_test2 / TEST23_ITERATIONS);

    // ==========================================================
    // [테스트 3] 모든 비트가 0인 키 스칼라 곱셈 (1,000번 반복 평균)
    // ==========================================================
    printf("--- Test 3 : All Bits are 0 (1,000 Iterations) ---\n");
    uint64_t sum_cycles_test3 = 0;

    for (int iter = 0; iter < TEST23_ITERATIONS; iter++) {
        start_cycle = xthal_get_ccount();
        curve25519_scalar_mult(&Projective_Result, Private_Key_All_Zeros, Base_X);
        end_cycle = xthal_get_ccount();
        
        sum_cycles_test3 += (end_cycle - start_cycle);

        // 💡 워치독 방지
        vTaskDelay(1);
    }

    get_affine_x(Final_Affine_X, &Projective_Result, P);
    print_256("Result_3", Final_Affine_X);
    printf("-> Average CPU Cycles (1,000 runs): %" PRIu64 "\n\n", sum_cycles_test3 / TEST23_ITERATIONS);

    printf("-> ESP32 Curve25519 Core Benchmark Complete!\n");

    // FreeRTOS Task는 작업 종료 시 반드시 스스로를 삭제해야 함
    vTaskDelete(NULL);
}

// ====================================================================
// [ESP-IDF 메인 엔트리 포인트]
// ====================================================================
void app_main(void) {
    // 대량의 벤치마크 루프와 변수를 다루기 위해 스택 크기를 32KB로 여유롭게 설정
    xTaskCreate(curve25519_benchmark_task, "curve25519_bench", 32768, NULL, 5, NULL);
}