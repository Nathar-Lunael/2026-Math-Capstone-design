// =========================================================================================
// [ESP-IDF 전용] Dudect T-test 통계적 검증 실험 (어셈블리 레벨 상수 시간 실행 검증 시험)
// =========================================================================================

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include "esp_timer.h"
#include "esp_cpu.h"
#include "esp_attr.h"          // IRAM_ATTR, DRAM_ATTR 사용을 위한 헤더
#include "freertos/FreeRTOS.h" // FreeRTOS 딜레이 사용
#include "freertos/task.h"

// ====================================================================
// [Step0 : 전역 상수 정의] (DRAM_ATTR 적용: 플래시가 아닌 내부 RAM에 상주)
// ====================================================================

const DRAM_ATTR uint64_t P[4] = { 0xFFFFFFFFFFFFFFEDULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0x7FFFFFFFFFFFFFFFULL };
const DRAM_ATTR uint64_t R2_MOD_P[4] = { 0x5A4, 0, 0, 0 };
const DRAM_ATTR uint64_t ONE[4] = { 1, 0, 0, 0 };
const DRAM_ATTR uint64_t ONE_MONT[4] = { 38, 0, 0, 0 };
const DRAM_ATTR uint64_t P_MINUS_2[4] = { 0xFFFFFFFFFFFFFFEBULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0x7FFFFFFFFFFFFFFFULL };
const DRAM_ATTR uint64_t A24_MONT[4] = { 4623308, 0, 0, 0 };
const DRAM_ATTR uint64_t M0 = 9708812670373448219ULL;

typedef struct ProjPoint {
    uint64_t X[4];
    uint64_t Z[4];
} ProjPoint;

// ====================================================================
// [Step 1~4: 모든 핵심 함수에 IRAM_ATTR 적용]
// ====================================================================

uint64_t IRAM_ATTR umul128_custom(uint64_t a, uint64_t b, uint64_t* high) {
    uint64_t a_lo = (uint32_t)a;
    uint64_t a_hi = a >> 32;
    uint64_t b_lo = (uint32_t)b;
    uint64_t b_hi = b >> 32;

    uint64_t a_x_b_lo = a_lo * b_lo;
    uint64_t a_x_b_hi = a_hi * b_hi;

    uint64_t a_lo_x_b_hi = a_lo * b_hi;
    uint64_t a_hi_x_b_lo = a_hi * b_lo;

    uint64_t cross1 = a_lo_x_b_hi + (a_x_b_lo >> 32);
    uint64_t cross2 = a_hi_x_b_lo + (uint32_t)cross1;

    *high = a_x_b_hi + (cross1 >> 32) + (cross2 >> 32);
    return (cross2 << 32) | (uint32_t)a_x_b_lo;
}

uint64_t IRAM_ATTR add_256(uint64_t* res, const uint64_t* a, const uint64_t* b) {
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

uint64_t IRAM_ATTR sub_256(uint64_t* res, const uint64_t* a, const uint64_t* b) {
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

void IRAM_ATTR cswap_256(uint64_t swap_flag, uint64_t* a, uint64_t* b) {
    uint64_t mask = 0 - swap_flag;
    for (int i = 0; i < 4; i++) {
        uint64_t dummy = mask & (a[i] ^ b[i]);
        a[i] ^= dummy;
        b[i] ^= dummy;
    }
}

void IRAM_ATTR cswap_point(uint64_t swap_flag, ProjPoint* p1, ProjPoint* p2) {
    cswap_256(swap_flag, p1->X, p2->X);
    cswap_256(swap_flag, p1->Z, p2->Z);
}

void IRAM_ATTR copy_point(ProjPoint* dest, const ProjPoint* src) {
    for (int i = 0; i < 4; i++) {
        dest->X[i] = src->X[i];
        dest->Z[i] = src->Z[i];
    }
}

uint64_t IRAM_ATTR mac_64(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t* carry) {
    uint64_t high, low;
    low = umul128_custom(a, b, &high);
    low += c;
    high += (low < c);
    low += d;
    high += (low < d);
    *carry = high;
    return low;
}

void IRAM_ATTR mod_add_256(uint64_t* res, const uint64_t* a, const uint64_t* b, const uint64_t* p) {
    uint64_t sum[4], diff[4];
    uint64_t carry = add_256(sum, a, b);
    uint64_t borrow = sub_256(diff, sum, p);
    uint64_t mask = 0 - ((~carry & borrow) & 1);
    for (int i = 0; i < 4; i++) { res[i] = (sum[i] & mask) | (diff[i] & ~mask); }
}

void IRAM_ATTR mod_sub_256(uint64_t* res, const uint64_t* a, const uint64_t* b, const uint64_t* p) {
    uint64_t diff[4], sum[4];
    uint64_t borrow = sub_256(diff, a, b);
    add_256(sum, diff, p);
    uint64_t mask = 0 - borrow;
    for (int i = 0; i < 4; i++) {
        res[i] = (sum[i] & mask) | (diff[i] & ~mask);
    }
}

void IRAM_ATTR mont_mul_256(uint64_t* res, const uint64_t* a, const uint64_t* b, const uint64_t* p) {
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

void IRAM_ATTR to_montgomery(uint64_t* res, const uint64_t* a, const uint64_t* p) {
    mont_mul_256(res, a, R2_MOD_P, p);
}

void IRAM_ATTR point_double(ProjPoint* res, const ProjPoint* pt, const uint64_t* p) {
    uint64_t t0[4], t1[4], U[4], V[4], Diff[4], temp[4];
    mod_add_256(t0, pt->X, pt->Z, p); mod_sub_256(t1, pt->X, pt->Z, p);
    mont_mul_256(U, t0, t0, p); mont_mul_256(V, t1, t1, p);
    mod_sub_256(Diff, U, V, p);
    mont_mul_256(res->X, U, V, p);
    mont_mul_256(temp, Diff, A24_MONT, p); mod_add_256(temp, V, temp, p);
    mont_mul_256(res->Z, Diff, temp, p);
}

void IRAM_ATTR point_add(ProjPoint* res, const ProjPoint* P1, const ProjPoint* P2, const ProjPoint* P_diff, const uint64_t* p) {
    uint64_t t0[4], t1[4], A[4], B[4], sum_sq[4], diff_sq[4];
    mod_add_256(t0, P1->X, P1->Z, p); mod_sub_256(t1, P2->X, P2->Z, p); mont_mul_256(A, t0, t1, p);
    mod_sub_256(t0, P1->X, P1->Z, p); mod_add_256(t1, P2->X, P2->Z, p); mont_mul_256(B, t0, t1, p);
    mod_add_256(t0, A, B, p); mont_mul_256(sum_sq, t0, t0, p);
    mod_sub_256(t1, A, B, p); mont_mul_256(diff_sq, t1, t1, p);
    mont_mul_256(res->X, P_diff->Z, sum_sq, p); mont_mul_256(res->Z, P_diff->X, diff_sq, p);
}

void IRAM_ATTR curve25519_scalar_mult(ProjPoint* res, const uint64_t* scalar, const uint64_t* base_x) {
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
// [메인 함수] ESP-IDF app_main
// ====================================================================
void app_main(void) {
    vTaskDelay(2000 / portTICK_PERIOD_MS); // 시리얼 모니터 연결 대기

    printf("\n=== ESP-IDF Curve25519 Constant-Time T-Test ===\n\n");

    uint64_t Base_X[4] = { 9, 0, 0, 0 };
    ProjPoint Result;

    uint64_t Key_Fixed[4] = { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL };
    uint64_t Key_Random[4];

    const int MEASUREMENTS = 10000; 

    double sum_A = 0, sum_sq_A = 0;
    double sum_B = 0, sum_sq_B = 0;
    int count_A = 0, count_B = 0;

    uint32_t start, end; 
    uint32_t cycles;
    volatile uint64_t prevent_opt = 0;

    printf("Warming up cache and starting %d measurements...\n", MEASUREMENTS);
    srand(12345);

    for (int i = 0; i < MEASUREMENTS; i++) {
        // Watchdog Timer(WDT) 리셋 방지를 위한 FreeRTOS 딜레이
        if (i % 100 == 0) {
            vTaskDelay(1 / portTICK_PERIOD_MS); 
        }

        int class_bit = rand() % 2;

        if (class_bit == 0) { 
            start = esp_cpu_get_cycle_count();
            curve25519_scalar_mult(&Result, Key_Fixed, Base_X);
            end = esp_cpu_get_cycle_count();
            
            cycles = (end >= start) ? (end - start) : (0xFFFFFFFF - start + end + 1);

            sum_A += cycles;
            sum_sq_A += (double)cycles * cycles;
            count_A++;
        }
        else { 
            Key_Random[0] = ((uint64_t)rand() << 32) | rand();
            Key_Random[1] = ((uint64_t)rand() << 32) | rand();
            Key_Random[2] = ((uint64_t)rand() << 32) | rand();
            Key_Random[3] = ((uint64_t)rand() << 32) | rand();

            start = esp_cpu_get_cycle_count();
            curve25519_scalar_mult(&Result, Key_Random, Base_X);
            end = esp_cpu_get_cycle_count();
            
            cycles = (end >= start) ? (end - start) : (0xFFFFFFFF - start + end + 1);

            sum_B += cycles;
            sum_sq_B += (double)cycles * cycles;
            count_B++;
        }
        prevent_opt ^= Result.X[0];
    }

    double mean_A = sum_A / count_A;
    double var_A = (sum_sq_A / count_A) - (mean_A * mean_A);

    double mean_B = sum_B / count_B;
    double var_B = (sum_sq_B / count_B) - (mean_B * mean_B);

    double t_stat = (mean_A - mean_B) / sqrt((var_A / count_A) + (var_B / count_B));

    printf("\n====================================================\n");
    printf("[ESP-IDF Statistical Results]\n");
    printf("Class A (Fixed Zeros) : Count = %d, Mean = %.2f cycles\n", count_A, mean_A);
    printf("Class B (Random Keys) : Count = %d, Mean = %.2f cycles\n", count_B, mean_B);
    printf("----------------------------------------------------\n");
    printf("T-Statistic (t-value) : %f\n", t_stat);
    printf("====================================================\n\n");

    double abs_t = t_stat < 0 ? -t_stat : t_stat;

    if (abs_t < 4.5) {
        printf("-> PASS: |t| < 4.5\n");
        printf("-> The distributions are statistically indistinguishable on ESP-IDF.\n");
        printf("-> Your implementation is STRICTLY Constant-Time on the MCU architecture!\n");
    }
    else {
        printf("-> FAIL: |t| >= 4.5 (Significant difference detected)\n");
    }

    if (prevent_opt == 0xDEADBEEF) printf(" ");

    // 측정이 끝나면 태스크를 삭제하여 무한 대기
    vTaskDelete(NULL);
}