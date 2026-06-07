// 몽고메리 사다리 알고리즘을 포함한 Curve25519 상의 ECC 코드.cpp 
// 리틀 엔디안 방식을 사용함. 

// ver 1.0 : 최초 완성본 (Step0~5 최초 병합 및 빌드 성공 확인)
// ver 1.1 : get_m0 함수를 제거하고 매직넘버 M0을 전역상수로 선언
// ver 1.2 : RFC7748 테스트벡터 통과 실험 중 발견한 오류 발견 및 수정 : 전역상수 A24_MONT[4] 를 기존 { 4623270, 0, 0, 0 } 에서 { 4623308, 0, 0, 0 }로 수정. 
// ver 1.3 : 결과 출력 창에 기준점으로 어떤 점을 사용하였는지도 출력하게끔 수정함.


#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <intrin.h> // MSVC 전용 내장 함수 헤더

// ====================================================================
// [Step0 : 전역 상수 정의]
// ====================================================================

const uint64_t P[4] = { 0xFFFFFFFFFFFFFFEDULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0x7FFFFFFFFFFFFFFFULL };               // P[4] : p = 2^{255} - 19 값. Fp 에서 소수 p의 그거 맞고, 이를 리틀엔디안 방식으로 저장한 것이다. (p = P[0] + P[1]*2^{64} + P[2]*2^{128} + P[3]*2^{192} 에서 P[0] = 2^{64}-19, P[1]=P[2]= 2^{64} -1, P[3]=2^{63}-1 )
const uint64_t R2_MOD_P[4] = { 0x5A4, 0, 0, 0 };                                                                                    // R2_MOD_P[4] : R^2 (mod p) where R = 2^{256}의 값. 일반 숫자를 연산 속도가 매우 빠른 몽고메리 형태(Montgomery Form) 로 변환할 때 사용되는 상수. to_montgomery 함수 내부에서 일반 숫자와 이 값을 몽고메리 곱셈하여 변환을 수행함. 
const uint64_t ONE[4] = {1, 0, 0, 0};                                                                                               // ONE[4] : Fp의 곱셈의 항등원 1. uint_t형 배열로 저장할 땐 1 = 1 + 0*2^{64} + 0*2^{128} + 0*2^{192}인 셈.
const uint64_t ONE_MONT[4] = { 38, 0, 0, 0 };                                                                                       // ONE_MONT[4]: 몽고메리 공간의 곱셈의 항등원 1*R (mod p). p = 2^{255} - 19 인걸 생각해보면 실제 값은 1*R (mod p) = 38이다. uint_t형 배열로 저장할 땐 38 = 38 + 0*2^{64} + 0*2^{128} + 0*2^{192}인 셈.
const uint64_t P_MINUS_2[4] = { 0xFFFFFFFFFFFFFFEBULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0x7FFFFFFFFFFFFFFFULL };       // P_MINUS_2[4] : p-2의 값. 즉,  2^{255} - 21. mod_inv_256 함수에서 a^{p-2} (mod p) 값 계산할 때 지수에 집어넣는 용도로 사용한다.
const uint64_t A24_MONT[4] = { 4623308, 0, 0, 0 };                                                                                  // A24_MONT[4] : 값 4623270. Curve25519:= y^2 = x^3 + 486662x^2 + x 에서 x^항의 계수 486662를 A라 하면, (A+2)/4 = . 이 Fp의 원소 (A+2)/4 = 를 몽고메리 형태로 변환하면 (A+2)/4 * R =  = . 이 값이 A24_MONT[4]이다. 
const uint64_t M0 = 9708812670373448219ULL;                                                                                         // p에 대한 매직넘버 M0 := -p^{-1} (mod R). 계산을 통해 얻는 실제 값은 M0 = 9708812670373448219.

typedef struct ProjPoint {                                                              // ProjPoint : ML알고리즘에선 아핀좌표계의 좌표(x, y)의 x좌표를 x := X/Z (mod p)의 형태로, 분리. 사영좌표계의 좌표 (X, Z)로 변환한 후. X와 Z를 각각 ϕ함수(몽고메리 영역변환함수, to_montgomery)에 넣어 몽고메리 형식 (ϕ(X), ϕ(Y)) = (XR mod p, YR mod p)으로 변환한다. 이 '몽고메리 형식으로 표현'된 사영공간상 좌표 (ϕ(X), ϕ(Y))를 구조체로 정의한다. 참고로 원래는 y도 따로 다뤄야하지만, ML의 논리상 y에 대한 논의가 없어도 x_Q를 구해낼 수 있다!
    uint64_t X[4];                                                                      // 사영좌표계의 점의 X좌표의 몽고메리 형식 ϕ(X). 256비트 정수
    uint64_t Z[4];                                                                      // 사영좌표계의 점의 Z좌표의 몽고메리 형식 ϕ(Z). 256비트 정수
} ProjPoint;

// ====================================================================
// [Step 1: 다중 정밀도 기본기]
// ====================================================================

uint64_t add_256(uint64_t* res, const uint64_t* a, const uint64_t* b) {                 // 상수시간 256비트 덧셈. 64비트 부호 없는 정수 a, b에 대해 res = a+b를 call by reference로, 올림수 carry 를 return으로 하여 나중에 사용함. 수학적으로 올바른 a+b의 값은 res와 carry를 조합해야 얻을 수 있음. 오버플로우 이슈로, 경우의 수가 좀 있다.
    uint64_t carry = 0;                                                                 // 일의자리.. 가 아니라 
    for (int i = 0; i < 4; i++) {
        uint64_t sum = a[i] + carry;                                                    // 이전 자릿수에서의 올림수 반영
        uint64_t carry1 = (sum < carry);                                                // 방금 연산으로 올림수가 발생carry를 더하다가 오버플로우가 났는지 확인 . 참고로 a = (b > c);라 하면 관계연산자 >의 결과(true 1, false 0)를 대입연산자 = 를 통해 a에 저장하는 것. 즉, b>c이면 a=1, b<=c이면 a=0.
        uint64_t total_sum = sum + b[i];
        uint64_t carry2 = (total_sum < b[i]);                                           // b[i]를 더하다가 오버플로우가 났는지 확인
        res[i] = total_sum;
        carry = carry1 | carry2;                                                        // 두 경우 중 한 번이라도 오버플로우가 났으면 carry는 1  (carry1과2의 비트 OR연산의 결과를 carry에 대입하는거다.)
    }
    return carry;                                                                       // 나중에 mod_add_256, mod_sub_256에서 
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

//  

void cswap_256(uint64_t swap_flag, uint64_t* a, uint64_t* b) {                      // 상수 시간 조건부 스왑(Constant-Time Conditional Swap) : 좌표의 성분을 서로 교환하는 함수. 부채널 공격 중 '타이밍 공격(Timing Attack)'을 방어하기 위해 if 문을 의도적으로 제거하고 오직 비트(Bit) 연산만으로 만들어졌다.
    uint64_t mask = 0 - swap_flag;                                                  // 매개변수 swap_flag는 0(안 바꿈) 또는 1(바꿈)입니다.
    for (int i = 0; i < 4; i++) {
        uint64_t dummy = mask & (a[i] ^ b[i]);                                      // 비트 XOR연산자 ^ : 숫자의 비트가 서로 다른 부분만 1로 표시하는 연산
        a[i] ^= dummy;
        b[i] ^= dummy;
    }
}

void cswap_point(uint64_t swap_flag, ProjPoint* p1, ProjPoint* p2) {                // 두 점의 좌표를 매개변수로 받아서 서로의 좌표를 교환하는 함수. 
    cswap_256(swap_flag, p1->X, p2->X);                                             // 위에서 정의한 cswap_256함수를 이용해 X좌표(256bit)끼리
    cswap_256(swap_flag, p1->Z, p2->Z);                                             // Z좌표끼리 스왑!  . cswap함수는 이 작업을 위해 정의한거다.
}

void copy_point(ProjPoint* dest, const ProjPoint* src) {                            // 좌표 복사 함수. 
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
    low = _umul128(a, b, &high);                                                         // _umul128 : 64비트 부호 없는 정수 두 개를 곱하여 128비트 결과를 생성. 상위 64비트는 포인터 인자로 반환하고, 하위 64비트는 함수 반환값으로 제공하여 오버플로 없이 큰 수의 곱셈을 효율적으로 처리. 여기서 상/하위란, 10진수 1234를 예로 들면 상위는 12, 하위는 34의 포지션.
    low += c;
    high += (low < c);
    low += d;                                                                            // 
    high += (low < d);
    *carry = high;
    return low;
}

void mod_add_256(uint64_t* res, const uint64_t* a, const uint64_t* b, const uint64_t* p) {     // 모듈러 덧셈 연산
    uint64_t sum[4], diff[4];                                                                  // 각각 오버플로우, 언더플로우가 발생했을 수 있는 덧셈, 뺄셈의 결과값(res) .. 을 저장해둘 지역변수
    uint64_t carry = add_256(sum, a, b);                                                       // 피연산자 : a, b add256의 return이 carry값.. 즉, 오퍼플로우 발생(1) 또는 안발생(0)임을 기억하자.
    uint64_t borrow = sub_256(diff, sum, p);                                                   // 피연산자 : sum, p where sum=a+b. sub256의 return이 carry값.. 즉, 언더플로우 발생(1) 또는 안발생(0)임을 기억하자.
    uint64_t mask = 0 - ((~carry & borrow) & 1);
    for (int i = 0; i < 4; i++) { res[i] = (sum[i] & mask) | (diff[i] & ~mask); }
}

void mod_sub_256(uint64_t* res, const uint64_t* a, const uint64_t* b, const uint64_t* p) {                  // 모듈러 뺄셈 연산
    uint64_t diff[4], sum[4];
    uint64_t borrow = sub_256(diff, a, b);
    add_256(sum, diff, p);
    uint64_t mask = 0 - borrow;
    for (int i = 0; i < 4; i++) {
        res[i] = (sum[i] & mask) | (diff[i] & ~mask);                                                        // P=Q|R : 
    }
}

void mont_mul_256(uint64_t* res, const uint64_t* a, const uint64_t* b, const uint64_t* p) {     // 몽고메리 곱셈 연산. 몽고메리 곱셈함수 MONT(X, Y) = XYR^{-1} (mod p) 을 구현한 것이다.. 모듈러 곱셈의 포지션. 
                                                                                                // 수학의 함수와 컴퓨터의 함수의 차이 : 컴퓨터의 함수는 매개변수가 수학으로 쳤을때 어떤 공간의 원소인지는 구분하는 능력이 없다! 그저 코딩된대로 계산할 뿐. 그리하여, 수학에선 엄밀히 구분하고 컴퓨터의 함수에선 인수로 무엇을 보내는지에 따라 코드는 같아도 역할이 달라질 수 있는거다! 예를 들어, a, b에 어떤 의도로 어떤 인수가 들어오느냐에 따라 mont_mul_256은 바
                                                                                                // a, b ∈ M     ⇒ mont_mul_256(res, a, b, p) : res에 MONT(a, b) = abR^{-1} (mod p) 값을 저장하는 함수. where MONT : M x M -> M           .. 의 역할과 기능을 한다! (계산 : input a, b ; output abR^{-1} (mod p) )                       -> 이걸 우리가 '몽고메리 곱셈'이라 부르는 연산이다!
                                                                                                // a∈M, b=R^2   ⇒ mont_mul_256(res, a, R2_MOD_P, p) : res에 ϕ(a) = aR (mod p) 값을 저장하는 함수. where ϕ: Fp -> M                       .. 의 역할과 기능을 한다! (계산 : input a, R^2 ; output aR^2R^{-1} (mod p) = aR (mod p) )      -> 이게 to_montgomery 함수이다!
                                                                                                // a∈M, b=1∈Fp ⇒ mont_mul_256(res, a, ONE, p) : res에 ϕ^{-1}(a) = MONT(a, 1) = aR^{-1} (mod p) 값을 저장하는 함수. where MONT = ϕ^{-1} : M -> Fp     .. 의 역할과 기능을 한다! (계신 : input a, 1 ; output aR^{-1} (mod p) )                       -> 이게 from_montgomery 함수이다!

    uint64_t t[5] = { 0 };
    uint64_t carry = 0;
    for (int i = 0; i < 4; i++) {                                                              // 여기서부터, T:=A'B'=ABR^2 '값' 계산. mod p에 대한게 아니라 그냥. 수학ver에선 
        carry = 0;
        for (int j = 0; j < 4; j++)                                 
            t[j] = mac_64(a[i], b[j], t[j], carry, &carry);                           
        t[4] += carry;                                                                           // 

        uint64_t m = t[0] * M0;                                                                  // 여기서부터 몽고메리 리덕션 과정. MONT에서 abR^{-1}을 구한, ϕ에서 aR을 구한, ϕ^{-1}에서 aR^{-1}을 계산한 후 mod p를 적용하는 단계임. 
        carry = 0;
        mac_64(m, p[0], t[0], carry, &carry);
        
        for (int j = 1; j < 4; j++)
            t[j - 1] = mac_64(m, p[j], t[j], carry, &carry);                                     // mac_64의 결과를 t[j]가 아니라 t[j-1]에 대입하는 것이 핵심! t가 64비트 정수형 배열이므로, 배열 한칸 아레에 저장한다는 것은 오른쪽으로 64비트 shift, 수학으로치면 R(=2^64)값으로 나누는 연산이다!
        
        uint64_t sum = t[4] + carry; t[3] = sum; t[4] = (sum < carry);
    }
    uint64_t temp[4];                                                                            // 여기서부터, 상수시간 조건부 보정 뺄셈. 수학으로 치면 [0, 2p-1]범위의 값을 갖게되는 몽고메리 곱셈의 결과에 mod p를 취하는 과정. 
    uint64_t borrow = sub_256(temp, t, p);
    uint64_t mask = 0 - borrow;
    for (int i = 0; i < 4; i++)
        res[i] = (t[i] & mask) | (temp[i] & ~mask);
}

void to_montgomery(uint64_t* res, const uint64_t* a, const uint64_t* p) {               // 동형사상 ϕ(x) = xR (mod p) 일반 숫자를 몽고메리 공간의 숫자로 변환                                                                 
    mont_mul_256(res, a, R2_MOD_P, p);                                                  // mont_mul_256(res, a, p) = res에 MONT(a, b) = abR^{-1} (mod p)임을 떠올려보면,  
                                                                                        // mont_mul_256(res, a, R2_MOD_P, p) = res에 MONT(a, R^2) = aR^2R^{-1} (mod p) = aR (mod p)를 저장하는 함수이다. 
}

void from_montgomery(uint64_t* res, const uint64_t* a_mont, const uint64_t* p) {       // 사상 ϕ^{-1}(X) = MONT(X, 1) =  XR^{-1} (mod p) 몽고메리 공간의 숫자를 일반 숫자로 변환
    mont_mul_256(res, a_mont, ONE, p);                                                 // 역시 mont_mul_256(res, a_mont, ONE, p) = res에 MONT(a, 1) = a*1*R^2 (mod p)
}

void mod_inv_256(uint64_t* res, const uint64_t* a, const uint64_t* p) {                              // mod_inv_256 : Fp의 어떤 원소 a에 대해 a^{-1} (mod p) 값을 구하는 함수. By 페르마의 소정리, a≠0, 1≤a≤p-1 ⇒ a^{p-2} ≡ a^{-1} (mod p)라는 사실을 이용해, 함수에서 실제로 계산해 얻어내는 값은 a^{p-2} (mod p)
                                                                                                     // 나눗셈 함수.. 무거운 연산! 그래서 Step4 몽고메리 사다리 알고리즘에선 사용하지 않고, Step5에서 affine x = XZ^{-1} (mod p) 를 계산하기위해  Z의 역원 Z^{-1}을 구해는데 딱 한번 쓰임.
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

void point_double(ProjPoint* res, const ProjPoint* pt, const uint64_t* p) {     // 타원곡선 점 두배 연산. 사영좌표계 점 R(X, Z)에 대해 2R((X^2-Z^2)^2, (4XZ)*{(X-Z)^2 + (A+2)/4 * 4XZ ) where A = 전역상수 A24_MONT 을 구하기. 
    uint64_t t0[4], t1[4], U[4], V[4], Diff[4], temp[4];
    mod_add_256(t0, pt->X, pt->Z, p); mod_sub_256(t1, pt->X, pt->Z, p);         // X+Z, X-Z 계산
    mont_mul_256(U, t0, t0, p); mont_mul_256(V, t1, t1, p);                     // U = (X+Z)^2, V = (X-Z)^2 계산
    mod_sub_256(Diff, U, V, p);                                                 // U-V (= 4XZ) 계산
    mont_mul_256(res->X, U, V, p);                                              // 최종 X좌표 X_2 = UV(=(X^2-Z^2)^2 )의 계산
    mont_mul_256(temp, Diff, A24_MONT, p); mod_add_256(temp, V, temp, p);       // 여기서부터, 최종 Z좌표 Z=(4XZ)*{(X-Z)^2 + (A+2)/4 * 4XZ } where A = 전역상수 A24_MONT 계산
    mont_mul_256(res->Z, Diff, temp, p);                                                
}

void point_add(ProjPoint* res, const ProjPoint* P1, const ProjPoint* P2, const ProjPoint* P_diff, const uint64_t* p) {          // 타원곡선 점 덧셈 연산. 사영좌표계 점 R0(X, Z) 과 R1(U, V)에 대해 (R0 + R1)( (기준점 P의 Z성분)*(UX-VZ)^2, (기준점 P이 X성분)(UZ-VX)^2) 계산. 몽고메리 사다리 알고리즘에 최적화된 차등 덧셈 방식을 사용중.
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
    k[0] &= 0xFFFFFFFFFFFFFFF8ULL; k[3] &= 0x7FFFFFFFFFFFFFFFULL; k[3] |= 0x4000000000000000ULL;        // 클램핑 : Curve25519의 꼬임보안의 핵심 중 두번째.. order 4에 대한 공격 방어.

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
    copy_point(res, &R0);                                                               // res에 R0을 저장 : 256번의 ML 알고리즘 반복 후의 점 R0이 바로 점 Q의 사영좌표 표현이다! 이 R0의 X성분과 Z성분으로 아핀좌표표현 x성분 x = XZ^{-1} (mod p) 값으로 변환하는건 이 함수에서 하진 않고 get_affine_x 함수에서 수행한다. 
}

// ====================================================================
// [Step 5: 최종 좌표 변환 (Affine X)] 
// ====================================================================
void get_affine_x(uint64_t* affine_x, const ProjPoint* pt_mont, const uint64_t* p) {
    uint64_t z_norm[4], z_inv[4];
    from_montgomery(z_norm, pt_mont->Z, p);
    mod_inv_256(z_inv, z_norm, p);                                  // Z^{-1} 값 구하기
    mont_mul_256(affine_x, pt_mont->X, z_inv, p);                   // x_Q <- XZ^{-1} mod p
}

// 결과 출력 편의 함수
void print_256(const char* name, const uint64_t* num) {
    printf("%s:\n%016llx %016llx %016llx %016llx\n", name, num[3], num[2], num[1], num[0]);
}

// ====================================================================
// [메인 함수] Curve25519 테스트
// ====================================================================
int main() {
    printf("=== Curve25519 Montgomery Ladder (Step 1 ~ 5) ===\n\n");

    // 1. 기준점 X (표준 값 9)                                                                          // Q=kP에서의 기준점 P ... 놀랍게도, gemini가 괜히 기준점을 9로 잡은게 아니었다. 
    uint64_t Base_X[4] = { 9, 0, 0, 0 };

    // 2. 사용자의 256비트 비밀키 (테스트용)
    uint64_t Private_Key[4] = { 0x123456789ABCDEF0ULL, 0x1111222233334444ULL, 0, 0 };                   // k값. 테스트용이다. 리틀엔디언 방식의 16진법 표현으로 000000000000000000000000000000001111222233334444123456789ABCDEF0 라는 값이라고 함.

    // 3. 스칼라 곱셈 실행 (Step 1 ~ Step 4 동작)
    ProjPoint Projective_Result;
    curve25519_scalar_mult(&Projective_Result, Private_Key, Base_X);

    // 4. 최종 좌표 복원 (Step 5 동작)
    uint64_t Final_Affine_X[4];                                                                         // Q=kP에서의 결과값 Q
    get_affine_x(Final_Affine_X, &Projective_Result, P);

    // 5. 완료.
    printf("[Base Point (Affine X Coordinate)]\n");
    print_256("Base_X", Base_X);
    printf("\n");

    printf("[Final Shared Secret (Affine X Coordinate)]\n");
    print_256("Result", Final_Affine_X);
    printf("\n-> Project Complete! The ECC core is fully functional.\n");

    return 0;
}
