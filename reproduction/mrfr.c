#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <time.h>

// 定义向量结构体
typedef struct {
    mpz_t x;//高精度整数
    mpz_t y;
} Vector;

// 初始化向量
void vector_init(Vector *v) {
    mpz_init(v->x);
    mpz_init(v->y);
}

// 清除向量
void vector_clear(Vector *v) {
    mpz_clear(v->x);
    mpz_clear(v->y);
}

// 打印向量
void vector_print(const Vector *v) {
    gmp_printf("(%Zd, %Zd)", v->x, v->y);
}

// 比较两个向量的∞范数
int vector_cmp(const Vector *a, const Vector *b) {
    mpz_t max_a, max_b;
    mpz_init(max_a);
    mpz_init(max_b);
    
    mpz_abs(max_a, a->x);//max_a=|a->x|
    if (mpz_cmp(max_a, a->y) < 0) {//若max_a>a->y 则返回1，若==，则0,若<,则-1
        mpz_abs(max_a, a->y);
    }
    
    mpz_abs(max_b, b->x);
    if (mpz_cmp(max_b, b->y) < 0) {
        mpz_abs(max_b, b->y);
    }
    
    int result = mpz_cmp(max_a, max_b);
    
    mpz_clear(max_a);
    mpz_clear(max_b);
    
    return result;
}

// (Improved RA)
void improved_ra(const char *sequence, int n, mpz_t p, mpz_t q) {
    // 初始化变量
    Vector alpha, beta;
    mpz_t z11, z12, z21, z22;
    mpz_t tmp1, tmp2, tmp3;
    
    vector_init(&alpha);
    vector_init(&beta);
    mpz_init(z11); mpz_init(z12); mpz_init(z21); mpz_init(z22);
    mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3);
    
    // 找到第一个非零位
    int k = 0;
    while (k < n && sequence[k] == '0') {
        k++;
    }
    
    // 初始化基础情况
    if (k <= 1) {
        mpz_set_ui(alpha.x, 1);//用一个unsigned long范围的正整数赋值给大整数alpha
        mpz_set_ui(alpha.y, 1);
        mpz_set_si(beta.x, -1);
        mpz_set_ui(beta.y, 1);
        
        mpz_set_ui(z11, 1);
        mpz_set_ui(z12, 0);
        mpz_set_ui(z21, 1);
        mpz_set_si(z22, -1);
    } else {
        // k>=2的情况
        mpz_set_ui(alpha.x, 0);
        mpz_set_ui(alpha.y, 2);
        
        mpz_ui_pow_ui(beta.x, 2, k-1);
        mpz_set_ui(beta.y, 1);
        
        mpz_set_ui(z11, 2);
        mpz_set_si(z12, -1);
        mpz_set_ui(z21, 1);
        mpz_set_ui(z22, 0);
    }
    
    // 处理剩余的位
    for (int i = k; i < n; i++) {
        int bit = sequence[i] - '0';
        
        // 计算条件
        mpz_mul_si(tmp1, z11, -bit);//tmp1=-ak*z11
        mpz_add(tmp1, tmp1, z12);
        mpz_mod_ui(tmp1, tmp1, 2);
        int cond1 = mpz_cmp_ui(tmp1, 0) == 0;
        
        mpz_mul_si(tmp1, z21, -bit);
        mpz_add(tmp1, tmp1, z22);
        mpz_mod_ui(tmp1, tmp1, 2);
        int cond2 = mpz_cmp_ui(tmp1, 0) == 0;
        
        if (cond1) {
            // 情况1: -ak*z11 + z12 ≡ 0 mod 2
            mpz_mul_si(tmp1, z11, -bit);
            mpz_add(tmp1, tmp1, z12);
            mpz_div_2exp(z12, tmp1, 1);
            
            // 计算2*beta
            Vector two_beta;
            vector_init(&two_beta);
            mpz_mul_2exp(two_beta.x, beta.x, 1);
            mpz_mul_2exp(two_beta.y, beta.y, 1);
            
            // 计算alpha - 2*beta
            Vector alpha_minus_2beta;
            vector_init(&alpha_minus_2beta);
            mpz_sub(alpha_minus_2beta.x, alpha.x, two_beta.x);
            mpz_sub(alpha_minus_2beta.y, alpha.y, two_beta.y);
            
            // 计算alpha + 2*beta
            Vector alpha_plus_2beta;
            vector_init(&alpha_plus_2beta);
            mpz_add(alpha_plus_2beta.x, alpha.x, two_beta.x);
            mpz_add(alpha_plus_2beta.y, alpha.y, two_beta.y);
            
            // 比较三种情况
            int cmp1 = vector_cmp(&two_beta, &alpha_minus_2beta);
            int cmp2 = vector_cmp(&two_beta, &alpha_plus_2beta);
            
            if (cmp1 <= 0 && cmp2 <= 0) {
                // 2*beta最小
                mpz_set(beta.x, two_beta.x);
                mpz_set(beta.y, two_beta.y);
                
                mpz_mul_si(tmp1, z21, -bit);
                mpz_add(z22, tmp1, z22);
                mpz_mul_2exp(z21, z21, 1);//
            } else {
                int cmp3 = vector_cmp(&alpha_minus_2beta, &alpha_plus_2beta);
                if (cmp3 <= 0) {
                    // alpha - 2*beta最小
                    mpz_set(beta.x, alpha_minus_2beta.x);
                    mpz_set(beta.y, alpha_minus_2beta.y);
                    
                    mpz_mul_si(tmp1, z21, bit);
                    mpz_add(tmp1, z12, tmp1);
                    mpz_sub(z22, tmp1, z22);
                    
                    mpz_mul_2exp(tmp1, z21, 1);
                    mpz_sub(z21, z11, tmp1);
                } else {
                    // alpha + 2*beta最小
                    mpz_set(beta.x, alpha_plus_2beta.x);
                    mpz_set(beta.y, alpha_plus_2beta.y);
                    
                    mpz_mul_si(tmp1, z21, bit);
                    mpz_sub(tmp1, z12, tmp1);
                    mpz_sub(z22, tmp1, z22);
                    
                    mpz_mul_2exp(tmp1, z21, 1);
                    mpz_add(z21, z11, tmp1);
                }
            }
            
            vector_clear(&two_beta);
            vector_clear(&alpha_minus_2beta);
            vector_clear(&alpha_plus_2beta);
        } else if (cond2) {
            // 情况2: -ak*z11 + z12 ≡ 1 mod 2 且 -ak*z21 + z22 ≡ 0 mod 2
            // 计算2*alpha
            Vector two_alpha;
            vector_init(&two_alpha);
            mpz_mul_2exp(two_alpha.x, alpha.x, 1);
            mpz_mul_2exp(two_alpha.y, alpha.y, 1);
            
            // 计算2*alpha - beta
            Vector two_alpha_minus_beta;
            vector_init(&two_alpha_minus_beta);
            mpz_sub(two_alpha_minus_beta.x, two_alpha.x, beta.x);
            mpz_sub(two_alpha_minus_beta.y, two_alpha.y, beta.y);
            
            // 计算2*alpha + beta
            Vector two_alpha_plus_beta;
            vector_init(&two_alpha_plus_beta);
            mpz_add(two_alpha_plus_beta.x, two_alpha.x, beta.x);
            mpz_add(two_alpha_plus_beta.y, two_alpha.y, beta.y);
            
            // 找出最小的两个向量
            Vector candidates[4] = {two_alpha, beta, two_alpha_minus_beta, two_alpha_plus_beta};
            int min1 = 0, min2 = 1;
            
            for (int j = 2; j < 4; j++) {
                if (vector_cmp(&candidates[j], &candidates[min1]) < 0) {
                    min2 = min1;
                    min1 = j;
                } else if (vector_cmp(&candidates[j], &candidates[min2]) < 0) {
                    min2 = j;
                }
            }
            
            // 更新alpha和beta
            Vector new_alpha = candidates[min1];
            Vector new_beta = candidates[min2];
            
            if (min1 == 0) {
                // 2*alpha最小
                mpz_set(alpha.x, two_alpha.x);
                mpz_set(alpha.y, two_alpha.y);
                
                mpz_mul_si(tmp1, z11, -bit);
                mpz_set(z12, tmp1);
                mpz_add(z12, z12, z12);
                mpz_mul_2exp(z11, z11, 1);
                
                mpz_mul_si(tmp1, z21, -bit);
                mpz_add(tmp1, tmp1, z22);
                mpz_div_2exp(z22, tmp1, 1);
            } else if (min1 == 1) {
                // beta最小
                // 交换alpha和beta
                Vector temp = alpha;
                alpha = beta;
                beta = temp;
                
                // 交换z12和z22
                mpz_swap(z12, z22);
                
                mpz_mul_si(tmp1, z11, -bit);
                mpz_add(z22, tmp1, z22);
                
                mpz_mul_si(tmp1, z21, -bit);
                mpz_add(tmp1, tmp1, z12);
                mpz_div_2exp(z12, tmp1, 1);
                
                // 交换z11和z21
                mpz_swap(z11, z21);
                mpz_mul_2exp(z21, z21, 1);
            } else if (min1 == 2) {
                // 2*alpha - beta最小
                // 交换alpha和beta
                Vector temp = alpha;
                alpha = beta;
                beta = temp;
                
                // 交换z12和z22
                mpz_swap(z12, z22);
                
                mpz_mul_si(tmp1, z21, -bit);
                mpz_add(tmp1, tmp1, z12);
                mpz_div_2exp(z12, tmp1, 1);
                
                mpz_mul_si(tmp1, z11, -bit);
                mpz_add(tmp2, tmp1, z22);
                mpz_sub(z22, tmp2, z12);
                
                // 交换z11和z21
                mpz_swap(z11, z21);
                mpz_mul_2exp(tmp1, z21, 1);
                mpz_sub(z21, tmp1, z11);
            } else {
                // 2*alpha + beta最小
                // 交换alpha和beta
                Vector temp = alpha;
                alpha = beta;
                beta = temp;
                
                // 交换z12和z22
                mpz_swap(z12, z22);
                
                mpz_mul_si(tmp1, z21, -bit);
                mpz_add(tmp1, tmp1, z12);
                mpz_div_2exp(z12, tmp1, 1);
                
                mpz_mul_si(tmp1, z11, -bit);
                mpz_add(tmp2, tmp1, z22);
                mpz_add(z22, tmp2, z12);
                
                // 交换z11和z21
                mpz_swap(z11, z21);
                mpz_mul_2exp(tmp1, z21, 1);
                mpz_add(z21, tmp1, z11);
            }
            
            vector_clear(&two_alpha);
            vector_clear(&two_alpha_minus_beta);
            vector_clear(&two_alpha_plus_beta);
        } else {
            // 情况3: -ak*z11 + z12 ≡ 1 mod 2 且 -ak*z21 + z22 ≡ 1 mod 2
            // 计算alpha - beta
            Vector alpha_minus_beta;
            vector_init(&alpha_minus_beta);
            mpz_sub(alpha_minus_beta.x, alpha.x, beta.x);
            mpz_sub(alpha_minus_beta.y, alpha.y, beta.y);
            
            // 计算alpha + beta
            Vector alpha_plus_beta;
            vector_init(&alpha_plus_beta);
            mpz_add(alpha_plus_beta.x, alpha.x, beta.x);
            mpz_add(alpha_plus_beta.y, alpha.y, beta.y);
            
            // 计算2*alpha
            Vector two_alpha;
            vector_init(&two_alpha);
            mpz_mul_2exp(two_alpha.x, alpha.x, 1);
            mpz_mul_2exp(two_alpha.y, alpha.y, 1);
            
            // 计算3*alpha - beta
            Vector three_alpha_minus_beta;
            vector_init(&three_alpha_minus_beta);
            mpz_mul_ui(tmp1, alpha.x, 3);
            mpz_sub(three_alpha_minus_beta.x, tmp1, beta.x);
            mpz_mul_ui(tmp1, alpha.y, 3);
            mpz_sub(three_alpha_minus_beta.y, tmp1, beta.y);
            
            // 计算3*alpha + beta
            Vector three_alpha_plus_beta;
            vector_init(&three_alpha_plus_beta);
            mpz_mul_ui(tmp1, alpha.x, 3);
            mpz_add(three_alpha_plus_beta.x, tmp1, beta.x);
            mpz_mul_ui(tmp1, alpha.y, 3);
            mpz_add(three_alpha_plus_beta.y, tmp1, beta.y);
            
            // 找出最小的两个向量
            Vector candidates[5] = {two_alpha, alpha_minus_beta, alpha_plus_beta, 
                                  three_alpha_minus_beta, three_alpha_plus_beta};
            int min1 = 0, min2 = 1;
            
            for (int j = 2; j < 5; j++) {
                if (vector_cmp(&candidates[j], &candidates[min1]) < 0) {
                    min2 = min1;
                    min1 = j;
                } else if (vector_cmp(&candidates[j], &candidates[min2]) < 0) {
                    min2 = j;
                }
            }
            
            // 更新alpha和beta
            Vector new_alpha = candidates[min1];
            Vector new_beta = candidates[min2];
            
            if (min1 == 0) {
                // 2*alpha最小
                mpz_set(alpha.x, two_alpha.x);
                mpz_set(alpha.y, two_alpha.y);
                
                mpz_mul_si(tmp1, z11, -bit);
                mpz_add(z12, tmp1, z12);
                mpz_mul_2exp(z11, z11, 1);
                
                mpz_mul_si(tmp1, z21, -bit);
                mpz_add(tmp1, tmp1, z22);
                mpz_div_2exp(z22, tmp1, 1);
            } else if (min1 == 1) {
                // alpha - beta最小
                mpz_set(alpha.x, alpha_minus_beta.x);
                mpz_set(alpha.y, alpha_minus_beta.y);
                
                mpz_sub(z11, z11, z21);
                mpz_add(tmp1, z11, z21);
                mpz_add(z21, tmp1, z21);
                
                mpz_sub(z12, z12, z22);
                mpz_add(tmp1, z12, z22);
                mpz_add(z22, tmp1, z22);
                
                mpz_mul_si(tmp1, z11, -bit);
                mpz_add(tmp1, tmp1, z12);
                mpz_div_2exp(z12, tmp1, 1);
                
                mpz_mul_si(tmp1, z21, -bit);
                mpz_add(tmp1, tmp1, z22);
                mpz_div_2exp(z22, tmp1, 1);
            } else if (min1 == 2) {
                // alpha + beta最小
                mpz_set(alpha.x, alpha_plus_beta.x);
                mpz_set(alpha.y, alpha_plus_beta.y);
                
                mpz_add(z11, z11, z21);
                mpz_sub(tmp1, z11, z21);
                mpz_sub(z21, tmp1, z21);
                
                mpz_add(z12, z12, z22);
                mpz_sub(tmp1, z12, z22);
                mpz_sub(z22, tmp1, z22);
                
                mpz_mul_si(tmp1, z11, -bit);
                mpz_add(tmp1, tmp1, z12);
                mpz_div_2exp(z12, tmp1, 1);
                
                mpz_mul_si(tmp1, z21, -bit);
                mpz_add(tmp1, tmp1, z22);
                mpz_div_2exp(z22, tmp1, 1);
            } else if (min1 == 3) {
                // 3*alpha - beta最小
                mpz_set(alpha.x, three_alpha_minus_beta.x);
                mpz_set(alpha.y, three_alpha_minus_beta.y);
                
                mpz_swap(z11, z21);
                mpz_sub(tmp1, z21, z11);
                mpz_set(z11, tmp1);
                mpz_mul_2exp(tmp1, z21, 1);
                mpz_add(z21, tmp1, z11);
                
                mpz_swap(z12, z22);
                mpz_sub(tmp1, z22, z12);
                mpz_set(z12, tmp1);
                mpz_mul_2exp(tmp1, z22, 1);
                mpz_add(z22, tmp1, z12);
                
                mpz_mul_si(tmp1, z11, -bit);
                mpz_add(tmp1, tmp1, z12);
                mpz_div_2exp(z12, tmp1, 1);
                
                mpz_mul_si(tmp1, z21, -bit);
                mpz_add(tmp1, tmp1, z22);
                mpz_div_2exp(z22, tmp1, 1);
            } else {
                // 3*alpha + beta最小
                mpz_set(alpha.x, three_alpha_plus_beta.x);
                mpz_set(alpha.y, three_alpha_plus_beta.y);
                
                mpz_swap(z11, z21);
                mpz_add(tmp1, z21, z11);
                mpz_set(z11, tmp1);
                mpz_mul_2exp(tmp1, z21, 1);
                mpz_add(z21, tmp1, z11);
                
                mpz_swap(z12, z22);
                mpz_add(tmp1, z22, z12);
                mpz_set(z12, tmp1);
                mpz_mul_2exp(tmp1, z22, 1);
                mpz_add(z22, tmp1, z12);
                
                mpz_mul_si(tmp1, z11, -bit);
                mpz_add(tmp1, tmp1, z12);
                mpz_div_2exp(z12, tmp1, 1);
                
                mpz_mul_si(tmp1, z21, -bit);
                mpz_add(tmp1, tmp1, z22);
                mpz_div_2exp(z22, tmp1, 1);
            }
            
            vector_clear(&alpha_minus_beta);
            vector_clear(&alpha_plus_beta);
            vector_clear(&two_alpha);
            vector_clear(&three_alpha_minus_beta);
            vector_clear(&three_alpha_plus_beta);
        }
    }
    
    // 检查哪个向量的第二个分量是奇数
    mpz_mod_ui(tmp1, alpha.y, 2);
    if (mpz_cmp_ui(tmp1, 1) == 0) {
        mpz_set(p, alpha.x);
        mpz_set(q, alpha.y);
    } else {
        mpz_set(p, beta.x);
        mpz_set(q, beta.y);
    }
    
    // 清理内存
    vector_clear(&alpha);
    vector_clear(&beta);
    mpz_clear(z11); mpz_clear(z12); mpz_clear(z21); mpz_clear(z22);
    mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3);
}

// 全局欧几里得算法(Glo-Euc)
void glo_euc(const char *sequence, int n, mpz_t p, mpz_t q) {
    // 初始化向量
    Vector a, b;
    vector_init(&a);
    vector_init(&b);
    
    // 计算Sn
    mpz_t Sn;
    mpz_init(Sn);
    mpz_set_ui(Sn, 0);
    
    for (int i = 0; i < n; i++) {
        if (sequence[i] == '1') {
            mpz_ui_pow_ui(a.x, 2, i);
            mpz_add(Sn, Sn, a.x);
        }
    }
    
    // 设置初始向量
    mpz_set(a.x, Sn);
    mpz_set_ui(a.y, 1);
    
    mpz_ui_pow_ui(b.x, 2, n);
    mpz_set_ui(b.y, 0);
    
    // 执行Lagrange约简算法
    while (1) {
        // 步骤1: 如果||a|| > ||b||，交换a和b
        if (vector_cmp(&a, &b) > 0) {
            Vector temp = a;
            a = b;
            b = temp;
        }
        
        // 步骤2: 如果||a-b|| > ||a+b||，设置b = -b
        Vector a_minus_b, a_plus_b;
        vector_init(&a_minus_b);
        vector_init(&a_plus_b);
        
        mpz_sub(a_minus_b.x, a.x, b.x);
        mpz_sub(a_minus_b.y, a.y, b.y);
        
        mpz_add(a_plus_b.x, a.x, b.x);
        mpz_add(a_plus_b.y, a.y, b.y);
        
        if (vector_cmp(&a_minus_b, &a_plus_b) > 0) {
            mpz_neg(b.x, b.x);
            mpz_neg(b.y, b.y);
        }
        
        vector_clear(&a_minus_b);
        vector_clear(&a_plus_b);
        
        // 步骤3: 如果||a|| <= ||b||，进入循环
        if (vector_cmp(&a, &b) <= 0) {
            break;
        }
        
        // 步骤4: 如果||a|| == ||b||，返回[a, -b]
        if (vector_cmp(&a, &b) == 0) {
            mpz_neg(b.x, b.x);
            mpz_neg(b.y, b.y);
            break;
        }
        
        // 步骤5: 设置a = b - a
        mpz_sub(a.x, b.x, a.x);
        mpz_sub(a.y, b.y, a.y);
    }
    
    // 检查哪个向量的第二个分量是奇数
    mpz_t tmp;
    mpz_init(tmp);
    mpz_mod_ui(tmp, a.y, 2);
    if (mpz_cmp_ui(tmp, 1) == 0) {
        mpz_set(p, a.x);
        mpz_set(q, a.y);
    } else {
        mpz_set(p, b.x);
        mpz_set(q, b.y);
    }
    
    // 清理内存
    vector_clear(&a);
    vector_clear(&b);
    mpz_clear(Sn);
    mpz_clear(tmp);
}

// 部分欧几里得算法(Par-Euc)
void par_euc(const char *sequence, int n, mpz_t p, mpz_t q) {
    // 初始化变量
    mpz_t r0, r1, u0, u1, v0, v1;
    mpz_t tmp, q_val, r, u2, v2;
    
    mpz_init(r0); mpz_init(r1); mpz_init(u0); mpz_init(u1); mpz_init(v0); mpz_init(v1);
    mpz_init(tmp); mpz_init(q_val); mpz_init(r); mpz_init(u2); mpz_init(v2);
    
    // 初始化r0 = 2^n, u0 = 1, v0 = 0
    mpz_ui_pow_ui(r0, 2, n);
    mpz_set_ui(u0, 1);
    mpz_set_ui(v0, 0);
    
    // 初始化r1 = Sn, u1 = 0, v1 = 1
    mpz_set_ui(r1, 0);
    for (int i = 0; i < n; i++) {
        if (sequence[i] == '1') {
            mpz_ui_pow_ui(tmp, 2, i);
            mpz_add(r1, r1, tmp);
        }
    }
    mpz_set_ui(u1, 0);
    mpz_set_ui(v1, 1);
    
    // 执行部分欧几里得算法
    while (mpz_cmp_ui(r1, 0) > 0) {
        // 计算商和余数
        mpz_fdiv_q(q_val, r0, r1);
        mpz_mod(r, r0, r1);
        
        // 更新u和v
        mpz_mul(tmp, q_val, u1);
        mpz_sub(u2, u0, tmp);
        
        mpz_mul(tmp, q_val, v1);
        mpz_sub(v2, v0, tmp);
        
        // 更新r, u, v
        mpz_set(r0, r1);
        mpz_set(u0, u1);
        mpz_set(v0, v1);
        
        mpz_set(r1, r);
        mpz_set(u1, u2);
        mpz_set(v1, v2);
        
        // 检查终止条件
        mpz_abs(tmp, v1);
        if (mpz_cmp(r1, tmp) <= 0) {
            break;
        }
    }
    
    // 检查v1是否为奇数
    mpz_mod_ui(tmp, v1, 2);
    if (mpz_cmp_ui(tmp, 1) == 0) {
        mpz_set(p, r1);
        mpz_set(q, v1);
    } else {
        // 需要进一步处理
        mpz_t x, y, v1_abs, r1_abs, v0_abs, r0_abs;
        mpz_init(x); mpz_init(y); mpz_init(v1_abs); mpz_init(r1_abs); 
        mpz_init(v0_abs); mpz_init(r0_abs);
        
        mpz_abs(v1_abs, v1);
        mpz_abs(r1_abs, r1);
        mpz_abs(v0_abs, v0);
        mpz_abs(r0_abs, r0);
        
        // 计算x = (|v1| - r1) / (v0 + r0)
        mpz_sub(tmp, v1_abs, r1_abs);
        mpz_add(x, v0_abs, r0_abs);
        mpz_fdiv_q(x, tmp, x);
        
        // 计算y = (r0 - |v0|) / (r1 + |v1|)
        mpz_sub(tmp, r0_abs, v0_abs);
        mpz_add(y, r1_abs, v1_abs);
        mpz_fdiv_q(y, tmp, y);
        
        // 计算候选向量
        Vector cand1, cand2;
        vector_init(&cand1);
        vector_init(&cand2);
        
        // 候选1: floor(x)
        mpz_set(tmp, x);
        mpz_mul(cand1.x, tmp, r0);
        mpz_sub(cand1.x, r1, cand1.x);
        mpz_mul(cand1.y, tmp, v0);
        mpz_sub(cand1.y, v1, cand1.y);
        
        // 候选2: ceil(x)
        mpz_add_ui(tmp, x, 1);
        mpz_mul(cand2.x, tmp, r0);
        mpz_sub(cand2.x, r1, cand2.x);
        mpz_mul(cand2.y, tmp, v0);
        mpz_sub(cand2.y, v1, cand2.y);
        
        // 选择较小的候选
        if (vector_cmp(&cand1, &cand2) <= 0) {
            mpz_set(p, cand1.x);
            mpz_set(q, cand1.y);
        } else {
            mpz_set(p, cand2.x);
            mpz_set(q, cand2.y);
        }
        
        // 清理内存
        vector_clear(&cand1);
        vector_clear(&cand2);
        mpz_clear(x); mpz_clear(y); mpz_clear(v1_abs); mpz_clear(r1_abs);
        mpz_clear(v0_abs); mpz_clear(r0_abs);
    }
    
    // 清理内存
    mpz_clear(r0); mpz_clear(r1); mpz_clear(u0); mpz_clear(u1); mpz_clear(v0); mpz_clear(v1);
    mpz_clear(tmp); mpz_clear(q_val); mpz_clear(r); mpz_clear(u2); mpz_clear(v2);
}

// 生成随机二进制序列
void generate_random_sequence(char *sequence, int n) {
    for (int i = 0; i < n; i++) {
        sequence[i] = (rand() % 2) ? '1' : '0';
    }
    sequence[n] = '\0';
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s <algorithm> <sequence_length>\n", argv[0]);
        printf("Algorithms: improved_ra, glo_euc, par_euc\n");
        return 1;
    }
    
    char *algorithm = argv[1];
    int n = atoi(argv[2]);
    
    if (n <= 0) {
        printf("Sequence length must be positive\n");
        return 1;
    }
    
    // 生成随机序列
    char *sequence = (char *)malloc(n + 1);
    generate_random_sequence(sequence, n);
    
    printf("Sequence: %s\n", sequence);
    
    // 初始化结果变量
    mpz_t p, q;
    mpz_init(p);
    mpz_init(q);
    
    // 记录开始时间
    clock_t start = clock();
    
    // 调用相应的算法
    if (strcmp(algorithm, "improved_ra") == 0) {
        improved_ra(sequence, n, p, q);
    } else if (strcmp(algorithm, "glo_euc") == 0) {
        glo_euc(sequence, n, p, q);
    } else if (strcmp(algorithm, "par_euc") == 0) {
        par_euc(sequence, n, p, q);
    } else {
        printf("Unknown algorithm: %s\n", algorithm);
        free(sequence);
        mpz_clear(p);
        mpz_clear(q);
        return 1;
    }
    
    // 记录结束时间
    clock_t end = clock();
    double time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    // 输出结果
    gmp_printf("Result: p = %Zd, q = %Zd\n", p, q);
    printf("Time used: %.6f seconds\n", time_used);
    
    // 清理内存
    free(sequence);
    mpz_clear(p);
    mpz_clear(q);
    
    return 0;
}
