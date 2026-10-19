// SPDX-FileCopyrightText: 2026 The project author
// SPDX-License-Identifier: Apache-2.0

#include "protocols.h"
#include "protocols_helper.h"
#include "protocol_constants.h"
#include "rng.h"
#include <assert.h>
#include <stdio.h>

void make_SecretKey(secretkey_t sk)
{
    for(int i = 0; i < P_LEN + M_LEN; i++)
    {
        do
        {
            randombytes((unsigned char *)&(sk[i]), sizeof(sk[i]));
            sk[i] &= secretkey_mask[i];
        } while (sk[i] > ITERATIONS * TORSION_ODD_POWERS[i]);
    }
}

void print_SecretKey(const secretkey_t sk)
{
    for(int i = 0; i < P_LEN + M_LEN - 1; i++)
    {
        printf("%d,", sk[i]);
    }
    printf("%d", sk[P_LEN + M_LEN - 1]);
    printf("\n");
}

void print_PublicKey(const orient_cycle_t *pk)
{
    print_cycle(pk);
    printf("\n");
}

void print_SecretAgreement(const fp2_t a)
{
    fp2_print(a);
    printf("\n");
}

void eval_action(orient_cycle_t* res, const orient_cycle_t* inp, const int *s)
{
    assert(res != inp);
    vertex_t R[2];
    
    for(int j = 0; j < CYCLE_LENGTH; j++)
    {
        copy_vertices(&(R[0]), &((*inp)[j]));
        copy_vertices_swap_kernel(&(R[1]), &((*inp)[j + 1 == CYCLE_LENGTH ? 0 : j + 1]));
        int mscalar_i, mscalar_bit_i;
        ec_point_t B24;

        // Isogenies with kernel in E[p+1]
        mscalar_i = 0;
        mscalar_bit_i = 0;
        for(int i = 0; i < P_LEN; i++)
        {
            for(int t = 1; t <= TORSION_ODD_POWERS[i]; t++)
            {
                digit_t b = compare_gt(t, s[i]);
                swap_vertices(&(R[0]), &(R[1]), b);
                // int b = compare_gt3(t, s[i]);
                ec_point_t K;
                xMULv2(&K, &((R[0]).Ps), &(mscalar_s[mscalar_i]), mscalar_s_bits[mscalar_bit_i], &((R[0]).A24));
                mscalar_i += mscalar_s_limbs[mscalar_bit_i];
                mscalar_bit_i++;
                kps(i, K, R[0].A24);
                xisog(&B24, i, R[0].A24);
                xeval(&(R[0].Ps), i, R[0].Ps, R[0].A24);
                xeval(&(R[0].Pt), i, R[0].Pt, R[0].A24);
                xeval(&(R[0].Qs), i, R[0].Qs, R[0].A24);
                xeval(&(R[0].Qt), i, R[0].Qt, R[0].A24);
                copy_point(&(R[0].A24), &B24);
                kps_clear(i);
                xMULv2(&(R[1].Ps), &(R[1].Ps), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &(R[1].A24));
                swap_vertices(&(R[0]), &(R[1]), b);
            }
        }

        // Isogenies with kernel in E[p-1]
        mscalar_i = 0;
        mscalar_bit_i = 0;
        for(int i = P_LEN; i < P_LEN + M_LEN; i++)
        {
            for(int t = 1; t <= TORSION_ODD_POWERS[i]; t++)
            {
                digit_t b = compare_gt(t, s[i]);
                swap_vertices(&(R[0]), &(R[1]), b);
                ec_point_t K;
                xMULv2(&K, &(R[0].Pt), &(mscalar_t[mscalar_i]), mscalar_t_bits[mscalar_bit_i], &(R[0].A24));
                mscalar_i += mscalar_t_limbs[mscalar_bit_i];
                mscalar_bit_i++;
                kps(i, K, R[0].A24);
                xisog(&B24, i, R[0].A24);
                xeval(&(R[0].Pt), i, R[0].Pt, R[0].A24);
                xeval(&(R[0].Qs), i, R[0].Qs, R[0].A24);
                xeval(&(R[0].Qt), i, R[0].Qt, R[0].A24);
                copy_point(&(R[0].A24), &B24);
                kps_clear(i);
                xMULv2(&(R[1].Pt), &(R[1].Pt), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &(R[1].A24));
                swap_vertices(&(R[0]), &(R[1]), b);
            }
        }
        

        #ifndef NDEBUG
        ec_curve_t E1, E2;
        A24_to_AC(&E1, &(R[0].A24));
        A24_to_AC(&E2, &(R[1].A24));
        fp2_t j1, j2;
        ec_j_inv(&j1, &E1);
        ec_j_inv(&j2, &E2);
        assert(fp2_is_equal(&j1, &j2));
        if(j < CYCLE_LENGTH - 1)
        {
            fp2_t t0, t1;
            fp2_copy(&t0, &(R[0].A24.z));
            fp2_inv(&t0);
            fp2_mul(&t0, &t0, &(R[0].A24.x));
            fp2_copy(&t1, &(R[1].A24.z));
            fp2_inv(&t1);
            fp2_mul(&t1, &t1, &(R[1].A24.x));
            assert(fp2_is_equal(&(t0), &(t1)));
        }
        #endif

        copy_point(&(R[0].Ps), &(R[1].Qs));
        copy_point(&(R[0].Pt), &(R[1].Qt));
        copy_vertices(&((*res)[j]), &(R[0]));
    }

    ec_isom_t correction;
    ec_curve_t head, tail;
    A24_to_AC(&head, &(R[1].A24));
    A24_to_AC(&tail, &(R[0].A24));
    ec_isomorphism(&correction, &head, &tail);
    ec_iso_eval(&((*res)[CYCLE_LENGTH - 1].Ps), &correction);
    ec_iso_eval(&((*res)[CYCLE_LENGTH - 1].Pt), &correction);
}

void eval_action_stra_del(orient_cycle_t* res, const orient_cycle_t* inp, const int *s)
{
    assert(res != inp);
    vertex_t R[2];
    
    for(int j = 0; j < CYCLE_LENGTH; j++)
    {
        copy_vertices(&(R[0]), &((*inp)[j]));
        copy_vertices_swap_kernel(&(R[1]), &((*inp)[j + 1 == CYCLE_LENGTH ? 0 : j + 1]));
        int strascalar_i/*, strascalar_bit_i*/;
        ec_point_t B24;

        // Isogenies with kernel in E[p+1]
        strascalar_i = 0;
        ec_point_t stack_kernel[2][stack_volume];
        int stack_scalar[stack_volume];
        size_t top = 0;
        copy_point(&(stack_kernel[0][top]), &(R[0].Ps));
        copy_point(&(stack_kernel[1][top]), &(R[1].Ps));
        stack_scalar[top] = 0;
        top++;
        int jj = 0, k = 0;
        for(int i = 0; i < P_LEN; i++)
        {
            for(int t = 1; t <= TORSION_ODD_POWERS[i]; t++)
            {
                digit_t b = compare_gt(t, s[i]);
                jj++;
                while(stack_scalar[top-1] < degree_s_num - jj)
                {
                    xMULv2(&(stack_kernel[0][top]), &(stack_kernel[0][top-1]), &(strascalar_s[strascalar_i]), strascalar_s_bits[k], &(R[0].A24));
                    xMULv2(&(stack_kernel[1][top]), &(stack_kernel[1][top-1]), &(strascalar_s[strascalar_i]), strascalar_s_bits[k], &(R[1].A24));
                    strascalar_i += strascalar_s_limbs[k];
                    stack_scalar[top] = stack_scalar[top-1] + strategy_s[k];
                    top++;
                    k++;
                    assert(top <= stack_volume);
                }
                for (int si = 0; si < top; si++)
                {
                    swap_points(&(stack_kernel[0][si]), &(stack_kernel[1][si]), b);
                }
                swap_vertices(&(R[0]), &(R[1]), b);
                top--;
                assert(top >= 0);
                kps(i, stack_kernel[0][top], R[0].A24);
                xisog(&B24, i, R[0].A24);
                for(int traverse_stack = 0; traverse_stack < top; traverse_stack++)
                {
                    xeval(&(stack_kernel[0][traverse_stack]), i, stack_kernel[0][traverse_stack], R[0].A24);
                }
                xeval(&(R[0].Pt), i, R[0].Pt, R[0].A24);
                xeval(&(R[0].Qs), i, R[0].Qs, R[0].A24);
                xeval(&(R[0].Qt), i, R[0].Qt, R[0].A24);
                copy_point(&(R[0].A24), &B24);
                kps_clear(i);
                for(int traverse_stack = 0; traverse_stack < top; traverse_stack++)
                {
                    xMULv2(&(stack_kernel[1][traverse_stack]), &(stack_kernel[1][traverse_stack]), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &(R[1].A24));
                }
                for (int si = 0; si < top; si++)
                {
                    swap_points(&(stack_kernel[0][si]), &(stack_kernel[1][si]), b);
                }
                swap_vertices(&(R[0]), &(R[1]), b);
            }
        }
        assert(top == 0);

        // Isogenies with kernel in E[p-1]
        strascalar_i = 0;
        top = 0;
        copy_point(&(stack_kernel[0][top]), &(R[0].Pt));
        copy_point(&(stack_kernel[1][top]), &(R[1].Pt));
        stack_scalar[top] = 0;
        top++;
        jj = 0; k = 0;
        for(int i = P_LEN; i < P_LEN + M_LEN; i++)
        {
            for(int t = 1; t <= TORSION_ODD_POWERS[i]; t++)
            {
                digit_t b = compare_gt(t, s[i]);
                jj++;
                while(stack_scalar[top-1] < degree_t_num - jj)
                {
                    xMULv2(&(stack_kernel[0][top]), &(stack_kernel[0][top-1]), &(strascalar_t[strascalar_i]), strascalar_t_bits[k], &(R[0].A24));
                    xMULv2(&(stack_kernel[1][top]), &(stack_kernel[1][top-1]), &(strascalar_t[strascalar_i]), strascalar_t_bits[k], &(R[1].A24));
                    strascalar_i += strascalar_t_limbs[k];
                    stack_scalar[top] = stack_scalar[top-1] + strategy_t[k];
                    top++;
                    k++;
                    assert(top <= stack_volume);
                }
                for (int si = 0; si < top; si++)
                {
                    swap_points(&(stack_kernel[0][si]), &(stack_kernel[1][si]), b);
                }
                swap_vertices(&(R[0]), &(R[1]), b);
                top--;
                assert(top >= 0);
                kps(i, stack_kernel[0][top], R[0].A24);
                xisog(&B24, i, R[0].A24);
                for(int traverse_stack = 0; traverse_stack < top; traverse_stack++)
                {
                    xeval(&(stack_kernel[0][traverse_stack]), i, stack_kernel[0][traverse_stack], R[0].A24);
                }
                xeval(&(R[0].Qs), i, R[0].Qs, R[0].A24);
                xeval(&(R[0].Qt), i, R[0].Qt, R[0].A24);
                copy_point(&(R[0].A24), &B24);
                kps_clear(i);
                for(int traverse_stack = 0; traverse_stack < top; traverse_stack++)
                {
                    xMULv2(&(stack_kernel[1][traverse_stack]), &(stack_kernel[1][traverse_stack]), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &(R[1].A24));
                }
                for (int si = 0; si < top; si++)
                {
                    swap_points(&(stack_kernel[0][si]), &(stack_kernel[1][si]), b);
                }
                swap_vertices(&(R[0]), &(R[1]), b);
            }
        }
        assert(top == 0);        

        #ifndef NDEBUG
        ec_curve_t E1, E2;
        A24_to_AC(&E1, &(R[0].A24));
        A24_to_AC(&E2, &(R[1].A24));
        fp2_t j1, j2;
        ec_j_inv(&j1, &E1);
        ec_j_inv(&j2, &E2);
        assert(fp2_is_equal(&j1, &j2));
        if(j < CYCLE_LENGTH - 1)
        {
            fp2_t t0, t1;
            fp2_copy(&t0, &(R[0].A24.z));
            fp2_inv(&t0);
            fp2_mul(&t0, &t0, &(R[0].A24.x));
            fp2_copy(&t1, &(R[1].A24.z));
            fp2_inv(&t1);
            fp2_mul(&t1, &t1, &(R[1].A24.x));
            assert(fp2_is_equal(&(t0), &(t1)));
        }
        #endif

        copy_point(&(R[0].Ps), &(R[1].Qs));
        copy_point(&(R[0].Pt), &(R[1].Qt));
        copy_vertices(&((*res)[j]), &R[0]);
    }

    ec_isom_t correction;
    ec_curve_t head, tail;
    A24_to_AC(&head, &(R[1].A24));
    A24_to_AC(&tail, &(R[0].A24));
    ec_isomorphism(&correction, &head, &tail);
    ec_iso_eval(&((*res)[CYCLE_LENGTH - 1].Ps), &correction);
    ec_iso_eval(&((*res)[CYCLE_LENGTH - 1].Pt), &correction);
}

void eval_action_stra(orient_cycle_t* res, const orient_cycle_t* inp, const int *s)
{
    assert(res != inp);
    vertex_t R[2];
    
    for(int j = 0; j < CYCLE_LENGTH; j++)
    {
        copy_vertices(&(R[0]), &((*inp)[j]));
        copy_vertices_swap_kernel(&(R[1]), &((*inp)[j + 1 == CYCLE_LENGTH ? 0 : j + 1]));
        int strascalar_i/*, strascalar_bit_i*/;
        ec_point_t B24;

        // Isogenies with kernel in E[p+1]
        strascalar_i = 0;
        ec_point_t stack_kernel[2][stack_volume];
        int stack_scalar[stack_volume];
        size_t top = 0;
        copy_point(&(stack_kernel[0][top]), &(R[0].Ps));
        copy_point(&(stack_kernel[1][top]), &(R[1].Ps));
        stack_scalar[top] = 0;
        top++;
        int jj = 0, k = 0;
        for(int i = 0; i < P_LEN; i++)
        {
            for(int t = 1; t <= TORSION_ODD_POWERS[i]; t++)
            {
                digit_t b = compare_gt(t, s[i]);
                jj++;
                for (int si = 0; si < top; si++)
                {
                    swap_points(&(stack_kernel[0][si]), &(stack_kernel[1][si]), b);
                }
                swap_vertices(&(R[0]), &(R[1]), b);
                while(stack_scalar[top-1] < degree_s_num - jj)
                {
                    stack_scalar[top] = stack_scalar[top-1] + strategy_s[k];
                    xMULv2(&(stack_kernel[0][top]), &(stack_kernel[0][top-1]), &(strascalar_s[strascalar_i]), strascalar_s_bits[k], &(R[0].A24));
                    if (stack_scalar[top] < degree_s_num - jj) // This conditional check does NOT leak any secret information
                    {
                        xMULv2(&(stack_kernel[1][top]), &(stack_kernel[1][top-1]), &(strascalar_s[strascalar_i]), strascalar_s_bits[k], &(R[1].A24));
                    }
                    strascalar_i += strascalar_s_limbs[k];
                    top++;
                    k++;
                    assert(top <= stack_volume);
                }
                top--;
                assert(top >= 0);
                kps(i, stack_kernel[0][top], R[0].A24);
                xisog(&B24, i, R[0].A24);
                for(int traverse_stack = 0; traverse_stack < top; traverse_stack++)
                {
                    xeval(&(stack_kernel[0][traverse_stack]), i, stack_kernel[0][traverse_stack], R[0].A24);
                }
                xeval(&(R[0].Pt), i, R[0].Pt, R[0].A24);
                xeval(&(R[0].Qs), i, R[0].Qs, R[0].A24);
                xeval(&(R[0].Qt), i, R[0].Qt, R[0].A24);
                copy_point(&(R[0].A24), &B24);
                kps_clear(i);
                for(int traverse_stack = 0; traverse_stack < top; traverse_stack++)
                {
                    xMULv2(&(stack_kernel[1][traverse_stack]), &(stack_kernel[1][traverse_stack]), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &(R[1].A24));
                }
                for (int si = 0; si < top; si++)
                {
                    swap_points(&(stack_kernel[0][si]), &(stack_kernel[1][si]), b);
                }
                swap_vertices(&(R[0]), &(R[1]), b);
            }
        }
        assert(top == 0);

        // Isogenies with kernel in E[p-1]
        strascalar_i = 0;
        top = 0;
        copy_point(&(stack_kernel[0][top]), &(R[0].Pt));
        copy_point(&(stack_kernel[1][top]), &(R[1].Pt));
        stack_scalar[top] = 0;
        top++;
        jj = 0; k = 0;
        for(int i = P_LEN; i < P_LEN + M_LEN; i++)
        {
            for(int t = 1; t <= TORSION_ODD_POWERS[i]; t++)
            {
                digit_t b = compare_gt(t, s[i]);
                jj++;
                for (int si = 0; si < top; si++)
                {
                    swap_points(&(stack_kernel[0][si]), &(stack_kernel[1][si]), b);
                }
                swap_vertices(&(R[0]), &(R[1]), b);
                while(stack_scalar[top-1] < degree_t_num - jj)
                {
                    xMULv2(&(stack_kernel[0][top]), &(stack_kernel[0][top-1]), &(strascalar_t[strascalar_i]), strascalar_t_bits[k], &(R[0].A24));
                    if (stack_scalar[top] < degree_s_num - jj) // This conditional check does NOT leak any secret information
                    {
                        xMULv2(&(stack_kernel[1][top]), &(stack_kernel[1][top-1]), &(strascalar_t[strascalar_i]), strascalar_t_bits[k], &(R[1].A24));
                    }
                    strascalar_i += strascalar_t_limbs[k];
                    stack_scalar[top] = stack_scalar[top-1] + strategy_t[k];
                    top++;
                    k++;
                    assert(top <= stack_volume);
                }
                top--;
                assert(top >= 0);
                kps(i, stack_kernel[0][top], R[0].A24);
                xisog(&B24, i, R[0].A24);
                for(int traverse_stack = 0; traverse_stack < top; traverse_stack++)
                {
                    xeval(&(stack_kernel[0][traverse_stack]), i, stack_kernel[0][traverse_stack], R[0].A24);
                }
                xeval(&(R[0].Qs), i, R[0].Qs, R[0].A24);
                xeval(&(R[0].Qt), i, R[0].Qt, R[0].A24);
                copy_point(&(R[0].A24), &B24);
                kps_clear(i);
                for(int traverse_stack = 0; traverse_stack < top; traverse_stack++)
                {
                    xMULv2(&(stack_kernel[1][traverse_stack]), &(stack_kernel[1][traverse_stack]), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &(R[1].A24));
                }
                for (int si = 0; si < top; si++)
                {
                    swap_points(&(stack_kernel[0][si]), &(stack_kernel[1][si]), b);
                }
                swap_vertices(&(R[0]), &(R[1]), b);
            }
        }
        assert(top == 0);        

        #ifndef NDEBUG
        ec_curve_t E1, E2;
        A24_to_AC(&E1, &(R[0].A24));
        A24_to_AC(&E2, &(R[1].A24));
        fp2_t j1, j2;
        ec_j_inv(&j1, &E1);
        ec_j_inv(&j2, &E2);
        assert(fp2_is_equal(&j1, &j2));
        if(j < CYCLE_LENGTH - 1)
        {
            fp2_t t0, t1;
            fp2_copy(&t0, &(R[0].A24.z));
            fp2_inv(&t0);
            fp2_mul(&t0, &t0, &(R[0].A24.x));
            fp2_copy(&t1, &(R[1].A24.z));
            fp2_inv(&t1);
            fp2_mul(&t1, &t1, &(R[1].A24.x));
            assert(fp2_is_equal(&(t0), &(t1)));
        }
        #endif

        copy_point(&(R[0].Ps), &(R[1].Qs));
        copy_point(&(R[0].Pt), &(R[1].Qt));
        copy_vertices(&((*res)[j]), &R[0]);
    }

    ec_isom_t correction;
    ec_curve_t head, tail;
    A24_to_AC(&head, &(R[1].A24));
    A24_to_AC(&tail, &(R[0].A24));
    ec_isomorphism(&correction, &head, &tail);
    ec_iso_eval(&((*res)[CYCLE_LENGTH - 1].Ps), &correction);
    ec_iso_eval(&((*res)[CYCLE_LENGTH - 1].Pt), &correction);
}

void eval_action_omit(orient_cycle_t* res, const orient_cycle_t* inp, const int *s, int omit)
{
    assert(res != inp);
    assert(omit > 0);
    assert(omit < CYCLE_LENGTH);

    vertex_t R[2];
    
    for(int j = 0; j < CYCLE_LENGTH - omit; j++)
    {
        copy_vertices(&(R[0]), &((*inp)[j]));
        copy_vertices_swap_kernel(&(R[1]), &((*inp)[j + 1]));
        int mscalar_i, mscalar_bit_i;
        ec_point_t B24;

        // Isogenies with kernel in E[p+1]
        mscalar_i = 0;
        mscalar_bit_i = 0;
        for(int i = 0; i < P_LEN; i++)
        {
            for(int t = 1; t <= TORSION_ODD_POWERS[i]; t++)
            {
                digit_t b = compare_gt(t, s[i]);
                swap_vertices(&(R[0]), &(R[1]), b);
                ec_point_t K;
                xMULv2(&K, &(R[0].Ps), &(mscalar_s[mscalar_i]), mscalar_s_bits[mscalar_bit_i], &(R[0].A24));
                mscalar_i += mscalar_s_limbs[mscalar_bit_i];
                mscalar_bit_i++;
                kps(i, K, R[0].A24);
                xisog(&B24, i, R[0].A24);
                xeval(&(R[0].Ps), i, R[0].Ps, R[0].A24);
                xeval(&(R[0].Pt), i, R[0].Pt, R[0].A24);
                xeval(&(R[0].Qs), i, R[0].Qs, R[0].A24);
                xeval(&(R[0].Qt), i, R[0].Qt, R[0].A24);
                copy_point(&(R[0].A24), &B24);
                kps_clear(i);
                xMULv2(&(R[1].Ps), &(R[1].Ps), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &(R[1].A24));
                swap_vertices(&(R[0]), &(R[1]), b);
            }
        }

        // Isogenies with kernel in E[p-1]
        mscalar_i = 0;
        mscalar_bit_i = 0;
        for(int i = P_LEN; i < P_LEN + M_LEN; i++)
        {
            for(int t = 1; t <= TORSION_ODD_POWERS[i]; t++)
            {
                digit_t b = compare_gt(t, s[i]);
                swap_vertices(&(R[0]), &(R[1]), b);
                ec_point_t K;
                xMULv2(&K, &(R[0].Pt), &(mscalar_t[mscalar_i]), mscalar_t_bits[mscalar_bit_i], &(R[0].A24));
                mscalar_i += mscalar_t_limbs[mscalar_bit_i];
                mscalar_bit_i++;
                kps(i, K, R[0].A24);
                xisog(&B24, i, R[0].A24);
                xeval(&(R[0].Pt), i, R[0].Pt, R[0].A24);
                xeval(&(R[0].Qs), i, R[0].Qs, R[0].A24);
                xeval(&(R[0].Qt), i, R[0].Qt, R[0].A24);
                copy_point(&(R[0].A24), &B24);
                kps_clear(i);
                xMULv2(&(R[1].Pt), &(R[1].Pt), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &(R[1].A24));
                swap_vertices(&(R[0]), &(R[1]), b);
            }
        }
        

        #ifndef NDEBUG
        ec_curve_t E1, E2;
        A24_to_AC(&E1, &(R[0].A24));
        A24_to_AC(&E2, &(R[1].A24));
        fp2_t j1, j2;
        ec_j_inv(&j1, &E1);
        ec_j_inv(&j2, &E2);
        assert(fp2_is_equal(&j1, &j2));
        if(j < CYCLE_LENGTH - 1)
        {
            fp2_t t0, t1;
            fp2_copy(&t0, &(R[0].A24.z));
            fp2_inv(&t0);
            fp2_mul(&t0, &t0, &(R[0].A24.x));
            fp2_copy(&t1, &(R[1].A24.z));
            fp2_inv(&t1);
            fp2_mul(&t1, &t1, &(R[1].A24.x));
            assert(fp2_is_equal(&(t0), &(t1)));
        }
        #endif

        copy_point(&(R[0].Ps), &(R[1].Qs));
        copy_point(&(R[0].Pt), &(R[1].Qt));
        copy_vertices(&((*res)[j]), &(R[0]));
    }
}

void eval_action_stra_omit(orient_cycle_t* res, const orient_cycle_t* inp, const int *s, int omit)
{
    assert(res != inp);
    assert(omit > 0);
    assert(omit < CYCLE_LENGTH);
    vertex_t R[2];
    
    for(int j = 0; j < CYCLE_LENGTH - omit; j++)
    {
        copy_vertices(&(R[0]), &((*inp)[j]));
        copy_vertices_swap_kernel(&(R[1]), &((*inp)[j + 1 == CYCLE_LENGTH ? 0 : j + 1]));
        int strascalar_i/*, strascalar_bit_i*/;
        ec_point_t B24;

        // Isogenies with kernel in E[p+1]
        strascalar_i = 0;
        ec_point_t stack_kernel[2][stack_volume];
        int stack_scalar[stack_volume];
        size_t top = 0;
        copy_point(&(stack_kernel[0][top]), &(R[0].Ps));
        copy_point(&(stack_kernel[1][top]), &(R[1].Ps));
        stack_scalar[top] = 0;
        top++;
        int jj = 0, k = 0;
        for(int i = 0; i < P_LEN; i++)
        {
            for(int t = 1; t <= TORSION_ODD_POWERS[i]; t++)
            {
                digit_t b = compare_gt(t, s[i]);
                jj++;
                for (int si = 0; si < top; si++)
                {
                    swap_points(&(stack_kernel[0][si]), &(stack_kernel[1][si]), b);
                }
                swap_vertices(&(R[0]), &(R[1]), b);
                while(stack_scalar[top-1] < degree_s_num - jj)
                {
                    stack_scalar[top] = stack_scalar[top-1] + strategy_s[k];
                    xMULv2(&(stack_kernel[0][top]), &(stack_kernel[0][top-1]), &(strascalar_s[strascalar_i]), strascalar_s_bits[k], &(R[0].A24));
                    if (stack_scalar[top] < degree_s_num - jj) // This conditional check does NOT leak any secret information
                    {
                        xMULv2(&(stack_kernel[1][top]), &(stack_kernel[1][top-1]), &(strascalar_s[strascalar_i]), strascalar_s_bits[k], &(R[1].A24));
                    }
                    strascalar_i += strascalar_s_limbs[k];
                    top++;
                    k++;
                    assert(top <= stack_volume);
                }
                top--;
                assert(top >= 0);
                kps(i, stack_kernel[0][top], R[0].A24);
                xisog(&B24, i, R[0].A24);
                for(int traverse_stack = 0; traverse_stack < top; traverse_stack++)
                {
                    xeval(&(stack_kernel[0][traverse_stack]), i, stack_kernel[0][traverse_stack], R[0].A24);
                }
                xeval(&(R[0].Pt), i, R[0].Pt, R[0].A24);
                xeval(&(R[0].Qs), i, R[0].Qs, R[0].A24);
                xeval(&(R[0].Qt), i, R[0].Qt, R[0].A24);
                copy_point(&(R[0].A24), &B24);
                kps_clear(i);
                for(int traverse_stack = 0; traverse_stack < top; traverse_stack++)
                {
                    xMULv2(&(stack_kernel[1][traverse_stack]), &(stack_kernel[1][traverse_stack]), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &(R[1].A24));
                }
                for (int si = 0; si < top; si++)
                {
                    swap_points(&(stack_kernel[0][si]), &(stack_kernel[1][si]), b);
                }
                swap_vertices(&(R[0]), &(R[1]), b);
            }
        }
        assert(top == 0);

        // Isogenies with kernel in E[p-1]
        strascalar_i = 0;
        top = 0;
        copy_point(&(stack_kernel[0][top]), &(R[0].Pt));
        copy_point(&(stack_kernel[1][top]), &(R[1].Pt));
        stack_scalar[top] = 0;
        top++;
        jj = 0; k = 0;
        for(int i = P_LEN; i < P_LEN + M_LEN; i++)
        {
            for(int t = 1; t <= TORSION_ODD_POWERS[i]; t++)
            {
                digit_t b = compare_gt(t, s[i]);
                jj++;
                for (int si = 0; si < top; si++)
                {
                    swap_points(&(stack_kernel[0][si]), &(stack_kernel[1][si]), b);
                }
                swap_vertices(&(R[0]), &(R[1]), b);
                while(stack_scalar[top-1] < degree_t_num - jj)
                {
                    stack_scalar[top] = stack_scalar[top-1] + strategy_t[k];
                    xMULv2(&(stack_kernel[0][top]), &(stack_kernel[0][top-1]), &(strascalar_t[strascalar_i]), strascalar_t_bits[k], &(R[0].A24));
                    if (stack_scalar[top] < degree_t_num - jj) // This conditional check does NOT leak any secret information
                    {
                        xMULv2(&(stack_kernel[1][top]), &(stack_kernel[1][top-1]), &(strascalar_t[strascalar_i]), strascalar_t_bits[k], &(R[1].A24));
                    }
                    strascalar_i += strascalar_t_limbs[k];
                    top++;
                    k++;
                    assert(top <= stack_volume);
                }
                top--;
                assert(top >= 0);
                kps(i, stack_kernel[0][top], R[0].A24);
                xisog(&B24, i, R[0].A24);
                for(int traverse_stack = 0; traverse_stack < top; traverse_stack++)
                {
                    xeval(&(stack_kernel[0][traverse_stack]), i, stack_kernel[0][traverse_stack], R[0].A24);
                }
                xeval(&(R[0].Qs), i, R[0].Qs, R[0].A24);
                xeval(&(R[0].Qt), i, R[0].Qt, R[0].A24);
                copy_point(&(R[0].A24), &B24);
                kps_clear(i);
                for(int traverse_stack = 0; traverse_stack < top; traverse_stack++)
                {
                    xMULv2(&(stack_kernel[1][traverse_stack]), &(stack_kernel[1][traverse_stack]), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &(R[1].A24));
                }
                for (int si = 0; si < top; si++)
                {
                    swap_points(&(stack_kernel[0][si]), &(stack_kernel[1][si]), b);
                }
                swap_vertices(&(R[0]), &(R[1]), b);
            }
        }
        assert(top == 0);        

        #ifndef NDEBUG
        ec_curve_t E1, E2;
        A24_to_AC(&E1, &(R[0].A24));
        A24_to_AC(&E2, &(R[1].A24));
        fp2_t j1, j2;
        ec_j_inv(&j1, &E1);
        ec_j_inv(&j2, &E2);
        assert(fp2_is_equal(&j1, &j2));
        if(j < CYCLE_LENGTH - 1)
        {
            fp2_t t0, t1;
            fp2_copy(&t0, &(R[0].A24.z));
            fp2_inv(&t0);
            fp2_mul(&t0, &t0, &(R[0].A24.x));
            fp2_copy(&t1, &(R[1].A24.z));
            fp2_inv(&t1);
            fp2_mul(&t1, &t1, &(R[1].A24.x));
            assert(fp2_is_equal(&(t0), &(t1)));
        }
        #endif

        copy_point(&(R[0].Ps), &(R[1].Qs));
        copy_point(&(R[0].Pt), &(R[1].Qt));
        copy_vertices(&((*res)[j]), &R[0]);
    }
}

void eval_action_id(orient_cycle_t* res, const orient_cycle_t* inp, const int *s)
{
    assert(res != inp);
    vertex_t R[2];
    fp2_t mult_sq_total_X = fp2_1;
    fp2_t mult_sq_total_Z = fp2_1;
    
    for(int j = 0; j < CYCLE_LENGTH; j++)
    {
        copy_vertices(&(R[0]), &((*inp)[j]));
        copy_vertices_swap_kernel(&(R[1]), &((*inp)[j + 1 == CYCLE_LENGTH ? 0 : j + 1]));
        int mscalar_i, mscalar_bit_i;
        ec_point_t B24;

        // Isogenies with kernel in E[p+1]
        fp2_t mult_sq_X = Ms_mont;
        fp2_t mult_sq_Z = fp2_1;
        mscalar_i = 0;
        mscalar_bit_i = 0;
        for(int i = 0; i < P_LEN; i++)
        {
            for(int t = 1; t <= TORSION_ODD_POWERS[i]; t++)
            {
                digit_t b = compare_gt(t, s[i]);
                swap_vertices(&(R[0]), &(R[1]), b);
                ec_point_t K;
                xMULv2(&K, &((R[0]).Ps), &(mscalar_s[mscalar_i]), mscalar_s_bits[mscalar_bit_i], &((R[0]).A24));
                mscalar_i += mscalar_s_limbs[mscalar_bit_i];
                mscalar_bit_i++;
                kps(i, K, R[0].A24);
                xisog(&B24, i, R[0].A24);
                fp2_t pi_X = fp2_1, pi_Z = fp2_1;
                ec_ker_mul(&pi_Z, &pi_X, &K, &(R[0].A24), TORSION_ODD_PRIMES[i]);
                fp2_sqr(&pi_X, &pi_X);
                fp2_sqr(&pi_Z, &pi_Z);
                fp_mul(pi_Z.re, pi_Z.re, ell_mont[i]);
                fp_mul(pi_Z.im, pi_Z.im, ell_mont[i]);
                fp2_swap(&pi_X, &pi_Z, b);
                fp2_mul(&mult_sq_X, &mult_sq_X, &pi_X);
                fp2_mul(&mult_sq_Z, &mult_sq_Z, &pi_Z);
                xeval(&(R[0].Ps), i, R[0].Ps, R[0].A24);
                xeval(&(R[0].Pt), i, R[0].Pt, R[0].A24);
                xeval(&(R[0].Qs), i, R[0].Qs, R[0].A24);
                xeval(&(R[0].Qt), i, R[0].Qt, R[0].A24);
                copy_point(&(R[0].A24), &B24);
                kps_clear(i);
                xMULv2(&(R[1].Ps), &(R[1].Ps), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &(R[1].A24));
                swap_vertices(&(R[0]), &(R[1]), b);
            }
        }
        fp2_mul(&mult_sq_total_X, &mult_sq_total_X, &mult_sq_X);
        fp2_mul(&mult_sq_total_Z, &mult_sq_total_Z, &mult_sq_Z);

        // Isogenies with kernel in E[p-1]
        mult_sq_X = Mt_mont;
        mult_sq_Z = fp2_1;
        mscalar_i = 0;
        mscalar_bit_i = 0;
        for(int i = P_LEN; i < P_LEN + M_LEN; i++)
        {
            for(int t = 1; t <= TORSION_ODD_POWERS[i]; t++)
            {
                digit_t b = compare_gt(t, s[i]);
                swap_vertices(&(R[0]), &(R[1]), b);
                ec_point_t K;
                xMULv2(&K, &(R[0].Pt), &(mscalar_t[mscalar_i]), mscalar_t_bits[mscalar_bit_i], &(R[0].A24));
                mscalar_i += mscalar_t_limbs[mscalar_bit_i];
                mscalar_bit_i++;
                kps(i, K, R[0].A24);
                xisog(&B24, i, R[0].A24);
                fp2_t pi_X = fp2_1, pi_Z = fp2_1;
                ec_ker_mul(&pi_Z, &pi_X, &K, &(R[0].A24), TORSION_ODD_PRIMES[i]);
                fp2_sqr(&pi_X, &pi_X);
                fp2_sqr(&pi_Z, &pi_Z);
                fp_mul(pi_Z.re, pi_Z.re, ell_mont[i]);
                fp_mul(pi_Z.im, pi_Z.im, ell_mont[i]);
                fp2_swap(&pi_X, &pi_Z, b);
                fp2_mul(&mult_sq_X, &mult_sq_X, &pi_X);
                fp2_mul(&mult_sq_Z, &mult_sq_Z, &pi_Z);
                xeval(&(R[0].Pt), i, R[0].Pt, R[0].A24);
                xeval(&(R[0].Qs), i, R[0].Qs, R[0].A24);
                xeval(&(R[0].Qt), i, R[0].Qt, R[0].A24);
                copy_point(&(R[0].A24), &B24);
                kps_clear(i);
                xMULv2(&(R[1].Pt), &(R[1].Pt), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &(R[1].A24));
                swap_vertices(&(R[0]), &(R[1]), b);
            }
        }
        fp2_mul(&mult_sq_total_X, &mult_sq_total_X, &mult_sq_X);
        fp2_mul(&mult_sq_total_Z, &mult_sq_total_Z, &mult_sq_Z);
        

        #ifndef NDEBUG
        ec_curve_t E1, E2;
        A24_to_AC(&E1, &(R[0].A24));
        A24_to_AC(&E2, &(R[1].A24));
        fp2_t j1, j2;
        ec_j_inv(&j1, &E1);
        ec_j_inv(&j2, &E2);
        assert(fp2_is_equal(&j1, &j2));
        if(j < CYCLE_LENGTH - 1)
        {
            fp2_t t0, t1;
            fp2_copy(&t0, &(R[0].A24.z));
            fp2_inv(&t0);
            fp2_mul(&t0, &t0, &(R[0].A24.x));
            fp2_copy(&t1, &(R[1].A24.z));
            fp2_inv(&t1);
            fp2_mul(&t1, &t1, &(R[1].A24.x));
            assert(fp2_is_equal(&(t0), &(t1)));
        }
        #endif

        copy_point(&(R[0].Ps), &(R[1].Qs));
        copy_point(&(R[0].Pt), &(R[1].Qt));
        copy_vertices(&((*res)[j]), &(R[0]));
    }

    ec_curve_t A_left, A_right;
    A24_to_AC(&A_left, &(R[1].A24));
    A24_to_AC(&A_right, &(R[0].A24));
    fp2_t u_sq_X, u_sq_Z, t1;
    fp2_t r_X, r_Z;
    fp2_copy(&u_sq_X, &mult_sq_total_X);
    fp2_mul(&u_sq_Z, &alpha_sq, &mult_sq_total_Z);
    fp2_mul(&r_X, &u_sq_X, &(A_left.C));
    fp2_mul(&r_X, &r_X, &(A_right.A));
    fp2_mul(&t1, &u_sq_Z, &(A_left.A));
    fp2_mul(&t1, &t1, &(A_right.C));
    fp2_sub(&r_X, &r_X, &t1);
    fp2_mul(&r_Z, &u_sq_Z, &(A_left.C));
    fp2_mul(&r_Z, &r_Z, &(A_right.C));
    fp2_add(&t1, &r_Z, &r_Z);
    fp2_add(&r_Z, &r_Z, &t1);
    fp2_mul(&(R[0].Ps.x), &(R[0].Ps.x), &r_Z);
    fp2_mul(&t1, &r_X, &(R[0].Ps.z));
    fp2_sub(&(R[0].Ps.x), &(R[0].Ps.x), &t1);
    fp2_mul(&(R[0].Ps.x), &(R[0].Ps.x), &u_sq_Z);
    fp2_mul(&(R[0].Ps.z), &u_sq_X, &(R[0].Ps.z));
    fp2_mul(&(R[0].Ps.z), &(R[0].Ps.z), &r_Z);
    
    fp2_mul(&(R[0].Pt.x), &(R[0].Pt.x), &r_Z);
    fp2_mul(&t1, &r_X, &(R[0].Pt.z));
    fp2_sub(&(R[0].Pt.x), &(R[0].Pt.x), &t1);
    fp2_mul(&(R[0].Pt.x), &(R[0].Pt.x), &u_sq_Z);
    fp2_mul(&(R[0].Pt.z), &u_sq_X, &(R[0].Pt.z));
    fp2_mul(&(R[0].Pt.z), &(R[0].Pt.z), &r_Z);
    copy_point(&((*res)[CYCLE_LENGTH - 1].Ps), &(R[0].Ps));
    copy_point(&((*res)[CYCLE_LENGTH - 1].Pt), &(R[0].Pt));
}

void eval_action_id_stra(orient_cycle_t* res, const orient_cycle_t* inp, const int *s)
{
    assert(res != inp);
    vertex_t R[2];
    fp2_t mult_sq_total_X = fp2_1;
    fp2_t mult_sq_total_Z = fp2_1;
    
    for(int j = 0; j < CYCLE_LENGTH; j++)
    {
        copy_vertices(&(R[0]), &((*inp)[j]));
        copy_vertices_swap_kernel(&(R[1]), &((*inp)[j + 1 == CYCLE_LENGTH ? 0 : j + 1]));
        int strascalar_i/*, strascalar_bit_i*/;
        ec_point_t B24;

        // Isogenies with kernel in E[p+1]
        fp2_t mult_sq_X = Ms_mont;
        fp2_t mult_sq_Z = fp2_1;
        strascalar_i = 0;
        ec_point_t stack_kernel[2][stack_volume];
        int stack_scalar[stack_volume];
        size_t top = 0;
        copy_point(&(stack_kernel[0][top]), &(R[0].Ps));
        copy_point(&(stack_kernel[1][top]), &(R[1].Ps));
        stack_scalar[top] = 0;
        top++;
        int jj = 0, k = 0;
        for(int i = 0; i < P_LEN; i++)
        {
            for(int t = 1; t <= TORSION_ODD_POWERS[i]; t++)
            {
                digit_t b = compare_gt(t, s[i]);
                jj++;
                for (int si = 0; si < top; si++)
                {
                    swap_points(&(stack_kernel[0][si]), &(stack_kernel[1][si]), b);
                }
                swap_vertices(&(R[0]), &(R[1]), b);
                while(stack_scalar[top-1] < degree_s_num - jj)
                {
                    stack_scalar[top] = stack_scalar[top-1] + strategy_s[k];
                    xMULv2(&(stack_kernel[0][top]), &(stack_kernel[0][top-1]), &(strascalar_s[strascalar_i]), strascalar_s_bits[k], &(R[0].A24));
                    if (stack_scalar[top] < degree_s_num - jj) // This conditional check does NOT leak any secret information
                    {
                        xMULv2(&(stack_kernel[1][top]), &(stack_kernel[1][top-1]), &(strascalar_s[strascalar_i]), strascalar_s_bits[k], &(R[1].A24));
                    }
                    strascalar_i += strascalar_s_limbs[k];
                    top++;
                    k++;
                    assert(top <= stack_volume);
                }
                top--;
                assert(top >= 0);
                kps(i, stack_kernel[0][top], R[0].A24);
                xisog(&B24, i, R[0].A24);
                fp2_t pi_X = fp2_1, pi_Z = fp2_1;
                ec_ker_mul(&pi_Z, &pi_X, &(stack_kernel[0][top]), &(R[0].A24), TORSION_ODD_PRIMES[i]);
                fp2_sqr(&pi_X, &pi_X);
                fp2_sqr(&pi_Z, &pi_Z);
                fp_mul(pi_Z.re, pi_Z.re, ell_mont[i]);
                fp_mul(pi_Z.im, pi_Z.im, ell_mont[i]);
                fp2_swap(&pi_X, &pi_Z, b);
                fp2_mul(&mult_sq_X, &mult_sq_X, &pi_X);
                fp2_mul(&mult_sq_Z, &mult_sq_Z, &pi_Z);
                for(int traverse_stack = 0; traverse_stack < top; traverse_stack++)
                {
                    xeval(&(stack_kernel[0][traverse_stack]), i, stack_kernel[0][traverse_stack], R[0].A24);
                }
                xeval(&(R[0].Pt), i, R[0].Pt, R[0].A24);
                xeval(&(R[0].Qs), i, R[0].Qs, R[0].A24);
                xeval(&(R[0].Qt), i, R[0].Qt, R[0].A24);
                copy_point(&(R[0].A24), &B24);
                kps_clear(i);
                for(int traverse_stack = 0; traverse_stack < top; traverse_stack++)
                {
                    xMULv2(&(stack_kernel[1][traverse_stack]), &(stack_kernel[1][traverse_stack]), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &(R[1].A24));
                }
                for (int si = 0; si < top; si++)
                {
                    swap_points(&(stack_kernel[0][si]), &(stack_kernel[1][si]), b);
                }
                swap_vertices(&(R[0]), &(R[1]), b);
            }
        }
        assert(top == 0);
        fp2_mul(&mult_sq_total_X, &mult_sq_total_X, &mult_sq_X);
        fp2_mul(&mult_sq_total_Z, &mult_sq_total_Z, &mult_sq_Z);

        // Isogenies with kernel in E[p-1]
        mult_sq_X = Mt_mont;
        mult_sq_Z = fp2_1;
        strascalar_i = 0;
        top = 0;
        copy_point(&(stack_kernel[0][top]), &(R[0].Pt));
        copy_point(&(stack_kernel[1][top]), &(R[1].Pt));
        stack_scalar[top] = 0;
        top++;
        jj = 0; k = 0;
        for(int i = P_LEN; i < P_LEN + M_LEN; i++)
        {
            for(int t = 1; t <= TORSION_ODD_POWERS[i]; t++)
            {
                digit_t b = compare_gt(t, s[i]);
                jj++;
                for (int si = 0; si < top; si++)
                {
                    swap_points(&(stack_kernel[0][si]), &(stack_kernel[1][si]), b);
                }
                swap_vertices(&(R[0]), &(R[1]), b);
                while(stack_scalar[top-1] < degree_t_num - jj)
                {
                    stack_scalar[top] = stack_scalar[top-1] + strategy_t[k];
                    xMULv2(&(stack_kernel[0][top]), &(stack_kernel[0][top-1]), &(strascalar_t[strascalar_i]), strascalar_t_bits[k], &(R[0].A24));
                    if (stack_scalar[top] < degree_t_num - jj) // This conditional check does NOT leak any secret information
                    {
                        xMULv2(&(stack_kernel[1][top]), &(stack_kernel[1][top-1]), &(strascalar_t[strascalar_i]), strascalar_t_bits[k], &(R[1].A24));
                    }
                    strascalar_i += strascalar_t_limbs[k];
                    top++;
                    k++;
                    assert(top <= stack_volume);
                }
                top--;
                assert(top >= 0);
                kps(i, stack_kernel[0][top], R[0].A24);
                xisog(&B24, i, R[0].A24);
                fp2_t pi_X = fp2_1, pi_Z = fp2_1;
                ec_ker_mul(&pi_Z, &pi_X, &(stack_kernel[0][top]), &(R[0].A24), TORSION_ODD_PRIMES[i]);
                fp2_sqr(&pi_X, &pi_X);
                fp2_sqr(&pi_Z, &pi_Z);
                fp_mul(pi_Z.re, pi_Z.re, ell_mont[i]);
                fp_mul(pi_Z.im, pi_Z.im, ell_mont[i]);
                fp2_swap(&pi_X, &pi_Z, b);
                fp2_mul(&mult_sq_X, &mult_sq_X, &pi_X);
                fp2_mul(&mult_sq_Z, &mult_sq_Z, &pi_Z);
                for(int traverse_stack = 0; traverse_stack < top; traverse_stack++)
                {
                    xeval(&(stack_kernel[0][traverse_stack]), i, stack_kernel[0][traverse_stack], R[0].A24);
                }
                xeval(&(R[0].Qs), i, R[0].Qs, R[0].A24);
                xeval(&(R[0].Qt), i, R[0].Qt, R[0].A24);
                copy_point(&(R[0].A24), &B24);
                kps_clear(i);
                for(int traverse_stack = 0; traverse_stack < top; traverse_stack++)
                {
                    xMULv2(&(stack_kernel[1][traverse_stack]), &(stack_kernel[1][traverse_stack]), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &(R[1].A24));
                }
                for (int si = 0; si < top; si++)
                {
                    swap_points(&(stack_kernel[0][si]), &(stack_kernel[1][si]), b);
                }
                swap_vertices(&(R[0]), &(R[1]), b);
            }
        }
        assert(top == 0);
        fp2_mul(&mult_sq_total_X, &mult_sq_total_X, &mult_sq_X);
        fp2_mul(&mult_sq_total_Z, &mult_sq_total_Z, &mult_sq_Z);

        #ifndef NDEBUG
        ec_curve_t E1, E2;
        A24_to_AC(&E1, &(R[0].A24));
        A24_to_AC(&E2, &(R[1].A24));
        fp2_t j1, j2;
        ec_j_inv(&j1, &E1);
        ec_j_inv(&j2, &E2);
        assert(fp2_is_equal(&j1, &j2));
        if(j < CYCLE_LENGTH - 1)
        {
            fp2_t t0, t1;
            fp2_copy(&t0, &(R[0].A24.z));
            fp2_inv(&t0);
            fp2_mul(&t0, &t0, &(R[0].A24.x));
            fp2_copy(&t1, &(R[1].A24.z));
            fp2_inv(&t1);
            fp2_mul(&t1, &t1, &(R[1].A24.x));
            assert(fp2_is_equal(&(t0), &(t1)));
        }
        #endif

        copy_point(&(R[0].Ps), &(R[1].Qs));
        copy_point(&(R[0].Pt), &(R[1].Qt));
        copy_vertices(&((*res)[j]), &R[0]);
    }

    ec_curve_t A_left, A_right;
    A24_to_AC(&A_left, &(R[1].A24));
    A24_to_AC(&A_right, &(R[0].A24));
    fp2_t u_sq_X, u_sq_Z, t1;
    fp2_t r_X, r_Z;
    fp2_copy(&u_sq_X, &mult_sq_total_X);
    fp2_mul(&u_sq_Z, &alpha_sq, &mult_sq_total_Z);
    fp2_mul(&r_X, &u_sq_X, &(A_left.C));
    fp2_mul(&r_X, &r_X, &(A_right.A));
    fp2_mul(&t1, &u_sq_Z, &(A_left.A));
    fp2_mul(&t1, &t1, &(A_right.C));
    fp2_sub(&r_X, &r_X, &t1);
    fp2_mul(&r_Z, &u_sq_Z, &(A_left.C));
    fp2_mul(&r_Z, &r_Z, &(A_right.C));
    fp2_add(&t1, &r_Z, &r_Z);
    fp2_add(&r_Z, &r_Z, &t1);
    fp2_mul(&(R[0].Ps.x), &(R[0].Ps.x), &r_Z);
    fp2_mul(&t1, &r_X, &(R[0].Ps.z));
    fp2_sub(&(R[0].Ps.x), &(R[0].Ps.x), &t1);
    fp2_mul(&(R[0].Ps.x), &(R[0].Ps.x), &u_sq_Z);
    fp2_mul(&(R[0].Ps.z), &u_sq_X, &(R[0].Ps.z));
    fp2_mul(&(R[0].Ps.z), &(R[0].Ps.z), &r_Z);
    
    fp2_mul(&(R[0].Pt.x), &(R[0].Pt.x), &r_Z);
    fp2_mul(&t1, &r_X, &(R[0].Pt.z));
    fp2_sub(&(R[0].Pt.x), &(R[0].Pt.x), &t1);
    fp2_mul(&(R[0].Pt.x), &(R[0].Pt.x), &u_sq_Z);
    fp2_mul(&(R[0].Pt.z), &u_sq_X, &(R[0].Pt.z));
    fp2_mul(&(R[0].Pt.z), &(R[0].Pt.z), &r_Z);
    copy_point(&((*res)[CYCLE_LENGTH - 1].Ps), &(R[0].Ps));
    copy_point(&((*res)[CYCLE_LENGTH - 1].Pt), &(R[0].Pt));
}

void eval_action_omit_id(orient_cycle_t* res, const orient_cycle_t* inp, const int *s, int omit)
{
    assert(res != inp);
    assert(omit > 0);
    assert(omit < CYCLE_LENGTH);
    // vertex_t R0, R1;
    vertex_t R[2];
    
    for(int j = 0; j < CYCLE_LENGTH - omit; j++)
    {
        copy_vertices(&(R[0]), &((*inp)[j]));
        copy_vertices_swap_kernel(&(R[1]), &((*inp)[j + 1]));
        int mscalar_i, mscalar_bit_i;
        ec_point_t B24;

        // Isogenies with kernel in E[p+1]
        mscalar_i = 0;
        mscalar_bit_i = 0;
        for(int i = 0; i < P_LEN; i++)
        {
            for(int t = 1; t <= TORSION_ODD_POWERS[i]; t++)
            {
                digit_t b = compare_gt(t, s[i]);
                swap_vertices(&(R[0]), &(R[1]), b);
                ec_point_t K;
                xMULv2(&K, &(R[0].Ps), &(mscalar_s[mscalar_i]), mscalar_s_bits[mscalar_bit_i], &(R[0].A24));
                mscalar_i += mscalar_s_limbs[mscalar_bit_i];
                mscalar_bit_i++;
                kps(i, K, R[0].A24);
                xisog(&B24, i, R[0].A24);
                xeval(&(R[0].Ps), i, R[0].Ps, R[0].A24);
                xeval(&(R[0].Pt), i, R[0].Pt, R[0].A24);
                xeval(&(R[0].Qs), i, R[0].Qs, R[0].A24);
                xeval(&(R[0].Qt), i, R[0].Qt, R[0].A24);
                copy_point(&(R[0].A24), &B24);
                kps_clear(i);
                xMULv2(&(R[1].Ps), &(R[1].Ps), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &(R[1].A24));
                swap_vertices(&(R[0]), &(R[1]), b);
            }
        }

        // Isogenies with kernel in E[p-1]
        mscalar_i = 0;
        mscalar_bit_i = 0;
        for(int i = P_LEN; i < P_LEN + M_LEN; i++)
        {
            for(int t = 1; t <= TORSION_ODD_POWERS[i]; t++)
            {
                digit_t b = compare_gt(t, s[i]);
                swap_vertices(&(R[0]), &(R[1]), b);
                ec_point_t K;
                xMULv2(&K, &(R[0].Pt), &(mscalar_t[mscalar_i]), mscalar_t_bits[mscalar_bit_i], &(R[0].A24));
                mscalar_i += mscalar_t_limbs[mscalar_bit_i];
                mscalar_bit_i++;
                kps(i, K, R[0].A24);
                xisog(&B24, i, R[0].A24);
                xeval(&(R[0].Pt), i, R[0].Pt, R[0].A24);
                xeval(&(R[0].Qs), i, R[0].Qs, R[0].A24);
                xeval(&(R[0].Qt), i, R[0].Qt, R[0].A24);
                copy_point(&(R[0].A24), &B24);
                kps_clear(i);
                xMULv2(&(R[1].Pt), &(R[1].Pt), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &(R[1].A24));
                swap_vertices(&(R[0]), &(R[1]), b);
            }
        }
        

        #ifndef NDEBUG
        ec_curve_t E1, E2;
        A24_to_AC(&E1, &(R[0].A24));
        A24_to_AC(&E2, &(R[1].A24));
        fp2_t j1, j2;
        ec_j_inv(&j1, &E1);
        ec_j_inv(&j2, &E2);
        assert(fp2_is_equal(&j1, &j2));
        if(j < CYCLE_LENGTH - 1)
        {
            fp2_t t0, t1;
            fp2_copy(&t0, &(R[0].A24.z));
            fp2_inv(&t0);
            fp2_mul(&t0, &t0, &(R[0].A24.x));
            fp2_copy(&t1, &(R[1].A24.z));
            fp2_inv(&t1);
            fp2_mul(&t1, &t1, &(R[1].A24.x));
            assert(fp2_is_equal(&(t0), &(t1)));
        }
        #endif

        copy_point(&(R[0].Ps), &(R[1].Qs));
        copy_point(&(R[0].Pt), &(R[1].Qt));
        copy_vertices(&((*res)[j]), &(R[0]));
    }
}

void osidhld_PublicKey(orient_cycle_t* out, secretkey_t s)
{
    secretkey_t copysk;
    memcpy(copysk, s, sizeof(secretkey_t));
    orient_cycle_t cache[2];
    size_t source = 0, destination = 1;
    eval_action(&(cache[source]), &base_cycle, copysk);
    for(int j = 0; j < P_LEN + M_LEN; j++)
    {
        copysk[j] -= TORSION_ODD_POWERS[j];
    }
    for (int i = 0; i < ITERATIONS - 2; i++)
    {
        eval_action(&(cache[destination]), &(cache[source]), copysk);
        for(int j = 0; j < P_LEN + M_LEN; j++)
        {
            copysk[j] -= TORSION_ODD_POWERS[j];
        }
        source ^= 1;
        destination ^= 1;
    }
    eval_action(out, &cache[source], copysk);
}

void osidhld_PublicKey_stra(orient_cycle_t* out, secretkey_t s)
{
    secretkey_t copysk;
    memcpy(copysk, s, sizeof(secretkey_t));
    orient_cycle_t cache[2];
    size_t source = 0, destination = 1;
    eval_action_stra(&(cache[source]), &base_cycle, copysk);
    for(int j = 0; j < P_LEN + M_LEN; j++)
    {
        copysk[j] -= TORSION_ODD_POWERS[j];
    }
    for (int i = 0; i < ITERATIONS - 2; i++)
    {
        eval_action_stra(&(cache[destination]), &(cache[source]), copysk);
        for(int j = 0; j < P_LEN + M_LEN; j++)
        {
            copysk[j] -= TORSION_ODD_POWERS[j];
        }
        source ^= 1;
        destination ^= 1;
    }
    eval_action_stra(out, &cache[source], copysk);
}

void osidhld_SecretAgreement(fp2_t* out, const orient_cycle_t* in, secretkey_t s)
{
    secretkey_t copysk;
    memcpy(copysk, s, sizeof(secretkey_t));
    orient_cycle_t cache[2];
    size_t source = 0, destination = 1;
    eval_action(&(cache[source]), in, copysk);
    for(int j = 0; j < P_LEN + M_LEN; j++)
    {
        copysk[j] -= TORSION_ODD_POWERS[j];
    }
    for (int i = 0; i < ITERATIONS - CYCLE_LENGTH; i++)
    {
        eval_action(&(cache[destination]), &(cache[source]), copysk);
        for(int j = 0; j < P_LEN + M_LEN; j++)
        {
            copysk[j] -= TORSION_ODD_POWERS[j];
        }
        source ^= 1;
        destination ^= 1;
    }
    for (int i = 1; i < CYCLE_LENGTH - 1; i++)
    {
        eval_action_omit(&(cache[destination]), &(cache[source]), copysk, i);
        for(int j = 0; j < P_LEN + M_LEN; j++)
        {
            copysk[j] -= TORSION_ODD_POWERS[j];
        }
        source ^= 1;
        destination ^= 1;
    }
    eval_action_omit(&(cache[destination]), &cache[source], copysk, CYCLE_LENGTH - 1);
    ec_curve_t result;
    A24_to_AC(&result, &(cache[destination][0].A24));
    ec_j_inv(out, &result);
}

void osidhld_SecretAgreement_stra(fp2_t* out, const orient_cycle_t* in, secretkey_t s)
{
    secretkey_t copysk;
    memcpy(copysk, s, sizeof(secretkey_t));
    orient_cycle_t cache[2];
    size_t source = 0, destination = 1;
    eval_action_stra(&(cache[source]), in, copysk);
    for(int j = 0; j < P_LEN + M_LEN; j++)
    {
        copysk[j] -= TORSION_ODD_POWERS[j];
    }
    for (int i = 0; i < ITERATIONS - CYCLE_LENGTH; i++)
    {
        eval_action_stra(&(cache[destination]), &(cache[source]), copysk);
        for(int j = 0; j < P_LEN + M_LEN; j++)
        {
            copysk[j] -= TORSION_ODD_POWERS[j];
        }
        source ^= 1;
        destination ^= 1;
    }
    for (int i = 1; i < CYCLE_LENGTH - 1; i++)
    {
        eval_action_stra_omit(&(cache[destination]), &(cache[source]), copysk, i);
        for(int j = 0; j < P_LEN + M_LEN; j++)
        {
            copysk[j] -= TORSION_ODD_POWERS[j];
        }
        source ^= 1;
        destination ^= 1;
    }
    eval_action_stra_omit(&(cache[destination]), &cache[source], copysk, CYCLE_LENGTH - 1);
    ec_curve_t result;
    A24_to_AC(&result, &(cache[destination][0].A24));
    ec_j_inv(out, &result);
}

void osidhld_PublicKey_id(orient_cycle_t* out, secretkey_t s)
{
    secretkey_t copysk;
    memcpy(copysk, s, sizeof(secretkey_t));
    orient_cycle_t cache[2];
    size_t source = 0, destination = 1;
    eval_action_id(&(cache[source]), &base_cycle, copysk);
    for(int j = 0; j < P_LEN + M_LEN; j++)
    {
        copysk[j] -= TORSION_ODD_POWERS[j];
    }
    for (int i = 0; i < ITERATIONS - 2; i++)
    {
        eval_action_id(&(cache[destination]), &(cache[source]), copysk);
        for(int j = 0; j < P_LEN + M_LEN; j++)
        {
            copysk[j] -= TORSION_ODD_POWERS[j];
        }
        source ^= 1;
        destination ^= 1;
    }
    eval_action_id(out, &cache[source], copysk);
}

void osidhld_SecretAgreement_id(fp2_t* out, const orient_cycle_t* in, secretkey_t s)
{
    secretkey_t copysk;
    memcpy(copysk, s, sizeof(secretkey_t));
    orient_cycle_t cache[2];
    size_t source = 0, destination = 1;
    eval_action_id(&(cache[source]), in, copysk);
    for(int j = 0; j < P_LEN + M_LEN; j++)
    {
        copysk[j] -= TORSION_ODD_POWERS[j];
    }
    for (int i = 0; i < ITERATIONS - CYCLE_LENGTH; i++)
    {
        eval_action_id(&(cache[destination]), &(cache[source]), copysk);
        for(int j = 0; j < P_LEN + M_LEN; j++)
        {
            copysk[j] -= TORSION_ODD_POWERS[j];
        }
        source ^= 1;
        destination ^= 1;
    }
    for (int i = 1; i < CYCLE_LENGTH - 1; i++)
    {
        eval_action_omit(&(cache[destination]), &(cache[source]), copysk, i);
        for(int j = 0; j < P_LEN + M_LEN; j++)
        {
            copysk[j] -= TORSION_ODD_POWERS[j];
        }
        source ^= 1;
        destination ^= 1;
    }
    eval_action_omit(&(cache[destination]), &cache[source], copysk, CYCLE_LENGTH - 1);
    ec_curve_t result;
    A24_to_AC(&result, &(cache[destination][0].A24));
    ec_j_inv(out, &result);
}

void osidhld_PublicKey_id_stra(orient_cycle_t* out, secretkey_t s)
{
    secretkey_t copysk;
    memcpy(copysk, s, sizeof(secretkey_t));
    orient_cycle_t cache[2];
    size_t source = 0, destination = 1;
    eval_action_id_stra(&(cache[source]), &base_cycle, copysk);
    for(int j = 0; j < P_LEN + M_LEN; j++)
    {
        copysk[j] -= TORSION_ODD_POWERS[j];
    }
    for (int i = 0; i < ITERATIONS - 2; i++)
    {
        eval_action_id_stra(&(cache[destination]), &(cache[source]), copysk);
        for(int j = 0; j < P_LEN + M_LEN; j++)
        {
            copysk[j] -= TORSION_ODD_POWERS[j];
        }
        source ^= 1;
        destination ^= 1;
    }
    eval_action_id_stra(out, &cache[source], copysk);
}

void osidhld_SecretAgreement_id_stra(fp2_t* out, const orient_cycle_t* in, secretkey_t s)
{
    secretkey_t copysk;
    memcpy(copysk, s, sizeof(secretkey_t));
    orient_cycle_t cache[2];
    size_t source = 0, destination = 1;
    eval_action_id_stra(&(cache[source]), in, copysk);
    for(int j = 0; j < P_LEN + M_LEN; j++)
    {
        copysk[j] -= TORSION_ODD_POWERS[j];
    }
    for (int i = 0; i < ITERATIONS - CYCLE_LENGTH; i++)
    {
        eval_action_id_stra(&(cache[destination]), &(cache[source]), copysk);
        for(int j = 0; j < P_LEN + M_LEN; j++)
        {
            copysk[j] -= TORSION_ODD_POWERS[j];
        }
        source ^= 1;
        destination ^= 1;
    }
    for (int i = 1; i < CYCLE_LENGTH - 1; i++)
    {
        eval_action_stra_omit(&(cache[destination]), &(cache[source]), copysk, i);
        for(int j = 0; j < P_LEN + M_LEN; j++)
        {
            copysk[j] -= TORSION_ODD_POWERS[j];
        }
        source ^= 1;
        destination ^= 1;
    }
    eval_action_stra_omit(&(cache[destination]), &cache[source], copysk, CYCLE_LENGTH - 1);
    ec_curve_t result;
    A24_to_AC(&result, &(cache[destination][0].A24));
    ec_j_inv(out, &result);
}
