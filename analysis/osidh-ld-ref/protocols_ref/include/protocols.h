// SPDX-FileCopyrightText: 2026 The project author
// SPDX-License-Identifier: Apache-2.0

#ifndef PROTOCOLS_H
#define PROTOCOLS_H

#include "ec.h"
#include "isog.h"

#include "protocol_setup.h"

typedef struct vertex_t {
    ec_point_t A24;
    ec_point_t Ps;
    ec_point_t Pt;
    ec_point_t Qs;
    ec_point_t Qt;
} vertex_t;

typedef vertex_t orient_cycle_t[CYCLE_LENGTH];

typedef int secretkey_t[P_LEN + M_LEN];

void make_SecretKey(secretkey_t sk);
void osidhld_PublicKey(orient_cycle_t* out, secretkey_t s);
void osidhld_PublicKey_stra(orient_cycle_t* out, secretkey_t s);
void osidhld_SecretAgreement(fp2_t* out, const orient_cycle_t* in, secretkey_t s);
void osidhld_SecretAgreement_stra(fp2_t* out, const orient_cycle_t* in, secretkey_t s);

void print_SecretKey(const secretkey_t sk);
void print_PublicKey(const orient_cycle_t *pk);
void print_SecretAgreement(const fp2_t a);

void osidhld_PublicKey_id(orient_cycle_t* out, secretkey_t s);
void osidhld_SecretAgreement_id(fp2_t* out, const orient_cycle_t* in, secretkey_t s);
void osidhld_PublicKey_id_stra(orient_cycle_t* out, secretkey_t s);
void osidhld_SecretAgreement_id_stra(fp2_t* out, const orient_cycle_t* in, secretkey_t s);

#endif