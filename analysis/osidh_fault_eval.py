#!/usr/bin/env python3
"""
osidh_fault_eval.py
===================
Evaluation / attack-counting harness for the disorientation fault against the
OSIDH-LD reference implementation of

    Wang, Lai, Lin, Zhao, "Efficient and Parallel Implementation of
    Isogeny-based Deterministic Group Actions", TCHES 2026(3), pp. 541-565.

Attack idea ("From if to Oops", Banegas-Trimoska): every eval_action variant
routes its single secret-dependent decision through

    digit_t b = compare_gt(t, s[i]);   //  b = 0xFF..FF  iff  t <= s[i]  (swap)
    swap_vertices(&R[0], &R[1], b);     //  pick the trail for this l_i-isogeny
        ... one l_i-isogeny step on R[0] ...
    swap_vertices(&R[0], &R[1], b);     //  swap back

Skipping compare_gt (or forcing b=0) on a step where t <= s[i] sends that one
l_i-isogeny down the WRONG orientation trail.  The faulty public key differs
from the honest one <=> the fault was "effective".  That effective/no-op bit is
a comparison oracle on the secret.

Encoding recovered from protocols.c + the driver loop:
  * secret  sk[i] in [0, ITERATIONS * T_i]   (make_SecretKey + mask + reject)
  * driver calls eval_action ITERATIONS times, doing  copysk[i] -= T_i  each round
  * inside a round the swap at (round r, prime i, step t) is genuine
        iff   t <= sk[i] - r*T_i   iff   sk[i] >= r*T_i + t.
  * as (r,t) range over 0<=r<ITERATIONS, 1<=t<=T_i, the threshold r*T_i+t
    takes every value in [1, N_i], N_i = ITERATIONS*T_i.
  => a single effective fault at the slot with threshold tau tests  sk[i] >= tau.
     Binary search over tau recovers sk[i] in ceil(log2(N_i+1)) effective faults.

This harness does NOT run the 256-bit Fp2 group action (impossible on the M4;
that is why target_osidh.c is a proxy).  It reproduces the *logical* fault
oracle exactly from the real parameter set and answers two questions:
  (A) idealised: how many EFFECTIVE faults recover the full key, and does the
      recovery reproduce the secret exactly (validation)?
  (B) realistic: given a per-attempt fault probability p (calibrated to the
      ChipWhisperer campaign), how many RAW glitch attempts does that cost?
"""

import argparse
import math
import random
from statistics import mean, pstdev

# ---------------------------------------------------------------------------
# Real parameters  (precomp/ref/ld256 + protocol_setup.h + ec_params.h)
# ---------------------------------------------------------------------------
ITERATIONS   = 35
CYCLE_LENGTH = 13
P_LEN        = 31
M_LEN        = 18
N_PRIMES     = P_LEN + M_LEN            # 49

# index-aligned layout (DIGIT_LEN 16/32 in torsion_constants.c; the 64-bit
# layout only permutes the primes, the multiset of powers is identical)
TORSION_ODD_PRIMES = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157,
    163, 167, 173, 179, 181, 191, 193, 199, 211, 229, 239, 241, 269]
TORSION_ODD_POWERS = [
    4, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

assert len(TORSION_ODD_PRIMES) == N_PRIMES == len(TORSION_ODD_POWERS)

# sk[i] in [0, N_i]
N = [ITERATIONS * T for T in TORSION_ODD_POWERS]          # per-prime max exponent
BITS = [math.ceil(math.log2(n + 1)) for n in N]           # binary-search depth


# ---------------------------------------------------------------------------
# Ground truth: a secret key, generated exactly like make_SecretKey()
# ---------------------------------------------------------------------------
def gen_secret(rng):
    # mask+reject in the C reduces to a uniform draw on [0, N_i]
    return [rng.randint(0, n) for n in N]


# ---------------------------------------------------------------------------
# (A) IDEALISED perfect oracle
#     oracle(i, tau) := (sk[i] >= tau)   -- one effective-vs-no-op fault
#     Count every oracle call as one effective fault query.
# ---------------------------------------------------------------------------
def bsearch_decisions(sk_i, n_i):
    """Binary search for sk_i in [0, n_i].
    Returns (recovered, n_ge, n_lt): #effective (output-changing) and
    #no-op (output-identical) faulted executions used."""
    lo, hi = 0, n_i
    n_ge = n_lt = 0
    while lo < hi:
        mid = (lo + hi + 1) // 2          # test threshold tau = mid in [1, n_i]
        if sk_i >= mid:                   # effective fault observed
            lo = mid
            n_ge += 1
        else:                             # no-op fault
            hi = mid - 1
            n_lt += 1
    return lo, n_ge, n_lt


def run_ideal(rng, trials):
    ok = 0
    tot_ge = tot_lt = 0
    for _ in range(trials):
        sk = gen_secret(rng)
        rec = []
        for i in range(N_PRIMES):
            v, g, l = bsearch_decisions(sk[i], N[i])
            rec.append(v)
            tot_ge += g
            tot_lt += l
        ok += (rec == sk)
    return ok, tot_ge / trials, tot_lt / trials


# ---------------------------------------------------------------------------
# (B) REALISTIC noisy oracle  (exact expectation, no sampling)
#     A physical glitch on a genuinely-swappable slot (sk[i] >= tau) is
#     effective with probability p; a non-swappable slot (sk[i] < tau) can
#     never become effective, so it always looks like a no-op.
#     Decision rule per threshold: glitch up to n_reps times.
#        - >=1 output change  -> "sk[i] >= tau"  (definitive, stop early)
#        - n_reps no-ops      -> "sk[i] <  tau"  (wrong for a true ">=" only if
#                                the glitch never fired: prob eps=(1-p)^n_reps)
#     Expected raw glitches:
#        true ">=" threshold : sum_{k=1..n} (1-p)^{k-1} = (1-(1-p)^n)/p
#        true "<"  threshold : n_reps  (always all reps, no early stop)
#     Averages over sk_i are EXACT: N_i is tiny, so enumerate every value.
# ---------------------------------------------------------------------------
def prime_decision_profile(n_i):
    """Average (#ge, #lt) binary-search decisions over sk_i uniform in [0,n_i]."""
    ge = lt = 0
    for v in range(n_i + 1):
        _, g, l = bsearch_decisions(v, n_i)
        ge += g
        lt += l
    m = n_i + 1
    return ge / m, lt / m


def expected_raw(p, n_reps, profiles):
    cost_ge = (1 - (1 - p) ** n_reps) / p        # geometric, capped at n_reps
    cost_lt = n_reps
    exp_ge = sum(g for g, _ in profiles)
    exp_lt = sum(l for _, l in profiles)
    return exp_ge * cost_ge + exp_lt * cost_lt, exp_ge, exp_lt


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--trials', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--p', type=float, nargs='+',
                    default=[0.00234, 0.01, 0.05, 0.10, 0.25],
                    help='per-attempt effective-fault probabilities to sweep')
    ap.add_argument('--target-error', type=float, default=2**-20,
                    help='max tolerated prob of a wrong "<" per threshold')
    args = ap.parse_args()
    rng = random.Random(args.seed)

    # ---- static structure ------------------------------------------------
    total_bits_entropy = sum(math.log2(n + 1) for n in N)
    total_bs_depth = sum(BITS)
    print("=" * 70)
    print("OSIDH-LD (TCHES 2026) disorientation-fault attack — evaluation")
    print("=" * 70)
    print(f"primes (P_LEN+M_LEN)          : {N_PRIMES}   (l = 3..269)")
    print(f"ITERATIONS / CYCLE_LENGTH     : {ITERATIONS} / {CYCLE_LENGTH}")
    print(f"per-prime exponent range N_i  : sk[i] in [0, ITERATIONS*T_i]")
    specials = [(TORSION_ODD_PRIMES[i], TORSION_ODD_POWERS[i], N[i], BITS[i])
                for i in range(N_PRIMES) if TORSION_ODD_POWERS[i] > 1]
    print(f"  T_i>1 primes                : "
          + ", ".join(f"l={l}:T={t},N={n},{b}b" for l, t, n, b in specials))
    print(f"  all other 47 primes         : T=1, N=35, {math.ceil(math.log2(36))}b each")
    print(f"secret entropy (sum log2(N_i+1)): {total_bits_entropy:.1f} bits"
          f"   (~NIST-I / 256-bit prime)")
    print(f"binary-search depth (sum bits) : {total_bs_depth}"
          f"   <- minimum effective faults, full key")
    print()

    # ---- (A) idealised ---------------------------------------------------
    ok, avg_ge, avg_lt = run_ideal(random.Random(args.seed), args.trials)
    avg_q = avg_ge + avg_lt
    print("-" * 70)
    print(f"(A) IDEALISED perfect oracle  [{args.trials} random keys]")
    print("-" * 70)
    print(f"exact key recoveries          : {ok}/{args.trials}  "
          f"({'VALIDATED' if ok == args.trials else 'FAILURES!'})")
    print(f"faulted executions / full key : mean {avg_q:.1f}  "
          f"(<= worst case sum ceil(log2(N_i+1)) = {total_bs_depth})")
    print(f"  of which effective (output changed): {avg_ge:.1f}")
    print(f"  of which no-op    (output identical): {avg_lt:.1f}")
    print(f"  = one faulted execution per {total_bits_entropy/avg_q:.2f} secret bits")
    print()

    # ---- (B) realistic (exact expectation) -------------------------------
    profiles = [prime_decision_profile(N[i]) for i in range(N_PRIMES)]
    print("-" * 70)
    print("(B) REALISTIC noisy oracle — expected RAW glitch attempts")
    print(f"    n_reps set so per-threshold miss eps=(1-p)^n_reps <= "
          f"{args.target_error:.1e}")
    print("-" * 70)
    exp_ge = sum(g for g, _ in profiles)
    exp_lt = sum(l for _, l in profiles)
    print(f"expected effective decisions/key : {exp_ge:.1f}   "
          f"no-op decisions/key : {exp_lt:.1f}")
    print()
    print(f"{'p':>8} | {'n_reps':>7} | {'raw glitches (exp)':>18} | "
          f"{'per-key err <=':>14}")
    for p in args.p:
        n_reps = max(1, math.ceil(math.log(args.target_error) / math.log(1 - p)))
        raw, eg, el = expected_raw(p, n_reps, profiles)
        eps = (1 - p) ** n_reps
        key_err = 1 - (1 - eps) ** eg      # any effective decision missed
        print(f"{p:>8.5f} | {n_reps:>7d} | {raw:>18,.0f} | {key_err:>14.1e}")
    print()
    print("Calibration: the ChipWhisperer campaign in the paper reports 3")
    print("effective skips in 1281 tries (p ~ 0.00234) against Houben; OSIDH-LD")
    print("exposes CYCLE_LENGTH*ITERATIONS = 455 physical occurrences of each")
    print("slot per action, so the effective p per action is expected higher.")
    print("'raw glitches' counts attempts on the fault-sensitive function only;")
    print("device resets / re-syncs add a constant 1/(1-p_reset) overhead.")


if __name__ == '__main__':
    main()
