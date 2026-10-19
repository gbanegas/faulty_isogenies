/*
 * fault_oracle_check.c
 * Real single-fault validation against the OSIDH-LD *reference* eval_action.
 *
 * We compute one honest orient_cycle = eval_action(base_cycle, s), then repeat
 * it while disorienting exactly one compare_gt call (index k within the run),
 * and compare the outputs. Prediction: the fault is EFFECTIVE (output changes)
 * iff the targeted (prime i, step t) satisfies t <= s[i] (a genuine swap).
 *
 * This exercises the actual field/isogeny code, not a model — it confirms the
 * effective/no-op oracle the counting harness assumes really exists.
 */
#include <stdio.h>
#include <string.h>
#include "protocols.h"
#include "protocols_helper.h"
#include "protocol_constants.h"

// fault-instrumentation globals (declared extern in protocols_helper.h)
long g_fault_skip = -1;
long g_fault_ctr  = 0;

extern void eval_action(orient_cycle_t* res, const orient_cycle_t* inp, const int *s);

// Map a compare_gt call index within the FIRST cycle slot (j=0) to (prime i, step t).
// eval_action order for one j: block E[p+1] i=0..P_LEN-1 then E[p-1] i=P_LEN..,
// each prime i contributing t=1..TORSION_ODD_POWERS[i] calls.
static int index_to_it(long k, int *pi, int *pt)
{
    long idx = 0;
    for (int i = 0; i < P_LEN + M_LEN; i++) {
        for (int t = 1; t <= (int)TORSION_ODD_POWERS[i]; t++) {
            if (idx == k) { *pi = i; *pt = t; return 1; }
            idx++;
        }
    }
    return 0; // k beyond one j-slice
}

static int run_faulted(orient_cycle_t *out, const int *s, long skip)
{
    g_fault_ctr = 0;
    g_fault_skip = skip;
    eval_action(out, &base_cycle, s);
    g_fault_skip = -1;
    return 0;
}

int main(void)
{
    // A fixed, in-range secret. Interesting choices at prime 0 (l=3, T=4):
    //   s[0] = 2  -> steps t=1,2 are genuine swaps; t=3,4 are no-ops.
    int s[P_LEN + M_LEN];
    for (int i = 0; i < P_LEN + M_LEN; i++) s[i] = 1;         // most primes: swap at t=1
    s[0] = 2;    // l=3: genuine for t<=2, no-op for t=3,4
    s[3] = 1;    // l=11 (T=2): genuine for t=1, no-op for t=2
    s[5] = 0;    // l=17: no swap at all (t=1 is a no-op)

    orient_cycle_t honest, faulted;
    run_faulted(&honest, s, -1);   // no fault

    printf("call |  l  | i  | t | s[i] | predict | observed | match\n");
    printf("-----+-----+----+---+------+---------+----------+------\n");
    int ok = 1, tested = 0;
    // probe the whole j=0 slice (53 calls) — every fault site of the first cycle
    for (long k = 0; ; k++) {
        int i, t;
        if (!index_to_it(k, &i, &t)) break;
        run_faulted(&faulted, s, k);
        int changed = (memcmp(&honest, &faulted, sizeof(orient_cycle_t)) != 0);
        int predict_eff = (t <= s[i]);      // genuine swap disoriented?
        int good = (changed == predict_eff);
        ok &= good; tested++;
        // print a few representative rows (prime 0, 3, 5 and a couple others)
        if (i <= 5 || t > 1) {
            printf("%4ld | %3llu | %2d | %d | %4d | %7s | %8s | %s\n",
                   k, (unsigned long long)TORSION_ODD_PRIMES[i], i, t, s[i],
                   predict_eff ? "EFFECT" : "no-op",
                   changed ? "changed" : "same",
                   good ? "ok" : "MISMATCH");
        }
    }
    printf("-----+-----+----+---+------+---------+----------+------\n");
    printf("tested %d fault sites (full j=0 cycle slice); all match prediction: %s\n",
           tested, ok ? "YES" : "NO");
    return ok ? 0 : 1;
}
