#!/usr/bin/env python3
"""
Toy OSIDH-LD with real isogenies, and the fault attack of
"From if to Oops: Turning a conditional into a key-recovery oracle".

What is implemented
  * Algorithm 2 (one iteration of the OSIDH-LD group action) step by step: a cycle of r
    tuples (E_j, P_j, Q_j); for every position j two lanes R0 = (E_j, P_j, Q_j) and
    R1 = (E_{j+1}, Q_{j+1}, P_{j+1}); for every prime l_i and inner step t = 1..e_i the
    mask b = compare_gt(t, s_eff) selects the lane by a swap before and after the step;
    kernel bookkeeping with the cofactor m; R0[1] <- R1[2] at the end; a lane check
    (the toy's IsoFix) that the two lanes reached the same curve.
  * The full action: B iterations with s_eff = clamp(s_i - u*e_i, 0, e_i), so s_i in [0, B*e_i].
  * The fault: the mask at one chosen (u, j, i, t) is forced to 0.
  * The attack of Sec. 4.3: two-stage binary search per component, with the oracle
    "faulted public key == honest public key" (an aborted run reads as a mismatch).

How it differs from the real OSIDH-LD (it is a toy)
  * The orientation is the Frobenius one over a small F_p (p = 4*M - 1, p = 3 mod 8): curves
    are Montgomery curves E_A over F_p, P_j is an F_p-point and Q_j a twist point of order
    M = prod l_i^e_i. The cycle closes because [prod l_i^e_i] has order r = 3 in this setting.
    The real scheme uses F_p^2 curves, a non-Frobenius orientation and r = 13.
  * Arithmetic is affine x-only (Velu formulas of Costello-Hisil) and not constant time.
  * The torsion powers mirror ld256: e = 4 for l = 3, e = 2 for l = 11, e = 1 otherwise.

Every result is cross-checked against an independent textbook CSIDH action that uses
random point sampling. Pure Python 3.8+, no dependencies; about a minute with defaults.
"""
import argparse, math, random, time

# ================================================================ parameters
ELLS = [3, 5, 7, 11, 13, 17, 23]
ES = [4, 1, 1, 2, 1, 1, 1]
M = math.prod(l**e for l, e in zip(ELLS, ES))
P = 4 * M - 1                                     # 6974587619, prime, = 3 mod 8
N = len(ELLS)
REF_RNG = random.Random(time.time())

# ================================================================ F_p and Montgomery x-only arithmetic
def inv(a):
    return pow(a, P - 2, P)

def is_square(a):
    return pow(a, (P - 1) // 2, P) == 1

def rhs(A, x):
    return (x * x * x + A * x * x + x) % P

def xdbl(X, Z, a24):
    t0 = (X + Z) * (X + Z) % P
    t1 = (X - Z) * (X - Z) % P
    t2 = (t0 - t1) % P
    return t0 * t1 % P, t2 * (t1 + a24 * t2) % P

def xadd(XP, ZP, XQ, ZQ, XD, ZD):
    u = (XP - ZP) * (XQ + ZQ) % P
    v = (XP + ZP) * (XQ - ZQ) % P
    return ZD * (u + v) * (u + v) % P, XD * (u - v) * (u - v) % P

def mul(A, k, x):
    """[k]x on E_A or its twist (x-only Montgomery ladder); None is the point at infinity."""
    if x is None or k == 0:
        return None
    a24 = (A + 2) * inv(4) % P
    X0, Z0, X1, Z1 = 1, 0, x, 1
    for bit in bin(k)[2:]:
        if bit == "1":
            X0, Z0 = xadd(X0, Z0, X1, Z1, x, 1)
            X1, Z1 = xdbl(X1, Z1, a24)
        else:
            X1, Z1 = xadd(X0, Z0, X1, Z1, x, 1)
            X0, Z0 = xdbl(X0, Z0, a24)
    return None if Z0 == 0 else X0 * inv(Z0) % P

def isogeny(A, xK, ell, pts):
    """odd-degree x-only Velu: codomain A' = pi^2 (A - 6 sum(x_i - 1/x_i)) and images of pts.
    A kernel point on E_A gives a positive step, one on the twist a negative step."""
    xs = [mul(A, k, xK) for k in range(1, (ell - 1) // 2 + 1)]
    assert None not in xs and mul(A, ell, xK) is None, "kernel point of wrong order"
    pi, sig = 1, 0
    for xi in xs:
        pi = pi * xi % P
        sig = (sig + xi - inv(xi)) % P
    A2 = pi * pi * (A - 6 * sig) % P
    out = []
    for x in pts:
        if x is None or x in xs:
            out.append(None)
            continue
        num, den = 1, 1
        for xi in xs:
            num = num * (x * xi - 1) % P
            den = den * (x - xi) % P
        out.append(x * num * num % P * inv(den * den % P) % P)
    return A2, out

# ================================================================ reference CSIDH action (independent)


def random_x(A, side, rng=REF_RNG):
    while True:
        x = rng.randrange(1, P)
        if is_square(rhs(A, x)) == (side > 0):
            return x

def reference_action(A, exps):
    """textbook CSIDH action prod [l_i]^exps_i * E_A with random point sampling."""
    e = list(exps)
    while any(e):
        for side in (1, -1):
            idx = [i for i in range(N) if e[i] * side > 0]
            if not idx:
                continue
            k = math.prod(ELLS[i] for i in idx)
            Pt = mul(A, (P + 1) // k, random_x(A, side))
            for i in idx:
                k //= ELLS[i]
                K = mul(A, k, Pt)
                if K is None:
                    continue
                A, (Pt,) = isogeny(A, K, ELLS[i], [Pt])
                e[i] -= side
    return A

def point_of_order_M(A, side, rng):
    """an F_p-point (side=+1) or twist point (side=-1) of exact order M."""
    while True:
        x = mul(A, (P + 1) // M, random_x(A, side, rng))
        if x is not None and all(mul(A, M // l, x) is not None for l in ELLS):
            return x

# ================================================================ Algorithm 2 with a fault hook
class LaneCheckFailure(Exception):
    def __init__(self, u, positions):
        super().__init__(f"lane check failed in iteration {u} at positions {positions}")
        self.u, self.positions = u, positions

def iteration(cycle, seff, u, fault=None):
    """One application of Algorithm 2. cycle = [(A_j, xP_j, xQ_j)], seff[i] in [0, e_i].
    fault = (u, j, i, t, kind) with kind 'mask' (b forced to 0), 'swap1' or 'swap2'
    (first or second conditional swap skipped). Returns (new cycle, failing positions)."""
    r = len(cycle)
    out, failed = [], []
    for j in range(r):
        R0 = list(cycle[j])                                       # (E_j, P_j, Q_j)
        A1, P1, Q1 = cycle[(j + 1) % r]
        R1 = [A1, Q1, P1]                                         # (E_{j+1}, Q_{j+1}, P_{j+1})
        m = M
        for i in range(N):
            for t in range(1, ES[i] + 1):
                b = 1 if t > seff[i] else 0                       # b <- compare_gt(t, s_i)
                kind = fault[4] if fault and fault[:4] == (u, j, i, t) else None
                if kind == "mask":
                    b = 0                                         # the fault
                if b and kind != "swap1":
                    R0, R1 = R1, R0                               # Cswap(R0, R1, b)
                m //= ELLS[i]
                K = mul(R0[0], m, R0[1])                          # K <- Mul(R0[0], R0[1], m)
                R0[0], (R0[1], R0[2]) = isogeny(R0[0], K, ELLS[i], [R0[1], R0[2]])
                R1[1] = mul(R1[0], ELLS[i], R1[1])                # R1[1] <- Mul(R1[0], R1[1], l_i)
                if b and kind != "swap2":
                    R0, R1 = R1, R0                               # Cswap(R0, R1, b)
        R0[1] = R1[2]                                             # R0[1] <- R1[2]
        if R0[0] != R1[0]:                                        # IsoFix: lanes must meet
            failed.append(j)
        out.append(tuple(R0))
    return out, failed

def group_action(cycle, s, fault=None):
    """B iterations of Algorithm 2; raises LaneCheckFailure like a build with assertions."""
    for u in range(B):
        seff = [min(max(s[i] - u * ES[i], 0), ES[i]) for i in range(N)]
        cycle, failed = iteration(cycle, seff, u, fault)
        if failed:
            raise LaneCheckFailure(u, failed)
    return cycle

def curves(cycle):
    return [A for A, _, _ in cycle]

# ================================================================ the attack (Sec. 4.3)
def oracle(base, honest, s_victim, u, i, t, j=0, stats=None):
    """1 if the faulted public key equals the honest one (no-op), 0 otherwise (abort = 0)."""
    if stats is not None:
        stats["queries"] += 1
    try:
        out = group_action(base, s_victim, fault=(u, j, i, t, "mask"))
    except LaneCheckFailure:
        return 0
    return 1 if curves(out) == honest else 0

def recover_component(base, honest, s_victim, i, stats):
    """stage 1: smallest u* in 0..B with O(u*, i, 1) = 0; stage 2 (e_i > 1): largest t with O(u*-1, i, t) = 1."""
    lo, hi = 0, B
    while lo < hi:
        mid = (lo + hi) // 2
        if oracle(base, honest, s_victim, mid, i, 1, stats=stats) == 0:
            hi = mid
        else:
            lo = mid + 1
    u = lo
    if u == 0:
        return 0
    if ES[i] == 1:
        return u
    lo, hi = 1, ES[i]
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if oracle(base, honest, s_victim, u - 1, i, mid, stats=stats) == 1:
            lo = mid
        else:
            hi = mid - 1
    return (u - 1) * ES[i] + lo

# ================================================================ main
def main():
    global B
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--seed", type=int, default=2026, help="seed for the base cycle and the keys")
    ap.add_argument("--B", type=int, default=5, help="number of iterations (s_i in [0, B*e_i])")
    ap.add_argument("--keys", type=int, default=5, help="number of random keys to attack")
    ap.add_argument("--swap-skips", action="store_true",
                    help="also scan faults that skip the first or the second conditional swap")
    args = ap.parse_args()
    if args.B < 3:
        ap.error("--B must be at least 3 (the demo needs a fault before the last iteration)")
    B = args.B
    rng = random.Random(args.seed)
    t0 = time.time()

    print("=" * 78)
    print("Toy OSIDH-LD (Algorithm 2) with real isogenies, and the mask fault attack")
    print("=" * 78)
    print(f"p = 4*M - 1 = {P} ({math.log2(P):.1f} bits), primes {ELLS}, torsion powers e_i = {ES}")
    print(f"B = {B} iterations, so s_i in [0, B*e_i]; key space "
          f"{math.prod(B * e + 1 for e in ES):,} = 2^{sum(math.log2(B * e + 1) for e in ES):.1f}")
    worst = sum(math.ceil(math.log2(B * e + 1)) for e in ES)
    print(f"information-theoretic worst case: sum ceil(log2(B*e_i + 1)) = {worst} oracle queries")

    # ---- public base cycle: E_1 = a0 * E_0, E_{j+1} = prod [l_i]^e_i * E_j, closes after r steps
    A = reference_action(0, [rng.randrange(-3, 4) for _ in range(N)])
    As = [A]
    while True:
        A = reference_action(A, ES)
        if A == As[0]:
            break
        As.append(A)
    r = len(As)
    base = [(Aj, point_of_order_M(Aj, +1, rng), point_of_order_M(Aj, -1, rng)) for Aj in As]
    print(f"base orientation cycle: r = {r} curves, A = {As}")

    # ---- 1. correctness of the deterministic action
    print("\n1. Correctness of Algorithm 2 (no fault)")
    while True:                         # demo key with an e_i = 1 component in [1, B-2]
        s = [rng.randrange(0, B * e + 1) for e in ES]
        if any(ES[k] == 1 and 1 <= s[k] <= B - 2 for k in range(N)):
            break
    honest_cycle = group_action(base, s)
    honest = curves(honest_cycle)
    ref = [reference_action(Aj, s) for Aj in As]
    print(f"   secret s = {s}")
    print(f"   output curves {honest}")
    print(f"   equal to prod [l_i]^s_i * E_j from the reference action: {honest == ref}")
    assert honest == ref

    # ---- 2. what one fault does (Sec. 4.3, 5.2, 6.1)
    print("\n2. Effect of forcing one mask to 0 (one e_i = 1 component, cycle position j = 1)")
    i = next(k for k in range(N) if ES[k] == 1 and 1 <= s[k] <= B - 2)
    j = 1 % r
    u_first = s[i]                      # e_i = 1: the mask at t = 1 is 1 from iteration s_i on
    for u, why in ((u_first - 1, "mask 0"),
                   (u_first, "mask 1, not the last iteration"),
                   (B - 1, "mask 1, last iteration")):
        try:
            out = group_action(base, s, fault=(u, j, i, 1, "mask"))
            changed = [k for k in range(r) if out[k][0] != honest[k]]
            msg = f"completes, curves changed at positions {changed}"
            if changed == [j]:
                shifted = reference_action(honest[j], [1 if k == i else 0 for k in range(N)])
                msg += f"; faulty E_j == [l_{ELLS[i]}] * honest E_j: {out[j][0] == shifted}"
        except LaneCheckFailure as e:
            msg = (f"lanes agree in iteration {u}; lane check fails in iteration {e.u} "
                   f"at positions {e.positions} (j-1, j = {(j - 1) % r}, {j})")
        print(f"   l = {ELLS[i]}, s_i = {s[i]}, fault at u = {u} [{why}]: {msg}")

    # ---- 3. oracle soundness over every fault position
    print("\n3. Every fault position (u, j, i, t) for this key")
    kinds = [("mask", "mask forced to 0 (the paper's model)")]
    if args.swap_skips:
        kinds += [("swap1", "first Cswap skipped"), ("swap2", "second Cswap skipped")]
    for kind, label in kinds:
        n_pos, sound, same_iter, wrong = 0, 0, 0, []
        for u in range(B):
            for jj in range(r):
                for ii in range(N):
                    for t in range(1, ES[ii] + 1):
                        n_pos += 1
                        mask = 1 if t > min(max(s[ii] - u * ES[ii], 0), ES[ii]) else 0
                        try:
                            out = group_action(base, s, fault=(u, jj, ii, t, kind))
                            o = 1 if curves(out) == honest else 0
                        except LaneCheckFailure as e:
                            o = 0
                            same_iter += (e.u == u)
                        if o == 1 - mask:
                            sound += 1
                        else:
                            wrong.append((u, jj, ELLS[ii], t))
        print(f"   {label}: {n_pos} positions; lanes disagree in the faulted iteration at "
              f"{same_iter}; output changed exactly when the mask was 1 at {sound}/{n_pos}")
        if wrong:
            print(f"      no visible change despite mask 1 at (u, j, l, t) = {wrong}: the last step of a "
                  f"position in the last iteration, where skipping the second swap only exchanges P' and Q'")

    # ---- 4. key recovery
    print(f"\n4. Key recovery with the two-stage binary search, {args.keys} random keys")
    total = 0
    for k in range(args.keys):
        sv = s if k == 0 else [rng.randrange(0, B * e + 1) for e in ES]
        hv = honest if k == 0 else curves(group_action(base, sv))
        stats = {"queries": 0}
        rec = [recover_component(base, hv, sv, ii, stats) for ii in range(N)]
        ok = rec == sv and curves(group_action(base, rec)) == hv
        total += stats["queries"]
        print(f"   key {k + 1}: s = {sv} -> recovered {rec}: {'correct' if ok else 'WRONG'}, "
              f"{stats['queries']} oracle queries (worst case {worst})")
        assert ok and stats["queries"] <= worst
    print(f"   average {total / args.keys:.1f} queries per key")
    print(f"\nTotal time {time.time() - t0:.1f} s")

if __name__ == "__main__":
    main()
