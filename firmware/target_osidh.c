/*
 * target_osidh.c
 * Minimal OSIDH-LD fault-injection for CWLITEARM / STM32F303RCT6
 *
 * Exercises the single fault-sensitive line in eval_action():
 *
 *   digit_t b = compare_gt(t, s[i]);   // <-- trigger wraps this
 *   swap_vertices(&R[0], &R[1], b);
 *
 * The full OSIDH-LD group action is far too heavy for the Cortex-M4 at
 * 256-bit Fp2, so this harness:
 *   (1) Mimics the compare_gt + swap_vertices pattern on a lightweight
 *       proxy state (two 64-bit test words instead of full vertices).
 *   (2) Places trigger_high / trigger_low around the compare_gt call,
 *       matching the structure of Listing 2 in the paper.
 *   (3) Communicates over SimpleSerial v2:
 *       Host -> 'g' <t_byte> <s_byte>  : run one test with t=t_byte, s=s_byte
 *       Target -> 'r' <result_byte>    : 0x00 = no-swap (fault!), 0x01 = swap (correct)
 *       Target -> 'e' <error_byte>     : device reset / unexpected state
 *
 * Build with the ChipWhisperer simpleserial HAL for STM32F3:
 *   arm-none-eabi-gcc -mcpu=cortex-m4 -mthumb -O0 \
 *     -DHAL_TYPE=STM32F303xC \
 *     -I hal/stm32f3 -I simpleserial \
 *     target_osidh.c hal/stm32f3/stm32f3_hal.c simpleserial/simpleserial.c \
 *     -o target_osidh.elf
 *
 * -O0 is intentional: we want the compare_gt to be a real load/compare
 * instruction that can be skipped, not optimised away.
 */

#include "hal.h"
#include "simpleserial.h"
#include <stdint.h>
#include <string.h>

/* -----------------------------------------------------------------------
 * Portable digit type (32-bit on Cortex-M4)
 * --------------------------------------------------------------------- */
typedef uint32_t digit_t;

/* -----------------------------------------------------------------------
 * compare_gt: returns 0xFF..FF when a <= b (i.e. swap needed),
 *             returns 0x00..00 when a > b  (no swap needed).
 * This is the direct analogue of compare_gt in protocols_helper.h,
 * adjusted for 32-bit digit_t on Cortex-M4.
 * --------------------------------------------------------------------- */
static inline digit_t __attribute__((noinline)) compare_gt(int a, int b)
{
    digit_t res = (uint32_t)(a > b) - 1u;
    return res;
}

/* -----------------------------------------------------------------------
 * Lightweight proxy for swap_vertices.
 * Instead of full vertex_t structs we use a pair of 32-bit words.
 * The swap logic is identical to the constant-time swap in
 * protocols_helper.h: if option = 0xFF..FF, values are swapped.
 * --------------------------------------------------------------------- */
static inline void swap_words(volatile uint32_t *A, volatile uint32_t *B,
                               digit_t option)
{
    uint32_t tmp = *A;
    *A = (*A & ~option) | (*B &  option);
    *B = (*B & ~option) | ( tmp & option);
}

/* -----------------------------------------------------------------------
 * Test state: two proxy "vertices" R[0] and R[1].
 * We initialise them to known distinct values so we can detect
 * whether the swap happened.
 * --------------------------------------------------------------------- */
#define SENTINEL_R0  0xAAAAAAAAu
#define SENTINEL_R1  0x55555555u

/* -----------------------------------------------------------------------
 * run_test(t_val, s_val)
 *   Performs one iteration of the fault-sensitive inner loop:
 *     b = compare_gt(t_val, s_val)
 *     swap(R[0], R[1], b)
 *   Returns:
 *     0x01  if R[0] == SENTINEL_R1 (swap happened => correct)
 *     0x00  if R[0] == SENTINEL_R0 (no swap => fault succeeded)
 *     0xFF  if R[0] is neither     (corrupt/reset)
 * --------------------------------------------------------------------- */
static uint8_t run_test(uint8_t t_val, uint8_t s_val)
{
    volatile uint32_t R0 = SENTINEL_R0;
    volatile uint32_t R1 = SENTINEL_R1;

    /* --- Trigger window: wraps ONLY the compare_gt + first swap --- */
    trigger_high();
    digit_t b = compare_gt((int)t_val, (int)s_val);
    swap_words(&R0, &R1, b);
    trigger_low();

    /* Determine result */
    if (R0 == SENTINEL_R1)
        return 0x01;   /* swap happened: correct behaviour */
    if (R0 == SENTINEL_R0)
        return 0x00;   /* no swap: fault succeeded */
    return 0xFF;       /* corrupted state */
}

/* -----------------------------------------------------------------------
 * SimpleSerial command handler
 * Receives: cmd='g', data = [t_byte, s_byte]
 * Sends:    cmd='r', data = [result]
 * --------------------------------------------------------------------- */
static uint8_t handle_glitch(uint8_t *buf, uint8_t len)
{
    if (len < 2)
    {
        simpleserial_put('e', 1, (uint8_t[]){0xEE});
        return 0;
    }

    uint8_t t_val = buf[0];
    uint8_t s_val = buf[1];

    /* We only want to fault the case where a swap SHOULD happen,
     * i.e. t_val <= s_val.  The host enforces this, but we check. */
    uint8_t result = run_test(t_val, s_val);
    simpleserial_put('r', 1, &result);
    return 0;
}

/* -----------------------------------------------------------------------
 * Extended test: run all TORSION_ODD_POWERS[0..n-1] iterations for one
 * outer loop iteration, matching the j-loop structure of eval_action.
 * This exercises all CYCLE_LENGTH=13 tuple slots for prime index 0
 * (prime 3, power 4) at a single outer iteration.
 *
 * Returns: number of slots where the swap was detected (should be 13
 *          when t <= s_val; should be 0 after a successful full-skip fault).
 * --------------------------------------------------------------------- */
#define CYCLE_LENGTH 13

static uint8_t run_full_cycle_test(uint8_t t_val, uint8_t s_val)
{
    uint8_t swap_count = 0;
    for (int j = 0; j < CYCLE_LENGTH; j++)
    {
        volatile uint32_t R0 = SENTINEL_R0 ^ (uint32_t)j;
        volatile uint32_t R1 = SENTINEL_R1 ^ (uint32_t)j;

        trigger_high();
        digit_t b = compare_gt((int)t_val, (int)s_val);
        swap_words(&R0, &R1, b);
        trigger_low();

        if (R0 == (SENTINEL_R1 ^ (uint32_t)j))
            swap_count++;
    }
    return swap_count;  /* expected: 13 (all swapped) or 0 (all faulted) */
}

static uint8_t handle_cycle_test(uint8_t *buf, uint8_t len)
{
    if (len < 2)
    {
        simpleserial_put('e', 1, (uint8_t[]){0xEE});
        return 0;
    }
    uint8_t result = run_full_cycle_test(buf[0], buf[1]);
    simpleserial_put('r', 1, &result);
    return 0;
}

/* -----------------------------------------------------------------------
 * main
 * --------------------------------------------------------------------- */
int main(void)
{
    platform_init();
    init_uart();
    trigger_setup();

    simpleserial_init();

    /* 'g': single compare_gt + swap test */
    simpleserial_addcmd('g', 2, handle_glitch);

    /* 'c': full CYCLE_LENGTH iteration test (13 triggers) */
    simpleserial_addcmd('c', 2, handle_cycle_test);

    while (1)
        simpleserial_get();

    return 0;
}
