// SPDX-FileCopyrightText: Copyright 2015 Google Inc. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

/*
 * This file is derived from cycleclock.h in the
 * benchmark project (https://github.com/google/benchmark),
 * which is licensed under the Apache-2.0 license.
 *
 * Modifications made:
 *   - Modified some macro definitions to make them self-contained
 */

// ----------------------------------------------------------------------
// CycleClock
//    A CycleClock tells you the current time in Cycles.  The "time"
//    is actually time since power-on.  This is like time() but doesn't
//    involve a system call and is much more precise.
//
// NOTE: Not all cpu/platform/kernel combinations guarantee that this
// clock increments at a constant rate or is synchronized across all logical
// cpus in a system.
//
// If you need the above guarantees, please consider using a different
// API. There are efforts to provide an interface which provides a millisecond
// granularity and implemented as a memory read. A memory read is generally
// cheaper than the CycleClock for many architectures.
//
// Also, in some out of order CPU implementations, the CycleClock is not
// serializing. So if you're trying to count at cycles granularity, your
// data might be inaccurate due to out of order instruction execution.
// ----------------------------------------------------------------------

#ifndef BENCHMARK_CYCLECLOCK_H_
#define BENCHMARK_CYCLECLOCK_H_

#include <stdint.h>

#if defined(__GNUC__) || defined(__clang__)
#define BENCHMARK_ALWAYS_INLINE __attribute__((always_inline))
#elif defined(_MSC_VER) && !defined(__clang__)
#define BENCHMARK_ALWAYS_INLINE __forceinline
#else
#define BENCHMARK_ALWAYS_INLINE
#endif

#if defined(__APPLE__)
  #include "TargetConditionals.h"
#endif

#if defined(TARGET_OS_MAC)
#include <mach/mach_time.h>
#endif
// For MSVC, we want to use '_asm rdtsc' when possible (since it works
// with even ancient MSVC compilers), and when not possible the
// __rdtsc intrinsic, declared in <intrin.h>.  Unfortunately, in some
// environments, <windows.h> and <intrin.h> have conflicting
// declarations of some other intrinsics, breaking compilation.
// Therefore, we simply declare __rdtsc ourselves. See also
// http://connect.microsoft.com/VisualStudio/feedback/details/262047
//
// Note that MSVC defines the x64 preprocessor macros when building
// for Arm64EC, despite it using Arm64 assembly instructions.
#if defined(_MSC_VER) && !defined(_M_IX86) && !defined(_M_ARM64) && \
    !defined(_M_ARM64EC)
extern "C" uint64_t __rdtsc();
#pragma intrinsic(__rdtsc)
#endif

#if !defined(_WIN32) || defined(__MINGW32__)
#include <sys/time.h>
#include <time.h>
#endif

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

inline BENCHMARK_ALWAYS_INLINE int64_t Now() {
#if defined(TARGET_OS_MAC)
  // this goes at the top because we need ALL Macs, regardless of
  // architecture, to return the number of "mach time units" that
  // have passed since startup.  See sysinfo.cc where
  // InitializeSystemInfo() sets the supposed cpu clock frequency of
  // macs to the number of mach time units per second, not actual
  // CPU clock frequency (which can change in the face of CPU
  // frequency scaling).  Also note that when the Mac sleeps, this
  // counter pauses; it does not continue counting, nor does it
  // reset to zero.
  return (int64_t)(mach_absolute_time());
#elif defined(__EMSCRIPTEN__)
  // this goes above x86-specific code because old versions of Emscripten
  // define __x86_64__, although they have nothing to do with it.
  return (int64_t)(emscripten_get_now() * 1e+6);
#elif defined(__i386__)
  int64_t ret;
  __asm__ volatile("rdtsc" : "=A"(ret));
  return ret;

// Note that Clang, like MSVC, defines the x64 preprocessor macros when building
// for Arm64EC, despite it using Arm64 assembly instructions.
#elif (defined(__x86_64__) || defined(__amd64__)) && !defined(__arm64ec__)
  uint64_t low, high;
  __asm__ volatile("rdtsc" : "=a"(low), "=d"(high));
  return (int64_t)((high << 32) | low);
#elif defined(__powerpc__) || defined(__ppc__)
  // This returns a time-base, which is not always precisely a cycle-count.
#if defined(__powerpc64__) || defined(__ppc64__)
  int64_t tb;
  asm volatile("mfspr %0, 268" : "=r"(tb));
  return tb;
#else
  uint32_t tbl, tbu0, tbu1;
  asm volatile(
      "mftbu %0\n"
      "mftb %1\n"
      "mftbu %2"
      : "=r"(tbu0), "=r"(tbl), "=r"(tbu1));
  tbl &= -(int32_t)(tbu0 == tbu1);
  // high 32 bits in tbu1; low 32 bits in tbl  (tbu0 is no longer needed)
  return ((uint64_t)(tbu1) << 32) | tbl;
#endif
#elif defined(__sparc__)
  int64_t tick;
  asm(".byte 0x83, 0x41, 0x00, 0x00");
  asm("mov   %%g1, %0" : "=r"(tick));
  return tick;
#elif defined(__ia64__)
  int64_t itc;
  asm("mov %0 = ar.itc" : "=r"(itc));
  return itc;
#elif defined(_MSC_VER) && defined(_M_IX86)
  // Older MSVC compilers (like 7.x) don't seem to support the
  // __rdtsc intrinsic properly, so I prefer to use _asm instead
  // when I know it will work.  Otherwise, I'll use __rdtsc and hope
  // the code is being compiled with a non-ancient compiler.
  _asm rdtsc
#elif defined(_MSC_VER) && (defined(_M_ARM64) || defined(_M_ARM64EC))
  // See // https://docs.microsoft.com/en-us/cpp/intrinsics/arm64-intrinsics
  // and https://reviews.llvm.org/D53115
  int64_t virtual_timer_value;
  virtual_timer_value = _ReadStatusReg(ARM64_CNTVCT);
  return virtual_timer_value;
#elif defined(_MSC_VER)
  return __rdtsc();
#elif defined(__native_client__)
  // Native Client validator on x86/x86-64 allows RDTSC instructions,
  // and this case is handled above. Native Client validator on ARM
  // rejects MRC instructions (used in the ARM-specific sequence below),
  // so we handle it here. Portable Native Client compiles to
  // architecture-agnostic bytecode, which doesn't provide any
  // cycle counter access mnemonics.

  // Native Client does not provide any API to access cycle counter.
  // Use clock_gettime(CLOCK_MONOTONIC, ...) instead of gettimeofday
  // because is provides nanosecond resolution (which is noticeable at
  // least for PNaCl modules running on x86 Mac & Linux).
  // Initialize to always return 0 if clock_gettime fails.
  struct timespec ts = {0, 0};
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (int64_t)(ts.tv_sec) * 1000000000 + ts.tv_nsec;
#elif defined(__aarch64__) || defined(__arm64ec__)
  // System timer of ARMv8 runs at a different frequency than the CPU's.
  // The frequency is fixed, typically in the range 1-50MHz.  It can be
  // read at CNTFRQ special register.  We assume the OS has set up
  // the virtual timer properly.
  int64_t virtual_timer_value;
  asm volatile("mrs %0, cntvct_el0" : "=r"(virtual_timer_value));
  return virtual_timer_value;
#elif defined(__ARM_ARCH)
  // V6 is the earliest arch that has a standard cyclecount
  // Native Client validator doesn't allow MRC instructions.
#if (__ARM_ARCH >= 6)
  uint32_t pmccntr;
  uint32_t pmuseren;
  uint32_t pmcntenset;
  // Read the user mode perf monitor counter access permissions.
  asm volatile("mrc p15, 0, %0, c9, c14, 0" : "=r"(pmuseren));
  if (pmuseren & 1) {  // Allows reading perfmon counters for user mode code.
    asm volatile("mrc p15, 0, %0, c9, c12, 1" : "=r"(pmcntenset));
    if (pmcntenset & 0x80000000ul) {  // Is it counting?
      asm volatile("mrc p15, 0, %0, c9, c13, 0" : "=r"(pmccntr));
      // The counter is set up to count every 64th cycle
      return (int64_t)(pmccntr) * 64;  // Should optimize to << 6
    }
  }
#endif
  struct timeval tv;
  gettimeofday(&tv, nullptr);
  return (int64_t)(tv.tv_sec) * 1000000 + tv.tv_usec;
#elif defined(__mips__) || defined(__m68k__)
  // mips apparently only allows rdtsc for superusers, so we fall
  // back to gettimeofday.  It's possible clock_gettime would be better.
  struct timeval tv;
  gettimeofday(&tv, nullptr);
  return (int64_t)(tv.tv_sec) * 1000000 + tv.tv_usec;
#elif defined(__loongarch__) || defined(__csky__)
  struct timeval tv;
  gettimeofday(&tv, nullptr);
  return (int64_t)(tv.tv_sec) * 1000000 + tv.tv_usec;
#elif defined(__s390__)  // Covers both s390 and s390x.
  // Return the CPU clock.
  uint64_t tsc;
#if defined(BENCHMARK_OS_ZOS)
  // z/OS HLASM syntax.
  asm(" stck %0" : "=m"(tsc) : : "cc");
#else
  // Linux on Z syntax.
  asm("stck %0" : "=Q"(tsc) : : "cc");
#endif
  return tsc;
#elif defined(__riscv)  // RISC-V
  // Use RDTIME (and RDTIMEH on riscv32).
  // RDCYCLE is a privileged instruction since Linux 6.6.
#if __riscv_xlen == 32
  uint32_t cycles_lo, cycles_hi0, cycles_hi1;
  // This asm also includes the PowerPC overflow handling strategy, as above.
  // Implemented in assembly because Clang insisted on branching.
  asm volatile(
      "rdtimeh %0\n"
      "rdtime %1\n"
      "rdtimeh %2\n"
      "sub %0, %0, %2\n"
      "seqz %0, %0\n"
      "sub %0, zero, %0\n"
      "and %1, %1, %0\n"
      : "=r"(cycles_hi0), "=r"(cycles_lo), "=r"(cycles_hi1));
  return (int64_t)(((uint64_t)(cycles_hi1) << 32) |
                              cycles_lo);
#else
  uint64_t cycles;
  asm volatile("rdtime %0" : "=r"(cycles));
  return (int64_t)(cycles);
#endif
#elif defined(__e2k__) || defined(__elbrus__)
  struct timeval tv;
  gettimeofday(&tv, nullptr);
  return (int64_t)(tv.tv_sec) * 1000000 + tv.tv_usec;
#elif defined(__hexagon__)
  uint64_t pcycle;
  asm volatile("%0 = C15:14" : "=r"(pcycle));
  return (int64_t)(pcycle);
#elif defined(__alpha__)
  // Alpha has a cycle counter, the PCC register, but it is an unsigned 32-bit
  // integer and thus wraps every ~4s, making using it for tick counts
  // unreliable beyond this time range.  The real-time clock is low-precision,
  // roughtly ~1ms, but it is the only option that can reasonable count
  // indefinitely.
  struct timeval tv;
  gettimeofday(&tv, nullptr);
  return (int64_t)(tv.tv_sec) * 1000000 + tv.tv_usec;
#elif defined(__hppa__) || defined(__linux__)
  // Fallback for all other architectures with a recent Linux kernel, e.g.:
  // HP PA-RISC provides a user-readable clock counter (cr16), but
  // it's not syncronized across CPUs and only 32-bit wide when programs
  // are built as 32-bit binaries.
  // Same for SH-4 and possibly others.
  // Use clock_gettime(CLOCK_MONOTONIC, ...) instead of gettimeofday
  // because is provides nanosecond resolution.
  // Initialize to always return 0 if clock_gettime fails.
  struct timespec ts = {0, 0};
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (int64_t)(ts.tv_sec) * 1000000000 + ts.tv_nsec;
#else
  // The soft failover to a generic implementation is automatic only for ARM.
  // For other platforms the developer is expected to make an attempt to create
  // a fast implementation and use generic version if nothing better is
  // available.
#error You need to define CycleTimer for your OS and CPU
#endif
}

#endif