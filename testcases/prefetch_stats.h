#ifndef PREFETCH_STATS_H
#define PREFETCH_STATS_H

#include <cstdio>
#include <sys/time.h>

#ifdef __APPLE__
#include <mach/mach.h>
#include <mach/task_info.h>
#endif

namespace prefetch_stats {

struct PrefetchStats {
  unsigned long prefetch_count = 0;
  unsigned long cache_misses_before = 0;
  unsigned long cache_misses_after = 0;
  double prefetch_overhead_us = 0.0;
};

static PrefetchStats global_stats;

inline double get_time_us() {
  struct timeval tv;
  gettimeofday(&tv, nullptr);
  return tv.tv_sec * 1e6 + tv.tv_usec;
}

inline void record_prefetch() { global_stats.prefetch_count++; }

inline void start_prefetch_timing() {
  global_stats.prefetch_overhead_us -= get_time_us();
}

inline void end_prefetch_timing() {
  global_stats.prefetch_overhead_us += get_time_us();
}

#ifdef __APPLE__
inline void get_memory_stats(unsigned long &resident_size,
                             unsigned long &page_faults) {
  struct mach_task_basic_info info;
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info,
                &count) == KERN_SUCCESS) {
    resident_size = info.resident_size;
  } else {
    resident_size = 0;
  }
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  page_faults = r.ru_majflt + r.ru_minflt;
}
#endif

inline void print_stats() {
  printf("\nPrefetch Statistics\n");
  printf("Prefetch calls:    %lu\n", global_stats.prefetch_count);
  printf("Prefetch overhead: %.2f ms\n",
         global_stats.prefetch_overhead_us / 1000.0);

#ifdef __APPLE__
  unsigned long resident, faults;
  get_memory_stats(resident, faults);
  printf("Resident memory:   %.2f MB\n", resident / (1024.0 * 1024.0));
  printf("Page faults:       %lu\n", faults);
#endif
}

struct ScopedTiming {
  ScopedTiming() { start_prefetch_timing(); }
  ~ScopedTiming() { end_prefetch_timing(); }
};

} // namespace prefetch_stats

#ifdef PROFILE_PREFETCH
#define PREFETCH_RECORD() prefetch_stats::record_prefetch()
#define PREFETCH_TIMING_SCOPE() prefetch_stats::ScopedTiming _pf_timing
#define PREFETCH_PRINT_STATS() prefetch_stats::print_stats()
#else
#define PREFETCH_RECORD() ((void)0)
#define PREFETCH_TIMING_SCOPE() ((void)0)
#define PREFETCH_PRINT_STATS() ((void)0)
#endif

#endif
