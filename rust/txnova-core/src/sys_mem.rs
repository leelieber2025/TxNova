//! Host probe and live BAM-worker scheduling.
//!
//! Pattern from scdm: `threads: 0` sets a CPU ceiling; live concurrency
//! follows `MemAvailable` and can rise after startup. Explicit `N` is a
//! ceiling, not a promise to ignore RAM. Does not share scdm extract slots.

use std::fs;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;
use std::time::Duration;

use crate::error::{CoreError, Result};

/// RAM held by one concurrent BAM sample worker (reader + PE pending).
pub const BAM_WORKER_BYTES: u64 = 2 * 1024 * 1024 * 1024;
/// Headroom for OS + Python + a neighbor job (e.g. scdm).
pub const RESERVE_BYTES: u64 = 4 * 1024 * 1024 * 1024;
/// Do not grow while free RAM is below this.
const CRITICAL_BYTES: u64 = 4 * 1024 * 1024 * 1024;
/// Extra BAM readers past this are disk-bound.
const BAM_WORKERS_CEILING: usize = 32;

const ACQUIRE_POLL: Duration = Duration::from_millis(50);

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WorkPlan {
    pub requested: usize,
    pub cpus: usize,
    pub usable: usize,
    pub available_bytes: Option<u64>,
    /// Concurrency *ceiling*. Live workers follow MemAvailable.
    pub bam_workers: usize,
    pub auto: bool,
}

pub fn available_memory_bytes() -> Option<u64> {
    #[cfg(target_os = "linux")]
    {
        parse_meminfo_available(&fs::read_to_string("/proc/meminfo").ok()?)
    }
    #[cfg(not(target_os = "linux"))]
    {
        None
    }
}

pub fn parse_meminfo_available(text: &str) -> Option<u64> {
    let mut mem_available_kb: Option<u64> = None;
    let mut mem_free_kb: Option<u64> = None;
    let mut buffers_kb: u64 = 0;
    let mut cached_kb: u64 = 0;

    for line in text.lines() {
        let mut parts = line.split_whitespace();
        let key = parts.next()?;
        let val: u64 = parts.next()?.parse().ok()?;
        match key {
            "MemAvailable:" => mem_available_kb = Some(val),
            "MemFree:" => mem_free_kb = Some(val),
            "Buffers:" => buffers_kb = val,
            "Cached:" => cached_kb = val,
            _ => {}
        }
    }

    let kb = mem_available_kb.or_else(|| {
        mem_free_kb.map(|free| free.saturating_add(buffers_kb).saturating_add(cached_kb))
    })?;
    Some(kb.saturating_mul(1024))
}

pub fn cpu_count() -> Option<usize> {
    thread::available_parallelism().ok().map(|n| n.get())
}

/// Leave one core for OS / the shell on anything that is not a tiny box.
pub fn usable_cpus(cpus: usize) -> usize {
    if cpus >= 4 {
        cpus - 1
    } else {
        cpus.max(1)
    }
}

/// Sample-parallel BAM workers from usable cores. Half the box: each worker
/// is one full BAM scan, not a CPU-bound inner loop.
pub fn auto_bam_from_cpu(usable: usize) -> usize {
    let n = match usable {
        0 | 1 => 1,
        2 | 3 => usable,
        _ => (usable / 2).max(2),
    };
    n.clamp(1, BAM_WORKERS_CEILING)
}

/// Ceiling only. `0` = CPU auto. Memory does not shrink the ceiling.
pub fn resolve_work_plan(
    requested: usize,
    cpu_limit: Option<usize>,
    available_bytes: Option<u64>,
    n_samples: usize,
) -> WorkPlan {
    let cpus = cpu_limit.filter(|n| *n > 0).unwrap_or(1);
    let usable = usable_cpus(cpus);
    let sample_cap = n_samples.max(1);
    let auto = requested == 0;
    let bam_workers = if auto {
        auto_bam_from_cpu(usable).min(sample_cap).max(1)
    } else {
        requested.max(1).min(sample_cap)
    };
    WorkPlan {
        requested,
        cpus,
        usable,
        available_bytes,
        bam_workers,
        auto,
    }
}

/// How many BAM workers are safe *right now*.
pub fn target_nproc_now(max_cap: usize, inflight: usize, available: Option<u64>) -> usize {
    let max_cap = max_cap.max(1);
    let Some(avail) = available else {
        return max_cap;
    };
    if avail < CRITICAL_BYTES {
        return inflight.max(1).min(max_cap);
    }
    let our_est = (inflight as u64).saturating_mul(BAM_WORKER_BYTES);
    let effective = avail.saturating_add(our_est);
    let usable = effective.saturating_sub(RESERVE_BYTES);
    let budget = (usable / BAM_WORKER_BYTES) as usize;
    budget.max(1).min(max_cap)
}

/// Runtime gate: at most `max_cap` BAM workers, shrinking under memory pressure.
pub struct AdaptiveLimiter {
    max_cap: usize,
    inflight: AtomicUsize,
}

/// RAII slot: release concurrency on drop.
pub struct AdaptiveSlot<'a> {
    lim: &'a AdaptiveLimiter,
}

impl Drop for AdaptiveSlot<'_> {
    fn drop(&mut self) {
        self.lim.release();
    }
}

impl AdaptiveLimiter {
    pub fn new(max_cap: usize) -> Self {
        Self {
            max_cap: max_cap.max(1),
            inflight: AtomicUsize::new(0),
        }
    }

    pub fn acquire(&self) -> AdaptiveSlot<'_> {
        loop {
            let avail = available_memory_bytes();
            let cur = self.inflight.load(Ordering::Acquire);
            let target = target_nproc_now(self.max_cap, cur, avail);
            if cur < target {
                if self
                    .inflight
                    .compare_exchange(cur, cur + 1, Ordering::AcqRel, Ordering::Acquire)
                    .is_ok()
                {
                    return AdaptiveSlot { lim: self };
                }
                continue;
            }
            thread::sleep(ACQUIRE_POLL);
        }
    }

    fn release(&self) {
        let prev = self.inflight.fetch_sub(1, Ordering::AcqRel);
        debug_assert!(prev > 0, "AdaptiveLimiter release without acquire");
    }
}

/// Ceiling for this call, then live limiter. Re-probes CPU each time.
pub fn bam_worker_cap(requested: usize, n_samples: usize) -> usize {
    resolve_work_plan(requested, cpu_count(), available_memory_bytes(), n_samples).bam_workers
}

/// Run `n` sample jobs under the live limiter. `work(i)` runs on sample index `i`.
pub fn run_bam_jobs<T, F>(requested: usize, n: usize, work: F) -> Result<Vec<T>>
where
    T: Send,
    F: Fn(usize) -> Result<T> + Sync,
{
    if n == 0 {
        return Ok(Vec::new());
    }
    let cap = bam_worker_cap(requested, n);
    let limiter = AdaptiveLimiter::new(cap);
    let mut out: Vec<Option<T>> = (0..n).map(|_| None).collect();
    std::thread::scope(|scope| -> Result<()> {
        let handles: Vec<_> = (0..n)
            .map(|i| {
                let work = &work;
                let limiter = &limiter;
                scope.spawn(move || {
                    let _slot = limiter.acquire();
                    work(i).map(|v| (i, v))
                })
            })
            .collect();
        for h in handles {
            let (i, v) = h
                .join()
                .unwrap_or_else(|_| Err(CoreError::fail("BAM worker panicked")))?;
            out[i] = Some(v);
        }
        Ok(())
    })?;
    Ok(out
        .into_iter()
        .map(|x| x.expect("worker slot filled"))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    const GIB: u64 = 1024 * 1024 * 1024;

    #[test]
    fn parse_memavailable() {
        let sample = "\
MemTotal:       32768000 kB
MemFree:         1000000 kB
MemAvailable:   20000000 kB
Buffers:          500000 kB
Cached:          4000000 kB
";
        assert_eq!(parse_meminfo_available(sample), Some(20_000_000 * 1024));
    }

    #[test]
    fn parse_fallback_without_memavailable() {
        let sample = "\
MemTotal:       32768000 kB
MemFree:         2000000 kB
Buffers:          500000 kB
Cached:          1500000 kB
";
        assert_eq!(parse_meminfo_available(sample), Some(4_000_000 * 1024));
    }

    #[test]
    fn usable_leaves_one_core_from_four() {
        assert_eq!(usable_cpus(1), 1);
        assert_eq!(usable_cpus(2), 2);
        assert_eq!(usable_cpus(3), 3);
        assert_eq!(usable_cpus(4), 3);
        assert_eq!(usable_cpus(8), 7);
        assert_eq!(usable_cpus(64), 63);
    }

    #[test]
    fn bam_workers_scale_with_machine() {
        assert_eq!(auto_bam_from_cpu(1), 1);
        assert_eq!(auto_bam_from_cpu(2), 2);
        assert_eq!(auto_bam_from_cpu(3), 3);
        assert_eq!(auto_bam_from_cpu(7), 3);
        assert_eq!(auto_bam_from_cpu(15), 7);
        assert_eq!(auto_bam_from_cpu(31), 15);
        assert_eq!(auto_bam_from_cpu(63), 31);
        assert_eq!(auto_bam_from_cpu(127), 32);
    }

    #[test]
    fn auto_ceiling_is_cpu_not_memory() {
        let tight = resolve_work_plan(0, Some(16), Some(6 * GIB), 8);
        assert!(tight.auto);
        assert_eq!(tight.usable, 15);
        assert_eq!(tight.bam_workers, 7);
        let fat = resolve_work_plan(0, Some(64), Some(256 * GIB), 24);
        assert_eq!(fat.bam_workers, 24);
        let io_cap = resolve_work_plan(0, Some(128), Some(512 * GIB), 48);
        assert_eq!(io_cap.bam_workers, 32);
    }

    #[test]
    fn explicit_n_is_ceiling() {
        let p = resolve_work_plan(6, Some(4), Some(4 * GIB), 8);
        assert!(!p.auto);
        assert_eq!(p.bam_workers, 6);
    }

    #[test]
    fn workers_never_exceed_sample_count() {
        let p = resolve_work_plan(0, Some(32), Some(128 * GIB), 2);
        assert_eq!(p.bam_workers, 2);
    }

    #[test]
    fn live_target_shrinks_when_ram_is_tight() {
        assert_eq!(target_nproc_now(8, 0, Some(6 * GIB)), 1);
        assert_eq!(target_nproc_now(8, 0, Some(20 * GIB)), 8);
        assert_eq!(target_nproc_now(8, 3, Some(3 * GIB)), 3);
        assert_eq!(target_nproc_now(8, 0, None), 8);
    }

    #[test]
    fn limiter_runs_all_jobs() {
        let out = run_bam_jobs(2, 5, |i| Ok(i * 10)).unwrap();
        assert_eq!(out, vec![0, 10, 20, 30, 40]);
    }
}
