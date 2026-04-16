#!/usr/bin/env python3
"""
Parse a periodic memory log like the one shown and emit:
  1) a CSV of memory samples
  2) an optional CSV aligning error-log timestamps to nearest memory sample
  3) a short summary on stdout

Usage:
  python parse_memlog.py --memlog mem.log --out mem_samples.csv
  python parse_memlog.py --memlog mem.log --error-log err.txt --out mem_samples.csv --events-out events_aligned.csv
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean, median
from typing import Optional, List, Tuple


TS_RE = re.compile(r"^───\s+(?P<ts>\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2})\s+")
MEM_RE = re.compile(
    r"^Mem:\s+(?P<total_gi>\d+(?:\.\d+)?)Gi\s+"
    r"(?P<used_gi>\d+(?:\.\d+)?)Gi\s+"
    r"(?P<free_gi>\d+(?:\.\d+)?)Gi\s+"
    r"(?P<shared>(?:\d+(?:\.\d+)?[KMGTPiB]*)|\d+(?:\.\d+)?)\s+"
    r"(?P<buff_cache_gi>\d+(?:\.\d+)?)Gi\s+"
    r"(?P<available_gi>\d+(?:\.\d+)?)Gi\s*$"
)
SWAP_RE = re.compile(
    r"^Swap:\s+(?P<total>(?:\d+(?:\.\d+)?[KMGTPiB]*)|\d+(?:\.\d+)?)\s+"
    r"(?P<used>(?:\d+(?:\.\d+)?[KMGTPiB]*)|\d+(?:\.\d+)?)\s+"
    r"(?P<free>(?:\d+(?:\.\d+)?[KMGTPiB]*)|\d+(?:\.\d+)?)\s*$"
)
PSI_RE = re.compile(
    r"^(?P<kind>some|full)\s+avg10=(?P<avg10>\d+\.\d+)\s+avg60=(?P<avg60>\d+\.\d+)\s+avg300=(?P<avg300>\d+\.\d+)\s+total=(?P<total>\d+)\s*$"
)
PID_HEADER_RE = re.compile(r"^\s*PID\s+RSS\s+COMMAND\s*$")
PID_LINE_RE = re.compile(r"^\s*(?P<pid>\d+)\s+(?P<rss_kb>\d+)\s+(?P<cmd>\S+)\s*$")
ERR_TS_RE = re.compile(r"^\[(?P<ts>\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2})\]")


@dataclass
class MemSample:
    ts: datetime
    total_gi: float
    used_gi: float
    free_gi: float
    shared: str
    buff_cache_gi: float
    available_gi: float
    mem_free_kb: Optional[int] = None
    mem_available_kb: Optional[int] = None
    anon_pages_kb: Optional[int] = None
    committed_as_kb: Optional[int] = None
    swap_total_kb: Optional[int] = None
    swap_free_kb: Optional[int] = None
    psi_some_avg10: Optional[float] = None
    psi_some_avg60: Optional[float] = None
    psi_some_avg300: Optional[float] = None
    psi_some_total: Optional[int] = None
    psi_full_avg10: Optional[float] = None
    psi_full_avg60: Optional[float] = None
    psi_full_avg300: Optional[float] = None
    psi_full_total: Optional[int] = None
    top_rss_kb: Optional[int] = None
    top_rss_pid: Optional[int] = None
    top_rss_cmd: Optional[str] = None

    @property
    def used_pct(self) -> float:
        return (self.used_gi / self.total_gi) * 100.0 if self.total_gi else 0.0

    @property
    def avail_pct(self) -> float:
        return (self.available_gi / self.total_gi) * 100.0 if self.total_gi else 0.0

    @property
    def pressure_flag(self) -> bool:
        # simple heuristic: any sustained memory pressure would show up here
        return bool((self.psi_some_avg10 or 0.0) > 0.0 or (self.psi_full_avg10 or 0.0) > 0.0)


def parse_int(s: str) -> int:
    return int(s.replace(",", "").strip())


def parse_ts(s: str) -> datetime:
    # Logs are UTC in your example.
    return datetime.strptime(s, "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)


def parse_memlog(path: Path) -> List[MemSample]:
    samples: List[MemSample] = []

    current: dict = {}
    in_pid_block = False
    pid_candidates: List[Tuple[int, int, str]] = []

    def flush_current():
        nonlocal current, in_pid_block, pid_candidates
        if "ts" not in current:
            current = {}
            pid_candidates = []
            in_pid_block = False
            return

        sample = MemSample(
            ts=current["ts"],
            total_gi=current.get("total_gi", 0.0),
            used_gi=current.get("used_gi", 0.0),
            free_gi=current.get("free_gi", 0.0),
            shared=current.get("shared", ""),
            buff_cache_gi=current.get("buff_cache_gi", 0.0),
            available_gi=current.get("available_gi", 0.0),
            mem_free_kb=current.get("mem_free_kb"),
            mem_available_kb=current.get("mem_available_kb"),
            anon_pages_kb=current.get("anon_pages_kb"),
            committed_as_kb=current.get("committed_as_kb"),
            swap_total_kb=current.get("swap_total_kb"),
            swap_free_kb=current.get("swap_free_kb"),
            psi_some_avg10=current.get("psi_some_avg10"),
            psi_some_avg60=current.get("psi_some_avg60"),
            psi_some_avg300=current.get("psi_some_avg300"),
            psi_some_total=current.get("psi_some_total"),
            psi_full_avg10=current.get("psi_full_avg10"),
            psi_full_avg60=current.get("psi_full_avg60"),
            psi_full_avg300=current.get("psi_full_avg300"),
            psi_full_total=current.get("psi_full_total"),
        )

        if pid_candidates:
            pid_candidates.sort(reverse=True)
            sample.top_rss_kb, sample.top_rss_pid, sample.top_rss_cmd = pid_candidates[0]

        samples.append(sample)
        current = {}
        pid_candidates = []
        in_pid_block = False

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")

            m = TS_RE.match(line)
            if m:
                flush_current()
                current["ts"] = parse_ts(m.group("ts"))
                continue

            if not current:
                continue

            m = MEM_RE.match(line)
            if m:
                current["total_gi"] = float(m.group("total_gi"))
                current["used_gi"] = float(m.group("used_gi"))
                current["free_gi"] = float(m.group("free_gi"))
                current["shared"] = m.group("shared")
                current["buff_cache_gi"] = float(m.group("buff_cache_gi"))
                current["available_gi"] = float(m.group("available_gi"))
                continue

            if line.startswith("MemTotal:"):
                parts = line.split()
                current["mem_total_kb"] = parse_int(parts[1])
                continue
            if line.startswith("MemFree:"):
                current["mem_free_kb"] = parse_int(line.split()[1])
                continue
            if line.startswith("MemAvailable:"):
                current["mem_available_kb"] = parse_int(line.split()[1])
                continue
            if line.startswith("AnonPages:"):
                current["anon_pages_kb"] = parse_int(line.split()[1])
                continue
            if line.startswith("Committed_AS:"):
                current["committed_as_kb"] = parse_int(line.split()[1])
                continue
            if line.startswith("SwapTotal:"):
                current["swap_total_kb"] = parse_int(line.split()[1])
                continue
            if line.startswith("SwapFree:"):
                current["swap_free_kb"] = parse_int(line.split()[1])
                continue

            m = PSI_RE.match(line)
            if m:
                kind = m.group("kind")
                current[f"psi_{kind}_avg10"] = float(m.group("avg10"))
                current[f"psi_{kind}_avg60"] = float(m.group("avg60"))
                current[f"psi_{kind}_avg300"] = float(m.group("avg300"))
                current[f"psi_{kind}_total"] = int(m.group("total"))
                continue

            if PID_HEADER_RE.match(line):
                in_pid_block = True
                continue

            if in_pid_block:
                if not line.strip():
                    in_pid_block = False
                    continue
                m = PID_LINE_RE.match(line)
                if m:
                    pid_candidates.append(
                        (int(m.group("rss_kb")), int(m.group("pid")), m.group("cmd"))
                    )
                continue

    flush_current()
    return samples


def parse_error_timestamps(path: Path) -> List[datetime]:
    out: List[datetime] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            m = ERR_TS_RE.match(raw)
            if m:
                out.append(parse_ts(m.group("ts")))
    return out


def nearest_sample(samples: List[MemSample], t: datetime) -> Optional[MemSample]:
    if not samples:
        return None
    # linear scan is fine for a few thousand samples; easy to read and robust
    best = min(samples, key=lambda s: abs((s.ts - t).total_seconds()))
    return best


def write_samples_csv(samples: List[MemSample], out_path: Path) -> None:
    fieldnames = list(asdict(samples[0]).keys()) + [
        "used_pct",
        "avail_pct",
        "pressure_flag",
    ] if samples else []

    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for s in samples:
            row = asdict(s)
            row["ts"] = s.ts.isoformat()
            row["used_pct"] = round(s.used_pct, 3)
            row["avail_pct"] = round(s.avail_pct, 3)
            row["pressure_flag"] = s.pressure_flag
            writer.writerow(row)


def write_events_csv(samples: List[MemSample], event_ts: List[datetime], out_path: Path) -> None:
    with out_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "event_ts",
            "sample_ts",
            "delta_seconds",
            "used_gi",
            "free_gi",
            "available_gi",
            "used_pct",
            "mem_available_kb",
            "committed_as_kb",
            "psi_some_avg10",
            "psi_full_avg10",
            "top_rss_kb",
            "top_rss_pid",
            "top_rss_cmd",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for ts in event_ts:
            s = nearest_sample(samples, ts)
            if s is None:
                continue
            writer.writerow({
                "event_ts": ts.isoformat(),
                "sample_ts": s.ts.isoformat(),
                "delta_seconds": round(abs((s.ts - ts).total_seconds()), 3),
                "used_gi": s.used_gi,
                "free_gi": s.free_gi,
                "available_gi": s.available_gi,
                "used_pct": round(s.used_pct, 3),
                "mem_available_kb": s.mem_available_kb,
                "committed_as_kb": s.committed_as_kb,
                "psi_some_avg10": s.psi_some_avg10,
                "psi_full_avg10": s.psi_full_avg10,
                "top_rss_kb": s.top_rss_kb,
                "top_rss_pid": s.top_rss_pid,
                "top_rss_cmd": s.top_rss_cmd,
            })


def print_summary(samples: List[MemSample]) -> None:
    if not samples:
        print("No samples parsed.")
        return

    used = [s.used_gi for s in samples]
    avail = [s.available_gi for s in samples]
    committed = [s.committed_as_kb for s in samples if s.committed_as_kb is not None]
    psi_some = [s.psi_some_avg10 for s in samples if s.psi_some_avg10 is not None]
    psi_full = [s.psi_full_avg10 for s in samples if s.psi_full_avg10 is not None]

    print(f"Samples: {len(samples)}")
    print(f"Time range: {samples[0].ts.isoformat()} -> {samples[-1].ts.isoformat()}")
    print(f"Used GiB:   min={min(used):.2f}  median={median(used):.2f}  max={max(used):.2f}  mean={mean(used):.2f}")
    print(f"Avail GiB:  min={min(avail):.2f}  median={median(avail):.2f}  max={max(avail):.2f}  mean={mean(avail):.2f}")
    if committed:
        print(f"Committed_AS kB: min={min(committed)}  median={int(median(committed))}  max={max(committed)}")
    if psi_some:
        print(f"PSI some avg10: max={max(psi_some):.2f}")
    if psi_full:
        print(f"PSI full avg10: max={max(psi_full):.2f}")

    peak = max(samples, key=lambda s: s.used_gi)
    print(
        "Peak sample: "
        f"{peak.ts.isoformat()} used={peak.used_gi:.2f}GiB "
        f"avail={peak.available_gi:.2f}GiB "
        f"committed={peak.committed_as_kb if peak.committed_as_kb is not None else 'NA'}kB "
        f"psi_some={peak.psi_some_avg10 if peak.psi_some_avg10 is not None else 'NA'} "
        f"psi_full={peak.psi_full_avg10 if peak.psi_full_avg10 is not None else 'NA'}"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--memlog", required=True, type=Path)
    ap.add_argument("--error-log", type=Path)
    ap.add_argument("--out", required=True, type=Path, help="CSV path for memory samples")
    ap.add_argument("--events-out", type=Path, help="CSV path for error timestamps aligned to nearest memory sample")
    args = ap.parse_args()

    samples = parse_memlog(args.memlog)
    if not samples:
        print("No memory samples found.")
        return 1

    write_samples_csv(samples, args.out)
    print_summary(samples)

    if args.error_log and args.events_out:
        event_ts = parse_error_timestamps(args.error_log)
        write_events_csv(samples, event_ts, args.events_out)
        print(f"Wrote aligned events CSV: {args.events_out}")

    print(f"Wrote memory samples CSV: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())