#!/usr/bin/env python3
"""Sequential launcher for seqtoid-pipelines short_read_mngs.

Runs each sample twice:
  1) with --use-diamond
  2) without --use-diamond

Logs are kept per sample/per mode:
  आउट dir/
    sample_label/
      diamond/
        out.txt
        err.txt
        command.sh
        run.json
      mmseqs/
        out.txt
        err.txt
        command.sh
        run.json

Sample input can come from a CSV/TSV sheet with columns:
  label,r1,r2
or from repeated --sample specs like:
  --sample label=/path/to/R1.fastq.gz,/path/to/R2.fastq.gz
  --sample label=/path/to/single.fastq.gz

The script is intentionally explicit: no directory scanning, no guessing.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class Sample:
    label: str
    r1: Path
    r2: Optional[Path] = None


def _clean_label(label: str) -> str:
    label = label.strip()
    if not label:
        raise ValueError("sample label cannot be empty")
    # Keep it filesystem-friendly.
    safe = []
    for ch in label:
        if ch.isalnum() or ch in ("-", "_", "."):
            safe.append(ch)
        else:
            safe.append("_")
    return "".join(safe)


def _parse_sample_spec(spec: str) -> Sample:
    """Parse LABEL=R1[,R2] or LABEL:R1[:R2]."""
    if "=" in spec:
        label, rest = spec.split("=", 1)
        parts = [p.strip() for p in rest.split(",") if p.strip()]
    elif ":" in spec:
        label, rest = spec.split(":", 1)
        parts = [p.strip() for p in rest.split(":") if p.strip()]
    else:
        raise ValueError(
            f"invalid --sample spec {spec!r}; use LABEL=R1[,R2] or LABEL:R1[:R2]"
        )

    label = _clean_label(label)
    if len(parts) == 1:
        return Sample(label=label, r1=Path(parts[0]).expanduser().resolve(), r2=None)
    if len(parts) == 2:
        return Sample(
            label=label,
            r1=Path(parts[0]).expanduser().resolve(),
            r2=Path(parts[1]).expanduser().resolve(),
        )
    raise ValueError(f"invalid --sample spec {spec!r}; expected 1 or 2 paths")


def _read_sample_sheet(path: Path) -> list[Sample]:
    if not path.exists():
        raise FileNotFoundError(f"sample sheet not found: {path}")

    text = path.read_text().splitlines()
    if not text:
        raise ValueError(f"sample sheet is empty: {path}")

    # Detect delimiter from the header line.
    header = next((line for line in text if line.strip() and not line.lstrip().startswith("#")), None)
    if header is None:
        raise ValueError(f"sample sheet has no usable rows: {path}")
    delimiter = "\t" if "\t" in header else ","

    samples: list[Sample] = []
    reader = csv.DictReader(
        (line for line in text if line.strip() and not line.lstrip().startswith("#")),
        delimiter=delimiter,
    )
    for row in reader:
        label = _clean_label((row.get("label") or row.get("sample") or "").strip())
        r1_raw = (row.get("r1") or row.get("R1") or "").strip()
        r2_raw = (row.get("r2") or row.get("R2") or "").strip()
        if not label or not r1_raw:
            raise ValueError(
                f"sample sheet row missing required label/r1 values: {row!r}"
            )
        samples.append(
            Sample(
                label=label,
                r1=Path(r1_raw).expanduser().resolve(),
                r2=Path(r2_raw).expanduser().resolve() if r2_raw else None,
            )
        )
    if not samples:
        raise ValueError(f"no samples parsed from {path}")
    return samples


def _add_optional_arg(cmd: list[str], flag: str, value: Optional[str]) -> None:
    if value is not None and value != "":
        cmd.extend([flag, value])


def build_command(
        runner: str,
        module: str,
        sample: Sample,
        *,
        threads: int,
        use_smt: bool,
        nvme_scratch: Path,
        use_diamond: bool,
        common_args: dict[str, Optional[str]],
        extra_args: list[str],
        out_dir: Path,
) -> list[str]:
    cmd: list[str] = [runner, "--module", module, "--threads", str(threads)]
    if use_smt:
        cmd.append("--use-smt")

    cmd.extend(["--nvme-scratch", str(nvme_scratch)])
    cmd.extend(["-o", str(out_dir)])

    if use_diamond:
        cmd.append("--use-diamond")

    cmd.extend(["-i", str(sample.r1)])
    if sample.r2 is not None:
        cmd.extend(["-I", str(sample.r2)])

    # Environment-specific pipeline arguments.
    for flag in (
            "kallisto_index",
            "ercc_bowtie2_index",
            "host_bowtie2_index",
            "human_bowtie2_index",
            "host_index",
            "taxid_lineages_db",
            "acc2taxid_db",
            "nt_db_size",
            "diamond_db",
            "nt",
            "nr",
            "nt_offset_db",
            "nr_offset_db",
            "nt_info_tab",
            "host_hisat2_index",
            "human_hisat2_index",
            "nt_split_dir",
            "mmseqs_db",
            "human_host",
            "seed",
            "max_subsample",
            "min_read_len",
            "max_read_len",
            "max_reads",
            "stall_threshold",
    ):
        val = common_args.get(flag)
        if val is None:
            continue
        cli_flag = "--" + flag.replace("_", "-")
        if isinstance(val, bool):
            if val:
                cmd.append(cli_flag)
        else:
            cmd.extend([cli_flag, str(val)])

    cmd.extend(extra_args)
    return cmd


def write_manifest(path: Path, *, cmd: list[str], rc: int, sample: Sample, use_diamond: bool) -> None:
    payload = {
        "sample": sample.label,
        "r1": str(sample.r1),
        "r2": str(sample.r2) if sample.r2 is not None else None,
        "use_diamond": use_diamond,
        "return_code": rc,
        "command": cmd,
    }
    path.write_text(json.dumps(payload, indent=2) + "\n")


def mode_tag(mode: str) -> str:
    return "diamond" if mode == "diamond" else "mmseqs2"


def rename_pipeline_output_dir(src_name: str, mode: str) -> str:
    tag = mode_tag(mode)
    if f"_{tag}_" in src_name or src_name.startswith(f"{tag}_") or src_name.endswith(f"_{tag}"):
        return src_name

    # Prefer to splice the backend tag immediately before the read marker.
    for needle in ("_R1.fastq", "_R2.fastq", "_R1", "_R2"):
        if needle in src_name:
            return src_name.replace(needle, f"_{tag}{needle}", 1)

    return f"{src_name}_{tag}"


def finalize_pipeline_output(pipeline_stage_dir: Path, run_dir: Path, mode: str) -> None:
    """Move the pipeline-generated output tree into a mode-tagged final location."""
    if not pipeline_stage_dir.exists():
        return

    # If the stage dir contains exactly one top-level directory, that is usually
    # the actual pipeline output tree. Move that tree instead of the stage wrapper.
    children = [p for p in pipeline_stage_dir.iterdir() if not p.name.startswith(".")]
    candidate = None
    if len(children) == 1 and children[0].is_dir():
        candidate = children[0]
    elif len(children) == 0:
        return
    else:
        candidate = pipeline_stage_dir

    if candidate is pipeline_stage_dir:
        dest_name = rename_pipeline_output_dir(pipeline_stage_dir.name, mode)
        dest = run_dir / dest_name
        if dest.exists():
            shutil.rmtree(dest)
        shutil.move(str(pipeline_stage_dir), str(dest))
        return

    dest_name = rename_pipeline_output_dir(candidate.name, mode)
    dest = run_dir / dest_name
    if dest.exists():
        shutil.rmtree(dest)
    shutil.move(str(candidate), str(dest))

    # Clean up the now-empty staging directory if possible.
    try:
        pipeline_stage_dir.rmdir()
    except OSError:
        pass


def run_one(
        *,
        runner: str,
        module: str,
        sample: Sample,
        out_root: Path,
        threads: int,
        use_smt: bool,
        nvme_scratch: Path,
        use_diamond: bool,
        common_args: dict[str, Optional[str]],
        extra_args: list[str],
        dry_run: bool,
) -> int:
    mode = "diamond" if use_diamond else "mmseqs"
    run_dir = out_root / sample.label / mode
    run_dir.mkdir(parents=True, exist_ok=True)

    cmd = build_command(
        runner,
        module,
        sample,
        threads=threads,
        use_smt=use_smt,
        nvme_scratch=nvme_scratch,
        use_diamond=use_diamond,
        common_args=common_args,
        extra_args=extra_args,
        out_dir=run_dir,
    )

    command_txt = run_dir / "command.sh"
    command_txt.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + shlex.join(cmd) + "\n")
    os.chmod(command_txt, 0o755)

    manifest = run_dir / "run.json"
    if dry_run:
        write_manifest(manifest, cmd=cmd, rc=0, sample=sample, use_diamond=use_diamond)
        print(f"[dry-run] {sample.label} {mode}: {shlex.join(cmd)}")
        return 0

    out_path = run_dir / "out.txt"
    err_path = run_dir / "err.txt"
    print(f"[{sample.label}][{mode}] starting")
    print(f"[{sample.label}][{mode}] {shlex.join(cmd)}")

    with out_path.open("w") as out_fh, err_path.open("w") as err_fh:
        proc = subprocess.run(cmd, stdout=out_fh, stderr=err_fh, text=True)

    write_manifest(manifest, cmd=cmd, rc=proc.returncode, sample=sample, use_diamond=use_diamond)

    if proc.returncode == 0:
        print(f"[{sample.label}][{mode}] finished OK")
    else:
        print(f"[{sample.label}][{mode}] failed with exit code {proc.returncode}", file=sys.stderr)

    return proc.returncode


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--runner", default="seqtoid-pipelines", help="pipeline executable name or path")
    p.add_argument("--module", default="short_read_mngs", help="pipeline module")
    p.add_argument("--out-root", required=True, help="directory where sample/mode subdirs will be created")
    p.add_argument("--threads", type=int, default=128)
    p.add_argument("--use-smt", action="store_true")
    p.add_argument("--nvme-scratch", required=True, type=Path)
    p.add_argument("--sample-sheet", type=Path, help="CSV/TSV with columns label,r1,r2")
    p.add_argument("--sample", action="append", default=[], help="LABEL=R1[,R2] or LABEL:R1[:R2]; repeatable")
    p.add_argument("--dry-run", action="store_true")

    # Required-ish inputs for your current pipeline example.
    p.add_argument("--kallisto-index")
    p.add_argument("--ercc-bowtie2-index")
    p.add_argument("--host-bowtie2-index")
    p.add_argument("--human-bowtie2-index")
    p.add_argument("--host-index")
    p.add_argument("--taxid-lineages-db")
    p.add_argument("--acc2taxid-db")
    p.add_argument("--nt-db-size")
    p.add_argument("--diamond-db")
    p.add_argument("--nt")
    p.add_argument("--nr")
    p.add_argument("--nt-offset-db")
    p.add_argument("--nr-offset-db")
    p.add_argument("--nt-info-tab")
    p.add_argument("--host-hisat2-index")
    p.add_argument("--human-hisat2-index")
    p.add_argument("--nt-split-dir")
    p.add_argument("--mmseqs-db")
    p.add_argument("--human-host", action="store_true")
    p.add_argument("--seed")
    p.add_argument("--max-subsample")
    p.add_argument("--min-read-len")
    p.add_argument("--max-read-len")
    p.add_argument("--max-reads")
    p.add_argument("--stall-threshold")

    p.add_argument(
        "--extra-arg",
        action="append",
        default=[],
        help="additional raw argument appended verbatim to every run; repeatable",
    )

    args = p.parse_args(argv)
    return args


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    out_root = Path(args.out_root).expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    nvme_scratch = args.nvme_scratch.expanduser().resolve()

    samples: list[Sample] = []
    if args.sample_sheet:
        samples.extend(_read_sample_sheet(args.sample_sheet.expanduser().resolve()))
    for spec in args.sample:
        samples.append(_parse_sample_spec(spec))

    if not samples:
        print("No samples provided. Use --sample-sheet or repeat --sample.", file=sys.stderr)
        return 2

    common_args: dict[str, Optional[str]] = {
        "kallisto_index": args.kallisto_index,
        "ercc_bowtie2_index": args.ercc_bowtie2_index,
        "host_bowtie2_index": args.host_bowtie2_index,
        "human_bowtie2_index": args.human_bowtie2_index,
        "host_index": args.host_index,
        "taxid_lineages_db": args.taxid_lineages_db,
        "acc2taxid_db": args.acc2taxid_db,
        "nt_db_size": args.nt_db_size,
        "diamond_db": args.diamond_db,
        "nt": args.nt,
        "nr": args.nr,
        "nt_offset_db": args.nt_offset_db,
        "nr_offset_db": args.nr_offset_db,
        "nt_info_tab": args.nt_info_tab,
        "host_hisat2_index": args.host_hisat2_index,
        "human_hisat2_index": args.human_hisat2_index,
        "nt_split_dir": args.nt_split_dir,
        "mmseqs_db": args.mmseqs_db,
        "human_host": args.human_host,
        "seed": args.seed,
        "max_subsample": args.max_subsample,
        "min_read_len": args.min_read_len,
        "max_read_len": args.max_read_len,
        "max_reads": args.max_reads,
        "stall_threshold": args.stall_threshold,
    }

    # Basic validation: only check the args that are commonly needed for the run.
    for key in (
            "kallisto_index",
            "ercc_bowtie2_index",
            "host_bowtie2_index",
            "host_hisat2_index",
            "taxid_lineages_db",
            "acc2taxid_db",
            "nt_db_size",
            "nt",
            "nr",
            "nt_offset_db",
            "nr_offset_db",
            "nt_split_dir",
    ):
        if not common_args.get(key):
            print(f"missing required pipeline arg: --{key.replace('_', '-')}", file=sys.stderr)
            return 2

    # Run Diamond first, then MMseqs.
    modes = [True, False]
    worst_rc = 0
    for sample in samples:
        # Verify inputs exist before kicking off any work.
        if not sample.r1.exists():
            print(f"missing R1 for {sample.label}: {sample.r1}", file=sys.stderr)
            worst_rc = 2
            continue
        if sample.r2 is not None and not sample.r2.exists():
            print(f"missing R2 for {sample.label}: {sample.r2}", file=sys.stderr)
            worst_rc = 2
            continue

        for use_diamond in modes:
            rc = run_one(
                runner=args.runner,
                module=args.module,
                sample=sample,
                out_root=out_root,
                threads=args.threads,
                use_smt=args.use_smt,
                nvme_scratch=nvme_scratch,
                use_diamond=use_diamond,
                common_args=common_args,
                extra_args=args.extra_arg,
                dry_run=args.dry_run,
            )
            worst_rc = max(worst_rc, rc)
            # Stop immediately on failure so runs stay sequential and easy to inspect.
            if rc != 0:
                return rc

    return worst_rc


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))