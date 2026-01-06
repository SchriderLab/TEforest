#!/usr/bin/env python3
"""
Generate paired alignment snapshots for each region in a BED,
expanding by 500 bp, from two BAMs, and stack the images.

Example
-------
$ python make_alignment_images.py \
    /nas/.../A1_TE_Annotation.filtered.bed \
    --path1 /path/to/first.bam \
    --path2 /path/to/second.bam \
    --workers 4
"""
import argparse
import os
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

def parse_args():
    p = argparse.ArgumentParser(description="Make stacked BamSnap images from two BAMs per BED region")
    p.add_argument("bed", type=Path,
                   help="Input BED file (4 cols: chr, start, end, ID, …, TE in col7)")
    p.add_argument("--path1", required=True, type=Path, help="First BAM file")
    p.add_argument("--path2", required=True, type=Path, help="Second BAM file")
    p.add_argument("-w", "--workers", type=int, default=os.cpu_count(),
                   help="Parallel workers (default: # cores)")
    return p.parse_args()

def read_bed(bed_path):
    regions = []
    with bed_path.open() as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                chrom, start, end, rid, *_rest = line.split()
                start, end = int(start), int(end)
                regions.append((chrom, max(0, start-500), end+500, rid))
    return regions

def make_commands(regions, bam1, bam2, out_dir):
    cmds = []
    for chrom, s, e, rid in regions:
        png1 = out_dir / f"{rid}.bam1.png"
        png2 = out_dir / f"{rid}.bam2.png"
        stacked = out_dir / f"{rid}.stacked.png"

        cmd1 = (
            f"bamsnap -bam {bam1} -pos {chrom}:{s}-{e} "
            f"-out {png1} -ref /nas/longleaf/home/adaigle/work/mcclintock_stuff/testing_truth_dataset/mcclintock_data/A1.fasta "
            "-bamplot coverage read -no_target_line -read_color_by interchrom"
        )
        cmd2 = cmd1.replace(str(bam1), str(bam2)).replace(str(png1), str(png2))
        # stack vertically via ImageMagick's `convert`
        cmd_stack = f"convert {png1} {png2} -append {stacked}"
        cmds.append((cmd1, cmd2, cmd_stack))
        # pass png1/png2 so we can delete them later
        cmds.append((cmd1, cmd2, cmd_stack, png1, png2))
    return cmds

def run_triple(cmds):
    failures = []
    for cmd1, cmd2, cmd_stack, png1, png2 in cmds:
        # 1) bam1
        if subprocess.run(cmd1, shell=True).returncode != 0:
            failures.append((cmd1, "exit != 0"))
            continue
        # 2) bam2
        if subprocess.run(cmd2, shell=True).returncode != 0:
            failures.append((cmd2, "exit != 0"))
            continue
        # 3) stack
        if subprocess.run(cmd_stack, shell=True).returncode != 0:
            failures.append((cmd_stack, "exit != 0"))
            continue
        # 4) cleanup intermediates
        try:
            png1.unlink()
            png2.unlink()
        except OSError as e:
            print(f"⚠️  Couldn’t delete intermediates {png1}, {png2}: {e}")
    return failures

def main():
    args = parse_args()

    out_dir = Path("/nas/longleaf/home/adaigle/work/mcclintock_stuff/"
                   "testing_truth_dataset/alignment_images")
    out_dir.mkdir(parents=True, exist_ok=True)

    if not args.bed.exists():
        raise FileNotFoundError(f"BED not found: {args.bed}")
    if not args.path1.exists() or not args.path2.exists():
        raise FileNotFoundError("One of the BAMs was not found")

    regions = read_bed(args.bed)
    cmds = make_commands(regions, args.path1, args.path2, out_dir)

    print(f"Launching {len(cmds)*3} commands ({len(cmds)} regions) with {args.workers} workers…")
    failures = []

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_triple, [triple]): triple for triple in cmds}
        for future in as_completed(futures):
            failed = future.result()
            if failed:
                failures.extend(failed)
                print(f"✘ Failure in region {futures[future]}")
            else:
                print(f"✔ Completed region {futures[future][0].split()[4]}")

    if failures:
        print(f"\n{len(failures)} command(s) failed; check above.")
    else:
        print("\nAll images generated successfully.")

if __name__ == "__main__":
    main()
