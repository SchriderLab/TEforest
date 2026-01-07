#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert TEforest BED calls to a single-sample SV VCF."
    )
    parser.add_argument("--bed", required=True, help="Input TEforest BED file.")
    parser.add_argument("--ref", required=True, help="Reference FASTA for REF base.")
    parser.add_argument("--output", required=True, help="Output VCF path.")
    parser.add_argument(
        "--sample",
        default=None,
        help="Override sample name (otherwise inferred from TE_string).",
    )
    return parser.parse_args()


def ensure_fai(ref_path: str) -> str:
    fai_path = ref_path + ".fai"
    if not os.path.exists(fai_path):
        subprocess.run(["samtools", "faidx", ref_path], check=True)
    return fai_path


class FastaIndex:
    def __init__(self, fasta_path: str):
        self.fasta_path = fasta_path
        self.fai_path = ensure_fai(fasta_path)
        self.index = {}
        with open(self.fai_path, "r", encoding="utf-8") as fh:
            for line in fh:
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                name = fields[0]
                length = int(fields[1])
                offset = int(fields[2])
                line_bases = int(fields[3])
                line_width = int(fields[4])
                self.index[name] = (length, offset, line_bases, line_width)
        self._fh = open(self.fasta_path, "rb")

    def get_base(self, chrom: str, pos1: int) -> str:
        if chrom not in self.index:
            return "N"
        length, offset, line_bases, line_width = self.index[chrom]
        if pos1 < 1 or pos1 > length:
            return "N"
        pos0 = pos1 - 1
        line = pos0 // line_bases
        col = pos0 % line_bases
        file_offset = offset + line * line_width + col
        self._fh.seek(file_offset)
        base = self._fh.read(1).decode("ascii", errors="ignore").upper()
        return base if base in {"A", "C", "G", "T", "N"} else "N"

    def close(self) -> None:
        self._fh.close()


def normalize_zyg(val: str) -> float:
    try:
        z = float(val)
    except (TypeError, ValueError):
        return 0.0
    if z >= 0.75:
        return 1.0
    if z >= 0.25:
        return 0.5
    return 0.0


def gt_from_zyg(zyg: float, teorig: str) -> str:
    if teorig == "non-reference":
        mapping = {1.0: "1/1", 0.5: "0/1", 0.0: "0/0"}
    elif teorig == "reference":
        mapping = {1.0: "0/0", 0.5: "0/1", 0.0: "1/1"}
    else:
        raise ValueError(f"Unknown TEORIG: {teorig}")
    return mapping[zyg]


def parse_te_string(te_string: str) -> dict:
    parts = te_string.split("|")
    if len(parts) < 7:
        raise ValueError(f"Invalid TE_string: {te_string}")
    return {
        "TEFAM": parts[0],
        "TEORIG": parts[1],
        "ZYG": parts[2],
        "SAMPLE": parts[3],
        "CALLER": parts[4],
        "EVID": parts[5],
        "IDNUM": parts[6],
    }


def vcf_header(sample: str, include_sample_info: bool) -> str:
    lines = [
        "##fileformat=VCFv4.2",
        "##source=TEforest",
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">',
        '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
        '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">',
        '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">',
        '##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info: NAME,START,END,POLARITY">',
        '##INFO=<ID=TEORIG,Number=1,Type=String,Description="TE origin (reference or non-reference)">',
        '##INFO=<ID=CALLER,Number=1,Type=String,Description="Calling method">',
        '##INFO=<ID=EVID,Number=1,Type=String,Description="Evidence tag">',
        '##INFO=<ID=IDNUM,Number=1,Type=Integer,Description="Record index from TEforest">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]
    if include_sample_info:
        lines.append(
            '##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name (only set when mixed samples are present)">'
        )
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample)
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    records = []
    samples = set()

    with open(args.bed, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0].lower() in {"seqnames", "chrom", "chr", "#chrom"}:
                continue
            if len(fields) < 6:
                continue
            chrom, start, end, te_string, _score, _strand = fields[:6]
            te = parse_te_string(te_string)
            samples.add(te["SAMPLE"])
            records.append((chrom, int(start), int(end), te))

    if not records:
        with open(args.output, "w", encoding="utf-8") as out:
            sample_name = args.sample or "SAMPLE"
            out.write(vcf_header(sample_name, include_sample_info=False))
        return

    if args.sample:
        sample_name = args.sample
        include_sample_info = False
    else:
        if len(samples) == 1:
            sample_name = next(iter(samples))
            include_sample_info = False
        else:
            sample_name = "SAMPLE"
            include_sample_info = True
            print(
                f"[bed_to_vcf] Multiple sample names found: {sorted(samples)}",
                file=sys.stderr,
            )

    fasta = FastaIndex(args.ref)
    with open(args.output, "w", encoding="utf-8") as out:
        out.write(vcf_header(sample_name, include_sample_info))
        for chrom, start, end, te in records:
            tefam = te["TEFAM"]
            teorig = te["TEORIG"]
            zyg = normalize_zyg(te["ZYG"])
            gt = gt_from_zyg(zyg, teorig)

            pos = start + 1  # BED is 0-based; VCF is 1-based.
            ref_base = fasta.get_base(chrom, pos)
            length = max(0, end - start)

            if teorig == "non-reference":
                svtype = "INS"
                alt = f"<INS:ME:{tefam}>"
                svlen = length
            elif teorig == "reference":
                svtype = "DEL"
                alt = f"<DEL:ME:{tefam}>"
                svlen = -length
            else:
                raise ValueError(f"Unknown TEORIG: {teorig}")

            try:
                idnum_int = int(te["IDNUM"])
                id_field = f"TEforest_{idnum_int}"
                idnum_value = str(idnum_int)
            except (TypeError, ValueError):
                id_field = "."
                idnum_value = "."

            info_parts = [
                f"SVTYPE={svtype}",
                f"END={end}",
                f"SVLEN={svlen}",
                "IMPRECISE",
                "CIPOS=0,0",
                "CIEND=0,0",
                f"MEINFO={tefam},.,.,.",
                f"TEORIG={teorig}",
                f"CALLER={te['CALLER']}",
                f"EVID={te['EVID']}",
                f"IDNUM={idnum_value}",
            ]
            if include_sample_info:
                info_parts.append(f"SAMPLE={te['SAMPLE']}")

            info_str = ";".join(info_parts)

            row = [
                chrom,
                str(pos),
                id_field,
                ref_base,
                alt,
                ".",
                "PASS",
                info_str,
                "GT",
                gt,
            ]
            out.write("\t".join(row) + "\n")

    fasta.close()


if __name__ == "__main__":
    main()
