import argparse
import multiprocessing as mp
import os
import re
import subprocess
import sys
from io import StringIO
from typing import List, Union, Dict, Tuple
import gc
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import norm
import logging

df = pd.DataFrame
arr = np.ndarray


class EmptyRegionError(Exception):
    """Raised when data accessed from BAM file is empty for a given region."""
    pass


class EmtpyFeatVecError(Exception):
    """Raised when a feature vector is empty but the BAM region was not."""
    pass


def compute_edge_intervals(start: int, end: int, edge_len: int = 500) -> List[Tuple[int,int]]:
    """
    Return a list of genomic (start,end) intervals to KEEP.
    - If region_len >= 2*edge_len: two edges of length edge_len.
    - Else: the whole region.
    """
    if end - start >= 2 * edge_len:
        return [(start, start + edge_len), (end - edge_len, end)]
    else:
        return [(start, end)]


class SVMaker:
    def __init__(self, bamfile, output_dir=None):
        """
        Class for creating feature vectors from short read data.

        Args:
            output_dir (str/pathlike): Directory to write output npy files to.
            bamfile (str/pathlike): Bamfile to access for alignment
        """
        self.output_dir = output_dir
        self.bamfile = bamfile
        # Will store the local insert size mean/SD after computing
        self.local_insert_mean = None
        self.local_insert_sd   = None

    def _get_bam_data(
        self, bamfile, chromosome: str, start: int, end: int, threads: int = 1
    ) -> Union[df, None]:
        """Gets data from BAM file, parses features that require it, and returns a DataFrame."""
        samtools_cmd = [
            "samtools",
            "view",
            "-@ " + str(threads),
            bamfile,
            chromosome,
        ]
        if (start and end) is not None:
            samtools_cmd[-1] = chromosome + ":" + str(int(start)) + "-" + str(int(end))

        # Count how many reads are in this slice
        num_reads = int(subprocess.check_output(samtools_cmd + ["-c"]).decode("utf-8"))
        if num_reads == 0:
            return None

        a = subprocess.Popen(" ".join(samtools_cmd), stdout=subprocess.PIPE, shell=True)
        b = StringIO(a.communicate()[0].decode("utf-8"))

        try:
            bam_df = pd.read_csv(
                b,
                sep="\t",
                usecols=list(range(11)),
                names=[
                    "ID",
                    "Bitflag",
                    "Chrom",
                    "Read_Start",
                    "Quality",
                    "CIGAR",
                    "Rnext",
                    "Pnext",
                    "Template_Length",
                    "Sequence",
                    "PHRED",
                ],
            )
        except:
            return None

        # Expand bitflags
        bitflag_df = df(
            self._parse_bitflag(list(bam_df["Bitflag"]), len(bam_df)),
            columns=[
                "Paired",
                "Proper_Pair",
                "Is_Read1_Unmapped",
                "Is_Read2_Unmapped",
                "Is_Read1_Rev_Comp",
                "Is_Read2_Rev_Comp",
                "Is_First_Read",
                "Is_Second_Read",
            ],
        )
        bam_df = pd.concat([bam_df, bitflag_df], axis=1)
        bam_df = self._filter_dataframe(bam_df)
        bam_df = self._classify_reads(bam_df)

        # Compute local insert stats (mean, SD) from "proper" pairs
        self._compute_local_insert_stats(bam_df)

        return bam_df

    def _compute_local_insert_stats(self, bam_df: df):
        """
        Computes the local insert size mean and SD for the region,
        using only reads that are properly paired (Proper_Pair=1),
        both mapped (Is_Read1_Unmapped=0, Is_Read2_Unmapped=0),
        and in correct orientation.
        """
        filtered = bam_df[
            (bam_df["Proper_Pair"] == 1)
            & (bam_df["Is_Read1_Unmapped"] == 0)
            & (bam_df["Is_Read2_Unmapped"] == 0)
            & (bam_df["Rnext"] == "=")
        ]

        vals = filtered["Template_Length"].abs()
        if len(vals) > 1:
            self.local_insert_mean = vals.mean()
            self.local_insert_sd   = vals.std()
        else:
            # Fallback if too few reads
            self.local_insert_mean = 300.0
            self.local_insert_sd   = 100.0

    def _filter_dataframe(self, region_data: df) -> df:
        region_data = region_data.drop_duplicates(ignore_index=True)
        region_data = region_data.dropna()
        region_data = region_data[region_data["CIGAR"] != "*"].reset_index(drop=True)
        return region_data

    def _classify_reads(self, bam_df: df) -> df:
        try:
            bam_df["Read_Length"] = bam_df["Sequence"].apply(len)
            read_len = stats.mode(np.abs(bam_df["Read_Length"]), keepdims=True).mode[0]
        except:
            return None

        bam_df["Split"] = np.where(
            np.abs(bam_df["Read_Length"]) > read_len + 5, 1, 0
        )
        bam_df["Orphan"] = np.where(bam_df["Rnext"] != "=", 1, 0)

        template_threshold_99 = abs(bam_df["Template_Length"]).quantile(0.99)
        template_threshold_02 = abs(bam_df["Template_Length"]).quantile(0.02)

        bam_df["Long_Insert"] = np.where(
            bam_df["Template_Length"] >= template_threshold_99, 1, 0
        )
        bam_df["Short_Insert"] = np.where(
            bam_df["Template_Length"] <= template_threshold_02, 1, 0
        )

        # Parallel/everted/orphan reads
        bam_df["Parallel_Read"] = np.where(
            (
                (bam_df["Is_Read1_Rev_Comp"] == 0)
                & (bam_df["Is_Read2_Rev_Comp"] == 0)
                & (bam_df["Is_Read1_Unmapped"] == 0)
                & (bam_df["Is_Read2_Unmapped"] == 0)
            )
            | (
                (bam_df["Is_Read1_Rev_Comp"] == 1)
                & (bam_df["Is_Read2_Rev_Comp"] == 1)
                & (bam_df["Is_Read1_Unmapped"] == 0)
                & (bam_df["Is_Read2_Unmapped"] == 0)
            ),
            1,
            0,
        )
        bam_df["Everted_Read"] = np.where(
            (
                (bam_df["Is_First_Read"] == 0)
                & ((bam_df["Rnext"] == "=") & (bam_df["Read_Start"] > bam_df["Pnext"]))
                & ((bam_df["Is_Read1_Rev_Comp"] == 0) & (bam_df["Is_Read2_Rev_Comp"] == 1))
            )
            | (
                (bam_df["Is_First_Read"] == 1)
                & ((bam_df["Rnext"] == "=") & (bam_df["Read_Start"] < bam_df["Pnext"]))
                & ((bam_df["Is_Read1_Rev_Comp"] == 0) & (bam_df["Is_Read2_Rev_Comp"] == 1))
            ),
            1,
            0,
        )
        bam_df["Orphan_Read"] = np.where(
            (bam_df["Orphan"] == 1)
            | ((bam_df["Is_Read1_Unmapped"] == 1) | (bam_df["Is_Read2_Unmapped"] == 1)),
            1,
            0,
        )
        return bam_df

    def _parse_bitflag(self, bitflag_list: List[str], num_reads: int) -> arr:
        bitflag_array = np.zeros((num_reads, 8))
        for i in range(len(bitflag_list)):
            bitflag_str = format(int(bitflag_list[i]), "#016b")[::-1][0:8]
            bitflag_array[i] = list(bitflag_str)
        return bitflag_array

    def _summarize_reads(self, reads_arr: arr) -> arr:
        return np.sum(reads_arr, axis=1) if reads_arr.ndim == 3 else reads_arr

    def _create_sum_stats(self, summ_arr: arr, read_counts: arr) -> arr:
        """
        Creates array of summary statistics for region.
        For each feature:
            Mean, Median, StdDev, IQR
        Uses empirical quantile resampling (100 points).
        """
        feat_vec = np.zeros(summ_arr.shape[0] * 4)
        read_counts_safe = np.where(read_counts == 0, 1, read_counts)

        for feat_ind in range(summ_arr.shape[0]):
            _stand_feat = summ_arr[feat_ind, :]
            if len(_stand_feat) == 0:
                return None

            # 1) Normalize by coverage
            normalized_feat = _stand_feat / read_counts_safe
            normalized_feat = np.nan_to_num(normalized_feat)

            # 2) Empirical CDF
            cumulative_coverage = np.cumsum(normalized_feat)
            if cumulative_coverage[-1] == 0:
                cumulative_coverage[-1] = 1
            empirical_cdf = cumulative_coverage / cumulative_coverage[-1]

            # 3) Resample onto fixed grid
            new_points = 100
            new_x = np.linspace(0, 1, new_points)
            resampled_feat = np.interp(new_x, empirical_cdf, normalized_feat)

            # 4) Summary stats
            vec_ind = feat_ind * 4
            feat_vec[vec_ind + 0] = np.mean(resampled_feat)
            feat_vec[vec_ind + 1] = np.median(resampled_feat)
            feat_vec[vec_ind + 2] = np.std(resampled_feat)
            feat_vec[vec_ind + 3] = stats.iqr(resampled_feat, interpolation="midpoint")

        return feat_vec

    def _set_outfile_name(self, region: pd.Series) -> str:
        outfile = "-".join(
            [
                region["Sample"],
                region["Chrom"],
                str(int(region["Ref_begin"])),
                str(int(region["Ref_end"])),
                str(int(region["Class"])),
                region["TE"],
            ]
        )
        return outfile

    def _set_reg_str(self, region: pd.Series):
        self.reg_str = "-".join(
            [
                region["Chrom"],
                str(region["Ref_begin"]),
                str(region["Ref_end"]),
            ]
        )

    def _check_create_dir(self, region: pd.Series):
        outdir = os.path.join(self.output_dir, region["Sample"])
        os.makedirs(outdir, exist_ok=True)

    def _save_featvec(self, outfile_name: str, region: pd.Series, feat_vec: arr):
        np.save(
            os.path.join(self.output_dir, region["Sample"], outfile_name + ".npy"),
            feat_vec,
        )

    def _save_readcount(self, outfile_name: str, region: pd.Series, read_counts: arr):
        np.save(
            os.path.join(self.output_dir, region["Sample"], outfile_name + "_readcount.npy"),
            read_counts,
        )

    def create_region(
        self,
        region: pd.Series,
        bamfile: str,
        savefile: bool = True,
        edge_len: int = 500,
    ):
        """
        Create a feature vector using only the edges:
          - If region_len >= 2*edge_len: only two  edge_len windows
          - Else: full region
        The center is dropped entirely.
        """
        self._set_reg_str(region)
        try:
            bam_df = self._get_bam_data(
                bamfile,
                region["Chrom"],
                region["Ref_begin"],
                region["Ref_end"],
            )
            if bam_df is None:
                raise EmptyRegionError

            reg_start = int(region["Ref_begin"])
            reg_end   = int(region["Ref_end"])
            start_diff = reg_start

            # -------------------------
            # compute intervals to keep
            keep_ivals = compute_edge_intervals(reg_start, reg_end, edge_len=edge_len)

            mapper = ReadMapper(start_diff, reg_start, reg_end)
            pileup_arr, read_counts = mapper.map_reads_with_intervals(
                bam_df, keep_ivals, compress=True
            )
            if pileup_arr is None:
                raise EmtpyFeatVecError

            summary_arr = self._summarize_reads(pileup_arr)
            feat_vec = self._create_sum_stats(summary_arr, read_counts)

            if savefile and feat_vec is not None:
                outfile = self._set_outfile_name(region)
                self._check_create_dir(region)
                self._save_featvec(outfile, region, feat_vec)
                self._save_readcount(outfile, region, read_counts)
                print(outfile + " (edges-only) made.")
            else:
                return feat_vec, read_counts

        except Exception:
            logging.warning(
                self.reg_str + " couldn't be loaded (edges-only).", exc_info=True
            )
            return None

    def update_region(
        self,
        region: pd.Series,
        bamfile: str,
        savefile: bool = True,
        edge_len: int = 500,
    ):
        """
        Update an existing feature vector by computing TE‑specific features
        using the same edge-only rule.
        """
        self._set_reg_str(region)
        try:
            bam_df = self._get_bam_data(
                bamfile,
                region["Chrom"],
                region["Ref_begin"],
                region["Ref_end"],
            )
            if bam_df is None:
                raise EmptyRegionError

            reg_start = int(region["Ref_begin"])
            reg_end   = int(region["Ref_end"])
            start_diff = int(bam_df["Read_Start"].min())

            # -------------------------
            keep_ivals = compute_edge_intervals(reg_start, reg_end, edge_len=edge_len)

            mapper = ReadMapper(start_diff, reg_start, reg_end)
            pileup_arr, read_counts = mapper.map_reads_with_intervals(
                bam_df, keep_ivals, compress=True
            )
            if pileup_arr is None:
                raise EmtpyFeatVecError

            summary_arr = self._summarize_reads(pileup_arr)
            feat_vec = self._create_sum_stats(summary_arr, read_counts)

            if savefile and feat_vec is not None:
                outfile = self._set_outfile_name(region)
                first_vec = np.load(
                    os.path.join(self.output_dir, region["Sample"], outfile + ".npy")
                )
                combined = np.concatenate((first_vec, feat_vec), axis=0)
                self._save_featvec(outfile, region, combined)
            else:
                return feat_vec

        except Exception:
            logging.warning(
                self.reg_str + " couldn't be updated (edges-only).", exc_info=True
            )
            return None


class ReadMapper:
    def __init__(self, start_diff: int, start: int, end: int):
        """
        Maps reads onto a 2D matrix: (features, bases_in_kept_intervals).
        """
        self.start_diff = start_diff
        self.start = start
        self.end   = end
        self.cigar_re = re.compile(r'(\d+)([MIDNSHP=X])')

        self.flag_feats = [
            "Paired", "Proper_Pair", "Is_Read1_Unmapped", "Is_Read2_Unmapped",
            "Is_Read1_Rev_Comp", "Is_Read2_Rev_Comp", "Is_First_Read",
            "Is_Second_Read", "Split", "Long_Insert", "Short_Insert",
            "Parallel_Read", "Everted_Read", "Orphan_Read"
        ]

    def map_reads_with_intervals(
        self,
        reads_df: df,
        keep_intervals: List[Tuple[int, int]],
        compress: bool = True,
    ) -> Union[None, Tuple[arr, np.ndarray]]:
        """
        Map reads but only into the union of 'keep_intervals'.
        If compress=True, the dest array length = sum(interval_lengths).
        """
        # 1) get region-relative starts
        df_rel = reads_df.copy()
        df_rel["Read_Start"] = df_rel["Read_Start"] - self.start
        df_rel = df_rel.reset_index(drop=True)

        # 2) compute relative intervals
        rel_ivals = [(s - self.start, e - self.start) for s, e in keep_intervals]
        # drop invalid
        rel_ivals = [(s, e) for s, e in rel_ivals if e > s]
        if not rel_ivals:
            return None

        seg_lens = [e - s for s, e in rel_ivals]
        region_len = self.end - self.start
        total_len = sum(seg_lens) if compress else region_len

        # 3) prepare accumulator
        pileup = np.zeros((22, total_len), dtype=np.int32)
        # compute offsets for each segment in dest
        offsets = []
        off = 0
        for L in seg_lens:
            offsets.append(off)
            off += L

        # 4) iterate reads
        for _, row in df_rel.iterrows():
            binary_cig = self._parse_cigar(row["CIGAR"])
            if binary_cig.ndim == 3:
                binary_cig = binary_cig[:, 0, :]
            read_len = binary_cig.shape[1]

            rstart = int(row["Read_Start"])
            rend   = rstart + read_len

            # for each kept segment
            for idx, (seg_s, seg_e) in enumerate(rel_ivals):
                o0 = max(rstart, seg_s)
                o1 = min(rend,   seg_e)
                if o1 <= o0:
                    continue

                c0 = o0 - rstart
                c1 = c0 + (o1 - o0)
                if compress:
                    d0 = offsets[idx] + (o0 - seg_s)
                else:
                    d0 = o0
                d1 = d0 + (o1 - o0)

                # CIGAR coverage (5 channels)
                pileup[0:5, d0:d1] += binary_cig[:, c0:c1].astype(np.int32, copy=False)
                # Flags
                for j, feat in enumerate(self.flag_feats, start=5):
                    if row[feat]:
                        pileup[j, d0:d1] += 1
                # Quality, template-length, read-count
                pileup[19, d0:d1] += int(row["Quality"])
                pileup[20, d0:d1] += abs(int(row["Template_Length"]))
                pileup[21, d0:d1] += 1

        read_counts = pileup[21, :].astype(np.int32, copy=False)
        return pileup[:21, :].astype(np.float32, copy=False), read_counts

    def map_reads(self, reads_df: df) -> Union[None, Tuple[arr, arr]]:
        # (left intact, if you ever need full-region fallback)
        reads_df["Read_Start"] = reads_df["Read_Start"] - self.start
        reads_df = reads_df.reset_index(drop=True)
        region_len = self.end - self.start

        if reads_df.shape[0] > 0:
            pileup_arr = np.zeros((22, region_len), dtype=np.int32)
            def _add_flags_inplace(acc, flag_feats, row, region_slice):
                for j, feat in enumerate(flag_feats, start=5):
                    if row[feat]:
                        acc[j, region_slice] += 1

            for _, row in reads_df.iterrows():
                binary_cig = self._parse_cigar(row["CIGAR"])
                if binary_cig.ndim == 3:
                    binary_cig = binary_cig[:, 0, :]
                start      = row["Read_Start"]
                end        = start + binary_cig.shape[1]

                r0 = max(0, start)
                r1 = min(end, region_len)
                if r1 <= r0:
                    continue

                c0 = r0 - start
                c1 = c0 + (r1 - r0)

                pileup_arr[0:5, r0:r1] += binary_cig[:, c0:c1]
                _add_flags_inplace(pileup_arr, self.flag_feats, row, slice(r0, r1))
                pileup_arr[19, r0:r1] += row["Quality"]
                pileup_arr[20, r0:r1] += abs(row["Template_Length"])
                pileup_arr[21, r0:r1] += 1

            read_counts_per_base = pileup_arr[21, :].astype(np.int32, copy=False)
            pileup_arr_no_count = pileup_arr[:21, :].astype(np.float32, copy=False)
            return pileup_arr_no_count, read_counts_per_base
        else:
            return None

    def _parse_cigar(self, cig_string: str) -> arr:
        cig_list = re.split("([a-zA-Z])", cig_string)[0:-1]
        cig_expanded = ""
        for i in range(0, len(cig_list), 2):
            length = int(cig_list[i])
            op     = cig_list[i+1]
            cig_expanded += length * op

        binary_cig_arr = np.zeros((5, 1, len(cig_expanded)))
        cig_opts = ["M", "D", "I", "S", "H"]
        for i, c in enumerate(cig_expanded):
            if c in cig_opts:
                binary_cig_arr[cig_opts.index(c), 0, i] = 1
        return binary_cig_arr


def region_generator(sub_regions, output_dir, bam_dir):
    for idx, region in sub_regions.iterrows():
        bamfile = os.path.join(bam_dir, region["Sample"] + ".bam")
        svm = SVMaker(bamfile, output_dir)
        svm.create_region(region, bamfile, True)
        if idx % 1000 == 0:
            gc.collect()


def te_specific_region_generator(sub_regions, output_dir, te_bam_dir):
    for idx, region in sub_regions.iterrows():
        bamfile = os.path.join(
            te_bam_dir, region["Sample"] + "/" + region["TE"] + "_to_ISO1.bam"
        )
        svm = SVMaker(bamfile, output_dir)
        svm.update_region(region, bamfile, True)
        if idx % 1000 == 0:
            gc.collect()


def fused_region_generator(sub_regions, output_dir, bam_dir, te_bam_dir):
    """
    One-pass generator (edges-only):
      1. build the base feature vector,
      2. build the TE-specific vector,
      3. save the combined result.
    """
    for idx, region in sub_regions.iterrows():
        bamfile = os.path.join(bam_dir, region["Sample"] + ".bam")
        te_bam  = os.path.join(
            te_bam_dir,
            region["Sample"],
            f'{region["TE"]}_to_ISO1.bam'
        )

        svm = SVMaker(bamfile, output_dir)
        result = svm.create_region(region, bamfile, savefile=False)
        if result is None:
            continue
        base_vec, read_counts = result

        te_vec = svm.update_region(region, te_bam, savefile=False)
        if te_vec is None:
            continue

        full_vec = np.concatenate((base_vec, te_vec), axis=0)
        outfile  = svm._set_outfile_name(region)
        svm._check_create_dir(region)

        np.save(
            os.path.join(output_dir, region["Sample"], outfile + ".npy"),
            full_vec,
        )
        np.save(
            os.path.join(output_dir, region["Sample"], outfile + "_readcount.npy"),
            read_counts,
        )

        if idx % 1000 == 0:
            gc.collect()


def parse_arguments() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i","--input_file",metavar="IN_FILE",required=False,dest="input_file",
        default="../test/testlocs.csv",
        help="5-column csv with colname header SAMPLE,CHROM,REF_BEGIN,REF_END,CLASS,TE"
    )
    parser.add_argument(
        "-od","--output_dir",metavar="OUTPUT_DIR",required=False,default=".",
        dest="output_dir",help="Directory to save npy files to. Defaults to cwd."
    )
    parser.add_argument(
        "-bd","--bam_dir",metavar="BAM_DIR",type=str,required=False,
        default="TrainingData/bamfiles/",dest="bam_dir",
        help="Path to BAM file directory. Filenames must match SAMPLE col."
    )
    parser.add_argument(
        "-tv","--train_val",metavar="TRAIN_VAL",type=str,required=False,default="train",
        dest="train_val",help="Store data in /train/ or /val/ subdirs."
    )
    parser.add_argument(
        "--stop_idx",metavar="STOP_INDEX",type=str,required=False,dest="stop_idx",
        help="Where to stop indexing the dataframe for parallelization."
    )
    parser.add_argument(
        "-tebd","--te_bam_dir",metavar="TE_BAM_DIR",type=str,required=False,
        default="TrainingData/bamfiles/",dest="te_bam_dir",
        help="Path to TE BAM file directory, TE name must match region TE col."
    )
    parser.add_argument(
        "-p", "--processes", metavar="N_PROC", type=int, default=1,
        dest="n_proc",
        help="How many worker processes to launch (Snakemake passes {threads})."
    )
    return parser.parse_args()


def main():
    ua = parse_arguments()
    regions = pd.read_csv(ua.input_file)
    N = max(1, int(ua.n_proc))
    print("Using", N, "processes")

    def chunkify(df, n):
        return np.array_split(df, n)

    if N == 1:
        region_generator(regions, ua.output_dir, ua.bam_dir)
        te_specific_region_generator(regions, ua.output_dir, ua.te_bam_dir)
    else:
        print("Spawning", N, "workers for fused_region_generator (edges-only)")
        with mp.Pool(N) as pool:
            pool.starmap(
                fused_region_generator,
                [
                    (chunk, ua.output_dir, ua.bam_dir, ua.te_bam_dir)
                    for chunk in chunkify(regions, N)
                ],
            )


if __name__ == "__main__":
    main()
