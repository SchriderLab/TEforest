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
            self._parse_bitflag(list(bam_df["Bitflag"]), num_reads),
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
        and in correct orientation (Rnext='=' for the same reference, etc. if desired).
        """
        # Filter for properly paired, uniquely mapped reads
        filtered = bam_df[
            (bam_df["Proper_Pair"] == 1)
            & (bam_df["Is_Read1_Unmapped"] == 0)
            & (bam_df["Is_Read2_Unmapped"] == 0)
            & (bam_df["Rnext"] == "=")  # same reference
        ]

        # Use absolute value of Template_Length to compute stats
        vals = filtered["Template_Length"].abs()
        if len(vals) > 1:
            self.local_insert_mean = vals.mean()
            self.local_insert_sd   = vals.std()
        else:
            # Fallback if no or too few reads match the filter
            self.local_insert_mean = 300.0
            self.local_insert_sd   = 100.0

    def _filter_dataframe(self, region_data: df) -> df:
        region_data = region_data.drop_duplicates(ignore_index=True)
        region_data = region_data.dropna()
        region_data = region_data[region_data["CIGAR"] != "*"].reset_index()
        return region_data

    def _classify_reads(self, bam_df: df) -> df:
        try:
            bam_df["Read_Length"] = bam_df["Sequence"].apply(len)
            read_len = stats.mode(np.abs(bam_df["Read_Length"]), keepdims=True).mode[0]
        except:
            return None

        # Mark reads that are significantly longer than the mode (possible split reads)
        bam_df["Split"] = np.where(
            np.abs(bam_df["Read_Length"]) > read_len + 5, 1, 0,
        )
        bam_df["Orphan"] = np.where(bam_df["Rnext"] != "=", 1, 0)

        # For reference
        template_threshold_99 = abs(bam_df["Template_Length"]).quantile(0.99)
        template_threshold_02 = abs(bam_df["Template_Length"]).quantile(0.02)

        bam_df["Long_Insert"] = np.where(
            bam_df["Template_Length"] >= template_threshold_99, 1, 0
        )
        bam_df["Short_Insert"] = np.where(
            bam_df["Template_Length"] <= template_threshold_02, 1, 0
        )

        # Parallel, everted, orphan read flags
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
        return np.sum(reads_arr, axis=1)

    def _create_sum_stats(self, summ_arr: arr, read_counts: arr) -> arr:
        """
        Creates array of summary statistics for region.
        For each feature:
            Mean, Median, StdDev, IQR
    
        Instead of Gaussian rescaling, we now use empirical quantile-based rescaling.
        This allows us to normalize feature distributions regardless of insert size SD or region length.
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
    
            # 2) Compute cumulative sum and normalize to get empirical quantile coordinates
            cumulative_coverage = np.cumsum(normalized_feat)
            if cumulative_coverage[-1] == 0:
                cumulative_coverage[-1] = 1  # avoid division by zero
            empirical_cdf = cumulative_coverage / cumulative_coverage[-1]
    
            # 3) Resample the feature onto a fixed coordinate grid
            new_points = 100
            new_x = np.linspace(0, 1, new_points)
            resampled_feat = np.interp(new_x, empirical_cdf, normalized_feat)
    
            # 4) Compute summary stats
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
        if not os.path.exists(outdir):
            print(outdir)
            os.makedirs(outdir)

    def _save_featvec(self, outfile_name: str, region: pd.Series, feat_vec: arr):
        np.save(
            os.path.join(self.output_dir, region["Sample"], outfile_name + ".npy"),
            feat_vec,
        )

    def _save_readcount(self, outfile_name: str, region: pd.Series, feat_vec: arr):
        np.save(
            os.path.join(self.output_dir, region["Sample"], outfile_name + "_readcount.npy"),
            feat_vec,
        )

    def create_region(self, region: pd.Series, bamfile: str, savefile: bool = True):
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

            mapper = ReadMapper(start_diff, reg_start, reg_end)
            pileup_arr, read_counts_per_base = mapper.map_reads(bam_df)
            if pileup_arr is None:
                raise EmtpyFeatVecError

            summary_arr = self._summarize_reads(pileup_arr)
            feat_vec = self._create_sum_stats(summary_arr, read_counts_per_base)

            if savefile and feat_vec is not None:
                outfile = self._set_outfile_name(region)
                if not os.path.exists(
                    os.path.join(self.output_dir, region["Sample"], outfile)
                ):
                    self._check_create_dir(region)
                    self._save_featvec(outfile, region, feat_vec)
                    self._save_readcount(outfile, region, read_counts_per_base)
                    print(outfile + " made.")
                else:
                    print(outfile + " already exists. Passing...")
            else:
                return feat_vec

        except Exception:
            logging.warning(
                self.reg_str + f" couldn't be loaded.",
                exc_info=True,
            )
            return None

    def update_region(self, region: pd.Series, bamfile: str, savefile: bool = True):
        """
        Instead of creating a fresh file, we load an existing feature vector and
        concatenate new features onto it (used for TE-specific updates).
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

            mapper = ReadMapper(start_diff, reg_start, reg_end)
            pileup_arr, read_counts_per_base = mapper.map_reads(bam_df)
            if pileup_arr is None:
                raise EmtpyFeatVecError

            summary_arr = self._summarize_reads(pileup_arr)
            outfile = self._set_outfile_name(region)
            # Load existing read counts
            read_counts_per_base = np.load(
                os.path.join(self.output_dir, region["Sample"], outfile + "_readcount.npy")
            )
            feat_vec = self._create_sum_stats(summary_arr, read_counts_per_base)

            if savefile and feat_vec is not None:
                first_feat_vec = np.load(
                    os.path.join(self.output_dir, region["Sample"], outfile + ".npy")
                )
                second_feat_vec = np.concatenate((first_feat_vec, feat_vec), axis=0)
                self._save_featvec(outfile, region, second_feat_vec)
            else:
                return feat_vec

        except Exception:
            logging.warning(
                self.reg_str + f" couldn't be loaded.",
                exc_info=True,
            )
            return None


class ReadMapper:
    def __init__(self, start_diff: int, start: int, end: int):
        """
        Maps reads onto a 3D matrix: (features, reads, bases).
        """
        self.start_diff = start_diff
        self.start = start
        self.end   = end

    def map_reads(self, reads_df: df) -> Union[None, Tuple[arr, arr]]:
        """
        Build a 3D matrix of shape (22, n_reads, region_length).
        Then sum across reads to get the per-base coverage for each feature.
        """
        reads_df["Read_Start"] = reads_df["Read_Start"] - self.start_diff
        reads_df = reads_df.reset_index()
        region_len = self.end - self.start

        if reads_df.shape[0] > 0:
            pileup_arr = np.zeros((22, reads_df.shape[0], region_len))
            max_chunk_index = region_len - 1

            for i in range(pileup_arr.shape[1]):
                read_start_ind = reads_df["Read_Start"][i]
                binary_cig_arr = self._parse_cigar(reads_df["CIGAR"][i])
                _read_inds = list(range(read_start_ind, read_start_ind + binary_cig_arr.shape[2]))
                read_inds = [r for r in _read_inds if r >= 0]

                if len(read_inds) == 0 or read_inds[0] >= max_chunk_index:
                    continue

                read_inds, pileup_arr = self._map_cigar(pileup_arr, i, read_inds, max_chunk_index, binary_cig_arr)
                pileup_arr = self._iterate_feats(reads_df, pileup_arr, i, read_inds)

                # Quality
                pileup_arr[19, i, read_inds] = reads_df["Quality"][i]
                # Template length
                pileup_arr[20, i, read_inds] = np.abs(reads_df["Template_Length"][i])
                # "Read count" feature
                pileup_arr[21, i, read_inds] = 1

            # read_counts = sum of last dimension
            read_counts_per_base = np.sum(pileup_arr[21, :, :], axis=0)
            # remove the "read count" row
            pileup_arr = pileup_arr[:-1, :, :]
            return pileup_arr, read_counts_per_base
        else:
            return None

    def _map_cigar(self, pileup_arr, read_idx, read_inds, max_chunk_index, binary_cig_arr):
        if read_inds[-1] >= max_chunk_index:
            read_inds = read_inds[: read_inds.index(max_chunk_index)]
            pileup_arr[0:5, read_idx : read_idx+1, read_inds] = binary_cig_arr[:, :, 0 : len(read_inds)]
        else:
            pileup_arr[0:5, read_idx : read_idx+1, read_inds] = binary_cig_arr
        return read_inds, pileup_arr

    def _parse_cigar(self, cig_string: str) -> arr:
        cig_list = re.split("([a-zA-Z])", cig_string)[0:-1]
        cig_expanded = ""
        for i in range(0, len(cig_list), 2):
            length = int(cig_list[i])
            op     = cig_list[i+1]
            cig_expanded += length * op

        # 5 possible ops: M, D, I, S, H
        binary_cig_arr = np.zeros((5, 1, len(cig_expanded)))
        cig_opts = ["M", "D", "I", "S", "H"]
        for i, c in enumerate(cig_expanded):
            if c in cig_opts:
                binary_cig_arr[cig_opts.index(c), 0, i] = 1
        return binary_cig_arr

    def _iterate_feats(self, reads_df, pileup_arr, read_idx, read_inds):
        featlist = [
            "Paired","Proper_Pair","Is_Read1_Unmapped","Is_Read2_Unmapped",
            "Is_Read1_Rev_Comp","Is_Read2_Rev_Comp","Is_First_Read","Is_Second_Read",
            "Split","Long_Insert","Short_Insert","Parallel_Read","Everted_Read","Orphan_Read"
        ]
        for feat in featlist:
            if reads_df[feat][read_idx] == 1:
                feat_idx = featlist.index(feat) + 5
                pileup_arr[feat_idx, read_idx, read_inds] = 1
        return pileup_arr


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

def parse_arguments() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i","--input_file",metavar="IN_FILE",required=False,dest="input_file",
        default="../test/testlocs.csv",
        help="5-column csv with colname header SAMPLE,CHROM,REF_BEGIN,REF_END,CLASS, TE"
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
    args = parser.parse_args()
    return args

def main() -> None:
    ua = parse_arguments()
    N = 1
    print("Using ", N, "processes")

    regions = pd.read_csv(ua.input_file, header=None)
    regions.columns = ["Sample","Chrom","Ref_begin","Ref_end","Class","TE"]
    print(len(regions), "files to generate data for.")
    print(regions.head())

    region_generator(regions, ua.output_dir, ua.bam_dir)
    te_specific_region_generator(regions, ua.output_dir, ua.te_bam_dir)

if __name__ == "__main__":
    main()
