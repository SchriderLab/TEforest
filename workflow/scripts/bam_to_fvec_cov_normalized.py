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

# from exceptions import *
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
            output_dir (str/pathlike): Directory to write output npy files to. Intermediate directory with sample (i.e. 'A4') will be created in this.
            bamfile (str/pathlike): Bamfile to access for alignment
        """
        self.output_dir = output_dir
        self.bamfile = bamfile

    def _get_bam_data(
        self, bamfile, chromosome: str, start: int, end: int, threads: int = 1
    ) -> Union[df, None]:
        """Gets data from BAM file, parses features that require it, fills in
        feature array, and then calculates coverage on a per-base step

        Args:
            chromosome (str): Chr identifier in BAM file
            start (int): Start coordinate for region
            end (int): End coordinate for region
            threads (int, optional): Threads to give to samtools. Defaults to 1.

        Returns:
            df: Dataframe containing the entirety of the SAM entry with added columns.
        """
        samtools_cmd = [
            "samtools",
            "view",
            "-@ " + str(threads),
            bamfile,
            chromosome,
        ]
        if (start and end) is not None:
            samtools_cmd[-1] = chromosome + ":" + str(int(start)) + "-" + str(int(end))

        # Get number of rows in bam file slice
        num_reads = int(subprocess.check_output(samtools_cmd + ["-c"]).decode("utf-8"))
        # if num_reads > 10000:
        #    logging.warning(
        #        self.reg_str + "has too many reads, subsetting to a random 10000."
        #    )
        #    samtools_cmd += "| shuf -n 10000"

        if num_reads == 0:
            # print("!" + samtools_cmd[-1] + " -- no reads in bam file.")
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

        # Add explicit bitflag cols
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

        # print(bam_df.head())

        return bam_df

    def _filter_dataframe(self, region_data: df) -> df:
        """
        Basic preprocessing to remove 0 quality reads and null CIGAR strings.

        Args:
            region_data (df): All read data from a given region.

        Returns:
            df: Cleaned up dataframe containing read data.
        """
        region_data = region_data.drop_duplicates(ignore_index=True)
        region_data = region_data.dropna()
        # removed this bc we need low mapping quality reads for repetitve regions
        # region_data = region_data[region_data["Quality"] > 0]
        region_data = region_data[region_data["CIGAR"] != "*"].reset_index()

        return region_data

    def _classify_reads(self, bam_df: df) -> df:
        """Creates a dict with flags for if a read falls into a category of interest.

        Args:
            read_slice (pd.Series): DataFrame containing all feature information extracted from BAM file

        Returns:
            dict: Dictionary with aspects of reads of interest with 0 or 1 values.
        """
        try:
            bam_df["Read_Length"] = bam_df["Sequence"].apply(len)
            read_len = stats.mode(np.abs(bam_df["Read_Length"]), keepdims=True).mode[0]
        except:
            return None

        # Probably a better way to define split-reads?
        bam_df["Split"] = np.where(
            np.abs(bam_df["Read_Length"]) > read_len + 5,
            1,
            0,
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
                & (
                    (bam_df["Is_Read1_Rev_Comp"] == 0)
                    & (bam_df["Is_Read2_Rev_Comp"] == 1)
                )
            )
            | (
                (bam_df["Is_First_Read"] == 1)
                & ((bam_df["Rnext"] == "=") & (bam_df["Read_Start"] < bam_df["Pnext"]))
                & (
                    (bam_df["Is_Read1_Rev_Comp"] == 0)
                    & (bam_df["Is_Read2_Rev_Comp"] == 1)
                )
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
        """Turns list of string bitflag values and turns them
            into a numpy array of expanded binary flags

        Args:
            bitflag_list (list[str]): List of bitflags in string format
            num_reads (int): Number of reads in the region

        Returns:
            bitflag_array: ndarray where each sub-array is a 1D array of binary flags
        """
        bitflag_array = np.zeros((num_reads, 8))
        for i in range(len(bitflag_list)):
            bitflag_str = format(int(bitflag_list[i]), "#016b")[::-1][0:8]
            bitflag_array[i] = list(bitflag_str)

        return bitflag_array

    def _summarize_reads(self, reads_arr: arr) -> arr:
        return np.sum(reads_arr, axis=1)

    def _create_sum_stats(self, summ_arr: arr, read_counts: arr) -> arr:
        """Creates array of summary statistics for region
        For each feature:
            Mean
            Median
            StdDev
            IQ Rangeh
        """
        feat_vec = np.zeros(summ_arr.shape[0] * 4)
        # Replace zeros in read_counts to avoid division by zero
        #in these cases all features should be 0 
        read_counts_safe = np.where(read_counts == 0, 1, read_counts)

        for feat_ind in range(0, summ_arr.shape[0]):
            _stand_feat = summ_arr[feat_ind, :]
            if len(_stand_feat) == 0:
                return None

            # Normalize by the number of reads at each base
            normalized_feat = _stand_feat / read_counts_safe
            normalized_feat = np.nan_to_num(normalized_feat)  # Handle division by zero or NaNs

            vec_ind = feat_ind * 4
            feat_vec[vec_ind + 0] = np.mean(normalized_feat)
            feat_vec[vec_ind + 1] = np.median(normalized_feat)
            feat_vec[vec_ind + 2] = np.std(normalized_feat)
            feat_vec[vec_ind + 3] = stats.iqr(normalized_feat, interpolation="midpoint")

        return feat_vec

    def _set_outfile_name(self, region: pd.Series) -> str:
        """Creates name of outfile given region details
        Args:
            region (pd tuple): One row of the regions dataframe
        Returns:
            outfile: Name of the new outfile
        """
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
        """Generates quick string representation of region for exception handling and info."""
        self.reg_str = f"-".join(
            [
                region["Chrom"],
                str(region["Ref_begin"]),
                str(region["Ref_end"]),
            ]
        )

    def _check_create_dir(self, region: pd.Series):
        """Checks for directory existing for file to save, creates if not."""
        if not os.path.exists(os.path.join(self.output_dir, region["Sample"])):
            print(os.path.join(self.output_dir, region["Sample"]))
            os.makedirs(os.path.join(self.output_dir, region["Sample"]))

    def _save_featvec(self, outfile_name: str, region: pd.Series, feat_vec: arr):
        """Saves feature vector to npy file."""
        np.save(
            os.path.join(self.output_dir, region["Sample"], outfile_name + ".npy"),
            feat_vec,
        )
    def _save_readcount(self, outfile_name: str, region: pd.Series, feat_vec: arr):
        """Saves readcount to npy file."""
        np.save(
            os.path.join(self.output_dir, region["Sample"], outfile_name + "_readcount.npy"),
            feat_vec,
        )

    def create_region(self, region: pd.Series, bamfile: str, savefile: bool = True):
        """Generates data for a given region according to inputs, for writing files.
        Args:
            region (pd.Df): Information about genomic region
            self.output_dir (str): Directory to output npy files to
            num_procs (int): Number of processes to give samtools for accessing BAM files
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
            reg_end = int(region["Ref_end"])
            start_diff = int(bam_df["Read_Start"].min())

            mapper = ReadMapper(start_diff, reg_start, reg_end)
            pileup_arr, read_counts_per_base = mapper.map_reads(bam_df)
            #print("Read counts per base:", read_counts_per_base)
            if pileup_arr is None:
                raise EmtpyFeatVecError

            summary_arr = self._summarize_reads(pileup_arr)
            # print("Summary arr", summary_arr)
            feat_vec = self._create_sum_stats(summary_arr, read_counts_per_base)
            # print("Feat vec", feat_vec)
            
            if savefile:
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
        """Generates data for a given region according to inputs, for writing files.
        Instead of saving to a file, this will load in the preexisting file, and
        concatenate the feature vectors together.
        Args:
            region (pd.Df): Information about genomic region
            self.output_dir (str): Directory to output npy files to
            num_procs (int): Number of processes to give samtools for accessing BAM files
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
            reg_end = int(region["Ref_end"])
            start_diff = int(bam_df["Read_Start"].min())

            mapper = ReadMapper(start_diff, reg_start, reg_end)
            pileup_arr, read_counts_per_base = mapper.map_reads(bam_df)
            if pileup_arr is None:
                raise EmtpyFeatVecError

            summary_arr = self._summarize_reads(pileup_arr)
            # print("Summary arr", summary_arr)
            outfile = self._set_outfile_name(region)
            read_counts_per_base = np.load(
                    os.path.join(self.output_dir, region["Sample"], outfile + "_readcount.npy")
                )
            feat_vec = self._create_sum_stats(summary_arr, read_counts_per_base)
            # print("Feat vec", feat_vec)

            if savefile:
                outfile = self._set_outfile_name(region)
                first_feat_vec = np.load(
                    os.path.join(self.output_dir, region["Sample"], outfile + ".npy")
                )
                second_feat_vec = np.concatenate((first_feat_vec, feat_vec), axis=0)
                # if not os.path.exists(
                #    os.path.join(self.output_dir, region["Sample"], outfile)
                # ):
                #    self._check_create_dir(region)
                self._save_featvec(outfile, region, second_feat_vec)
                #    print(outfile + " made.")

                # else:
                #    print(outfile + " already exists. Passing...")
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
        Maps reads onto multidimensional array, segregates this process from the overall feature vector generation.

        Args:
            start_diff (int): Number of bases to subtract from all start regions to get relative 0
            start (int): Actual start coordiante
            end (int): Actual end coordinate
        """
        self.start_diff = start_diff
        self.start = start
        self.end = end

    def map_reads(self, reads_df: df) -> Union[None, arr, arr]:
        """
        Map reads onto 3D matrix of size (features, reads, bases).

        Args:
            reads_df (df): Annotated reads data from BAM and preprocessing


        Returns:
            Union[None, arr]: 3D matrix of sparse read feature maps
        """
        # fmt: off

        reads_df.loc[: , "Read_Start"] = reads_df.loc[: , "Read_Start"] - self.start_diff

        reads_df = reads_df.reset_index()
        if reads_df.shape[0] > 0:
            # Feats, reads, bp
            pileup_arr = np.zeros((22 , reads_df.shape[0] , self.end - self.start))
            max_chunk_index = pileup_arr.shape[2] - 1

            for i in range(pileup_arr.shape[1]):
                # Draw CIGAR string values on row
                read_start_ind = reads_df["Read_Start"][i]
                binary_cig_arr = self._parse_cigar(reads_df["CIGAR"][i])
                _read_inds = list(range(read_start_ind , read_start_ind + binary_cig_arr.shape[2]))
                read_inds =[i for i in _read_inds if i >= 0]

                # Nothing to map
                if read_inds[0] >= max_chunk_index:
                    continue

                read_inds , pileup_arr = self._map_cigar(pileup_arr , i , read_inds , max_chunk_index , binary_cig_arr)
                pileup_arr = self._iterate_feats(reads_df , pileup_arr , i , read_inds)
                #pileup_arr = self._draw_inserts(reads_df, pileup_arr, i, read_inds, max_chunk_index)

                # Quality scores, will just scale these post-generation instead of binary vectors
                pileup_arr[19 , i , read_inds] = reads_df["Quality"][i]

                # Template length as a separate one?
                pileup_arr[20 , i , read_inds] = np.abs(reads_df["Template_Length"][i])
                # Add new feature: read count (each read at each base gets a value of '1')
                pileup_arr[21, i, read_inds] = 1
            # Count the number of reads mapping to each base
            read_counts_per_base = np.sum(pileup_arr[21, :, :], axis=0)
            #remove read counts from features, we don't want it
            pileup_arr = pileup_arr[:-1, :, :]
            return pileup_arr, read_counts_per_base

        else:
            return None

        # fmt: on

    def _map_cigar(
        self, pileup_arr, read_idx, read_inds, max_chunk_index, binary_cig_arr
    ) -> Tuple[List[int], arr]:
        """
        Handle reads laying partially on region by only slicing read to max index.
        Unlike templates these can't lay over the start since the first index is at the first read.
        Also returns read_inds, in case they needed trimmed, to be used for downstream processes.

        Args:
            pileup_arr (arr): Feature by bases map of genomic region.
            read_idx (int): Base-pair index for accessing pileup_arr
            read_inds (arr): 1D array of BP indices where read is laying on genomic region.
            max_chunk_idndex (int): Maximum region cutoff, everything trimmed to this location.
            binary_cig_arr (arr): CIGAR string mapped onto a OHE CIGAR by BP matrix.

        Returns:
            arr: Genomic feature map with CIGAR feature mapped on it.
        """

        if read_inds[-1] >= max_chunk_index:
            # Shorten indices to match
            read_inds = read_inds[0 : read_inds.index(max_chunk_index)]
            pileup_arr[0:5, read_idx : read_idx + 1, read_inds] = binary_cig_arr[
                :, :, 0 : len(read_inds)
            ]

        else:
            pileup_arr[0:5, read_idx : read_idx + 1, read_inds] = binary_cig_arr

        return read_inds, pileup_arr

    def _parse_cigar(self, cig_string: str) -> arr:
        """Generates a list of greyscale values that represent different CIGAR string entries per-base.

        Args:
            cig_string (str): CIGAR String, in form of "54M25I" etc.

        Returns:
            binary_cig_arr: Numpy array with 1 dim per CIGAR entry, 1 if present at location 0 if not.
        """
        # Split CIGAR into list separated by numbers, get matching sites
        # as the start/end depending on the soft clipping status
        cig_list = re.split("([a-zA-Z])", cig_string)[
            0:-1
        ]  # Create list for easier parsing of numbers

        cig_expanded = ""
        _inds = list(range(len(cig_list)))[::2]
        for i in _inds:
            cig_expanded += int(cig_list[i]) * cig_list[i + 1]

        # Make list of binary values for each opt
        binary_cig_arr = np.zeros((5, 1, len(cig_expanded)))
        cig_opts = ["M", "D", "I", "S", "H"]
        for i in range(len(cig_opts)):
            for j in range(len(cig_expanded)):
                if cig_expanded[j] == cig_opts[i]:
                    binary_cig_arr[i, 0:1, j] = 1

        return binary_cig_arr

    def _iterate_feats(self, reads_df, pileup_arr, read_idx, read_inds) -> arr:
        """
        Iterates through bitflag in SAM entry and maps onto feature vectors.

        Args:
            reads_df (df): All read information for given region.
            pileup_arr (arr): Feature by bases map of genomic region.
            read_idx (int): Base-pair index for accessing pileup_arr
            read_inds (arr): 1D array of BP indices where read is laying on genomic region.

        Returns:
            arr: Genomic feature map with CIGAR feature mapped on it.
        """
        featlist = [
            "Paired",
            "Proper_Pair",
            "Is_Read1_Unmapped",
            "Is_Read2_Unmapped",
            "Is_Read1_Rev_Comp",
            "Is_Read2_Rev_Comp",
            "Is_First_Read",
            "Is_Second_Read",
            "Split",
            "Long_Insert",
            "Short_Insert",
            "Parallel_Read",
            "Everted_Read",
            "Orphan_Read",
        ]
        for feat in featlist:
            if reads_df[feat][read_idx] == 1:
                pileup_arr[featlist.index(feat) + 5, read_idx, read_inds] = 1

        return pileup_arr

    def _draw_inserts(
        self, reads_df, pileup_arr, read_idx, read_inds, max_chunk_index
    ) -> arr:
        """
        Draw insert between reads on row.
        Similar check as above, where values are cut off at end/beginning of region if they run over.

        Vals are 0.5 because you get one half of template from each read direction.
        Perfect paired map will result in value of 1 for that bp.

        Args:
            reads_df (df): All read information for given region.
            pileup_arr (arr): Feature by bases map of genomic region.
            read_idx (int): Base-pair index for accessing pileup_arr
            read_inds (arr): 1D array of BP indices where read is laying on genomic region.
            max_chunk_idndex (int): Maximum region cutoff, everything trimmed to this location.

        Returns:
            arr: Genomic feature map with CIGAR feature mapped on it.
        """

        if reads_df["Template_Length"][read_idx] >= 0:
            template_start = read_inds[-1] + 1
            template_end = read_inds[-1] + 1 + reads_df["Template_Length"][read_idx]

            if template_end > max_chunk_index:
                pileup_arr[19, read_idx, template_start:max_chunk_index] = 0.5
            else:
                pileup_arr[19, read_idx, template_start:max_chunk_index] = 0.5

        else:
            template_start = read_inds[0] - np.abs(
                reads_df["Template_Length"][read_idx]
            )
            template_end = read_inds[0]

            if template_start <= 0:
                pileup_arr[19, read_idx, 0:template_end] = 0.5
            else:
                pileup_arr[19, read_idx, template_start:template_end] = 0.5

        return pileup_arr


def region_generator(sub_regions, output_dir, bam_dir):
    """
    Run this if running this script separately on a dataset.
    Otherwise create an SVMaker object and call create_mapper() for a single region.
    """

    for idx, region in sub_regions.iterrows():
        # print(region)
        bamfile = os.path.join(bam_dir, region["Sample"] + ".bam")
        # this was used to access mcclintock bams, which add a _1.
        # bamfile = os.path.join(bam_dir, region["Sample"] + "_1.sorted.bam")
        # this line was used to get the final heterozygote file
        # in the future the appropriate sample name should be used
        # bamfile = os.path.join(bam_dir, "final_product.sort.bam")
        svm = SVMaker(bamfile, output_dir)

        svm.create_region(region, bamfile, True)

        if idx % 1000 == 0:
            gc.collect()


def te_specific_region_generator(sub_regions, output_dir, te_bam_dir):
    """
    Run this if running this script separately on a dataset.
    Otherwise create an SVMaker object and call create_mapper() for a single region.
    """

    for idx, region in sub_regions.iterrows():
        # print(region)
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
        "-i",
        "--input_file",
        metavar="IN_FILE",
        required=False,
        dest="input_file",
        help='5-column csv with colname header SAMPLE,CHROM,REF_BEGIN,REF_END,CLASS \
                                where SAMPLE is the name of the BAM file being accessed (without the ".bam"), \
                                can have multiple SAMPLES in the same input file, but have all \
                                bam files in the same directory. CHROM is the name of the chromosome \
                                exactly as it is accessed through the BAM file. REF_BEGIN and REF_END \
                                are fairly self-explanatory, but should be accessible in the BAM or will \
                                be ignored. Sample file can be found in the /examples/ directory.',
        default="../test/testlocs.csv",
    )

    parser.add_argument(
        "-od",
        "--output_dir",
        metavar="OUTPUT_DIR",
        required=False,
        default=".",
        dest="output_dir",
        help="Directory to save npy files to. Defaults to cwd.",
    )

    parser.add_argument(
        "-bd",
        "--bam_dir",
        metavar="BAM_DIR",
        type=str,
        required=False,
        default="TrainingData/bamfiles/",
        dest="bam_dir",
        help="Path to BAM file directory. BAM filenames must match the entries in SAMPLE \
                            column in input-file.csv. Defaults to /data/.",
    )

    parser.add_argument(
        "-tv",
        "--train_val",
        metavar="TRAIN_VAL",
        type=str,
        required=False,
        default="train",
        dest="train_val",
        help="Whether to store data in the /train/ or /val/ yolo subdirectories.",
    )

    parser.add_argument(
        "--stop_idx",
        metavar="STOP_INDEX",
        type=str,
        required=False,
        dest="stop_idx",
        help="Where to stop indexing the pandas dataframe for SLURM array parallelization.",
    )

    parser.add_argument(
        "-tebd",
        "--te_bam_dir",
        metavar="TE_BAM_DIR",
        type=str,
        required=False,
        default="TrainingData/bamfiles/",
        dest="te_bam_dir",
        help="Path to TE BAM file directory. BAM filenames must match the entries in SAMPLE \
                            column in input-file.csv. TE specific bam files must follow \
                            structure TEname_to_ISO1.bam",
    )
    args = parser.parse_args()

    return args


def main() -> None:
    ua = parse_arguments()

    N = mp.cpu_count() - 1 or 1
    print("Using ", N, "processes")

    regions = pd.read_csv(ua.input_file, header=None)  # , sep="\t")
    regions.columns = ["Sample", "Chrom", "Ref_begin", "Ref_end", "Class", "TE"]
    print(len(regions), "files to generate data for.")
    print(regions.head())

    # n = 5000
    # region_chunks = [regions[i : i + n] for i in range(0, regions.shape[0], n)]
    region_generator(regions, ua.output_dir, ua.bam_dir)
    te_specific_region_generator(regions, ua.output_dir, ua.te_bam_dir)

    # p = mp.Pool(processes=1)
    # p.starmap(
    #    region_generator,
    #    [
    #        (region, od, bd)
    #        for region in region_chunks
    #        for od in [ua.output_dir]
    #        for bd in [ua.bam_dir]
    #    ],
    # )
    #
    # p.close()
    # p.join()


if __name__ == "__main__":
    main()
