# distutils: extra_compile_args = ["-O3"]
# cython: language_level=3
import sys
import MistLogger as logging
import argparse
from AnnoR1 import *
import numpy as np
import subprocess
import re
import csv
import os
from Bio import SeqIO
import math

def cli(cli_args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--R1',
                        required=True,
                        help='Read 1 FASTQ file from the split_fastqs node.')
    parser.add_argument('--cutadapt_log',
                        required=True,
                        help='cutadapt log file for retrive total read count.')
    args = parser.parse_args(cli_args)
    sample = os.path.basename(args.R1).split('_R1')[0]
    args_dict = dict(sample=sample, **args.__dict__)
    #logging.info('Running with options: {}'.format(args_dict))
    return args_dict


@utils.node_timer
def annotate_r1(R1, cutadapt_log, sample):
    # read metadata
    checkfastq_lib_name = get_library_name(R1).split("_")
    # this check is needed for the functional tests
    if len(checkfastq_lib_name) == 1:
        library = checkfastq_lib_name[0]
    else:
        library = checkfastq_lib_name[:-1]
        library = "_".join(library)
        
    #automatically detect Rhapsody bead version
    bead_version = get_bead_version(R1)

    #bead_version = run_metadata["Bead_Version"][library]
    #assay = run_metadata["Assay"]

    # Get total read count of this sample from filter metrics
    #sample_filter_metrics = []
    #for filter_metrics_file in filter_metrics.split(","):
    #    if library in filter_metrics_file:
    #        sample_filter_metrics.append(filter_metrics_file)
    #sample_total_read_count = 0
    #for filter_metrics_file in sample_filter_metrics:
    #    with gzip.open(filter_metrics_file,'rt') as f:
    #        first_line = f.readline()
    #        read_count = int(first_line.split(",")[0])
    #        sample_total_read_count += read_count
    
    # Get total read count of this sample from cutadapt filter metrics
    sample_total_read_count = np.loadtxt(cutadapt_log, dtype='int64', delimiter='\t', skiprows=1, usecols=[2])
    sample_total_read_count = np.log(int(sample_total_read_count))
    
    
    logging.info("Running Annotate R1 with bead version: " + bead_version)
    #logging.info("Total number of reads in this sample across all splits is " + str(int(round(np.exp(sample_total_read_count)))))
    AnnoR1_object = AnnoR1(R1, sample, library, bead_version, sample_total_read_count)
    #logging.info("Calling cell labels")
    AnnoR1_object.call_cell_labels()
    #logging.info("Annotated all {} reads of R1".format(AnnoR1_object.misc_functions.stat_counters["total_num_reads"]))
    intermediate1_file = "{}_Annotation_R1_intermediate1.csv".format(sample)
    os.remove(intermediate1_file)
    intermediate2_file = "{}_Annotation_R1_intermediate2.csv".format(sample)
    os.remove(intermediate2_file)
    #os.remove(R1)

def sufficient_poly_t_check(poly_t_seq, minimum_poly_t_seq_len=7):
    # type: (str, int) -> bool
    """ensure capture of polyadenylation site, if read is long enough

    (From the handbook) Following the UMI, a poly(T) tail, the polyadenylation [poly(A)] complement of an mRNA
    molecule, is expected. Each read with a valid cell label is kept for further consideration
    only if >= 6 out of 8 bases after UMI are found to be Ts.

    If the read contains fewer than eight bases in the polyadenylation region, skip this test

    Args:
        poly_t_seq: Read 1 DNA sequence, following the cell label
        minimum_poly_t_seq_len: the minimum length of the poly-T sequence on which to run the poly-t-check

    Returns:
        True for read passes or skips; False for failures

    Usage:
        >>> poly_t_start_position = 60
        >>> s = "GCGCAATCAACTGGCCTGCGACCGACAAGAGGTAGCGGTGACGAAGGGTCAGCGTAATTTTTTTTTTTTTTTTTTT"
        >>> sufficient_poly_t_check(s[poly_t_start_position:])
        True
        >>> s = "AGAACTTCCACTGGCCTGCGACAGAAATCGGGTAGCGGTGACACAGGGAGGGATCAATAATTTTTTTTTTCCATTG"
        >>> sufficient_poly_t_check(s[poly_t_start_position:])
        True
        >>> s = "TAGCTTGTAACTGGCCTTTTTTTTTTTTTTTTTTTTTTTAGAACAGGGTAATTTTTTTTTTGTAAAAACAACTGTG"
        >>> sufficient_poly_t_check(s[poly_t_start_position:])
        False
    """
    if minimum_poly_t_seq_len <= len(poly_t_seq):
        poly_t_check_passing = (6 <= poly_t_seq[:8].count('T'))
    else:
        poly_t_check_passing = True
    return poly_t_check_passing


def get_bead_version(R1):
    key1Dict = dict(enumerate(cell_keys.v2_cell_key1[:96]))
    key1Dict = dict(zip(key1Dict.values(), key1Dict.keys()))

    key2Dict = dict(enumerate(cell_keys.v2_cell_key2[:96]))
    key2Dict = dict(zip(key2Dict.values(), key2Dict.keys()))

    key3Dict = dict(enumerate(cell_keys.v2_cell_key3[:96]))
    key3Dict = dict(zip(key3Dict.values(), key3Dict.keys()))

    # Bead Versions:
    # V1: Original V1 Bead
    # V1P: Original V1 Pretzelball Bead = (V1) or (V1 + 20) for V1 <= V2_Lnk2StartEnd; (V1) or (V1 + 19) for V1 > V1_Lnk2StartEnd; if V1 == V1_Lnk2StartEnd, V1_Lnk2StartEnd[1] -= 1
    # V2_L1: Ligation V2 Bead
    # V2_L1P: Ligation V2 Pretzelball Bead (V2_L1) or (V2_L1 + 20)
    # V2_L3: Ligation V2 Bead with 0-3 fixed insert in front of the oligo (V2_L1 + (0 or 1 or 2 or 3))
    # V2_L3P:  Ligation V2 Pretzelball Bead with 0-3 fixed insert in front of the oligo (V2_L1 + (0 or 1 or 2 or 3)) or (V2_L1 + (0 or 1 or 2 or 3) + 20)

    # V1 cell label locations
    V1_CL1StartEnd = [0, 9]
    V1_Lnk1StartEnd = [9, 21]
    V1_CL2StartEnd = [21, 30]
    V1_Lnk2StartEnd = [30, 43]
    V1_CL3StartEnd = [43, 52]
    V1_UMIStartEnd = [52, 60]
    V1_PolyTStartEnd = [60, 68]

    # V2 Bead label locations
    V2_CL1StartEnd = [0, 9]
    V2_Lnk1StartEnd = [9, 13]
    V2_CL2StartEnd = [13, 22]
    V2_Lnk2StartEnd = [22, 26]
    V2_CL3StartEnd = [26, 35]
    V2_UMIStartEnd = [35, 43]
    V2_PolyTStartEnd = [43, 51]

    # Bead Version Recognition by Cell Labels 2 and 3

    # Counting sequences that belong to each version of the beads. (One sequence can belong to several versions.)
    BeadVersion_Census = {"V1": 1, "V1P": 1, "V2_L1": 1, "V2_L1P": 1, "V2_L3": 1, "V2_L3P": 1,
                          "V2_Op3": 1}  # Start with 1 so we never divide by zero

    with utils.quick_gzip_open(R1) as r1_file:
        records = SeqIO.parse(r1_file, 'fastq')
        countReads = 0  # Count number of valid reads
        for inx, record in enumerate(records, start=1):
            countReads += 1
            if inx == 100000:
                break

    num_skip = 100
    if countReads < 100000:
        num_skip = math.floor(countReads/1000)

    with utils.quick_gzip_open(R1) as r1_file:
        records = SeqIO.parse(r1_file, 'fastq')
        countValid = 0  # Count number of valid reads
        for inx, record in enumerate(records, start=1):
            if num_skip > 1 and inx % num_skip != 0:
                continue
            if countValid > 1000:
                break

            r1_seq = str(record.seq)
            # V1
            if (getSequenceByStartEnd(r1_seq, V1_CL2StartEnd) in key2Dict and
                    getSequenceByStartEnd(r1_seq, V1_CL3StartEnd) in key3Dict):
                BeadVersion_Census["V1"] += 1
                countValid += 1

            # V1P
            if (getSequenceByStartEnd(r1_seq, [x + 20 for x in V1_CL2StartEnd]) in key2Dict and
                    getSequenceByStartEnd(r1_seq, [x + 19 for x in V1_CL3StartEnd]) in key3Dict):
                BeadVersion_Census["V1P"] += 1
                countValid += 1

            # V2_L1
            if (getSequenceByStartEnd(r1_seq, V2_CL2StartEnd) in key2Dict and
                    getSequenceByStartEnd(r1_seq, V2_CL3StartEnd) in key3Dict):
                BeadVersion_Census["V2_L1"] += 1
                BeadVersion_Census["V2_L3"] += 1
                countValid += 1

            # V2_L1P
            if (getSequenceByStartEnd(r1_seq, [x + 20 for x in V2_CL2StartEnd]) in key2Dict and
                    getSequenceByStartEnd(r1_seq, [x + 20 for x in V2_CL3StartEnd]) in key3Dict):
                BeadVersion_Census["V2_L1P"] += 1
                BeadVersion_Census["V2_L3P"] += 1
                BeadVersion_Census["V2_Op3"] += 1
                countValid += 1

            # V2_L3
            if (
                    (getSequenceByStartEnd(r1_seq, [x + 1 for x in V2_CL2StartEnd]) in key2Dict and
                     getSequenceByStartEnd(r1_seq, [x + 1 for x in V2_CL3StartEnd]) in key3Dict)
                    or
                    (getSequenceByStartEnd(r1_seq, [x + 2 for x in V2_CL2StartEnd]) in key2Dict and
                     getSequenceByStartEnd(r1_seq, [x + 2 for x in V2_CL3StartEnd]) in key3Dict)
                    or
                    (getSequenceByStartEnd(r1_seq, [x + 3 for x in V2_CL2StartEnd]) in key2Dict and
                     getSequenceByStartEnd(r1_seq, [x + 3 for x in V2_CL3StartEnd]) in key3Dict)
            ):
                BeadVersion_Census["V2_L3"] += 1
                countValid += 1

            # V2_L3P
            if (
                    (getSequenceByStartEnd(r1_seq, [x + 21 for x in V2_CL2StartEnd]) in key2Dict and
                     getSequenceByStartEnd(r1_seq, [x + 21 for x in V2_CL3StartEnd]) in key3Dict)
                    or
                    (getSequenceByStartEnd(r1_seq, [x + 22 for x in V2_CL2StartEnd]) in key2Dict and
                     getSequenceByStartEnd(r1_seq, [x + 22 for x in V2_CL3StartEnd]) in key3Dict)
                    or
                    (getSequenceByStartEnd(r1_seq, [x + 23 for x in V2_CL2StartEnd]) in key2Dict and
                     getSequenceByStartEnd(r1_seq, [x + 23 for x in V2_CL3StartEnd]) in key3Dict)
            ):
                BeadVersion_Census["V2_L3P"] += 1
                BeadVersion_Census["V2_Op3"] += 1
                countValid += 1

            # V2_Op3
            if (
                    (getSequenceByStartEnd(r1_seq, [x + 24 for x in V2_CL2StartEnd]) in key2Dict and
                     getSequenceByStartEnd(r1_seq, [x + 24 for x in V2_CL3StartEnd]) in key3Dict)
                    or
                    (getSequenceByStartEnd(r1_seq, [x + 25 for x in V2_CL2StartEnd]) in key2Dict and
                     getSequenceByStartEnd(r1_seq, [x + 25 for x in V2_CL3StartEnd]) in key3Dict)
                    or
                    (getSequenceByStartEnd(r1_seq, [x + 26 for x in V2_CL2StartEnd]) in key2Dict and
                     getSequenceByStartEnd(r1_seq, [x + 26 for x in V2_CL3StartEnd]) in key3Dict)
                    or
                    (getSequenceByStartEnd(r1_seq, [x + 27 for x in V2_CL2StartEnd]) in key2Dict and
                     getSequenceByStartEnd(r1_seq, [x + 27 for x in V2_CL3StartEnd]) in key3Dict)
                    or
                    (getSequenceByStartEnd(r1_seq, [x + 28 for x in V2_CL2StartEnd]) in key2Dict and
                     getSequenceByStartEnd(r1_seq, [x + 28 for x in V2_CL3StartEnd]) in key3Dict)
            ):
                BeadVersion_Census["V2_Op3"] += 1
                countValid += 1

    # Check if file is empty
    bead_version = "V1"
    if countValid > 1000:
        # Identify bead version based on counts
        # Check V1 or V2
        sumV1 = BeadVersion_Census["V1"] + BeadVersion_Census["V1P"]
        sumV2 = BeadVersion_Census["V2_L3"] + BeadVersion_Census[
            "V2_Op3"]  # Sum of these two includes all inserts from 0-3 and 20-28
        if sumV1 > sumV2:
            # V1 Bead
            if BeadVersion_Census["V1"] > BeadVersion_Census["V1P"]:
                bead_version = "V1"
            else:
                bead_version = "V1P"
        else:
            # V2 Bead
            # Is it a 5 prime bead
            if BeadVersion_Census["V2_Op3"] > BeadVersion_Census["V2_L3"]:
                # It is a 5 prime bead
                # Does it have the PhiX Option 3 inserts (extra 0, 3, 5 bases that causes total insert length to range from 24-28)?
                if BeadVersion_Census["V2_L3P"] / BeadVersion_Census["V2_Op3"] < 0.8:
                    # This is the V2_Op3 bead, as > ~20% of the reads have the opt 3 inserts.
                    bead_version = "V2_Op3"
                else:
                    # It does not have the PhiX Option 3 inserts.
                    # Does it have the V2 inserts?
                    if BeadVersion_Census["V2_L1P"] / BeadVersion_Census["V2_L3P"] < 0.8:
                        # This is the V2_L3P bead, > ~20% of the reads have the V2 inserts.
                        bead_version = "V2_L3P"
                    else:
                        bead_version = "V2_L1P"

            else:
                # Not a 5 prime bead
                if BeadVersion_Census["V2_L1"] / BeadVersion_Census["V2_L3"] < 0.8:
                    # This is the V2_L3 bead, as ~ more than 20% of the reads have the V2 inserts.
                    bead_version = "V2_L3"
                else:
                    bead_version = "V2_L1"

    elif countValid > 0:
        # Number of valid read count is low (number of input sequence is very low). Just use the type with the highest
        bead_version = max(BeadVersion_Census, key=BeadVersion_Census.get)
    else:
        raise Exception('Zero read pairs with valid cell labels in this library. Please check FASTQ files and make sure they are from a Rhapsody library.  Pipeline cannot proceed.')

    return bead_version


def get_library_name(filename):
    """
    Tries to get the library_name from the filename.

    Args:
        filename: a path to a fastq file
    Returns:
        library_name: the library name for the fastq file
    """
    base_name = os.path.basename(filename)

    # remove any leading periods that would
    # result in the files being hidden
    base_name = base_name.lstrip(".")

    # the library name will be the entire filename up to the illumina
    # basespace elements if they exist or the entire filename
    # up to the first period (usually for the file extension)
    # if the illumina basespace elements do not exist
    # LIB_NAME_REGEX = '^(.*?)(_S[0-9]*)?(_L[0-9]*)?(_R[1|2].*)?\.(.*?)\.(.*)$'
    LIB_NAME_REGEX = re.compile(r'^(.*?)'
                            r'(_S[0-9]*)?'
                            r'(_L[0-9]*)?'
                            r'(_R[1|2].*)?'
                            r'\.(.*?)'
                            r'\.(.*)$')
    reg_matches = re.findall(LIB_NAME_REGEX, base_name)

    # if the regex fails, use a default name
    if len(reg_matches) == 0:
        library_name = "SampleName"
    elif len(reg_matches[0]) == 0:
        library_name = "SampleName"
    # if the first group is an empty string, use a default name
    elif reg_matches[0][0] == "":
        library_name = "SampleName"
    else:
        library_name = reg_matches[0][0]

    # if the resulting library name doesn't have any
    # letters or numbers, use the default name
    ALPHA_NUMERIC_REGEX = re.compile(r'[a-zA-Z0-9]+')
    if re.search(ALPHA_NUMERIC_REGEX, library_name) is None:
        library_name = "SampleName"

    # if the library name is empty or just part of the
    # file extension, then return a default sample name
    if library_name.lower() in ["", "fastq", "fq", "gz"]:
        library_name = "SampleName"

    if library_name == "SampleName":
        # the library_name does not have to be
        # unique but it cannot be empty
        logging.warning("Could not determine a sample name "
                        "for the file '{}'. Using the default: "
                        "'SampleName'.".format(filename))

    return library_name


def getSequenceByStartEnd(Sequence,StartEnd):
    return Sequence[StartEnd[0]:StartEnd[1]]

def package_main():
    """entry point for pipeline"""
    main()


def main(cli_args=None):
    """entry point for testing"""
    return annotate_r1(**cli(cli_args))


if __name__ == '__main__':
    package_main()
