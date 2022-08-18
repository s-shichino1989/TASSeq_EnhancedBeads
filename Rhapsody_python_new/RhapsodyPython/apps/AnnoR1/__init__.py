# distutils: extra_compile_args = ["-O3"]
# cython: language_level=3
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import utils
import csv
from AnnoR1._trim_primer import trim_primer
import os
import numpy as np
import copy
import collections
import json
import xopen

class AnnoR1:
    from ._bead_specific_settings import bead_specific_settings
    from ._cell_key_hash_table import cell_key_hash_table
    from ._error_rate_table import read1_error_rate_table
    from ._call_cell_key import cell_key_caller
    from ._misc_functions import misc_functions

    def __init__(self, R1file="", sample="", library="", bead_version="", sample_total_read_count=1):
        # Inputs
        self.R1file = R1file
        self.sample = sample
        self.library = library
        self.bead_version = bead_version
        self.sample_total_read_count = sample_total_read_count

        # Setting thresholds:
        # Baseline LD/bp threshold for rapid filtering of cell key candidates
        self.LDpBP_threshold = 1 / 3

        # In the third pass, if a read has >= 50% chance of being real, call it.
        prob_valid_read = 0.5
        self.prob_valid_read_threshold = np.log(1 - np.exp(np.log(1-prob_valid_read) / np.exp(self.sample_total_read_count)))
        prob_valid_part = 0.5
        self.prob_valid_part_threshold = np.log(1 - np.exp(np.log(1-prob_valid_part) / np.exp(self.sample_total_read_count)))
        self.prob_cell_threshold = 0.5

        # Initialize output file name
        self.intermediate1_file = "{}_Annotation_R1_intermediate1.csv".format(sample)
        if os.path.exists(self.intermediate1_file):
            os.remove(self.intermediate1_file)
        self.intermediate2_file = "{}_Annotation_R1_intermediate2.csv".format(sample)
        if os.path.exists(self.intermediate2_file):
            os.remove(self.intermediate2_file)
        self.output_file = "{}_Annotation_R1.csv.gz".format(sample)
        if os.path.exists(self.output_file):
            os.remove(self.output_file)

        # Annotation Settings
        self.bead = self.bead_specific_settings(bead_version)

        # Miscellaneous functions
        self.misc_functions = self.misc_functions(self)

        # Cell Key Hash Tables
        self.hash_tables = self.cell_key_hash_table(self)

        # Error rate table
        self.error_rate_tables = self.read1_error_rate_table(self)

        # Candidate information passing from first pass to second pass
        self.ins_lnk_from_perfect_hash_table = collections.deque([])
        self.max_num_inslnk_stored = 3000000

        # Candidate information passing from second pass to third pass
        self.candidates_from_hash_table = collections.deque([])
        self.max_num_candidates_stored = 400000

        np.random.seed(seed=7865)

    def call_cell_labels(self):
        # Trim off 5' primer if it is a 5' bead
        if self.bead_version in ['V1P', 'V2_L1P', 'V2_L3P', 'V2_Op3']:
            trimmed_file = "{}_trimmed.fastq.gz".format(self.sample)
            trim_primer(self.R1file, trimmed_file)
            self.R1file = trimmed_file

        # Call perfect cell labels
        self.pass1_call_perfect_cell_labels()

        # Remove low occurrence ins/lnk combinations
        self.misc_functions.filter_ins_lnk_combinations()

        # Call cell labels with 1 LD hash tables
        self.pass2_call_cell_labels_with_hash_tables()

        # Estimate final error rate
        self.error_rate_tables.update_log_prob_tables()

        # Run error-rate-based cell label calling
        self.pass3_call_cell_labels_with_error_correction()

        self.error_rate_tables.output_error_count_tables()
        self.misc_functions.output_stat_counters()

    def pass1_call_perfect_cell_labels(self):
        # Pass-1 Call perfect cell labels
        with open(self.intermediate1_file, 'wt') as out_file, utils.quick_gzip_open(self.R1file) as r1_file:
            r1_fgi = FastqGeneralIterator(r1_file)
            count = 0
            inslnk_stored = 0
            write_buffer = ""
            for r1_record in r1_fgi:
                count += 1
                if count == 1000000:
                    self.misc_functions.filter_ins_lnk_combinations()
                r1_title, r1_seq, r1_qual = r1_record

                # Create cell key caller object
                this_cell_key_caller = self.cell_key_caller(self, r1_seq)

                # Call cell keys
                this_cell_key_caller.call_perfect_cell_keys()

                if this_cell_key_caller.accepted_cell_key:
                    # Cell key called
                    self.misc_functions.update_counts(this_cell_key_caller)
                    self.misc_functions.update_linker_stats(this_cell_key_caller)
                    write_buffer += ",".join([this_cell_key_caller.cell_key,"0-0-0",this_cell_key_caller.umi]) + "\n"
                else:
                    # Read with error, pass it to second pass
                    if inslnk_stored < self.max_num_inslnk_stored:
                        self.ins_lnk_from_perfect_hash_table.append(this_cell_key_caller.potential_ins_lnk_combinations)
                        inslnk_stored += 1
                        write_buffer += "?" + r1_seq[:(self.bead.UMI_start_end_pos[1] + len(self.bead.inserts[-1]) + 4)] + "\n"
                    else:
                        ins_lnk_string = json.dumps(this_cell_key_caller.potential_ins_lnk_combinations)
                        write_buffer += "?" + ins_lnk_string + "?" + r1_seq[:(self.bead.UMI_start_end_pos[1] + len(self.bead.inserts[-1]) + 4)] + "\n"
                if count % 10000 == 0:
                    out_file.write(write_buffer)
                    write_buffer = ""
            out_file.write(write_buffer)
            self.misc_functions.stat_counters['total_num_reads'] = count

    def pass2_call_cell_labels_with_hash_tables(self):
        # Pass-2 Call cell labels with 1-LD hash tables
        self.misc_functions.static_cell_read_count = copy.deepcopy(self.misc_functions.cell_counter)
        num_candidate_stored = 0
        with open(self.intermediate2_file, 'wt') as out_file,  open(self.intermediate1_file, 'rt') as int_file:
            while True:
                lines = int_file.readlines(10000)
                if not lines:
                    break
                write_buffer = ""
                for line in lines:
                    if line[0] != "?":
                        # A perfect cell label, write out directly
                        write_buffer += line
                        continue

                    if len(self.ins_lnk_from_perfect_hash_table) > 0:
                        r1_seq = line[1:].strip("\n")
                        # Create cell key caller object
                        this_cell_key_caller = self.cell_key_caller(self, r1_seq)
                        # Get potential insert and linker identities from first pass
                        this_cell_key_caller.potential_ins_lnk_combinations = self.ins_lnk_from_perfect_hash_table.popleft()
                    else:
                        ins_lnk_string, r1_seq = line[1:].strip("\n").split("?")
                        # Create cell key caller object
                        this_cell_key_caller = self.cell_key_caller(self, r1_seq)
                        # Get potential insert and linker identities from first pass
                        this_cell_key_caller.potential_ins_lnk_combinations = json.loads(ins_lnk_string,
                                                                                         object_hook=self.int_conversion)

                    # Call cell label
                    this_cell_key_caller.call_cell_keys_with_1LD_hash_tables()

                    # Check insert linker quality
                    this_cell_key_caller.get_potential_ins_lnk_combinations()
                    self.misc_functions.update_linker_stats(this_cell_key_caller)

                    if this_cell_key_caller.accepted_cell_key:
                        # Cell key called
                        # Update read count data
                        self.misc_functions.update_counts_second_pass(this_cell_key_caller)
                        # Update error count tables
                        self.error_rate_tables.update_error_count_table(this_cell_key_caller.insert_idx,
                                                                        this_cell_key_caller.linker_idx,
                                                                        this_cell_key_caller.editops)

                        # Write to file
                        error_str = self.misc_functions.editops_to_output(this_cell_key_caller.insert_idx,
                                                                  this_cell_key_caller.editops)
                        write_buffer += ",".join([
                            this_cell_key_caller.cell_key,
                            error_str,
                            this_cell_key_caller.umi
                        ]) + "\n"
                    else:
                        # Read not called
                        # Update read count data
                        self.misc_functions.update_counts_second_pass2(this_cell_key_caller)
                        if this_cell_key_caller.illegible_insert_linker:
                            # If insert/linker cannot be identified, call read as illegible
                            self.misc_functions.update_counts(this_cell_key_caller)
                            write_buffer += ",".join([
                                "x-x-x",
                                "x-x-x",
                                ""
                            ]) + "\n"
                        else:
                            # Pass read to third pass
                            if num_candidate_stored < self.max_num_candidates_stored:
                                self.candidates_from_hash_table.append([this_cell_key_caller.potential_ins_lnk_combinations, this_cell_key_caller.candidates])
                                num_candidate_stored += 1
                                write_buffer += "?" + r1_seq + "\n"
                            else:
                                ins_lnk_string = json.dumps(this_cell_key_caller.potential_ins_lnk_combinations)
                                candidate_string = json.dumps(this_cell_key_caller.candidates)
                                write_buffer += "?" + "?".join([ins_lnk_string, candidate_string, r1_seq]) + "\n"
                out_file.write(write_buffer)

        self.misc_functions.stat_counters["p2_cellular_read"] = self.misc_functions.stat_counters["cellular_read"]

    def pass3_call_cell_labels_with_error_correction(self):
        self.misc_functions.static_cell_read_count = copy.deepcopy(self.misc_functions.cell_counter)
        with xopen.xopen(self.output_file, 'wt', threads=2, compresslevel=1) as out_file, open(
                self.intermediate2_file, 'rt') as int_file:
            while True:
                lines = int_file.readlines(10000)
                if not lines:
                    break
                write_buffer = ""
                for line in lines:
                    if line[0] != "?":
                        # A called cell label, write out directly
                        write_buffer += line
                        continue

                    # Load candidate information from second pass
                    if len(self.candidates_from_hash_table) > 0:
                        # Read the line
                        r1_seq = line[1:].strip("\n")

                        # Create cell key caller object
                        this_cell_key_caller = self.cell_key_caller(self, r1_seq)

                        # Load information from pass 2
                        this_cell_key_caller.potential_ins_lnk_combinations, this_cell_key_caller.candidates = self.candidates_from_hash_table.popleft()
                    else:
                        # Read the line
                        ins_lnk_string, candidate_string, r1_seq = line[1:].strip("\n").split("?")

                        # Create cell key caller object
                        this_cell_key_caller = self.cell_key_caller(self, r1_seq)

                        # Load information from pass 2
                        this_cell_key_caller.potential_ins_lnk_combinations = json.loads(ins_lnk_string,
                                                                                         object_hook=self.int_conversion)
                        this_cell_key_caller.candidates = json.loads(candidate_string, object_hook=self.int_conversion)
                        if "ck1s" not in this_cell_key_caller.candidates:
                            this_cell_key_caller.candidates["ck1s"] = {}
                        if "ck2s" not in this_cell_key_caller.candidates:
                            this_cell_key_caller.candidates["ck2s"] = {}
                        if "ck3s" not in this_cell_key_caller.candidates:
                            this_cell_key_caller.candidates["ck3s"] = {}
                        if "ck1_indel_adjusts" not in this_cell_key_caller.candidates:
                            this_cell_key_caller.candidates["ck1_indel_adjusts"] = {}
                        if "ck2_indel_adjusts" not in this_cell_key_caller.candidates:
                            this_cell_key_caller.candidates["ck2_indel_adjusts"] = {}


                    # Call cell keys
                    this_cell_key_caller.call_cell_keys_with_error_correction()

                    # Update read counts
                    self.misc_functions.update_counts(this_cell_key_caller)

                    if this_cell_key_caller.accepted_cell_key:
                        # Update error count table
                        self.error_rate_tables.update_error_count_table(this_cell_key_caller.insert_idx,
                                                                        this_cell_key_caller.linker_idx,
                                                                        this_cell_key_caller.editops)
                        # Write to file
                        write_buffer += ",".join([
                            this_cell_key_caller.cell_key,
                            self.misc_functions.editops_to_output(this_cell_key_caller.insert_idx,
                                                                  this_cell_key_caller.editops),
                            this_cell_key_caller.umi
                        ]) + "\n"
                    else:
                        # Write illegible read to file
                        write_buffer += ",".join([
                            "x-x-x",
                            "x-x-x",
                            ""
                        ]) + "\n"

                out_file.write(write_buffer)

    def int_conversion(self, obj):
        output = {}
        for k, v in obj.items():
            try:
                k = int(k)
            except ValueError:
                pass
            if isinstance(v, dict):
                v = self.int_conversion(v)
            else:
                if not isinstance(v, list):
                    try:
                        v = int(v)
                    except ValueError:
                        pass
            output[k] = v
        return output
