# distutils: extra_compile_args = ["-O3"]
# cython: language_level=3
from AnnoR1 import cell_keys
import copy
#import Levenshtein
import MistLogger as logging
import numpy as np
import math
from rapidfuzz.distance import Levenshtein as lv

class cell_key_hash_table:
    def __init__(self, AnnoR1):
        # Class for obtaining Levenshtein editops between read 1 and truth
        # Part 1's output format: [cell key index, Levenshtein editops, [indel adjustment for cell key 2, insert index, linker set index]]
        # Part 2's output format: [cell key index, Levenshtein editops, [indel adjustment for cell key 3, linker set index]]
        # Part 3's output format: [cell key index, Levenshtein editops, [indel adjustment for UMI, linker set index]]
        self.AnnoR1 = AnnoR1

        # Hash tables for cell keys, without other sections of read1. (Key: cell key. Value: cell key index)
        key1Dict = dict(enumerate(cell_keys.v2_cell_key1, 1))
        self.cell_keys_1 = dict(zip(key1Dict.values(), key1Dict.keys()))
        key2Dict = dict(enumerate(cell_keys.v2_cell_key2, 1))
        self.cell_keys_2 = dict(zip(key2Dict.values(), key2Dict.keys()))
        key3Dict = dict(enumerate(cell_keys.v2_cell_key3, 1))
        self.cell_keys_3 = dict(zip(key3Dict.values(), key3Dict.keys()))

        # Lists for storing part locations, overhang lengths and LD thresholds for fast collection of potential candidates
        self.fast_check_ins_lnk_LD_threshold = []
        self.fast_collect_part1_data = []
        self.fast_collect_part2_data = []
        self.fast_collect_part3_data = []
        self.__initialize_fast_check_ins_lnk_LD_threshold()
        self.__initialize_fast_part1_collection_data()
        self.__initialize_fast_part2_collection_data()
        self.__initialize_fast_part3_collection_data()

        # Lists for storing data for calculating Levenshtein editops of candidate cell keys
        self.part1_LD_cal_data = []
        self.part2_LD_cal_data = []
        self.part3_LD_cal_data = []
        self.__initialize_part1_locations_and_sequences_for_editops()
        self.__initialize_part2_locations_and_sequences_for_editops()
        self.__initialize_part3_locations_and_sequences_for_editops()

        # Each part's positions on the read
        self.part1_ck_start_ends = np.zeros((len(self.AnnoR1.bead.inserts), 2)).astype(int)
        self.part1_1del_ck_start_ends = np.zeros((len(self.AnnoR1.bead.inserts), 2)).astype(int)
        self.part1_1ins_ck_start_ends = np.zeros((len(self.AnnoR1.bead.inserts), 2)).astype(int)

        self.part2_ck_start_ends = np.zeros((len(self.AnnoR1.bead.inserts), 2)).astype(int)
        self.part2_1del_ck_start_ends = np.zeros((len(self.AnnoR1.bead.inserts), 2)).astype(int)
        self.part2_1ins_ck_start_ends = np.zeros((len(self.AnnoR1.bead.inserts), 2)).astype(int)

        self.part3_ck_start_ends = np.zeros((len(self.AnnoR1.bead.inserts), 2)).astype(int)
        self.part3_1del_ck_start_ends = np.zeros((len(self.AnnoR1.bead.inserts), 2)).astype(int)
        self.part3_1ins_ck_start_ends = np.zeros((len(self.AnnoR1.bead.inserts), 2)).astype(int)

        # Hash tables for cell label calling. Keys are sequences, values are cell label indices
        self.perfect_cell_keys_1 = [ [{} for lnk in range(len(self.AnnoR1.bead.linker_sets))] for ins in range(len(self.AnnoR1.bead.inserts)) ]
        self.perfect_cell_keys_2 = [{} for lnk in range(len(self.AnnoR1.bead.linker_sets))]
        self.perfect_cell_keys_3 = [{} for lnk in range(len(self.AnnoR1.bead.linker_sets))]

        self.one_sub_cell_keys_1 = [ [{} for lnk in range(len(self.AnnoR1.bead.linker_sets))] for ins in range(len(self.AnnoR1.bead.inserts)) ]
        self.one_sub_cell_keys_2 = [{} for lnk in range(len(self.AnnoR1.bead.linker_sets))]
        self.one_sub_cell_keys_3 = [{} for lnk in range(len(self.AnnoR1.bead.linker_sets))]

        self.one_del_cell_keys_1 = [ [{} for lnk in range(len(self.AnnoR1.bead.linker_sets))] for ins in range(len(self.AnnoR1.bead.inserts)) ]
        self.one_del_cell_keys_2 = [{} for lnk in range(len(self.AnnoR1.bead.linker_sets))]
        self.one_del_cell_keys_3 = [{} for lnk in range(len(self.AnnoR1.bead.linker_sets))]

        self.one_ins_cell_keys_1 = [ [{} for lnk in range(len(self.AnnoR1.bead.linker_sets))] for ins in range(len(self.AnnoR1.bead.inserts)) ]
        self.one_ins_cell_keys_2 = [{} for lnk in range(len(self.AnnoR1.bead.linker_sets))]
        self.one_ins_cell_keys_3 = [{} for lnk in range(len(self.AnnoR1.bead.linker_sets))]

        logging.info("Construct cell label calling dictionaries")
        self.__construct_cell_key_1_hash_table()
        self.__construct_cell_key_2_hash_table()
        self.__construct_cell_key_3_hash_table()


    def __initialize_fast_check_ins_lnk_LD_threshold(self):
        self.fast_check_ins_lnk_LD_threshold = np.zeros(
            (len(self.AnnoR1.bead.inserts), len(self.AnnoR1.bead.linker_sets))).astype(int)
        for ins in range(len(self.AnnoR1.bead.inserts)):
            for lnk in range(len(self.AnnoR1.bead.linker_sets)):
                ins_lnk_length = len(self.AnnoR1.bead.inserts[ins]) + len(self.AnnoR1.bead.linker_sets[lnk][0]) + len(
                    self.AnnoR1.bead.linker_sets[lnk][1])
                self.fast_check_ins_lnk_LD_threshold[ins][lnk] = ins_lnk_length * self.AnnoR1.LDpBP_threshold

    def __initialize_fast_part1_collection_data(self):
        fast_part1_start_ends = np.zeros((len(self.AnnoR1.bead.inserts), 2)).astype(int)
        overhang_lengths = np.zeros((len(self.AnnoR1.bead.inserts)))
        for ins in range(len(self.AnnoR1.bead.inserts)):
            CL1_Start_End_ins = [x + len(self.AnnoR1.bead.inserts[ins]) for x in self.AnnoR1.bead.CL1_Start_End]
            if CL1_Start_End_ins[0] > 0:
                fast_part1_start_ends[ins][0] = CL1_Start_End_ins[0] - 1
                overhang_lengths[ins] += 1
            fast_part1_start_ends[ins][1] = CL1_Start_End_ins[1] + 1
            overhang_lengths[ins] += 1

        LD_threshold = self.AnnoR1.LDpBP_threshold * (
                self.AnnoR1.bead.CL1_Start_End[1] - self.AnnoR1.bead.CL1_Start_End[0])
        self.fast_collect_part1_data = [fast_part1_start_ends, overhang_lengths, LD_threshold]

    def __initialize_fast_part2_collection_data(self):
        fast_part2_start_ends = np.zeros((len(self.AnnoR1.bead.inserts), 2)).astype(int)
        overhang_length = 2
        for ins in range(len(self.AnnoR1.bead.inserts)):
            CL2_Start_End_ins = [x + len(self.AnnoR1.bead.inserts[ins]) for x in self.AnnoR1.bead.CL2_Start_End]
            fast_part2_start_ends[ins][0] = CL2_Start_End_ins[0] - 1
            fast_part2_start_ends[ins][1] = CL2_Start_End_ins[1] + 1

        LD_threshold = self.AnnoR1.LDpBP_threshold * (
                self.AnnoR1.bead.CL2_Start_End[1] - self.AnnoR1.bead.CL2_Start_End[0])
        self.fast_collect_part2_data = [fast_part2_start_ends, overhang_length, LD_threshold]

    def __initialize_fast_part3_collection_data(self):
        fast_part3_start_ends = np.zeros((len(self.AnnoR1.bead.inserts), 2)).astype(int)
        overhang_length = 2
        for ins in range(len(self.AnnoR1.bead.inserts)):
            CL3_Start_End_ins = [x + len(self.AnnoR1.bead.inserts[ins]) for x in self.AnnoR1.bead.CL3_Start_End]
            fast_part3_start_ends[ins][0] = CL3_Start_End_ins[0] - 1
            fast_part3_start_ends[ins][1] = CL3_Start_End_ins[1] + 1

        LD_threshold = self.AnnoR1.LDpBP_threshold * (
                self.AnnoR1.bead.CL3_Start_End[1] - self.AnnoR1.bead.CL3_Start_End[0])
        self.fast_collect_part3_data = [fast_part3_start_ends, overhang_length, LD_threshold]

    def __initialize_part1_locations_and_sequences_for_editops(self):
        # Initialize variables
        # lists for storing truth sequences (so that we don't need to reconstruct them repeatedly)
        part1_seq = [[[] for i in range(len(self.AnnoR1.bead.linker_sets))] for j in
                     range(len(self.AnnoR1.bead.inserts))]

        # End location
        part1_start_ends = np.zeros((len(self.AnnoR1.bead.inserts), 2)).astype(int)

        # LD thresholds
        part1_ld_threshold = [[] for i in range(len(self.AnnoR1.bead.inserts))]

        # Number of bases
        part1_len = [0 for i in range(len(self.AnnoR1.bead.inserts))]

        # Define start end positions and calculate thresholds
        lnk1_len = self.AnnoR1.bead.Lnk1_Start_End[1] - self.AnnoR1.bead.Lnk1_Start_End[0]
        lnk1_overhang_len = math.ceil(lnk1_len / 2)
        for ins in range(len(self.AnnoR1.bead.inserts)):
            part1_start_ends[ins][1] = len(self.AnnoR1.bead.inserts[ins]) + self.AnnoR1.bead.Lnk1_Start_End[1]
            part1_len[ins] = part1_start_ends[ins][1] - lnk1_overhang_len
            part1_ld_threshold[ins] = int(round(part1_len[ins] * self.AnnoR1.LDpBP_threshold))

        for ins in range(len(self.AnnoR1.bead.inserts)):
            for lnk_idx in range(len(self.AnnoR1.bead.linker_sets)):
                for cell_key_1 in self.cell_keys_1:
                    # Define sequences for each part
                    part1_seq[ins][lnk_idx].append(
                        self.AnnoR1.bead.inserts[ins] + cell_key_1 + self.AnnoR1.bead.linker_sets[lnk_idx][0])

        part1_seq = np.asarray(part1_seq)
        part1_ld_threshold = np.asarray(part1_ld_threshold)
        self.part1_LD_cal_data = (
            part1_start_ends,
            part1_ld_threshold,
            part1_seq,
            part1_len
        )

    def __initialize_part2_locations_and_sequences_for_editops(self):
        # Initialize variables
        # lists for storing truth sequences (so that we don't need to reconstruct them repeatedly)
        part2_seq = [[] for i in range(len(self.AnnoR1.bead.linker_sets))]

        # Start end location
        part2_start_ends = np.zeros((len(self.AnnoR1.bead.inserts), 2)).astype(int)

        # LD threshold
        part2_ld_threshold = 0

        # Number of bases
        part2_len = 0

        # Define start end positions and calculate thresholds
        lnk1_len = self.AnnoR1.bead.Lnk1_Start_End[1] - self.AnnoR1.bead.Lnk1_Start_End[0]
        lnk2_len = self.AnnoR1.bead.Lnk2_Start_End[1] - self.AnnoR1.bead.Lnk2_Start_End[0]
        lnk1_overhang_len = math.floor(lnk1_len / 2)
        lnk2_overhang_len = math.ceil(lnk2_len / 2)
        for ins in range(len(self.AnnoR1.bead.inserts)):
            part2_start_ends[ins][0] = len(self.AnnoR1.bead.inserts[ins]) + self.AnnoR1.bead.Lnk1_Start_End[0]
            part2_start_ends[ins][1] = len(self.AnnoR1.bead.inserts[ins]) + self.AnnoR1.bead.Lnk2_Start_End[1]

        lnk2_overhang_start = part2_start_ends[0][1] - part2_start_ends[0][0] - lnk2_overhang_len
        part2_len = part2_start_ends[0][1] - part2_start_ends[0][0] - lnk1_overhang_len - lnk2_overhang_len
        part2_ld_threshold = int(round(part2_len * self.AnnoR1.LDpBP_threshold))

        for lnk_idx in range(len(self.AnnoR1.bead.linker_sets)):
            # Define sequences for each part
            for cell_key_2 in self.cell_keys_2:
                part2_seq[lnk_idx].append(
                    self.AnnoR1.bead.linker_sets[lnk_idx][0] + cell_key_2 + self.AnnoR1.bead.linker_sets[lnk_idx][1])

        part2_seq = np.asarray(part2_seq)
        self.part2_LD_cal_data = (
            part2_start_ends,
            part2_ld_threshold,
            part2_seq,
            part2_len,
            lnk1_overhang_len,
            lnk2_overhang_start
        )

    def __initialize_part3_locations_and_sequences_for_editops(self):
        # Initialize variables
        # lists for storing truth sequences (so that we don't need to reconstruct them repeatedly)
        part3_seq = [[] for i in range(len(self.AnnoR1.bead.linker_sets))]

        # Start end location
        part3_start_ends = np.zeros((len(self.AnnoR1.bead.inserts), 2)).astype(int)

        # LD threshold
        part3_ld_threshold = 0

        # Number of bases
        part3_len = 0

        # Define start end positions and calculate thresholds
        lnk2_len = self.AnnoR1.bead.Lnk2_Start_End[1] - self.AnnoR1.bead.Lnk2_Start_End[0]
        lnk2_overhang_len = math.floor(lnk2_len / 2)
        for ins in range(len(self.AnnoR1.bead.inserts)):
            part3_start_ends[ins][0] = len(self.AnnoR1.bead.inserts[ins]) + self.AnnoR1.bead.Lnk2_Start_End[0]
            part3_start_ends[ins][1] = len(self.AnnoR1.bead.inserts[ins]) + self.AnnoR1.bead.CL3_Start_End[1]

        # Define thresholds
        part3_len = part3_start_ends[0][1] - part3_start_ends[0][0] - lnk2_overhang_len
        part3_ld_threshold = int(round(part3_len * self.AnnoR1.LDpBP_threshold))

        for lnk_idx in range(len(self.AnnoR1.bead.linker_sets)):
            # Define sequences for each part
            for cell_key_3 in self.cell_keys_3:
                part3_seq[lnk_idx].append(self.AnnoR1.bead.linker_sets[lnk_idx][1] + cell_key_3)

        part3_seq = np.asarray(part3_seq)

        self.part3_LD_cal_data = (
            part3_start_ends,
            part3_ld_threshold,
            part3_seq,
            part3_len,
            lnk2_overhang_len
        )

    def __construct_cell_key_1_hash_table(self):
        part1_start_ends, part1_ld_thresholds, part1_seqs, part1_lens = self.part1_LD_cal_data
        for ins in range(len(self.AnnoR1.bead.inserts)):
            for lnk in range(len(self.AnnoR1.bead.linker_sets)):
                part1_len = part1_lens[ins]
                # Part 1 location on read 1, by hash table types
                self.part1_ck_start_ends[ins][1] = part1_len
                self.part1_1del_ck_start_ends[ins][1] = part1_len + 1
                self.part1_1ins_ck_start_ends[ins][1] = part1_len + 3
                for ck in range(len(part1_seqs[ins][lnk])):
                    part1_seq = part1_seqs[ins][lnk][ck]
                    perfect_seq = part1_seq[:part1_len]
                    if perfect_seq in self.perfect_cell_keys_1[ins][lnk]:
                        logging.error(f"Error: Duplicate cell key sequence in CL1, insert {ins}, linker set {lnk}.")
                        logging.error(f"Duplicated sequence: {perfect_seq}")
                        logging.error(f"Cell key index and editops in dictionary: {self.perfect_cell_keys_1[ins][lnk][perfect_seq]}")
                        logging.error(f"Duplicated cell key index and editops: {[ck+1, [], 0]}")
                        raise
                    else:
                        self.perfect_cell_keys_1[ins][lnk][perfect_seq] = [[str(ck+1),[], 0]] # sequence is key, cell key index and levenshtein editops are values

                    # Put 1-substitution sequences into hash table
                    mutated_sequences, mutation_information = self.__get_one_sub_sequences(perfect_seq, ck + 1)
                    for ms_idx in range(len(mutated_sequences)):
                        if mutated_sequences[ms_idx] in self.one_sub_cell_keys_1[ins][lnk]:
                            self.one_sub_cell_keys_1[ins][lnk][mutated_sequences[ms_idx]] += [mutation_information[ms_idx]]
                        else:
                            self.one_sub_cell_keys_1[ins][lnk][mutated_sequences[ms_idx]] = [mutation_information[ms_idx]]

                    # Put 1-deletion sequences into hash table
                    perfect_seq = part1_seq[:(part1_len + 2)]
                    mutated_sequences, mutation_information = self.__get_one_del_sequences(perfect_seq, ck + 1, 0, len(perfect_seq)-2)
                    for ms_idx in range(len(mutated_sequences)):
                        if mutated_sequences[ms_idx] in self.one_del_cell_keys_1[ins][lnk]:
                            if mutation_information[ms_idx][0] in self.one_del_cell_keys_1[ins][lnk][mutated_sequences[ms_idx]][:][0]:
                                continue
                            self.one_del_cell_keys_1[ins][lnk][mutated_sequences[ms_idx]] += [mutation_information[ms_idx]]
                        else:
                            self.one_del_cell_keys_1[ins][lnk][mutated_sequences[ms_idx]] = [mutation_information[ms_idx]]

                    # Put 1-insertion sequences into hash table
                    mutated_sequences, mutation_information = self.__get_one_ins_sequences(perfect_seq, ck + 1, 0, len(perfect_seq)-2)
                    for ms_idx in range(len(mutated_sequences)):
                        if mutated_sequences[ms_idx] in self.one_ins_cell_keys_1[ins][lnk]:
                            if mutation_information[ms_idx][0] in self.one_ins_cell_keys_1[ins][lnk][mutated_sequences[ms_idx]][:][0]:
                                continue
                            self.one_ins_cell_keys_1[ins][lnk][mutated_sequences[ms_idx]] += [mutation_information[ms_idx]]
                        else:
                            self.one_ins_cell_keys_1[ins][lnk][mutated_sequences[ms_idx]] = [mutation_information[ms_idx]]

    def __construct_cell_key_2_hash_table(self):
        part2_start_ends, \
        part2_ld_threshold, \
        part2_seqs, \
        part2_len, \
        lnk1_overhang_len, \
        lnk2_overhang_stsart \
            = self.part2_LD_cal_data
        # Part 2 location on read 1, by hash table types
        for ins in range(len(self.AnnoR1.bead.inserts)):
            part2_start_end = part2_start_ends[ins]
            self.part2_ck_start_ends[ins][0] = part2_start_end[0] + lnk1_overhang_len
            self.part2_1del_ck_start_ends[ins][0] = part2_start_end[0] + lnk1_overhang_len - 2
            self.part2_1ins_ck_start_ends[ins][0] = part2_start_end[0] + lnk1_overhang_len - 2
            self.part2_ck_start_ends[ins][1] = part2_start_end[0] + lnk1_overhang_len + part2_len
            self.part2_1del_ck_start_ends[ins][1] = part2_start_end[0] + lnk1_overhang_len + part2_len + 1
            self.part2_1ins_ck_start_ends[ins][1] = part2_start_end[0] + lnk1_overhang_len + part2_len + 3
        part2_end = part2_len + lnk1_overhang_len
        for lnk in range(len(self.AnnoR1.bead.linker_sets)):
            for ck in range(len(cell_keys.v2_cell_key2)):
                part2_seq = part2_seqs[lnk][ck]
                perfect_seq = part2_seq[lnk1_overhang_len:part2_end]
                if perfect_seq in self.perfect_cell_keys_2[lnk]:
                    logging.error(f"Error: Duplicate cell key sequence in CL2, linker set {lnk}.")
                    logging.error(f"Perfect sequence: {perfect_seq}")
                    logging.error(f"Duplicated sequence: {perfect_seq}")
                    logging.error(f"Cell key index and editops in dictionary: {self.perfect_cell_keys_2[lnk][perfect_seq]}")
                    logging.error(f"Duplicated cell key index and editops: {[ck+1, [], 0]}")
                    raise
                else:
                    self.perfect_cell_keys_2[lnk][perfect_seq] = [[str(ck+1),[], 0]] # sequence is key, cell key index and levenshtein editops are values

                # Put 1-substitution sequences into hash table
                mutated_sequences, mutation_information = self.__get_one_sub_sequences(perfect_seq, ck+1)
                for ms_idx in range(len(mutated_sequences)):
                    if mutated_sequences[ms_idx] in self.one_sub_cell_keys_2[lnk]:
                        self.one_sub_cell_keys_2[lnk][mutated_sequences[ms_idx]] += [mutation_information[ms_idx]]
                    else:
                        self.one_sub_cell_keys_2[lnk][mutated_sequences[ms_idx]] = [mutation_information[ms_idx]]

                # Put 1-deletion sequences into hash table
                perfect_seq = part2_seq[(lnk1_overhang_len - 2):(part2_end + 2)]
                mutated_sequences, mutation_information = self.__get_one_del_sequences(perfect_seq, ck + 1, 2, len(perfect_seq)-2)
                for ms_idx in range(len(mutated_sequences)):
                    if mutated_sequences[ms_idx] in self.one_del_cell_keys_2[lnk]:
                        if mutation_information[ms_idx][0] in self.one_del_cell_keys_2[lnk][mutated_sequences[ms_idx]][:][0]:
                            continue
                        self.one_del_cell_keys_2[lnk][mutated_sequences[ms_idx]] += [mutation_information[ms_idx]]
                    else:
                        self.one_del_cell_keys_2[lnk][mutated_sequences[ms_idx]] = [mutation_information[ms_idx]]

                # Put 1-insertion sequences into hash table
                perfect_seq = part2_seq[(lnk1_overhang_len - 2):(part2_end + 2)]
                mutated_sequences, mutation_information = self.__get_one_ins_sequences(perfect_seq, ck + 1, 2, len(perfect_seq) - 2)
                for ms_idx in range(len(mutated_sequences)):
                    if mutated_sequences[ms_idx] in self.one_ins_cell_keys_2[lnk]:
                        if mutation_information[ms_idx][0] in self.one_ins_cell_keys_2[lnk][mutated_sequences[ms_idx]][:][0]:
                            continue
                        self.one_ins_cell_keys_2[lnk][mutated_sequences[ms_idx]] += [mutation_information[ms_idx]]
                    else:
                        self.one_ins_cell_keys_2[lnk][mutated_sequences[ms_idx]] = [mutation_information[ms_idx]]

    def __construct_cell_key_3_hash_table(self):
        part3_start_ends, \
        part3_ld_threshold, \
        part3_seqs, \
        part3_len, \
        lnk2_overhang_len \
            = self.part3_LD_cal_data
        for ins in range(len(self.AnnoR1.bead.inserts)):
            part3_start_end = part3_start_ends[ins]
            self.part3_ck_start_ends[ins][0] = part3_start_end[0] + lnk2_overhang_len
            self.part3_1del_ck_start_ends[ins][0] = part3_start_end[0] + lnk2_overhang_len - 2
            self.part3_1ins_ck_start_ends[ins][0] = part3_start_end[0] + lnk2_overhang_len - 2
            self.part3_ck_start_ends[ins][1] = part3_start_end[1]
            self.part3_1del_ck_start_ends[ins][1] = part3_start_end[1] - 1
            self.part3_1ins_ck_start_ends[ins][1] = part3_start_end[1] + 1
        for lnk in range(len(self.AnnoR1.bead.linker_sets)):
            for ck in range(len(cell_keys.v2_cell_key3)):
                part3_seq = part3_seqs[lnk][ck]
                perfect_seq = part3_seq[lnk2_overhang_len:]
                if perfect_seq not in self.perfect_cell_keys_3[lnk]:
                    self.perfect_cell_keys_3[lnk][perfect_seq] = [[str(ck+1),[], 0]] # sequence is key, cell key index and levenshtein editops are values
                else:
                    logging.error(f"Error: Duplicate cell key sequence in CL3, linker set {lnk}.")
                    logging.error(f"Duplicated sequence: {perfect_seq}")
                    logging.error(f"Cell key index and editops in dictionary: {self.perfect_cell_keys_3[lnk][perfect_seq]}")
                    logging.error(f"Duplicated cell key index and editops: {[ck+1, [], 0]}")
                    raise

                # Put 1-substitution sequences into hash table
                mutated_sequences, mutation_information = self.__get_one_sub_sequences(perfect_seq, ck + 1)
                for ms_idx in range(len(mutated_sequences)):
                    if mutated_sequences[ms_idx] in self.one_sub_cell_keys_3[lnk]:
                        self.one_sub_cell_keys_3[lnk][mutated_sequences[ms_idx]] += [mutation_information[ms_idx]]
                    else:
                        # Update mutation location to reflect its position in the read
                        self.one_sub_cell_keys_3[lnk][mutated_sequences[ms_idx]] = [mutation_information[ms_idx]]

                # Put 1-deletion sequences into hash table
                perfect_seq = part3_seq[(lnk2_overhang_len - 2):]
                mutated_sequences, mutation_information = self.__get_one_del_sequences(perfect_seq, ck + 1, 2, len(perfect_seq) - 2)
                for ms_idx in range(len(mutated_sequences)):
                    if mutated_sequences[ms_idx] in self.one_del_cell_keys_3[lnk]:
                        if mutation_information[ms_idx][0] in self.one_del_cell_keys_3[lnk][mutated_sequences[ms_idx]][:][0]:
                            continue
                        self.one_del_cell_keys_3[lnk][mutated_sequences[ms_idx]] += [mutation_information[ms_idx]]
                    else:
                        self.one_del_cell_keys_3[lnk][mutated_sequences[ms_idx]] = [mutation_information[ms_idx]]

                # Put 1-insertion sequences into hash table
                mutated_sequences, mutation_information = self.__get_one_ins_sequences(perfect_seq, ck + 1, 2, len(perfect_seq) - 2)
                for ms_idx in range(len(mutated_sequences)):
                    if mutated_sequences[ms_idx] in self.one_ins_cell_keys_3[lnk]:
                        if mutation_information[ms_idx][0] in self.one_ins_cell_keys_3[lnk][mutated_sequences[ms_idx]][:][0]:
                            continue
                        self.one_ins_cell_keys_3[lnk][mutated_sequences[ms_idx]] += [mutation_information[ms_idx]]
                    else:
                        self.one_ins_cell_keys_3[lnk][mutated_sequences[ms_idx]] = [mutation_information[ms_idx]]

    @staticmethod
    def __generate_UMIs(UMI_length):
        # Generate all combinations of nucleotides given UMI length.
        if UMI_length <= 0:
            return [""]
        bases = ["A", "T", "C", "G", "N"]
        UMIs = ["A", "T", "C", "G", "N"]
        for i in range(UMI_length - 1):
            tmp_UMI = []
            for j in bases:
                for k in UMIs:
                    tmp_UMI.append(k + j)

            UMIs = tmp_UMI
        return UMIs

    @staticmethod
    def __get_one_sub_sequences(original_sequence, cell_key_index):
        """
        Given a sequence, output all possible 1 substitution mutations from it.
        Input:
            original_sequence: the sequence, where 1 substitution mutations will be introduced into.
            cell_key_index: cell key index of the original sequence
        Output:
            mutated_sequences: a list of sequences that are 1 substitution away from the original_sequence.
            mutation_information: a list of editops objects, in the same order as mutated sequences.
              editops objects:
                Format: [cell key index, [[type of mutation, mutation location]]]
                    type of mutation: int, 0: substitution, 1: deletion, 2: insertion
                    mutation location: location on the original sequence.
                Note: it has the same format as Levenshtein.editops output, except that
                    the type of mutation is recorded in int, and
                    there is no 3rd element in the list.
        """
        bases = ["A", "T", "C", "G", "N"]
        mutated_sequences = []
        mutation_information = []
        for mutation_location in range(len(original_sequence)):
            for base in bases:
                mutated_sequence = copy.deepcopy(original_sequence)
                if mutated_sequence[mutation_location] == base:
                    continue

                if mutation_location != len(mutated_sequence) - 1:
                    mutated_sequence = mutated_sequence[:mutation_location] + base + mutated_sequence[(mutation_location + 1):]
                else:
                    mutated_sequence = mutated_sequence[:mutation_location] + base
                mutation_info = [0, mutation_location]
                # above: 0 indicates substitution (1 is deletion and 2 is insertion, not used here)
                mutation_information.append([str(cell_key_index), [mutation_info], 0])
                mutated_sequences.append(mutated_sequence)

        return mutated_sequences, mutation_information

    @staticmethod
    def __get_one_del_sequences(original_sequence, cell_key_index, min_idx, max_idx):
        """
        Given a sequence, output all possible 1 deletion mutations from it.
        Input:
            original_sequence: the sequence, where 1 substitution mutations will be introduced into.
            cell_key_index: cell key index of the original sequence
            min_idx, max_idx: only add mutation within this min max range
        Output:
            mutated_sequences: a list of sequences that are 1 deletion away from the original_sequence.
            mutation_information: a list of editops objects, in the same order as mutated sequences.
              editops objects:
                Format: [cell key index, [[type of mutation, mutation location]]]
                    type of mutation: int, 0: substitution, 1: deletion, 2: insertion
                    mutation location: location on the original sequence.
                Note: it has the same format as Levenshtein.editops output, except that
                    the type of mutation is recorded in int, and
                    there is no 3rd element in the list.
        """
        mutated_sequences = []
        mutation_information = []
        for mutation_location in range(min_idx, max_idx):
            mutated_sequence = copy.deepcopy(original_sequence)
            if mutation_location == 0:
                mutated_sequence = mutated_sequence[1:]
            elif mutation_location == len(original_sequence) - 1:
                mutated_sequence = mutated_sequence[0:mutation_location]
            else:
                mutated_sequence = mutated_sequence[0:mutation_location] + mutated_sequence[(mutation_location+1):]
            mutation_info = [1, mutation_location - min_idx]
            # above: 0 indicates substitution 1 is deletion and 2 is insertion
            mutation_information.append([str(cell_key_index), [mutation_info], -1])
            mutated_sequences.append(mutated_sequence)

        return mutated_sequences, mutation_information

    @staticmethod
    def __get_one_ins_sequences(original_sequence, cell_key_index, min_idx, max_idx):
        """
        Given a sequence, output all possible 1 insertion mutations from it.
        Input:
            original_sequence: the sequence, where 1 insertion mutations will be introduced into.
            cell_key_index: cell key index of the original sequence
            min_idx, max_idx: only add mutation within this min max range
        Output:
            mutated_sequences: a list of sequences that are 1 insertion away from the original_sequence.
            mutation_information: a list of editops objects, in the same order as mutated sequences.
              editops objects:
                Format: [cell key index, [[type of mutation, mutation location]]]
                    type of mutation: int, 0: substitution, 1: deletion, 2: insertion
                    mutation location: location on the original sequence.
                Note: it has the same format as Levenshtein.editops output, except that
                    the type of mutation is recorded in int, and
                    there is no 3rd element in the list.
        """
        bases = ["A", "T", "C", "G", "N"]
        mutated_sequences = []
        mutation_information = []
        for mutation_location in range(min_idx, max_idx):
            for base in bases:
                mutated_sequence = copy.deepcopy(original_sequence)
                if mutation_location == 0:
                    mutated_sequence = base + mutated_sequence
                elif mutation_location == len(mutated_sequence) - 1:
                    mutated_sequence = mutated_sequence + base
                else:
                    mutated_sequence = mutated_sequence[:mutation_location] + base + mutated_sequence[mutation_location:]
                mutation_info = [2, mutation_location - min_idx]
                # above: 0 indicates substitution 1 is deletion and 2 is insertion
                mutation_information.append([str(cell_key_index), [mutation_info], 1])
                mutated_sequences.append(mutated_sequence)

        return mutated_sequences, mutation_information

    def calculate_part1_editops(self, sequence, ins, lnk, candidate_idx):
        # Get information about part 1
        part1_start_ends, part1_ld_thresholds, part1_seqs, part1_lens = self.AnnoR1.hash_tables.part1_LD_cal_data

        # Part 1's sequence from the read
        part1_seq_from_read1 = sequence[part1_start_ends[ins][0]:part1_start_ends[ins][1]]
        # Proposed truth part 1 sequence
        part1_seq_truth = part1_seqs[ins][lnk][candidate_idx - 1]
        # Length of part 1
        part1_len = part1_lens[ins]
        # Baseline part 1 LD threshold
        part1_ld_threshold = part1_ld_thresholds[ins]

        # Calculate Levenshtein editops between read 1 and the proposed truth
        #editops = Levenshtein.editops(part1_seq_truth, part1_seq_from_read1)
        editops = lv.editops(part1_seq_truth, part1_seq_from_read1).as_list()
        
        LD = len(editops)

        # Trim overhangs of part 1 from editops
        for LD in reversed(range(len(editops))):
            if editops[LD][1] < part1_len:
                LD += 1
                break

        if LD <= part1_ld_threshold:
            final_editops = editops[:LD]
            part1_prob = self.AnnoR1.error_rate_tables.get_part1_prob(final_editops, ins, lnk)
            if part1_prob < self.AnnoR1.prob_valid_part_threshold:
                # If the probability of the read belonging to the current candidate is < 50%. Skip.
                return -1
            # Calculate indel adjustments caused by this candidate's error
            ck1_indel_adjust = 0
            for i in range(len(final_editops)):
                if final_editops[i][0] == 1:
                    ck1_indel_adjust -= 1
                elif final_editops[i][0] == 2:
                    ck1_indel_adjust += 1
            return [[str(candidate_idx), final_editops, ck1_indel_adjust]]

        return -1

    def calculate_part2_editops(self, sequence, ck1_indel_adjust, lnk, candidate_idx):
        # Get information about part 2
        part2_start_ends, \
        part2_ld_threshold, \
        part2_seqs, \
        part2_len, \
        lnk1_overhang_len, \
        lnk2_overhang_start \
            = self.AnnoR1.hash_tables.part2_LD_cal_data

        # Part 2 from the read and the truth candidate sequence
        part2_seq_from_read1 = sequence[(part2_start_ends[0][0] + ck1_indel_adjust):(part2_start_ends[0][1] + ck1_indel_adjust)]
        part2_seq_truth = part2_seqs[lnk][candidate_idx - 1]

        # Calculate Levenshtein editops between read 1 and the proposed truth
        #editops = Levenshtein.editops(part2_seq_truth, part2_seq_from_read1)
        editops = lv.editops(part2_seq_truth, part2_seq_from_read1).as_list()

        if len(editops) == 0:
            return [[str(candidate_idx), [], 0]]

        # Trim head and tail overhangs
        i_head = 0
        i_end = len(editops)
        while i_head < len(editops):
            if editops[i_head][1] >= lnk1_overhang_len:
                break
            i_head += 1
        for i_end in reversed(range(len(editops))):
            if editops[i_end][1] < lnk2_overhang_start:
                i_end += 1
                break

        if i_end - i_head <= part2_ld_threshold:
            final_editops = [(editops[i][0], editops[i][1] - lnk1_overhang_len) for i in range(i_head, i_end)]
            part2_prob = self.AnnoR1.error_rate_tables.get_part2_prob(final_editops, lnk)
            if part2_prob < self.AnnoR1.prob_valid_part_threshold:
                # If the probability of the read belonging to the current candidate is < 50%. Skip.
                return -1
            # Calculate indel adjustments caused by this candidate's error
            ck2_indel_adjust = 0
            for i in range(len(final_editops)):
                if final_editops[i][0] == 1:
                    ck2_indel_adjust -= 1
                elif final_editops[i][0] == 2:
                    ck2_indel_adjust += 1
            return [[str(candidate_idx), final_editops, ck2_indel_adjust]]

        return -1

    def calculate_part3_editops(self, sequence, ck2_indel_adjust, lnk, candidate_idx):
        part3_start_ends, \
        part3_ld_threshold, \
        part3_seqs, \
        part3_len, \
        lnk2_overhang_len \
            = self.AnnoR1.hash_tables.part3_LD_cal_data

        # Part 3 sequences from the read and the candidate
        part3_seq_from_read1 = sequence[(part3_start_ends[0][0] + ck2_indel_adjust):(part3_start_ends[0][1] + ck2_indel_adjust)]
        part3_seq_truth = part3_seqs[lnk][candidate_idx - 1]

        # Calculate Levenshtein editops between read 1 and the proposed truth
        #editops = Levenshtein.editops(part3_seq_truth, part3_seq_from_read1)
        editops = lv.editops(part3_seq_truth, part3_seq_from_read1).as_list()

        # Trim overhang
        i_head = 0
        while i_head < len(editops):
            if editops[i_head][1] >= lnk2_overhang_len:
                break
            i_head += 1
        editops = editops[i_head:]

        if len(editops) == 0:
            return [[str(candidate_idx), [], 0]]

        # Get the range of end of read 1 (or umi starting) location based on the initial editops indels
        umi_loc_tentative_ins = 0
        umi_loc_tentative_del = 0
        for ed_idx in reversed(range(len(editops))):
            if editops[ed_idx][0] == "insert":
                umi_loc_tentative_ins -= 1
            elif editops[ed_idx][0] == "delete":
                umi_loc_tentative_del += 1

        if umi_loc_tentative_ins == 0 and umi_loc_tentative_del == 0:
            if len(editops) <= part3_ld_threshold:
                final_editops = [(editops[i][0], editops[i][1] - lnk2_overhang_len) for i in range(len(editops))]
                part3_prob = self.AnnoR1.error_rate_tables.get_part3_prob(final_editops, lnk)
                if part3_prob < self.AnnoR1.prob_valid_part_threshold:
                    # If the probability of the read belonging to the current candidate is < 50%. Skip.
                    return -1
                return [[str(candidate_idx), final_editops, 0]]
            else:
                return -1

        # Fine tune end of read 1 (and umi starting) location, so that LD is minimized.
        editops = []
        current_best_LD = 99
        final_umi_loc = 0
        for umi_test in reversed(range(umi_loc_tentative_ins, umi_loc_tentative_del)):
            tmp_editops = lv.editops(part3_seq_truth,
                                              sequence[(part3_start_ends[0][0] + ck2_indel_adjust):(part3_start_ends[0][1] + ck2_indel_adjust + umi_test)]).as_list()
            #tmp_editops = Levenshtein.editops(part3_seq_truth,
            #                                  sequence[(part3_start_ends[0][0] + ck2_indel_adjust):(part3_start_ends[0][1] + ck2_indel_adjust + umi_test)])
            i_head = 0
            while i_head < len(tmp_editops):
                if tmp_editops[i_head][1] >= lnk2_overhang_len:
                    break
                i_head += 1
            tmp_editops = tmp_editops[i_head:]
            if len(tmp_editops) != 0 and tmp_editops[-1][1] >= len(part3_seq_truth):
                continue
            if len(tmp_editops) < current_best_LD:
                current_best_LD = len(tmp_editops)
                editops = tmp_editops
                final_umi_loc = umi_test

        # Adjust indices in editops to reflect real locations on read 1
        if current_best_LD <= part3_ld_threshold:
            final_editops = [(editops[i][0], editops[i][1] - lnk2_overhang_len) for i in range(len(editops))]
            part3_prob = self.AnnoR1.error_rate_tables.get_part3_prob(final_editops, lnk)
            if part3_prob < self.AnnoR1.prob_valid_part_threshold:
                # If the probability of the read belonging to the current candidate is < 50%. Skip.
                return -1
            return [[str(candidate_idx), final_editops, final_umi_loc ]]

        return -1
