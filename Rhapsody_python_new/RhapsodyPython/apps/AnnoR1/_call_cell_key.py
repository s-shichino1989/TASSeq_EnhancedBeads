# distutils: extra_compile_args = ["-O3"]
# cython: language_level=3
#import Levenshtein
from AnnoR1 import cell_keys
import numpy as np
from rapidfuzz.distance import Levenshtein
from polyleven import levenshtein

class cell_key_caller:
    def __init__(self, AnnoR1, sequence):
        self.sequence = sequence
        self.AnnoR1 = AnnoR1

        # tmp storage for potential candidates
        # Candidates (ck1s, ck2s and ck3s) are stored in a tree structure:
        # The first layer is the indel adjustments from the previous part,
        # the second layer is the linker set index.
        # The ckx_indel_adjusts objects save the tree of available indel adjustments and
        # linker set indices from their respective parts (ck1 or ck2); so that we can search for the next
        # part's candidates (parts 2 and 3) using them.
        self.candidates = dict(ck1s={}, ck1_indel_adjusts={},
                               ck2s={}, ck2_indel_adjusts={},
                               ck3s={})

        # Insert/Linker quality information
        self.unique_insert_linker: bool = False
        self.ambiguous_insert_linker: bool = False
        self.illegible_insert_linker: bool = False
        self.potential_ins_lnk_combinations = {}

        # Current best candidate information
        self.cell_key: str = ""
        self.insert_idx: int = 0
        self.linker_idx: int = 0
        self.editops: list = []
        self.umi: str = ""
        self.accepted_cell_key: bool = False

        # Second pass candidates' LD tracking
        self.LD1_found = False
        self.LD2_found = False

    def call_perfect_cell_keys(self):
        for ins in self.AnnoR1.misc_functions.all_ins_lnk_combinations:
            for lnk in self.AnnoR1.misc_functions.all_ins_lnk_combinations[ins]:
                # Call part 1 cell key
                part1 = self.sequence[self.AnnoR1.hash_tables.part1_ck_start_ends[ins][0]:self.AnnoR1.hash_tables.part1_ck_start_ends[ins][1]]
                if part1 not in self.AnnoR1.hash_tables.perfect_cell_keys_1[ins][lnk]:
                    continue
                # If cell key 1 has a perfect match, make it known to the second pass if the algorithm ever gets there.
                self.__update_potential_ins_lnk_combinations(self.potential_ins_lnk_combinations, ins, lnk)

                # Call part2 cell key
                part2 = self.sequence[self.AnnoR1.hash_tables.part2_ck_start_ends[ins][0]:self.AnnoR1.hash_tables.part2_ck_start_ends[ins][1]]
                if part2 not in self.AnnoR1.hash_tables.perfect_cell_keys_2[lnk]:
                    continue

                # Call part 3 cell key
                part3 = self.sequence[self.AnnoR1.hash_tables.part3_ck_start_ends[ins][0]:self.AnnoR1.hash_tables.part3_ck_start_ends[ins][1]]
                if part3 not in self.AnnoR1.hash_tables.perfect_cell_keys_3[lnk]:
                    continue

                # Cell key is found, record information and return
                ck1, ck1_mutation, _ = self.AnnoR1.hash_tables.perfect_cell_keys_1[ins][lnk][part1][0]
                ck2, ck2_mutation, _ = self.AnnoR1.hash_tables.perfect_cell_keys_2[lnk][part2][0]
                ck3, ck3_mutation, _ = self.AnnoR1.hash_tables.perfect_cell_keys_3[lnk][part3][0]

                self.cell_key = ",".join([ck1, ck2, ck3])
                self.insert_idx = ins
                self.linker_idx = lnk
                self.editops = []
                self.umi = self.sequence[(ins + self.AnnoR1.bead.UMI_start_end_pos[0]):(ins + self.AnnoR1.bead.UMI_start_end_pos[1])]
                self.accepted_cell_key = True
                self.unique_insert_linker = True
                return

    def call_cell_keys_with_1LD_hash_tables(self):
        # Objects for storing the best candidate
        best_cell_key = "x,x,x"
        best_cell_key_LD = 99
        best_read_count = 0
        editops = []
        insert_idx = 0
        linker_idx = 0
        umi_loc = 0
        total_candidate_read_count = 0

        # Obtain potential insert/linker combinations
        if self.potential_ins_lnk_combinations:
            ins_lnk_combs = self.potential_ins_lnk_combinations
        else:
            ins_lnk_combs = self.AnnoR1.misc_functions.all_ins_lnk_combinations

        # Get part1 candidates from given insert/linker combinations
        self.part1_hash_table_candidates(ins_lnk_combs)
        if not self.candidates["ck1s"]:
            return

        # Get part 2 candidates from given cell key 1 indel adjustment and linker combinations
        self.part2_hash_table_candidates()
        if not self.candidates["ck2s"]:
            return

        # Get part 3 candidates from given cell key 2 indel adjustment and linker combinations
        self.part3_hash_table_candidates()
        if not self.candidates["ck3s"]:
            return

        # check all combinations of ck1, ck2 and ck3 for the best cell candidate
        for ins in self.candidates["ck1s"]:
            for lnk in self.candidates["ck1s"][ins]:
                for ck1_info in self.candidates["ck1s"][ins][lnk]:
                    ck1_indel_adjust= ck1_info[2]
                    # Indel adjustment from the current ck1 candidate
                    ck1_indel_adjust += ins
                    if ck1_indel_adjust not in self.candidates["ck2s"]:
                        continue
                    if lnk not in self.candidates["ck2s"][ck1_indel_adjust]:
                        continue
                    # Cell key 1 candidate identity
                    ck1 = ck1_info[0]
                    # Cell key 1 candidate Levenshtein editops
                    ck1_editops = ck1_info[1]
                    for ck2_info in self.candidates["ck2s"][ck1_indel_adjust][lnk]:
                        ck2_indel_adjust = ck2_info[2]
                        # Indel adjustment from the current ck2 candidate
                        ck2_indel_adjust += ck1_indel_adjust
                        if ck2_indel_adjust not in self.candidates["ck3s"]:
                            continue
                        if lnk not in self.candidates["ck3s"][ck2_indel_adjust]:
                            continue
                        # Cell key 2 candidate identity
                        ck2 = ck2_info[0]
                        # Cell key 2 candidate Levenshtein editops
                        ck2_editops = ck2_info[1]
                        for ck3_info in self.candidates["ck3s"][ck2_indel_adjust][lnk]:
                            ck3_indel_adjust = ck3_info[2]
                            # Indel adjustment from the current ck3 candidate. Need this for UMI location
                            ck3_indel_adjust += ck2_indel_adjust
                            # Cell key 3 candidate identity
                            ck3 = ck3_info[0]
                            # Cell key 3 candidate Levenshtein editops
                            ck3_editops = ck3_info[1]
                            # Candidate cell label
                            this_cell = ",".join([ck1,ck2,ck3])
                            # Candidate cell's LD
                            LD = len(ck1_editops) + len(ck2_editops) + len(ck3_editops)

                            # Get candidate cell's read count
                            if this_cell in self.AnnoR1.misc_functions.static_cell_read_count:
                                cell_read_count = self.AnnoR1.misc_functions.static_cell_read_count[this_cell]
                                if LD == 1:
                                    self.LD1_found = True
                                elif LD == 2:
                                    self.LD2_found = True
                            else:
                                cell_read_count = 0.5
                            if best_cell_key != this_cell:
                                total_candidate_read_count += cell_read_count

                            # If no candidate can have both 1) lowest LD and 2) highest cell read count
                            # We will not have a best candidate.
                            if cell_read_count < best_read_count:
                                if LD < best_cell_key_LD:
                                    best_cell_key = "x,x,x"
                                    best_cell_key_LD = LD
                            elif cell_read_count > best_read_count:
                                if LD <= best_cell_key_LD:
                                    # The current candidate has the lowest LD and highest cell read count. Remember it.
                                    best_cell_key = this_cell
                                    best_cell_key_LD = LD
                                    best_read_count = cell_read_count
                                    insert_idx = ins
                                    linker_idx = lnk
                                    umi_loc = ins + ck3_indel_adjust
                                    editops = [ck1_editops, ck2_editops, ck3_editops]
                                else:
                                    best_cell_key = "x,x,x"
                                    best_read_count = cell_read_count
                            elif LD < best_cell_key_LD:
                                # The current candidate has the lowest LD and highest cell read count. Remember it.
                                best_cell_key = this_cell
                                best_cell_key_LD = LD
                                best_read_count = cell_read_count
                                insert_idx = ins
                                linker_idx = lnk
                                umi_loc = ins + ck3_indel_adjust
                                editops = [ck1_editops, ck2_editops, ck3_editops]
                            elif LD == best_cell_key_LD and best_cell_key != this_cell:
                                if best_cell_key != this_cell:
                                    best_cell_key = "x,x,x"
                                else:
                                    # Edge case caused by ambiguous error
                                    if editops[0] != ck1_editops:
                                        if editops[0][0] > ck1_editops[0]:
                                            editops[0] = ck1_editops
                                    if editops[1] != ck2_editops:
                                        if editops[1][0] > ck2_editops[0]:
                                            editops[1] = ck2_editops
                                    if editops[2] != ck3_editops:
                                        if editops[2][0] > ck3_editops[0]:
                                            editops[2] = ck3_editops

        if "x" not in best_cell_key and best_read_count / total_candidate_read_count > 0.98:
            # The best candidate is accepted.
            self.cell_key = best_cell_key
            self.insert_idx = insert_idx
            self.linker_idx = linker_idx
            # Calculate ck2 and ck3 mutation locations based on insert identity
            ck2_editops = [
                [editops[1][i][0], editops[1][i][1] + self.AnnoR1.hash_tables.part2_ck_start_ends[insert_idx][0]] \
                for i in range(len(editops[1]))
            ]
            ck3_editops = [
                [editops[2][i][0], editops[2][i][1] + self.AnnoR1.hash_tables.part3_ck_start_ends[insert_idx][0]] \
                for i in range(len(editops[2]))
            ]
            self.editops = editops[0] + ck2_editops + ck3_editops
            UMI_start = insert_idx + umi_loc + self.AnnoR1.bead.UMI_start_end_pos[0]
            UMI_end = insert_idx + umi_loc + self.AnnoR1.bead.UMI_start_end_pos[1]
            self.umi = self.sequence[UMI_start:UMI_end]
            self.accepted_cell_key = True

    ''' This function finds part 1 candidates using hash tables'''
    def part1_hash_table_candidates(self, ins_lnk_combos):
        for ins in ins_lnk_combos:
            part1 = self.sequence[self.AnnoR1.hash_tables.part1_ck_start_ends[ins][0]:self.AnnoR1.hash_tables.part1_ck_start_ends[ins][1]]
            part1_del = ""
            part1_ins = ""
            for lnk in ins_lnk_combos[ins]:
                # Perfect hash table
                if part1 in self.AnnoR1.hash_tables.perfect_cell_keys_1[ins][lnk]:
                    ck1_info = self.AnnoR1.hash_tables.perfect_cell_keys_1[ins][lnk][part1]
                    self.__update_candidate_indel_info(self.candidates["ck1_indel_adjusts"], ins, lnk)
                    self.__update_candidate_ck_info(self.candidates["ck1s"], ins, lnk, ck1_info)
                    self.__update_potential_ins_lnk_combinations(self.potential_ins_lnk_combinations, ins, lnk)
                    continue
                # 1-substitution hash table
                if part1 in self.AnnoR1.hash_tables.one_sub_cell_keys_1[ins][lnk]:
                    ck1_info = self.AnnoR1.hash_tables.one_sub_cell_keys_1[ins][lnk][part1]
                    self.__update_candidate_indel_info(self.candidates["ck1_indel_adjusts"], ins, lnk)
                    self.__update_candidate_ck_info(self.candidates["ck1s"], ins, lnk, ck1_info)
                    self.__update_potential_ins_lnk_combinations(self.potential_ins_lnk_combinations, ins, lnk)

                # 1-deletion hash table
                if part1_del == "":
                    part1_del = self.sequence[self.AnnoR1.hash_tables.part1_1del_ck_start_ends[ins][0]:self.AnnoR1.hash_tables.part1_1del_ck_start_ends[ins][1]]
                if part1_del in self.AnnoR1.hash_tables.one_del_cell_keys_1[ins][lnk]:
                    ck1_info = self.AnnoR1.hash_tables.one_del_cell_keys_1[ins][lnk][part1_del]
                    this_indel_adjust = ins - 1
                    self.__update_candidate_indel_info(self.candidates["ck1_indel_adjusts"], this_indel_adjust, lnk)
                    self.__update_candidate_ck_info(self.candidates["ck1s"], ins, lnk, ck1_info)
                    self.__update_potential_ins_lnk_combinations(self.potential_ins_lnk_combinations, ins, lnk)

                # 1-insertion hash table
                if part1_ins == "":
                    part1_ins = self.sequence[self.AnnoR1.hash_tables.part1_1ins_ck_start_ends[ins][0]:self.AnnoR1.hash_tables.part1_1ins_ck_start_ends[ins][1]]
                if part1_ins in self.AnnoR1.hash_tables.one_ins_cell_keys_1[ins][lnk]:
                    ck1_info = self.AnnoR1.hash_tables.one_ins_cell_keys_1[ins][lnk][part1_ins]
                    this_indel_adjust = ins + 1
                    self.__update_candidate_indel_info(self.candidates["ck1_indel_adjusts"], this_indel_adjust, lnk)
                    self.__update_candidate_ck_info(self.candidates["ck1s"], ins, lnk, ck1_info)
                    self.__update_potential_ins_lnk_combinations(self.potential_ins_lnk_combinations, ins, lnk)


    ''' This function finds part 2 candidates using hash tables'''
    def part2_hash_table_candidates(self):
        for ck1_indel_adjust in self.candidates["ck1_indel_adjusts"]:
            start = ck1_indel_adjust + self.AnnoR1.hash_tables.part2_ck_start_ends[0][0]
            end = ck1_indel_adjust + self.AnnoR1.hash_tables.part2_ck_start_ends[0][1]
            part2 = self.sequence[start:end]
            part2_del = ""
            part2_ins = ""
            for lnk in self.candidates["ck1_indel_adjusts"][ck1_indel_adjust]:
                # Perfect read
                if part2 in self.AnnoR1.hash_tables.perfect_cell_keys_2[lnk]:
                    ck2_info = self.AnnoR1.hash_tables.perfect_cell_keys_2[lnk][part2]
                    self.__update_candidate_indel_info(self.candidates["ck2_indel_adjusts"], ck1_indel_adjust, lnk)
                    self.__update_candidate_ck_info(self.candidates["ck2s"], ck1_indel_adjust, lnk, ck2_info)
                    continue

                # Substitution
                if part2 in self.AnnoR1.hash_tables.one_sub_cell_keys_2[lnk]:
                    ck2_info = self.AnnoR1.hash_tables.one_sub_cell_keys_2[lnk][part2]
                    self.__update_candidate_indel_info(self.candidates["ck2_indel_adjusts"], ck1_indel_adjust, lnk)
                    self.__update_candidate_ck_info(self.candidates["ck2s"], ck1_indel_adjust, lnk, ck2_info)

                # Deletion
                if part2_del == "":
                    start = ck1_indel_adjust + self.AnnoR1.hash_tables.part2_1del_ck_start_ends[0][0]
                    end = ck1_indel_adjust + self.AnnoR1.hash_tables.part2_1del_ck_start_ends[0][1]
                    part2_del = self.sequence[start:end]
                if part2_del in self.AnnoR1.hash_tables.one_del_cell_keys_2[lnk]:
                    ck2_info = self.AnnoR1.hash_tables.one_del_cell_keys_2[lnk][part2_del]
                    this_indel_adjust = ck1_indel_adjust - 1
                    self.__update_candidate_indel_info(self.candidates["ck2_indel_adjusts"], this_indel_adjust, lnk)
                    self.__update_candidate_ck_info(self.candidates["ck2s"], ck1_indel_adjust, lnk, ck2_info)

                # Insertion
                if part2_ins == "":
                    start = ck1_indel_adjust + self.AnnoR1.hash_tables.part2_1ins_ck_start_ends[0][0]
                    end = ck1_indel_adjust + self.AnnoR1.hash_tables.part2_1ins_ck_start_ends[0][1]
                    part2_ins = self.sequence[start:end]
                if part2_ins in self.AnnoR1.hash_tables.one_ins_cell_keys_2[lnk]:
                    ck2_info = self.AnnoR1.hash_tables.one_ins_cell_keys_2[lnk][part2_ins]
                    this_indel_adjust = ck1_indel_adjust + 1
                    self.__update_candidate_indel_info(self.candidates["ck2_indel_adjusts"], this_indel_adjust, lnk)
                    self.__update_candidate_ck_info(self.candidates["ck2s"], ck1_indel_adjust, lnk, ck2_info)

    ''' This function finds part 3 candidates using hash tables'''
    def part3_hash_table_candidates(self):
        for ck2_indel_adjust in self.candidates["ck2_indel_adjusts"]:
            start = ck2_indel_adjust + self.AnnoR1.hash_tables.part3_ck_start_ends[0][0]
            end = ck2_indel_adjust + self.AnnoR1.hash_tables.part3_ck_start_ends[0][1]
            part3 = self.sequence[start:end]
            part3_del = ""
            part3_ins = ""
            for lnk in self.candidates["ck2_indel_adjusts"][ck2_indel_adjust]:
                # Perfect read
                if part3 in self.AnnoR1.hash_tables.perfect_cell_keys_3[lnk]:
                    ck3_info = self.AnnoR1.hash_tables.perfect_cell_keys_3[lnk][part3]
                    self.__update_candidate_ck_info(self.candidates["ck3s"], ck2_indel_adjust, lnk, ck3_info)
                    continue

                # Substitution
                if part3 in self.AnnoR1.hash_tables.one_sub_cell_keys_3[lnk]:
                    ck3_info = self.AnnoR1.hash_tables.one_sub_cell_keys_3[lnk][part3]
                    self.__update_candidate_ck_info(self.candidates["ck3s"], ck2_indel_adjust, lnk, ck3_info)

                # Deletion
                if part3_del == "":
                    start = ck2_indel_adjust + self.AnnoR1.hash_tables.part3_1del_ck_start_ends[0][0]
                    end = ck2_indel_adjust + self.AnnoR1.hash_tables.part3_1del_ck_start_ends[0][1]
                    part3_del = self.sequence[start:end]
                if part3_del in self.AnnoR1.hash_tables.one_del_cell_keys_3[lnk]:
                    ck3_info = self.AnnoR1.hash_tables.one_del_cell_keys_3[lnk][part3_del]
                    self.__update_candidate_ck_info(self.candidates["ck3s"], ck2_indel_adjust, lnk, ck3_info)

                # Insertion
                if part3_ins == "":
                    start = ck2_indel_adjust + self.AnnoR1.hash_tables.part3_1ins_ck_start_ends[0][0]
                    end = ck2_indel_adjust + self.AnnoR1.hash_tables.part3_1ins_ck_start_ends[0][1]
                    part3_ins = self.sequence[start:end]
                if part3_ins in self.AnnoR1.hash_tables.one_ins_cell_keys_3[lnk]:
                    ck3_info = self.AnnoR1.hash_tables.one_ins_cell_keys_3[lnk][part3_ins]
                    self.__update_candidate_ck_info(self.candidates["ck3s"], ck2_indel_adjust, lnk, ck3_info)

    def get_potential_ins_lnk_combinations(self):
        """
        Check potential insert size and linker identity of a given read
        output:
            2D list: each row is a insert+linker combination, column 1 is insert index, 2 is linker index

        Algorithm:
            If there are one or more perfect matches, output all perfect matches.
            If there is only one perfect match, output that match.
            If there are no perfect matches, output all matches that are: 1) <= 1 LD away from the lowest LD match that is
              below the threshold 2) <= the second lowest LD/bp match that is below the threshold
        """
        # Calculate LD and LD/bp for all insert and linker combinations
        best_ins_lnk_LD = []
        best_ins_lnk_LDpbp = []
        ins_lnk_combinations = []
        for ins in self.AnnoR1.misc_functions.all_ins_lnk_combinations:
            for lnk in self.AnnoR1.misc_functions.all_ins_lnk_combinations[ins]:
                ins_start = 0
                ins_end = ins
                lnk1_startend = [lnk1_SE + ins for lnk1_SE in self.AnnoR1.bead.Lnk1_Start_End]
                lnk2_startend = [lnk2_SE + ins for lnk2_SE in self.AnnoR1.bead.Lnk2_Start_End]
                ins_lnk_length = ins + len(self.AnnoR1.bead.linker_sets[lnk][0]) + len(
                    self.AnnoR1.bead.linker_sets[lnk][1])
                #LD = Levenshtein.distance(self.sequence[ins_start:ins_end], self.AnnoR1.bead.inserts[ins]) + \
                #     Levenshtein.distance(self.sequence[lnk1_startend[0]:lnk1_startend[1]],
                #                          self.AnnoR1.bead.linker_sets[lnk][0]) + \
                #     Levenshtein.distance(self.sequence[lnk2_startend[0]:lnk2_startend[1]],
                #                          self.AnnoR1.bead.linker_sets[lnk][1])
                LD = levenshtein(self.sequence[ins_start:ins_end], self.AnnoR1.bead.inserts[ins]) + \
                     levenshtein(self.sequence[lnk1_startend[0]:lnk1_startend[1]],
                                          self.AnnoR1.bead.linker_sets[lnk][0]) + \
                     levenshtein(self.sequence[lnk2_startend[0]:lnk2_startend[1]],
                                          self.AnnoR1.bead.linker_sets[lnk][1])
                if LD <= self.AnnoR1.hash_tables.fast_check_ins_lnk_LD_threshold[ins][lnk]:
                    best_ins_lnk_LD.append(LD)
                    best_ins_lnk_LDpbp.append(LD / ins_lnk_length)
                else:
                    best_ins_lnk_LD.append(ins_lnk_length)
                    best_ins_lnk_LDpbp.append(1)
                ins_lnk_combinations.append([ins, lnk])

        # Sort
        sorted_LD_idx = np.argsort(best_ins_lnk_LD)
        sorted_LDpbp_idx = np.argsort(best_ins_lnk_LDpbp)

        # No good insert/linker detected. Bad read.
        if best_ins_lnk_LDpbp[sorted_LDpbp_idx[0]] == 1:
            self.illegible_insert_linker = True
            return

        # Get potential ins lnk combinations by the lowest LD and lowest LD + 1
        lowest_LD = best_ins_lnk_LD[sorted_LD_idx[0]]

        if lowest_LD == 0:
            acceptable_LD = 0
        else:
            acceptable_LD = lowest_LD + 1

        potential_ins_lnk_combinations = [ins_lnk_combinations[sorted_LD_idx[i]] for i in range(len(sorted_LD_idx)) \
                                          if best_ins_lnk_LD[sorted_LD_idx[i]] <= acceptable_LD]

        # Get potential ins lnk combinations by LD/bp and the second lowest LD/bp
        lowest_LDpbp = best_ins_lnk_LDpbp[sorted_LDpbp_idx[0]]
        second_lowest_LDpbp = -1
        for i in range(len(sorted_LDpbp_idx)):
            if best_ins_lnk_LDpbp[sorted_LDpbp_idx[i]] == lowest_LDpbp:
                if ins_lnk_combinations[sorted_LDpbp_idx[i]] not in potential_ins_lnk_combinations:
                    potential_ins_lnk_combinations.append(ins_lnk_combinations[sorted_LDpbp_idx[i]])
            elif second_lowest_LDpbp == -1:
                second_lowest_LDpbp = best_ins_lnk_LDpbp[sorted_LDpbp_idx[i]]
                if second_lowest_LDpbp == 1:
                    break
                if ins_lnk_combinations[sorted_LDpbp_idx[i]] not in potential_ins_lnk_combinations:
                    potential_ins_lnk_combinations.append(ins_lnk_combinations[sorted_LDpbp_idx[i]])
            elif second_lowest_LDpbp == best_ins_lnk_LDpbp[sorted_LDpbp_idx[i]]:
                if ins_lnk_combinations[sorted_LDpbp_idx[i]] not in potential_ins_lnk_combinations:
                    potential_ins_lnk_combinations.append(ins_lnk_combinations[sorted_LDpbp_idx[i]])
            else:
                break

        # Found potential insert/linker combinations, add them to the class object
        if len(potential_ins_lnk_combinations) > 0:
            for ins, lnk in potential_ins_lnk_combinations:
                self.__update_potential_ins_lnk_combinations(self.potential_ins_lnk_combinations, ins, lnk)

        num_combinations = 0
        for ins in self.potential_ins_lnk_combinations:
            for lnk in self.potential_ins_lnk_combinations[ins]:
                num_combinations += 1

        # Insert/linker quality metrics
        if num_combinations == 0:
            self.illegible_insert_linker = True
            return
        elif num_combinations == 1:
            self.unique_insert_linker = True
        else:
            self.ambiguous_insert_linker = True

    def call_cell_keys_with_error_correction(self):
        #if self.potential_ins_lnk_combinations:
        ins_lnk_combs = self.potential_ins_lnk_combinations
        #else:
        #    ins_lnk_combs = self.AnnoR1.misc_functions.all_ins_lnk_combinations
        # Part 1
        self.part1_strcomp_candidates(ins_lnk_combs)
        if len(self.candidates["ck1s"]) == 0:
            return
        # Part 2
        self.part2_strcomp_candidates()
        if not self.candidates["ck2s"]:
            return

        # Part 3
        self.part3_strcomp_candidates()
        if not self.candidates["ck3s"]:
            return

        # Objects for saving the best candidate
        best_prob: float = float("-inf")
        best_cell_key = "x,x,x"
        best_ins: int = 0
        best_lnk: int = 0
        best_editops = []
        best_umi_loc: int = 0
        sum_all_prob = 0

        # Identity the best candidate
        for ins in self.candidates["ck1s"]:
            for lnk in self.candidates["ck1s"][ins]:
                for ck1_info in self.candidates["ck1s"][ins][lnk]:
                    # Indel adjustments for part 2 caused by part 1
                    ck1_indel_adjust = ck1_info[2]
                    ck1_indel_adjust += ins
                    if ck1_indel_adjust not in self.candidates["ck2s"]:
                        continue
                    if lnk not in self.candidates["ck2s"][ck1_indel_adjust]:
                        continue
                    # Part 1 candidate cell key index
                    ck1 = ck1_info[0]
                    # This candidate's Levenshtein editops
                    ck1_editops = ck1_info[1]
                    # Initialize perfect read probability
                    perfect_prob = self.AnnoR1.error_rate_tables.log_prob_perfect_read[ins][lnk]
                    # Update this read's probability with editops
                    ck1_prob = self.AnnoR1.error_rate_tables.update_cl_prob(perfect_prob, ck1_editops, ins, lnk)
                    for ck2_info in self.candidates["ck2s"][ck1_indel_adjust][lnk]:
                        ck2_indel_adjust = ck2_info[2]
                        ck2_indel_adjust += ck1_indel_adjust
                        if ck2_indel_adjust not in self.candidates["ck3s"]:
                            continue
                        if lnk not in self.candidates["ck3s"][ck2_indel_adjust]:
                            continue
                        ck2 = ck2_info[0]
                        ck2_editops = ck2_info[1]
                        ck2_editops = [
                            [ck2_editops[i][0], ck2_editops[i][1] + self.AnnoR1.hash_tables.part2_ck_start_ends[ins][0]]
                            for i in range(len(ck2_editops))
                        ]
                        ck2_prob = self.AnnoR1.error_rate_tables.update_cl_prob(ck1_prob, ck2_editops, ins, lnk)
                        for ck3_info in self.candidates["ck3s"][ck2_indel_adjust][lnk]:
                            ck3_indel_adjust = ck3_info[2]
                            ck3_indel_adjust += ck2_indel_adjust
                            ck3_editops = ck3_info[1]
                            ck3_editops = [
                                [ck3_editops[i][0], ck3_editops[i][1] + self.AnnoR1.hash_tables.part3_ck_start_ends[ins][0]]
                                for i in range(len(ck3_editops))
                            ]
                            ck3_prob = self.AnnoR1.error_rate_tables.update_cl_prob(ck2_prob, ck3_editops, ins, lnk)
                            if ck3_prob < self.AnnoR1.prob_valid_read_threshold:
                                # If the probability of getting this candidate's error rate is < 50%, skip
                                continue

                            # Probability of getting a read from this candidate cell
                            ck3 = ck3_info[0]
                            this_cell = ",".join([ck1, ck2, ck3])
                            if this_cell in self.AnnoR1.misc_functions.static_cell_read_count:
                                prob_cell_and_error = ck3_prob + np.log(self.AnnoR1.misc_functions.static_cell_read_count[this_cell])
                            else:
                                prob_cell_and_error = ck3_prob

                            sum_all_prob += np.exp(ck3_prob)

                            if prob_cell_and_error > best_prob:
                                best_prob = prob_cell_and_error
                                best_cell_key = this_cell
                                best_ins: int = ins
                                best_lnk: int = lnk
                                best_editops = ck1_editops + ck2_editops + ck3_editops
                                best_umi_loc: int = ck3_indel_adjust

        if best_prob == float('-inf'):
            # No best candidate
            return

        # Probability of this read belonging to the best candidate cell
        best_candidate_cell_probability = np.exp(best_prob)/sum_all_prob
        if best_candidate_cell_probability < self.AnnoR1.prob_cell_threshold:
            # If the probability of this read belonging to the best candidate is < 50%, skip.
            return

        self.accepted_cell_key = True
        self.cell_key = best_cell_key
        self.insert_idx = best_ins
        self.linker_idx = best_lnk
        self.editops = best_editops
        final_umi_location = [x + best_umi_loc for x in self.AnnoR1.bead.UMI_start_end_pos]
        self.umi = self.sequence[final_umi_location[0]:final_umi_location[1]]

    ''' This function finds part 1 candidates using string comparison against all 96 truth cell keys'''
    def part1_strcomp_candidates(self, ins_lnk_combo):
        for ins in ins_lnk_combo:
            # Fast filter of candidates
            fast_collected_candidates = self.__fast_collect_part1_candidates(ins)
            for lnk in ins_lnk_combo[ins]:
                if ins in self.candidates["ck1s"]:
                    if lnk in self.candidates["ck1s"][ins]:
                        continue
                for candidate_idx in fast_collected_candidates:
                    # Check Levenshtein editops of each fast candidate and keep the potential candidates
                    ck1_info = self.AnnoR1.hash_tables.calculate_part1_editops(self.sequence, ins, lnk, candidate_idx)
                    if ck1_info == -1:
                        continue
                    ck1_indel_adjust = ck1_info[0][2] + ins
                    self.__update_candidate_indel_info(self.candidates["ck1_indel_adjusts"], ck1_indel_adjust, lnk)
                    self.__update_candidate_ck_info(self.candidates["ck1s"], ins, lnk, ck1_info)

    def __fast_collect_part1_candidates(self, ins):
        part1_start_ends, overhang_lengths, LD_threshold = self.AnnoR1.hash_tables.fast_collect_part1_data
        part1 = self.sequence[part1_start_ends[ins][0]:part1_start_ends[ins][1]]
        candidate_list = []
        for ck1 in range(len(cell_keys.v2_cell_key1)):
            #LD = Levenshtein.distance(cell_keys.v2_cell_key1[ck1], part1) - overhang_lengths[ins]
            LD = levenshtein(cell_keys.v2_cell_key1[ck1], part1) - overhang_lengths[ins]
            if LD <= LD_threshold:
                candidate_list.append(ck1+1)
        return candidate_list

    ''' This function finds part 2 candidates using string comparison against all 96 truth cell keys'''
    def part2_strcomp_candidates(self):
        for ck1_indel_adjust in self.candidates["ck1_indel_adjusts"]:
            # Fast filter of candidates
            fast_collected_candidates = self.__fast_collect_part2_candidates(ck1_indel_adjust)
            for lnk in self.candidates["ck1_indel_adjusts"][ck1_indel_adjust]:
                if ck1_indel_adjust in self.candidates["ck2s"]:
                    if lnk in self.candidates["ck2s"][ck1_indel_adjust]:
                        continue
                for candidate_idx in fast_collected_candidates:
                    # Check Levenshtein editops of each fast candidate and keep the potential candidates
                    ck2_info = self.AnnoR1.hash_tables.calculate_part2_editops(self.sequence, ck1_indel_adjust, lnk, candidate_idx)
                    if ck2_info == -1:
                        continue
                    ck2_indel_adjust = ck2_info[0][2] + ck1_indel_adjust
                    self.__update_candidate_indel_info(self.candidates["ck2_indel_adjusts"], ck2_indel_adjust, lnk)
                    self.__update_candidate_ck_info(self.candidates["ck2s"], ck1_indel_adjust, lnk, ck2_info)
                    if len(ck2_info[0][1]) == 0:
                        break

    def __fast_collect_part2_candidates(self, ck1_indel_adjust):
        part2_start_ends, overhang_length, LD_threshold = self.AnnoR1.hash_tables.fast_collect_part2_data
        part2 = self.sequence[(part2_start_ends[0][0] + ck1_indel_adjust):(part2_start_ends[0][1] + ck1_indel_adjust)]
        candidate_list = []
        for ck2 in range(len(cell_keys.v2_cell_key2)):
            #LD = Levenshtein.distance(cell_keys.v2_cell_key2[ck2], part2) - overhang_length
            LD = levenshtein(cell_keys.v2_cell_key2[ck2], part2) - overhang_length
            if LD <= LD_threshold:
                candidate_list.append(ck2+1)
        return candidate_list

    ''' This function finds part 3 candidates using string comparison against all 96 truth cell keys'''
    def part3_strcomp_candidates(self):
        for ck2_indel_adjust in self.candidates["ck2_indel_adjusts"]:
            # Fast filter of candidates
            fast_collected_candidates = self.__fast_collect_part3_candidates(ck2_indel_adjust)
            for lnk in self.candidates["ck2_indel_adjusts"][ck2_indel_adjust]:
                if ck2_indel_adjust in self.candidates["ck3s"]:
                    if lnk in self.candidates["ck3s"][ck2_indel_adjust]:
                        continue
                for candidate_idx in fast_collected_candidates:
                    # Check Levenshtein editops of each fast candidate and keep the potential candidates
                    ck3_info = self.AnnoR1.hash_tables.calculate_part3_editops(self.sequence, ck2_indel_adjust, lnk, candidate_idx)
                    if ck3_info == -1:
                        continue
                    self.__update_candidate_ck_info(self.candidates["ck3s"], ck2_indel_adjust, lnk, ck3_info)
                    if len(ck3_info[0][1]) == 0:
                        break

    def __fast_collect_part3_candidates(self, ck2_indel_adjust):
        part3_start_ends, overhang_length, LD_threshold = self.AnnoR1.hash_tables.fast_collect_part3_data
        part3 =self.sequence[(part3_start_ends[0][0] + ck2_indel_adjust):(part3_start_ends[0][1] + ck2_indel_adjust)]
        candidate_list = []
        for ck3 in range(len(cell_keys.v2_cell_key3)):
            #LD = Levenshtein.distance(cell_keys.v2_cell_key3[ck3], part3) - overhang_length
            LD = levenshtein(cell_keys.v2_cell_key3[ck3], part3) - overhang_length
            if LD <= LD_threshold:
                candidate_list.append(ck3+1)
        return candidate_list

    '''This function saves cell key information into a dictionary '''
    @staticmethod
    def __update_candidate_ck_info(ck_info_dict, indel_adjust, lnk, ck_info_raw):
        ck_info = ck_info_raw.copy()
        if indel_adjust not in ck_info_dict:
            ck_info_dict[indel_adjust] = {lnk: ck_info}
        elif lnk not in ck_info_dict[indel_adjust] or len(ck_info[0][1]) == 0:
            ck_info_dict[indel_adjust][lnk] = ck_info
        else:
            ck_info_dict[indel_adjust][lnk] += ck_info

    '''This function adds all indel adjustment and linker identity of a given part into ck_indel_info_dict'''
    @staticmethod
    def __update_candidate_indel_info(ck_indel_info_dict, indel_adjust, lnk):
        if indel_adjust not in ck_indel_info_dict:
            ck_indel_info_dict[indel_adjust] = {lnk:True}
        elif lnk not in ck_indel_info_dict[indel_adjust]:
            ck_indel_info_dict[indel_adjust][lnk] = True

    @staticmethod
    def __update_potential_ins_lnk_combinations(potential_ins_lnk_combinations, ins, lnk):
        if ins not in potential_ins_lnk_combinations:
            potential_ins_lnk_combinations[ins] = {lnk:True}
        elif lnk not in potential_ins_lnk_combinations[ins]:
            potential_ins_lnk_combinations[ins][lnk] = True
