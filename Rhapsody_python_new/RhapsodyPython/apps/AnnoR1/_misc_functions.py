# distutils: extra_compile_args = ["-O3"]
# cython: language_level=3
import numpy as np
import json

class misc_functions:
    # This class contains functions that support the other classes. It contains:
    # 1) hash table for translating a base's location into its (V1, not V2 parts) part number.
    # 2) stat_counters for read 1 quality stats
    # 3) Available insert/linker combinations in this sample
    def __init__(self, AnnoR1):
        self.AnnoR1 = AnnoR1

        # list for fast translation of base location into cell label part index (0, 1 or 2)
        self.base_to_part = []
        self.__initialize_base_to_part()
        
        # stat_counters
        self.stat_counters = dict(total_num_reads=0, cellular_read=0, p2_cellular_read=0,
                                  estimated_cellular_read=0,
                             read_count_breakdown=[[0 for _ in range(len(self.AnnoR1.bead.linker_sets))]
                                                   for _ in range(len(self.AnnoR1.bead.inserts))],
                             perfect_CL_count_breakdown=[[0 for _ in range(len(self.AnnoR1.bead.linker_sets))]
                                                         for _ in range(len(self.AnnoR1.bead.inserts))],
                             illegible_insert_linker_count=0, 
                             ambiguous_insert_linker_count=0,
                             unique_insert_linker_count=0)


        self.num_reads_1LD = [[1 for _ in range(len(self.AnnoR1.bead.linker_sets))]
                                                         for _ in range(len(self.AnnoR1.bead.inserts))]
        self.num_reads_2LD = [[1 for _ in range(len(self.AnnoR1.bead.linker_sets))]
                              for _ in range(len(self.AnnoR1.bead.inserts))]
        self.additional_LD1_reads = 0
        self.additional_LD2_reads = 0

        self.cell_counter = {}
        # Static cell read count for deterministic cell label calling
        self.static_cell_read_count = {}

        self.all_ins_lnk_combinations = {}
        self.__generate_all_ins_lnk_comb()

    def filter_ins_lnk_combinations(self):
        # Remove redundant ins lnk combinations from consideration
        
        if self.stat_counters['cellular_read'] == 0:
            print('Warning: no valid cell labels at this point.  Keeping all existing ins lnk combinations')
            return

        new_ins_lnk_combinations = {}
        for ins in self.all_ins_lnk_combinations:
            for lnk in self.all_ins_lnk_combinations[ins]:
                portion_of_reads = self.stat_counters['read_count_breakdown'][ins][lnk]/self.stat_counters['cellular_read']
                if portion_of_reads < 0.01:
                    continue
                if ins not in new_ins_lnk_combinations:
                    new_ins_lnk_combinations[ins] = [lnk]
                elif lnk not in new_ins_lnk_combinations[ins]:
                    new_ins_lnk_combinations[ins].append(lnk)
        self.all_ins_lnk_combinations = new_ins_lnk_combinations

    def __generate_all_ins_lnk_comb(self):
        for ins in range(len(self.AnnoR1.bead.inserts)):
            self.all_ins_lnk_combinations[ins] = []
            for lnk in range(len(self.AnnoR1.bead.linker_sets)):
                self.all_ins_lnk_combinations[ins].append(lnk)

    def __initialize_base_to_part(self):
        self.base_to_part = np.zeros((len(self.AnnoR1.bead.inserts), self.AnnoR1.bead.UMI_start_end_pos[1] + len(self.AnnoR1.bead.inserts[-1]))).astype(int)
        for ins in range(len(self.AnnoR1.bead.inserts)):
            for pos in range(self.AnnoR1.bead.CL1_Start_End[0] + len(self.AnnoR1.bead.inserts[ins]), self.AnnoR1.bead.Lnk1_Start_End[1] + len(self.AnnoR1.bead.inserts[ins])):
                self.base_to_part[ins][pos] = 0
            for pos in range(self.AnnoR1.bead.CL2_Start_End[0] + len(self.AnnoR1.bead.inserts[ins]), self.AnnoR1.bead.Lnk2_Start_End[1] + len(self.AnnoR1.bead.inserts[ins])):
                self.base_to_part[ins][pos] = 1
            for pos in range(self.AnnoR1.bead.CL3_Start_End[0] + len(self.AnnoR1.bead.inserts[ins]), self.AnnoR1.bead.UMI_start_end_pos[1] + len(self.AnnoR1.bead.inserts[ins])):
                self.base_to_part[ins][pos] = 2
            self.base_to_part[ins] = np.asarray(self.base_to_part[ins])

    def editops_to_output(self, insert, editops):
        LD_count = [0,0,0]
        for editop in editops:
            LD_count[self.base_to_part[insert][editop[1]]] += 1
        return "-".join([str(LD) for LD in LD_count])

    def update_linker_stats(self, this_cell_key_caller):
        if this_cell_key_caller.unique_insert_linker:
            self.stat_counters["unique_insert_linker_count"] += 1
        elif this_cell_key_caller.illegible_insert_linker:
            self.stat_counters["illegible_insert_linker_count"] += 1
        elif this_cell_key_caller.ambiguous_insert_linker:
            self.stat_counters["ambiguous_insert_linker_count"] += 1

    def update_counts(self, this_cell_key_caller):
        ins = this_cell_key_caller.insert_idx
        lnk = this_cell_key_caller.linker_idx
        if this_cell_key_caller.accepted_cell_key:
            self.stat_counters["cellular_read"] += 1
            self.stat_counters["read_count_breakdown"][ins][lnk] += 1
            if len(this_cell_key_caller.editops) == 0:
                self.stat_counters["perfect_CL_count_breakdown"][ins][lnk] += 1
            if this_cell_key_caller.cell_key in self.cell_counter:
                self.cell_counter[this_cell_key_caller.cell_key] += 1
            else:
                self.cell_counter[this_cell_key_caller.cell_key] = 1

    def update_counts_second_pass(self, this_cell_key_caller):
        # Update function for reads that are called.
        ins = this_cell_key_caller.insert_idx
        lnk = this_cell_key_caller.linker_idx
        if this_cell_key_caller.accepted_cell_key:
            self.stat_counters["cellular_read"] += 1
            self.stat_counters["read_count_breakdown"][ins][lnk] += 1
            if len(this_cell_key_caller.editops) == 0:
                self.stat_counters["perfect_CL_count_breakdown"][ins][lnk] += 1
            elif len(this_cell_key_caller.editops) == 1:
                self.num_reads_1LD[ins][lnk] += 1
            elif len(this_cell_key_caller.editops) == 2:
                self.num_reads_2LD[ins][lnk] += 1
            if this_cell_key_caller.cell_key in self.cell_counter:
                self.cell_counter[this_cell_key_caller.cell_key] += 1
            else:
                self.cell_counter[this_cell_key_caller.cell_key] = 1

    def update_counts_second_pass2(self, this_cell_key_caller):
        # Update function for eads that are not called
        if this_cell_key_caller.LD1_found:
            self.additional_LD1_reads += 1
        elif this_cell_key_caller.LD2_found:
            self.additional_LD2_reads += 1

    def output_stat_counters(self):
        output_file = "{}_R1_read_count_breakdown.json".format(self.AnnoR1.library)
        with open(output_file, "w") as outfile:
            json.dump(self.stat_counters, outfile)
