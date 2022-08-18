# distutils: extra_compile_args = ["-O3"]
# cython: language_level=3
from AnnoR1 import cell_keys

class bead_specific_settings:
    def __init__(self, bead_version):
        # Bead Specific settings
        # Part locations in sequence to idx dictionaries
        self.bead_version = bead_version
        # Section locations
        self.CL1_Start_End,  self.CL2_Start_End,  self.CL3_Start_End     = None, None, None
        self.Lnk1_Start_End, self.Lnk2_Start_End, self.UMI_start_end_pos = None, None, None
        # insert sequences
        self.inserts = ['']
        # linker sequences
        self.linker_sets = []
        # Cell label 1 range for each insert size
        self.cl1_range = []
        self.__bead_specific_settings(bead_version)

    def __bead_specific_settings(self, bead_version):
        if bead_version == "V1" or bead_version == "V1P":
            self.CL1_Start_End = (0, 9)
            self.Lnk1_Start_End = (9, 21)
            self.CL2_Start_End = (21, 30)
            self.Lnk2_Start_End = (30, 43)
            self.CL3_Start_End = (43, 52)
            self.UMI_start_end_pos = (52, 60)
            self.linker_sets = [[cell_keys.v1_linker1, cell_keys.v1_linker2]]
            self.inserts = ['']
            self.cl1_range = [range(0, 96)]
        else:
            self.CL1_Start_End = (0, 9)
            self.Lnk1_Start_End = (9, 13)
            self.CL2_Start_End = (13, 22)
            self.Lnk2_Start_End = (22, 26)
            self.CL3_Start_End = (26, 35)
            self.UMI_start_end_pos = (35, 43)
            self.linker_sets = [[cell_keys.v2_linker1, cell_keys.v2_linker2],[cell_keys.v2_5p_linker1, cell_keys.v2_5p_linker2]]
            self.inserts = cell_keys.v2_inserts # The last insert in the list must be the longest one.
            self.cl1_range = [range(0, 96), range(24, 48), range(48, 72), range(72, 96)]
