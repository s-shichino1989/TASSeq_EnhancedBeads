# distutils: extra_compile_args = ["-O3"]
# cython: language_level=3
import collections
import MistLogger as logging
import numpy as np

class read1_error_rate_table:
    def __init__(self, AnnoR1):
        self.AnnoR1 = AnnoR1

        # Format: self.__error_count_table[insert size, linker set index, error type (0: substitution, 1: deletion, 2: insertion)), base location]
        self.__error_count_table = np.zeros((len(self.AnnoR1.bead.inserts), len(self.AnnoR1.bead.linker_sets), 3,
                                             self.AnnoR1.bead.CL3_Start_End[1] + len(self.AnnoR1.bead.inserts[-1])))
        self.__estimated_error_count = np.zeros((len(self.AnnoR1.bead.inserts), len(self.AnnoR1.bead.linker_sets), 3,
                                                self.AnnoR1.bead.CL3_Start_End[1] + len(self.AnnoR1.bead.inserts[-1])))

        self.mutation_type_to_idx = {"replace": 0, "delete": 1, "insert": 2, 0: 0, 1: 1, 2: 2}

        # Log probability of getting a perfect read
        self.log_prob_perfect_read = np.zeros((len(self.AnnoR1.bead.inserts), len(self.AnnoR1.bead.linker_sets)))

        # Log probability of getting a perfect part1, 2 or 3 on a read
        self.log_prob_perfect_part1 = np.zeros((len(self.AnnoR1.bead.inserts), len(self.AnnoR1.bead.linker_sets)))
        self.log_prob_perfect_part2 = np.zeros(len(self.AnnoR1.bead.linker_sets))
        self.log_prob_perfect_part3 = np.zeros(len(self.AnnoR1.bead.linker_sets))

        # Pre-calculated variables for updating perfect read probability into error'ed read probability
        # self.__log_prob_error_rate_update_table[insert size, linker set index, error type (0: substitution, 1: deletion, 2: insertion)), base location]
        self.__log_prob_error_rate_update_table = np.zeros((len(self.AnnoR1.bead.inserts),
                                                            len(self.AnnoR1.bead.linker_sets), 3,
                                                            self.AnnoR1.bead.CL3_Start_End[1] + len(
                                                                self.AnnoR1.bead.inserts[-1])))
        self.__log_prob_error_rate_update_table_part2 = []
        self.__log_prob_error_rate_update_table_part3 = []

        # Estimation of final cell label calling stats
        self.estimated_cellular_reads_by_inslnk = np.zeros(
            (len(self.AnnoR1.bead.inserts), len(self.AnnoR1.bead.linker_sets)))
        self.estimated_cellular_reads = 0

    def estimate_final_error_rate(self):
        """
        Estimate final error rates.
        Data: Number of perfect cell labels, 1/2 error cell labels from second pass, error count from second pass.
        Task: Estimate the true error count of the entire data set.
        """
        # Initiate estimate error count object
        self.__estimated_error_count = self.__error_count_table.copy()
        self.__estimated_error_count += 1

        # Total number of legible reads in the data
        total_num_reads = self.AnnoR1.misc_functions.stat_counters['total_num_reads'] - \
                          self.AnnoR1.misc_functions.stat_counters["illegible_insert_linker_count"]
        # Number of reads that have a called cell label up to the second pass
        total_num_called_reads = np.sum(self.AnnoR1.misc_functions.stat_counters['cellular_read'])
        for ins in self.AnnoR1.misc_functions.all_ins_lnk_combinations:
            for lnk in self.AnnoR1.misc_functions.all_ins_lnk_combinations[ins]:
                # Number of cell-label bases on read 1
                num_base = len(self.AnnoR1.bead.inserts[ins]) + self.AnnoR1.bead.CL3_Start_End[1]
                # Number of perfect reads in this insert/linker combination
                num_reads_without_errors_ins_lnk = self.AnnoR1.misc_functions.stat_counters["perfect_CL_count_breakdown"][ins][lnk]
                # Number of 1 LD reads in this insert/linker combination
                num_reads_with_1_error_ins_lnk = self.AnnoR1.misc_functions.num_reads_1LD[ins][lnk]
                # Number of 2 LD reads in this insert/linker combination
                num_reads_with_2_error_ins_lnk = self.AnnoR1.misc_functions.num_reads_2LD[ins][lnk]
                # Percentage of reads that belong to this insert/linker combination
                ins_lnk_read_ratio = self.AnnoR1.misc_functions.stat_counters['read_count_breakdown'][ins][
                                         lnk] / np.sum(self.AnnoR1.misc_functions.stat_counters['read_count_breakdown'])
                # Incorporate potential 1 LD read that were not called in the second pass into the variable
                num_reads_with_1_error_ins_lnk += self.AnnoR1.misc_functions.additional_LD1_reads * ins_lnk_read_ratio
                # Incorporate potential 2 LD read that were not called in the second pass into the variable
                num_reads_with_2_error_ins_lnk += self.AnnoR1.misc_functions.additional_LD2_reads * ins_lnk_read_ratio
                # Total number of reads with cell label called up to the second pass in this insert/linker combination
                num_reads_ins_lnk = self.AnnoR1.misc_functions.stat_counters['read_count_breakdown'][ins][lnk]

                # Total occurrence of error on each position in read 1
                error_count = np.asarray(
                    [np.sum(self.__estimated_error_count[ins, lnk, :, base]) for base in range(num_base)])
                # Estimated total number of reads in this insert/linker combination
                est_total_read_count_ins_lnk = (num_reads_ins_lnk / total_num_called_reads) * total_num_reads

                # Estimate total number of cellular read and the final error count multiplier
                est_num_valid_reads, error_count_multiplier = self.__estimate_num_valid_reads_and_error_rate(
                    error_count, est_total_read_count_ins_lnk,
                    num_reads_without_errors_ins_lnk, num_reads_with_1_error_ins_lnk,
                    num_reads_with_2_error_ins_lnk, num_base, ins)

                # Store them in class variables
                self.estimated_cellular_reads += est_num_valid_reads
                self.estimated_cellular_reads_by_inslnk[ins][lnk] = est_num_valid_reads

                # Update error rates to reflect the final estimated error rate
                self.__estimated_error_count[ins][lnk] *= error_count_multiplier

        self.__estimated_error_count = np.asarray(self.__estimated_error_count)
        self.AnnoR1.misc_functions.stat_counters["estimated_cellular_read"] = self.estimated_cellular_reads

    def __get_2_err_read_count(self, updated_error_count, est_num_valid_reads, num_base, ins):
        log_perfect_count = np.log(est_num_valid_reads - updated_error_count)
        log_error_count = np.log(updated_error_count)
        d_valid_read = (num_base - 1) * np.log(est_num_valid_reads)
        sum_log_perfect_count = np.sum(log_perfect_count)
        est_num_2LD_reads = 0
        for base1 in range(self.AnnoR1.hash_tables.part1_ck_start_ends[ins][1]):
            part1_error_rate = log_error_count[base1] + sum_log_perfect_count - log_perfect_count[base1] - d_valid_read
            for base2 in range(self.AnnoR1.hash_tables.part2_ck_start_ends[ins][0], num_base):
                est_num_2LD_reads += np.exp(part1_error_rate + log_error_count[base2] - log_perfect_count[base2])
        for base1 in range(self.AnnoR1.hash_tables.part2_ck_start_ends[ins][0],
                           self.AnnoR1.hash_tables.part2_ck_start_ends[ins][1]):
            part1_error_rate = log_error_count[base1] + sum_log_perfect_count - log_perfect_count[base1] - d_valid_read
            for base2 in range(self.AnnoR1.hash_tables.part3_ck_start_ends[ins][0],
                               self.AnnoR1.hash_tables.part3_ck_start_ends[ins][1]):
                est_num_2LD_reads += np.exp(part1_error_rate + log_error_count[base2] - log_perfect_count[base2])
        return est_num_2LD_reads

    @staticmethod
    def __get_1_err_read_count(updated_error_count, est_num_valid_reads, num_base):
        log_perfect_count = np.log(est_num_valid_reads - updated_error_count)
        log_error_count = np.log(updated_error_count)
        d_valid_read = (num_base - 1) * np.log(est_num_valid_reads)
        sum_log_perfect_count = np.sum(log_perfect_count)
        est_num_1LD_reads = np.sum(np.exp(
            [log_error_count[base] + sum_log_perfect_count - log_perfect_count[base] - d_valid_read for base in
             range(num_base)]))
        return est_num_1LD_reads

    @staticmethod
    def __get_perfect_read_count(updated_error_count, est_num_valid_reads, num_base):
        log_perfect_count = np.log(est_num_valid_reads - updated_error_count)
        d_valid_read = (num_base - 1) * np.log(est_num_valid_reads)
        return np.exp(np.sum(log_perfect_count) - d_valid_read)

    def __estimate_num_valid_reads_and_error_rate(self, error_count, est_total_read_count, num_reads_without_errors,
                                   num_reads_with_1_error, num_reads_with_2_error, num_base, ins):
        """
        Using the equations of getting the number of perfect read counts and 1&2 LD read counts, we estimate the number
        of valid reads and the true error rate.
        True error rate will be calculated by the recorded error rate * error rate multiplier.

        Let Z_0T = the true number of perfect reads
        Z_0E = the estimated number of perfect reads based on the estimated number of valid reads and estimated error rate
        Z_1T = the true number of 1 LD raeds
        Z_1E = the estimated number of 1 LD reads based on the estimated number of valid reads and estimated error rate
        Z_2T = the true number of 2 LD raeds in the second pass
        Z_2E = the estimated number of 2 LD reads based on the estimated number of valid reads and estimated error rate
        Z = sum(abs(Z_iT - Z_iE)  for i in range(4))

        X = estimated number of valid reads

        Y_E = vector of estimated number of error for each base
        Y_T = True error count for each base up to the second pass of the algorithm
        Y = Y_E/Y_T

        Since Z_iE = f_i(X,Y), Z = f(X,Y).
        We perform stochastic gradient descent to minimize Z, which gives us the estimated number of valid reads (X) and
        estimated final error rate (Y_E = Y * Y_T) as a result.
        """
        # Z
        Z_0T = num_reads_without_errors
        Z_1T = num_reads_with_1_error
        Z_2T = num_reads_with_2_error
        Z_1w = Z_0T / Z_1T
        Z_2w = Z_0T / Z_2T
        # X
        min_X = Z_0T + Z_1T + Z_2T
        max_X = est_total_read_count
        X = np.random.uniform(min_X, max_X)
        delta_X = (max_X - min_X) / 1000
        if delta_X < 1:
            delta_X = 1
        h_x = delta_X
        # Y
        Y = 1.5
        Y_T = error_count
        delta_Y = 0.0005 / np.max(Y_T)
        h_y = delta_Y

        # Sanity check
        if max_X < min_X:
            logging.warning(
                "Number of called cell label is greater than the number of total reads in estimate_num_valid_reads.")
            return X, Y

        i = 0
        previous_Zs = collections.deque(np.zeros(20))
        Y_E = Y_T * Y
        Z_0E = self.__get_perfect_read_count(Y_E, X, num_base)
        Z_1E = self.__get_1_err_read_count(Y_E, X, num_base)
        Z_2E = self.__get_2_err_read_count(Y_E, X, num_base, ins)
        Z = abs(Z_0T - Z_0E) + Z_1w * abs(Z_1T - Z_1E) + Z_2w * abs(Z_2T - Z_2E)
        lowest_Z = Z
        lowest_X = X
        lowest_Y = Y
        restart_count = 0
        max_iteration = 300
        delta_X_raw = delta_X
        delta_Y_raw = delta_Y
        h_x_raw = h_x
        h_y_raw = h_y
        while abs(1 - np.mean(previous_Zs) / Z) > 0.0001 and restart_count < 10:
            # update parameters
            increment_size_adjustment = np.power(1 - (i / max_iteration),0.75)
            i += 1
            delta_X = delta_X_raw * increment_size_adjustment
            h_x = h_x_raw * increment_size_adjustment
            delta_Y = delta_Y_raw * increment_size_adjustment
            h_y = h_y_raw * increment_size_adjustment

            # dZ/dY
            dY_E = Y_T * (Y + h_y)
            Z_0E = self.__get_perfect_read_count(dY_E, X, num_base)
            Z_1E = self.__get_1_err_read_count(dY_E, X, num_base)
            Z_2E = self.__get_2_err_read_count(dY_E, X, num_base, ins)
            dZ = abs(Z_0T - Z_0E) + Z_1w * abs(Z_1T - Z_1E) + Z_2w * abs(Z_2T - Z_2E) - Z
            dZ_dY = dZ / h_y
            y_step = dZ_dY * delta_Y
            y_step_stdev = abs(y_step)

            # dZ/dX
            Y_E = Y_T * Y
            dX = X + h_x
            Z_0E = self.__get_perfect_read_count(Y_E, dX, num_base)
            Z_1E = self.__get_1_err_read_count(Y_E, dX, num_base)
            Z_2E = self.__get_2_err_read_count(Y_E, dX, num_base, ins)
            dZ = abs(Z_0T - Z_0E) + Z_1w * abs(Z_1T - Z_1E) + Z_2w * abs(Z_2T - Z_2E) - Z
            dZ_dX = dZ / h_x
            x_step = dZ_dX * delta_X
            x_step_stdev = abs(x_step)

            X -= np.random.normal(x_step, x_step_stdev)
            Y -= np.random.normal(y_step, y_step_stdev)

            # Sanity checks
            if X < min_X:
                X = np.random.uniform(min_X, max_X)
                Y = np.random.uniform(1.2, 3)
                restart_count += 1
                i = 0
            elif X > max_X:
                X = np.random.uniform(min_X, max_X)
                Y = np.random.uniform(1.2, 3)
                restart_count += 1
                i = 0
            if Y < 0.8:
                X = np.random.uniform(min_X, max_X)
                Y = np.random.uniform(1.2, 3)
                restart_count += 1
                i = 0
            if i == max_iteration:
                X = np.random.uniform(min_X, max_X)
                Y = np.random.uniform(1.2, 3)
                restart_count += 1
                i = 0

            dY_E = Y_T * Y
            Z_0E = self.__get_perfect_read_count(dY_E, X, num_base)
            Z_1E = self.__get_1_err_read_count(dY_E, X, num_base)
            Z_2E = self.__get_2_err_read_count(dY_E, X, num_base, ins)
            Z = abs(Z_0T - Z_0E) + Z_1w * abs(Z_1T - Z_1E) + Z_2w * abs(Z_2T - Z_2E)
            _ = previous_Zs.popleft()
            previous_Zs.append(Z)
            if Z < lowest_Z:
                lowest_X = X
                lowest_Y = Y
                lowest_Z = Z

        return lowest_X, lowest_Y

    def update_log_prob_tables(self):
        """
        Using the estimated final error rate, calculate the probability of a perfect read/part
        Also, calculate the error rate update tables for speedy calculation of the probability of an error'ed read
        """
        # Estimate the final error rate
        self.estimate_final_error_rate()

        # Position and length of each part on read 1
        part1_start_ends = self.AnnoR1.hash_tables.part1_ck_start_ends
        part2_start_ends = self.AnnoR1.hash_tables.part2_ck_start_ends
        part3_start_ends = self.AnnoR1.hash_tables.part3_ck_start_ends
        part2_len = part2_start_ends[0][1] - part2_start_ends[0][0]
        part3_len = part3_start_ends[0][1] - part3_start_ends[0][0]

        # Initialize part 2 and 3 error rate update table
        self.__log_prob_error_rate_update_table_part2 = np.zeros((len(self.AnnoR1.bead.linker_sets), 3, part2_len))
        self.__log_prob_error_rate_update_table_part3 = np.zeros((len(self.AnnoR1.bead.linker_sets), 3, part3_len))

        for ins in self.AnnoR1.misc_functions.all_ins_lnk_combinations:
            for lnk in self.AnnoR1.misc_functions.all_ins_lnk_combinations[ins]:
                # Number of bases in this insert/linker combination
                num_bases = len(self.AnnoR1.bead.inserts[ins]) + self.AnnoR1.bead.CL3_Start_End[1]
                # Log number of reads without any error on each base
                errorless_log_counts = np.log([self.estimated_cellular_reads_by_inslnk[ins][lnk] - np.sum(
                    self.__estimated_error_count[ins, lnk, :, base]) for base in range(num_bases)])
                # Perfect read probability denominator
                self.log_prob_perfect_read[ins, lnk] = - num_bases * np.log(
                    self.estimated_cellular_reads_by_inslnk[ins][lnk])
                part1_len = part1_start_ends[ins][1] - part1_start_ends[ins][0]
                self.log_prob_perfect_part1[ins, lnk] = - part1_len * np.log(
                    self.estimated_cellular_reads_by_inslnk[ins][lnk])
                log_prob_perfect_part2 = - part2_len * np.log(self.estimated_cellular_reads_by_inslnk[ins][lnk])
                log_prob_perfect_part3 = - part3_len * np.log(self.estimated_cellular_reads_by_inslnk[ins][lnk])

                # Incorporate perfect read probability numerators onto their respective denominators
                self.log_prob_perfect_read[ins, lnk] += np.sum(errorless_log_counts)
                self.log_prob_perfect_part1[ins, lnk] += np.sum(
                    errorless_log_counts[part1_start_ends[ins][0]:part1_start_ends[ins][1]])
                log_prob_perfect_part2 += np.sum(
                    errorless_log_counts[part2_start_ends[ins][0]:part2_start_ends[ins][1]])
                log_prob_perfect_part3 += np.sum(
                    errorless_log_counts[part3_start_ends[ins][0]:part3_start_ends[ins][1]])
                if self.log_prob_perfect_part2[lnk] == 0 or self.log_prob_perfect_part2[lnk] < log_prob_perfect_part2:
                    self.log_prob_perfect_part2[lnk] = log_prob_perfect_part2
                if self.log_prob_perfect_part3[lnk] == 0 or self.log_prob_perfect_part3[lnk] < log_prob_perfect_part3:
                    self.log_prob_perfect_part3[lnk] = log_prob_perfect_part3

                for err_type in range(3):
                    # Construct error rate update tables
                    for base in range(num_bases):
                        self.__log_prob_error_rate_update_table[ins, lnk, err_type, base] = \
                            np.log(self.__estimated_error_count[ins, lnk, err_type, base]) - errorless_log_counts[base]

                    # Use the most lenient error rates (across all insert sizes)
                    # for part 2 and part 3 error rate update tables.
                    for base in range(part2_start_ends[ins][0], part2_start_ends[ins][1]):
                        base_on_part2 = base - part2_start_ends[ins][0]
                        error_rate_from_table = self.__log_prob_error_rate_update_table[ins, lnk, err_type, base]
                        error_rate_from_part2 = self.__log_prob_error_rate_update_table_part2[
                            lnk, err_type, base_on_part2]
                        if error_rate_from_table > error_rate_from_part2 or error_rate_from_part2 == 0:
                            self.__log_prob_error_rate_update_table_part2[
                                lnk, err_type, base_on_part2] = error_rate_from_table

                    for base in range(part3_start_ends[ins][0], part3_start_ends[ins][1]):
                        base_on_part3 = base - part3_start_ends[ins][0]
                        error_rate_from_table = self.__log_prob_error_rate_update_table[ins, lnk, err_type, base]
                        error_rate_from_part3 = self.__log_prob_error_rate_update_table_part3[
                            lnk, err_type, base_on_part3]
                        if error_rate_from_table > error_rate_from_part3 or error_rate_from_part3 == 0:
                            self.__log_prob_error_rate_update_table_part3[
                                lnk, err_type, base_on_part3] = error_rate_from_table

    def update_cl_prob(self, cl_prob, editops, ins, lnk):
        """
        This function updates the current cell label probability (cl_prob) using the error rate update table,
        based on the errors in editops.
        """
        for editop in editops:
            cl_prob += self.__log_prob_error_rate_update_table[
                ins, lnk, self.mutation_type_to_idx[editop[0]], editop[1]]
        return cl_prob

    def get_part1_prob(self, editops, ins, lnk):
        """This function gets part 1's probability based on editops"""
        part1_prob = self.log_prob_perfect_part1[ins][lnk]
        for editop in editops:
            part1_prob += self.__log_prob_error_rate_update_table[
                ins, lnk, self.mutation_type_to_idx[editop[0]], editop[1]]
        return part1_prob

    def get_part2_prob(self, editops, lnk):
        """This function gets part 2's probability based on editops"""
        part2_prob = self.log_prob_perfect_part2[lnk]
        for editop in editops:
            part2_prob += self.__log_prob_error_rate_update_table_part2[
                lnk, self.mutation_type_to_idx[editop[0]], editop[1]]
        return part2_prob

    def get_part3_prob(self, editops, lnk):
        """This function gets part 3's probability based on editops"""
        part3_prob = self.log_prob_perfect_part3[lnk]
        for editop in editops:
            part3_prob += self.__log_prob_error_rate_update_table_part2[
                lnk, self.mutation_type_to_idx[editop[0]], editop[1]]
        return part3_prob

    def get_error_count_table(self):
        return self.__error_count_table

    def update_error_count_table(self, insert, linker_set, editops):
        for editop in editops:
            self.__error_count_table[insert][linker_set][self.mutation_type_to_idx[editop[0]], editop[1]] += 1

    def output_error_count_tables(self):
        output_er_table = "{}_R1_error_count_table.npy".format(self.AnnoR1.library)
        np.save(output_er_table, self.__error_count_table)
