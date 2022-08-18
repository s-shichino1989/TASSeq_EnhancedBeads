# distutils: extra_compile_args = ["-O3"]
# cython: language_level=3
from MistShellUtils import shell_command_log_stderr

def trim_primer(R1, trimmed_r1):
    def cutadaptCmd(r1, trimmed_r1):
        runCmd = f"cutadapt -g ACAGGAAACTCATGGTGCGT -O 16 -o {trimmed_r1} {r1}"
        return runCmd

    runCmd = cutadaptCmd(R1, trimmed_r1)
    shell_command_log_stderr(runCmd, shell=True, outLog=False, outReturn=False)
