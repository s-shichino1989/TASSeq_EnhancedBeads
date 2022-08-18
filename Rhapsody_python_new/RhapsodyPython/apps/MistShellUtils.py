# distutils: extra_compile_args = ["-O3"]
# cython: language_level=3
import subprocess
import logging


def shell_command_log_stderr(cmd, **kwargs):
    """run a shell command using the subprocess.call api, but redirect output to log"""
    outLog = kwargs.get('outLog', True)
    if "outLog" in kwargs.keys():
        del kwargs["outLog"]
    outReturn = kwargs.get('outReturn', False)
    if "outReturn" in kwargs.keys():
        del kwargs["outReturn"]
    logging.debug('Running the following: `{}`'.format(cmd))
    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        **kwargs
    )
    stdoutLines = []
    for line in iter(p.stdout.readline, ''):
        stdoutLines.append(line.strip())
    p.stdout.close()

    if outLog:
        logging.debug("\n".join(stdoutLines))

    return_code = p.wait()
    if return_code != 0:
        logging.debug("An error occurred:\n" + "\n".join(stdoutLines))
        raise subprocess.CalledProcessError(returncode=return_code, cmd=cmd)
    if outReturn:
        return stdoutLines
