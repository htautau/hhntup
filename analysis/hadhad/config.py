import os, errno
import getpass


PLOTS_DIR = './plots'


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise


def plots_dir(script):

    script = os.path.basename(script)
    script = os.path.splitext(script)[0]
    dir = os.path.join(PLOTS_DIR, script, getpass.getuser())
    if not os.path.exists(dir):
        mkdir_p(dir)
    return dir
