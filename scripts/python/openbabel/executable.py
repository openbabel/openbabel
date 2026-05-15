# -*- coding: utf-8 -*-
import os
import subprocess
import sys


ROOT_DIR = os.path.dirname(__file__)


def _program(name, args):
    return subprocess.call([os.path.join(ROOT_DIR, "bin", name)] + args, close_fds=False)


def obabel():
    suffix = '.exe' if os.name == 'nt' else ''
    raise SystemExit(_program('obabel' + suffix, sys.argv[1:]))
