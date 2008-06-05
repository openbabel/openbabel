import sys
if sys.platform.find("linux") != -1:
    try:
        import dl
    except ImportError:
        import DLFCN as dl
    sys.setdlopenflags(sys.getdlopenflags() | dl.RTLD_GLOBAL)

import obcore
import obconversion
import obtemplate
