"""
Set BABEL_DATADIR to the location where the
openbabel_data was installed
"""

import os, sys
import _winreg

if len(sys.argv)==2 and sys.argv[1]=="-install":
    # Connect to the registry
    registry = _winreg.ConnectRegistry(None,_winreg.HKEY_LOCAL_MACHINE)
    # Open Environment key for writing
    environment =_winreg.OpenKey(registry,
       r"SYSTEM\CurrentControlSet\Control\Session Manager\Environment",
       0,_winreg.KEY_ALL_ACCESS)
    # Set the value of BABEL_DATADIR
    datadir = os.path.join(sys.prefix, "Lib", "site-packages", "openbabel_data")
    _winreg.SetValueEx(environment, "BABEL_DATADIR", 0, _winreg.REG_EXPAND_SZ,
       datadir)
    _winreg.CloseKey(environment)
    _winreg.CloseKey(registry)

    print "BABEL_DATADIR is set to %s" % datadir
    print
    print "You will need to reboot before the openbabel module"
    print "can access the new value of BABEL_DATADIR. However,"
    print "you can start using the module right away if you wish."
