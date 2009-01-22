import glob
try:
    import set
except ImportError:
    from sets import Set as set

extns = ['txt', 'par', 'ff', 'prm']

def get_data_files():
    # All .txt, .par, .prm and .ff files in ../../../data
    datafiles = []
    for ext in extns:
        files = [x.split("\\")[-1].lower() for x in
                 glob.glob("..\\..\\data\\*." + ext)]
        datafiles.extend(files)
    datafiles.remove("CMakeLists.txt".lower())
    return datafiles

def get_installed_files(script):
    installedfiles = []
    for line in open(script, "r"):
        if line.lstrip().startswith("File "):
            if line.find("/oname")>=0:
                temp = line.split("=")[1].split(" ")[0]
            else:
                temp = line.rstrip().split()[-1].split("\\")[-1]
            installedfiles.append(temp.lower())
    installedfiles.append("Uninstall.exe".lower())
    installedfiles.sort()
    return installedfiles

def get_uninstalled_files(script):
    uninstalledfiles = []
    for line in open(script, "r"):
        if line.lstrip().startswith("Delete \"$INSTDIR"):
            temp = line.rstrip().split("\\")[-1][:-1]
            uninstalledfiles.append(temp.lower())
    uninstalledfiles.sort()
    return uninstalledfiles

if __name__ == "__main__":
    script = "NSISScriptToCreateInstallerOBF.nsi"
##    script = "CreateInstallerForLanguages.nsi"    
    installedfiles = set(get_installed_files(script))
    uninstalledfiles = set(get_uninstalled_files(script))
    datafiles = set(get_data_files())
    print
    print "%d files are installed, and %d files are uninstalled." % (
           len(installedfiles), len(uninstalledfiles))
    print
    print "The following files are installed but not uninstalled:", installedfiles - uninstalledfiles
    print
    print "The following files are uninstalled but were never installed:", uninstalledfiles - installedfiles
    print
    print "The following files should be installed...but are not:", datafiles - installedfiles
    print
    installeddatafiles = set([x for x in installedfiles if x.split(".")[-1] in extns])
    print "The following .txt, .par, .prm or .ff files are installed but are not datafiles (they may still be legit):", installeddatafiles - datafiles
    raw_input()
