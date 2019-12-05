#!/usr/bin/env python
import sys, wget, os
for line in open(sys.argv[1]):
  pdb,lig = line.rstrip().split()
  print(pdb,lig)
  fname = wget.download('https://www.rcsb.org/pdb/download/downloadLigandFiles.do?ligandIdList=%s&structIdList=%s&instanceType=all&excludeUnobserved=false&includeHydrogens=false'%(lig,pdb),bar=None)
  os.rename(fname, '%s_%s.sdf'%(pdb,lig)) 
