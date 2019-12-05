#!/usr/bin/env python3

'''Go through all the sdfs in the current directory and output
duplicates (so they can be removed).  Looks for exact duplicates.'''

import sys,glob

seen = set()
for fname in glob.glob('*.sdf'):
  contents = open(fname).read()
  if contents in seen:
    print(fname)
  else:
    seen.add(contents)
