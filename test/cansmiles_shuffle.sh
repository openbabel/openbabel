#!/bin/sh
#======================================================================
# FILE:		cansmiles_shuffle.sh
# AUTHOR:	Craig A. James
# DESCRIPTION:
#	Tests SMILES canonicalization.  For each SMILES of input:
#
#	  * Writes the SMILES (repeats it) N times
#
#	  * Runs the N repeated SMILES through "obabel ... -o smi -xC", (the
#	    "anti-canonicalizer" option that generates randomly-ordered
#	    SMILES) to generate N different but equivalent SMILES.
#
#	  * The N random SMILES are canonicalized using "obabel ... -o can",
#	    which in theory should result in N identical SMILES.
#
#	  * The N canonical smiles are run through "sort -u" (sort-unique),
#	    which should produce a file with exactly one SMILES.
#
#	  * The final output file's lines are counted.  If there is more
#	    than one line (more than one SMILES), it is reported as
#	    an error.
#
#======================================================================	
# Copyright (c) 2009, Craig A. James, eMolecules Inc.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 2 of the License.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#======================================================================

niter=$1
echo "$0: $niter shuffles per molecule";

n=0
errors=0
tmpfile=/tmp/babel$$.smi
canfile=/tmp/babel$$.can

unset xpg_echo

while read -r smiles ; do

  # Generate N random SMILES for this molecule
  n=`expr $n + 1`
  (
   i=0
   while test $i -lt $niter ; do
      /bin/echo $smiles
      i=`expr $i + 1`
   done
  ) | obabel -i smi -o smi -xC >$tmpfile 2>/dev/null 

  # Canonicalize them: They should all come out the same
  obabel $tmpfile -o can 2>/dev/null | sort -u >$canfile

  # We should only have one SMILES now.  Any more (or zero) is an error.
  count=`wc --lines <$canfile`
  if [ "$count" -ne "1" ] ; then
    echo ""
    echo "error: Record $n: Got $count SMILES for one molecule:"
    cat $canfile
    errors=`expr $errors + 1`
  fi

  echo -n "."
done

echo ""
echo "$n SMILES tested, $errors errors."
rm $tmpfile $canfile

