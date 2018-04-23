#!/bin/bash

#test that an informative error message reporting needed environment variable
#is produced if plugins can't be found
export BABEL_LIBDIR=
echo "CCC" | obabel -ismi -ocan |& grep BABEL_LIBDIR >& /dev/null
if [ $? -ne 0 ]
then
 echo "FAIL"
fi
