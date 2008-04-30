#!/bin/sh

OSTYPE=`uname -s`

AMFLAGS="--add-missing"
if test "$OSTYPE" = "IRIX" -o "$OSTYPE" = "IRIX64"; then
   AMFLAGS=$AMFLAGS" --include-deps";
fi

echo "Running aclocal$AMVER"
rm -fr autom4te.cache
aclocal
echo "Running libtoolize"
libtoolize --force --copy
echo "Running autoheader"
autoheader
echo "Running automake$AMVER"
automake $AMFLAGS
echo "Running autoconf"
autoconf

echo "======================================"
echo "Now you are ready to run './configure'"
echo "======================================"
