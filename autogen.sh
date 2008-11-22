#!/bin/sh

OSTYPE=`uname -s`

AMFLAGS="--foreign --add-missing"
if test "$OSTYPE" = "IRIX" -o "$OSTYPE" = "IRIX64"; then
   AMFLAGS=$AMFLAGS" --include-deps";
fi

rm -fr autom4te.cache
echo "Running libtoolize"
libtoolize --force --copy
echo "Running aclocal"
aclocal --force
echo "Running autoheader"
autoheader --force
echo "Running automake"
automake $AMFLAGS
echo "Running autoconf"
autoconf
echo "Running autoreconf"
autoreconf

echo "======================================"
echo "Now you are ready to run './configure'"
echo "======================================"
