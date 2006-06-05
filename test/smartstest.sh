#!/bin/sh
#
# This is just a small shell wrapper to be sure that we find the correct libs
# since 'make check' is likely run before 'make install'
#
unset BABEL_LIBDIR
BABEL_LIBDIR="../src/formats/.libs"
export BABEL_LIBDIR

./smartstest
