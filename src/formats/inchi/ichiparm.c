/*
 * International Union of Pure and Applied Chemistry (IUPAC)
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.01
 * July 21, 2006
 * Developed at NIST
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
/* #include <varargs.h> */
#include <errno.h>
#include <limits.h>

#include "mode.h"       /* moved from below, suggestion by David Mosenkis */

#ifndef INCHI_ANSI_ONLY
#include <conio.h>
#endif

#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichidrp.h"
#include "ichierr.h"
#include "ichimain.h"
#include "extr_ct.h"

#ifdef INCHI_LIB
#include "ichi_lib.h"
#endif

#include "ichicomp.h"

#if( ADD_CMLPP == 1 )
#include "readcml.hpp"
#include "debug.h"
#endif

#include "ichiparm.h"
