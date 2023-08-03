/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.04
 * September 9, 2011
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC).
 * Originally developed at NIST. Modifications and additions by IUPAC 
 * and the InChI Trust.
 *
 * IUPAC/InChI-Trust Licence for the International Chemical Identifier (InChI) 
 * Software version 1.0.
 * Copyright (C) IUPAC and InChI Trust Limited
 * 
 * This library is free software; you can redistribute it and/or modify it under the 
 * terms of the IUPAC/InChI Trust Licence for the International Chemical Identifier 
 * (InChI) Software version 1.0; either version 1.0 of the License, or 
 * (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 * See the IUPAC/InChI Trust Licence for the International Chemical Identifier (InChI) 
 * Software version 1.0 for more details.
 * 
 * You should have received a copy of the IUPAC/InChI Trust Licence for the 
 * International Chemical Identifier (InChI) Software version 1.0 along with 
 * this library; if not, write to:
 * 
 * The InChI Trust
 * c/o FIZ CHEMIE Berlin
 * Franklinstrasse 11
 * 10587 Berlin
 * GERMANY
 * 
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

#ifndef COMPILE_ANSI_ONLY
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

#ifdef TARGET_LIB_FOR_WINCHI
#include "ichi_lib.h"
#endif

#include "ichicomp.h"

#if ( ADD_CMLPP == 1 )
#include "readcml.hpp"
#include "debug.h"
#endif

#include "ichi_io.h"
#include "ichiparm.h"
