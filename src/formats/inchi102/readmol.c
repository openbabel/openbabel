/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.02
 * October 31, 2008
 * Developed at NIST
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-license.php
 */


#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>

#include "mode.h"

#include "incomdef.h"   /*^^^ comdef.h renamed to incomdef.h (avoid collision with MS comdef.h) */
#include "util.h"
#include "readmol.h"
#include "ichi_io.h"
#include "ichicomp.h"
#include "inpdef.h"


#include "lreadmol.h"
