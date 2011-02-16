/* admissable_c.h
 * =========================================================================
 * This file is part of the pregraphs project
 *
 * Copyright (C) 2010-2011 Universiteit Gent
 *
 * Author: Nicolas Van Cleemput
 * In collaboration with Gunnar Brinkmann

 *  * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * A copy of the GNU General Public License can be found in the file
 * LICENSE.txt provided with this distribution. This license can also
 * be found on the GNU website at http://www.gnu.org/licenses/gpl.html.
 *
 * If you did not receive a copy of the GNU General Public License along
 * with this program, contact the lead developer, or write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 */

/*
 * File:   admissable_c.h
 * Author: nvcleemp
 *
 * Created on June 9, 2010, 9:52 AM
 */

#ifndef _ADMISSABLE_C_H
#define	_ADMISSABLE_C_H

/******************Includes**********************/

#ifndef MAXN
#define MAXN 64
#endif

#include "util.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <ctype.h>
#include <string.h>

/******************Defines*******************************/

#define TRUE 1
#define FALSE 0

typedef int boolean;

struct _3regpregraph {
    int order;

    int adjList[MAXN][3];

    unsigned long adjMatrix[MAXN];
    //bit vector for each vertex

    int degree[MAXN];
};

typedef struct _3regpregraph PREGRAPH;

/******************Global Variables**********************/

//FILE *outputFile = NULL; //NULL == standard out

short defaultEndian = LITTLE_ENDIAN; // defines which endian should be used while exporting pregraph code

boolean debugOutput = FALSE;

/******************Methods*******************************/


#endif	/* _ADMISSABLE_C_H */

