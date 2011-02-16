/* has3edgecolouring.h
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
 * File:   has3edgecolouring.h
 * Author: nvcleemp
 *
 * Created on August 10, 2010, 11:20 AM
 */

#ifndef _HAS3EDGECOLOURING_H
#define	_HAS3EDGECOLOURING_H

/******************Includes**********************/

#ifndef MAXN
#define MAXN 64
#endif

#include "util.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>

/******************Defines*******************************/

#define TRUE 1
#define FALSE 0

typedef int boolean;

struct _3regpregraph {
    int order;

    int adjList[MAXN][3];

    //unsigned long adjMatrix[MAXN];
    //bit vector for each vertex

    int degree[MAXN];
};

typedef struct _3regpregraph PREGRAPH;

/******************Global Variables**********************/

//FILE *outputFile = NULL; //NULL == standard out

short defaultEndian = LITTLE_ENDIAN; // defines which endian should be used while exporting pregraph code

boolean debugOutput = FALSE;

/******************Methods*******************************/

/* Store the colouring of the current graph
 * This goes to 4 for efficiency reason, only used until 3 at most
 */
int colours[MAXN][4];

int coloursAroundVertex[MAXN];

int neighbourToIndexMapping[MAXN][MAXN];

int coloursCopy[MAXN][4];
int coloursAroundVertexCopy[MAXN];
int neighbourToIndexMappingCopy[MAXN][MAXN];

#define MAXVAL INT_MAX - 1
static int markvalue_edges = MAXVAL;
unsigned int marks_edges[MAXN][4]; //goes to 4 for efficiency reason, only used until 3 at most
#define RESETMARKS_EDGES {int mki, mkj; if ((markvalue_edges += 1) > MAXVAL) \
      { markvalue_edges = 1; for(mki=0;mki<MAXN;++mki) for(mkj=0;mkj<3;++mkj) marks_edges[mki][mkj]=0;}}
#define MARK_EDGES(v, w) marks_edges[v][w] = markvalue_edges
#define UNMARK_EDGES(v, w) marks_edges[v][w] = markvalue_edges - 1
#define ISMARKED_EDGES(v, w) (marks_edges[v][w] == markvalue_edges)

#endif	/* _HAS3EDGECOLOURING_H */

