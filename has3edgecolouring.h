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

