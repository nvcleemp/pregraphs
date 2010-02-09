/* 
 * File:   printmulticode.h
 * Author: nvcleemp
 *
 * Created on February 8, 2010, 1:26 PM
 */

#ifndef _PRINTMULTICODE_H
#define	_PRINTMULTICODE_H

#ifndef MAXN
#define MAXN 50
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

struct _multigraph {
    int order;
    int multiEdgeCount; //TODO: necessary?

    //an adjacency list of the graph
    int adjList[3*MAXN];

    //the degrees of the vertices
    int degree[MAXN];

    //multiedge stores the target of the multiedge for vertices of degree 2
    //i.e. if degree[v]=2, then the edge (v,multiedge[v]) is a multiedge
    int multiedge[MAXN];
};

typedef struct _multigraph MULTIGRAPH;


/******************Global Variables**********************/

short defaultEndian = LITTLE_ENDIAN; // defines which endian should be used while exporting pregraph code

boolean debugOutput = FALSE;


#endif	/* _PRINTMULTICODE_H */

