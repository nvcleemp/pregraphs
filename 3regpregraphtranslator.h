/* 
 * File:   3regpregraphtranslator.h
 * Author: nvcleemp
 *
 * Created on February 4, 2010, 9:51 AM
 */

#ifndef _3REGPREGRAPHSTRANSLATOR_H
#define	_3REGPREGRAPHSTRANSLATOR_H

/******************Includes**********************/

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

struct _primpregraph {
    int order;
    int degree1Count;
    int multiEdgeCount; //TODO: necessary?

    //an adjacency list of the graph
    int adjList[3*MAXN];

    //the degrees of the vertices
    int degree[MAXN];

    //multiedge stores the target of the multiedge for vertices of degree 2
    //i.e. if degree[v]=2, then the edge (v,multiedge[v]) is a multiedge
    int multiedge[MAXN];
};

typedef struct _primpregraph PRIMPREGRAPH;

struct _pregraph {
    int order;

    PRIMPREGRAPH *ppgraph;

    boolean semiEdgeVertices[MAXN];
};

typedef struct _pregraph PREGRAPH;

/******************Global Variables**********************/

//FILE *outputFile = NULL; //NULL == standard out

short defaultEndian = LITTLE_ENDIAN; // defines which endian should be used while exporting pregraph code

boolean debugOutput = FALSE;

/******************Methods*******************************/


#endif	/* _3REGPREGRAPHSTRANSLATOR_H */

