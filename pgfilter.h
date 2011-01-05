/*
 * File:   pgfilter.h
 * Author: nvcleemp
 *
 * Created on January 5, 2011, 9:52 AM
 */

#ifndef _PGFILTER_H
#define	_PGFILTER_H

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


#endif	/* _PGFILTER_H */

