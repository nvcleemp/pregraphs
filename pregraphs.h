/* pregraphs.h
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
 * File:   pregraphs.h
 * Author: nvcleemp
 *
 * Created on December 4, 2009, 2:18 PM
 */

#ifndef _PREGRAPHS_H
#define	_PREGRAPHS_H

/******************Includes**********************/

//#define MAXN 50 //needs to be defined before nauty.h is included !

#include "util.h"
#include "nauty/nauty.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <ctype.h>
#include <signal.h>

#ifdef _TEST
    #include "tests.h"
#endif

/******************Defines*******************************/

#ifdef _DEBUG

#define DEBUGPPGRAPHPRINT(ppgraph) {\
                                        fprintf(stderr, "==========================\n");\
                                        fprintf(stderr, "%s:%u %s:\n", __FILE__, __LINE__, #ppgraph);\
                                        fprintf(stderr, "Order   : %d\n", (ppgraph)->order);\
                                        fprintf(stderr, "# deg 1 : %d\n", (ppgraph)->degree1Count);\
                                        fprintf(stderr, "# multi : %d\n", (ppgraph)->multiEdgeCount);\
                                        int debugppgraphprintcounter;\
                                        for(debugppgraphprintcounter=0;debugppgraphprintcounter<(ppgraph)->order;debugppgraphprintcounter++){\
                                            fprintf(stderr, "%2d) ", debugppgraphprintcounter);\
                                            int debugppgraphprintcounter2;\
                                            for(debugppgraphprintcounter2=0;debugppgraphprintcounter2<(ppgraph)->degree[debugppgraphprintcounter];debugppgraphprintcounter2++){\
                                                fprintf(stderr, "%2d ", (ppgraph)->adjList[3*debugppgraphprintcounter + debugppgraphprintcounter2]);\
                                            }\
                                            if((ppgraph)->degree[debugppgraphprintcounter]==2){\
                                                fprintf(stderr, "[%2d]", (ppgraph)->multiedge[debugppgraphprintcounter]);\
                                            }\
                                            fprintf(stderr, "\n");\
                                        }\
                                        fprintf(stderr, "==========================\n");\
                                        fprintf(stderr, "Nauty:\n");\
                                        for(debugppgraphprintcounter=0;debugppgraphprintcounter<(ppgraph)->order;debugppgraphprintcounter++){\
                                            set *debugGraphRow;\
                                            debugGraphRow = GRAPHROW((ppgraph)->ulgraph, debugppgraphprintcounter, MAXM);\
                                            fprintf(stderr, "%2d) ", debugppgraphprintcounter);\
                                            int debugppgraphprintcounter2;\
                                            for(debugppgraphprintcounter2=-1;(debugppgraphprintcounter2 = nextelement(debugGraphRow, MAXM, debugppgraphprintcounter2)) >=0;){\
                                                fprintf(stderr, "%2d ", debugppgraphprintcounter2);\
                                            }\
                                            fprintf(stderr, "\n");\
                                        }\
                                        fprintf(stderr, "==========================\n");\
                                        fflush(stderr);\
                                    }

#define DEBUGBRIDGESPRINT(ppgraph) {\
                                        fprintf(stderr, "=========Bridges========\n");\
                                        fprintf(stderr, "%s:%u %s:\n", __FILE__, __LINE__, #ppgraph);\
                                        int debugBridgesPrintCounter;\
                                        for(debugBridgesPrintCounter=0;debugBridgesPrintCounter<(ppgraph)->bridgeCount;debugBridgesPrintCounter++){\
                                            fprintf(stderr, "%2d) ", debugBridgesPrintCounter);\
                                            fprintf(stderr, "%2d ", (ppgraph)->bridges[debugBridgesPrintCounter][0]);\
                                            fprintf(stderr, "%2d ", (ppgraph)->bridges[debugBridgesPrintCounter][1]);\
                                            fprintf(stderr, "\n");\
                                        }\
                                        fprintf(stderr, "=========Bridges========\n");\
                                   }

#else

#define DEBUGPPGRAPHPRINT(ppgraph)

#define DEBUGBRIDGESPRINT(ppgraph)

#endif

#ifdef _CONSISTCHECK

    boolean doConsistencyCheck = FALSE;

    #define ENABLECONSISTENCYCHECK { doConsistencyCheck = TRUE; }

    #define DISABLECONSISTENCYCHECK { doConsistencyCheck = FALSE; }

    #define CHECKCONSISTENCY(ppgraph) {\
                                        if(doConsistencyCheck){\
                                            int consistencyCheckI, consistencyCheckJ;\
                                            for(consistencyCheckI = 0; consistencyCheckI < ppgraph->order; consistencyCheckI++){\
                                                for(consistencyCheckJ = 0; consistencyCheckJ < ppgraph->degree[consistencyCheckI]; consistencyCheckJ++){\
                                                    if(!areAdjacent(ppgraph, ppgraph->adjList[consistencyCheckI*3 + consistencyCheckJ], consistencyCheckI)){\
                                                        fprintf(stderr, "%s:%u Consistency check failed:\n", __FILE__, __LINE__);\
                                                        fprintf(stderr, "   %d and %d are not adjacent.\n", ppgraph->adjList[consistencyCheckI*3 + consistencyCheckJ], consistencyCheckI);\
                                                        DEBUGPPGRAPHPRINT(ppgraph)\
                                                        show_stackframe();\
                                                    }\
                                                }\
                                            }\
                                        }\
                                    }

#else

    #define ENABLECONSISTENCYCHECK
    #define DISABLECONSISTENCYCHECK
    #define CHECKCONSISTENCY(ppgraph)

#endif

#ifdef _PROFILING
    #define _PROFILING_DEG1
    #define _PROFILING_OPERATION1_SAME_COLOUR
    #define _PROFILING_DEG2
#endif

#define REG 3

//For the timing
#define time_factor sysconf(_SC_CLK_TCK)

struct _primpregraph {
    int order;
    int degree1Count;
    int multiEdgeCount;

    //an adjacency list of the graph
    int adjList[3*MAXN];

    //the degrees of the vertices
    int degree[MAXN];

    //multiedge stores the target of the multiedge for vertices of degree 2
    //i.e. if degree[v]=2, then the edge (v,multiedge[v]) is a multiedge
    int multiedge[MAXN];

    //the underlying graph in nauty format
    graph ulgraph[MAXN*MAXM];

    //a list of the bridges in this graph
    //MAXN-1 is an upperbound for the number of bridges, because a spanning
    //tree contains all bridges and order - 1 edges. This upperbound is
    //sharp because a tree is a possible graph we get.
    //currently this information is only correct as long as we're working on
    //degree 1 operations. As soon as we start with degree 2 operations this
    //information is no longer updated.
    //the bridges need to be stored so that the largest vertex comes first.
    //Because operation 1.1 can remove bridges and this is currently not
    //checked, this list only contains a possible bridges and some edges
    //might actually not be a bridge
    int bridges[MAXN-1][2];

    int bridgeCount;

    int bridgePosition[MAXN][MAXN];
};

typedef struct _primpregraph PRIMPREGRAPH;

struct _pregraph {
    int order;

    PRIMPREGRAPH *ppgraph;

    set *semiEdgeVertices;
};

typedef struct _pregraph PREGRAPH;

typedef int VERTEXPAIR[2];

#define NAUTY_WORKSIZE 50 * MAXM

unsigned int marksArray[MAXN];
unsigned int markValue = UINT_MAX;

#define MARKED(i) (marksArray[i] == markValue)
#define SET_MARK(i) (marksArray[i] = markValue)
#define RESET_MARK {if(markValue == UINT_MAX){int i;for(i=0;i<MAXN;i++){marksArray[i]=0;};markValue=1;}else{markValue++;}}

//the radius of the neighbourhood that is used as colour for degree 1 vertices
#define DEG1_DISTANCE_COLOUR_VALUE 4

/******************Global Variables**********************/
int vertexCount;
int currentVertexCount;
int minVertexCount;
int maxVertexCount;
unsigned long long structureCount;
unsigned long long primitivesCount;

boolean logStatistics = FALSE;
//used for statistics
unsigned long long *graphsWithLoopsCount;
unsigned long long *graphsWithSemiEdgesCount;
unsigned long long *graphsWithMultiEdgesCount;
unsigned long long graphsWithOnlyLoopsCount;
unsigned long long graphsWithOnlySemiEdgesCount;
unsigned long long graphsWithOnlyMultiEdgesCount;
unsigned long long simplegraphsCount;

#ifdef _PROFILING_DEG1

unsigned long long degree1Operation1Total;
unsigned long long degree1Operation1Canonical;
unsigned long long degree1Operation2Total;
unsigned long long degree1Operation2Canonical;

unsigned long long canonicalDegree1Calls;
unsigned long long canonicalDegree1BecauseOnlyOneVertexOfDegree1;
unsigned long long canonicalDegree1TrivialRemainsTrivial;
unsigned long long canonicalDegree1BridgeFixed;
unsigned long long canonicalDegree1NotBecauseNotSmallestColour;
unsigned long long canonicalDegree1BecauseOnlyOneMinimumColour;
unsigned long long canonicalDegree1YesWithNauty;
unsigned long long canonicalDegree1NoWithNauty;

int canonicalDegree1PossibleColoursCount;
unsigned long long *canonicalDegree1MinimumColourFrequency;

int *canonicalDegree1PartitionCountFrequency;

unsigned long long canonicalDegree1Degree3PartitionCount[MAXN+1];

unsigned long long *canonicalDegree1Degree1PartitionSize;
unsigned long long *canonicalDegree1Degree1PartitionCount;
unsigned long long canonicalDegree1Degree3Neighbours0PartitionSize;
unsigned long long canonicalDegree1Degree3Neighbours0PartitionCount;
unsigned long long *canonicalDegree1Degree3Neighbours1PartitionSize;
unsigned long long *canonicalDegree1Degree3Neighbours1PartitionCount;
unsigned long long *canonicalDegree1Degree3Neighbours2PartitionSize;
unsigned long long *canonicalDegree1Degree3Neighbours2PartitionCount;
#endif

#ifdef _PROFILING_OPERATION1_SAME_COLOUR
unsigned long long degree1Operation1TooFewVertices = 0;
unsigned long long degree1Operation1OperationTried = 0;
unsigned long long degree1Operation1HaveSameColour = 0;
unsigned long long degree1Operation1HaveDifferentColour = 0;
unsigned long long degree1Operation1SameColourSolved = 0;
unsigned long long degree1Operation1SolvedWithPaths = 0;
unsigned long long degree1Operation1NotSolvedWithPaths = 0;
unsigned long long degree1Operation1Degree1Counts[MAXN];
#endif

int degree1OperationsDepth = 0;
int degree1OperationsDepthMaximum = 0;
int degree2OperationsDepth = 0;
int degree2OperationsDepthMaximum = 0;

VERTEXPAIR *globalMultiEdgeList;
int *globalMultiEdgeListSize;
int *globalMultiEdgeOrbits;
int *globalMultiEdgeOrbitCount;

char outputType = 'n'; //defaults to no output
char *outputFile = NULL; //NULL == standard out

short endian = LITTLE_ENDIAN; // defines which endian should be used while exporting pregraph code

boolean onlyPrimitives = FALSE;
boolean allowLoops = FALSE; /* TRUE if loops are allowed*/
boolean allowMultiEdges = FALSE; /* TRUE if multi-edges are allowed*/
boolean allowSemiEdges = FALSE; /* TRUE if semi-edges are allowed*/
boolean onlyColourable = FALSE; /* TRUE if only 3-edge-colourable pregraphs are allowed */
boolean onlyBipartite = FALSE; /* TRUE if only bipartite pregraphs are allowed */

boolean operation11Disabled = FALSE;
boolean operation12Disabled = FALSE;
boolean operation21Disabled = FALSE;
boolean operation22Disabled = FALSE;
boolean operation23Disabled = FALSE;

boolean noRejections = FALSE;

boolean moduloEnabled = FALSE;
int moduloRest;
int moduloMod;
unsigned long long int splitPointCount = 0;
int splitDepth = 0;

/* Provide space for the generators at each recursion depth (maximum depth = MAXN + 1)
 * There are at most n<=MAXN generators in a graph with n vertices and the length
 * of each generator is n.
 */
permutation automorphismGroupGenerators[MAXN+1][MAXN][MAXN];
int numberOfGenerators[MAXN+1];

/* Keep track of the bridges that are changed in each step.
 */
int changedBridges[MAXN+1][2];
int numberOfChangedBridges[MAXN+1];

/* Variables for nauty */
int nautyLabelling[MAXN], nautyPtn[MAXN];
static DEFAULTOPTIONS_GRAPH(nautyOptions);
statsblk nautyStats;
setword nautyWorkspace[50 * MAXM];
graph canonicalGraph[MAXN * MAXM];

/* Store the colouring of the current graph
 * This goes to 4 for efficiency reason, only used until 3 at most
 */
int colours[MAXN][4];

/* Store the vertex colouring of the current graph
 */
int vertexColours[MAXN];
#define BLACK 0;
#define WHITE 1;

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

/******************Methods*******************************/

char writePregraphCode(FILE *f, PREGRAPH *ppgraph);
char writePregraphTable(FILE *f, PREGRAPH *ppgraph);
char writePrimpregraphCode(FILE *f, PRIMPREGRAPH *ppgraph);
char writePrimpregraphTable(FILE *f, PRIMPREGRAPH *ppgraph);

inline boolean areAdjacent(PRIMPREGRAPH *ppgraph, int u, int v);

//these operations provide the actual implementation of the different operations
//they don't perform any checking of the validity of performing this method
void apply_deg1_operation1(PRIMPREGRAPH *ppgraph, int u, int v);
void apply_deg1_operation2(PRIMPREGRAPH *ppgraph, int u, int v);
void apply_deg2_operation1(PRIMPREGRAPH *ppgraph, int u, int v);
void apply_deg2_operation2(PRIMPREGRAPH *ppgraph, int u, int v);
void apply_deg2_operation3(PRIMPREGRAPH *ppgraph, int u, int v);

void revert_deg1_operation1(PRIMPREGRAPH *ppgraph, int u, int v);
void revert_deg1_operation2(PRIMPREGRAPH *ppgraph, int u, int v);
void revert_deg2_operation1(PRIMPREGRAPH *ppgraph, int u, int v);
void revert_deg2_operation2(PRIMPREGRAPH *ppgraph, int u, int v);
void revert_deg2_operation3(PRIMPREGRAPH *ppgraph, int u, int v);

void get_deg1_pairs(PRIMPREGRAPH *ppgraph, VERTEXPAIR *vertexPairList, int *vertexPairListSize);
void get_single_edges(PRIMPREGRAPH *ppgraph, VERTEXPAIR *vertexPairList, int *vertexPairListSize);
void get_multi_edges(PRIMPREGRAPH *ppgraph, VERTEXPAIR *vertexPairList, int *vertexPairListSize);

//union-find algorithm
void unionElements(int *forest, int *treeSizes, int *numberOfComponents, int element1, int element2);
int findRootOfElement(int *forest, int element);

int nextDegree1Vertex(int current, PRIMPREGRAPH *ppgraph);
void determine_possible_sets_of_degree1_vertices(set *tempSet, set *vertexSetList, int* currentListPosition, int maximumSetSize, int currentSetSize, int currentSetElement, PRIMPREGRAPH *ppgraph, int skippedVertices);

void handle_pregraph_result(PREGRAPH *pregraph);
void handle_primpregraph_result(PRIMPREGRAPH *ppgraph);
void handle_deg1_operation_result(PRIMPREGRAPH *ppgraph);
void handle_deg2_operation_result(PRIMPREGRAPH *ppgraph, boolean multiEdgesDetermined);
void handle_deg1_operation1(PRIMPREGRAPH *ppgraph);
void handle_deg1_operation2(PRIMPREGRAPH *ppgraph);
void handle_deg2_operation1(PRIMPREGRAPH *ppgraph);
void handle_deg2_operation2(PRIMPREGRAPH *ppgraph, boolean *multiEdgesDetermined);
void handle_deg2_operation3(PRIMPREGRAPH *ppgraph, boolean *multiEdgesDetermined);
void do_deg1_operations(PRIMPREGRAPH *ppgraph);
void do_deg2_operations(PRIMPREGRAPH *ppgraph, boolean multiEdgesDetermined);
void grow(PRIMPREGRAPH *ppgraph);
void growWithoutDeg1Operations(PRIMPREGRAPH *ppgraph);
void handle_snarkhunter_result(unsigned char snarkhunter_graph[MAXN][REG + 1], int order);
void start();
void startFromFile(FILE *inputFile);

char read_old_or_new(FILE *f, boolean bignum, int endian, unsigned short *number);
char read_2byte_number(FILE *f, unsigned short *n, int endian);
char readPregraphCode(FILE *f, PRIMPREGRAPH *ppgraph, int endian);

void writeThetaGraph();

void construct_K2(PRIMPREGRAPH *ppgraph);
void construct_C4(PRIMPREGRAPH *ppgraph);
void construct_K3_with_spike(PRIMPREGRAPH *ppgraph);
void construct_K3_3(PRIMPREGRAPH *ppgraph);

void saveGenerators(int count, permutation perm[], nvector orbits[],
        int numorbits, int stabvertex, int n);

void call_snarkhunter(int n, int g, void (*userproc) (unsigned char (*)[REG + 1], int));

void initInfo();
void logInfo(PREGRAPH *pregraph);
void printInfo();
#endif	/* _PREGRAPHS_H */

