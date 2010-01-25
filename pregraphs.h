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

#else

#define DEBUGPPGRAPHPRINT(ppgraph)

#endif

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

    //the underlying graph in nauty format
    graph ulgraph[MAXN*MAXM];
};

typedef struct _primpregraph PRIMPREGRAPH;

struct _pregraph {
    int order;

    PRIMPREGRAPH *ppgraph;

    set *semiEdgeVertices;
};

typedef struct _pregraph PREGRAPH;

typedef int VERTEXPAIR[2];

#define WORKSIZE 50 * MAXM

/******************Global Variables**********************/
int vertexCount;
int currentVertexCount;
int minVertexCount;
int maxVertexCount;
long structureCount;

char outputType = 'n'; //defaults to no output
char *outputFile = NULL; //NULL == standard out

short endian = LITTLE_ENDIAN; // defines which endian should be used while exporting pregraph code

boolean allowLoops = FALSE; /* TRUE if loops are allowed*/
boolean allowMultiEdges = FALSE; /* TRUE if multi-edges are allowed*/
boolean allowSemiEdges = FALSE; /* TRUE if semi-edges are allowed*/

permutation automorphismGroupGenerators[MAXN][MAXN];
int numberOfGenerators;

/* Variables for nauty */
int nautyLabelling[MAXN], nautyPtn[MAXN];
static DEFAULTOPTIONS_GRAPH(nautyOptions);
statsblk nautyStats;
setword workspace[50 * MAXM];
graph canonicalGraph[MAXN * MAXM];

/******************Methods*******************************/
        
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
void determine_possible_sets_of_degree1_vertices(set *tempSet, set *vertexSetList, int* currentListPosition, int maximumSetSize, int currentSetSize, int currentSetElement, PRIMPREGRAPH *ppgraph);

void handle_pregraph_result(PREGRAPH *pregraph);
void handle_primpregraph_result(PRIMPREGRAPH *ppgraph);
void handle_deg1_operation_result(PRIMPREGRAPH *ppgraph);
void handle_deg2_operation_result(PRIMPREGRAPH *ppgraph,
        VERTEXPAIR *multiEdgeList, int multiEdgeListSize, int *multiEdgeOrbits, int multiEdgeOrbitCount);
void handle_deg1_operation1(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators);
void handle_deg1_operation2(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators);
void handle_deg2_operation1(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators);
void handle_deg2_operation2(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators,
        VERTEXPAIR **oldMultiEdgeList, int *oldMultiEdgeListSize, int **oldMultiEdgeOrbits, int *oldMultiEdgeOrbitCount);
void handle_deg2_operation3(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators,
        VERTEXPAIR **oldMultiEdgeList, int *oldMultiEdgeListSize, int **oldMultiEdgeOrbits, int *oldMultiEdgeOrbitCount);
void do_deg1_operations(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators);
void do_deg2_operations(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators,
        VERTEXPAIR *multiEdgeList, int multiEdgeListSize, int *multiEdgeOrbits, int multiEdgeOrbitCount);
void grow(PRIMPREGRAPH *ppgraph);
void handle_3_regular_result(graph *g);
void start();

void writeFatK2();

void construct_K2(PRIMPREGRAPH *ppgraph);
void construct_C4(PRIMPREGRAPH *ppgraph);
void construct_K3_with_spike(PRIMPREGRAPH *ppgraph);


void saveGenerators(int count, permutation perm[], nvector orbits[],
        int numorbits, int stabvertex, int n);
void copyGenerators(permutation (*copy)[MAXN][MAXN], int n);
#endif	/* _PREGRAPHS_H */

