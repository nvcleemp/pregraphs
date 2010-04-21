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

#else

#define DEBUGPPGRAPHPRINT(ppgraph)

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

#define NAUTY_WORKSIZE 50 * MAXM

unsigned int marksArray[MAXN];
unsigned int markValue = UINT_MAX;

#define MARKED(i) (marksArray[i] == markValue)
#define SET_MARK(i) (marksArray[i] = markValue)
#define RESET_MARK {if(markValue == UINT_MAX){int i;for(i=0;i<MAXN;i++){marksArray[i]=0;};markValue=1;}else{markValue++;}}

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

/* Variables for nauty */
int nautyLabelling[MAXN], nautyPtn[MAXN];
static DEFAULTOPTIONS_GRAPH(nautyOptions);
statsblk nautyStats;
setword nautyWorkspace[50 * MAXM];
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
void handle_3_regular_result(int *adjacencyList);
void start();
void startFromFile(FILE *inputFile);

char read_old_or_new(FILE *f, boolean bignum, int endian, unsigned short *number);
char read_2byte_number(FILE *f, unsigned short *n, int endian);
char readPregraphCode(FILE *f, PRIMPREGRAPH *ppgraph, int endian);

void writeFatK2();

void construct_K2(PRIMPREGRAPH *ppgraph);
void construct_C4(PRIMPREGRAPH *ppgraph);
void construct_K3_with_spike(PRIMPREGRAPH *ppgraph);
void construct_K3_3(PRIMPREGRAPH *ppgraph);

void saveGenerators(int count, permutation perm[], nvector orbits[],
        int numorbits, int stabvertex, int n);

void init_irreducible_graphs(int n);

void initInfo();
void logInfo(PREGRAPH *pregraph);
void printInfo();
#endif	/* _PREGRAPHS_H */

