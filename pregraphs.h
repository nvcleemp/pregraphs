/* 
 * File:   pregraphs.h
 * Author: nvcleemp
 *
 * Created on December 4, 2009, 2:18 PM
 */

#ifndef _PREGRAPHS_H
#define	_PREGRAPHS_H

/******************Includes**********************/

#define MAXN 50 //needs to be defined before nauty.h is included !

#include "util.h"
#include "nauty/nauty.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#ifdef _TEST
    #include "tests.h"
#endif

/******************Defines*******************************/

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

    //the graph in nauty format
    graph graph;
};

typedef struct _primpregraph PRIMPREGRAPH;

struct _pregraph {
    int order;

    PRIMPREGRAPH *ppgraph;

    set *semiEdgeVertices;
};

typedef struct _pregraph PREGRAPH;

typedef int VERTEXPAIR[2];

/******************Global Variables**********************/
int vertexCount;
int currentVertexCount;
int minVertexCount;
int maxVertexCount;

boolean onlyCount = FALSE; /* TRUE if the graphs don't need to be outputted*/
boolean allowLoops = FALSE; /* TRUE if loops are allowed*/
boolean allowMultiEdges = FALSE; /* TRUE if multi-edges are allowed*/
boolean allowSemiEdges = FALSE; /* TRUE if semi-edges are allowed*/

permutation generators[MAXN][MAXN];
int number_of_generators;

/* Variables for nauty */
int lab[MAXN], ptn[MAXN];
static DEFAULTOPTIONS_GRAPH(options);
statsblk stats;
setword workspace[50 * MAXM];

/******************Methods*******************************/
        
inline boolean areAdjacent(PRIMPREGRAPH *ppgraph, int u, int v);

//these operations provide the actual implementation of the different operations
//they don't perform any checking of the validity of performing this method
void apply_deg1_operation1(PRIMPREGRAPH *ppgraph, int u, int v);
void apply_deg1_operation2(PRIMPREGRAPH *ppgraph, int u, int v);
void apply_deg2_operation1(PRIMPREGRAPH *ppgraph, int u, int v);
void apply_deg2_operation2(PRIMPREGRAPH *ppgraph, int u, int v);
void apply_deg2_operation3(PRIMPREGRAPH *ppgraph, int u, int v);


void get_deg1_pairs(PRIMPREGRAPH *ppgraph, VERTEXPAIR *vertexPairList, int *vertexPairListSize);
void get_single_edges(PRIMPREGRAPH *ppgraph, VERTEXPAIR *vertexPairList, int *vertexPairListSize);
void get_multi_edges(PRIMPREGRAPH *ppgraph, VERTEXPAIR *vertexPairList, int *vertexPairListSize);

//union-find algorithm
void union_elements(int *forest, int *treeSizes, int *numberOfComponents, int element1, int element2);
int find_root_of_element(int *forest, int element);

int nextDegree1Vertex(int current, PRIMPREGRAPH *ppgraph);
void determine_possible_sets_of_degree1_vertices(set *tempSet, set *vertexSetList, int* currentListPosition, int maximumSetSize, int currentSetSize, int currentSetElement, PRIMPREGRAPH *ppgraph);

void handle_pregraph_result(PREGRAPH *pregraph);
void handle_primpregraph_result(PRIMPREGRAPH *ppgraph);
void handle_deg1_operation_result(PRIMPREGRAPH *ppgraph);
void handle_deg1_operation1(PRIMPREGRAPH *ppgraph);
void handle_deg1_operation2(PRIMPREGRAPH *ppgraph);
void handle_deg2_operation1(PRIMPREGRAPH *ppgraph);
void handle_deg2_operation2(PRIMPREGRAPH *ppgraph);
void handle_deg2_operation3(PRIMPREGRAPH *ppgraph);
void do_deg1_operations(PRIMPREGRAPH *ppgraph);
void do_deg2_operations(PRIMPREGRAPH *ppgraph);
void grow(PRIMPREGRAPH *ppgraph);
void start();

void construct_K2(PRIMPREGRAPH *ppgraph);
void construct_C4(PRIMPREGRAPH *ppgraph);


void save_generators(int count, permutation perm[], nvector orbits[],
        int numorbits, int stabvertex, int n);
#endif	/* _PREGRAPHS_H */

