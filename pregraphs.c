/*
 * File:   pregraphs.c
 * Author: nvcleemp
 *
 * Created on December 8, 2009, 2:32 PM
 */

//#define _TEST
//#define _DEBUG
//#define _CONSISTCHECK

#include "pregraphs.h"

#include <stdio.h>
#include <signal.h>
#include <execinfo.h>

void show_stackframe() {
  void *trace[16];
  char **messages = (char **)NULL;
  int i, trace_size = 0;

  trace_size = backtrace(trace, 16);
  messages = backtrace_symbols(trace, trace_size);
  fprintf(stderr, "[bt] Execution path:\n");
  for (i=0; i<trace_size; ++i){
	fprintf(stderr, "[bt] %s\n", messages[i]);
  }
}

inline boolean areAdjacent(PRIMPREGRAPH *ppgraph, int u, int v){
    int i;
    for (i = 0; i < ppgraph->degree[u]; i++) {
        if(ppgraph->adjList[u*3+i]==v)
            return TRUE;
    }
    return FALSE;
}

/*
 * Performs a depth-first search of the graph looking for the vertex find.
 * Returns TRUE if the vertex was reached and FALSE if it wasn't reached.
 */
boolean DFSearch(PRIMPREGRAPH *ppgraph, int current, int find, set *visited){
    int i;
    for (i = 0; i < ppgraph->degree[current]; i++) {
        if(!ISELEMENT(visited, ppgraph->adjList[current*3+i])){
            if(ppgraph->adjList[current*3+i]==find) return TRUE;

            ADDELEMENT(visited, ppgraph->adjList[current*3+i]);
            if(DFSearch(ppgraph, ppgraph->adjList[current*3+i], find, visited))
                return TRUE;
        }
    }
    return FALSE;
}

/*
 * Returns true if the edge uv is a bridge
 */
boolean isBridge(PRIMPREGRAPH *ppgraph, int u, int v){
    DEBUGMSG("Start isBridge")
    DEBUGASSERT(areAdjacent(ppgraph, u, v))

    if(ppgraph->degree[u]==1 || ppgraph->degree[v]==1){DEBUGMSG("End isBridge") return TRUE;}

    set visited[MAXM];
    EMPTYSET(visited, MAXM);
    ADDELEMENT(visited, u);
    int i;
    for (i = 0; i < ppgraph->degree[u]; i++) {
        if(ppgraph->adjList[u*3+i]!=v && !ISELEMENT(visited, ppgraph->adjList[u*3+i])){
            ADDELEMENT(visited, ppgraph->adjList[u*3+i]);
            if(DFSearch(ppgraph, ppgraph->adjList[u*3+i], v, visited)){
                DEBUGMSG("End isBridge")
                return FALSE;
            }
        }
    }
    DEBUGMSG("End isBridge")
    return TRUE;
}
//-----------------------------------------------------------------

/*                         o v
 *   u     v               |
 *   o     o               o u
 *  _|_   _|_            _/ \_
 * /   \_/   \   ===>   / \_/ \
 * \_________/          \_____/
 */
void apply_deg1_operation1(PRIMPREGRAPH *ppgraph, int u, int v){
    DEBUGMSG("Start apply_deg1_operation1")
    DEBUGDUMP(u, "%d")
    DEBUGDUMP(v, "%d")
    DEBUGASSERT(ppgraph->degree[u]==1 && ppgraph->degree[v]==1)
    int t = ppgraph->adjList[v*3]; //original neighbour of v
    ppgraph->degree[u]=3;
    ppgraph->adjList[u*3+1]=v;
    ppgraph->adjList[u*3+2]=t;
    int i=0;
    while(ppgraph->adjList[t*3+i]!=v) i++; //t and v are adjacent so will stop before i == 3
    ppgraph->adjList[t*3+i]=u;
    ppgraph->adjList[v*3]=u;
    ppgraph->degree1Count--;

    set *gu, *gv, *gt;
    gu = GRAPHROW(ppgraph->ulgraph, u, MAXM);
    gv = GRAPHROW(ppgraph->ulgraph, v, MAXM);
    gt = GRAPHROW(ppgraph->ulgraph, t, MAXM);
    ADDELEMENT(gu, v);
    DELELEMENT(gv,t);
    ADDELEMENT(gv,u);
    ADDELEMENT(gu,t);
    DELELEMENT(gt,v);
    ADDELEMENT(gt,u);
    CHECKCONSISTENCY(ppgraph)
    DEBUGMSG("End apply_deg1_operation1")
}

void revert_deg1_operation1(PRIMPREGRAPH *ppgraph, int u, int v){
    DEBUGMSG("Start revert_deg1_operation1")
    DEBUGDUMP(u, "%d")
    DEBUGDUMP(v, "%d")
    //the original neighbour of u is the first element in the adjList
    //so no need to change anything except the degree
    int t = ppgraph->adjList[u*3+2];
    ppgraph->degree[u]=1;
    ppgraph->adjList[v*3]=t;
    int i=0;
    while(ppgraph->adjList[t*3+i]!=u) i++; //u and t are adjacent so will stop before i == 3
    ppgraph->adjList[t*3+i]=v;
    ppgraph->degree1Count++;

    set *gu, *gv, *gt;
    gu = GRAPHROW(ppgraph->ulgraph, u, MAXM);
    gv = GRAPHROW(ppgraph->ulgraph, v, MAXM);
    gt = GRAPHROW(ppgraph->ulgraph, t, MAXM);
    DELELEMENT(gv, u);
    ADDELEMENT(gv, t);
    DELELEMENT(gu, v);
    DELELEMENT(gu, t);
    ADDELEMENT(gt, v);
    DELELEMENT(gt, u);
    CHECKCONSISTENCY(ppgraph)
    DEBUGMSG("End revert_deg1_operation1")
}

//-----------------------------------------------------------------

/*                        o t
 *                        |
 *                        o s
 *  _______            __/ \__
 * / u\ /v \   ===>   / u\ /v \
 * \__/ \__/          \__/ \__/
 */
void apply_deg1_operation2(PRIMPREGRAPH *ppgraph, int u, int v){
    DEBUGMSG("Start apply_deg1_operation2")
    DEBUGDUMP(u, "%d")
    DEBUGDUMP(v, "%d")
    DEBUGASSERT(areAdjacent(ppgraph, u, v))
    int s, t, i;
    s = ppgraph->order;
    t = s + 1;
    i=0;
    while(ppgraph->adjList[u*3+i]!=v) i++; //u and v are adjacent so will stop before i == 3
    ppgraph->adjList[u*3+i]=s;
    i=0;
    while(ppgraph->adjList[v*3+i]!=u) i++; //u and v are adjacent so will stop before i == 3
    ppgraph->adjList[v*3+i]=s;

    ppgraph->degree[s]=3;
    ppgraph->adjList[s*3]=u;
    ppgraph->adjList[s*3+1]=v;
    ppgraph->adjList[s*3+2]=t;

    ppgraph->degree[t]=1;
    ppgraph->adjList[t*3]=s;

    ppgraph->degree1Count++;
    ppgraph->order+=2;

    set *gu, *gv, *gs, *gt;
    gu = GRAPHROW(ppgraph->ulgraph, u, MAXM);
    gv = GRAPHROW(ppgraph->ulgraph, v, MAXM);
    gs = GRAPHROW(ppgraph->ulgraph, s, MAXM);
    gt = GRAPHROW(ppgraph->ulgraph, t, MAXM);
    EMPTYSET(gs, MAXM);
    EMPTYSET(gt, MAXM);
    DELELEMENT(gu, v);
    DELELEMENT(gv, u);
    ADDELEMENT(gu, s);
    ADDELEMENT(gv, s);
    ADDELEMENT(gs, u);
    ADDELEMENT(gs, v);
    ADDELEMENT(gs, t);
    ADDELEMENT(gt, s);
    CHECKCONSISTENCY(ppgraph)
    DEBUGMSG("End apply_deg1_operation2")
}

void revert_deg1_operation2(PRIMPREGRAPH *ppgraph, int u, int v){
    DEBUGMSG("Start revert_deg1_operation2")
    DEBUGDUMP(u, "%d")
    DEBUGDUMP(v, "%d")
    int s, i;
    s = ppgraph->order-2;
    i=0;
    while(ppgraph->adjList[u*3+i]!=s) i++;
    ppgraph->adjList[u*3+i]=v;
    i=0;
    while(ppgraph->adjList[v*3+i]!=s) i++;
    ppgraph->adjList[v*3+i]=u;

    ppgraph->degree1Count--;
    ppgraph->order-=2;

    set *gu, *gv;
    gu = GRAPHROW(ppgraph->ulgraph, u, MAXM);
    gv = GRAPHROW(ppgraph->ulgraph, v, MAXM);
    DELELEMENT(gu, s);
    DELELEMENT(gv, s);
    ADDELEMENT(gu, v);
    ADDELEMENT(gv, u);
    CHECKCONSISTENCY(ppgraph)
    DEBUGMSG("End revert_deg1_operation2")
}

//-----------------------------------------------------------------

/*
 *  __u      v__            __u  s  t  v__
 * /  o------o  \          /  o--o--o--o  \
 * |   \____/   |   ===>   |   \______/   |
 * \____________/          \______________/
 */
void apply_deg2_operation1(PRIMPREGRAPH *ppgraph, int u, int v){
    DEBUGMSG("Start apply_deg2_operation1")
    DEBUGDUMP(u, "%d")
    DEBUGDUMP(v, "%d")
    DEBUGASSERT(areAdjacent(ppgraph, u, v))
    DEBUGASSERT(ppgraph->degree[u]!=2 || ppgraph->multiedge[u]!=v)

    int s, t, i;
    s = ppgraph->order;
    t = s + 1;
    i=0;
    while(ppgraph->adjList[u*3+i]!=v) i++; //u and v are adjacent so will stop before i == 3
    ppgraph->adjList[u*3+i]=s;
    i=0;
    while(ppgraph->adjList[v*3+i]!=u) i++; //u and v are adjacent so will stop before i == 3
    ppgraph->adjList[v*3+i]=t;

    ppgraph->degree[s]=2;
    ppgraph->adjList[s*3]=u;
    ppgraph->adjList[s*3+1]=t;
    ppgraph->multiedge[s]=t;

    ppgraph->degree[t]=2;
    ppgraph->adjList[t*3]=v;
    ppgraph->adjList[t*3+1]=s;
    ppgraph->multiedge[t]=s;

    ppgraph->multiEdgeCount++;
    ppgraph->order+=2;

    set *gu, *gv, *gs, *gt;
    gu = GRAPHROW(ppgraph->ulgraph, u, MAXM);
    gv = GRAPHROW(ppgraph->ulgraph, v, MAXM);
    gs = GRAPHROW(ppgraph->ulgraph, s, MAXM);
    gt = GRAPHROW(ppgraph->ulgraph, t, MAXM);
    EMPTYSET(gs, MAXM);
    EMPTYSET(gt, MAXM);
    DELELEMENT(gu, v);
    DELELEMENT(gv, u);
    ADDELEMENT(gu, s);
    ADDELEMENT(gs, u);
    ADDELEMENT(gs, t);
    ADDELEMENT(gv, t);
    ADDELEMENT(gt, v);
    ADDELEMENT(gt, s);
    CHECKCONSISTENCY(ppgraph)
    DEBUGMSG("End apply_deg2_operation1")
}

void revert_deg2_operation1(PRIMPREGRAPH *ppgraph, int u, int v){
    DEBUGMSG("Start revert_deg2_operation1")
    DEBUGDUMP(u, "%d")
    DEBUGDUMP(v, "%d")
    int s, t, i;
    s = ppgraph->order-2;
    t = s + 1;
    i=0;
    while(ppgraph->adjList[u*3+i]!=s) i++;
    DEBUGASSERT(i<3)
    ppgraph->adjList[u*3+i]=v;
    i=0;
    while(ppgraph->adjList[v*3+i]!=t) i++;
    DEBUGASSERT(i<3)
    ppgraph->adjList[v*3+i]=u;

    ppgraph->multiEdgeCount--;
    ppgraph->order-=2;

    set *gu, *gv;
    gu = GRAPHROW(ppgraph->ulgraph, u, MAXM);
    gv = GRAPHROW(ppgraph->ulgraph, v, MAXM);
    DELELEMENT(gu, s);
    DELELEMENT(gv, t);
    ADDELEMENT(gu, v);
    ADDELEMENT(gv, u);
    CHECKCONSISTENCY(ppgraph)
    DEBUGMSG("End revert_deg2_operation1")
}

//-----------------------------------------------------------------

/*                              s  t
 *                              o--o
 *  __  u  v  __            __ u|  |v __
 * /  \-o--o-/  \          /  \-o--o-/  \
 * |   \____/   |   ===>   |   \____/   |
 * \____________/          \____________/
 */
void apply_deg2_operation2(PRIMPREGRAPH *ppgraph, int u, int v){
    DEBUGMSG("Start apply_deg2_operation2")
    DEBUGDUMP(u, "%d")
    DEBUGDUMP(v, "%d")
    DEBUGASSERT(areAdjacent(ppgraph, u, v))
    DEBUGASSERT(ppgraph->degree[u]==2 && ppgraph->multiedge[u]==v)
    DEBUGASSERT(ppgraph->degree[v]==2 && ppgraph->multiedge[v]==u)

    int s, t;
    s = ppgraph->order;
    t = s + 1;
    ppgraph->degree[u]=3;
    ppgraph->adjList[u*3+2]=s;
    ppgraph->degree[v]=3;
    ppgraph->adjList[v*3+2]=t;

    ppgraph->degree[s]=2;
    ppgraph->adjList[s*3]=u;
    ppgraph->adjList[s*3+1]=t;
    ppgraph->multiedge[s]=t;

    ppgraph->degree[t]=2;
    ppgraph->adjList[t*3]=v;
    ppgraph->adjList[t*3+1]=s;
    ppgraph->multiedge[t]=s;

    ppgraph->order+=2;

    set *gu, *gv, *gs, *gt;
    gu = GRAPHROW(ppgraph->ulgraph, u, MAXM);
    gv = GRAPHROW(ppgraph->ulgraph, v, MAXM);
    gs = GRAPHROW(ppgraph->ulgraph, s, MAXM);
    gt = GRAPHROW(ppgraph->ulgraph, t, MAXM);
    EMPTYSET(gs, MAXM);
    EMPTYSET(gt, MAXM);
    ADDELEMENT(gu, s);
    ADDELEMENT(gs, u);
    ADDELEMENT(gs, t);
    ADDELEMENT(gv, t);
    ADDELEMENT(gt, v);
    ADDELEMENT(gt, s);
    CHECKCONSISTENCY(ppgraph)
    DEBUGMSG("End apply_deg2_operation2")
}

void revert_deg2_operation2(PRIMPREGRAPH *ppgraph, int u, int v){
    DEBUGMSG("Start revert_deg2_operation2")
    DEBUGDUMP(u, "%d")
    DEBUGDUMP(v, "%d")
    int s, t;
    s = ppgraph->order-2;
    t = s + 1;

    //s (resp t) are placed last in the adjlist of u (resp v)
    //so we can just change the degree of u and v
    ppgraph->degree[u]=2;
    ppgraph->degree[v]=2;
    ppgraph->multiedge[u]=v;
    ppgraph->multiedge[v]=u;
    ppgraph->order-=2;

    set *gu, *gv;
    gu = GRAPHROW(ppgraph->ulgraph, u, MAXM);
    gv = GRAPHROW(ppgraph->ulgraph, v, MAXM);
    DELELEMENT(gu, s);
    DELELEMENT(gv, t);
    CHECKCONSISTENCY(ppgraph)
    DEBUGMSG("End revert_deg2_operation2")
}

//-----------------------------------------------------------------

/*                             s   t
 *                             o---o
 *                              \ /
 *                               o v
 *  __  u  v  __            __  u|   __
 * /  \-o--o-/  \          /  \--o--/  \
 * |   \____/   |   ===>   |   \___/   |
 * \____________/          \___________/
 */
void apply_deg2_operation3(PRIMPREGRAPH *ppgraph, int u, int v){
    DEBUGMSG("Start apply_deg2_operation3")
    DEBUGDUMP(u, "%d")
    DEBUGDUMP(v, "%d")
    DEBUGASSERT(areAdjacent(ppgraph, u, v))
    DEBUGASSERT(ppgraph->degree[u]==2 && ppgraph->multiedge[u]==v)
    DEBUGASSERT(ppgraph->degree[v]==2 && ppgraph->multiedge[v]==u)

    int s, t, i, x;
    s = ppgraph->order;
    t = s + 1;

    i=0;
    while(ppgraph->adjList[v*3+i]==u) i++; //v has degree 2, so this will stop
    DEBUGASSERT(i<2)
    x = ppgraph->adjList[v*3+i]; //other neighbour of v
    ppgraph->adjList[u*3+2]=x;
    ppgraph->degree[u]=3;

    ppgraph->degree[v]=3;
    ppgraph->adjList[v*3+i]=s;
    ppgraph->adjList[v*3+2]=t;

    i=0;
    while(ppgraph->adjList[x*3+i]!=v) i++; //x has degree 3, so this will stop
    DEBUGASSERT(i<3)
    ppgraph->adjList[x*3+i]=u;

    ppgraph->degree[s]=2;
    ppgraph->adjList[s*3]=v;
    ppgraph->adjList[s*3+1]=t;
    ppgraph->multiedge[s]=t;

    ppgraph->degree[t]=2;
    ppgraph->adjList[t*3]=v;
    ppgraph->adjList[t*3+1]=s;
    ppgraph->multiedge[t]=s;

    ppgraph->order+=2;

    set *gu, *gv, *gs, *gt, *gx;
    gu = GRAPHROW(ppgraph->ulgraph, u, MAXM);
    gv = GRAPHROW(ppgraph->ulgraph, v, MAXM);
    gs = GRAPHROW(ppgraph->ulgraph, s, MAXM);
    gt = GRAPHROW(ppgraph->ulgraph, t, MAXM);
    gx = GRAPHROW(ppgraph->ulgraph, x, MAXM);
    EMPTYSET(gs, MAXM);
    EMPTYSET(gt, MAXM);
    DELELEMENT(gv, x);
    DELELEMENT(gx, v);
    ADDELEMENT(gx, u);
    ADDELEMENT(gu, x);
    ADDELEMENT(gv, s);
    ADDELEMENT(gv, t);
    ADDELEMENT(gs, t);
    ADDELEMENT(gs, v);
    ADDELEMENT(gt, v);
    ADDELEMENT(gt, s);

    CHECKCONSISTENCY(ppgraph)
    DEBUGMSG("End apply_deg2_operation3")
}

void revert_deg2_operation3(PRIMPREGRAPH *ppgraph, int u, int v){
    DEBUGMSG("Start revert_deg2_operation3")
    DEBUGDUMP(u, "%d")
    DEBUGDUMP(v, "%d")
    int s, t, i, x;

    s = ppgraph->order-2;
    t = s + 1;
    x = ppgraph->adjList[u*3+2];

    i=0;
    while(ppgraph->adjList[v*3+i]!=s) i++;
    DEBUGASSERT(i<2)
    ppgraph->adjList[v*3+i]=x;

    i=0;
    while(ppgraph->adjList[x*3+i]!=u) i++;
    DEBUGASSERT(i<3)
    ppgraph->adjList[x*3+i]=v;

    ppgraph->degree[u]=2;
    ppgraph->degree[v]=2;
    ppgraph->multiedge[u]=v;
    ppgraph->multiedge[v]=u;
    ppgraph->order-=2;

    set *gu, *gv, *gx;
    gu = GRAPHROW(ppgraph->ulgraph, u, MAXM);
    gv = GRAPHROW(ppgraph->ulgraph, v, MAXM);
    gx = GRAPHROW(ppgraph->ulgraph, x, MAXM);
    DELELEMENT(gv, s);
    DELELEMENT(gv, t);
    ADDELEMENT(gv, x);
    DELELEMENT(gu, x);
    ADDELEMENT(gx, v);
    DELELEMENT(gx, u);
    CHECKCONSISTENCY(ppgraph)
    DEBUGMSG("End revert_deg2_operation3")
}

//-----------------------------------------------------------------

void get_deg1_pairs(PRIMPREGRAPH *ppgraph, VERTEXPAIR *vertexPairList, int *vertexPairListSize){
    int i, j;
    *vertexPairListSize = 0;
    for(i=0;i<ppgraph->order-1;i++){
        if(ppgraph->degree[i]==1){
            for(j=i+1; j<ppgraph->order;j++){
                if(ppgraph->degree[j]==1 && ppgraph->adjList[i*3]!=j && ppgraph->adjList[i*3]!=ppgraph->adjList[j*3]){
                    vertexPairList[*vertexPairListSize][0]=i;
                    vertexPairList[*vertexPairListSize][1]=j;
                    (*vertexPairListSize)++;
                }
            }
        }
    }
}

void get_single_edges(PRIMPREGRAPH *ppgraph, VERTEXPAIR *vertexPairList, int *vertexPairListSize){
    int i, j;
    *vertexPairListSize = 0;
    for(i=0;i<ppgraph->order-1;i++){
        for(j=0; j<ppgraph->degree[i];j++){
            if(ppgraph->adjList[i*3+j]>i &&
                    !(ppgraph->degree[i]==2 && ppgraph->multiedge[i]==ppgraph->adjList[i*3+j])){
                vertexPairList[*vertexPairListSize][0]=i;
                vertexPairList[*vertexPairListSize][1]=ppgraph->adjList[i*3+j];
                (*vertexPairListSize)++;
            }
        }
    }

}

void get_multi_edges(PRIMPREGRAPH *ppgraph, VERTEXPAIR *vertexPairList, int *vertexPairListSize){
    int i;
    *vertexPairListSize = 0;
    for(i=0;i<ppgraph->order-1;i++){
        if(ppgraph->degree[i]==2 && ppgraph->multiedge[i]>i){
            vertexPairList[*vertexPairListSize][0]=i;
            vertexPairList[*vertexPairListSize][1]=ppgraph->multiedge[i];
            (*vertexPairListSize)++;
        }
    }

}

/*
 * determines the orbits of a given list of vertices pairs. The pairs must be so that the image of a pair in the list
 * is also in the list.
 *
 * In the end vertexPairOrbits[i] will contain the index of the canonical representant of the orbit to which vertexPairList[i]
 * belongs and orbitCount will contain the number of orbits.
 * The first two parameters are read-only, the last two are write-only.
 */
void determine_vertex_pairs_orbits(VERTEXPAIR *vertexPairList, int vertexPairListSize, int *vertexPairOrbits, int *orbitCount, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators){
    DEBUGMSG("Start determine_vertex_pairs_orbits")
    DEBUG2DARRAYDUMP(vertexPairList, vertexPairListSize, 2, "%d")

    int i, j, k, temp;
    int orbitSize[vertexPairListSize];

    //initialization of the variables
    for(i=0; i<vertexPairListSize; i++){
        vertexPairOrbits[i]=i;
        orbitSize[i]=1;
    }
    *orbitCount=vertexPairListSize;

    permutation *permutation;
    VERTEXPAIR pair;
    DEBUGDUMP(currentNumberOfGenerators, "%d")
    for(i = 0; i < currentNumberOfGenerators; i++) {
        permutation = (*currentGenerators)[i];
        DEBUGARRAYDUMP(permutation, currentVertexCount, "%d")

        for(j = 0; j<vertexPairListSize; j++){
            //apply permutation to current vertex pair
            pair[0] = permutation[vertexPairList[j][0]];
            pair[1] = permutation[vertexPairList[j][1]];

            //canonical form of vertex pair
            if(pair[0]>pair[1]){
                temp = pair[1];
                pair[1] = pair[0];
                pair[0] = temp;
            }

            //search the pair in the list
            for(k = 0; k<vertexPairListSize; k++){
                if(pair[0] == vertexPairList[k][0] && pair[1] == vertexPairList[k][1]){
                    unionElements(vertexPairOrbits, orbitSize, orbitCount, j, k);
                    break; //the list of pairs doesn't contain any duplicates so we can stop
                }
            }
        }
    }

    //make sure that each element is connected to its root
    for(i = 0; i < vertexPairListSize; i++){
        findRootOfElement(vertexPairOrbits, i);
    }

    DEBUGARRAYDUMP(vertexPairOrbits, vertexPairListSize, "%d")

    DEBUGMSG("End determine_vertex_pairs_orbits")
}

/*
 * determines the orbits of a given list of vertices sets. The sets must be so that the image of a set in the list
 * is also in the list.
 *
 * In the end vertexSetOrbits[i] will contain the index of the canonical representant of the orbit to which vertexSetList[i]
 * belongs and orbitCount will contain the number of orbits.
 * The first two parameters are read-only, the last two are write-only.
 */
void determine_vertex_sets_orbits(set *vertexSetList, int vertexSetListSize, int *vertexSetOrbits, int *orbitCount){
    //it is assumed that all sets have an equal number of elements
    int i, j, k, l;
    int orbitSize[vertexSetListSize];

    //initialization of the variables
    for(i=0; i<vertexSetListSize; i++){
        vertexSetOrbits[i]=i;
        orbitSize[i]=1;
    }
    *orbitCount=vertexSetListSize;

    permutation *permutation;
    set set[MAXM];
    for(i = 0; i < numberOfGenerators; i++) {
        //the generators were stored in the global variable generators by the method save_generators
        permutation = automorphismGroupGenerators[i];
        DEBUGARRAYDUMP(permutation, currentVertexCount, "%d")

        for(j = 0; j<vertexSetListSize; j++){
            //apply permutation to current vertex pair
            EMPTYSET(set, MAXM);
            for(l=-1; (l = nextelement(vertexSetList+j*MAXM, MAXM, l)) >=0;){
                ADDELEMENT(set, permutation[l]);
            }

            //search the pair in the list
            for(k = 0; k<vertexSetListSize; k++){
                //TODO: is there a better way to check whether two sets are equal?
                l = nextelement(set, MAXM, -1);
                while(l>=0 && ISELEMENT(vertexSetList+k*MAXM, l)){
                    l = nextelement(set, MAXM, l);
                }
                if(l==-1){
                    unionElements(vertexSetOrbits, orbitSize, orbitCount, j, k);
                    break; //the list of sets doesn't contain any duplicates so we can stop
                }
            }
        }
    }
}

void unionElements(int *forest, int *treeSizes, int *numberOfComponents, int element1, int element2){
    int root1 = findRootOfElement(forest, element1);
    int root2 = findRootOfElement(forest, element2);

    DEBUGMSG("Union:")
    DEBUGDUMP(element1, "%d")
    DEBUGDUMP(element2, "%d")
    DEBUGDUMP(root1, "%d")
    DEBUGDUMP(root2, "%d")

    if(root1==root2) return;

    if(treeSizes[root1]<treeSizes[root2]){
        forest[root1]=root2;
        treeSizes[root2]+=treeSizes[root1];
    } else {
        forest[root2]=root1;
        treeSizes[root1]+=treeSizes[root2];
    }
    (*numberOfComponents)--;
}

int findRootOfElement(int *forest, int element) {
    //find with path-compression
    if(element!=forest[element]){
        forest[element]=findRootOfElement(forest, forest[element]);
    }
    return forest[element];
}

//--------------------OUTPUT------------------------------------

char writePregraphTable(FILE *f, PREGRAPH *pregraph) {
    fprintf(f, "==============================\n");
    fprintf(f, "|  Graph number: %20llu  |\n", structureCount);
    fprintf(f, "|  Number of vertices: %14d  |\n", pregraph->order);
    fprintf(f, "==============================\n");

    unsigned short i, j;
    PRIMPREGRAPH *ppgraph = pregraph->ppgraph;
    int primPregraph2Pregraph[ppgraph->order];
    j = 1; //vertices are labeled starting from 1
    for (i = 0; i < ppgraph->order; i++) {
        if (ISELEMENT(pregraph->semiEdgeVertices, i)) {
            primPregraph2Pregraph[i] = pregraph->order + 1;
        } else {
            primPregraph2Pregraph[i] = j++;
        }
    }
    for (i = 0; i < ppgraph->order; i++) {
        if (primPregraph2Pregraph[i] != pregraph->order + 1) { //don't write vertices that correspond to semi-edges
            fprintf(f, "|%4d ||", primPregraph2Pregraph[i]);
            for (j = 0; j < ppgraph->degree[i]; j++) {
                if(primPregraph2Pregraph[ppgraph->adjList[i * 3 + j]] == pregraph->order + 1){
                    fprintf(f, "    S |");
                } else {
                    fprintf(f, " %4d |", primPregraph2Pregraph[ppgraph->adjList[i * 3 + j]]);
                }
            }
            if (j == 1) {
                //add loop
                fprintf(f, "    L |      |");
            } else if (j == 2) {
                //add multi-edge
                fprintf(f, " %4d |", primPregraph2Pregraph[ppgraph->multiedge[i]]);
            }
            fprintf(f,"|\n");
        }
    }
    fprintf(f, "==============================\n");
    fprintf(f,"\n");
    return (ferror(f) ? 2 : 1);
}

char write_2byte_number(FILE *f, unsigned short n, short writeEndian) {
    if (writeEndian == BIG_ENDIAN) {
        fprintf(f, "%c%c", n / 256, n % 256);
    } else {
        fprintf(f, "%c%c", n % 256, n / 256);
    }
    return (ferror(f) ? 2 : 1);
}

char writePregraphCode(FILE *f, PREGRAPH *pregraph) {
    unsigned short i, j;
    PRIMPREGRAPH *ppgraph = pregraph->ppgraph;
    unsigned short primPregraph2Pregraph[ppgraph->order];
    j = 1; //vertices are labeled starting from 1
    for (i = 0; i < ppgraph->order; i++) {
        if (ISELEMENT(pregraph->semiEdgeVertices, i)) {
            primPregraph2Pregraph[i] = pregraph->order + 1;
        } else {
            primPregraph2Pregraph[i] = j++;
        }
    }
    if (structureCount == 1) { //if first graph
        fprintf(f, ">>pregraph_code %s<<", (endian == LITTLE_ENDIAN ? "le" : "be"));
    }
    if (pregraph->order + 1 <= UCHAR_MAX) {
        fprintf(f, "%c", (unsigned char) pregraph->order);
    } else {
        fprintf(f, "%c", 0);
        /* big graph */
        if (write_2byte_number(f, (unsigned short) pregraph->order, endian) == 2) {
            return (2);
        }
    }
    for (i = 0; i < ppgraph->order; i++) {
        if (primPregraph2Pregraph[i] != pregraph->order + 1) { //don't write vertices that correspond to semi-edges
            for (j = 0; j < ppgraph->degree[i]; j++) {
                if(primPregraph2Pregraph[ppgraph->adjList[i * 3 + j]]>primPregraph2Pregraph[i]){
                    //only include adjacency information for vertices with a larger index
                    if (pregraph->order + 1 <= UCHAR_MAX) {
                        fprintf(f, "%c", (unsigned char) primPregraph2Pregraph[ppgraph->adjList[i * 3 + j]]);
                    } else {
                        if (write_2byte_number(f, primPregraph2Pregraph[ppgraph->adjList[i * 3 + j]], endian) == 2) {
                            return (2);
                        }
                    }
                }
            }
            if (j == 1) {
                //add loop
                if (pregraph->order + 1 <= UCHAR_MAX) {
                    fprintf(f, "%c", (unsigned char) primPregraph2Pregraph[i]);
                } else {
                    if (write_2byte_number(f, primPregraph2Pregraph[i], endian) == 2) {
                        return (2);
                    }
                }
            } else if (j == 2 && primPregraph2Pregraph[ppgraph->multiedge[i]]>primPregraph2Pregraph[i]) {
                //add multi-edge if it has a larger index
                if (pregraph->order + 1 <= UCHAR_MAX) {
                    fprintf(f, "%c", (unsigned char) primPregraph2Pregraph[ppgraph->multiedge[i]]);
                } else {
                    if (write_2byte_number(f, primPregraph2Pregraph[ppgraph->multiedge[i]], endian) == 2) {
                        return (2);
                    }
                }
            }
            //closing 0
            if (pregraph->order + 1 <= UCHAR_MAX) {
                fprintf(f, "%c", 0);
            } else {
                if (write_2byte_number(f, 0, endian) == 2) {
                    return (2);
                }
            }
        }
    }
    return (ferror(f) ? 2 : 1);
}

char writePrimpregraphTable(FILE *f, PRIMPREGRAPH *ppgraph) {
    fprintf(f, "==============================\n");
    fprintf(f, "|  Graph number: %20llu  |\n", primitivesCount);
    fprintf(f, "|  Number of vertices: %14d  |\n", ppgraph->order);
    fprintf(f, "|  Number deg1 vertices: %12d  |\n", ppgraph->degree1Count);
    fprintf(f, "|  Number multi-edges: %14d  |\n", ppgraph->multiEdgeCount);
    fprintf(f, "==============================\n");

    unsigned short i, j;
    for (i = 0; i < ppgraph->order; i++) {
        fprintf(f, "|%4d ||", i+1);
        for (j = 0; j < ppgraph->degree[i]; j++) {
            fprintf(f, " %4d |", ppgraph->adjList[i * 3 + j] + 1);
        }
        if(j==2){
            fprintf(f, "(%4d)|", ppgraph->multiedge[i] + 1);
        }
        fprintf(f,"|\n");
    }
    fprintf(f, "==============================\n");
    fprintf(f,"\n");
    return (ferror(f) ? 2 : 1);
}

char writePrimpregraphCode(FILE *f, PRIMPREGRAPH *ppgraph) {
    unsigned short i, j;
    if (structureCount == 1) { //if first graph
        fprintf(f, ">>multi_code %s<<", (endian == LITTLE_ENDIAN ? "le" : "be"));
    }
    if (ppgraph->order <= UCHAR_MAX) {
        fprintf(f, "%c", (unsigned char) ppgraph->order);
    } else {
        fprintf(f, "%c", 0);
        /* big graph */
        if (write_2byte_number(f, (unsigned short) ppgraph->order, endian) == 2) {
            return (2);
        }
    }
    for (i = 0; i < ppgraph->order - 1; i++) {
            for (j = 0; j < ppgraph->degree[i]; j++) {
                if(ppgraph->adjList[i * 3 + j]>i){
                    if (ppgraph->order + 1 <= UCHAR_MAX) {
                        fprintf(f, "%c", (unsigned char) ppgraph->adjList[i * 3 + j] + 1);
                    } else {
                        if (write_2byte_number(f, ppgraph->adjList[i * 3 + j] + 1, endian) == 2) {
                            return (2);
                        }
                    }
                }
            }
            //closing 0
            if (ppgraph->order <= UCHAR_MAX) {
                fprintf(f, "%c", 0);
            } else {
                if (write_2byte_number(f, 0, endian) == 2) {
                    return (2);
                }
            }
    }
    return (ferror(f) ? 2 : 1);
}

/*
void write_pregraph(PREGRAPH *pregraph, FILE *file){
    fprintf(file, "%d ", pregraph->order);
    int i, j;
    for (i = 0; i < pregraph->order; i++) {
        for (j = 0; j < (pregraph->ppgraph)->degree[i]; j++) {
        }
        fprintf(file, "0 ", pregraph->order);
    }
}
*/

void handle_pregraph_result(PREGRAPH *pregraph){
    DEBUGMSG("Start handle_pregraph_result")
    structureCount++;
    DEBUGDUMP(structureCount, "%ld")
    if(outputType != 'n'){
        FILE *file = stdout;
        if(outputFile != NULL){
            file = fopen(outputFile, "a");
            DEBUGMSG("Opened file")
        }

        if(outputType == 'h'){
            if (writePregraphTable(file, pregraph) == 2) {
                fprintf(stderr, "Error while writing graph %llu\n", structureCount);
                exit(1);
            }
        } else if(outputType == 'c'){
            if (writePregraphCode(file, pregraph) == 2) {
                fprintf(stderr, "Error while writing graph %llu\n", structureCount);
                exit(1);
            }
        }

        if(outputFile != NULL){
            fclose(file);
            DEBUGMSG("Closed file")
        }
    }
    if(logStatistics) logInfo(pregraph);
    DEBUGMSG("End handle_pregraph_result")
}

void doPrimpregraphExport(PRIMPREGRAPH *ppgraph){
    if(outputType != 'n'){
        FILE *file = stdout;
        if(outputFile != NULL){
            file = fopen(outputFile, "a");
            DEBUGMSG("Opened file")
        }

        if(outputType == 'h'){
            if (writePrimpregraphTable(file, ppgraph) == 2) {
                fprintf(stderr, "Error while writing pregraph primitive %llu\n", primitivesCount);
                exit(1);
            }
        } else if(outputType == 'c'){
            if (writePrimpregraphCode(file, ppgraph) == 2) {
                fprintf(stderr, "Error while writing pregraph primitive %llu\n", primitivesCount);
                exit(1);
            }
        }

        if(outputFile != NULL){
            fclose(file);
            DEBUGMSG("Closed file")
        }
    }
}

/*
 * Handles a result in the form of a primitive pregraph. If both semi-edges and loops are
 * allowed all the different ways to select these are taken and passed on to handle_pregraph_result.
 */
void handle_primpregraph_result(PRIMPREGRAPH *ppgraph){
    DEBUGMSG("Start handle_primpregraph_result")
    DEBUGPPGRAPHPRINT(ppgraph)
    DEBUGDUMP(ppgraph->order, "%d")
    DEBUGDUMP(vertexCount, "%d")
    DEBUGASSERT(ppgraph->order >= vertexCount)
    DEBUGASSERT(allowSemiEdges || vertexCount == ppgraph->order)

    int semiEdgeCount = ppgraph->order - vertexCount; DEBUGDUMP(semiEdgeCount, "%d")
    int degree1Count = ppgraph->degree1Count; DEBUGDUMP(degree1Count, "%d")
    DEBUGASSERT(semiEdgeCount <= degree1Count)
    int loopCount = degree1Count - semiEdgeCount;

    if(!allowLoops && loopCount>0){
        DEBUGMSG("End handle_primpregraph_result")
        return;
    }

    primitivesCount++;
    if(onlyPrimitives){
        doPrimpregraphExport(ppgraph);
        DEBUGMSG("End handle_primpregraph_result")
        return;
    }

    //determine up to automorphism all the ways to select semiEdgeCount vertices
    //of degree 1 by using union-find
    int listSize, i;
    listSize = 1;
    if(loopCount>0){
        //otherwise listSize would be reset to 0 if semiEdgeCount == degree1Count
        for(i=1; i<=semiEdgeCount;i++){
            listSize = listSize*(degree1Count - i + 1)/i;
        }
    }
    set vertexSetList[MAXM*listSize];
    for(i=0; i<listSize;i++){
        EMPTYSET(vertexSetList+i*MAXM,MAXM);
    }

    if(semiEdgeCount>0){
        int position = 0;
        set tempSet[MAXM];
        EMPTYSET(tempSet, MAXM);
        determine_possible_sets_of_degree1_vertices(tempSet, vertexSetList, &position, semiEdgeCount, 0, nextDegree1Vertex(-1, ppgraph), ppgraph, 0);
        DEBUGDUMP(position, "%d")
        DEBUGDUMP(listSize, "%d")
        DEBUGASSERT(position==listSize)
        #ifdef _DEBUG
        // print the sets
        for(i=0;i<listSize;i++){
            fprintf(stderr, "%s:%u set %d= [", __FILE__, __LINE__, i);
            int l;
            for(l=-1; (l = nextelement(vertexSetList + i*MAXM, MAXM, l)) >=0;){
                fprintf(stderr, "%d ", l);
            }
            fprintf(stderr, "]\n");
        }
        #endif

    } //else: listsize is already 1 and that set is already empty

    int orbitCount;
    int orbits[listSize];
    determine_vertex_sets_orbits(vertexSetList, listSize, orbits, &orbitCount);

    //output pregraph
    PREGRAPH pregraph;
    pregraph.order = vertexCount;
    pregraph.ppgraph = ppgraph;
    for(i=0; i<listSize; i++){
        if(orbits[i]==i){
            pregraph.semiEdgeVertices = vertexSetList + i*MAXM;
            handle_pregraph_result(&pregraph);
        }
    }
    DEBUGMSG("End handle_primpregraph_result")
}

int nextDegree1Vertex(int current, PRIMPREGRAPH *ppgraph){
    current++;
    while(current<ppgraph->order && ppgraph->degree[current]!=1) current++;
    if(current==ppgraph->order)
        return -1;
    else
        return current;
}

void determine_possible_sets_of_degree1_vertices
(set *tempSet, set *vertexSetList, int* currentListPosition, int maximumSetSize,
        int currentSetSize, int currentSetElement, PRIMPREGRAPH *ppgraph,
        int skippedVertices){
    DEBUGMSG("Start determine_possible_sets_of_degree1_vertices")
    if(ppgraph->degree1Count-skippedVertices < maximumSetSize){
        DEBUGMSG("End determine_possible_sets_of_degree1_vertices")
        return;
    }
    if(currentSetElement!=-1){
        ADDELEMENT(tempSet, currentSetElement);
        DEBUGDUMP(currentSetElement, "%d added")
        if(currentSetSize + 1 == maximumSetSize){
            //add to list
            int i;
            for(i=0;i<MAXM;i++){
                vertexSetList[(*currentListPosition)*MAXM+i] = tempSet[i];
            }
            (*currentListPosition)++;
        } else {
            determine_possible_sets_of_degree1_vertices
                    (tempSet, vertexSetList, currentListPosition, maximumSetSize,
                    currentSetSize + 1, nextDegree1Vertex(currentSetElement, ppgraph),
                    ppgraph, skippedVertices);
        }
        DELELEMENT(tempSet, currentSetElement);
        DEBUGDUMP(currentSetElement, "%d removed")
        determine_possible_sets_of_degree1_vertices
               (tempSet, vertexSetList, currentListPosition, maximumSetSize,
                currentSetSize, nextDegree1Vertex(currentSetElement, ppgraph),
                ppgraph, skippedVertices+1);
    }
    DEBUGMSG("End determine_possible_sets_of_degree1_vertices")
}

/*
 * Handles the result of a degree 1 operation. If the graph has a valid size it is
 * passed on to handle_primpregraph_result, then we continue with degree 1 operations
 * and degree 2 operations.
 */
void handle_deg1_operation_result(PRIMPREGRAPH *ppgraph){
    DEBUGMSG("Start handle_deg1_operation_result")
    degree1OperationsDepth++;
    if(degree1OperationsDepth>degree1OperationsDepthMaximum) degree1OperationsDepthMaximum = degree1OperationsDepth;
    //if correct number of vertices
    if(ppgraph->order >= minVertexCount && ppgraph->order<=maxVertexCount && ppgraph->order - vertexCount <= ppgraph->degree1Count)
        handle_primpregraph_result(ppgraph);

    if(ppgraph->order - ppgraph->degree1Count + 1 > vertexCount){
        //+1 because deg1 operation always adds 1 vertex of degree 3 and degree 2 operation always add 2 vertices (of degree 2 or 3)
        //too many degree 3 and 2 vertices
        DEBUGMSG("End handle_deg1_operation_result")
        degree1OperationsDepth--;
        return;
    }

    permutation currentGenerators[MAXN][MAXN]; //TODO: can't we make this smaller because we now the size at this point
    int currentNumberOfGenerators = numberOfGenerators;
    copyGenerators(&currentGenerators, ppgraph->order);

    #ifdef _DEBUG
    //check that the generators were copied correctly
    int i, j;
    for(i=0; i<numberOfGenerators; i++){
        for (j = 0; j < ppgraph->order; j++) {
            DEBUGASSERT(currentGenerators[i][j]==automorphismGroupGenerators[i][j])
        }
    }
    DEBUGDUMP(currentNumberOfGenerators, "%d")
    DEBUG2DARRAYDUMP(currentGenerators, currentNumberOfGenerators, ppgraph->order, "%d")
    #endif

    do_deg1_operations(ppgraph, &currentGenerators, currentNumberOfGenerators); //when this returns &ppgraph is unchanged

    if(allowMultiEdges && ppgraph->order - ppgraph->degree1Count + 2 <= vertexCount){
        do_deg2_operations(ppgraph, &currentGenerators, currentNumberOfGenerators, FALSE);
    }
    degree1OperationsDepth--;
    DEBUGMSG("End handle_deg1_operation_result")
}

/*
 * Handles the result of a degree 2 operation. If the graph has a valid size it is
 * passed on to handle_primpregraph_result, then we continue with degree 2 operations.
 * Degree 1 operations are no longer possible.
 */
void handle_deg2_operation_result(PRIMPREGRAPH *ppgraph, boolean multiEdgesDetermined){
    DEBUGMSG("Start handle_deg2_operation_result")
    degree2OperationsDepth++;
    if(degree2OperationsDepth>degree2OperationsDepthMaximum) degree2OperationsDepthMaximum = degree2OperationsDepth;
    //if correct number of vertices
    if(ppgraph->order >= minVertexCount && ppgraph->order<=maxVertexCount && ppgraph->order - vertexCount <= ppgraph->degree1Count)
        handle_primpregraph_result(ppgraph);

    if(ppgraph->order - ppgraph->degree1Count + 2 > vertexCount){
        //+2 because degree 2 operations always add 2 vertices (of degree 2 or 3)
        //too many degree 3 and 2 vertices
        DEBUGMSG("End handle_deg2_operation_result")
        degree2OperationsDepth--;
        return;
    }

    permutation currentGenerators[MAXN][MAXN]; //TODO: can't we make this smaller because we now the size at this point
    int currentNumberOfGenerators = numberOfGenerators;
    copyGenerators(&currentGenerators, ppgraph->order);

    #ifdef _DEBUG
    //check that the generators were copied correctly
    int i, j;
    for(i=0; i<numberOfGenerators; i++){
        for (j = 0; j < ppgraph->order; j++) {
            DEBUGASSERT(currentGenerators[i][j]==automorphismGroupGenerators[i][j])
        }
    }
    DEBUG2DARRAYDUMP(currentGenerators, numberOfGenerators, ppgraph->order, "%d")
    #endif

    do_deg2_operations(ppgraph, &currentGenerators, currentNumberOfGenerators, multiEdgesDetermined);
    degree2OperationsDepth--;
    DEBUGMSG("End handle_deg2_operation_result")
}

/**
 * Returns a vector with the number of vertices at a distance d <= maxdistance
 * of v. Uses BFS. Vertices at distance 0,1 and 2 aren't included.
 */
void vertexNeighbourhoodSizeVector(PRIMPREGRAPH *ppgraph, int v, int maxDistance, int *vector){
    static int upperboundDistance[MAXN];
    static int queue[MAXN], queueHead, queueTail;
    int i;

    for(i=0; i<maxDistance-3+1; i++) vector[i]=0;

    RESET_MARK

    queue[0] = v;
    queueHead = -1;
    queueTail = 0;
    upperboundDistance[v] = 0;
    SET_MARK(v);

    while (queueHead < queueTail) {
        queueHead++;
        int currentVertex = queue[queueHead];
        int currentDistance = upperboundDistance[currentVertex];
        if(currentDistance >= 3)
            vector[currentDistance-3]++;
        if(currentDistance < maxDistance){
            for(i = 0; i < ppgraph->degree[currentVertex]; i++){
                if(MARKED(ppgraph->adjList[currentVertex*3+i])){
                    //we already have calculated a 'distance' for this vertex
                    //only check if we need to improve it
                    if(upperboundDistance[ppgraph->adjList[currentVertex*3+i]] >
                            currentDistance + 1){
                        upperboundDistance[ppgraph->adjList[currentVertex*3+i]]
                                = currentDistance + 1;
                    }
                } else {
                    //first time that we see this vertex
                    upperboundDistance[ppgraph->adjList[currentVertex*3+i]]
                            = currentDistance + 1;
                    //mark the vertex and put it in the queue
                    SET_MARK(ppgraph->adjList[currentVertex*3+i]);
                    queueTail++;
                    queue[queueTail]=ppgraph->adjList[currentVertex*3+i];
                }
            }
        }

    }
}

/**
 * Colours each vertex of degree 1 with a vector with the number of vertices
 * at a distance d <= maxdistance of that vertex.
 * Returns a vertex with the smallest 'colour' used.
 */
int colourDegree1VertexNeighbourhoodSizeVector(PRIMPREGRAPH *ppgraph, int maxDistance, int* minimumCount, int* representants){
    int i, j, minimum = -1;
    int vectorSize = maxDistance - 3 + 1;
    int colours[MAXN*vectorSize];
    for(i = 0; i < ppgraph->order; i++){
        if(ppgraph->degree[i]==1){
            vertexNeighbourhoodSizeVector(ppgraph, i, maxDistance, colours+i*vectorSize);
            DEBUGDUMP(i, "%d")
            DEBUGARRAYDUMP(colours, vectorSize, "%d")
            if(minimum==-1){
                minimum = i;
                representants[i]=i;
                (*minimumCount)=1;
            } else {
                j = 0;
                while(j<vectorSize &&
                        colours[minimum*vectorSize+j]==colours[i*vectorSize+j])
                    j++;
                if(j!=vectorSize &&
                        colours[minimum*vectorSize+j]>colours[i*vectorSize+j]){
                    minimum = i;
                    representants[i]=i;
                    (*minimumCount)=1;
                } else if(j==vectorSize){
                    (*minimumCount)++;
                    representants[i]=minimum;
                } else {
                    representants[i]=-1;
                }
            }
        }
    }
    DEBUGDUMP(minimum, "%d")
    return minimum;
}

/**
 * Returns the number of vertices at a distance d <= maxdistance of v.
 * Uses BFS.
 */
int vertexNeighbourhoodSize(PRIMPREGRAPH *ppgraph, int v, int maxDistance){
    int neighbourhoodSize = 0;
    static int upperboundDistance[MAXN];
    static int queue[MAXN], queueHead, queueTail;
    int i;

    RESET_MARK

    queue[0] = v;
    queueHead = -1;
    queueTail = 0;
    upperboundDistance[v] = 0;
    SET_MARK(v);

    while (queueHead < queueTail) {
        neighbourhoodSize++;
        queueHead++;
        int currentVertex = queue[queueHead];
        int currentDistance = upperboundDistance[currentVertex];
        if(currentDistance < maxDistance){
            for(i = 0; i < ppgraph->degree[currentVertex]; i++){
                if(MARKED(ppgraph->adjList[currentVertex*3+i])){
                    //we already have calculated a 'distance' for this vertex
                    //only check if we need to improve it
                    if(upperboundDistance[ppgraph->adjList[currentVertex*3+i]] >
                            currentDistance + 1){
                        upperboundDistance[ppgraph->adjList[currentVertex*3+i]]
                                = currentDistance + 1;
                    }
                } else {
                    //first time that we see this vertex
                    upperboundDistance[ppgraph->adjList[currentVertex*3+i]]
                            = currentDistance + 1;
                    //mark the vertex and put it in the queue
                    SET_MARK(ppgraph->adjList[currentVertex*3+i]);
                    queueTail++;
                    queue[queueTail]=ppgraph->adjList[currentVertex*3+i];
                }
            }
        }

    }

    return neighbourhoodSize;
}

/**
 * Colours each vertex of degree 1 with the number of vertices
 * at a distance d <= maxdistance of that vertex.
 * Returns the smallest 'colour' used.
 */
int colourDegree1VertexNeighbourhoodSize(PRIMPREGRAPH *ppgraph, int maxDistance, int* colours){
    int i, minimum = MAXN;
    for(i = 0; i < ppgraph->order; i++){
        if(ppgraph->degree[i]==1){
            colours[i] = vertexNeighbourhoodSize(ppgraph, i, maxDistance);
            DEBUGDUMP(i, "%d")
            DEBUGDUMP(colours[i], "%d")
            if(colours[i] < minimum)
                minimum = colours[i];
        }
    }
    DEBUGDUMP(minimum, "%d")
    return minimum;
}

boolean isCanonicalDegree1Edge(PRIMPREGRAPH *ppgraph, int v){
    DEBUGASSERT(ppgraph->degree[v]==1)
    DEBUGDUMP(v,"%d")
    int i, j;

    if(ppgraph->degree1Count==1){
        //only one degree 1 vertex: garantueed to be canonical
        //call nauty and return true
        int vertexOrbits[ppgraph->order];
        DEBUGMSG("Start nauty")
        numberOfGenerators = 0; //reset the generators
        nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, vertexOrbits, &nautyOptions, &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order, canonicalGraph);
        DEBUGMSG("End nauty")
        DEBUGARRAYDUMP(vertexOrbits, ppgraph->order, "%d")
        return TRUE;
    }

    //first check colors (which should be cheaper than always calling nauty)
    int degree1Vertices[ppgraph->degree1Count];
    int colours[MAXN];
    j = 0;
    for (i = 0; i < ppgraph->order; i++){
        if(ppgraph->degree[i] == 1){
            degree1Vertices[j] = i;
            j++;
        }
    }
    DEBUGASSERT(j==ppgraph->degree1Count)

    int minimumColour = colourDegree1VertexNeighbourhoodSize(ppgraph, 4, colours);

    //if v hasn't got the smallest colour, then it isn't canonical
    if(minimumColour != colours[v]) return FALSE;

    int minimumColourCount = 0;
    for (i = 0; i < ppgraph->degree1Count; i++){
        if(minimumColour == colours[degree1Vertices[i]]){
            minimumColourCount++;
        }
    }

    if(minimumColourCount==1){
        //only one degree 1 vertex with minimal colour, i.e. v is canonical
        //call nauty and return true
        int vertexOrbits[ppgraph->order];
        DEBUGMSG("Start nauty")
        numberOfGenerators = 0; //reset the generators
        nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, vertexOrbits, &nautyOptions, &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order, canonicalGraph);
        DEBUGMSG("End nauty")
        DEBUGARRAYDUMP(vertexOrbits, ppgraph->order, "%d")
        return TRUE;
    }

    //just call nauty and we'll be sure whether it is canonical

    int vertexOrbits[ppgraph->order];
    DEBUGMSG("Start nauty")
    numberOfGenerators = 0; //reset the generators
    nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, vertexOrbits, &nautyOptions, &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order, canonicalGraph);
    DEBUGMSG("End nauty")
    DEBUGARRAYDUMP(vertexOrbits, ppgraph->order, "%d")

    DEBUGARRAYDUMP(nautyLabelling, ppgraph->order, "%d")
    int reverseLabelling[ppgraph->order];
    for (i = 0; i < ppgraph->order; i++) {
        reverseLabelling[nautyLabelling[i]]=i;
    }
    DEBUGARRAYDUMP(reverseLabelling, ppgraph->order, "%d")
    int smallestLabelOrbitV = reverseLabelling[v]; //i.e. the smallest label of a vertex in the orbit of vertex v
    int smallestOtherDegree1Label = ppgraph->order;
    for (i = 0; i < ppgraph->order; i++) {
        if(ppgraph->degree[i]==1){
            if(vertexOrbits[i]==vertexOrbits[v]){
                if(reverseLabelling[i]<smallestLabelOrbitV) smallestLabelOrbitV = reverseLabelling[i];
            } else if(colours[i] == minimumColour){
                if(reverseLabelling[i]<smallestOtherDegree1Label) smallestOtherDegree1Label = reverseLabelling[i];
            }
        }
    }
    DEBUGDUMP(smallestLabelOrbitV, "%d")
    DEBUGDUMP(smallestOtherDegree1Label, "%d")
    if(noRejections) return TRUE;
    return (smallestLabelOrbitV < smallestOtherDegree1Label);
}

/*
 * Unused because more expansive than the time that was saved
 *
boolean isCanonicalDegree1Edge(PRIMPREGRAPH *ppgraph, int v){
    DEBUGASSERT(ppgraph->degree[v]==1)
    DEBUGDUMP(v,"%d")
    int i, j;

    if(ppgraph->degree1Count==1){
        //only one degree 1 vertex: garantueed to be canonical
        //call nauty and return true
        int vertexOrbits[ppgraph->order];
        DEBUGMSG("Start nauty")
        numberOfGenerators = 0; //reset the generators
        nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, vertexOrbits, &nautyOptions, &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order, canonicalGraph);
        DEBUGMSG("End nauty")
        DEBUGARRAYDUMP(vertexOrbits, ppgraph->order, "%d")
        return TRUE;
    }

    //first check colors (which should be cheaper than always calling nauty)
    int degree1Vertices[ppgraph->degree1Count];
    int representants[MAXN];
    int minimumColourCount;
    j = 0;
    for (i = 0; i < ppgraph->order; i++){
        if(ppgraph->degree[i] == 1){
            degree1Vertices[j] = i;
            j++;
        }
    }
    DEBUGASSERT(j==ppgraph->degree1Count)

    int minimumColourVertex = colourDegree1VertexNeighbourhoodSizeVector(ppgraph, 5, &minimumColourCount, representants);

    //if v hasn't got the smallest colour, then it isn't canonical
    if(representants[v]!=minimumColourVertex) return FALSE;


    if(minimumColourCount==1){
        //only one degree 1 vertex with minimal colour, i.e. v is canonical
        //call nauty and return true
        int vertexOrbits[ppgraph->order];
        DEBUGMSG("Start nauty")
        numberOfGenerators = 0; //reset the generators
        nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, vertexOrbits, &nautyOptions, &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order, canonicalGraph);
        DEBUGMSG("End nauty")
        DEBUGARRAYDUMP(vertexOrbits, ppgraph->order, "%d")
        return TRUE;
    }

    //just call nauty and we'll be sure whether it is canonical

    int vertexOrbits[ppgraph->order];
    DEBUGMSG("Start nauty")
    numberOfGenerators = 0; //reset the generators
    nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, vertexOrbits, &nautyOptions, &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order, canonicalGraph);
    DEBUGMSG("End nauty")
    DEBUGARRAYDUMP(vertexOrbits, ppgraph->order, "%d")

    DEBUGARRAYDUMP(nautyLabelling, ppgraph->order, "%d")
    int reverseLabelling[ppgraph->order];
    for (i = 0; i < ppgraph->order; i++) {
        reverseLabelling[nautyLabelling[i]]=i;
    }
    DEBUGARRAYDUMP(reverseLabelling, ppgraph->order, "%d")
    int smallestLabelOrbitV = reverseLabelling[v]; //i.e. the smallest label of a vertex in the orbit of vertex v
    int smallestOtherDegree1Label = ppgraph->order;
    for (i = 0; i < ppgraph->order; i++) {
        if(ppgraph->degree[i]==1){
            if(vertexOrbits[i]==vertexOrbits[v]){
                if(reverseLabelling[i]<smallestLabelOrbitV) smallestLabelOrbitV = reverseLabelling[i];
            } else if(representants[i] == minimumColourVertex){
                if(reverseLabelling[i]<smallestOtherDegree1Label) smallestOtherDegree1Label = reverseLabelling[i];
            }
        }
    }
    DEBUGDUMP(smallestLabelOrbitV, "%d")
    DEBUGDUMP(smallestOtherDegree1Label, "%d")
    if(noRejections) return TRUE;
    return (smallestLabelOrbitV < smallestOtherDegree1Label);
}
*/

/*
 * Handles the first degree 1 operation, i.e. find all the degree 1 pairs, determine the orbits
 * apply the operation for each pair, handle the result and then revert the operation.
 */
void handle_deg1_operation1(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators){
    if(operation11Disabled) return;
    DEBUGMSG("Start handle_deg1_operation1")
    DEBUG2DARRAYDUMP((*currentGenerators), currentNumberOfGenerators, ppgraph->order, "%d")
    int maxSize = ppgraph->degree1Count*ppgraph->degree1Count/2;
    VERTEXPAIR deg1PairList[maxSize]; //initialize an array that is large enough to hold all the degree 1 pairs
    int listSize;
    get_deg1_pairs(ppgraph, deg1PairList, &listSize);

    int orbitCount;
    int orbits[listSize];
    determine_vertex_pairs_orbits(deg1PairList, listSize, orbits, &orbitCount, currentGenerators, currentNumberOfGenerators);

    int i;
    for (i = 0; i < listSize; i++) {
        if(orbits[i]==i){
            DEBUGPPGRAPHPRINT(ppgraph)
            apply_deg1_operation1(ppgraph, deg1PairList[i][0], deg1PairList[i][1]);

            //the only deg 1 vertex after this operation is v. This is a valid action
            //if v belongs to the first orbit of degree 1 vertices

            if(isCanonicalDegree1Edge(ppgraph, deg1PairList[i][1])){
                //v belongs to the orbit of degree 1 vertices with the smallest representant
                //Therefore this graph was created from the correct parent.
                handle_deg1_operation_result(ppgraph);
            }

            revert_deg1_operation1(ppgraph, deg1PairList[i][0], deg1PairList[i][1]);
            DEBUGPPGRAPHPRINT(ppgraph)
        }
    }
    DEBUGMSG("End handle_deg1_operation1")
}

void handle_deg1_operation2(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators){
    if(operation12Disabled) return;
    DEBUGMSG("Start handle_deg1_operation2")
    DEBUG2DARRAYDUMP((*currentGenerators), currentNumberOfGenerators, ppgraph->order, "%d")
    int maxSize = ppgraph->order*3/2-ppgraph->degree1Count; //this upper bound is not tight (it is tight in case of no degree 2 vertices?)
    VERTEXPAIR edgeList[maxSize]; //initialize an array that is large enough to hold all single edges
    int listSize;
    get_single_edges(ppgraph, edgeList, &listSize);

    int orbitCount;
    int orbits[listSize];
    determine_vertex_pairs_orbits(edgeList, listSize, orbits, &orbitCount, currentGenerators, currentNumberOfGenerators);
    //TODO: the calculation above is done both for degree 1 operation 2 and degree 2 operation 1: avoid duplicating this work!!!
    DEBUGARRAYDUMP(orbits, orbitCount, "%d")

    int i;
    for (i = 0; i < listSize; i++) {
        if(orbits[i]==i){
            //TODO: maybe only enumerate the bridges?
            if(isBridge(ppgraph, edgeList[i][0], edgeList[i][1])){
                DEBUGPPGRAPHPRINT(ppgraph)
                apply_deg1_operation2(ppgraph, edgeList[i][0], edgeList[i][1]);

                //the new deg 1 vertex after this operation is t. This is a valid action
                //if t belongs to the first orbit of degree 1 vertices

                if(isCanonicalDegree1Edge(ppgraph, ppgraph->order-1)){
                    //t belongs to the orbit of degree 1 vertices with the smallest representant
                    //Therefore this graph was created from the correct parent.
                    handle_deg1_operation_result(ppgraph);
                }

                revert_deg1_operation2(ppgraph, edgeList[i][0], edgeList[i][1]);
                DEBUGPPGRAPHPRINT(ppgraph)
            }
        }
    }
    DEBUGMSG("End handle_deg1_operation2")
}

/*
 * Returns TRUE if the edge (v1, v2) belongs to the orbit with the canonically smallest label.
 * Also stores the list of multi-edges, its size, its orbits and the number of orbits.
 */
boolean isCanonicalMultiEdge(PRIMPREGRAPH *ppgraph, int v1, int v2, boolean *multiEdgesDetermined){
    DEBUGMSG("Start isCanonicalMultiEdge")
    //TODO: optimize if multiEdgeCount == 1
    if(v2<v1){
        int temp = v1;
        v1 = v2;
        v2 = temp;
    }
    int orbits[ppgraph->order + ppgraph->multiEdgeCount]; //the graph is enlarged so we have to provide a large enough array

    VERTEXPAIR *multiEdgeList = globalMultiEdgeList + ((degree2OperationsDepth+1)*HALFFLOOR(vertexCount));
    int *multiEdgeListSize = globalMultiEdgeListSize + degree2OperationsDepth + 1;

    get_multi_edges(ppgraph, multiEdgeList, multiEdgeListSize);

    DEBUGASSERT(*multiEdgeListSize == ppgraph->multiEdgeCount)

    DEBUGMSG("Start nauty")
    numberOfGenerators = 0; //reset the generators
    nautyOptions.defaultptn = FALSE; //use colourings for multigraphs

    int j;
    for(j = 0; j < ppgraph->order + ppgraph->multiEdgeCount; j++){
        nautyLabelling[j]=j;
        nautyPtn[j]=(j!=ppgraph->order-1);
    }
    for(j = 0; j < ppgraph->multiEdgeCount; j++){
        int n1 = multiEdgeList[j][0];
        int n2 = multiEdgeList[j][1];
        set *multiEdgeVertex, *v1, *v2;
        multiEdgeVertex = GRAPHROW(ppgraph->ulgraph, ppgraph->order + j, MAXM);
        v1 = GRAPHROW(ppgraph->ulgraph, n1, MAXM);
        v2 = GRAPHROW(ppgraph->ulgraph, n2, MAXM);
        EMPTYSET(multiEdgeVertex, MAXM);
        DELELEMENT(v1, n2);
        DELELEMENT(v2, n1);
        ADDELEMENT(v1, ppgraph->order + j);
        ADDELEMENT(v2, ppgraph->order + j);
        ADDELEMENT(multiEdgeVertex, n1);
        ADDELEMENT(multiEdgeVertex, n2);
    }

    nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, orbits, &nautyOptions, &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order + ppgraph->multiEdgeCount, canonicalGraph);
    nautyOptions.defaultptn = TRUE;

    //restore original graph
    for(j = 0; j < ppgraph->multiEdgeCount; j++){
        int n1 = multiEdgeList[j][0];
        int n2 = multiEdgeList[j][1];
        set *v1, *v2;
        v1 = GRAPHROW(ppgraph->ulgraph, n1, MAXM);
        v2 = GRAPHROW(ppgraph->ulgraph, n2, MAXM);
        ADDELEMENT(v1, n2);
        ADDELEMENT(v2, n1);
        DELELEMENT(v1, ppgraph->order + j);
        DELELEMENT(v2, ppgraph->order + j);
    }

    DEBUGMSG("End nauty")

    int *multiEdgeOrbits = globalMultiEdgeOrbits + ((degree2OperationsDepth+1)*HALFFLOOR(vertexCount));
    int *multiEdgeOrbitCount = globalMultiEdgeOrbitCount + degree2OperationsDepth + 1;

    determine_vertex_pairs_orbits(multiEdgeList, *multiEdgeListSize, multiEdgeOrbits, multiEdgeOrbitCount, &automorphismGroupGenerators, numberOfGenerators);

    int i = 0;
    while(i<*multiEdgeListSize && !(multiEdgeList[i][0]==v1 && multiEdgeList[i][1]==v2)) i++;
    DEBUGASSERT(i<*multiEdgeListSize)
    int newEdgeOrbit = multiEdgeOrbits[i]; //contains the number of the orbit of the edge (v1, v2)

    DEBUGARRAYDUMP(nautyLabelling, ppgraph->order, "%d")
    DEBUGARRAYDUMP(nautyPtn, ppgraph->order, "%d")
    int reverseLabelling[ppgraph->order];
    for (i = 0; i < ppgraph->order; i++) {
        reverseLabelling[nautyLabelling[i]]=i;
    }
    DEBUGARRAYDUMP(reverseLabelling, ppgraph->order, "%d")

    VERTEXPAIR smallestRepresentantNewEdge;
    VERTEXPAIR smallestOtherMultiEdge;
    VERTEXPAIR currentEdge;
    smallestRepresentantNewEdge[0]=smallestRepresentantNewEdge[1]=ppgraph->order;
    smallestOtherMultiEdge[0]=smallestOtherMultiEdge[1]=ppgraph->order;

    for (i = 0; i < *multiEdgeListSize; i++) {
        currentEdge[0]=reverseLabelling[multiEdgeList[i][0]];
        currentEdge[1]=reverseLabelling[multiEdgeList[i][1]];
        if(currentEdge[1]<currentEdge[0]){
            int temp = currentEdge[0];
            currentEdge[0] = currentEdge[1];
            currentEdge[1] = temp;
        }
        //check if this is a valid edge
        /*
         *    v
         *    o       ___
         *    |\t   s/   \
         *    | o---o     o---
         *    |/     \___/
         *    o
         *    u
         *
         *  This type of edge on the left side can't be constructed and therefore can't count as valid multi-edge for the construction
         */
        int neighbourV, neighbourU, t, s;
        boolean validEdge = TRUE;
        neighbourV = (ppgraph->adjList[(multiEdgeList[i][0])*3] == multiEdgeList[i][1]) ? ppgraph->adjList[(multiEdgeList[i][0])*3 + 1] : ppgraph->adjList[(multiEdgeList[i][0])*3];
        neighbourU = (ppgraph->adjList[(multiEdgeList[i][1])*3] == multiEdgeList[i][0]) ? ppgraph->adjList[(multiEdgeList[i][1])*3 + 1] : ppgraph->adjList[(multiEdgeList[i][1])*3];
        if(neighbourU == neighbourV){
            t = neighbourU;
            if(ppgraph->adjList[t*3]!=multiEdgeList[i][0] && ppgraph->adjList[t*3]!=multiEdgeList[i][1])
                s = ppgraph->adjList[t*3];
            else if(ppgraph->adjList[t*3+1]!=multiEdgeList[i][0] && ppgraph->adjList[t*3+1]!=multiEdgeList[i][1])
                s = ppgraph->adjList[t*3+1];
            else
                s = ppgraph->adjList[t*3+2];
            validEdge = (ppgraph->degree[s]!=2);
        }

        if(validEdge){
            if(multiEdgeOrbits[i]==newEdgeOrbit){
                if(currentEdge[0]<smallestRepresentantNewEdge[0] ||
                        (currentEdge[0]==smallestRepresentantNewEdge[0] && currentEdge[1] < smallestRepresentantNewEdge[1])){
                    smallestRepresentantNewEdge[0]=currentEdge[0];
                    smallestRepresentantNewEdge[1]=currentEdge[1];
                }
            } else {
                if(currentEdge[0]<smallestOtherMultiEdge[0] ||
                        (currentEdge[0]==smallestOtherMultiEdge[0] && currentEdge[1] < smallestOtherMultiEdge[1])){
                    smallestOtherMultiEdge[0]=currentEdge[0];
                    smallestOtherMultiEdge[1]=currentEdge[1];
                }
            }
        }
    }
    DEBUGDUMP(smallestRepresentantNewEdge[0], "%d")
    DEBUGDUMP(smallestRepresentantNewEdge[1], "%d")
    DEBUGDUMP(smallestOtherMultiEdge[0], "%d")
    DEBUGDUMP(smallestOtherMultiEdge[1], "%d")

    DEBUGMSG("End isCanonicalMultiEdge")
    *multiEdgesDetermined = TRUE;
    if(noRejections) return TRUE;
    return smallestRepresentantNewEdge[0] < smallestOtherMultiEdge[0] ||
            (smallestRepresentantNewEdge[0] == smallestOtherMultiEdge[0] && smallestRepresentantNewEdge[1] < smallestOtherMultiEdge[1]);
}

void handle_deg2_operation1(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators){
    if(operation21Disabled) return;
    DEBUGMSG("Start handle_deg2_operation1")
    DEBUG2DARRAYDUMP((*currentGenerators), currentNumberOfGenerators, ppgraph->order, "%d")
    int maxSize = ppgraph->order*3/2-ppgraph->degree1Count; //this upper bound is not tight (it is tight in case of no degree 2 vertices?)
    VERTEXPAIR edgeList[maxSize]; //initialize an array that is large enough to hold all single edges
    int listSize;
    get_single_edges(ppgraph, edgeList, &listSize);

    int orbitCount;
    int orbits[listSize];
    determine_vertex_pairs_orbits(edgeList, listSize, orbits, &orbitCount, currentGenerators, currentNumberOfGenerators);

    int i;
    for (i = 0; i < listSize; i++) {
        if(orbits[i]==i){
            DEBUGPPGRAPHPRINT(ppgraph)
            apply_deg2_operation1(ppgraph, edgeList[i][0], edgeList[i][1]);

            boolean newMultiEdgesDetermined = FALSE;
            if(isCanonicalMultiEdge(ppgraph, ppgraph->order-2, ppgraph->order-1, &newMultiEdgesDetermined))
                handle_deg2_operation_result(ppgraph, newMultiEdgesDetermined);

            revert_deg2_operation1(ppgraph, edgeList[i][0], edgeList[i][1]);
            DEBUGPPGRAPHPRINT(ppgraph)
        }
    }
    DEBUGMSG("End handle_deg2_operation1")
}

void handle_deg2_operation2(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators,
        boolean *multiEdgesDetermined){
    if(operation22Disabled) return;
    DEBUGMSG("Start handle_deg2_operation2")
    DEBUG2DARRAYDUMP((*currentGenerators), currentNumberOfGenerators, ppgraph->order, "%d")

    VERTEXPAIR *oldMultiEdgeList = globalMultiEdgeList + (degree2OperationsDepth*HALFFLOOR(vertexCount));
    int *oldMultiEdgeListSize = globalMultiEdgeListSize + degree2OperationsDepth;
    int *oldMultiEdgeOrbits = globalMultiEdgeOrbits + (degree2OperationsDepth*HALFFLOOR(vertexCount));
    int *oldMultiEdgeOrbitCount = globalMultiEdgeOrbitCount + degree2OperationsDepth;

    if(!(*multiEdgesDetermined)){
        //if the multi-edges have't been determined, do it now and store the result
        get_multi_edges(ppgraph, oldMultiEdgeList, oldMultiEdgeListSize);

        DEBUGASSERT(*oldMultiEdgeListSize == ppgraph->multiEdgeCount)

        DEBUGMSG("Start nauty")
        int orbits[ppgraph->order + ppgraph->multiEdgeCount]; //the graph is enlarged so we have to provide a large enough array
        numberOfGenerators = 0; //reset the generators
        nautyOptions.defaultptn = FALSE; //use colourings for multigraphs

        int j;
        for(j = 0; j < ppgraph->order + ppgraph->multiEdgeCount; j++){
            nautyLabelling[j]=j;
            nautyPtn[j]=(j!=ppgraph->order-1);
        }
        for(j = 0; j < ppgraph->multiEdgeCount; j++){
            int n1 = oldMultiEdgeList[j][0];
            int n2 = oldMultiEdgeList[j][1];
            set *multiEdgeVertex, *v1, *v2;
            multiEdgeVertex = GRAPHROW(ppgraph->ulgraph, ppgraph->order + j, MAXM);
            v1 = GRAPHROW(ppgraph->ulgraph, n1, MAXM);
            v2 = GRAPHROW(ppgraph->ulgraph, n2, MAXM);
            EMPTYSET(multiEdgeVertex, MAXM);
            DELELEMENT(v1, n2);
            DELELEMENT(v2, n1);
            ADDELEMENT(v1, ppgraph->order + j);
            ADDELEMENT(v2, ppgraph->order + j);
            ADDELEMENT(multiEdgeVertex, n1);
            ADDELEMENT(multiEdgeVertex, n2);
        }

        nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, orbits, &nautyOptions, &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order + ppgraph->multiEdgeCount, canonicalGraph);
        nautyOptions.defaultptn = TRUE;

        //restore original graph
        for(j = 0; j < ppgraph->multiEdgeCount; j++){
            int n1 = oldMultiEdgeList[j][0];
            int n2 = oldMultiEdgeList[j][1];
            set *v1, *v2;
            v1 = GRAPHROW(ppgraph->ulgraph, n1, MAXM);
            v2 = GRAPHROW(ppgraph->ulgraph, n2, MAXM);
            ADDELEMENT(v1, n2);
            ADDELEMENT(v2, n1);
            DELELEMENT(v1, ppgraph->order + j);
            DELELEMENT(v2, ppgraph->order + j);
        }

        DEBUGMSG("End nauty")


        determine_vertex_pairs_orbits(oldMultiEdgeList, *oldMultiEdgeListSize, oldMultiEdgeOrbits, oldMultiEdgeOrbitCount, currentGenerators, currentNumberOfGenerators);
        *multiEdgesDetermined = TRUE;
    }
    int i;
    for (i = 0; i < *oldMultiEdgeListSize; i++) {
        if(oldMultiEdgeOrbits[i]==i){
            DEBUGPPGRAPHPRINT(ppgraph)
            apply_deg2_operation2(ppgraph, oldMultiEdgeList[i][0], oldMultiEdgeList[i][1]);

            boolean newMultiEdgesDetermined = FALSE;
            if(isCanonicalMultiEdge(ppgraph, ppgraph->order-2, ppgraph->order-1, &newMultiEdgesDetermined))
                handle_deg2_operation_result(ppgraph, newMultiEdgesDetermined);

           revert_deg2_operation2(ppgraph, oldMultiEdgeList[i][0], oldMultiEdgeList[i][1]);
            DEBUGPPGRAPHPRINT(ppgraph)
        }
    }
    DEBUGMSG("End handle_deg2_operation2")
}

void handle_deg2_operation3(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators,
        boolean *multiEdgesDetermined){
    if(operation23Disabled) return;
    DEBUGMSG("Start handle_deg2_operation3")
    DEBUG2DARRAYDUMP((*currentGenerators), currentNumberOfGenerators, ppgraph->order, "%d")

    VERTEXPAIR *oldMultiEdgeList = globalMultiEdgeList + (degree2OperationsDepth*HALFFLOOR(vertexCount));
    int *oldMultiEdgeListSize = globalMultiEdgeListSize + degree2OperationsDepth;
    int *oldMultiEdgeOrbits = globalMultiEdgeOrbits + (degree2OperationsDepth*HALFFLOOR(vertexCount));
    int *oldMultiEdgeOrbitCount = globalMultiEdgeOrbitCount + degree2OperationsDepth;

    if(*multiEdgesDetermined){
         //if the multi-edges have't been determined, do it now and store the result
        get_multi_edges(ppgraph, oldMultiEdgeList, oldMultiEdgeListSize);

        DEBUGASSERT(*oldMultiEdgeListSize == ppgraph->multiEdgeCount)

        DEBUGMSG("Start nauty")
        int orbits[ppgraph->order + ppgraph->multiEdgeCount]; //the graph is enlarged so we have to provide a large enough array
        numberOfGenerators = 0; //reset the generators
        nautyOptions.defaultptn = FALSE; //use colourings for multigraphs

        int j;
        for(j = 0; j < ppgraph->order + ppgraph->multiEdgeCount; j++){
            nautyLabelling[j]=j;
            nautyPtn[j]=(j!=ppgraph->order-1);
        }
        for(j = 0; j < ppgraph->multiEdgeCount; j++){
            int n1 = oldMultiEdgeList[j][0];
            int n2 = oldMultiEdgeList[j][1];
            set *multiEdgeVertex, *v1, *v2;
            multiEdgeVertex = GRAPHROW(ppgraph->ulgraph, ppgraph->order + j, MAXM);
            v1 = GRAPHROW(ppgraph->ulgraph, n1, MAXM);
            v2 = GRAPHROW(ppgraph->ulgraph, n2, MAXM);
            EMPTYSET(multiEdgeVertex, MAXM);
            DELELEMENT(v1, n2);
            DELELEMENT(v2, n1);
            ADDELEMENT(v1, ppgraph->order + j);
            ADDELEMENT(v2, ppgraph->order + j);
            ADDELEMENT(multiEdgeVertex, n1);
            ADDELEMENT(multiEdgeVertex, n2);
        }

        nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, orbits, &nautyOptions, &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order + ppgraph->multiEdgeCount, canonicalGraph);
        nautyOptions.defaultptn = TRUE;

        //restore original graph
        for(j = 0; j < ppgraph->multiEdgeCount; j++){
            int n1 = oldMultiEdgeList[j][0];
            int n2 = oldMultiEdgeList[j][1];
            set *v1, *v2;
            v1 = GRAPHROW(ppgraph->ulgraph, n1, MAXM);
            v2 = GRAPHROW(ppgraph->ulgraph, n2, MAXM);
            ADDELEMENT(v1, n2);
            ADDELEMENT(v2, n1);
            DELELEMENT(v1, ppgraph->order + j);
            DELELEMENT(v2, ppgraph->order + j);
        }

        DEBUGMSG("End nauty")

        determine_vertex_pairs_orbits(oldMultiEdgeList, *oldMultiEdgeListSize, oldMultiEdgeOrbits, oldMultiEdgeOrbitCount, currentGenerators, currentNumberOfGenerators);
        *multiEdgesDetermined = TRUE;
    }
    int i;
    for (i = 0; i < *oldMultiEdgeListSize; i++) {
        if(oldMultiEdgeOrbits[i]==i){
            //check if the second neighbour of u and v aren't equal
            int x, y;
            x = (ppgraph->adjList[oldMultiEdgeList[i][0]*3] == oldMultiEdgeList[i][1]) ?
                ppgraph->adjList[oldMultiEdgeList[i][0]*3+1] :
                ppgraph->adjList[oldMultiEdgeList[i][0]*3];
            y = (ppgraph->adjList[oldMultiEdgeList[i][1]*3] == oldMultiEdgeList[i][0]) ?
                ppgraph->adjList[oldMultiEdgeList[i][1]*3+1] :
                ppgraph->adjList[oldMultiEdgeList[i][1]*3];
            if(x!=y){
                DEBUGPPGRAPHPRINT(ppgraph)
                apply_deg2_operation3(ppgraph, oldMultiEdgeList[i][0], oldMultiEdgeList[i][1]);

                boolean newMultiEdgesDetermined = FALSE;
                if(isCanonicalMultiEdge(ppgraph, ppgraph->order-2, ppgraph->order-1, &newMultiEdgesDetermined))
                    handle_deg2_operation_result(ppgraph, newMultiEdgesDetermined);

                revert_deg2_operation3(ppgraph, oldMultiEdgeList[i][0], oldMultiEdgeList[i][1]);
                DEBUGPPGRAPHPRINT(ppgraph)
            }
        }
    }
    DEBUGMSG("End handle_deg2_operation3")
}

/*
 * Performs the different degree 1 operations. When this method returns &ppgraph
 * will be unchanged.
 */
void do_deg1_operations(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators){
    DEBUGMSG("Start do_deg1_operations")
    DEBUGASSERT(allowLoops || allowSemiEdges)
    if(ppgraph->order<=maxVertexCount) handle_deg1_operation1(ppgraph, currentGenerators, currentNumberOfGenerators);
    if(ppgraph->order<=maxVertexCount-2) handle_deg1_operation2(ppgraph, currentGenerators, currentNumberOfGenerators);
    DEBUGMSG("End do_deg1_operations")
}

/*
 * Performs the different degree 2 operations. When this method returns &ppgraph
 * will be unchanged.
 */
void do_deg2_operations(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators, boolean multiEdgesDetermined){
    DEBUGMSG("Start do_deg2_operations")
    DEBUGASSERT(allowMultiEdges)

    boolean stillPossible = ppgraph->order<=maxVertexCount-2;
    if(stillPossible && !allowSemiEdges){ //the degree 2 operations preserve the parity
        stillPossible = (ppgraph->order%2 == vertexCount%2);
    }
    if(stillPossible){ //if we turn all the degree 1 vertices into semi-edges and we still have too much vertices, then we can stop
        stillPossible = ppgraph->order - ppgraph->degree1Count + 2 <= vertexCount;
    }
    if(stillPossible){
        boolean newMultiEdgesDetermined = multiEdgesDetermined;
        handle_deg2_operation1(ppgraph, currentGenerators, currentNumberOfGenerators);
        //the orbits of the multi-edges are already determined, so we pass them on
        handle_deg2_operation2(ppgraph, currentGenerators, currentNumberOfGenerators, &newMultiEdgesDetermined);
        handle_deg2_operation3(ppgraph, currentGenerators, currentNumberOfGenerators, &newMultiEdgesDetermined);
    }
    DEBUGMSG("End do_deg2_operations")
}

void grow(PRIMPREGRAPH *ppgraph){
    DEBUGMSG("Start grow")
    int orbits[ppgraph->order];
    DEBUGMSG("Start nauty")
    numberOfGenerators = 0; //reset the generators
    //there are no multiedges at this point!
    nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, orbits, &nautyOptions, &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order, canonicalGraph);
    DEBUGMSG("End nauty")
    //the generators for these start graphs need to be calculated
    permutation currentGenerators[MAXN][MAXN]; //TODO: can't we make this smaller because we now the size at this point
    int currentNumberOfGenerators = numberOfGenerators;
    copyGenerators(&currentGenerators, ppgraph->order);

    #ifdef _DEBUG
    //check that the generators were copied correctly
    int i, j;
    for(i=0; i<numberOfGenerators; i++){
        for (j = 0; j < ppgraph->order; j++) {
            DEBUGASSERT(currentGenerators[i][j]==automorphismGroupGenerators[i][j])
        }
    }
    DEBUG2DARRAYDUMP(currentGenerators, numberOfGenerators, ppgraph->order, "%d")
    #endif

    if(ppgraph->order >= minVertexCount && ppgraph->order<=maxVertexCount && ppgraph->order - vertexCount <= ppgraph->degree1Count)
        handle_primpregraph_result(ppgraph);


    if(allowLoops || allowSemiEdges){
        do_deg1_operations(ppgraph, &currentGenerators, currentNumberOfGenerators);
    }
    if(allowMultiEdges){
        do_deg2_operations(ppgraph, &currentGenerators, currentNumberOfGenerators, FALSE);
    }
    DEBUGMSG("End grow")
}

void growWithoutDeg1Operations(PRIMPREGRAPH *ppgraph){
    DEBUGMSG("Start growWithoutDeg1Operations")
    int orbits[ppgraph->order];

    VERTEXPAIR *multiEdgeList = globalMultiEdgeList; //depth = 0
    int *multiEdgeListSize = globalMultiEdgeListSize;
    get_multi_edges(ppgraph, multiEdgeList, multiEdgeListSize);

    DEBUGASSERT(*multiEdgeListSize == ppgraph->multiEdgeCount)

    DEBUGMSG("Start nauty")
    numberOfGenerators = 0; //reset the generators
    nautyOptions.defaultptn = FALSE; //use colourings for multigraphs

    int k;
    for(k = 0; k < ppgraph->order + ppgraph->multiEdgeCount; k++){
        nautyLabelling[k]=k;
        nautyPtn[k]=(k!=ppgraph->order-1);
    }
    for(k = 0; k < ppgraph->multiEdgeCount; k++){
        int n1 = multiEdgeList[k][0];
        int n2 = multiEdgeList[k][1];
        set *multiEdgeVertex, *v1, *v2;
        multiEdgeVertex = GRAPHROW(ppgraph->ulgraph, ppgraph->order + k, MAXM);
        v1 = GRAPHROW(ppgraph->ulgraph, n1, MAXM);
        v2 = GRAPHROW(ppgraph->ulgraph, n2, MAXM);
        EMPTYSET(multiEdgeVertex, MAXM);
        DELELEMENT(v1, n2);
        DELELEMENT(v2, n1);
        ADDELEMENT(v1, ppgraph->order + k);
        ADDELEMENT(v2, ppgraph->order + k);
        ADDELEMENT(multiEdgeVertex, n1);
        ADDELEMENT(multiEdgeVertex, n2);
    }

    nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, orbits, &nautyOptions, &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order + ppgraph->multiEdgeCount, canonicalGraph);

    nautyOptions.defaultptn = TRUE;

    //restore original graph
    for(k = 0; k < ppgraph->multiEdgeCount; k++){
        int n1 = multiEdgeList[k][0];
        int n2 = multiEdgeList[k][1];
        set *v1, *v2;
        v1 = GRAPHROW(ppgraph->ulgraph, n1, MAXM);
        v2 = GRAPHROW(ppgraph->ulgraph, n2, MAXM);
        ADDELEMENT(v1, n2);
        ADDELEMENT(v2, n1);
        DELELEMENT(v1, ppgraph->order + k);
        DELELEMENT(v2, ppgraph->order + k);
    }

    DEBUGMSG("End nauty")
    //the generators for these start graphs need to be calculated
    permutation currentGenerators[MAXN][MAXN]; //TODO: can't we make this smaller because we now the size at this point
    int currentNumberOfGenerators = numberOfGenerators;
    copyGenerators(&currentGenerators, ppgraph->order);

    int *multiEdgeOrbits = globalMultiEdgeOrbits;
    int *multiEdgeOrbitCount = globalMultiEdgeOrbitCount;
    determine_vertex_pairs_orbits(multiEdgeList, *multiEdgeListSize, multiEdgeOrbits, multiEdgeOrbitCount, &currentGenerators, currentNumberOfGenerators);

    #ifdef _DEBUG
    //check that the generators were copied correctly
    int i, j;
    for(i=0; i<numberOfGenerators; i++){
        for (j = 0; j < ppgraph->order; j++) {
            DEBUGASSERT(currentGenerators[i][j]==automorphismGroupGenerators[i][j])
        }
    }
    DEBUG2DARRAYDUMP(currentGenerators, numberOfGenerators, ppgraph->order, "%d")
    #endif

    if(ppgraph->order >= minVertexCount && ppgraph->order<=maxVertexCount && ppgraph->order - vertexCount <= ppgraph->degree1Count)
        handle_primpregraph_result(ppgraph);

    if(allowMultiEdges){
        do_deg2_operations(ppgraph, &currentGenerators, currentNumberOfGenerators, TRUE);
    }
    DEBUGMSG("End growWithoutDeg1Operations")
}

static int current3RegOrder;
static PRIMPREGRAPH *currentPpgraph;

/*
 * Handles the input from the external graph by translating the graph into a pregraph primitive
 * and feeding it to the grow method.
 */
void handle_3_regular_result(int *adjacencyList){
    DEBUGMSG("Start handle_3_regular_result")
    PRIMPREGRAPH *ppgraph = currentPpgraph;
    ppgraph->order = current3RegOrder;
    ppgraph->degree1Count = 0;
    ppgraph->multiEdgeCount = 0;

    int i, j;

    for(i=0; i<current3RegOrder; i++){
        ppgraph->degree[i]=3;

        set *v;
        v = GRAPHROW(ppgraph->ulgraph, i, MAXM);
        EMPTYSET(v, MAXM);
        j=0;
        for(j=0; j<3; j++){
            ppgraph->adjList[i*3+j] = adjacencyList[i*3+j];
            ADDELEMENT(v, adjacencyList[i*3+j]);
        }
    }

    DEBUGPPGRAPHPRINT(ppgraph)

    grow(ppgraph);
    DEBUGMSG("End handle_3_regular_result")
}

void start(){
    DEBUGMSG("Start start")
    structureCount=0;
    primitivesCount=0;
    if(!allowSemiEdges){
        minVertexCount = maxVertexCount = vertexCount;
    } else {
        minVertexCount = vertexCount;
        maxVertexCount = 2*vertexCount+2;
    }
    DEBUGDUMP(minVertexCount, "%d")
    DEBUGDUMP(maxVertexCount, "%d")
    PRIMPREGRAPH ppgraph;
    if(allowLoops || allowSemiEdges){
        construct_K2(&ppgraph);
        grow(&ppgraph);
    }
    if(allowMultiEdges){
        construct_C4(&ppgraph);
        growWithoutDeg1Operations(&ppgraph);
    }
    if((allowLoops || allowSemiEdges) && allowMultiEdges){
        construct_K3_with_spike(&ppgraph);
        growWithoutDeg1Operations(&ppgraph);
    }

    int i;
    currentPpgraph = &ppgraph;
    for(i = 4; i <= vertexCount; i+=2){//TODO: is this the correct upperbound for i
        #ifdef _DEBUG
        fprintf(stderr, "Starting snarkhunter for %d vertices\n", i);
        #endif
        current3RegOrder = i;
        //TODO: call into snarkhunter
        init_irreducible_graphs(i);
    }

    if(allowMultiEdges && vertexCount == 2){
        writeFatK2();
    }
    DEBUGMSG("End start")
}

/* For debugging purposes*/
void start2(){
    DEBUGMSG("Start start2")
    structureCount=0;
    primitivesCount=0;
    if(!allowSemiEdges){
        minVertexCount = maxVertexCount = vertexCount;
    } else {
        minVertexCount = vertexCount;
        maxVertexCount = 2*vertexCount+2;
    }
    DEBUGDUMP(minVertexCount, "%d")
    DEBUGDUMP(maxVertexCount, "%d")
    PRIMPREGRAPH ppgraph;
    construct_K3_3(&ppgraph);
    grow(&ppgraph);
    DEBUGMSG("End start2")
}

void startFromFile(FILE *inputFile){
    DEBUGMSG("Start startFromFile")
    int i, endian = LITTLE_ENDIAN;
    unsigned char formatIdentifier;

    //initialize some globals
    structureCount=0;
    primitivesCount=0;
    if(!allowSemiEdges){
        minVertexCount = maxVertexCount = vertexCount;
    } else {
        minVertexCount = vertexCount;
        maxVertexCount = 2*vertexCount+2;
    }
    DEBUGDUMP(minVertexCount, "%d")
    DEBUGDUMP(maxVertexCount, "%d")

    //First we read the header and determine the file format
    DEBUGMSG("")
    unsigned char c[1];
    do {
        if (fread(c, sizeof (unsigned char), 1, inputFile) == 0) {
            fprintf(stderr, "Error while reading input file: aborting!\n");
            exit(EXIT_FAILURE);
        }
    } while (isspace(*c));
    DEBUGMSG("")
    if(*c!='>'){
        fprintf(stderr, "First character was %c (%d).\n", *c, *c);
        fprintf(stderr, "File doesn't start with a header: aborting!\n");
        exit(EXIT_FAILURE);
    }
    DEBUGMSG("")
    if (fread(c, sizeof (unsigned char), 1, inputFile) == 0) {
        fprintf(stderr, "Error while reading input file: aborting!\n");
        exit(EXIT_FAILURE);
    }
    DEBUGMSG("")
    if(*c!='>'){
        fprintf(stderr, "Second character was %c(%d).\n", *c, *c);
        fprintf(stderr, "File doesn't start with a header: aborting!\n");
        exit(EXIT_FAILURE);
    }
    DEBUGMSG("")
    if (fread(c, sizeof (unsigned char), 1, inputFile) == 0) {
        fprintf(stderr, "Error while reading input file: aborting!\n");
        exit(EXIT_FAILURE);
    }
    DEBUGMSG("")
    if(*c=='m'){
        //should be multicode
        unsigned char format[10];
        format[0]='m';
        for(i = 1; i < 10; i++){
            if (fread(format + i, sizeof (unsigned char), 1, inputFile) == 0) {
                fprintf(stderr, "Error while reading input file: aborting!\n");
                exit(EXIT_FAILURE);
            }
        }
        if (strncmp((char *)  format, "multi_code", 10) != 0) {
            DEBUGARRAYDUMP(format, 10, "%c")
            fprintf(stderr, "Input file is in an illegal file format (expected multi_code): aborting!\n");
            exit(EXIT_FAILURE);
        }
        formatIdentifier = 'm';
    } else if(*c=='p'){
        //should be pregraphcode
        unsigned char format[13];
        format[0]='p';
        for(i = 1; i < 13; i++){
            if (fread(format + i, sizeof (unsigned char), 1, inputFile) == 0) {
                fprintf(stderr, "Error while reading input file: aborting!\n");
                exit(EXIT_FAILURE);
            }
        }
        if (strncmp((char *)format, "pregraph_code", 13) != 0) {
            DEBUGARRAYDUMP(format, 13, "%c")
            fprintf(stderr, "Input file is in an illegal file format (expected pregraph_code): aborting!\n");
            exit(EXIT_FAILURE);
        }
        formatIdentifier = 'p';

    } else {
        //not supported format
        fprintf(stderr, "Third character was %c(%d).\n", *c, *c);
        fprintf(stderr, "Input file is in an illegal file format: aborting!\n");
        exit(EXIT_FAILURE);
    }
    DEBUGMSG("")

    //next we read the endian of the file
    unsigned char nextTwo[2];
    do {
        if (fread(&nextTwo[0], sizeof (unsigned char), 1, inputFile) == 0) {
            fprintf(stderr, "Error while reading input file: aborting!\n");
            exit(EXIT_FAILURE);
        }
    } while (isspace(nextTwo[0]));
    if (fread(&nextTwo[1], sizeof (unsigned char), 1, inputFile) == 0) {
        fprintf(stderr, "Error while reading input file: aborting!\n");
        exit(EXIT_FAILURE);
    }
    //endian defaults to le
    if (strncmp((char *) nextTwo, "le", 2) == 0) {
        endian = LITTLE_ENDIAN;
    } else if (strncmp((char *) & nextTwo[0], "be", 2) == 0) {
        endian = BIG_ENDIAN;
    }

    //finally we read on to the end of the header
    while (strncmp((char *) nextTwo, "<<", 2) != 0) {
        nextTwo[0] = nextTwo[1];
        if (fread(&nextTwo[1], sizeof (unsigned char), 1, inputFile) == 0) {
            fprintf(stderr, "Error while reading input file: aborting!\n");
            exit(EXIT_FAILURE);
        }
    }

    fprintf(stderr, "Reading file in %s format with %s endian.\n",
            formatIdentifier=='p' ? "pregraph_code" : "multi_code",
            endian == LITTLE_ENDIAN ? "little" : "big");

    //=====================Finished reading header========================

    if(formatIdentifier=='p'){
        PRIMPREGRAPH ppgraph;
        while(readPregraphCode(inputFile, &ppgraph, endian)==(char)1){
            if((!allowLoops && !allowSemiEdges) && ppgraph.degree1Count!=0){
                fprintf(stderr, "Input graph contains vertex of degree 1, but loops and semi-edges aren't allowed.\n");
            }
            if(!allowMultiEdges && ppgraph.multiEdgeCount!=0){
                fprintf(stderr, "Input graph contains vertices of degree 2, but multi-edges aren't allowed.\n");
            }

            if(ppgraph.multiEdgeCount==0){
                grow(&ppgraph);
            } else {
                growWithoutDeg1Operations(&ppgraph);
            }
        }
    } else if(formatIdentifier=='m'){
        fprintf(stderr, "Multi_code not yet implemented.\n");
    }

    DEBUGMSG("End startFromFile")
}

//----------------------Begin input methods------------------------------

char read_2byte_number(FILE *f, unsigned short *n, int endian) {
    unsigned char c[2];
    if (fread(&c[0], sizeof (unsigned char), 2, f) < 2) {
        return (2);
    }
    if (endian == BIG_ENDIAN) {
        *n = c[0]*256 + c[1];
    } else {
        *n = c[1]*256 + c[0];
    }
    return (1);
}

char read_old_or_new(FILE *f, boolean bignum, int endian, unsigned short *number) {
    unsigned char k;
    if (bignum) {
        if (read_2byte_number(f, number, endian) == 2) {
            return (2);
        }
    } else {
        if (fread(&k, sizeof (unsigned char), 1, f) == 0) {
            return (2);
        }
        *number = (unsigned short) k;
    }
    return (1);
}

/*
 * Reads the pregraph primitive of a given graph in pregraphcode.
 * Returns 2 or larger number if an error occurred. Returns 0 if
 * EOF was reached. Returns 1 if all went OK.
 */
char readPregraphCode(FILE *f, PRIMPREGRAPH *ppgraph, int endian) {
    DEBUGMSG("Start readPregraphCode")
    int i, j, n, dummyVertex;
    unsigned short signum, number;
    if (read_old_or_new(f, FALSE, endian, &signum) == 2) {
        DEBUGMSG("End readPregraphCode")
        return (feof(f) ? 0 : 2);
    }
    //if the code starts with a zero, all the entries are two bytes
    //else the number we just read was the order of the graph
    if (signum == 0) {
        if (read_old_or_new(f, TRUE, endian, &number) == 2) {
            DEBUGMSG("End readPregraphCode")
            return (2);
        }
    } else {
        number = signum;
    }

    if ((n = (int) number) > MAXN) {
        DEBUGMSG("End readPregraphCode")
        return (3);
    }

    //initialize the pregraph
    ppgraph->order=n; //will be increased when we find semi-edges
    for(i = 0; i < MAXN; i++){
        ppgraph->degree[i]=0;
    }
    ppgraph->degree1Count=0;
    ppgraph->multiEdgeCount=0;

    i = 1;
    dummyVertex = n;
    while (i <= n) {
        if (read_old_or_new(f, signum == 0, endian, &number) == 2) {
            DEBUGMSG("End readPregraphCode")
            return (2);
        }
        DEBUGDUMP(i, "%d")
        DEBUGDUMP(number, "%d")
        if (number != 0) {
            if(number == i){
                DEBUGMSG("loop")
                //we found a loop
                ppgraph->degree1Count++;
            } else if (number == n+1){
                DEBUGMSG("semi-edge")
                //we found a semi-edge
                ppgraph->degree1Count++;
                ppgraph->order++;
                ppgraph->adjList[(i-1)*3 + ppgraph->degree[i-1]] = dummyVertex;
                ppgraph->degree[i-1]++;
                ppgraph->adjList[dummyVertex*3 + 0] = i-1;
                ppgraph->degree[dummyVertex]=1;
                dummyVertex++;
            } else {
                //first check to see if we have a multi-edge
                j=0;
                while(j<ppgraph->degree[i-1] && ppgraph->adjList[(i-1)*3 + j] != (number-1)) j++;
                DEBUGDUMP(j, "%d")
                DEBUGDUMP(ppgraph->degree[i-1], "%d")
                if(j != ppgraph->degree[i-1]){
                    DEBUGMSG("multi-edge")
                    //we found a multi-edge (do not increase the degree)
                    ppgraph->multiEdgeCount++;
                    ppgraph->multiedge[i-1]=number-1;
                    ppgraph->multiedge[number-1]=i-1;
                } else {
                    DEBUGMSG("simple edge")
                    ppgraph->adjList[(i-1)*3 + ppgraph->degree[i-1]] = number-1;
                    ppgraph->adjList[(number-1)*3 + ppgraph->degree[number-1]] = i-1;
                    ppgraph->degree[i-1]++;
                    ppgraph->degree[number-1]++;
                }
            }
        } else {
            DEBUGMSG("=========next vertex=========")
            i++;
        }
    }

    for(i = 0; i < ppgraph->order; i++){
        set *vertex;
        vertex = GRAPHROW(ppgraph->ulgraph, i, MAXM);
        EMPTYSET(vertex, MAXM);
        for(j = 0; j < ppgraph->degree[i]; j++){
            ADDELEMENT(vertex, ppgraph->adjList[i*3+j]);
        }
    }
    DEBUGMSG("End readPregraphCode")
    return (1);
}

//--------------------Begin extra output graph---------------------------

/*   __
 *  /  \
 * o----o
 *  \__/
 *
 */
void writeFatK2(){
    structureCount++;
    FILE *file = stdout;
    if(outputFile != NULL){
        file = fopen(outputFile, "a");
        DEBUGMSG("Opened file")
    }

    if(outputType == 'h'){
        fprintf(file, "========================================\n");
        fprintf(file, "|  Graph number: %20llu  |\n", structureCount);
        fprintf(file, "|  Number of vertices:              2  |\n");
        fprintf(file, "========================================\n");
        fprintf(file, "|   1 ||    2 |    2 |    2 ||\n");
        fprintf(file, "|   2 ||    1 |    1 |    1 ||\n");
        fprintf(file, "==============================\n\n");
    } else if(outputType == 'p'){
        if (structureCount == 1) { //if first graph
            fprintf(file, ">>pregraph_code %s<<", (endian == LITTLE_ENDIAN ? "le" : "be"));
        }
        fprintf(file, "%c", 2);
        fprintf(file, "%c", 2);
        fprintf(file, "%c", 2);
        fprintf(file, "%c", 2);
        fprintf(file, "%c", 1);
        fprintf(file, "%c", 1);
        fprintf(file, "%c", 1);
        fprintf(file, "%c", 0);
    }

    if(outputFile != NULL){
        fclose(file);
        DEBUGMSG("Closed file")
    }

    if(logStatistics){
        graphsWithOnlyMultiEdgesCount++;
        graphsWithMultiEdgesCount[1]++;
        graphsWithLoopsCount[0]++;
        graphsWithSemiEdgesCount[0]++;
    }
}

//---------------------End extra output graph----------------------------

//------------------------Begin start graphs-----------------------------

/*
 * o----o
 *
 */
void construct_K2(PRIMPREGRAPH *ppgraph){
    ppgraph->order = 2;
    ppgraph->degree1Count = 2;
    ppgraph->multiEdgeCount = 0;
    ppgraph->adjList[0*3] = 1;
    ppgraph->adjList[1*3] = 0;
    ppgraph->degree[0]=1;
    ppgraph->degree[1]=1;

    set *g0, *g1;
    g0 = GRAPHROW(ppgraph->ulgraph, 0, MAXM);
    g1 = GRAPHROW(ppgraph->ulgraph, 1, MAXM);
    EMPTYSET(g0, MAXM);
    EMPTYSET(g1, MAXM);
    ADDELEMENT(g0, 1);
    ADDELEMENT(g1, 0);
    DEBUGMSG("Created K2")
}

/*
 * o---o
 * |   |
 * |   |
 * o---o
 *
 */
void construct_C4(PRIMPREGRAPH *ppgraph){
    ppgraph->order = 4;
    ppgraph->degree1Count = 0;
    ppgraph->multiEdgeCount = 2;
    ppgraph->adjList[0*3] = 3;
    ppgraph->adjList[0*3+1] = 1;
    ppgraph->adjList[1*3] = 0;
    ppgraph->adjList[1*3+1] = 2;
    ppgraph->adjList[2*3] = 1;
    ppgraph->adjList[2*3+1] = 3;
    ppgraph->adjList[3*3] = 2;
    ppgraph->adjList[3*3+1] = 0;
    ppgraph->degree[0]=2;
    ppgraph->degree[1]=2;
    ppgraph->degree[2]=2;
    ppgraph->degree[3]=2;
    ppgraph->multiedge[0]=1;
    ppgraph->multiedge[1]=0;
    ppgraph->multiedge[2]=3;
    ppgraph->multiedge[3]=2;

    set *g0, *g1, *g2, *g3;
    g0 = GRAPHROW(ppgraph->ulgraph, 0, MAXM);
    g1 = GRAPHROW(ppgraph->ulgraph, 1, MAXM);
    g2 = GRAPHROW(ppgraph->ulgraph, 2, MAXM);
    g3 = GRAPHROW(ppgraph->ulgraph, 3, MAXM);
    EMPTYSET(g0, MAXM);
    EMPTYSET(g1, MAXM);
    EMPTYSET(g2, MAXM);
    EMPTYSET(g3, MAXM);
    ADDELEMENT(g0, 3); ADDELEMENT(g0, 1);
    ADDELEMENT(g1, 0); ADDELEMENT(g1, 2);
    ADDELEMENT(g2, 1); ADDELEMENT(g2, 3);
    ADDELEMENT(g3, 2); ADDELEMENT(g3, 0);
    DEBUGMSG("Created C4")
}

/*    o
 *    |
 *    o
 *   / \
 *  o---o
 */
void construct_K3_with_spike(PRIMPREGRAPH *ppgraph){
    ppgraph->order = 4;
    ppgraph->degree1Count = 1;
    ppgraph->multiEdgeCount = 1;
    ppgraph->adjList[0*3] = 1;
    ppgraph->adjList[0*3+1] = 2;
    ppgraph->adjList[0*3+2] = 3;
    ppgraph->adjList[1*3] = 0;
    ppgraph->adjList[1*3+1] = 2;
    ppgraph->adjList[2*3] = 1;
    ppgraph->adjList[2*3+1] = 0;
    ppgraph->adjList[3*3] = 0;
    ppgraph->degree[0]=3;
    ppgraph->degree[1]=2;
    ppgraph->degree[2]=2;
    ppgraph->degree[3]=1;
    ppgraph->multiedge[1]=2;
    ppgraph->multiedge[2]=1;

    set *g0, *g1, *g2, *g3;
    g0 = GRAPHROW(ppgraph->ulgraph, 0, MAXM);
    g1 = GRAPHROW(ppgraph->ulgraph, 1, MAXM);
    g2 = GRAPHROW(ppgraph->ulgraph, 2, MAXM);
    g3 = GRAPHROW(ppgraph->ulgraph, 3, MAXM);
    EMPTYSET(g0, MAXM);
    EMPTYSET(g1, MAXM);
    EMPTYSET(g2, MAXM);
    EMPTYSET(g3, MAXM);
    ADDELEMENT(g0, 1); ADDELEMENT(g0, 2); ADDELEMENT(g0, 3);
    ADDELEMENT(g1, 0); ADDELEMENT(g1, 2);
    ADDELEMENT(g2, 1); ADDELEMENT(g2, 0);
    ADDELEMENT(g3, 0);
    DEBUGMSG("Created K3 with spike")
}

/*  For testing purposes
 */
void construct_K3_3(PRIMPREGRAPH *ppgraph){
    ppgraph->order = 6;
    ppgraph->degree1Count = 0;
    ppgraph->multiEdgeCount = 0;
    ppgraph->adjList[0*3] = 3;
    ppgraph->adjList[0*3+1] = 4;
    ppgraph->adjList[0*3+2] = 5;
    ppgraph->adjList[1*3] = 3;
    ppgraph->adjList[1*3+1] = 4;
    ppgraph->adjList[1*3+2] = 5;
    ppgraph->adjList[2*3] = 3;
    ppgraph->adjList[2*3+1] = 4;
    ppgraph->adjList[2*3+2] = 5;
    ppgraph->adjList[3*3] = 0;
    ppgraph->adjList[3*3+1] = 1;
    ppgraph->adjList[3*3+2] = 2;
    ppgraph->adjList[4*3] = 0;
    ppgraph->adjList[4*3+1] = 1;
    ppgraph->adjList[4*3+2] = 2;
    ppgraph->adjList[5*3] = 0;
    ppgraph->adjList[5*3+1] = 1;
    ppgraph->adjList[5*3+2] = 2;
    ppgraph->degree[0]=3;
    ppgraph->degree[1]=3;
    ppgraph->degree[2]=3;
    ppgraph->degree[3]=3;
    ppgraph->degree[4]=3;
    ppgraph->degree[5]=3;

    set *g0, *g1, *g2, *g3, *g4, *g5;
    g0 = GRAPHROW(ppgraph->ulgraph, 0, MAXM);
    g1 = GRAPHROW(ppgraph->ulgraph, 1, MAXM);
    g2 = GRAPHROW(ppgraph->ulgraph, 2, MAXM);
    g3 = GRAPHROW(ppgraph->ulgraph, 3, MAXM);
    g4 = GRAPHROW(ppgraph->ulgraph, 4, MAXM);
    g5 = GRAPHROW(ppgraph->ulgraph, 5, MAXM);
    EMPTYSET(g0, MAXM);
    EMPTYSET(g1, MAXM);
    EMPTYSET(g2, MAXM);
    EMPTYSET(g3, MAXM);
    EMPTYSET(g4, MAXM);
    EMPTYSET(g5, MAXM);
    ADDELEMENT(g0, 3); ADDELEMENT(g0, 4); ADDELEMENT(g0, 5);
    ADDELEMENT(g1, 3); ADDELEMENT(g1, 4); ADDELEMENT(g1, 5);
    ADDELEMENT(g2, 3); ADDELEMENT(g2, 4); ADDELEMENT(g2, 5);
    ADDELEMENT(g3, 0); ADDELEMENT(g3, 1); ADDELEMENT(g3, 2);
    ADDELEMENT(g4, 0); ADDELEMENT(g4, 1); ADDELEMENT(g4, 2);
    ADDELEMENT(g5, 0); ADDELEMENT(g5, 1); ADDELEMENT(g5, 2);
    DEBUGMSG("Created K3,3")
}

//-------------------------End start graphs------------------------------

//----------------------Begin Nauty interaction--------------------------

void initNautyOptions() {
    //TODO also options without getcanon?
    nautyOptions.getcanon = TRUE;
    nautyOptions.userautomproc = saveGenerators;
    #ifdef _DEBUG
    nautyOptions.writeautoms = TRUE;
    //options.writemarkers = TRUE;
    nautyOptions.outfile = stderr;
    #endif
}

void saveGenerators(int count, permutation perm[], nvector orbits[],
        int numorbits, int stabvertex, int n) {
    memcpy(automorphismGroupGenerators + numberOfGenerators, perm, sizeof(permutation) * n);

    numberOfGenerators++;
}

void copyGenerators(permutation (*copy)[MAXN][MAXN], int n) {
    int i;
    for(i=0; i<numberOfGenerators; i++){
        memcpy((*copy) + i, automorphismGroupGenerators + i, sizeof(permutation) * n);
    }
}

//------------------------End Nauty interaction-------------------------

/*
print a usage message. name is the name of the current program.
 */
void usage(char *name) {
    fprintf(stderr, "Usage: %s [options] n\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

/*
print a help message. name is the name of the current program.
 */
void help(char *name) {
    fprintf(stderr, "The program %s calculates canonical pregraphs of a given order n.\n", name);
    fprintf(stderr, "Usage: %s [options] n \n\n", name);
    fprintf(stderr, "Valid options:\n");
    fprintf(stderr, "  -h          : Print this help and return.\n");
    fprintf(stderr, "  -i          : Causes %s to print extra info about the generated structures.\n", name);
    fprintf(stderr, "  -L          : Allow loops.\n");
    fprintf(stderr, "  -S          : Allow semi-edges.\n");
    fprintf(stderr, "  -M          : Allow multi-edges.\n");
    fprintf(stderr, "  -P          : Only generate the corresponding pregraph primitives.\n");
    fprintf(stderr, "  -f file     : Specifies the output file. If absent, the output is written to standard out.\n");
    fprintf(stderr, "  -o c        : Specifies the export format where c is one of\n");
    fprintf(stderr, "                c    pregraph code (or multicode if -P is used)\n");
    fprintf(stderr, "                h    human-readable output in tabular format\n");
    fprintf(stderr, "                n    no output: only count (default)\n");
    fprintf(stderr, "  -D #        : Disable some operations:\n");
    fprintf(stderr, "                1    Disable operation 1.1\n");
    fprintf(stderr, "                2    Disable operation 1.2\n");
    fprintf(stderr, "                3    Disable operation 2.1\n");
    fprintf(stderr, "                4    Disable operation 2.2\n");
    fprintf(stderr, "                5    Disable operation 2.3\n");
    fprintf(stderr, "  -d #        : Disable some operations:\n");
    fprintf(stderr, "                1    Disable operations for degree 1\n");
    fprintf(stderr, "                2    Disable operations for degree 2\n");
    fprintf(stderr, "  -X          : No operation will be discarded as being not-canonical. This will cause\n");
    fprintf(stderr, "                isomorphic graphs to be constructed. (This option is for debugging purposes.)\n");
}

void initInfo(){
    int i;
    //initialize arrays and variables
    graphsWithLoopsCount = (unsigned long long *)malloc(sizeof(unsigned long long)*(vertexCount + 1));
    for(i = 0; i <= vertexCount; i++) graphsWithLoopsCount[i]=0;

    graphsWithSemiEdgesCount = (unsigned long long *)malloc(sizeof(unsigned long long)*(vertexCount+2 + 1));
    for(i = 0; i <= vertexCount + 2; i++) graphsWithSemiEdgesCount[i]=0;

    graphsWithMultiEdgesCount = (unsigned long long *)malloc(sizeof(unsigned long long)*(vertexCount/2 + 1));
    for(i = 0; i <= vertexCount/2; i++) graphsWithMultiEdgesCount[i]=0;

    graphsWithOnlyLoopsCount = 0;
    graphsWithOnlySemiEdgesCount = 0;
    graphsWithOnlyMultiEdgesCount = 0;
    simplegraphsCount = 0;
}

void logInfo(PREGRAPH *pregraph){
    int multiEdgeCount, loopCount, semiEdgeCount;
    multiEdgeCount = pregraph->ppgraph->multiEdgeCount;
    semiEdgeCount = pregraph->ppgraph->order - pregraph->order;
    loopCount = pregraph->ppgraph->degree1Count - semiEdgeCount;
    DEBUGASSERT(multiEdgeCount<=vertexCount/2)
    DEBUGASSERT(semiEdgeCount<=vertexCount+2)
    DEBUGASSERT(loopCount<=vertexCount)

    graphsWithLoopsCount[loopCount]++;
    graphsWithSemiEdgesCount[semiEdgeCount]++;
    graphsWithMultiEdgesCount[multiEdgeCount]++;

    if(multiEdgeCount == 0 && semiEdgeCount == 0 && loopCount == 0) simplegraphsCount++;
    if(multiEdgeCount == 0 && semiEdgeCount == 0 && loopCount != 0) graphsWithOnlyLoopsCount++;
    if(multiEdgeCount == 0 && semiEdgeCount != 0 && loopCount == 0) graphsWithOnlySemiEdgesCount++;
    if(multiEdgeCount != 0 && semiEdgeCount == 0 && loopCount == 0) graphsWithOnlyMultiEdgesCount++;
}

void printInfo(){
    if(!onlyPrimitives){
        int i;
        boolean printSpacer = FALSE;
        for(i = 0; i <= vertexCount; i++){
            if(graphsWithLoopsCount[i]!=0){
                fprintf(stderr, "Generated %llu graphs with %d loop%s.\n", graphsWithLoopsCount[i], i, i==1 ? (char *)"" : (char *)"s");
                printSpacer = TRUE;
            }
        }
        if(printSpacer){
            fprintf(stderr, "\n");
            printSpacer = FALSE;
        }
        for(i = 0; i <= vertexCount + 2; i++){
            if(graphsWithSemiEdgesCount[i]!=0){
                fprintf(stderr, "Generated %llu graphs with %d semi-edge%s.\n", graphsWithSemiEdgesCount[i], i, i==1 ? (char *)"" : (char *)"s");
                printSpacer = TRUE;
            }
        }
        if(printSpacer){
            fprintf(stderr, "\n");
            printSpacer = FALSE;
        }
        for(i = 0; i <= vertexCount/2; i++){
            if(graphsWithMultiEdgesCount[i]!=0){
                fprintf(stderr, "Generated %llu graphs with %d multi-edge%s.\n", graphsWithMultiEdgesCount[i], i, i==1 ? (char *)"" : (char *)"s");
                printSpacer = TRUE;
            }
        }
        if(printSpacer){
            fprintf(stderr, "\n");
            printSpacer = FALSE;
        }

        fprintf(stderr, "Generated %llu simple graph%s.\n", simplegraphsCount, simplegraphsCount==1 ? (char *)"" : (char *)"s");
        fprintf(stderr, "Generated %llu graph%s with only loops (and at least one loop).\n", graphsWithOnlyLoopsCount, graphsWithOnlyLoopsCount==1 ? (char *)"" : (char *)"s");
        fprintf(stderr, "Generated %llu graph%s with only semi-edges (and at least one semi-edge).\n", graphsWithOnlySemiEdgesCount, graphsWithOnlySemiEdgesCount==1 ? (char *)"" : (char *)"s");
        fprintf(stderr, "Generated %llu graph%s with only multi-edges (and at least one multi-edge).\n", graphsWithOnlyMultiEdgesCount, graphsWithOnlyMultiEdgesCount==1 ? (char *)"" : (char *)"s");
    }
    fprintf(stderr, "\nGenerated %llu pregraph primitive%s.\n", primitivesCount, primitivesCount==1 ? (char *)"" : (char *)"s");

    fprintf(stderr, "\nDegree 1 operations maximum recursion depth: %d.\n", degree1OperationsDepthMaximum);
    fprintf(stderr, "Degree 2 operations maximum recursion depth: %d.\n", degree2OperationsDepthMaximum);
}

#ifdef PREGRAPH_NO_MAIN
    #define PREGRAPH_MAIN_FUNCTION pregraphnomain
#else
    #define PREGRAPH_MAIN_FUNCTION main
#endif
/*
 *
 */
int PREGRAPH_MAIN_FUNCTION(int argc, char** argv) {

    /*=========== commandline parsing ===========*/

    int c;
    char *name = argv[0];
    char *disabled;
    boolean fromFile = FALSE;
    char *inputFileName = NULL;
    FILE *inputFile = stdin;

    while ((c = getopt(argc, argv, "LSMPXf:F:o:D:d:Ihi")) != -1) {
        switch (c) {
            case 'L': //(defaults to FALSE)
                allowLoops = TRUE;
                break;
            case 'S': //(defaults to FALSE)
                allowSemiEdges = TRUE;
                break;
            case 'M': //(defaults to FALSE)
                allowMultiEdges = TRUE;
                break;
            case 'P': //(defaults to FALSE)
                onlyPrimitives = TRUE;
                break;
            case 'X': //(defaults to FALSE)
                noRejections = TRUE;
                break;
            case 'f': //(defaults to stdout)
                outputFile = optarg;
                break;
            case 'F': //(defaults to stdin)
                inputFileName = optarg;
                break;
            case 'o':
                outputType = optarg[0];
                switch (outputType) {
                    case 'n': //no output (default)
                    case 'c': //pregraph code or multicode
                    case 'h': //human-readable
                        break;
                    default:
                        fprintf(stderr, "Illegal output format %c.\n", c);
                        usage(name);
                        return 1;
                }
                break;
            case 'D':
                //disable certain operations
                disabled = optarg;
                while(*disabled != 0){
                    switch (*disabled) {
                        case '1': //disable operation 1.1
                            operation11Disabled = TRUE;
                            break;
                        case '2': //disable operation 1.2
                            operation12Disabled = TRUE;
                            break;
                        case '3': //disable operation 2.1
                            operation21Disabled = TRUE;
                            break;
                        case '4': //disable operation 2.2
                            operation22Disabled = TRUE;
                            break;
                        case '5': //disable operation 2.3
                            operation23Disabled = TRUE;
                            break;
                        default:
                            fprintf(stderr, "Illegal parameter %c for option -D.\n", *disabled);
                            usage(name);
                            return 1;
                    }
                    disabled++;
                }
                break;
            case 'd':
                //disable certain operations
                disabled = optarg;
                while(*disabled != 0){
                    switch (*disabled) {
                        case '1': //disable operations for degree 1
                            operation11Disabled = TRUE;
                            operation12Disabled = TRUE;
                            break;
                        case '2': //disable operations for degree 2
                            operation21Disabled = TRUE;
                            operation22Disabled = TRUE;
                            operation23Disabled = TRUE;
                            break;
                        default:
                            fprintf(stderr, "Illegal parameter %c for option -d.\n", *disabled);
                            usage(name);
                            return 1;
                    }
                    disabled++;
                }
                break;
            case 'I':
                fromFile = TRUE;
                break;
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case 'i':
                logStatistics = TRUE;
                break;
            default:
                fprintf(stderr, "Illegal option %c.\n", c);
                usage(name);
                return EXIT_FAILURE;
        }
    }

    // check the non-option arguments
    if (argc - optind != 1) {
        usage(name);
        return EXIT_FAILURE;
    }

    //parse the order
    vertexCount = strtol(argv[optind], NULL, 10);
    DEBUGDUMP(vertexCount, "%d")

    if(vertexCount > MAXN || (allowSemiEdges && 2*vertexCount + 2 > MAXN)
            || (allowMultiEdges && vertexCount + vertexCount/2 > MAXN)) {
        fprintf(stderr, "%s needs to be recompiled to support graphs of order %d and with the given option.\n", name, vertexCount);
        return EXIT_FAILURE;
    }

    /*=========== initialization ===========*/

    if(logStatistics) initInfo();

    if(outputFile != NULL){
        FILE *file = fopen(outputFile, "r");
        if(file != NULL){
            fprintf(stderr, "File %s already exists: aborting!\n", outputFile);
            return EXIT_FAILURE;
        }
    }

    if(fromFile && inputFileName != NULL){
        inputFile = fopen(inputFileName, "r");
        if(inputFile == NULL){
            fprintf(stderr, "File %s doesn't exist: aborting!\n", inputFileName);
            return EXIT_FAILURE;
        }
    }

    #ifdef _DEBUG
        fprintf(stderr, "%s:%u MAXN: %d and MAXM: %d \n", __FILE__, __LINE__, MAXN, MAXM);
        fprintf(stderr, "%s:%u Wordsize: %d \n", __FILE__, __LINE__, WORDSIZE);
        fprintf(stderr, "%s:%u Default workspacesize %d \n", __FILE__, __LINE__, NAUTY_WORKSIZE);
    #endif

    nauty_check(WORDSIZE, 1, 30, NAUTYVERSIONID);

    initNautyOptions();

    //create the memory used for storing the orbits of multi-edges

    //first we calculate the theoretical recursion depth of the degree 2 operations
    int depth = 0;
    if(allowMultiEdges){
        if(allowSemiEdges){
            depth = vertexCount%2==0 ? vertexCount/2 : (vertexCount-1)/2;
        } else if(allowLoops){
            depth = vertexCount/2-1;
        } else {
            depth = vertexCount/2-2;
        }
        if(depth < 0) depth = 0;
    }

    //create the different arrays and store the pointers
    //TODO: currently *2 because depth was wrong. It was too small.
    VERTEXPAIR multiEdgeList[2*depth*HALFFLOOR(vertexCount)];
    int multiEdgeListSize[2*depth];
    int multiEdgeOrbits[2*depth*HALFFLOOR(vertexCount)];
    int multiEdgeOrbitCount[2*depth];

    globalMultiEdgeList = multiEdgeList;
    globalMultiEdgeListSize = multiEdgeListSize;
    globalMultiEdgeOrbits = multiEdgeOrbits;
    globalMultiEdgeOrbitCount = multiEdgeOrbitCount;

    /*=========== generation ===========*/

    fprintf(stderr, "Generating %s%s with %d %s%s%s%s.\n", onlyPrimitives ? (char *)"pregraph primitives of " : (char *)"" , allowMultiEdges ? (char *)"multigraphs" : (char *)"simple graphs",
            vertexCount, vertexCount==1 ? (char *)"vertex" : (char *)"vertices",
            allowLoops && allowSemiEdges ? (char *)", " : (allowLoops ? (char *)" and " : (char *)""),
            allowLoops ? (char *)"loops" : (char *)"", allowSemiEdges ? (char *)" and semi-edges" : (char *)"");

    if(!allowSemiEdges && vertexCount%2==1){
        structureCount = primitivesCount = 0;
    } else {
        if(fromFile){
            startFromFile(inputFile);
        } else {
            start();
        }
    }

    if(!onlyPrimitives)
        fprintf(stderr, "Found %llu pregraph%s with %d %s.\n", structureCount, structureCount==1 ? (char *)"" : (char *)"s",
                                                            vertexCount, vertexCount==1 ? (char *)"vertex" : (char *)"vertices");

    if(logStatistics || onlyPrimitives ) printInfo();

    return EXIT_SUCCESS;
}

