/*
 * File:   pregraphs.c
 * Author: nvcleemp
 *
 * Created on December 8, 2009, 2:32 PM
 */

//#define _TEST
//#define _DEBUG

#include "pregraphs.h"

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
    fprintf(f, "|  Graph number: %10ld  |\n", structureCount);
    fprintf(f, "|  Number of vertices: %4d  |\n", pregraph->order);
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
                if (pregraph->order + 1 <= UCHAR_MAX) {
                    fprintf(f, "%c", (unsigned char) primPregraph2Pregraph[ppgraph->adjList[i * 3 + j]]);
                } else {
                    if (write_2byte_number(f, primPregraph2Pregraph[ppgraph->adjList[i * 3 + j]], endian) == 2) {
                        return (2);
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
            } else if (j == 2) {
                //add multi-edge
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
                fprintf(stderr, "Error while writing graph %ld\n", structureCount);
                exit(1);
            }
        } else if(outputType == 'p'){
            if (writePregraphCode(file, pregraph) == 2) {
                fprintf(stderr, "Error while writing graph %ld\n", structureCount);
                exit(1);
            }
        }

        if(outputFile != NULL){
            fclose(file);
            DEBUGMSG("Closed file")
        }
    }
    DEBUGMSG("End handle_pregraph_result")
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
        determine_possible_sets_of_degree1_vertices(tempSet, vertexSetList, &position, semiEdgeCount, 0, nextDegree1Vertex(-1, ppgraph), ppgraph);
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
        int currentSetSize, int currentSetElement, PRIMPREGRAPH *ppgraph){
    DEBUGMSG("Start determine_possible_sets_of_degree1_vertices")
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
                    currentSetSize + 1, nextDegree1Vertex(currentSetElement, ppgraph), ppgraph);
        }
        DELELEMENT(tempSet, currentSetElement);
        DEBUGDUMP(currentSetElement, "%d removed")
        determine_possible_sets_of_degree1_vertices
               (tempSet, vertexSetList, currentListPosition, maximumSetSize,
                currentSetSize, nextDegree1Vertex(currentSetElement, ppgraph), ppgraph);
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
    //if correct number of vertices
    if(ppgraph->order >= minVertexCount && ppgraph->order<=maxVertexCount && ppgraph->order - vertexCount <= ppgraph->degree1Count)
        handle_primpregraph_result(ppgraph);

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

    if(allowMultiEdges){
        do_deg2_operations(ppgraph, &currentGenerators, currentNumberOfGenerators, NULL, 0, NULL, 0);
    }
    DEBUGMSG("End handle_deg1_operation_result")
}

/*
 * Handles the result of a degree 2 operation. If the graph has a valid size it is
 * passed on to handle_primpregraph_result, then we continue with degree 2 operations.
 * Degree 1 operations are no longer possible.
 */
void handle_deg2_operation_result(PRIMPREGRAPH *ppgraph,
        VERTEXPAIR *multiEdgeList, int multiEdgeListSize, int *multiEdgeOrbits, int multiEdgeOrbitCount){
    DEBUGMSG("Start handle_deg2_operation_result")
    //if correct number of vertices
    if(ppgraph->order >= minVertexCount && ppgraph->order<=maxVertexCount && ppgraph->order - vertexCount <= ppgraph->degree1Count)
        handle_primpregraph_result(ppgraph);

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

    do_deg2_operations(ppgraph, &currentGenerators, currentNumberOfGenerators, multiEdgeList, multiEdgeListSize, multiEdgeOrbits, multiEdgeOrbitCount);
    DEBUGMSG("End handle_deg2_operation_result")
}

boolean isCanonicalDegree1Edge(PRIMPREGRAPH *ppgraph, int v, int *vertexOrbits){
    DEBUGASSERT(ppgraph->degree[v]==1)
    DEBUGDUMP(v,"%d")
    DEBUGARRAYDUMP(nautyLabelling, ppgraph->order, "%d")
    int reverseLabelling[ppgraph->order];
    int i;
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
            } else {
                if(reverseLabelling[i]<smallestOtherDegree1Label) smallestOtherDegree1Label = reverseLabelling[i];
            }
        }
    }
    DEBUGDUMP(smallestLabelOrbitV, "%d")
    DEBUGDUMP(smallestOtherDegree1Label, "%d")
    return (smallestLabelOrbitV < smallestOtherDegree1Label);
}

/*
 * Handles the first degree 1 operation, i.e. find all the degree 1 pairs, determine the orbits
 * apply the operation for each pair, handle the result and then revert the operation.
 */
void handle_deg1_operation1(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators){
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
            int vOrbits[ppgraph->order];
            DEBUGMSG("Start nauty")
            numberOfGenerators = 0; //reset the generators
            nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, vOrbits, &nautyOptions, &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order, canonicalGraph);
            DEBUGMSG("End nauty")
            DEBUGARRAYDUMP(vOrbits, ppgraph->order, "%d")
            
            if(isCanonicalDegree1Edge(ppgraph, deg1PairList[i][1], vOrbits)){
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
                int vOrbits[ppgraph->order];
                DEBUGMSG("Start nauty")
                numberOfGenerators = 0; //reset the generators
                nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, vOrbits, &nautyOptions, &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order, canonicalGraph);
                DEBUGMSG("End nauty")
                DEBUGARRAYDUMP(vOrbits, ppgraph->order, "%d")

                if(isCanonicalDegree1Edge(ppgraph, ppgraph->order-1, vOrbits)){
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
 * Also stores the list of multi-edges, its size, its orbits and the number of orbits. The values
 * multiEdgeList and multiEdgeOrbits will be allocated with malloc. The size and count need to be
 * allocated on the stack.
 */
boolean isCanonicalMultiEdge(PRIMPREGRAPH *ppgraph, int v1, int v2,
        VERTEXPAIR **multiEdgeList, int *multiEdgeListSize, int **multiEdgeOrbits, int *multiEdgeOrbitCount){
    DEBUGMSG("Start isCanonicalMultiEdge")
    //TODO: optimize if multiEdgeCount == 1
    if(v2<v1){
        int temp = v1;
        v1 = v2;
        v2 = temp;
    }
    int orbits[ppgraph->order];

    int maxSize = ppgraph->multiEdgeCount;
    *multiEdgeList = malloc(sizeof(VERTEXPAIR)*maxSize); //initialize an array that is large enough to hold all multi-edges
    get_multi_edges(ppgraph, *multiEdgeList, multiEdgeListSize);

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
        int n1 = (*multiEdgeList)[j][0];
        int n2 = (*multiEdgeList)[j][1];
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
        int n1 = (*multiEdgeList)[j][0];
        int n2 = (*multiEdgeList)[j][1];
        set *v1, *v2;
        v1 = GRAPHROW(ppgraph->ulgraph, n1, MAXM);
        v2 = GRAPHROW(ppgraph->ulgraph, n2, MAXM);
        ADDELEMENT(v1, n2);
        ADDELEMENT(v2, n1);
        DELELEMENT(v1, ppgraph->order + j);
        DELELEMENT(v2, ppgraph->order + j);
    }

    DEBUGMSG("End nauty")

    *multiEdgeOrbits = malloc(sizeof(int)*maxSize);
    determine_vertex_pairs_orbits(*multiEdgeList, *multiEdgeListSize, *multiEdgeOrbits, multiEdgeOrbitCount, &automorphismGroupGenerators, numberOfGenerators);

    int i = 0;
    while(i<*multiEdgeListSize && !((*multiEdgeList)[i][0]==v1 && (*multiEdgeList)[i][1]==v2)) i++;
    DEBUGASSERT(i<*multiEdgeListSize)
    int newEdgeOrbit = (*multiEdgeOrbits)[i]; //contains the number of the orbit of the edge (v1, v2)

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
        currentEdge[0]=reverseLabelling[(*multiEdgeList)[i][0]];
        currentEdge[1]=reverseLabelling[(*multiEdgeList)[i][1]];
        if(currentEdge[1]<currentEdge[0]){
            int temp = currentEdge[0];
            currentEdge[0] = currentEdge[1];
            currentEdge[1] = temp;
        }
        if((*multiEdgeOrbits)[i]==newEdgeOrbit){
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

    DEBUGMSG("End isCanonicalMultiEdge")
    return smallestRepresentantNewEdge[0] < smallestOtherMultiEdge[0] ||
            (smallestRepresentantNewEdge[0] == smallestOtherMultiEdge[0] && smallestRepresentantNewEdge[1] < smallestOtherMultiEdge[1]);
}

void handle_deg2_operation1(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators){
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

            VERTEXPAIR *multiEdgeList;
            int multiEdgeListSize;
            int *multiEdgeOrbits;
            int multiEdgeOrbitCount;

            if(isCanonicalMultiEdge(ppgraph, ppgraph->order-2, ppgraph->order-1, &multiEdgeList, &multiEdgeListSize, &multiEdgeOrbits, &multiEdgeOrbitCount))
                handle_deg2_operation_result(ppgraph, multiEdgeList, multiEdgeListSize, multiEdgeOrbits, multiEdgeOrbitCount);

            //multiEdgeList and multiEdgeOrbits are freed in do_deg2_operations

            revert_deg2_operation1(ppgraph, edgeList[i][0], edgeList[i][1]);
            DEBUGPPGRAPHPRINT(ppgraph)
        }
    }
    DEBUGMSG("End handle_deg2_operation1")
}

void handle_deg2_operation2(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators,
        VERTEXPAIR **oldMultiEdgeList, int *oldMultiEdgeListSize, int **oldMultiEdgeOrbits, int *oldMultiEdgeOrbitCount){
    DEBUGMSG("Start handle_deg2_operation2")
    DEBUG2DARRAYDUMP((*currentGenerators), currentNumberOfGenerators, ppgraph->order, "%d")
    DEBUGDUMP(*oldMultiEdgeList, "%p")
    #ifdef _DEBUG
    if(*oldMultiEdgeList!=NULL){
        DEBUGDUMP(*oldMultiEdgeListSize, "%d")
        DEBUG2DARRAYDUMP((*oldMultiEdgeList), (*oldMultiEdgeListSize), 2, "%d")
        DEBUGARRAYDUMP((*oldMultiEdgeOrbits), (*oldMultiEdgeListSize), "%d")
    }
    #endif
    if(*oldMultiEdgeList==NULL){
        //if the multi-edges have't been determined, do it now and store the result
        int maxSize = ppgraph->multiEdgeCount;
        *oldMultiEdgeList = malloc(sizeof(VERTEXPAIR)*maxSize);//initialize an array that is large enough to hold all multi-edges


        get_multi_edges(ppgraph, *oldMultiEdgeList, oldMultiEdgeListSize);

        DEBUGASSERT(*oldMultiEdgeListSize == ppgraph->multiEdgeCount)

        *oldMultiEdgeOrbits = malloc(sizeof(int)*maxSize);
        determine_vertex_pairs_orbits(*oldMultiEdgeList, *oldMultiEdgeListSize, *oldMultiEdgeOrbits, oldMultiEdgeOrbitCount, currentGenerators, currentNumberOfGenerators);
    }
    int i;
    for (i = 0; i < *oldMultiEdgeListSize; i++) {
        if((*oldMultiEdgeOrbits)[i]==i){
            DEBUGPPGRAPHPRINT(ppgraph)
            apply_deg2_operation2(ppgraph, (*oldMultiEdgeList)[i][0], (*oldMultiEdgeList)[i][1]);

            VERTEXPAIR *multiEdgeList;
            int multiEdgeListSize;
            int *multiEdgeOrbits;
            int multiEdgeOrbitCount;

            if(isCanonicalMultiEdge(ppgraph, ppgraph->order-2, ppgraph->order-1, &multiEdgeList, &multiEdgeListSize, &multiEdgeOrbits, &multiEdgeOrbitCount))
                handle_deg2_operation_result(ppgraph, multiEdgeList, multiEdgeListSize, multiEdgeOrbits, multiEdgeOrbitCount);

            //multiEdgeList and multiEdgeOrbits are freed in do_deg2_operations

            revert_deg2_operation2(ppgraph, (*oldMultiEdgeList)[i][0], (*oldMultiEdgeList)[i][1]);
            DEBUGPPGRAPHPRINT(ppgraph)
        }
    }
    DEBUGMSG("End handle_deg2_operation2")
}

void handle_deg2_operation3(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators,
        VERTEXPAIR **oldMultiEdgeList, int *oldMultiEdgeListSize, int **oldMultiEdgeOrbits, int *oldMultiEdgeOrbitCount){
    DEBUGMSG("Start handle_deg2_operation3")
    DEBUG2DARRAYDUMP((*currentGenerators), currentNumberOfGenerators, ppgraph->order, "%d")
    DEBUGDUMP(*oldMultiEdgeList, "%p")
    #ifdef _DEBUG
    if(*oldMultiEdgeList!=NULL){
        DEBUGDUMP(*oldMultiEdgeListSize, "%d")
        DEBUG2DARRAYDUMP((*oldMultiEdgeList), (*oldMultiEdgeListSize), 2, "%d")
        DEBUGARRAYDUMP((*oldMultiEdgeOrbits), (*oldMultiEdgeListSize), "%d")
    }
    #endif
    if(*oldMultiEdgeList==NULL){
        //if the multi-edges have't been determined, do it now and store the result
        int maxSize = ppgraph->multiEdgeCount;
        *oldMultiEdgeList = malloc(sizeof(VERTEXPAIR)*maxSize);//initialize an array that is large enough to hold all multi-edges

        get_multi_edges(ppgraph, *oldMultiEdgeList, oldMultiEdgeListSize);

        DEBUGASSERT(*oldMultiEdgeListSize == ppgraph->multiEdgeCount)

        *oldMultiEdgeOrbits = malloc(sizeof(int)*maxSize);
        determine_vertex_pairs_orbits(*oldMultiEdgeList, *oldMultiEdgeListSize, *oldMultiEdgeOrbits, oldMultiEdgeOrbitCount, currentGenerators, currentNumberOfGenerators);
    }
    int i;
    for (i = 0; i < *oldMultiEdgeListSize; i++) {
        if((*oldMultiEdgeOrbits)[i]==i){
            //check if the second neighbour of u and v aren't equal
            int x, y;
            x = (ppgraph->adjList[(*oldMultiEdgeList)[i][0]*3] == (*oldMultiEdgeList)[i][1]) ?
                ppgraph->adjList[(*oldMultiEdgeList)[i][0]*3+1] :
                ppgraph->adjList[(*oldMultiEdgeList)[i][0]*3];
            y = (ppgraph->adjList[(*oldMultiEdgeList)[i][1]*3] == (*oldMultiEdgeList)[i][0]) ?
                ppgraph->adjList[(*oldMultiEdgeList)[i][1]*3+1] :
                ppgraph->adjList[(*oldMultiEdgeList)[i][1]*3];
            if(x!=y){
                DEBUGPPGRAPHPRINT(ppgraph)
                apply_deg2_operation3(ppgraph, (*oldMultiEdgeList)[i][0], (*oldMultiEdgeList)[i][1]);

                VERTEXPAIR *multiEdgeList;
                int multiEdgeListSize;
                int *multiEdgeOrbits;
                int multiEdgeOrbitCount;

                if(isCanonicalMultiEdge(ppgraph, ppgraph->order-2, ppgraph->order-1, &multiEdgeList, &multiEdgeListSize, &multiEdgeOrbits, &multiEdgeOrbitCount))
                    handle_deg2_operation_result(ppgraph, multiEdgeList, multiEdgeListSize, multiEdgeOrbits, multiEdgeOrbitCount);

                //multiEdgeList and multiEdgeOrbits are freed in do_deg2_operations

                revert_deg2_operation3(ppgraph, (*oldMultiEdgeList)[i][0], (*oldMultiEdgeList)[i][1]);
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
void do_deg2_operations(PRIMPREGRAPH *ppgraph, permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators,
        VERTEXPAIR *multiEdgeList, int multiEdgeListSize, int *multiEdgeOrbits, int multiEdgeOrbitCount){
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
        handle_deg2_operation1(ppgraph, currentGenerators, currentNumberOfGenerators);
        //the orbits of the multi-edges are already determined, so we pass them on
        handle_deg2_operation2(ppgraph, currentGenerators, currentNumberOfGenerators, &multiEdgeList, &multiEdgeListSize, &multiEdgeOrbits, &multiEdgeOrbitCount);
        handle_deg2_operation3(ppgraph, currentGenerators, currentNumberOfGenerators, &multiEdgeList, &multiEdgeListSize, &multiEdgeOrbits, &multiEdgeOrbitCount);
        //free the multiEdgeList and multiEdgeOrbits before returning
        free(multiEdgeList);
        free(multiEdgeOrbits);
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
        do_deg2_operations(ppgraph, &currentGenerators, currentNumberOfGenerators, NULL, 0, NULL, 0);
    }
    DEBUGMSG("End grow")
}

void growWithoutDeg1Operations(PRIMPREGRAPH *ppgraph){
    DEBUGMSG("Start growWithoutDeg1Operations")
    int orbits[ppgraph->order];

    VERTEXPAIR *multiEdgeList = malloc(sizeof(VERTEXPAIR)*(ppgraph->multiEdgeCount)); //initialize an array that is large enough to hold all multi-edges
    int multiEdgeListSize;
    get_multi_edges(ppgraph, multiEdgeList, &multiEdgeListSize);

    DEBUGASSERT(multiEdgeListSize == ppgraph->multiEdgeCount)

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

    int *multiEdgeOrbits = malloc(sizeof(int)*(ppgraph->multiEdgeCount));
    int multiEdgeOrbitCount;
    determine_vertex_pairs_orbits(multiEdgeList, multiEdgeListSize, multiEdgeOrbits, &multiEdgeOrbitCount, &currentGenerators, currentNumberOfGenerators);

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
        do_deg2_operations(ppgraph, &currentGenerators, currentNumberOfGenerators, multiEdgeList, multiEdgeListSize, multiEdgeOrbits, multiEdgeOrbitCount);
    }
    DEBUGMSG("End growWithoutDeg1Operations")
}

static int current3RegOrder;
static PRIMPREGRAPH *currentPpgraph;

/*
 * Handles the input from the external graph by translating the graph into a pregraph primitive
 * and feeding it to the grow method.
 */
void handle_3_regular_result(graph *g){
    DEBUGMSG("Start handle_3_regular_result")
    PRIMPREGRAPH *ppgraph = currentPpgraph;
    ppgraph->order = current3RegOrder;
    ppgraph->degree1Count = 0;
    ppgraph->multiEdgeCount = 0;

    int i, l, j;

    for(i=0; i<current3RegOrder; i++){
        ppgraph->degree[i]=3;

        set *vNew, *vOld;
        vNew = GRAPHROW(ppgraph->ulgraph, i, MAXM);
        vOld = GRAPHROW(g, i, MAXM);
        EMPTYSET(vNew, MAXM);
        j=0;
        for(l=-1; (l = nextelement(vOld, MAXM, l)) >=0;){
            ppgraph->adjList[i*3+(j++)] = l;
            ADDELEMENT(vNew, l);
        }
    }

    DEBUGPPGRAPHPRINT(ppgraph)

    grow(ppgraph);
    DEBUGMSG("End handle_3_regular_result")
}

void start(){
    DEBUGMSG("Start start")
    structureCount=0;
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
        fprintf(file, "==============================\n");
        fprintf(file, "|  Graph number: %10ld  |\n", structureCount);
        fprintf(file, "|  Number of vertices:    2  |\n");
        fprintf(file, "==============================\n");
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
    fprintf(stderr, "  -L          : Allow loops.\n");
    fprintf(stderr, "  -S          : Allow semi-edges.\n");
    fprintf(stderr, "  -M          : Allow multi-edges.\n");
    fprintf(stderr, "  -f file     : Specifies the output file. If absent, the output is written to standard out.\n");
    fprintf(stderr, "  -o c        : Specifies the export format where c is one of\n");
    fprintf(stderr, "                p    planar code\n");
    fprintf(stderr, "                h    human-readable output in tabular format\n");
    fprintf(stderr, "                n    no output: only count (default)\n");
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

    while ((c = getopt(argc, argv, "LSMf:o:h")) != -1) {
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
            case 'f': //(defaults to stdout)
                outputFile = optarg;
                break;
            case 'o':
                outputType = optarg[0];
                switch (outputType) {
                    case 'n': //no output (default)
                    case 'p': //pregraph code
                    case 'h': //human-readable
                        break;
                    default:
                        usage(name);
                        return 1;
                }
                break;
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            default:
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

    /*=========== initialization ===========*/

    if(outputFile != NULL){
        FILE *file = fopen(outputFile, "r");
        if(file != NULL){
            fprintf(stderr, "File %s already exists: aborting!\n", outputFile);
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

    /*=========== generation ===========*/

    start();

    fprintf(stderr, "Found %ld structure%s with %d %s.\n", structureCount, structureCount==1 ? (char *)"" : (char *)"s",
                                                            vertexCount, vertexCount==1 ? (char *)"vertex" : (char *)"vertices");

    return EXIT_SUCCESS;
}

