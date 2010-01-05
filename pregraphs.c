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

//-----------------------------------------------------------------

/*                         o v
 *   u     v               |
 *   o     o               o u
 *  _|_   _|_            _/ \_
 * /   \_/   \   ===>   / \_/ \
 * \_________/          \_____/
 */
void apply_deg1_operation1(PRIMPREGRAPH *ppgraph, int u, int v){
    DEBUGASSERT(ppgraph->degree[u]==1 && ppgraph->degree[v]==1);
    ppgraph->degree[u]=3;
    ppgraph->adjList[u*3+1]=v;
    ppgraph->adjList[u*3+2]=ppgraph->adjList[v*3];
    ppgraph->adjList[v*3]=u;
    ppgraph->degree1Count--;

    set *gu, *gv;
    gu = GRAPHROW(ppgraph->graph, u, MAXM);
    gv = GRAPHROW(ppgraph->graph, v, MAXM);
    ADDELEMENT(gu, v);
    DELELEMENT(gv,ppgraph->adjList[v*3]);
    ADDELEMENT(gv,u);
    ADDELEMENT(gu,ppgraph->adjList[v*3]);
}

void revert_deg1_operation1(PRIMPREGRAPH *ppgraph, int u, int v){
    //the original neighbour of u is the first element in the adjList
    //so no need to change anything except the degree
    ppgraph->degree[u]=1;
    ppgraph->adjList[v*3]=ppgraph->adjList[u*3+2];
    
    set *gu, *gv;
    gu = GRAPHROW(ppgraph->graph, u, MAXM);
    gv = GRAPHROW(ppgraph->graph, v, MAXM);
    DELELEMENT(gv, u);
    ADDELEMENT(gv, ppgraph->adjList[u*3+2]);
    DELELEMENT(gu, v);
    DELELEMENT(gu, ppgraph->adjList[u*3+2]);
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
    DEBUGASSERT(areAdjacent(ppgraph, u, v));
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
    gu = GRAPHROW(ppgraph->graph, u, MAXM);
    gv = GRAPHROW(ppgraph->graph, v, MAXM);
    gs = GRAPHROW(ppgraph->graph, s, MAXM);
    gt = GRAPHROW(ppgraph->graph, t, MAXM);
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
}

void revert_deg1_operation2(PRIMPREGRAPH *ppgraph, int u, int v){
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

    set *gu, *gv, *gs, *gt;
    gu = GRAPHROW(ppgraph->graph, u, MAXM);
    gv = GRAPHROW(ppgraph->graph, v, MAXM);
    DELELEMENT(gu, s);
    DELELEMENT(gv, s);
    ADDELEMENT(gu, v);
    ADDELEMENT(gv, u);
}

//-----------------------------------------------------------------

/*
 *  __u      v__            __u  s  t  v__
 * /  o------o  \          /  o--o--o--o  \
 * |   \____/   |   ===>   |   \______/   |
 * \____________/          \______________/
 */
void apply_deg2_operation1(PRIMPREGRAPH *ppgraph, int u, int v){
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
    gu = GRAPHROW(ppgraph->graph, u, MAXM);
    gv = GRAPHROW(ppgraph->graph, v, MAXM);
    gs = GRAPHROW(ppgraph->graph, s, MAXM);
    gt = GRAPHROW(ppgraph->graph, t, MAXM);
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
}

void revert_deg2_operation1(PRIMPREGRAPH *ppgraph, int u, int v){
    int s, t, i;
    s = ppgraph->order-2;
    t = s + 1;
    i=0;
    while(ppgraph->adjList[u*3+i]!=s) i++;
    ppgraph->adjList[u*3+i]=v;
    i=0;
    while(ppgraph->adjList[v*3+i]!=t) i++;
    ppgraph->adjList[v*3+i]=u;

    ppgraph->multiEdgeCount--;
    ppgraph->order-=2;

    set *gu, *gv;
    gu = GRAPHROW(ppgraph->graph, u, MAXM);
    gv = GRAPHROW(ppgraph->graph, v, MAXM);
    DELELEMENT(gu, s);
    DELELEMENT(gv, t);
    ADDELEMENT(gu, v);
    ADDELEMENT(gv, u);
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
    DEBUGASSERT(areAdjacent(ppgraph, u, v))
    DEBUGASSERT(ppgraph->degree[u]==2 || ppgraph->multiedge[u]==v)
    DEBUGASSERT(ppgraph->degree[v]==2 || ppgraph->multiedge[v]==u)

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
    gu = GRAPHROW(ppgraph->graph, u, MAXM);
    gv = GRAPHROW(ppgraph->graph, v, MAXM);
    gs = GRAPHROW(ppgraph->graph, s, MAXM);
    gt = GRAPHROW(ppgraph->graph, t, MAXM);
    EMPTYSET(gs, MAXM);
    EMPTYSET(gt, MAXM);
    ADDELEMENT(gu, s);
    ADDELEMENT(gs, u);
    ADDELEMENT(gs, t);
    ADDELEMENT(gv, t);
    ADDELEMENT(gt, v);
    ADDELEMENT(gt, s);
}

void revert_deg2_operation2(PRIMPREGRAPH *ppgraph, int u, int v){
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
    DEBUGASSERT(areAdjacent(ppgraph, u, v))
    DEBUGASSERT(ppgraph->degree[u]==2 || ppgraph->multiedge[u]==v)
    DEBUGASSERT(ppgraph->degree[v]==2 || ppgraph->multiedge[v]==u)

    int s, t, i;
    s = ppgraph->order;
    t = s + 1;
    i=0;
    while(ppgraph->adjList[v*3+i]==u) i++; //v has degree 2, so this will stop
    ppgraph->adjList[u*3+2]=ppgraph->adjList[v*3+i];
    ppgraph->degree[u]=3;

    ppgraph->degree[v]=3;
    ppgraph->adjList[v*3]=u;
    ppgraph->adjList[v*3+1]=s;
    ppgraph->adjList[v*3+2]=t;

    ppgraph->degree[s]=2;
    ppgraph->adjList[s*3]=v;
    ppgraph->adjList[s*3+1]=t;
    ppgraph->multiedge[s]=t;

    ppgraph->degree[t]=2;
    ppgraph->adjList[t*3]=v;
    ppgraph->adjList[t*3+1]=s;
    ppgraph->multiedge[t]=s;

    ppgraph->order+=2;

    set *gu, *gv, *gs, *gt;
    gu = GRAPHROW(ppgraph->graph, u, MAXM);
    gv = GRAPHROW(ppgraph->graph, v, MAXM);
    gs = GRAPHROW(ppgraph->graph, s, MAXM);
    gt = GRAPHROW(ppgraph->graph, t, MAXM);
    EMPTYSET(gs, MAXM);
    EMPTYSET(gt, MAXM);
    int j = nextelement(gv, MAXM, -1);
    if(j==u) j = nextelement(gv, MAXM, j);
    DELELEMENT(gv, j);
    ADDELEMENT(gu, j);
    ADDELEMENT(gv, s);
    ADDELEMENT(gv, t);
    ADDELEMENT(gs, t);
    ADDELEMENT(gs, v);
    ADDELEMENT(gt, v);
    ADDELEMENT(gt, s);

}

void revert_deg2_operation3(PRIMPREGRAPH *ppgraph, int u, int v){
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
void determine_vertex_pairs_orbits(VERTEXPAIR *vertexPairList, int vertexPairListSize, int *vertexPairOrbits, int *orbitCount){
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
    for(i = 0; i < number_of_generators; i++) {
        //the generators were stored in the global variable generators by the method save_generators
        permutation = generators[i];
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
                    union_elements(vertexPairOrbits, orbitSize, orbitCount, j, k);
                    break; //the list of pairs doesn't contain any duplicates so we can stop
                }
            }
        }
    }

    DEBUGARRAYDUMP(vertexPairOrbits, vertexPairListSize, "%d")


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
    set set;
    for(i = 0; i < number_of_generators; i++) {
        //the generators were stored in the global variable generators by the method save_generators
        permutation = generators[i];
        DEBUGARRAYDUMP(permutation, currentVertexCount, "%d")

        for(j = 0; j<vertexSetListSize; j++){
            //apply permutation to current vertex pair
            EMPTYSET(set, MAXM);
            for(l=-1; (l = nextelement(vertexSetList+j, MAXM, l)) >=0;){
                ADDELEMENT(&set, permutation[l]);
            }

            //search the pair in the list
            for(k = 0; k<vertexSetListSize; k++){
                //TODO: is there a better way to check whether two sets are equal?
                l = nextelement(&set, MAXM, -1);
                while(l>=0 && ISELEMENT(vertexSetList+k, l)){
                    l = nextelement(&set, MAXM, l);
                }
                if(l==-1){
                    union_elements(vertexSetOrbits, orbitSize, orbitCount, j, k);
                    break; //the list of sets doesn't contain any duplicates so we can stop
                }
            }
        }
    }
}

void union_elements(int *forest, int *treeSizes, int *numberOfComponents, int element1, int element2){
    int root1 = find_root_of_element(forest, element1);
    int root2 = find_root_of_element(forest, element2);

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

int find_root_of_element(int *forest, int element) {
    //find with path-compression
    if(element!=forest[element]){
        forest[element]=find_root_of_element(forest, forest[element]);
    }
    return forest[element];
}

//--------------------------------------------------------------------

void handle_pregraph_result(PREGRAPH *pregraph){

}

/*
 * Handles a result in the form of a primitive pregraph. If both semi-edges and loops are
 * allowed all the different ways to select these are taken and passed on to handle_pregraph_result.
 */
void handle_primpregraph_result(PRIMPREGRAPH *ppgraph){
    DEBUGASSERT(ppgraph->order >= vertexCount)
    DEBUGASSERT(allowSemiEdges || vertexCount == ppgraph->order)
    
    int semiEdgeCount = ppgraph->order - vertexCount;
    int degree1Count = ppgraph->degree1Count;
    DEBUGASSERT(semiEdgeCount <= degree1Count)
    int loopCount = degree1Count - semiEdgeCount;

    //determine up to automorphism all the ways to select semiEdgeCount vertices
    //of degree 1 by using union-find
    int listSize, i;
    listSize = 1;
    for(i=1; i<=semiEdgeCount;i++){
        listSize = listSize*(loopCount - 1 + i)/i;
    }
    set vertexSetList[listSize];
    for(i=0; i<listSize;i++){
        EMPTYSET(vertexSetList+i,MAXM);
    }
    set tempSet;
    int position = 0;
    determine_possible_sets_of_degree1_vertices(&tempSet, vertexSetList, &position, semiEdgeCount, 0, nextDegree1Vertex(-1, ppgraph), ppgraph);
    DEBUGASSERT(position==listSize)

    int orbitCount;
    int orbits[listSize];
    determine_vertex_sets_orbits(vertexSetList, listSize, orbits, &orbitCount);
    
    //output pregraph
    PREGRAPH pregraph;
    pregraph.order = vertexCount;
    pregraph.ppgraph = ppgraph;
    for(i=0; i<listSize; i++){
        if(orbits[i]==i){
            pregraph.semiEdgeVertices = vertexSetList + i;
            handle_pregraph_result(&pregraph);
        }
    }
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
    ADDELEMENT(tempSet, currentSetElement);
    if(currentSetSize + 1 == maximumSetSize){
        //add to list
        vertexSetList[*currentListPosition] = *tempSet;
        (*currentListPosition)++;
    } else {
        determine_possible_sets_of_degree1_vertices
                (tempSet, vertexSetList, currentListPosition, maximumSetSize,
                currentSetSize + 1, nextDegree1Vertex(currentSetElement, ppgraph), ppgraph);
    }
    DELELEMENT(tempSet, currentSetElement);
    if(currentSetSize + ppgraph->degree1Count - currentSetElement > maximumSetSize){
        determine_possible_sets_of_degree1_vertices
                (tempSet, vertexSetList, currentListPosition, maximumSetSize,
                currentSetSize, nextDegree1Vertex(currentSetElement, ppgraph), ppgraph);
    }
}

/*
 * Handles the result of a degree 1 operation. If the graph has a valid size it is
 * passed on to handle_primpregraph_result, then we continue with degree 1 operations
 * and degree 2 operations.
 */
void handle_deg1_operation_result(PRIMPREGRAPH *ppgraph){
    //if correct number of vertices
    if(ppgraph->order >= minVertexCount && ppgraph->order<=maxVertexCount)
        handle_primpregraph_result(ppgraph);

    do_deg1_operations(ppgraph); //when this returns &ppgraph is unchanged

    if(allowMultiEdges){
        do_deg2_operations(ppgraph);
    }
}

/*
 * Handles the result of a degree 2 operation. If the graph has a valid size it is
 * passed on to handle_primpregraph_result, then we continue with degree 2 operations.
 * Degree 1 operations are no longer possible.
 */
void handle_deg2_operation_result(PRIMPREGRAPH *ppgraph){
    //if correct number of vertices
    if(ppgraph->order >= minVertexCount && ppgraph->order<=maxVertexCount)
        handle_primpregraph_result(ppgraph);

    do_deg2_operations(ppgraph);
}

/*
 * Handles the first degree 1 operation, i.e. find all the degree 1 pairs, determine the orbits
 * apply the operation for each pair, handle the result and then revert the operation.
 */
void handle_deg1_operation1(PRIMPREGRAPH *ppgraph){
    VERTEXPAIR deg1PairList[0];
    int listSize;
    get_deg1_pairs(ppgraph, deg1PairList, &listSize);

    int orbitCount;
    int orbits[listSize];
    determine_vertex_pairs_orbits(deg1PairList, listSize, orbits, &orbitCount);

    int i;
    for (i = 0; i < listSize; i++) {
        if(orbits[i]==i){
            apply_deg1_operation1(ppgraph, deg1PairList[i][0], deg1PairList[i][1]);
            //check if this was a valid action

            handle_deg1_operation_result(ppgraph);

            revert_deg1_operation1(ppgraph, deg1PairList[i][0], deg1PairList[i][1]);
        }
    }

}

void handle_deg1_operation2(PRIMPREGRAPH *ppgraph){

}

void handle_deg2_operation1(PRIMPREGRAPH *ppgraph){

}

void handle_deg2_operation2(PRIMPREGRAPH *ppgraph){

}

void handle_deg2_operation3(PRIMPREGRAPH *ppgraph){

}

/*
 * Performs the different degree 1 operations. When this method returns &ppgraph
 * will be unchanged.
 */
void do_deg1_operations(PRIMPREGRAPH *ppgraph){
    DEBUGASSERT(allowLoops || allowSemiEdges)
    if(ppgraph->order<=maxVertexCount) handle_deg1_operation1(ppgraph);
    if(ppgraph->order<=maxVertexCount-2) handle_deg1_operation2(ppgraph);
}

/*
 * Performs the different degree 2 operations. When this method returns &ppgraph
 * will be unchanged.
 */
void do_deg2_operations(PRIMPREGRAPH *ppgraph){
    DEBUGASSERT(allowMultiEdges)
    if(ppgraph->order<=maxVertexCount-2){
        handle_deg2_operation1(ppgraph);
        handle_deg2_operation2(ppgraph);
        handle_deg2_operation3(ppgraph);
    }
}

void grow(PRIMPREGRAPH *ppgraph){
    if(allowLoops || allowSemiEdges){
        do_deg1_operations(ppgraph);
    }
    if(allowMultiEdges){
        do_deg2_operations(ppgraph);
    }

}

void start(){
    if(!allowSemiEdges){
        minVertexCount = maxVertexCount = vertexCount;
    } else {
        minVertexCount = vertexCount;
        maxVertexCount = 2*vertexCount; //TODO: better and correct(?) upperbound
    }
    PRIMPREGRAPH ppgraph;
    if(allowLoops || allowSemiEdges){
        construct_K2(&ppgraph);
        grow(&ppgraph);
    }
    if(allowMultiEdges){
        construct_C4(&ppgraph);
        grow(&ppgraph);
    }
    if((allowLoops || allowSemiEdges) && allowMultiEdges){
        construct_K3_with_spike(&ppgraph);
        grow(&ppgraph);
    }
    //TODO: start generation of cubic graphs
}

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
    g0 = GRAPHROW(ppgraph->graph, 0, MAXM);
    g1 = GRAPHROW(ppgraph->graph, 1, MAXM);
    EMPTYSET(g0, MAXM);
    EMPTYSET(g1, MAXM);
    ADDELEMENT(g0, 1);
    ADDELEMENT(g1, 0);
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
    g0 = GRAPHROW(ppgraph->graph, 0, MAXM);
    g1 = GRAPHROW(ppgraph->graph, 1, MAXM);
    g2 = GRAPHROW(ppgraph->graph, 2, MAXM);
    g3 = GRAPHROW(ppgraph->graph, 3, MAXM);
    EMPTYSET(g0, MAXM);
    EMPTYSET(g1, MAXM);
    EMPTYSET(g2, MAXM);
    EMPTYSET(g3, MAXM);
    ADDELEMENT(g0, 3); ADDELEMENT(g0, 1);
    ADDELEMENT(g1, 0); ADDELEMENT(g1, 2);
    ADDELEMENT(g2, 1); ADDELEMENT(g2, 3);
    ADDELEMENT(g3, 2); ADDELEMENT(g3, 0);
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
    g0 = GRAPHROW(ppgraph->graph, 0, MAXM);
    g1 = GRAPHROW(ppgraph->graph, 1, MAXM);
    g2 = GRAPHROW(ppgraph->graph, 2, MAXM);
    g3 = GRAPHROW(ppgraph->graph, 3, MAXM);
    EMPTYSET(g0, MAXM);
    EMPTYSET(g1, MAXM);
    EMPTYSET(g2, MAXM);
    EMPTYSET(g3, MAXM);
    ADDELEMENT(g0, 1); ADDELEMENT(g0, 2); ADDELEMENT(g0, 3);
    ADDELEMENT(g1, 0); ADDELEMENT(g1, 2);
    ADDELEMENT(g2, 1); ADDELEMENT(g2, 0);
    ADDELEMENT(g3, 0);
}

//-------------------------End start graphs------------------------------

//----------------------Begin Nauty interaction--------------------------

void init_nauty_options() {
    //TODO also options without getcanon?
    options.getcanon = TRUE;
    options.userautomproc = save_generators;
    options.writeautoms = TRUE;
    //options.writemarkers = TRUE;
}

void save_generators(int count, permutation perm[], nvector orbits[],
        int numorbits, int stabvertex, int n) {
    memcpy(generators + number_of_generators, perm, sizeof(permutation) * n);

    number_of_generators++;
}

//------------------------End Nauty interaction-------------------------

#ifdef NO_MAIN
    #define MAIN_FUNCTION nomain
#else
    #define MAIN_FUNCTION main
#endif
/*
 * 
 */
int MAIN_FUNCTION(int argc, char** argv) {

    return (EXIT_SUCCESS);
}

