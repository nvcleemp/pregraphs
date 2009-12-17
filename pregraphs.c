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
    int i = 0;
    while(i<ppgraph->order && !ISELEMENT(gv,i)) i++;
    if(i<ppgraph->order){
        DELELEMENT(gv,i);
        ADDELEMENT(gv,u);
        ADDELEMENT(gu,i);
    } else {
        DEBUGMSG("Something went horribly wrong!")
        exit(EXIT_FAILURE);
    }
}

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

void handle_primpregraph_result(PRIMPREGRAPH *ppgraph){
    DEBUGASSERT(ppgraph->order >= vertexCount)
    DEBUGASSERT(allowSemiEdges || vertexCount == ppgraph->order)
    
    int semiEdgeCount = ppgraph->order - vertexCount;
    DEBUGASSERT(semiEdgeCount <= ppgraph->degree1Count)
    int loopCount = ppgraph->degree1Count - semiEdgeCount;
    //determine up to automorphism all the ways to select semiEdgeCount vertices
    //of degree 1 by using union-find

    //output pregraph
}

void handle_deg1_operation_result(PRIMPREGRAPH *ppgraph){
    //if correct number of vertices
    handle_primpregraph_result(ppgraph);

    //if still more vertices possible
    do_deg1_operations(ppgraph);

    if(allowMultiEdges){
        do_deg2_operations(ppgraph);
    }
}

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
            handle_deg1_operation1(ppgraph);
        }
    }

}

void do_deg1_operations(PRIMPREGRAPH *ppgraph){

}

void do_deg2_operations(PRIMPREGRAPH *ppgraph){

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
        //
        grow(&ppgraph);
    }
    //TODO: start generation of cubic graphs
}

//------------------------Begin start graphs-----------------------------

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

