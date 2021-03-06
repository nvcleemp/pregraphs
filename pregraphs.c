/* pregraphs.c
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
 * File:   pregraphs.c
 * Author: nvcleemp
 *
 * Created on December 8, 2009, 2:32 PM
 */

//#define _TEST
//#define _DEBUG
//#define _CONSISTCHECK
//#define _PROFILING

#include "pregraphs.h"

#include <execinfo.h>
#include <sys/times.h>

void determine_vertex_pairs_orbits(VERTEXPAIR *vertexPairList, int vertexPairListSize,
        int *vertexPairOrbits, int *orbitCount, permutation (*currentGenerators)[MAXN][MAXN] ,
        int currentNumberOfGenerators);

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

//------------------Start colouring methods--------------------------------
void initIsColourable(PRIMPREGRAPH *ppgraph) {
    RESETMARKS_EDGES;
    int i, j;
    for (i = 0; i < ppgraph->order; i++) {
        coloursAroundVertex[i] = 0;
        for(j = 0; j < ppgraph->degree[i]; j++) {
            neighbourToIndexMapping[i][ppgraph->adjList[3*i+j]] = j;
        }
    }

}

void copyColourInformation(int order) {
    int i, j;
    for (i = 0; i < order; i++) {
        coloursAroundVertexCopy[i] = coloursAroundVertex[i];
        for(j = 0; j < order; j++) {
            neighbourToIndexMappingCopy[i][j] = neighbourToIndexMapping[i][j];
        }
        for(j = 0; j < 3; j++) {
            coloursCopy[i][j] = colours[i][j];
        }
    }
}

void restoreColourInformation(int order) {
    int i, j;
    for (i = 0; i < order; i++) {
        coloursAroundVertex[i] = coloursAroundVertexCopy[i];
        for(j = 0; j < order; j++) {
            neighbourToIndexMapping[i][j] = neighbourToIndexMappingCopy[i][j];
        }
        for(j = 0; j < 3; j++) {
            colours[i][j] = coloursCopy[i][j];
        }
    }
}

void printColouring(FILE *f, PRIMPREGRAPH *ppgraph){
    int i, j;
    for(i = 0; i < ppgraph->order; i++){
        fprintf(f, "%2d) ", i+1);
        for(j = 0; j < ppgraph->degree[i]; j++){
            if(ISMARKED_EDGES(i, j)){
                fprintf(f, "%d ", colours[i][j]);
            } else {
                fprintf(f, "0 ");
            }
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
}

void printMapping(FILE *f, PRIMPREGRAPH *ppgraph){
    int i,j;
    fprintf(f, "    ");
    for(i = 0; i < ppgraph->order; i++){
        fprintf(f, "%2d ", i+1);
    }
    fprintf(f, "\n");
    for(i = 0; i < ppgraph->order; i++){
        fprintf(f, "%2d) ", i+1);
        for(j = 0; j < ppgraph->order; j++){
            fprintf(f, "%2d ", neighbourToIndexMapping[i][j]);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
}

int isConflictingColouring(int currentVertex, int colour) {
    int i;
    for(i = 0; i < 3; i++) {
        if(ISMARKED_EDGES(currentVertex, i) && colours[currentVertex][i] == colour)
            return 1;
    }
    return 0;

}


inline void determineAvailableColours(int usedColour, int *availableColours) {
    switch(usedColour) {
        case 1:
            availableColours[0] = 2;
            availableColours[1] = 3;
            break;
        case 2:
            availableColours[0] = 1;
            availableColours[1] = 3;
            break;
        case 3:
            availableColours[0] = 1;
            availableColours[1] = 2;
            break;
        default:
            fprintf(stdout, "Invalid previous colour");
            exit(0);
    }
}

inline int determineMissingColour(int sumColours){
    switch(sumColours) {
        case 3:
            return 3;
            break;
        case 4:
            return 2;
            break;
        case 5:
            return 1;
            break;
        default:
            ERRORMSG("Invalid sum of colours");
    }
}

void determineUncolouredVertex(int vertex, int *uncolouredVertex, int *missingColour, PRIMPREGRAPH *ppgraph) {
    DEBUGASSERT(coloursAroundVertex[vertex] == 2);

    int i;
    int sumColours = 0;
    for(i = 0; i < ppgraph->degree[vertex]; i++) {
        if(!ISMARKED_EDGES(vertex, i)) {
            *uncolouredVertex = ppgraph->adjList[3*vertex+i];
        } else {
            sumColours += colours[vertex][i];
        }
    }

    *missingColour = determineMissingColour(sumColours);
    
}

void unmarkColours(int nonfree_labelled[][2], int nonfree_labelled_size) {
    int i;
    int vertex0, vertex1;
    for(i = 0; i < nonfree_labelled_size; i++) {
        vertex0 = nonfree_labelled[i][0];
        vertex1 = nonfree_labelled[i][1];
        UNMARK_EDGES(vertex0, neighbourToIndexMapping[vertex0][vertex1]);
        UNMARK_EDGES(vertex1, neighbourToIndexMapping[vertex1][vertex0]);
        coloursAroundVertex[vertex0]--;
        coloursAroundVertex[vertex1]--;
    }
}

/*
 * Starts from currentvertex and sets colours that are fixed by the current colours.
 * The edges that are coloured are stored, so this can be rolled back in case of a conflict.
 */
boolean propagateFixedColours(int currentVertex, int nonfree_labelled[][2], int *nonfree_labelled_size, PRIMPREGRAPH *ppgraph) {
    int uncolouredVertex = -1, missingColour = -1; //values assigned to avoid compiler warnings
    //while the colour is fixed for currentVertex
    while(coloursAroundVertex[currentVertex] == 2) {
        //find the colour of the remaining edge
        determineUncolouredVertex(currentVertex, &uncolouredVertex, &missingColour, ppgraph);
        //check that this colour gives no conflicts
        if(!isConflictingColouring(uncolouredVertex, missingColour)) {
            int indexUncolouredVertex = neighbourToIndexMapping[currentVertex][uncolouredVertex];
            int indexCurrentVertex = neighbourToIndexMapping[uncolouredVertex][currentVertex];
            colours[currentVertex][indexUncolouredVertex] = missingColour;
            colours[uncolouredVertex][indexCurrentVertex] = missingColour;

            //DEBUGASSERT(missing_colour > 0 && missing_colour < 4);

            MARK_EDGES(currentVertex, indexUncolouredVertex);
            MARK_EDGES(uncolouredVertex, indexCurrentVertex);
            coloursAroundVertex[currentVertex] = 3;
            coloursAroundVertex[uncolouredVertex]++;

            nonfree_labelled[*nonfree_labelled_size][0] = currentVertex;
            nonfree_labelled[*nonfree_labelled_size][1] = uncolouredVertex;
            (*nonfree_labelled_size)++;

            currentVertex = uncolouredVertex;
        } else {
            //in case of conflicts: remove colours and return FALSE
            unmarkColours(nonfree_labelled, *nonfree_labelled_size);
            return FALSE;
        }
    }
    return TRUE;

}

/*
 * For this method we assume that the graph is a simple graph.
 */
int tryExtendingColouring(int numberOfColouredEdges, int numberOfEdges, PRIMPREGRAPH *ppgraph) {
    if(numberOfColouredEdges != numberOfEdges) {
        int currentVertex;

        for(currentVertex = 0; currentVertex < ppgraph->order; currentVertex++) {
            if(coloursAroundVertex[currentVertex] == 1 && ppgraph->degree[currentVertex] != 1) {
                break;
            }
            DEBUGASSERT(coloursAroundVertex[currentVertex] != 2)
        }
        DEBUGASSERT(currentVertex < ppgraph->order);

        int usedColour = -1; //the colour already used at this vertex
                             //value assigned to avoid compiler warnings
        int i;
        for(i = 0; i < ppgraph->degree[currentVertex]; i++) {
            if(ISMARKED_EDGES(currentVertex, i)) {
                usedColour = colours[currentVertex][i];
                break;
            }
        }
        DEBUGASSERT(i < ppgraph->degree[currentVertex]);

        int availableVertices[2];
        int indexAvailableVertex0 = (i + 1) % 3;
        int indexAvailableVertex1 = (i + 2) % 3;
        availableVertices[0] = ppgraph->adjList[3*currentVertex+indexAvailableVertex0];
        availableVertices[1] = ppgraph->adjList[3*currentVertex+indexAvailableVertex1];

        int availableColours[2];
        determineAvailableColours(usedColour, availableColours);

        DEBUGASSERT((availableColours[0] > 0 && availableColours[0] < 4) && (availableColours[1] > 0 && availableColours[1] < 4));

        int nonfree_labelled[numberOfEdges - numberOfColouredEdges][2];
        int j;
        for(i = 0; i < 2; i++) {
            if(isConflictingColouring(availableVertices[0], availableColours[i]) ||
                    isConflictingColouring(availableVertices[1], availableColours[(i + 1) % 2])) {
                continue;
            }

            int indexCurrentVertex0 = neighbourToIndexMapping[availableVertices[0]][currentVertex];
            int indexCurrentVertex1 = neighbourToIndexMapping[availableVertices[1]][currentVertex];

            colours[availableVertices[0]][indexCurrentVertex0] = availableColours[i];
            colours[availableVertices[1]][indexCurrentVertex1] = availableColours[(i + 1) % 2];

            MARK_EDGES(availableVertices[0], indexCurrentVertex0);
            MARK_EDGES(availableVertices[1], indexCurrentVertex1);

            coloursAroundVertex[availableVertices[0]]++;
            coloursAroundVertex[availableVertices[1]]++;

            int nonfree_labelled_size = 0;
            boolean abort = FALSE;
            for(j = 0; j < 2;j++) {
                if(!propagateFixedColours(availableVertices[j], nonfree_labelled, &nonfree_labelled_size, ppgraph)) {
                    DEBUGASSERT(nonfree_labelled_size <= numberOfEdges - numberOfColouredEdges)
                    abort = TRUE;
                    break;
                }
                DEBUGASSERT(nonfree_labelled_size <= numberOfEdges - numberOfColouredEdges)
            }

            if(!abort) {
                colours[currentVertex][indexAvailableVertex0] = availableColours[i];
                colours[currentVertex][indexAvailableVertex1] = availableColours[(i + 1) % 2];
                MARK_EDGES(currentVertex, indexAvailableVertex0);
                MARK_EDGES(currentVertex, indexAvailableVertex1);
                coloursAroundVertex[currentVertex] = 3;

                if(tryExtendingColouring(numberOfColouredEdges + nonfree_labelled_size + 2, numberOfEdges, ppgraph)) {
                    return 1;
                } else {
                    unmarkColours(nonfree_labelled, nonfree_labelled_size);
                }
                UNMARK_EDGES(currentVertex, indexAvailableVertex0);
                UNMARK_EDGES(currentVertex, indexAvailableVertex1);
                coloursAroundVertex[currentVertex] = 1;
            }

            UNMARK_EDGES(availableVertices[0], indexCurrentVertex0);
            UNMARK_EDGES(availableVertices[1], indexCurrentVertex1);

            //number_of_colours_snarks[current_vertex] = 1;
            coloursAroundVertex[availableVertices[0]]--;
            coloursAroundVertex[availableVertices[1]]--;
        }
        return 0;
    } else {
        return 1;
    }
}

/*
 * This method only works for simple graphs!
 * Stores the colouring in the array
 */
boolean isColourableGraph(PRIMPREGRAPH *ppgraph) {
    if(ppgraph->bridgeCount > 0 && ppgraph->degree1Count==0)
        return FALSE;
    if(ppgraph->degree1Count == 1)
        return FALSE;

    initIsColourable(ppgraph);

    int currentVertex = 0;
    int i, neighbour, currentIndex;
    for(i = 0; i < ppgraph->degree[currentVertex]; i++) {
        colours[currentVertex][i] = i + 1;
        neighbour = ppgraph->adjList[3*currentVertex+i];
        currentIndex = neighbourToIndexMapping[neighbour][currentVertex];
        colours[neighbour][currentIndex] = i + 1;

        MARK_EDGES(currentVertex, i);
        MARK_EDGES(neighbour, currentIndex);
        coloursAroundVertex[neighbour] = 1;
    }
    coloursAroundVertex[currentVertex] = ppgraph->degree[currentVertex];
    int numberOfEdges = (3*(ppgraph->order-ppgraph->degree1Count) + ppgraph->degree1Count)/2;
    return tryExtendingColouring(ppgraph->degree[currentVertex], numberOfEdges, ppgraph);
}

/*
 * This method only works for simple graphs!
 * Does BFS to check whether graph is bipartite
 */
boolean isBipartiteGraph(PRIMPREGRAPH *ppgraph) {
    int i, queueHead, queueTail;
    int queue[MAXN];
    boolean queued[MAXN];
    for(i=0; i<ppgraph->order; i++) {
        queued[i] = FALSE;
        vertexColours[i]=-1;
    }

    queueHead = 0;
    queueTail = 1;
    queue[0] = 0;
    queued[0] = TRUE;

    while(queueTail > queueHead){
        int currentVertex = queue[queueHead++];
        int numberOfColouredNeighbours = 0, sumOfColours = 0;

        //put all neighbours in queue
        for(i=0; i<ppgraph->degree[currentVertex]; i++){
            if(!queued[ppgraph->adjList[3*currentVertex+i]]){
                queued[ppgraph->adjList[3*currentVertex+i]]=TRUE;
                queue[queueTail++] = ppgraph->adjList[3*currentVertex+i];
            } else if (vertexColours[ppgraph->adjList[3*currentVertex+i]]>-1){
                numberOfColouredNeighbours++;
                sumOfColours+=vertexColours[ppgraph->adjList[3*currentVertex+i]];
            }
        }

        if(sumOfColours==0){
            //either all coloured neighbours have colour 0 or no coloured neighbours
            vertexColours[currentVertex] = 1;
        } else if(sumOfColours==numberOfColouredNeighbours){
            //all coloured neighbours have colour 1
            vertexColours[currentVertex] = 0;
        } else {
            //conflicting colours
            return FALSE;
        }
    }

    return TRUE;
}

int getNeighbourAtColour(PRIMPREGRAPH *ppgraph, int vertex, int colour){
    int i;
    for(i = 0; i < ppgraph->degree[vertex]; i++){
        if(colours[vertex][i]==colour){
            return ppgraph->adjList[3*vertex+i];
        }
    }
    return -1;
}

inline void setColour(PRIMPREGRAPH *ppgraph, int v1, int v2, int colour){
    int i=0;
    while(i<ppgraph->degree[v1] && ppgraph->adjList[3*v1+i]!=v2) i++;

    if(i<ppgraph->degree[v1]){
        colours[v1][i]=colour;
    } else {
        ERRORMSG("v2 is not a neighbour of v1")
    }
    i=0;
    while(i<ppgraph->degree[v2] && ppgraph->adjList[3*v2+i]!=v1) i++;

    if(i<ppgraph->degree[v2]){
        colours[v2][i]=colour;
    } else {
        ERRORMSG("v1 is not a neighbour of v2")
    }
}

void switchColours(PRIMPREGRAPH *ppgraph, int vertex, int firstColour, int secondColour){
    int previousVertex;
    int currentVertex = vertex;
    int nextVertex = getNeighbourAtColour(ppgraph, currentVertex, firstColour);
    int currentColour = firstColour;

    while(nextVertex>=0){
        previousVertex = currentVertex;
        currentVertex = nextVertex;
        currentColour = currentColour == firstColour ? secondColour : firstColour;
        nextVertex = getNeighbourAtColour(ppgraph, currentVertex, currentColour);
        setColour(ppgraph, previousVertex, currentVertex, currentColour);
    }
}

/*
 * This method only works for simple graphs!
 * Tries to find a colouring such that the two edges incident with the two edges of degree 1
 * have a different colouring
 */
boolean canBeColouredDifferently(PRIMPREGRAPH *ppgraph, int v1, int v2) {
    if(ppgraph->degree1Count == 2)
        return FALSE;
    if(ppgraph->degree[v1]!=1 || ppgraph->degree[v2]!=1){
        ERRORMSG("Not two vertices of degree 1")
    }

    initIsColourable(ppgraph);

    colours[v1][0]=1;
    colours[v2][0]=2;

    int v1Neighbour = ppgraph->adjList[3*v1+0];
    int v2Neighbour = ppgraph->adjList[3*v2+0];
    int v1Index = neighbourToIndexMapping[v1Neighbour][v1];
    int v2Index = neighbourToIndexMapping[v2Neighbour][v2];

    colours[v1Neighbour][v1Index] = 1;
    colours[v2Neighbour][v2Index] = 2;

    coloursAroundVertex[v1] = 1;
    coloursAroundVertex[v2] = 1;
    coloursAroundVertex[v1Neighbour]++;
    coloursAroundVertex[v2Neighbour]++;
    MARK_EDGES(v1, 0);
    MARK_EDGES(v2, 0);
    MARK_EDGES(v1Neighbour, v1Index);
    MARK_EDGES(v2Neighbour, v2Index);

    int numberOfEdges = (3*(ppgraph->order-ppgraph->degree1Count) + ppgraph->degree1Count)/2;
    return tryExtendingColouring(2, numberOfEdges, ppgraph);
}

/*
 * Tries to switch the colours colour1 and colour2 along a Kempe path starting in start without
 * changing the colour of the edge incident with illegalEnd. Returns TRUE when succesfull
 * and FALSE when not. Only in the case of TRUE the colouring will have changed.
 */
boolean trySwitchingColours(PRIMPREGRAPH *ppgraph, int start, int illegalEnd, int colour1, int colour2){
    int current = start;
    int next = ppgraph->adjList[3*current+0];
    int currentColour = colour1;
    int nextColour = colour2;
    while(ppgraph->degree[next]!=1){
	    int temp = currentColour;
        currentColour = nextColour;
        nextColour = temp;
	    current = next;
		next = getNeighbourAtColour(ppgraph, current, currentColour);
    }

    if(next==illegalEnd){
        return FALSE;
    } else {
        //switch the colours
        switchColours(ppgraph, start, colour1, colour2);
        return TRUE;
    }
}

/*
 * This method only works for simple graphs!
 * Tries to find a colouring such that the two edges incident with the two edges of degree 1
 * have a different colouring. First tries to switch the colours along the two Kempe paths.
 * If this is unsuccesful just calls canBeColouredDifferently, but first stores
 * the current colouring. If the answer is false, it restores the old colouring.
 */
boolean tryAlternateColouring(PRIMPREGRAPH *ppgraph, int v1, int v2){
    //find path starting from v1
    if(trySwitchingColours(ppgraph, v1, v2, colours[v1][0], colours[v1][0]%3 + 1)){
        #ifdef _PROFILING_OPERATION1_SAME_COLOUR
        degree1Operation1SolvedWithPaths++;
        #endif
        return TRUE;
    }
    if(trySwitchingColours(ppgraph, v1, v2, colours[v1][0], (colours[v1][0]%3 + 1)%3 + 1)){
        #ifdef _PROFILING_OPERATION1_SAME_COLOUR
        degree1Operation1SolvedWithPaths++;
        #endif
        return TRUE;
    }

    //copy the colours
    int i, j;
    for(i = 0; i < ppgraph->order; i++){
        for(j = 0; j < ppgraph->degree[i]; j++){
            coloursCopy[i][j]=colours[i][j];
        }
        coloursAroundVertexCopy[i]=coloursAroundVertex[i];
        //TODO: do we need to guarantee that coloursAroundVertex
        //stays up to date? It is probably only used in colouring methods
    }

    boolean alternateColouringExists = canBeColouredDifferently(ppgraph, v1, v2);

    if(!alternateColouringExists){
        for(i = 0; i < ppgraph->order; i++){
            for(j = 0; j < ppgraph->degree[i]; j++){
                colours[i][j]=coloursCopy[i][j];
            }
            coloursAroundVertex[i]=coloursAroundVertexCopy[i];
        }
    } else {
        #ifdef _PROFILING_OPERATION1_SAME_COLOUR
        degree1Operation1NotSolvedWithPaths++;
        #endif
    }

    return alternateColouringExists;
}

inline int getColour(PRIMPREGRAPH *ppgraph, int v1, int v2){
    int i=0;
    while(i<ppgraph->degree[v1] && ppgraph->adjList[3*v1+i]!=v2) i++;

    if(i<ppgraph->degree[v1]){
        return colours[v1][i];
    } else {
        ERRORMSG("v2 is not a neighbour of v1")
    }
}

/*
 * Performs BFS to switch vertex colours in the subgraph corresponding to the
 * component containing vertex after vertexNot has been removed
 */
void switchVertexColours(PRIMPREGRAPH *ppgraph, int vertex, int vertexNot){
    int i, queueHead, queueTail;
    int queue[MAXN];
    boolean queued[MAXN];
    for(i=0; i<ppgraph->order; i++) {
        queued[i] = FALSE;
    }

    queueHead = 0;
    queueTail = 0;
    queued[vertex] = TRUE;
    queued[vertexNot] = TRUE;

    //switch colour of vertex
    vertexColours[vertex] = 1 - vertexColours[vertex];
    
    //add all neighbours of vertex except vertexNot to the queue
    for(i=0; i<ppgraph->degree[vertex]; i++){
        if(ppgraph->adjList[3*vertex+i]!=vertexNot){
            queued[ppgraph->adjList[3*vertex+i]]=TRUE;
            queue[queueTail++] = ppgraph->adjList[3*vertex+i];
        }
    }

    while(queueTail > queueHead){
        int currentVertex = queue[queueHead++];

        //switch colour of vertex
        vertexColours[currentVertex] = 1 - vertexColours[currentVertex];

        //put all neighbours in queue
        for(i=0; i<ppgraph->degree[currentVertex]; i++){
            if(!queued[ppgraph->adjList[3*currentVertex+i]]){
                queued[ppgraph->adjList[3*currentVertex+i]]=TRUE;
                queue[queueTail++] = ppgraph->adjList[3*currentVertex+i];
            }
        }


    }
}

//-------------------End colouring methods---------------------------------

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

inline int findBridge(PRIMPREGRAPH *ppgraph, int u, int v){
    DEBUGMSG("Start findBridge")
    DEBUGDUMP(ppgraph->bridgeCount, "%d")
    DEBUGDUMP(u, "%d")
    DEBUGDUMP(v, "%d")
    /*
    int min = u < v ? u : v;
    int max = u < v ? v : u;
    int i = 0;
    while(i < ppgraph->bridgeCount &&
            (ppgraph->bridges[i][0]!=max || ppgraph->bridges[i][1]!=min)){
        i++;
    }
    DEBUGASSERT(i < ppgraph->bridgeCount)
    if(i >= ppgraph->bridgeCount){
        fprintf(stderr, "ERROR\n");
    }*/
    int i = ppgraph->bridgePosition[u][v];
    DEBUGDUMP(i, "%d")
    DEBUGMSG("End findBridge")
    return i;
}

inline void removeBridge(PRIMPREGRAPH *ppgraph, int index){
    DEBUGMSG("Start removeBridge")
    DEBUGDUMP(index, "%d")
    DEBUGDUMP(ppgraph->bridgeCount, "%d")
    ppgraph->bridgeCount--;
    if(ppgraph->bridgeCount > 0 && index < ppgraph->bridgeCount) {
        ppgraph->bridges[index][0] = ppgraph->bridges[ppgraph->bridgeCount][0];
        ppgraph->bridges[index][1] = ppgraph->bridges[ppgraph->bridgeCount][1];
        ppgraph->bridgePosition[ppgraph->bridges[index][0]][ppgraph->bridges[index][1]] = index;
        ppgraph->bridgePosition[ppgraph->bridges[index][1]][ppgraph->bridges[index][0]] = index;
    }
    DEBUGDUMP(ppgraph->bridgeCount, "%d")
    DEBUGMSG("End removeBridge")
}

inline void removeBridgesAtEnd(PRIMPREGRAPH *ppgraph, int count){
    DEBUGMSG("Start removeBridgesAtEnd")
    DEBUGDUMP(count, "%d")
    DEBUGDUMP(ppgraph->bridgeCount, "%d")
    ppgraph->bridgeCount-=count;
    DEBUGDUMP(ppgraph->bridgeCount, "%d")
    DEBUGMSG("End removeBridgesAtEnd")
}

//insert the bridge at position index, and move the current bridge at that
//position to the end
inline void insertBridge(PRIMPREGRAPH *ppgraph, int index, int v1, int v2){
    DEBUGMSG("Start insertBridge")
    DEBUGDUMP(v1, "%d")
    DEBUGDUMP(v2, "%d")
    DEBUGDUMP(index, "%d")
    DEBUGDUMP(ppgraph->bridgeCount, "%d")
    int min = v1 < v2 ? v1 : v2;
    int max = v1 < v2 ? v2 : v1;
    if(index != ppgraph->bridgeCount) {
        ppgraph->bridges[ppgraph->bridgeCount][0] = ppgraph->bridges[index][0];
        ppgraph->bridges[ppgraph->bridgeCount][1] = ppgraph->bridges[index][1];
        ppgraph->bridgePosition[ppgraph->bridges[ppgraph->bridgeCount][0]][ppgraph->bridges[ppgraph->bridgeCount][1]] =
                ppgraph->bridgeCount;
        ppgraph->bridgePosition[ppgraph->bridges[ppgraph->bridgeCount][1]][ppgraph->bridges[ppgraph->bridgeCount][0]] =
                ppgraph->bridgeCount;
    }
    ppgraph->bridges[index][0] = max;
    ppgraph->bridges[index][1] = min;
    ppgraph->bridgePosition[min][max] = index;
    ppgraph->bridgePosition[max][min] = index;
    ppgraph->bridgeCount++;
    DEBUGDUMP(ppgraph->bridgeCount, "%d")
    DEBUGMSG("End insertBridge")
}

inline void appendBridge(PRIMPREGRAPH *ppgraph, int v1, int v2){
    DEBUGMSG("Start appendBridge")
    DEBUGDUMP(v1, "%d")
    DEBUGDUMP(v2, "%d")
    DEBUGDUMP(ppgraph->bridgeCount, "%d")
    int min = v1 < v2 ? v1 : v2;
    int max = v1 < v2 ? v2 : v1;
    ppgraph->bridges[ppgraph->bridgeCount][0] = max;
    ppgraph->bridges[ppgraph->bridgeCount][1] = min;
    ppgraph->bridgePosition[min][max] = ppgraph->bridgeCount;
    ppgraph->bridgePosition[max][min] = ppgraph->bridgeCount;
    ppgraph->bridgeCount++;
    DEBUGDUMP(ppgraph->bridgeCount, "%d")
    DEBUGMSG("End appendBridge")
}

int bridgesDFSNumber[MAXN];
int bridgesDFSReach[MAXN];
int bridgesDFSNextLabel;

void bridgesDFS(PRIMPREGRAPH *ppgraph, int v, int parent){
    int i, b;
    bridgesDFSNumber[v]=bridgesDFSReach[v]=bridgesDFSNextLabel;
    bridgesDFSNextLabel++;
    for(i = 0; i < ppgraph->degree[v]; i++){
        b = ppgraph->adjList[v*3 + i];
        if(bridgesDFSNumber[b] == MAXN){
            bridgesDFS(ppgraph, b, v);

            if(bridgesDFSReach[b] < bridgesDFSReach[v]){
                bridgesDFSReach[v] = bridgesDFSReach[b];
            }
            if(bridgesDFSReach[b] > bridgesDFSNumber[v]){
                ppgraph->bridges[ppgraph->bridgeCount][0] = MAX(v, b);
                ppgraph->bridges[ppgraph->bridgeCount][1] = MIN(v, b);
                ppgraph->bridgePosition[b][v] = ppgraph->bridgeCount;
                ppgraph->bridgePosition[v][b] = ppgraph->bridgeCount;
                ppgraph->bridgeCount++;
            }
        } else if(b != parent && bridgesDFSNumber[b] < bridgesDFSReach[v]){
            bridgesDFSReach[v] = bridgesDFSNumber[b];
        }
    }
}

inline void determineBridges(PRIMPREGRAPH *ppgraph){
    DEBUGMSG("Start determineBridges")
    ppgraph->bridgeCount = 0;

    bridgesDFSNextLabel = 0;
    int i;
    for(i = 0; i < MAXN; i++) bridgesDFSNumber[i] = MAXN;

    bridgesDFS(ppgraph, 0, 0);

    DEBUGMSG("End determineBridges")
}

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
    DEBUGASSERT(!onlyColourable || colours[u][0] != colours[v][0])
    DEBUGASSERT(!onlyBipartite || vertexColours[u][0] == vertexColours[v][0])
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

    //adjust colouring
    if(onlyColourable){
        int cu = colours[u][0];
        int cv = colours[v][0];
        colours[u][1] = determineMissingColour(cu+cv);
        colours[u][2] = cv;
        colours[v][0] = determineMissingColour(cu+cv);
        coloursAroundVertex[u]=3;
    }

    //adjust vertex colouring
    if(onlyBipartite){
        vertexColours[v] = 1 - vertexColours[u];
    }

    //keep track of bridges
    int oldBridge1 = findBridge(ppgraph, u, ppgraph->adjList[u*3]);

    changedBridges[degree1OperationsDepth][0] = oldBridge1;
    removeBridge(ppgraph, oldBridge1);

    int oldBridge2 = findBridge(ppgraph, v, t);
    changedBridges[degree1OperationsDepth][1] = oldBridge2;
    numberOfChangedBridges[degree1OperationsDepth] = 2;
    removeBridge(ppgraph, oldBridge2);
    appendBridge(ppgraph, u, v);

    DEBUGPPGRAPHPRINT(ppgraph)
    DEBUGBRIDGESPRINT(ppgraph)

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

    //restore colouring
    if(onlyColourable){
        colours[v][0]=colours[u][2];
        coloursAroundVertex[u]=1;
    }

    //restore vertex colouring
    if(onlyBipartite){
        vertexColours[v] = vertexColours[u];
    }

    //keep track of bridges
    removeBridgesAtEnd(ppgraph, 1);
    insertBridge(ppgraph, changedBridges[degree1OperationsDepth][1], v, t);
    insertBridge(ppgraph, changedBridges[degree1OperationsDepth][0], u, ppgraph->adjList[u*3]);

    DEBUGPPGRAPHPRINT(ppgraph)
    DEBUGBRIDGESPRINT(ppgraph)

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
    int oldBridge = findBridge(ppgraph, u, v);
    int bridgeColour = -1;  //value assigned to avoid compiler warnings
    if(onlyColourable){
        bridgeColour = getColour(ppgraph, u, v);
    }
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

    //update colouring
    if(onlyColourable){
        int switchedColour = bridgeColour%3 + 1;
        int remainingColour = switchedColour%3 + 1;

        switchColours(ppgraph, u, switchedColour, bridgeColour);
        setColour(ppgraph, u, s, switchedColour);
        colours[s][0]=switchedColour;
        colours[s][1]=bridgeColour;
        colours[s][2]=remainingColour;
        colours[t][0]=remainingColour;

        coloursAroundVertex[s]=3;
        coloursAroundVertex[t]=1;
    }

    //update vertex colouring
    if(onlyBipartite){
        vertexColours[s]=vertexColours[u];
        vertexColours[t]=vertexColours[v];
        switchVertexColours(ppgraph, u, s);
    }

    //keep track of bridges
    changedBridges[degree1OperationsDepth][0] = oldBridge;
    numberOfChangedBridges[degree1OperationsDepth] = 1;
    removeBridge(ppgraph, oldBridge);
    appendBridge(ppgraph, u, s);
    appendBridge(ppgraph, v, s);
    appendBridge(ppgraph, t, s);

    DEBUGPPGRAPHPRINT(ppgraph)
    DEBUGBRIDGESPRINT(ppgraph)

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

    //restore colouring
    if(onlyColourable){
        int bridgeColour = colours[s][1];
        int switchedColour = colours[s][0];
        switchColours(ppgraph, u, bridgeColour, switchedColour);
        setColour(ppgraph, u, v, bridgeColour);
    }

    //restore vertex colouring
    if(onlyBipartite){
        switchVertexColours(ppgraph, u, v);
    }

    //keep track of bridges
    removeBridgesAtEnd(ppgraph, 3);
    insertBridge(ppgraph, changedBridges[degree1OperationsDepth][0], u, v);

    DEBUGPPGRAPHPRINT(ppgraph)
    DEBUGBRIDGESPRINT(ppgraph)

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

void get_bridges(PRIMPREGRAPH *ppgraph, VERTEXPAIR *vertexPairList, int *vertexPairListSize){
    int i;
    *vertexPairListSize = ppgraph->bridgeCount;
    for(i=0;i<ppgraph->bridgeCount;i++){
        vertexPairList[i][0]=ppgraph->bridges[i][1];
        vertexPairList[i][1]=ppgraph->bridges[i][0];
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
void determine_vertex_pairs_orbits(VERTEXPAIR *vertexPairList, int vertexPairListSize, int *vertexPairOrbits, int *orbitCount,
        permutation (*currentGenerators)[MAXN][MAXN] , int currentNumberOfGenerators){
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

    if(currentNumberOfGenerators==0){
        DEBUGARRAYDUMP(vertexPairOrbits, vertexPairListSize, "%d")
        DEBUGMSG("End determine_vertex_pairs_orbits")
        return;
    }

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
 * determines the orbits of a given list of vertices pairs. The pairs must be so that the image of a pair in the list
 * is also in the list.
 *
 * In the end vertexPairOrbits[i] will contain the index of the canonical representant of the orbit to which vertexPairList[i]
 * belongs and orbitCount will contain the number of orbits.
 * The first two parameters are read-only, the last two are write-only.
 */
void determine_vertex_pairs_orbits_and_sizes(VERTEXPAIR *vertexPairList, int vertexPairListSize, 
        int *vertexPairOrbits, int *orbitCount, permutation (*currentGenerators)[MAXN][MAXN] ,
        int currentNumberOfGenerators, int *orbitSize){
    DEBUGMSG("Start determine_vertex_pairs_orbits")
    DEBUG2DARRAYDUMP(vertexPairList, vertexPairListSize, 2, "%d")

    int i, j, k, temp;

    //initialization of the variables
    for(i=0; i<vertexPairListSize; i++){
        vertexPairOrbits[i]=i;
        orbitSize[i]=1;
    }
    *orbitCount=vertexPairListSize;

    if(currentNumberOfGenerators==0){
        DEBUGARRAYDUMP(vertexPairOrbits, vertexPairListSize, "%d")
        DEBUGMSG("End determine_vertex_pairs_orbits")
        return;
    }

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
    for(i = 0; i < numberOfGenerators[degree1OperationsDepth + degree2OperationsDepth]; i++) {
        //the generators were stored in the global variable generators by the method save_generators
        permutation = automorphismGroupGenerators[degree1OperationsDepth + degree2OperationsDepth][i];
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

//--------------------FILTERS-----------------------------------

int getThirdNeighbour(FILTERPREGRAPH *pregraph, int v, int n1, int n2){
    int i;
    for(i=0;i<3;i++){
	if((pregraph->adjList[v][i]==n1 && pregraph->adjList[v][(i+1)%3]==n2) ||
	    (pregraph->adjList[v][i]==n2 && pregraph->adjList[v][(i+1)%3]==n1)){
	    return pregraph->adjList[v][(i+2)%3];
	}
    }
    fprintf(stderr, "Error in 3rd neighbour: %d, %d, %d\n", v, n1, n2);
    exit(1);
}

//returns -1 if v, n1 and n2 are not in a square
int getSquare(FILTERPREGRAPH *pregraph, int v, int n1, int n2){
    int i,j;

    //first check for semi-edges
    if(n1 == pregraph->order || n2 == pregraph->order) return -1;

    //then check for multi-edges
    if(n1 == n2) return -1;

    for(i=0; i<3; i++){
	for(j=0; j<3; j++){
	    if(pregraph->adjList[n1][i]==pregraph->adjList[n2][j] &&
                    pregraph->adjList[n1][i]!=v &&
                    pregraph->adjList[n1][i]!=pregraph->order){
		return pregraph->adjList[n1][i];
	    }
	}
    }

    return -1;
}

//needs the square in the form
//      v1     v3
//       o----o--...
//       |    |
//       |    |
//       o----o--...
//      v2    v4
void finishChain(FILTERPREGRAPH *pregraph, boolean *removed, int v1, int v2, int v3,
        int v4, boolean *admissable, unsigned long long *adjMatrix){
    int prevUp, prevDown, currentUp, currentDown, nextUp, nextDown;

    prevUp = v1;
    prevDown = v2;
    currentUp = v3;
    currentDown = v4;
    removed[prevUp]=TRUE;
    removed[prevDown]=TRUE;

    while(TRUE){
	removed[currentUp]=TRUE;
	removed[currentDown]=TRUE;

	nextUp = getThirdNeighbour(pregraph, currentUp, currentDown, prevUp);
	nextDown = getThirdNeighbour(pregraph, currentDown, currentUp, prevDown);
        if(nextUp==nextDown && nextUp!=pregraph->order){
            *admissable=FALSE;
            return;
        }
	if(nextUp==pregraph->order || nextDown==pregraph->order){
	    //end of chain is reached because semi-edge is found
	    return;
	}
	if(nextUp==currentDown){
	    //end of chain is reached because multi-edge is found
	    return;
	}
	if(removed[nextUp] || removed[nextDown]){
	    //end of chain is reached because next vertices are already removed
	    return;
	}

	if(adjMatrix[nextUp] & (1<<nextDown)){
	    //still in the chain of square
	    prevUp = currentUp;
	    prevDown = currentDown;
	    currentUp = nextUp;
	    currentDown = nextDown;
	} else {
	    //end of chain reached because nextUp and nextDown aren't adjacent
	    return;
	}

    }
}

//needs the square in the form
//      v1     v3
//       o----o--...
//       |    |
//       |    |
//       o----o--...
//      v2    v4
void finishChainWithout(FILTERPREGRAPH *pregraph, boolean *removed, int v1, int v2,
        int v3, int v4, unsigned long long *adjMatrix){
    int prevUp, prevDown, currentUp, currentDown, nextUp, nextDown;

    prevUp = v1;
    prevDown = v2;
    currentUp = v3;
    currentDown = v4;
    removed[prevUp]=TRUE;
    removed[prevDown]=TRUE;

    while(TRUE){
	removed[currentUp]=TRUE;
	removed[currentDown]=TRUE;

	nextUp = getThirdNeighbour(pregraph, currentUp, currentDown, prevUp);
	nextDown = getThirdNeighbour(pregraph, currentDown, currentUp, prevDown);
	if(nextUp==pregraph->order || nextDown==pregraph->order){
	    //end of chain is reached because semi-edge is found
	    return;
	}
	if(nextUp==currentDown){
	    //end of chain is reached because multi-edge is found
	    return;
	}
	if(removed[nextUp] || removed[nextDown]){
	    //end of chain is reached because next vertices are already removed
	    return;
	}

	if(adjMatrix[nextUp] & (1<<nextDown)){
	    //still in the chain of square
	    prevUp = currentUp;
	    prevDown = currentDown;
	    currentUp = nextUp;
	    currentDown = nextDown;
	} else {
	    //end of chain reached because nextUp and nextDown aren't adjacent
	    return;
	}

    }
}

//needs the square in the form
//      v1     v2
//       o----o
//       |    |
//       |    |
//       o----o
//      v3    v4
void detectAndRemoveChain(FILTERPREGRAPH *pregraph, boolean *removed, int v1, int v2,
        int v3, int v4, boolean *admissable, unsigned long long *adjMatrix){
    int nrOfSquares = 1;
    int side1_a, side1_b, side2_a, side2_b;
    int n1, n2, n3, n4;
    int semiEdge = pregraph->order;

    //determine third neighbours for corner vertices
    n1 = getThirdNeighbour(pregraph, v1, v2, v3);
    n2 = getThirdNeighbour(pregraph, v2, v1, v4);
    n3 = getThirdNeighbour(pregraph, v3, v1, v4);
    n4 = getThirdNeighbour(pregraph, v4, v2, v3);

    //check for multi-edge
    if(n1==v2){
        finishChain(pregraph, removed, v1, v2, v3, v4, admissable, adjMatrix);
	return;
    } else if (n1==v3){
	finishChain(pregraph, removed, v1, v3, v2, v4, admissable, adjMatrix);
	return;
    } else if(n4==v2){
        finishChain(pregraph, removed, v4, v2, v3, v1, admissable, adjMatrix);
	return;
    } else if (n4==v3){
        finishChain(pregraph, removed, v4, v3, v2, v1, admissable, adjMatrix);
	return;
    }

    //check for semi-edges
    if(n1 == semiEdge && n2==semiEdge){
        finishChain(pregraph, removed, v1, v2, v3, v4, admissable, adjMatrix);
	return;
    } else if (n1==semiEdge && n3==semiEdge){
	finishChain(pregraph, removed, v1, v3, v2, v4, admissable, adjMatrix);
	return;
    } else if(n4== semiEdge && n2==semiEdge){
        finishChain(pregraph, removed, v4, v2, v3, v1, admissable, adjMatrix);
	return;
    } else if (n4==semiEdge && n3==semiEdge){
        finishChain(pregraph, removed, v4, v3, v2, v1, admissable, adjMatrix);
	return;
    }

    //check for semi-edges at opposite vertices
    //we know that there aren't any at neighbouring sites
    if((n1==semiEdge && n4==semiEdge) || (n2==semiEdge && n3==semiEdge)){
        removed[v1]=TRUE;
        removed[v2]=TRUE;
        removed[v3]=TRUE;
        removed[v4]=TRUE;
	return;
    }

    if(n1==semiEdge){
        //in this case we know that n2,n3 and n4 are't semi-edge vertices
        if(adjMatrix[n4] & (1<<n2)){
            //    v1   v2
            //   -o----o-
            //    |    |
            //    |    |
            //   -o----o-
            //    v3   v4
            side1_a = v1;
            side1_b = v3;
            side2_a = v2;
            side2_b = v4;
        } else {
            //    v1   v3
            //   -o----o-
            //    |    |
            //    |    |
            //   -o----o-
            //    v2   v4
            side1_a = v1;
            side1_b = v2;
            side2_a = v3;
            side2_b = v4;
        }
    } else if(n4==semiEdge){
        //in this case we know that n2,n3 and n1 are't semi-edge vertices
        if(adjMatrix[n1] & (1<<n2)){
            //    v1   v3
            //   -o----o-
            //    |    |
            //    |    |
            //   -o----o-
            //    v2   v4
            side1_a = v3;
            side1_b = v4;
            side2_a = v1;
            side2_b = v2;
        } else {
            //    v1   v2
            //   -o----o-
            //    |    |
            //    |    |
            //   -o----o-
            //    v3   v4
            side1_a = v2;
            side1_b = v4;
            side2_a = v1;
            side2_b = v3;
        }
    } else {
        //determine sides of chains
        if(adjMatrix[n1] & (1<<n2)){
            //    v1   v3
            //   -o----o-
            //    |    |
            //    |    |
            //   -o----o-
            //    v2   v4
            side1_a = v1;
            side1_b = v2;
            side2_a = v3;
            side2_b = v4;
        } else if(adjMatrix[n1] & (1<<n3)){
            //    v1   v2
            //   -o----o-
            //    |    |
            //    |    |
            //   -o----o-
            //    v3   v4
            side1_a = v1;
            side1_b = v3;
            side2_a = v2;
            side2_b = v4;
        } else if(adjMatrix[n4] & (1<<n2)){
            //    v1   v2
            //   -o----o-
            //    |    |
            //    |    |
            //   -o----o-
            //    v3   v4
            side1_a = v1;
            side1_b = v3;
            side2_a = v2;
            side2_b = v4;
        } else {
            //    v1   v3
            //   -o----o-
            //    |    |
            //    |    |
            //   -o----o-
            //    v2   v4
            side1_a = v1;
            side1_b = v2;
            side2_a = v3;
            side2_b = v4;
        }
    }

    //propagate chain of squares firts to side1 then to side2
    int prevUp, prevDown, currentUp, currentDown, nextUp, nextDown;

    prevUp = side2_a;
    prevDown = side2_b;
    currentUp = side1_a;
    currentDown = side1_b;
    removed[prevUp]=TRUE;
    removed[prevDown]=TRUE;

    int firstUp, firstDown;
    firstUp = prevUp;
    firstDown = prevDown;

    boolean inChain = TRUE;
    while(inChain){
	removed[currentUp]=TRUE;
	removed[currentDown]=TRUE;

	nextUp = getThirdNeighbour(pregraph, currentUp, currentDown, prevUp);
	nextDown = getThirdNeighbour(pregraph, currentDown, currentUp, prevDown);
	if(nextUp==semiEdge && nextDown==semiEdge){
	    //end of chain is reached because 2 semi-edges are found
	    //finish chain in opposite direction
	    finishChain(pregraph, removed, side1_a, side1_b, side2_a, side2_b,
                    admissable, adjMatrix);
	    return;
	}
	if(nextUp==currentDown){
	    //end of chain is reached because multi-edge is found
	    //finish chain in opposite direction
	    finishChain(pregraph, removed, side1_a, side1_b, side2_a, side2_b,
                    admissable, adjMatrix);
	    return;
	}
	if((nextUp!=semiEdge && removed[nextUp]) || (nextDown!=semiEdge && removed[nextDown])){
	    //end of chain is reached because next vertices are already removed
	    //finish chain in opposite direction
            if(nextUp == firstUp || nextUp == firstDown || nextDown == firstUp || nextDown == firstDown){
                //no need to finish chain
                *admissable = nrOfSquares % 2;
                return;
            } else {
                inChain = FALSE;
            }
            /*
	    finishChainWithout(pregraph, removed, side1_a, side1_b, side2_a,
                    side2_b);
            fprintf(stdout, "Squares: %d\n", nrOfSquares);
	    *admissable = nrOfSquares % 2;
	    return;
            */
            //inChain = FALSE;
	}

	if(nextUp==semiEdge || nextDown==semiEdge){
            //exactly one semi-edge is found, so end of chain is reached

            inChain = FALSE;
        } else if(adjMatrix[nextUp] & (1<<nextDown)){
	    //still in the chain of square
	    prevUp = currentUp;
	    prevDown = currentDown;
	    currentUp = nextUp;
	    currentDown = nextDown;
	    nrOfSquares++;
	} else {
	    //end of chain reached because nextUp and nextDown aren't adjacent
	    inChain = FALSE;
	}

    }
    
    int firstEndUp = currentUp;
    int firstEndDown = currentDown;

    //if we get here then one side of the chain is ended by an edge with
    //an ordinary vertex
    //now check the other side of the chain
    prevUp = side1_a;
    prevDown = side1_b;
    currentUp = side2_a;
    currentDown = side2_b;

    inChain = TRUE;
    while(inChain){
	removed[currentUp]=TRUE;
	removed[currentDown]=TRUE;

	nextUp = getThirdNeighbour(pregraph, currentUp, currentDown, prevUp);
	nextDown = getThirdNeighbour(pregraph, currentDown, currentUp, prevDown);
	if(nextUp==semiEdge && nextDown==semiEdge){
	    //end of chain is reached because 2 semi-edges are found
	    return;
	}
	if(nextUp==currentDown){
	    //end of chain is reached because multi-edge is found
	    return;
	}
	if((nextUp!=semiEdge && removed[nextUp]) || (nextDown!=semiEdge && removed[nextDown])){
	    //end of chain is reached because next vertices are already removed
	    *admissable = nrOfSquares % 2;
	    return;
	}

	if(nextUp==semiEdge || nextDown==semiEdge){
            //exactly one semi-edge is found, so end of chain is reached
            *admissable = nrOfSquares % 2;
            inChain = FALSE;
        } else if(adjMatrix[nextUp] & (1<<nextDown)){
	    //still in the chain of square
	    prevUp = currentUp;
	    prevDown = currentDown;
	    currentUp = nextUp;
	    currentDown = nextDown;
	    nrOfSquares++;
	} else {
	    //end of chain reached because nextUp and nextDown aren't adjacent
            if(nrOfSquares == 2 && (adjMatrix[nextUp] & (1<<firstEndUp)) && (adjMatrix[nextDown] & (1<<firstEndDown))){
                removed[nextUp] = TRUE;
                removed[nextDown] = TRUE;
                *admissable = TRUE;
            } else {
                *admissable = nrOfSquares % 2;
            }
	    inChain = FALSE;
	}
    }
}

void removeChainsOfSquares(FILTERPREGRAPH *pregraph, boolean *removed,
        boolean *admissable, unsigned long long *adjMatrix){
    int i, j;

    for(i=0; i<pregraph->order; i++){
	if(!removed[i]){
	    for(j=0;j<3;j++){
		int n1 = pregraph->adjList[i][(j+1)%3];
		int n2 = pregraph->adjList[i][(j+2)%3];
		if(!removed[n1] && !removed[n2] && n1!=pregraph->order &&
                        n2!=pregraph->order){
		    int oppositeCorner = getSquare(pregraph, i, n1, n2);
		    if(oppositeCorner!=-1 && oppositeCorner != pregraph->order){
//                        fprintf(stdout, "%d, %d, %d, %d\n\n", i+1, n1+1, n2+1, oppositeCorner+1);
			detectAndRemoveChain(pregraph, removed, i, n1, n2,
                                oppositeCorner, admissable, adjMatrix);
			if(!(*admissable)) return;
//                        writePregraphTableDebug(stdout, pregraph, removed);
			j=3; //exit the for
		    }
		}
	    }
	}
    }

}

//checks whether vertex v is incident with a multi-edge and removes that edge if so.
//assumes that v is not incident with two or more semi-edges
void checkAndRemoveMultiEdge(FILTERPREGRAPH *pregraph, boolean *removed, int v,
        boolean *admissable, unsigned long long *adjMatrix){
    if(pregraph->adjList[3*v+0]==pregraph->adjList[3*v+1]){
	removed[v]=TRUE;
	removed[pregraph->adjList[v][0]]=TRUE;
	*admissable = !(adjMatrix[pregraph->adjList[v][0]] &
                (1<<pregraph->adjList[v][2]));
    } else if(pregraph->adjList[v][0]==pregraph->adjList[v][2]){
	removed[v]=TRUE;
	removed[pregraph->adjList[v][0]]=TRUE;
	*admissable = !(adjMatrix[pregraph->adjList[v][0]] &
                (1<<pregraph->adjList[v][1]));
    } else if(pregraph->adjList[v][1]==pregraph->adjList[v][2]){
	removed[v]=TRUE;
	removed[pregraph->adjList[v][1]]=TRUE;
	*admissable = !(adjMatrix[pregraph->adjList[v][1]] &
                (1<<pregraph->adjList[v][0]));
    }
}

void removeMultiEdges(FILTERPREGRAPH *pregraph, boolean *removed, boolean *admissable,
                       unsigned long long *adjMatrix){
    int i,j,semiEdgeCount;
    for(i=0;i<pregraph->order;i++){
	if(!removed[i]){
	    semiEdgeCount=0;
	    for(j=0;j<3;j++){
		if(pregraph->adjList[i][j]==pregraph->order){
		    semiEdgeCount++;
		}
	    }
	    if(semiEdgeCount<2){
		checkAndRemoveMultiEdge(pregraph, removed, i, admissable, adjMatrix);
		if(!(*admissable)) return;
	    }
	}
    }
}

//assumes that no multi-edges are left
boolean checkRemainder(FILTERPREGRAPH *pregraph, boolean *removed){
    int i, j, vertex, degree[MAXN], neighbours[MAXN];
    #ifdef _DEBUG
    for(i=0; i<pregraph->order; i++){
        fprintf(stderr, "%d) %d, %d, %d   %s\n", i+1, pregraph->adjList[i][0] + 1,
                pregraph->adjList[i][1] + 1, pregraph->adjList[i][2] + 1,
                removed[i] ? "X" : "");
    }
    #endif

    //first check to see if there are no vertices of degree 1
    //We only cut at edges which have to have colour 1. If we removed
    //two edges at a vertex without removing the vertex, there were
    //conflicting colours at that vertex. (TODO: verify!!!)
    //We also check the number of 'real' neighbours: this should be less than 3
    //If a vertex has degree 2, than it should be incident with at least
    //one semi-edge
    for(i=0; i<pregraph->order; i++){
        if(!removed[i]){
            degree[i] = 0;
            neighbours[i] = 0;
            for(j=0; j<3; j++){
                if(pregraph->adjList[i][j]==pregraph->order ||
                        !removed[pregraph->adjList[i][j]]) degree[i]++;
                if(pregraph->adjList[i][j]!=pregraph->order &&
                        !removed[pregraph->adjList[i][j]]) neighbours[i]++;
            }
            if(degree[i]<2)
                return FALSE;
            if(neighbours[i]==3)
                return FALSE;
            if(degree[i]==2 && neighbours[i]==2)
                return FALSE;
        }
    }

    #ifdef _DEBUG
    fprintf(stderr, "    d  n\n");
    for(i=0; i<pregraph->order; i++){
        fprintf(stderr, "%d)  %d  %d\n", i, degree[i], neighbours[i]);
    }
    fprintf(stderr, "\n");
    #endif



    boolean visited[MAXN];
    int stack[3*MAXN], parent[MAXN];
    int stackSize = 0;

    for(i=0; i<pregraph->order; i++) visited[i] = FALSE;

    for(i=0; i<pregraph->order; i++){
	if(!(removed[i] || visited[i])){
	    stack[0]=i;
	    stackSize = 1;
	    parent[i]=-1;
            visited[i]=TRUE;
	    while(stackSize){
		vertex = stack[--stackSize];
                for(j=0;j<3;j++){
                    if(pregraph->adjList[vertex][j]!=parent[vertex] &&
                            pregraph->adjList[vertex][j]!=pregraph->order){
                        if(visited[pregraph->adjList[vertex][j]]) return FALSE;
                        if(!removed[pregraph->adjList[vertex][j]]){
                            stack[stackSize]=pregraph->adjList[vertex][j];
                            stackSize++;
                            parent[pregraph->adjList[vertex][j]]=vertex;
                            visited[pregraph->adjList[vertex][j]]=TRUE;
                        }
                    }
		}
	    }
	}
    }

    for(i=0; i<pregraph->order; i++){
        if(!removed[i] && degree[i]==2 && neighbours[i]==1){
            int currentVertex = i;
            int nextVertex = -1;  //value assigned to avoid compiler warnings
            int previousVertex;
            int vertexCount = 2;

            for(j=0;j<3;j++){
                if((pregraph->adjList[i][j]==pregraph->order &&
                        removed[pregraph->adjList[i][(j+1)%3]]) ||
                    (removed[pregraph->adjList[i][j]] &&
                        pregraph->adjList[i][(j+1)%3]==pregraph->order)){
                    nextVertex = pregraph->adjList[i][(j+2)%3];
                }
            }

            while(!(degree[nextVertex]==2 || (neighbours[nextVertex]==1 &&
                    degree[nextVertex]==3))){
                 previousVertex = currentVertex;
                 currentVertex = nextVertex;
                 nextVertex = getThirdNeighbour(pregraph, currentVertex,
                         previousVertex, pregraph->order);
                 vertexCount++;
            }
            if(degree[nextVertex]!=3 && vertexCount%2) return FALSE;
        }
    }

    return TRUE;
}

boolean existsVertexWithMultipleSemiEdges(FILTERPREGRAPH *pregraph){
    int i, j, semiEdgeCount;
    for(i=0; i<pregraph->order; i++){
	semiEdgeCount=0;
	for(j=0; j<3; j++){
	    if(pregraph->adjList[i][j]==pregraph->order) semiEdgeCount++;
	}
	if(semiEdgeCount>1) return TRUE;
    }
    return FALSE;
}

/**
 * Checks to see whether this pregraph has a 2-factor where each component is
 * the quotient of a 4-cycle.
 */
boolean isAdmissablePregraph(PRIMPREGRAPH *ppgraph){
    int i, j;
    boolean removed[MAXN];
    int old2new[MAXN];
    int new2old[MAXN];
    for(j=0; j<MAXN; j++){
        removed[j]=FALSE;
        old2new[j]=1111111111;
        new2old[j]=1111111111;
    }

    for(i=0,j=0; i<ppgraph->order; i++){
        if(ppgraph->degree[i]!=1){
            old2new[i]=j;
            new2old[j]=i;
            removed[j] = FALSE;
            j++;
        } else {
            old2new[i] = ppgraph->order - ppgraph->degree1Count;
        }
    }

    unsigned long long adjMatrix[MAXN];
    FILTERPREGRAPH fpregraph;
    fpregraph.order=ppgraph->order - ppgraph->degree1Count;

    for(i = 0; i < fpregraph.order; i++){
        adjMatrix[i]=0;
        if(ppgraph->degree[new2old[i]]==3){
            for(j=0; j<3; j++){
                fpregraph.adjList[i][j]=old2new[ppgraph->adjList[(new2old[i])*3+j]];
                adjMatrix[i]=adjMatrix[i] | (1<<fpregraph.adjList[i][j]);
            }
        } else {
            for(j=0; j<2; j++){
                fpregraph.adjList[i][j]=old2new[ppgraph->adjList[(new2old[i])*3+j]];
                adjMatrix[i]=adjMatrix[i] | (1<<fpregraph.adjList[i][j]);
            }
            fpregraph.adjList[i][2]=old2new[ppgraph->multiedge[new2old[i]]];
        }
    }

    //if all vertices are adjacent with at least one semi-edge, then the graph
    //has an admissable colouring if there are an even number of vertices
    i=0;
    while(i<fpregraph.order && (adjMatrix[i] & (1<<(fpregraph.order))))
        i++;
    if(i==fpregraph.order) return !(fpregraph.order%2) ||
            existsVertexWithMultipleSemiEdges(&fpregraph);

    boolean admissable = TRUE;

    removeChainsOfSquares(&fpregraph, removed, &admissable, adjMatrix);

    if(!admissable) return FALSE;

    removeMultiEdges(&fpregraph, removed, &admissable, adjMatrix);

    if(!admissable) return FALSE;

    return checkRemainder(&fpregraph, removed);
}

//needs the square in the form
//      v1     v2
//       o----o
//       |    |
//       |    |
//       o----o
//      v3    v4
void detectAndRemoveChainForC4(FILTERPREGRAPH *pregraph, boolean *removed, 
        int v1, int v2, int v3, int v4, boolean *illegalConfiguration,
        unsigned long *adjMatrix){
    int nrOfSquares = 1;
    int side1_a, side1_b, side2_a, side2_b;
    int n1, n2, n3, n4;
    int semiEdge = pregraph->order;
    //propagate chain of squares first to side1 then to side2
    int prevUp, prevDown, currentUp = -1, currentDown = -1, nextUp, nextDown;
    //values assigned to avoid compiler warnings
    boolean inChain;

    //determine third neighbours for corner vertices
    n1 = getThirdNeighbour(pregraph, v1, v2, v3);
    n2 = getThirdNeighbour(pregraph, v2, v1, v4);
    n3 = getThirdNeighbour(pregraph, v3, v1, v4);
    n4 = getThirdNeighbour(pregraph, v4, v2, v3);
    //fprintf(stdout, "%d %d %d %d\n",n1,n2,n3,n4);

    //check for semi-edges at opposite vertices
    if((n1==semiEdge && n4==semiEdge) || (n2==semiEdge && n3==semiEdge)){
        removed[v1]=TRUE;
        removed[v2]=TRUE;
        removed[v3]=TRUE;
        removed[v4]=TRUE;
	return;
    }

    //now check if one side is already closed

    //check for multi-edge
    if(n1==v2){
	//    v1   v3
	//    o----o-
	//   /|    |
	//   \|    |
	//    o----o-
	//    v2   v4
	side1_a = v1;
	side1_b = v2;
	side2_a = v3;
	side2_b = v4;
        removed[v1]=TRUE;
        removed[v2]=TRUE;
    } else if (n1==v3){
	//    v1   v2
	//    o----o-
	//   /|    |
	//   \|    |
	//    o----o-
	//    v3   v4
	side1_a = v1;
	side1_b = v3;
	side2_a = v2;
	side2_b = v4;
        removed[v1]=TRUE;
        removed[v3]=TRUE;
    } else if(n4==v2){
	//    v4   v3
	//    o----o-
	//   /|    |
	//   \|    |
	//    o----o-
	//    v2   v1
	side1_a = v4;
	side1_b = v2;
	side2_a = v3;
	side2_b = v1;
        removed[v4]=TRUE;
        removed[v2]=TRUE;
    } else if (n4==v3){
	//    v4   v2
	//    o----o-
	//   /|    |
	//   \|    |
	//    o----o-
	//    v3   v1
	side1_a = v4;
	side1_b = v3;
	side2_a = v2;
	side2_b = v1;
        removed[v4]=TRUE;
        removed[v3]=TRUE;
    } else

    //check for semi-edges
    if(n1 == semiEdge && n2==semiEdge){
	//    v1   v3
	//   _o----o-
	//    |    |
	//   _|    |
	//    o----o-
	//    v2   v4
	side1_a = v1;
	side1_b = v2;
	side2_a = v3;
	side2_b = v4;
        removed[v1]=TRUE;
        removed[v2]=TRUE;
    } else if (n1==semiEdge && n3==semiEdge){
	//    v1   v2
	//   _o----o-
	//    |    |
	//   _|    |
	//    o----o-
	//    v3   v4
	side1_a = v1;
	side1_b = v3;
	side2_a = v2;
	side2_b = v4;
        removed[v1]=TRUE;
        removed[v3]=TRUE;
    } else if(n4== semiEdge && n2==semiEdge){
	//    v4   v3
	//   _o----o-
	//    |    |
	//   _|    |
	//    o----o-
	//    v2   v1
	side1_a = v4;
	side1_b = v2;
	side2_a = v3;
	side2_b = v1;
        removed[v4]=TRUE;
        removed[v2]=TRUE;
    } else if (n4==semiEdge && n3==semiEdge){
	//    v4   v2
	//   _o----o-
	//    |    |
	//   _|    |
	//    o----o-
	//    v3   v1
	side1_a = v4;
	side1_b = v3;
	side2_a = v2;
	side2_b = v1;
        removed[v4]=TRUE;
        removed[v3]=TRUE;
    } else

    //both ends are open
    {
    if(n1==semiEdge){
        //in this case we know that n2,n3 and n4 are't semi-edge vertices
        if(adjMatrix[n4] & (1<<n2)){
            //    v1   v2
            //   -o----o-
            //    |    |
            //    |    |
            //   -o----o-
            //    v3   v4
            side1_a = v1;
            side1_b = v3;
            side2_a = v2;
            side2_b = v4;
        } else {
            //    v1   v3
            //   -o----o-
            //    |    |
            //    |    |
            //   -o----o-
            //    v2   v4
            side1_a = v1;
            side1_b = v2;
            side2_a = v3;
            side2_b = v4;
        }
    } else if(n4==semiEdge){
        //in this case we know that n2,n3 and n1 are't semi-edge vertices
        if(adjMatrix[n1] & (1<<n2)){
            //    v1   v3
            //   -o----o-
            //    |    |
            //    |    |
            //   -o----o-
            //    v2   v4
            side1_a = v3;
            side1_b = v4;
            side2_a = v1;
            side2_b = v2;
        } else {
            //    v1   v2
            //   -o----o-
            //    |    |
            //    |    |
            //   -o----o-
            //    v3   v4
            side1_a = v2;
            side1_b = v4;
            side2_a = v1;
            side2_b = v3;
        }
    } else {
        //determine sides of chains
        if(adjMatrix[n1] & (1<<n2)){
            //    v1   v3
            //   -o----o-
            //    |    |
            //    |    |
            //   -o----o-
            //    v2   v4
            side1_a = v1;
            side1_b = v2;
            side2_a = v3;
            side2_b = v4;
        } else if(adjMatrix[n1] & (1<<n3)){
            //    v1   v2
            //   -o----o-
            //    |    |
            //    |    |
            //   -o----o-
            //    v3   v4
            side1_a = v1;
            side1_b = v3;
            side2_a = v2;
            side2_b = v4;
        } else if(adjMatrix[n4] & (1<<n2)){
            //    v1   v2
            //   -o----o-
            //    |    |
            //    |    |
            //   -o----o-
            //    v3   v4
            side1_a = v1;
            side1_b = v3;
            side2_a = v2;
            side2_b = v4;
        } else {
            //    v1   v3
            //   -o----o-
            //    |    |
            //    |    |
            //   -o----o-
            //    v2   v4
            side1_a = v1;
            side1_b = v2;
            side2_a = v3;
            side2_b = v4;
        }
    }

    prevUp = side2_a;
    prevDown = side2_b;
    currentUp = side1_a;
    currentDown = side1_b;
    removed[prevUp]=TRUE;
    removed[prevDown]=TRUE;
    //fprintf(stdout, "next side: %d %d %d %d\n", prevUp, prevDown, currentUp, currentDown);

    int firstUp, firstDown;
    firstUp = prevUp;
    firstDown = prevDown;

    inChain = TRUE;
    while(inChain){
	removed[currentUp]=TRUE;
	removed[currentDown]=TRUE;

	nextUp = getThirdNeighbour(pregraph, currentUp, currentDown, prevUp);
	nextDown = getThirdNeighbour(pregraph, currentDown, currentUp, prevDown);
	if(nextUp==pregraph->order || nextDown==pregraph->order){
	    inChain = FALSE;
	} else if(nextUp==currentDown){
	    inChain = FALSE;
	} else if(removed[nextUp] || removed[nextDown]){
            if(nextUp == firstUp || nextUp == firstDown || nextDown == firstUp || nextDown == firstDown){
                //no need to finish chain
                *illegalConfiguration = !(nrOfSquares % 2);
                return;
            } else {
                inChain = FALSE;
            }
	    //end of chain is reached because next vertices are already removed
	} else if(adjMatrix[nextUp] & (1<<nextDown)){
	    //still in the chain of square
	    prevUp = currentUp;
	    prevDown = currentDown;
	    currentUp = nextUp;
	    currentDown = nextDown;
	    nrOfSquares++;
	} else {
	        //end of chain reached because nextUp and nextDown aren't adjacent
                inChain = FALSE;
        }

    }
    }
    //first end is finished
    int firstEndUp = currentUp;
    int firstEndDown = currentDown;

    //now check the other side of the chain
    prevUp = side1_a;
    prevDown = side1_b;
    currentUp = side2_a;
    currentDown = side2_b;
    //fprintf(stdout, "next side: %d %d %d %d\n", prevUp, prevDown, currentUp, currentDown);

    inChain = TRUE;
    while(inChain){
	removed[currentUp]=TRUE;
	removed[currentDown]=TRUE;

	nextUp = getThirdNeighbour(pregraph, currentUp, currentDown, prevUp);
	nextDown = getThirdNeighbour(pregraph, currentDown, currentUp, prevDown);
        //fprintf(stdout, "%d %d\n", nextUp, nextDown);
	if(nextUp==pregraph->order || nextDown==pregraph->order){
	    //end of chain is reached because at least 1 semi-edge is found
	    *illegalConfiguration = !(nrOfSquares % 2);
	    return;
	}
	if(nextUp==currentDown){
	    //end of chain is reached because multi-edge is found
	    *illegalConfiguration = !(nrOfSquares % 2);
	    return;
	}
	if(removed[nextUp] || removed[nextDown]){
	    //end of chain is reached because next vertices are already removed
	    *illegalConfiguration = !(nrOfSquares % 2);
	    return;
	}

	if(adjMatrix[nextUp] & (1<<nextDown)){
	    //still in the chain of square
	    prevUp = currentUp;
	    prevDown = currentDown;
	    currentUp = nextUp;
	    currentDown = nextDown;
	    nrOfSquares++;
	} else {
	    //end of chain reached because nextUp and nextDown aren't adjacent
            if(nrOfSquares == 2 && (adjMatrix[nextUp] & (1<<firstEndUp)) && (adjMatrix[nextDown] & (1<<firstEndDown))){
                removed[nextUp] = TRUE;
                removed[nextDown] = TRUE;
                *illegalConfiguration = FALSE;
            } else {
                *illegalConfiguration = !(nrOfSquares % 2);
            }
	    inChain = FALSE;
	}
    }
}

void removeChainsOfSquaresForC4(FILTERPREGRAPH *pregraph, boolean *removed, 
                      boolean *illegalConfiguration, unsigned long *adjMatrix){
    int i, j;

    for(i=0; i<pregraph->order; i++){
	if(!removed[i]){
	    for(j=0;j<3;j++){
		int n1 = pregraph->adjList[i][(j+1)%3];
		int n2 = pregraph->adjList[i][(j+2)%3];
		if(!removed[n1] && !removed[n2]){
		    int oppositeCorner = getSquare(pregraph, i, n1, n2);
		    //fprintf(stdout, "get square: %d %d ==> %d\n", i, j, oppositeCorner);
		    if(oppositeCorner!=-1){
			detectAndRemoveChainForC4(pregraph, removed, i, n1, n2, oppositeCorner, illegalConfiguration, adjMatrix);
			if((*illegalConfiguration)) return;
			j=3; //exit the for
		    }
		}
	    }
	    if(!removed[i]){
		//vertex i is not removed so it wasn contained in a square
		*illegalConfiguration= TRUE;
		return;
	    }
	}
    }

}

/**
 * Checks to see whether this pregraph has a 2-factor where each component is
 * a 4-cycle.
 */
boolean hasC4Cover(PRIMPREGRAPH *ppgraph){
    int i, j;
    boolean removed[MAXN];
    int old2new[MAXN];
    int new2old[MAXN];

    for(i=0,j=0; i<ppgraph->order; i++){
        if(ppgraph->degree[i]!=1){
            old2new[i]=j;
            new2old[j]=i;
            removed[j] = FALSE;
            j++;
        } else {
            old2new[i] = ppgraph->order - ppgraph->degree1Count;
        }
    }

    unsigned long adjMatrix[MAXN];
    FILTERPREGRAPH fpregraph;
    fpregraph.order=ppgraph->order - ppgraph->degree1Count;

    for(i = 0; i < fpregraph.order; i++){
        adjMatrix[i]=0;
        if(ppgraph->degree[new2old[i]]==3){
            for(j=0; j<3; j++){
                fpregraph.adjList[i][j]=old2new[ppgraph->adjList[(new2old[i])*3+j]];
                adjMatrix[i]=adjMatrix[i] | (1<<fpregraph.adjList[i][j]);
            }
        } else {
            for(j=0; j<2; j++){
                fpregraph.adjList[i][j]=old2new[ppgraph->adjList[(new2old[i])*3+j]];
                adjMatrix[i]=adjMatrix[i] | (1<<fpregraph.adjList[i][j]);
            }
            fpregraph.adjList[i][2]=old2new[ppgraph->multiedge[new2old[i]]];
        }
    }

    boolean illegalConfiguration = FALSE;

    removeChainsOfSquaresForC4(&fpregraph, removed, &illegalConfiguration, adjMatrix);

    if(illegalConfiguration) return FALSE;

    for(i=0; i<fpregraph.order; i++)
	if(!removed[i])
	    return FALSE;

    return TRUE;
}

//--------------------OUTPUT------------------------------------

char writePregraphTable(FILE *f, PREGRAPH *pregraph) {
    fprintf(f, "==============================\n");
    fprintf(f, "|  Graph number: %10llu  |\n", structureCount);
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
    DEBUGMSG("Checking for filters and applying")
    if(onlyAdmissable){
        if(!isAdmissablePregraph(pregraph->ppgraph)){
            DEBUGMSG("Pregraph refused by filter")
            return;
        }
        DEBUGMSG("Pregraph accepted by filter")
    } else if(onlyC4Coverable){
        if(!hasC4Cover(pregraph->ppgraph)){
            DEBUGMSG("Pregraph refused by filter")
            return;
        }
        DEBUGMSG("Pregraph accepted by filter")
    } else {
        DEBUGMSG("No filter set")
    }
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
    /* Handle splitting of construction */
    if(moduloEnabled && splitDepth>(degree1OperationsDepth+degree2OperationsDepth)){
        splitPointCount++;
        if(splitPointCount%moduloMod != moduloRest) {
            return;
        }
    }
    /**/
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
        if(loopCount==0){
            //all the degree 1 vertices correspond to semi-edges
            DEBUGASSERT(listSize==1)
            int vertex = -1;
            while((vertex = nextDegree1Vertex(vertex, ppgraph))!=-1){
                ADDELEMENT(vertexSetList, vertex);
            }
        } else {
            int position = 0;
            set tempSet[MAXM];
            EMPTYSET(tempSet, MAXM);
            determine_possible_sets_of_degree1_vertices(tempSet, vertexSetList,
                    &position, semiEdgeCount, 0, nextDegree1Vertex(-1, ppgraph), ppgraph, 0);
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
        }
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
    /* Handle splitting of construction */
    if(moduloEnabled && splitDepth==(degree1OperationsDepth)){
        splitPointCount++;
        if(splitPointCount%moduloMod != moduloRest) {
            degree1OperationsDepth--;
            return;
        }
    }
    /**/

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

    do_deg1_operations(ppgraph); //when this returns &ppgraph is unchanged

    if(allowMultiEdges && ppgraph->order - ppgraph->degree1Count + 2 <= vertexCount){
        do_deg2_operations(ppgraph, FALSE);
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
    /* Handle splitting of construction */
    if(moduloEnabled && splitDepth==(degree1OperationsDepth+degree2OperationsDepth)){
        splitPointCount++;
        if(splitPointCount%moduloMod != moduloRest) {
            degree2OperationsDepth--;
            return;
        }
    }
    /**/

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

    do_deg2_operations(ppgraph, multiEdgesDetermined);
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

boolean isCanonicalDegree1Edge(PRIMPREGRAPH *ppgraph, int v, boolean groupMayBeCopied){
    //groupMayBeCopied is TRUE only if the operation is 1.2 and the bridge was fixed by the automorphism group
    DEBUGASSERT(ppgraph->degree[v]==1)
    DEBUGDUMP(v,"%d")
    int i, j;

    //set to false if we reach the end of the generation and the group doesn't need to be calculated for the next round
    //in case we have both loops and semi-edges the group is always needed to create all sets.
    boolean groupstillNeeded = (allowLoops && allowSemiEdges) || (ppgraph->order - ppgraph->degree1Count + 1 <= vertexCount);

    #ifdef _PROFILING_DEG1
        canonicalDegree1Calls++;
    #endif

    if(ppgraph->degree1Count==1){
        #ifdef _PROFILING_DEG1
            canonicalDegree1BecauseOnlyOneVertexOfDegree1++;
        #endif
        //only one degree 1 vertex: garantueed to be canonical
        if(groupstillNeeded){
            if(numberOfGenerators[degree1OperationsDepth + degree2OperationsDepth]==0){
                //group was trivial and remains trivial: no need to call nauty
                numberOfGenerators[degree1OperationsDepth + degree2OperationsDepth + 1] = 0;
                #ifdef _PROFILING_DEG1
                    canonicalDegree1TrivialRemainsTrivial++;
                #endif
            } else if(groupMayBeCopied){
                #ifdef _PROFILING_DEG1
                    canonicalDegree1BridgeFixed++;
                #endif
                    numberOfGenerators[degree1OperationsDepth + degree2OperationsDepth + 1]=
                            numberOfGenerators[degree1OperationsDepth + degree2OperationsDepth];
                    for(i = 0; i < numberOfGenerators[degree1OperationsDepth + degree2OperationsDepth]; i++){
                        memcpy(automorphismGroupGenerators[degree1OperationsDepth + degree2OperationsDepth + 1] + i,
                               automorphismGroupGenerators[degree1OperationsDepth + degree2OperationsDepth] + i,
                               sizeof(permutation) * ppgraph->order);
                        automorphismGroupGenerators[degree1OperationsDepth + degree2OperationsDepth + 1][i][v-1] = v-1;
                        automorphismGroupGenerators[degree1OperationsDepth + degree2OperationsDepth + 1][i][v] = v;
                    }
            } else {
                //call nauty and return true
                int vertexOrbits[ppgraph->order];
                DEBUGMSG("Start nauty")
                numberOfGenerators[degree1OperationsDepth + degree2OperationsDepth + 1] = 0; //reset the generators
                nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, vertexOrbits, &nautyOptions, &nautyStats, nautyWorkspace,
                        NAUTY_WORKSIZE, MAXM, ppgraph->order, canonicalGraph);
                DEBUGMSG("End nauty")
                DEBUGARRAYDUMP(vertexOrbits, ppgraph->order, "%d")
            }
        }
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

    int minimumColour = colourDegree1VertexNeighbourhoodSize(ppgraph, DEG1_DISTANCE_COLOUR_VALUE, colours);

    //if v hasn't got the smallest colour, then it isn't canonical
    if(minimumColour != colours[v]) {
        #ifdef _PROFILING_DEG1
            canonicalDegree1NotBecauseNotSmallestColour++;
        #endif
        return FALSE;
    }

    int minimumColourCount = 0;
    for (i = 0; i < ppgraph->degree1Count; i++){
        if(minimumColour == colours[degree1Vertices[i]]){
            minimumColourCount++;
        }
    }

    if(minimumColourCount==1){
        //only one degree 1 vertex with minimal colour, i.e. v is canonical
        //call nauty and return true
        #ifdef _PROFILING_DEG1
            canonicalDegree1BecauseOnlyOneMinimumColour++;
        #endif
        if(groupstillNeeded){
            int vertexOrbits[ppgraph->order];
            DEBUGMSG("Start nauty")
            numberOfGenerators[degree1OperationsDepth + degree2OperationsDepth + 1] = 0; //reset the generators
            nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, vertexOrbits, &nautyOptions,
                    &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order, canonicalGraph);
            DEBUGMSG("End nauty")
            DEBUGARRAYDUMP(vertexOrbits, ppgraph->order, "%d")
        }
        return TRUE;
    }

    #ifdef _PROFILING_DEG1
        canonicalDegree1MinimumColourFrequency[minimumColour*MAXN + minimumColourCount]++;
    #endif

    //just call nauty and we'll be sure whether it is canonical

    //we provide nauty with the colouring we already obtained
    //first we set complete the partioning to 1
    for(j = 0; j < ppgraph->order; j++){
        nautyPtn[j]=1;
    }
    //next we provide a grouping of the degree 1 vertices based on their colour
    int currentColour = minimumColour;
    int labellingIndex = 0;
    int partitionCount = 0;
    while(labellingIndex < ppgraph->degree1Count){
        for(j=0; j < ppgraph->degree1Count; j++){
            if(colours[degree1Vertices[j]] == currentColour){
                if(labellingIndex == 0 || nautyPtn[labellingIndex - 1] == 0){
                    partitionCount++;
                }
                nautyLabelling[labellingIndex++] = degree1Vertices[j];
            }
        }
        nautyPtn[labellingIndex - 1]=0;
        currentColour++;
    }

    //determine colours for degree 3 vertices
    //given a degree 3 vertex v we have the following possibilities:
    // - v has 0 degree 1 vertices => d3colour[v]=0
    // - v has 1 degree 1 vertex w => d3colour[v]=d1colour[w] (>0)
    // - v has 2 degree 1 vertices w1 and w2 with d1colour[w1]=c1 and d1colour[w2]=c2 =>
    //      d3colour[v] = max(c1,c2) + canonicalDegree1PossibleColoursCount*min(c1,c2)
    //   because c1,c2 < canonicalDegree1PossibleColoursCount this will be different
    //   for each different pair c1 and c2
    // - v has 3 degree 1 vertices => there is only one degree 3 vertex and colouring
    //   makes no difference
    int d3Colours[MAXN];
    for(j=0; j < ppgraph->order; j++){
        d3Colours[j]=0;
    }
    for(j=0; j < ppgraph->degree1Count; j++){
        d3Colours[ppgraph->adjList[degree1Vertices[j]*3]] =
            (d3Colours[ppgraph->adjList[degree1Vertices[j]*3]] ? (1<<DEG1_DISTANCE_COLOUR_VALUE) + 1 : 0)
            + colours[degree1Vertices[j]];
    }

    int currentPartitionStart = labellingIndex;
    int lastColour = -1;
    int currentd3Colour = -1;

    #ifdef _PROFILING_DEG1
        int degree3PartitionCount = 0;
    #endif

    while(labellingIndex < ppgraph->order){
        for(j = 0; j < ppgraph->order; j++){
            if(ppgraph->degree[j]==3){
                if(currentd3Colour == lastColour){
                    if(d3Colours[j] > currentd3Colour){
                        currentd3Colour = d3Colours[j];
                        nautyLabelling[labellingIndex++] = j;
                    }
                } else if(d3Colours[j] < currentd3Colour && d3Colours[j] > lastColour){
                    currentd3Colour = d3Colours[j];
                    labellingIndex = currentPartitionStart;
                    nautyLabelling[labellingIndex++] = j;
                } else if(d3Colours[j] == currentd3Colour){
                    nautyLabelling[labellingIndex++] = j;
                }
            }
        }
        nautyPtn[labellingIndex - 1]=0;
        partitionCount++;
        currentPartitionStart = labellingIndex;
        lastColour = currentd3Colour;
        #ifdef _PROFILING_DEG1
            degree3PartitionCount++;
        #endif
        /*
        fprintf(stderr, "index: %d\n", labellingIndex);
        fprintf(stderr, "order: %d\n", ppgraph->order);
        fprintf(stderr, "current: %d\n", lastColour);
        fprintf(stderr, "last: %d\ncolouring :", currentd3Colour);
        for(i = 0; i < ppgraph->order; i++) fprintf(stderr, "%3d ", d3Colours[i]);
        fprintf(stderr, "\nlabelling :");
        for(i = 0; i < ppgraph->order; i++) fprintf(stderr, "%3d ", nautyLabelling[i]);
        fprintf(stderr, "\ndegree    :");
        for(i = 0; i < ppgraph->order; i++) fprintf(stderr, "%3d ", ppgraph->degree[i]);
        fprintf(stderr, "\n");
        */
    }

    DEBUGASSERT(labellingIndex == ppgraph->order)
    #ifdef _PROFILING_DEG1
    int lastPartitionEnd = -1;
    for (j = 0; j < ppgraph->order; j++){
        if(nautyPtn[j]==0){
            int size = j - lastPartitionEnd;
            if(ppgraph->degree[j]==1){
                int position = colours[j];
                canonicalDegree1Degree1PartitionCount[position]++;
                canonicalDegree1Degree1PartitionSize[position]+=size;
            } else {
                int position = d3Colours[j];
                if(position==0){
                    canonicalDegree1Degree3Neighbours0PartitionCount++;
                    canonicalDegree1Degree3Neighbours0PartitionSize+=size;
                } else if(position <= canonicalDegree1PossibleColoursCount){
                    canonicalDegree1Degree3Neighbours1PartitionCount[position]++;
                    canonicalDegree1Degree3Neighbours1PartitionSize[position]+=size;
                } else {
                    position-=(canonicalDegree1PossibleColoursCount+1);
                    canonicalDegree1Degree3Neighbours2PartitionCount[position]++;
                    canonicalDegree1Degree3Neighbours2PartitionSize[position]+=size;
                }
            }
            lastPartitionEnd = j;
        }
    }
        canonicalDegree1PartitionCountFrequency[partitionCount]++;
        canonicalDegree1Degree3PartitionCount[degree3PartitionCount]++;
    #endif

    nautyOptions.defaultptn = FALSE;
//    nautyOptions.invarproc = twopaths;
    int vertexOrbits[ppgraph->order];
    DEBUGMSG("Start nauty")
    numberOfGenerators[degree1OperationsDepth + degree2OperationsDepth + 1] = 0; //reset the generators
    nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, vertexOrbits, &nautyOptions,
            &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order, canonicalGraph);
    DEBUGMSG("End nauty")
    DEBUGARRAYDUMP(vertexOrbits, ppgraph->order, "%d")
//    nautyOptions.invarproc = NULL;
    nautyOptions.defaultptn = TRUE;

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
    #ifdef _PROFILING_DEG1
        if(smallestLabelOrbitV < smallestOtherDegree1Label)
            canonicalDegree1YesWithNauty++;
        else
            canonicalDegree1NoWithNauty++;
    #endif
    #ifdef _DEBUG_MODIFY_GENERATION
    if(noRejections) return TRUE;
    #endif
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
        nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, vertexOrbits,
                 &nautyOptions, &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order, canonicalGraph);
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
        nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, vertexOrbits, &nautyOptions,
              &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order, canonicalGraph);
        DEBUGMSG("End nauty")
        DEBUGARRAYDUMP(vertexOrbits, ppgraph->order, "%d")
        return TRUE;
    }

    //just call nauty and we'll be sure whether it is canonical

    int vertexOrbits[ppgraph->order];
    DEBUGMSG("Start nauty")
    numberOfGenerators = 0; //reset the generators
    nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, vertexOrbits, &nautyOptions,
            &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order, canonicalGraph);
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
void handle_deg1_operation1(PRIMPREGRAPH *ppgraph){
    #ifdef _PROFILING_OPERATION1_SAME_COLOUR
    degree1Operation1Degree1Counts[ppgraph->degree1Count]++;
    #endif
    #ifdef _DEBUG_MODIFY_GENERATION
    if(operation11Disabled) return;
    #endif
    DEBUGMSG("Start handle_deg1_operation1")
    DEBUG2DARRAYDUMP(automorphismGroupGenerators[degree1OperationsDepth + degree2OperationsDepth],
            numberOfGenerators[degree1OperationsDepth + degree2OperationsDepth], ppgraph->order, "%d")
    if(onlyColourable && ppgraph->degree1Count<3){
        #ifdef _PROFILING_OPERATION1_SAME_COLOUR
        degree1Operation1TooFewVertices++;
        #endif
        //garantueed to lead to a non-3-edge-colourable graph
        DEBUGMSG("End handle_deg1_operation1")
        return;
    }
    int maxSize = ppgraph->degree1Count*ppgraph->degree1Count/2;
    VERTEXPAIR deg1PairList[maxSize]; //initialize an array that is large enough to hold all the degree 1 pairs
    int listSize;
    get_deg1_pairs(ppgraph, deg1PairList, &listSize);

    int orbitCount;
    int orbits[listSize];
    determine_vertex_pairs_orbits(deg1PairList, listSize, orbits, &orbitCount,
            automorphismGroupGenerators + (degree1OperationsDepth + degree2OperationsDepth),
            numberOfGenerators[degree1OperationsDepth + degree2OperationsDepth]);

    int i;
    for (i = 0; i < listSize; i++) {
        if(orbits[i]==i){
            #ifdef _PROFILING_OPERATION1_SAME_COLOUR
                degree1Operation1OperationTried++;
            #endif
            DEBUGPPGRAPHPRINT(ppgraph)
            #ifdef _PROFILING_DEG1
                degree1Operation1Total++;
            #endif
            if(onlyColourable && colours[deg1PairList[i][0]][0]==colours[deg1PairList[i][1]][0]){
                #ifdef _PROFILING_OPERATION1_SAME_COLOUR
                    degree1Operation1HaveSameColour++;
                #endif
                //check that this graph can be coloured so that these two vertices of degree
                //1 are incident with two edges with a different colour
                if(!tryAlternateColouring(ppgraph, deg1PairList[i][0], deg1PairList[i][1])){
                    continue; //don't apply operation on these vertices
                }
                #ifdef _PROFILING_OPERATION1_SAME_COLOUR
                else {
                    degree1Operation1SameColourSolved++;
                }
                #endif
            }
            #ifdef _PROFILING_OPERATION1_SAME_COLOUR
            else {
                degree1Operation1HaveDifferentColour++;
            }
            #endif
            if(onlyBipartite && vertexColours[deg1PairList[i][0]]!=vertexColours[deg1PairList[i][1]]){
                continue; //don't apply operation on these vertices
            }
            apply_deg1_operation1(ppgraph, deg1PairList[i][0], deg1PairList[i][1]);

            //the only deg 1 vertex after this operation is v. This is a valid action
            //if v belongs to the first orbit of degree 1 vertices

            if(isCanonicalDegree1Edge(ppgraph, deg1PairList[i][1], FALSE)){
                #ifdef _PROFILING_DEG1
                    degree1Operation1Canonical++;
                #endif
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

void handle_deg1_operation2(PRIMPREGRAPH *ppgraph){
    #ifdef _DEBUG_MODIFY_GENERATION
    if(operation12Disabled) return;
    #endif
    if(ppgraph->bridgeCount==0) return;
    DEBUGMSG("Start handle_deg1_operation2")
    DEBUG2DARRAYDUMP(automorphismGroupGenerators[degree1OperationsDepth + degree2OperationsDepth],
            numberOfGenerators[degree1OperationsDepth + degree2OperationsDepth], ppgraph->order, "%d")
    //int maxSize = ppgraph->order*3/2-ppgraph->degree1Count;
            //this upper bound is not tight (it is tight in case of no degree 2 vertices?)
    VERTEXPAIR edgeList[ppgraph->bridgeCount]; //initialize an array that is large enough to hold all single edges
    int listSize;
    get_bridges(ppgraph, edgeList, &listSize);

    int orbitCount;
    int orbits[listSize];
    int orbitSize[listSize];
    determine_vertex_pairs_orbits_and_sizes(edgeList, listSize, orbits, &orbitCount,
            automorphismGroupGenerators+(degree1OperationsDepth + degree2OperationsDepth),
            numberOfGenerators[degree1OperationsDepth + degree2OperationsDepth], orbitSize);

    DEBUGARRAYDUMP(orbits, orbitCount, "%d")
    DEBUGARRAYDUMP(orbitSize, orbitCount, "%d")

    int i;
    for (i = 0; i < listSize; i++) {
        if(orbits[i]==i){
            //TODO: maybe only enumerate the bridges?
            if(isBridge(ppgraph, edgeList[i][0], edgeList[i][1])){
                #ifdef _PROFILING_DEG1
                    degree1Operation2Total++;
                #endif
                DEBUGPPGRAPHPRINT(ppgraph)
                apply_deg1_operation2(ppgraph, edgeList[i][0], edgeList[i][1]);

                //the new deg 1 vertex after this operation is t. This is a valid action
                //if t belongs to the first orbit of degree 1 vertices

                if(isCanonicalDegree1Edge(ppgraph, ppgraph->order-1, orbitSize[i]==1)){
                    #ifdef _PROFILING_DEG1
                        degree1Operation2Canonical++;
                    #endif
                    //t belongs to the orbit of degree 1 vertices with the smallest representant
                    //Therefore this graph was created from the correct parent.

                    //TODO: if onlyColourable -> switch colours along a path on one side of the bridge
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
    numberOfGenerators[degree1OperationsDepth + degree2OperationsDepth + 1] = 0; //reset the generators
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

    nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, orbits, &nautyOptions,
            &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order + ppgraph->multiEdgeCount, canonicalGraph);
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

    determine_vertex_pairs_orbits(multiEdgeList, *multiEdgeListSize, multiEdgeOrbits,
            multiEdgeOrbitCount, automorphismGroupGenerators + (degree1OperationsDepth + degree2OperationsDepth + 1),
            numberOfGenerators[degree1OperationsDepth + degree2OperationsDepth + 1]);

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
        neighbourV = (ppgraph->adjList[(multiEdgeList[i][0])*3] == multiEdgeList[i][1]) ?
            ppgraph->adjList[(multiEdgeList[i][0])*3 + 1] : ppgraph->adjList[(multiEdgeList[i][0])*3];
        neighbourU = (ppgraph->adjList[(multiEdgeList[i][1])*3] == multiEdgeList[i][0]) ?
            ppgraph->adjList[(multiEdgeList[i][1])*3 + 1] : ppgraph->adjList[(multiEdgeList[i][1])*3];
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
    #ifdef _DEBUG_MODIFY_GENERATION
    if(noRejections) return TRUE;
    #endif
    return smallestRepresentantNewEdge[0] < smallestOtherMultiEdge[0] ||
            (smallestRepresentantNewEdge[0] == smallestOtherMultiEdge[0] && smallestRepresentantNewEdge[1] < smallestOtherMultiEdge[1]);
}

void handle_deg2_operation1(PRIMPREGRAPH *ppgraph){
    #ifdef _DEBUG_MODIFY_GENERATION
    if(operation21Disabled) return;
    #endif
    DEBUGMSG("Start handle_deg2_operation1")
    DEBUG2DARRAYDUMP(automorphismGroupGenerators[degree1OperationsDepth+degree2OperationsDepth],
            numberOfGenerators[degree1OperationsDepth+degree2OperationsDepth], ppgraph->order, "%d")
    int maxSize = ppgraph->order*3/2-ppgraph->degree1Count; //this upper bound is not tight (it is tight in case of no degree 2 vertices?)
    VERTEXPAIR edgeList[maxSize]; //initialize an array that is large enough to hold all single edges
    int listSize;
    get_single_edges(ppgraph, edgeList, &listSize);

    int orbitCount;
    int orbits[listSize];
    determine_vertex_pairs_orbits(edgeList, listSize, orbits, &orbitCount, 
            automorphismGroupGenerators+(degree1OperationsDepth+degree2OperationsDepth),
            numberOfGenerators[degree1OperationsDepth+degree2OperationsDepth]);

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

void handle_deg2_operation2(PRIMPREGRAPH *ppgraph, boolean *multiEdgesDetermined){
    #ifdef _DEBUG_MODIFY_GENERATION
    if(operation22Disabled) return;
    #endif
    DEBUGMSG("Start handle_deg2_operation2")
    DEBUG2DARRAYDUMP(automorphismGroupGenerators[degree1OperationsDepth+degree2OperationsDepth],
            numberOfGenerators[degree1OperationsDepth+degree2OperationsDepth], ppgraph->order, "%d")

    VERTEXPAIR *oldMultiEdgeList = globalMultiEdgeList + (degree2OperationsDepth*HALFFLOOR(vertexCount));
    int *oldMultiEdgeListSize = globalMultiEdgeListSize + degree2OperationsDepth;
    int *oldMultiEdgeOrbits = globalMultiEdgeOrbits + (degree2OperationsDepth*HALFFLOOR(vertexCount));
    int *oldMultiEdgeOrbitCount = globalMultiEdgeOrbitCount + degree2OperationsDepth;

    if(!(*multiEdgesDetermined)){
        //if the multi-edges have't been determined, do it now and store the result
        get_multi_edges(ppgraph, oldMultiEdgeList, oldMultiEdgeListSize);

        DEBUGASSERT(*oldMultiEdgeListSize == ppgraph->multiEdgeCount)

        determine_vertex_pairs_orbits(oldMultiEdgeList, *oldMultiEdgeListSize, oldMultiEdgeOrbits,
                oldMultiEdgeOrbitCount, automorphismGroupGenerators+(degree1OperationsDepth+degree2OperationsDepth),
                numberOfGenerators[degree1OperationsDepth+degree2OperationsDepth]);
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

void handle_deg2_operation3(PRIMPREGRAPH *ppgraph, boolean *multiEdgesDetermined){
    #ifdef _DEBUG_MODIFY_GENERATION
    if(operation23Disabled) return;
    #endif
    DEBUGMSG("Start handle_deg2_operation3")
    DEBUG2DARRAYDUMP(automorphismGroupGenerators[degree1OperationsDepth+degree2OperationsDepth],
            numberOfGenerators[degree1OperationsDepth+degree2OperationsDepth], ppgraph->order, "%d")

    VERTEXPAIR *oldMultiEdgeList = globalMultiEdgeList + (degree2OperationsDepth*HALFFLOOR(vertexCount));
    int *oldMultiEdgeListSize = globalMultiEdgeListSize + degree2OperationsDepth;
    int *oldMultiEdgeOrbits = globalMultiEdgeOrbits + (degree2OperationsDepth*HALFFLOOR(vertexCount));
    int *oldMultiEdgeOrbitCount = globalMultiEdgeOrbitCount + degree2OperationsDepth;

    if(*multiEdgesDetermined){
         //if the multi-edges have't been determined, do it now and store the result
        get_multi_edges(ppgraph, oldMultiEdgeList, oldMultiEdgeListSize);

        DEBUGASSERT(*oldMultiEdgeListSize == ppgraph->multiEdgeCount)

        determine_vertex_pairs_orbits(oldMultiEdgeList, *oldMultiEdgeListSize, oldMultiEdgeOrbits, 
                oldMultiEdgeOrbitCount, automorphismGroupGenerators+(degree1OperationsDepth+degree2OperationsDepth),
                numberOfGenerators[degree1OperationsDepth+degree2OperationsDepth]);
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
void do_deg1_operations(PRIMPREGRAPH *ppgraph){
    DEBUGMSG("Start do_deg1_operations")
    DEBUGASSERT(allowLoops || allowSemiEdges)
    if(ppgraph->order<=maxVertexCount) handle_deg1_operation1(ppgraph);
    if(ppgraph->order<=maxVertexCount-2) handle_deg1_operation2(ppgraph);
    DEBUGMSG("End do_deg1_operations")
}

/*
 * Performs the different degree 2 operations. When this method returns &ppgraph
 * will be unchanged.
 */
void do_deg2_operations(PRIMPREGRAPH *ppgraph, boolean multiEdgesDetermined){
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
        handle_deg2_operation1(ppgraph);
        //the orbits of the multi-edges are already determined, so we pass them on
        handle_deg2_operation2(ppgraph, &newMultiEdgesDetermined);
        if(!onlyColourable && !onlyBipartite){
            //don't do operation 2.3 in case we want 3-edge-colourable pregraphs
            //don't do operation 2.3 in case we want bipartite pregraphs
            handle_deg2_operation3(ppgraph, &newMultiEdgesDetermined);
        }
    }
    DEBUGMSG("End do_deg2_operations")
}

/*
 * Start the generation process with the given graph. At this points the bridges need to be calculated.
 */
void grow(PRIMPREGRAPH *ppgraph){
    /* Handle splitting of construction */
    if(splitDepth==0){
        splitPointCount++;
        if(moduloEnabled && (splitPointCount%moduloMod != moduloRest)) {
            return;
        }
    }
    /**/

    DEBUGMSG("Start grow")
    int orbits[ppgraph->order];
    DEBUGMSG("Start nauty")
    numberOfGenerators[0] = 0; //reset the generators
    //there are no multiedges at this point!
    degree1OperationsDepth=-1;
    nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, orbits, &nautyOptions,
            &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order, canonicalGraph);
    degree1OperationsDepth=0;
    DEBUGMSG("End nauty")
    //the generators for these start graphs need to be calculated

    if(ppgraph->order >= minVertexCount && ppgraph->order<=maxVertexCount && ppgraph->order - vertexCount <= ppgraph->degree1Count)
        handle_primpregraph_result(ppgraph);


    if(allowLoops || allowSemiEdges){
        determineBridges(ppgraph);
        do_deg1_operations(ppgraph);
    }
    if(allowMultiEdges){
        do_deg2_operations(ppgraph, FALSE);
    }
    DEBUGMSG("End grow")
}

void growWithoutDeg1Operations(PRIMPREGRAPH *ppgraph){
    /* Handle splitting of construction */
    if(splitDepth==0){
        splitPointCount++;
        if(moduloEnabled && (splitPointCount%moduloMod != moduloRest)) {
            return;
        }
    }
    /**/

    DEBUGMSG("Start growWithoutDeg1Operations")
    int orbits[ppgraph->order];

    VERTEXPAIR *multiEdgeList = globalMultiEdgeList; //depth = 0
    int *multiEdgeListSize = globalMultiEdgeListSize;
    get_multi_edges(ppgraph, multiEdgeList, multiEdgeListSize);

    DEBUGASSERT(*multiEdgeListSize == ppgraph->multiEdgeCount)

    DEBUGMSG("Start nauty")
    numberOfGenerators[0] = 0; //reset the generators
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

    degree1OperationsDepth=-1;
    nauty(ppgraph->ulgraph, nautyLabelling, nautyPtn, NULL, orbits, &nautyOptions,
            &nautyStats, nautyWorkspace, NAUTY_WORKSIZE, MAXM, ppgraph->order + ppgraph->multiEdgeCount, canonicalGraph);
    degree1OperationsDepth=0;

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

    int *multiEdgeOrbits = globalMultiEdgeOrbits;
    int *multiEdgeOrbitCount = globalMultiEdgeOrbitCount;
    determine_vertex_pairs_orbits(multiEdgeList, *multiEdgeListSize, multiEdgeOrbits, 
            multiEdgeOrbitCount, automorphismGroupGenerators+(degree1OperationsDepth+degree2OperationsDepth),
            numberOfGenerators[degree1OperationsDepth+degree2OperationsDepth]);

    if(ppgraph->order >= minVertexCount && ppgraph->order<=maxVertexCount && ppgraph->order - vertexCount <= ppgraph->degree1Count)
        handle_primpregraph_result(ppgraph);

    if(allowMultiEdges){
        do_deg2_operations(ppgraph, TRUE);
    }
    DEBUGMSG("End growWithoutDeg1Operations")
}

static PRIMPREGRAPH *currentPpgraph;

/*
 * Handles the input from the external graph by translating the graph into a pregraph primitive
 * and feeding it to the grow method.
 */
 void handle_snarkhunter_result(unsigned char snarkhunter_graph[MAXN][REG + 1], int order){
    DEBUGMSG("Start handle_snarkhunter_result")
    PRIMPREGRAPH *ppgraph = currentPpgraph;
    ppgraph->order = order;
    ppgraph->degree1Count = 0;
    ppgraph->multiEdgeCount = 0;
    ppgraph->bridgeCount = 0; //will be set when necessary

    int i, j;

    for(i=0; i<order; i++){
        ppgraph->degree[i]=3;

        set *v;
        v = GRAPHROW(ppgraph->ulgraph, i, MAXM);
        EMPTYSET(v, MAXM);
        j=0;
        for(j=0; j<3; j++){
            ppgraph->adjList[i*3+j] = snarkhunter_graph[i][j];
            ADDELEMENT(v, snarkhunter_graph[i][j]);
        }
    }

    DEBUGPPGRAPHPRINT(ppgraph)

    if(onlyColourable){
        DEBUGMSG("Checking colourability")
        if(!isColourableGraph(ppgraph)){
            DEBUGMSG("Cubic graph is not 3-edge-colourable")
            DEBUGMSG("End handle_snarkhunter_result")
            return;
        }
    }

    if(onlyBipartite){
        DEBUGMSG("Checking vertex-colourability")
        if(!isBipartiteGraph(ppgraph)){
            DEBUGMSG("Cubic graph is not bipartite")
            DEBUGMSG("End handle_snarkhunter_result")
            return;
        }
    }

    if(ppgraph->order == vertexCount)
        handle_primpregraph_result(ppgraph);
	else
		grow(ppgraph);

    DEBUGMSG("End handle_snarkhunter_result")

}

/*
 * Handles the input from the external graph by translating the graph into a pregraph primitive
 * and feeding it to the grow method.
 */
 void handle_minibaum_result(unsigned char minibaum_graph[knoten+1][reg], int order){
    DEBUGMSG("Start handle_minibaum_result")
    PRIMPREGRAPH *ppgraph = currentPpgraph;
    ppgraph->order = order;
    ppgraph->degree1Count = 0;
    ppgraph->multiEdgeCount = 0;
    ppgraph->bridgeCount = 0; //will be set when necessary

    int i, j;

    for(i=0; i<order; i++){
        ppgraph->degree[i]=3;

        set *v;
        v = GRAPHROW(ppgraph->ulgraph, i, MAXM);
        EMPTYSET(v, MAXM);
        j=0;
        for(j=0; j<3; j++){
            ppgraph->adjList[i*3+j] = minibaum_graph[i+1][j]-1;
            ADDELEMENT(v, minibaum_graph[i+1][j]-1);
        }
    }

    DEBUGPPGRAPHPRINT(ppgraph)

    if(onlyColourable){
        DEBUGMSG("Checking colourability")
        if(!isColourableGraph(ppgraph)){
            DEBUGMSG("Cubic graph is not 3-edge-colourable")
            DEBUGMSG("End handle_minibaum_result")
            return;
        }
    }

    if(onlyBipartite){
        DEBUGMSG("Checking vertex-colourability")
        //still need to call this method because it will set the colours which we need
        if(!isBipartiteGraph(ppgraph)){
            DEBUGMSG("Cubic graph is not bipartite")
            DEBUGMSG("End handle_minibaum_result")
            return;
        }
    }

    if(ppgraph->order == vertexCount)
        handle_primpregraph_result(ppgraph);
    else
        grow(ppgraph);

    DEBUGMSG("End handle_minibaum_result")

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
        if(!onlyColourable && !onlyBipartite){
            construct_K3_with_spike(&ppgraph);
            growWithoutDeg1Operations(&ppgraph);
        }
    }

    int i;
    currentPpgraph = &ppgraph;
    int cubicGeneratorStart = 4;
    if(!allowLoops && !allowSemiEdges && !allowMultiEdges){
        cubicGeneratorStart = vertexCount;
    }
    if(!onlyColourable || vertexCount%2!=1){
        if(is_minibaum_available(vertexCount) && onlyBipartite){
            #ifdef _DEBUG
            fprintf(stderr, "Starting minibaum for %d vertices\n", vertexCount);
            #endif
            call_minibaum(vertexCount, onlyBipartite, vertexCount!=cubicGeneratorStart);
        } else {
            for(i = cubicGeneratorStart; i <= vertexCount; i+=2){//TODO: is this the correct upperbound for i
                #ifdef _DEBUG
                fprintf(stderr, "Starting snarkhunter for %d vertices\n", i);
                #endif
                call_snarkhunter(i, 3, *handle_snarkhunter_result);
            }
        }
    }
    if(allowMultiEdges && vertexCount == 2){
        //is 3-edge-colourable
        //is bipartite
        //is admissable
        writeThetaGraph();
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
void writeThetaGraph(){
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
    } else if(outputType == 'c'){
        if (structureCount == 1) { //if first graph
            fprintf(file, ">>pregraph_code %s<<", (endian == LITTLE_ENDIAN ? "le" : "be"));
        }
        fprintf(file, "%c", 2);
        fprintf(file, "%c", 2);
        fprintf(file, "%c", 2);
        fprintf(file, "%c", 2);
        fprintf(file, "%c", 0);
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
    if(onlyColourable){
        colours[0][0]=1;
        colours[1][0]=1;
    }
    if(onlyBipartite){
        vertexColours[0]=BLACK;
        vertexColours[1]=WHITE;
    }
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
 *
 * not 3-edge-colourable
 * not bipartite
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
    //depth + 1, because we always call nauty for the new graph before the depth is increased
    memcpy(
            automorphismGroupGenerators[degree1OperationsDepth +
                                        degree2OperationsDepth + 1]
            + numberOfGenerators[degree1OperationsDepth +
                                 degree2OperationsDepth + 1]
           , perm, sizeof(permutation) * n);

    numberOfGenerators[degree1OperationsDepth + degree2OperationsDepth + 1]++;
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
    fprintf(stderr, "  -l          : Print generation limits and return.\n");
    fprintf(stderr, "  -i          : Causes %s to print extra info about the generated\n", name);
    fprintf(stderr, "                structures.\n");
    fprintf(stderr, "  -L          : Allow loops.\n");
    fprintf(stderr, "  -S          : Allow semi-edges.\n");
    fprintf(stderr, "  -M          : Allow multi-edges.\n");
    fprintf(stderr, "  -C          : Only generate 3-edge-colourable pregraphs.\n");
    fprintf(stderr, "  -B          : Only generate bipartite pregraphs.\n");
    fprintf(stderr, "  -q          : Only generate pregraphs that have a 2-factor where each\n");
    fprintf(stderr, "                component is the quotient of a 4-cycle.\n");
    fprintf(stderr, "  -4          : Only generate pregraphs that have a 2-factor where each\n");
    fprintf(stderr, "                component is a 4-cycle.\n");
    fprintf(stderr, "  -P          : Only generate the corresponding pregraph primitives.\n");
    fprintf(stderr, "  -I          : Only start the generation from the files provided by the input\n");
    fprintf(stderr, "                file (see -F).\n");
    fprintf(stderr, "  -f file     : Specifies the output file. If absent, the output is written to\n");
    fprintf(stderr, "                standard out.\n");
    fprintf(stderr, "  -F file     : Specifies the input file. If absent, the input is taken from\n");
    fprintf(stderr, "                standard in.\n");
    fprintf(stderr, "                This option is only used if -I is used.\n");
    fprintf(stderr, "  -o c        : Specifies the export format where c is one of\n");
    fprintf(stderr, "                c    pregraph code (or multicode if -P is used)\n");
    fprintf(stderr, "                h    human-readable output in tabular format\n");
    fprintf(stderr, "                n    no output: only count (default)\n");
#ifdef _DEBUG_MODIFY_GENERATION
    fprintf(stderr, "  -D #        : Disable some operations:\n");
    fprintf(stderr, "                1    Disable operation 1.1\n");
    fprintf(stderr, "                2    Disable operation 1.2\n");
    fprintf(stderr, "                3    Disable operation 2.1\n");
    fprintf(stderr, "                4    Disable operation 2.2\n");
    fprintf(stderr, "                5    Disable operation 2.3\n");
    fprintf(stderr, "  -d #        : Disable some operations:\n");
    fprintf(stderr, "                1    Disable operations for degree 1\n");
    fprintf(stderr, "                2    Disable operations for degree 2\n");
    fprintf(stderr, "  -X          : No operation will be discarded as being not-canonical. This will\n");
    fprintf(stderr, "                cause isomorphic graphs to be constructed. (This option is for\n");
    fprintf(stderr, "                debugging purposes.)\n");
#endif
    fprintf(stderr, "  -m r:m[:d]  : Split the generation into several parts. This basically means\n");
    fprintf(stderr, "                that at depth d a counter will be kept and the program will\n");
    fprintf(stderr, "                only continue beyond this point if the counter mod m is equal\n");
    fprintf(stderr, "                to r. Special measures are taken if some graphs are already\n");
    fprintf(stderr, "                outputted before depth d. The default for d is 0.\n");
}

void printLimits(char *name) {
    fprintf(stderr, "The program %s was compiled with MAXN set to %d.\n", name, MAXN);
    fprintf(stderr, "This imposes the following bounds on the generated structures:\n\n");
    fprintf(stderr, "simple cubic graphs (no options)                   : max %d vertices.\n", MAXN);
    fprintf(stderr, "simple cubic graphs with loops (-L)                : max %d vertices.\n", MAXN);
    fprintf(stderr, "simple cubic graphs with semi-edges (-S)           : max %d vertices.\n", (MAXN-2)/2);
    fprintf(stderr, "cubic multigraphs (-M)                             : max %d vertices.\n", 2*MAXN/3);
    fprintf(stderr, "cubic multigraphs with loops (-LM)                 : max %d vertices.\n", 2*MAXN/3);
    fprintf(stderr, "cubic multigraphs with semi-edges (-SM)            : max %d vertices.\n", (MAXN-2)/2);
    fprintf(stderr, "cubic multigraphs with loops and semi-edges (-LSM) : max %d vertices.\n\n", (MAXN-2)/2);
    fprintf(stderr, "The options -C, -B, -q and -4 do not change the maximum number of vertices.\n");
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

    #ifdef _PROFILING_DEG1
        degree1Operation1Total = 0;
        degree1Operation1Canonical = 0;
        degree1Operation2Total = 0;
        degree1Operation2Canonical = 0;

        canonicalDegree1Calls = 0;
        canonicalDegree1BecauseOnlyOneVertexOfDegree1 = 0;
        canonicalDegree1TrivialRemainsTrivial = 0;
        canonicalDegree1BridgeFixed = 0;
        canonicalDegree1NotBecauseNotSmallestColour = 0;
        canonicalDegree1BecauseOnlyOneMinimumColour = 0;
        canonicalDegree1YesWithNauty = 0;
        canonicalDegree1NoWithNauty = 0;

        canonicalDegree1PossibleColoursCount = (1<<DEG1_DISTANCE_COLOUR_VALUE) + 1;
        canonicalDegree1MinimumColourFrequency = (unsigned long long *)malloc(sizeof(unsigned long long)*(canonicalDegree1PossibleColoursCount * MAXN));
        for(i = 0; i < canonicalDegree1PossibleColoursCount * MAXN; i++) canonicalDegree1MinimumColourFrequency[i]=0;

        canonicalDegree1PartitionCountFrequency = (int *)malloc(sizeof(int)*(canonicalDegree1PossibleColoursCount + 2));
        for(i = 0; i < canonicalDegree1PossibleColoursCount + 2; i++) canonicalDegree1PartitionCountFrequency[i]=0;

        for(i = 0; i <= MAXN; i++) canonicalDegree1Degree3PartitionCount[i]=0;

        canonicalDegree1Degree1PartitionCount = (unsigned long long *)malloc(sizeof(unsigned long long)*(canonicalDegree1PossibleColoursCount+2));
        canonicalDegree1Degree1PartitionSize = (unsigned long long *)malloc(sizeof(unsigned long long)*(canonicalDegree1PossibleColoursCount+2));
        canonicalDegree1Degree3Neighbours1PartitionCount = (unsigned long long *)malloc(sizeof(unsigned long long)*(canonicalDegree1PossibleColoursCount+2));
        canonicalDegree1Degree3Neighbours1PartitionSize = (unsigned long long *)malloc(sizeof(unsigned long long)*(canonicalDegree1PossibleColoursCount+2));
        canonicalDegree1Degree3Neighbours2PartitionCount = (unsigned long long *)malloc(sizeof(unsigned long long)*(canonicalDegree1PossibleColoursCount+2));
        canonicalDegree1Degree3Neighbours2PartitionSize = (unsigned long long *)malloc(sizeof(unsigned long long)*(canonicalDegree1PossibleColoursCount+2));

        for(i = 0; i < canonicalDegree1PossibleColoursCount + 2; i++){
            canonicalDegree1Degree1PartitionCount[i] = 0;
            canonicalDegree1Degree1PartitionSize[i] = 0;
            canonicalDegree1Degree3Neighbours1PartitionCount[i] = 0;
            canonicalDegree1Degree3Neighbours1PartitionSize[i] = 0;
            canonicalDegree1Degree3Neighbours2PartitionCount[i] = 0;
            canonicalDegree1Degree3Neighbours2PartitionSize[i] = 0;
        }
        canonicalDegree1Degree3Neighbours0PartitionCount = 0;
        canonicalDegree1Degree3Neighbours0PartitionSize = 0;
    #endif

    #ifdef _PROFILING_OPERATION1_SAME_COLOUR
    for(i=0;i<MAXN;i++){
        degree1Operation1Degree1Counts[i]=0;
    }
    #endif

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
        fprintf(stderr, "Generated %llu graph%s with only loops (and at least one loop).\n",
                graphsWithOnlyLoopsCount, graphsWithOnlyLoopsCount==1 ? (char *)"" : (char *)"s");
        fprintf(stderr, "Generated %llu graph%s with only semi-edges (and at least one semi-edge).\n",
                graphsWithOnlySemiEdgesCount, graphsWithOnlySemiEdgesCount==1 ? (char *)"" : (char *)"s");
        fprintf(stderr, "Generated %llu graph%s with only multi-edges (and at least one multi-edge).\n",
                graphsWithOnlyMultiEdgesCount, graphsWithOnlyMultiEdgesCount==1 ? (char *)"" : (char *)"s");
    }
    fprintf(stderr, "\nGenerated %llu pregraph primitive%s.\n", primitivesCount, primitivesCount==1 ? (char *)"" : (char *)"s");

    fprintf(stderr, "\nDegree 1 operations maximum recursion depth: %d.\n", degree1OperationsDepthMaximum);
    fprintf(stderr, "Degree 2 operations maximum recursion depth: %d.\n", degree2OperationsDepthMaximum);

    #ifdef _PROFILING_DEG1

    fprintf(stderr, "\n");
    fprintf(stderr, "Degree 1 operations\n");
    fprintf(stderr, "===================\n");
    fprintf(stderr, "Calls to canonicalDegree1 : %llu\n", canonicalDegree1Calls);
    fprintf(stderr, "Found canonical because only one vertex of degree 1  : %llu\n", canonicalDegree1BecauseOnlyOneVertexOfDegree1);
    fprintf(stderr, "   - with a trivial group : %llu\n", canonicalDegree1TrivialRemainsTrivial);
    fprintf(stderr, "   - group may be copied  : %llu\n", canonicalDegree1BridgeFixed);
    fprintf(stderr, "Found not canonical because not minimum colour       : %llu\n", canonicalDegree1NotBecauseNotSmallestColour);
    fprintf(stderr, "Found canonical because only one with minimum colour : %llu\n", canonicalDegree1BecauseOnlyOneMinimumColour);
    fprintf(stderr, "Found canonical after nauty was called               : %llu\n", canonicalDegree1YesWithNauty);
    fprintf(stderr, "Found not canonical after nauty was called           : %llu\n", canonicalDegree1NoWithNauty);
    fprintf(stderr, "\n");
    fprintf(stderr, "┌───────────────────────┬────────────┬────────────┬────────────┐\n");
    fprintf(stderr, "│ DEGREE 1              │   Total    │  Accepted  │  Rejected  │\n");
    fprintf(stderr, "├───────────────────────┼────────────┼────────────┼────────────┤\n");
    fprintf(stderr, "│ Decided without nauty │ %10llu │ %10llu │ %10llu │\n",
            canonicalDegree1BecauseOnlyOneVertexOfDegree1 + canonicalDegree1BecauseOnlyOneMinimumColour + canonicalDegree1NotBecauseNotSmallestColour,
            canonicalDegree1BecauseOnlyOneVertexOfDegree1 + canonicalDegree1BecauseOnlyOneMinimumColour,
            canonicalDegree1NotBecauseNotSmallestColour);
    fprintf(stderr, "│ Decided with nauty    │ %10llu │ %10llu │ %10llu │\n",
            canonicalDegree1YesWithNauty + canonicalDegree1NoWithNauty,
            canonicalDegree1YesWithNauty,
            canonicalDegree1NoWithNauty);
    fprintf(stderr, "├───────────────────────┼────────────┼────────────┼────────────┤\n");
    fprintf(stderr, "│ Operation 1           │ %10llu │ %10llu │ %10llu │\n",
            degree1Operation1Total,
            degree1Operation1Canonical,
            degree1Operation1Total - degree1Operation1Canonical);
    fprintf(stderr, "│ Operation 2           │ %10llu │ %10llu │ %10llu │\n",
            degree1Operation2Total,
            degree1Operation2Canonical,
            degree1Operation2Total - degree1Operation2Canonical);
    fprintf(stderr, "└───────────────────────┴────────────┴────────────┴────────────┘\n");

    //determine bounding box for colour frequency table
    int i, j;
    int minVertex, maxVertex, minColour, maxColour;
    minVertex = MAXN;
    maxVertex = 0;
    minColour = canonicalDegree1PossibleColoursCount;
    maxColour = 0;
    for(i = 0; i < canonicalDegree1PossibleColoursCount; i++){
        for(j = 0; j < MAXN; j++){
            if(canonicalDegree1MinimumColourFrequency[i*MAXN+j]>0){
                if(minColour > i){
                    minColour = i;
                }
                if(maxColour < i){
                    maxColour = i;
                }
                if(minVertex > j){
                    minVertex = j;
                }
                if(maxVertex < j){
                    maxVertex = j;
                }
            }
        }
    }

    fprintf(stderr, "\nThe following table shows the number of times the first colour was minimal, but still nauty was needed\n");
    fprintf(stderr, "because there were multiple vertices that shared that colour. The rows correspond to the colours and the\n");
    fprintf(stderr, "columns correspond with the number of vertices that shared that colour.\n");
    fprintf(stderr, "┌───────────────────┬");
    for(j = minVertex; j < maxVertex; j++){
        fprintf(stderr, "────────────┬");
    }
    fprintf(stderr, "────────────┐\n");
    fprintf(stderr, "│ Colour \\ vertices │");
    for(j = minVertex; j <= maxVertex; j++){
        fprintf(stderr, " %10d │", j);
    }
    fprintf(stderr, "\n");

    for(i = minColour; i <= maxColour; i++){
        fprintf(stderr, "├───────────────────┼");
        for(j = minVertex; j < maxVertex; j++){
            fprintf(stderr, "────────────┼");
        }
        fprintf(stderr, "────────────┤\n");
        fprintf(stderr, "│ %17d │", i);
        for(j = minVertex; j <= maxVertex; j++){
            fprintf(stderr, " %10llu │", canonicalDegree1MinimumColourFrequency[i*MAXN+j]);
        }
        fprintf(stderr, "\n");
    }

    fprintf(stderr, "└───────────────────┴");
    for(j = minVertex; j < maxVertex; j++){
        fprintf(stderr, "────────────┴");
    }
    fprintf(stderr, "────────────┘\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Size of partitions provided to nauty\n");
    for(j = 0; j < canonicalDegree1PossibleColoursCount + 2; j++){
        if(canonicalDegree1PartitionCountFrequency[j])
            fprintf(stderr, "%10d time%s a partition of size %d\n",
                canonicalDegree1PartitionCountFrequency[j],
                canonicalDegree1PartitionCountFrequency[j]==1 ? (char *)" " : (char *)"s",
                j);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "Size of partitions for degree 3 vertices\n");
    for(j = 0; j <= MAXN; j++){
        if(canonicalDegree1Degree3PartitionCount[j])
            fprintf(stderr, "%10llu time%s a partition of size %d\n",
                canonicalDegree1Degree3PartitionCount[j],
                canonicalDegree1Degree3PartitionCount[j]==1 ? (char *)" " : (char *)"s",
                j);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "Average size of partition classes of partitions provided to nauty\n");
    fprintf(stderr, "   - Degree 1 vertices\n");
    for(j = 0; j < canonicalDegree1PossibleColoursCount + 2; j++){
        if(canonicalDegree1Degree1PartitionCount[j]){
            fprintf(stderr, "         colour %3d: %f\n", j, canonicalDegree1Degree1PartitionSize[j]*1.0/canonicalDegree1Degree1PartitionCount[j]);
        }
    }
    fprintf(stderr, "   - Degree 3 vertices\n");
    if(canonicalDegree1Degree3Neighbours0PartitionCount){
        fprintf(stderr, "      * with 0 degree 1 neighbours\n");
        fprintf(stderr, "         colour   0: %f\n", canonicalDegree1Degree3Neighbours0PartitionSize*1.0/canonicalDegree1Degree3Neighbours0PartitionCount);
    }
    fprintf(stderr, "      * with 1 degree 1 neighbour\n");
    for(j = 0; j < canonicalDegree1PossibleColoursCount + 2; j++){
        if(canonicalDegree1Degree3Neighbours1PartitionCount[j]){
            fprintf(stderr, "        colour %3d: %f\n", j, canonicalDegree1Degree3Neighbours1PartitionSize[j]*1.0/canonicalDegree1Degree3Neighbours1PartitionCount[j]);
        }
    }
    fprintf(stderr, "      * with 2 degree 1 neighbours\n");
    for(j = 0; j < canonicalDegree1PossibleColoursCount + 2; j++){
        if(canonicalDegree1Degree3Neighbours2PartitionCount[j]){
            fprintf(stderr, "         colour %3d: %f\n", j, canonicalDegree1Degree3Neighbours2PartitionSize[j]*1.0/canonicalDegree1Degree3Neighbours2PartitionCount[j]);
        }
    }
    #endif

    #ifdef _PROFILING_OPERATION1_SAME_COLOUR
    fprintf(stderr, "\n\nStatistics about the first degree 1 operation.\n");
    fprintf(stderr, "Too few vertices to apply operation        : %llu\n\n", degree1Operation1TooFewVertices);
    fprintf(stderr, "Number of times the operation was tried    : %llu\n", degree1Operation1OperationTried);
    fprintf(stderr, "             - with different colour       : %llu\n", degree1Operation1HaveDifferentColour);
    fprintf(stderr, "             - with same colour            : %llu   (%.2f)\n\n", degree1Operation1HaveSameColour, (double)degree1Operation1HaveSameColour/degree1Operation1OperationTried);
    fprintf(stderr, "Number of times the colours weren't solved : %llu\n", degree1Operation1HaveSameColour - degree1Operation1SameColourSolved);
    fprintf(stderr, "Number of times the colours were solved    : %llu\n", degree1Operation1SameColourSolved);
    fprintf(stderr, "             - with paths                  : %llu   (%.2f)\n", degree1Operation1SolvedWithPaths, (double)degree1Operation1SolvedWithPaths/degree1Operation1SameColourSolved);
    fprintf(stderr, "             - without paths               : %llu\n\n", degree1Operation1NotSolvedWithPaths);

    fprintf(stderr, "Overview of the number of vertices of degree 1 that a graph had when this operation was called.\n");
    int k;
    for(k=2;k<vertexCount;k++){
        fprintf(stderr, "%10llu time%s called with %d vertices of degree 1\n", degree1Operation1Degree1Counts[k], degree1Operation1Degree1Counts[k]==1 ? (char *)"" : (char *)"s", k);
    }
    #endif

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
    #ifdef _DEBUG_MODIFY_GENERATION
    char *disabled;
    #endif
    char *moduloString;
    boolean fromFile = FALSE;
    char *inputFileName = NULL;
    FILE *inputFile = stdin;

    while ((c = getopt(argc, argv, "LSMPXCBf:F:o:D:d:Ihim:4ql")) != -1) {
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
#ifdef _DEBUG_MODIFY_GENERATION
            case 'X': //(defaults to FALSE)
                noRejections = TRUE;
                break;
#endif
            case 'C': //(defaults to FALSE)
                onlyColourable = TRUE;
                break;
            case 'B': //(defaults to FALSE)
                onlyBipartite = TRUE;
                break;
            case '4': //(defaults to FALSE)
                onlyC4Coverable = TRUE;
                break;
            case 'q': //(defaults to FALSE)
                onlyAdmissable = TRUE;
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
#ifdef _DEBUG_MODIFY_GENERATION
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
#endif
            case 'I':
                fromFile = TRUE;
                break;
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case 'l':
                printLimits(name);
                return EXIT_SUCCESS;
            case 'i':
                logStatistics = TRUE;
                break;
            case 'm':
                //modulo
                moduloEnabled = TRUE;
                moduloString = optarg;
                moduloRest = atoi(moduloString);
                moduloString = strchr(moduloString, ':');
                if(moduloString==NULL){
                    fprintf(stderr, "Illegal format for modulo.\n");
                    usage(name);
                    return EXIT_FAILURE;
                }
                moduloMod = atoi(moduloString+1);
                if (moduloRest >= moduloMod) {
                    fprintf(stderr, "Illegal format for modulo: rest must be smaller than mod.\n");
                    usage(name);
                    return EXIT_FAILURE;
                }
                if (moduloRest < 0) {
                    fprintf(stderr, "Illegal format for modulo: rest must be positive.\n");
                    usage(name);
                    return EXIT_FAILURE;
                }
                /**/
                moduloString = strchr(moduloString+1, ':');
                if(moduloString!=NULL){
                    splitDepth = atoi(moduloString+1);
                }
                /**/
                break;
            default:
                fprintf(stderr, "Illegal option %c.\n", c);
                usage(name);
                return EXIT_FAILURE;
        }
    }

    if(onlyC4Coverable) onlyAdmissable = FALSE;

    // check the non-option arguments
    if (argc - optind != 1) {
        usage(name);
        return EXIT_FAILURE;
    }

    //parse the order
    vertexCount = strtol(argv[optind], NULL, 10);
    DEBUGDUMP(vertexCount, "%d")

    if(vertexCount > MAXN){
        fprintf(stderr, "%s needs to be recompiled with a larger value for MAXN to support\ngraphs of order %d and with the given option.\n", name, vertexCount);
        fprintf(stderr, "MAXN needs to be at least %d for these parameters.\n", vertexCount);
        return EXIT_FAILURE;
    } else if (allowSemiEdges && 2*vertexCount + 2 > MAXN) {
        fprintf(stderr, "%s needs to be recompiled with a larger value for MAXN to support\ngraphs of order %d and with the given option.\n", name, vertexCount);
        fprintf(stderr, "MAXN needs to be at least %d for these parameters.\n", 2*vertexCount + 2);
        return EXIT_FAILURE;
    } else if (allowMultiEdges && vertexCount + vertexCount/2 > MAXN) {
        fprintf(stderr, "%s needs to be recompiled with a larger value for MAXN to support\ngraphs of order %d and with the given option.\n", name, vertexCount);
        fprintf(stderr, "MAXN needs to be at least %d for these parameters.\n", vertexCount + vertexCount/2);
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

    fprintf(stderr, "Generating %s%s%s%s%s with %d %s%s%s%s%s%s.\n",
            onlyPrimitives ? (char *)"pregraph primitives of " : (char *)"" ,
            onlyColourable ? (char *)"3-edge-colourable " : (char *)"",
            onlyColourable && onlyBipartite ? (char *)", " : (char *)"",
            onlyBipartite ? (char *)"bipartite " : (char *)"",
            allowMultiEdges ? (char *)"multigraphs" : (char *)"simple graphs",
            vertexCount, vertexCount==1 ? (char *)"vertex" : (char *)"vertices",
            allowLoops && allowSemiEdges ? (char *)", " : (allowLoops ? (char *)" and " : (char *)""),
            allowLoops ? (char *)"loops" : (char *)"", allowSemiEdges ? (char *)" and semi-edges" : (char *)"",
            onlyAdmissable ? (char *)" and filtering graphs that have a 2-factor where each component is a quotient of a 4-cycle" : (char *)"",
            onlyC4Coverable ? (char *)" and filtering graphs that have a 2-factor where each component is a 4-cycle" : (char *)"");
    if(moduloEnabled){
        fprintf(stderr, "Only generating part %d of %d (Splitting at depth %d).\n", moduloRest+1, moduloMod, splitDepth);
    }

    if(onlyAdmissable || onlyC4Coverable){
        onlyColourable=allowMultiEdges;
    }
    boolean onlySimpleGraphs = !allowSemiEdges && !allowLoops && !allowMultiEdges;
    
    struct tms TMS;
    unsigned int oldtime = 0;

    if(!allowSemiEdges && vertexCount%2==1){
        structureCount = primitivesCount = 0;
    } else if(onlyColourable && allowLoops){
        structureCount = primitivesCount = 0;
    } else if(onlyBipartite && allowLoops){
        structureCount = primitivesCount = 0;
    } else if(onlyAdmissable && allowLoops){
        structureCount = primitivesCount = 0;
    } else if(onlyC4Coverable && allowLoops){
        structureCount = primitivesCount = 0;
    } else if(onlyC4Coverable && vertexCount%4!=0){
        structureCount = primitivesCount = 0;
    } else if((onlySimpleGraphs && onlyAdmissable) && vertexCount%4!=0){
        structureCount = primitivesCount = 0;
    } else if(onlySimpleGraphs && vertexCount < 4){
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

    times(&TMS);
    unsigned int savetime = oldtime + (unsigned int) TMS.tms_utime;
    fprintf(stderr, "CPU time: %.1f seconds.\n", (double) savetime / time_factor);

    return EXIT_SUCCESS;
}

