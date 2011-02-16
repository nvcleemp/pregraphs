/* has3edgecolouring.c
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
 * File:   has3edgecolouring.c
 * Author: nvcleemp
 *
 * Created on August 10, 2010, 11:20 AM
 *
 * Checks whether the given pregraph has a 3-edge-colouring. The only accepted pregraphs are
 * multigraphs with semi-edges.
 */
//#define _DEBUG

#include "has3edgecolouring.h"

//=========================Methods copied from gconv.c=========================

/****************************READ_TO_END_OF_HEADER************************/

char read_to_end_of_header(FILE *f, unsigned char *c) {
    while (strncmp((char *) c, "<<", 2) != 0) {
        c[0] = c[1];
        if (fread(&c[1], sizeof (unsigned char), 1, f) == 0) {
            return (2);
        }
    }
    return (1);
}

/***************************READ_HEADER_ENTRY*****************************/
/*  This function finds out whether the string "s" is the next entry in
   the header. If it is, it returns 1, if it is not, it returns 0, if
   If an error occurrs, it returns 2.
   Whitespaces in the header are skipped.                               */

char read_header_entry(FILE *f, char *s) {
    unsigned int pos = 0;
    char c;
    while (s[pos]) { /* not end of string */
        if (fread(&c, sizeof (char), 1, f) == 0) {
            return (2);
        }
        if (c != ' ') { /* skip whitespace */
            if (s[pos++] != c) {
                if (ungetc(c, f) == EOF) {
                    return (2);
                } else {
                    return (0);
                }
                /* last character must be read again because it might be a part
                   of the header end "<<" */
            }
        }
    }
    return (1);
}

/**********************WRITE_2BYTE_NUMBER***************************/
/*  This procedure takes care of the endian                        */

char write_2byte_number(FILE *f, unsigned short n, int endian) {
    if (endian == BIG_ENDIAN) {
        fprintf(f, "%c%c", n / 256, n % 256);
    } else {
        fprintf(f, "%c%c", n % 256, n / 256);
    }
    return (ferror(f) ? 2 : 1);
}

/***********************READ_2BYTE_NUMBER***************************/
/*  This procedure takes care of the endian                        */

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

/************************READ_ENDIAN**************************************/
/*  Reads in the endian type (le/be), if existing, and skips to the end  */
/*  of the header                                                 	 */

char read_endian(FILE *f, int *endian) {
    unsigned char c[2];
    do {
        if (fread(&c[0], sizeof (unsigned char), 1, f) == 0) {
            return (2);
        }
    } while (isspace(c[0]));
    if (fread(&c[1], sizeof (unsigned char), 1, f) == 0) {
        return (2);
    }
    if (strncmp((char *) & c[0], "le", 2) == 0) {
        *endian = LITTLE_ENDIAN;
    } else if (strncmp((char *) & c[0], "be", 2) == 0) {
        *endian = BIG_ENDIAN;
    } else {
        *endian = defaultEndian;
    }
    if (read_to_end_of_header(f, &c[0]) == 2) {
        return (2);
    }
    return (1);
}

//======================End Methods copied from gconv.c=========================

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
    //fprintf(outputFile, "%d ", *number);
    return (1);
}

//modified version: only multi-edges and semi-edges (no loops)
char readPregraphCodeNoHeader(FILE *f, PREGRAPH *pregraph, int endian) {
    int i, n;
    unsigned short signum, number;
    if (read_old_or_new(f, FALSE, endian, &signum) == 2) {
        return (feof(f) ? 0 : 2);
    }
    //if the code starts with a zero, all the entries are two bytes
    //else the number we just read was the order of the graph
    if (signum == 0) {
        if (read_old_or_new(f, TRUE, endian, &number) == 2) {
            return (2);
        }
    } else {
        number = signum;
    }

    if ((n = (int) number) > MAXN) {
        return (3);
    }

    //initialize the pregraph
    pregraph->order=n;
    for(i = 0; i < MAXN; i++){
	pregraph->degree[i]=0;
    }

    i = 0;
    while (i < n) {
        if (read_old_or_new(f, signum == 0, endian, &number) == 2) {
            return (2);
        }
        //DEBUGDUMP(i, "%d")
        //DEBUGDUMP(number-1, "%d")
        if (number != 0) {
            if (number == n+1){
                //DEBUGMSG("semi-edge")
                //we found a semi-edge
                pregraph->adjList[i][pregraph->degree[i]] = n;
		pregraph->degree[i]++;
            } else {
                pregraph->adjList[i][pregraph->degree[i]] = number-1;
                pregraph->adjList[number-1][pregraph->degree[number-1]] = i;
		pregraph->degree[i]++;
		pregraph->degree[number-1]++;
            }
        } else {
            //DEBUGMSG("=========next vertex=========")
            i++;
        }
    }
    DEBUGMSG("Graph read")
    return (1);
}

char readPregraphCode(FILE *f, PREGRAPH *pregraph, int *endian, unsigned long count) {
    if (count == 0) {
        if (read_endian(f, endian) == 2) {
            DEBUGMSG("Read error")
            return (2);
        }
    }
    return (readPregraphCodeNoHeader(f, pregraph, *endian));
}

int writePregraphCode(FILE *f, PREGRAPH *pregraph, int endian, unsigned long structureCount){
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
    int i, j;
    for(i=0; i<pregraph->order; i++){
        for(j=0; j<3; j++){
            if(pregraph->adjList[i][j]>=i){
                if (pregraph->order + 1 <= UCHAR_MAX) {
                    fprintf(f, "%c", (unsigned char) (pregraph->adjList[i][j]+1));
                } else {
                    if (write_2byte_number(f, (unsigned short) (pregraph->adjList[i][j]+1), endian) == 2) {
                        return (2);
                    }
                }
            }
        }
        if (pregraph->order + 1 <= UCHAR_MAX) {
            fprintf(f, "%c", (unsigned char) 0);
        } else {
            if (write_2byte_number(f, (unsigned short) 0, endian) == 2) {
                return (2);
            }
        }
    }
    return 1;
}

void writePregraphTable(FILE *f, PREGRAPH *pregraph){
    int i, j;
    for(i=0; i<pregraph->order; i++){
        fprintf(f, "%d) ", i+1);
        for(j=0; j<pregraph->degree[i]; j++){
            fprintf(f, "%d ", pregraph->adjList[i][j]+1);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
}

//======================================================================================================
//======================================================================================================

//------------------Start colouring methods--------------------------------
void printColours(FILE *f, PREGRAPH *pregraph){
    int i, j;
    for(i = 0; i<pregraph->order; i++){
        fprintf(f, "%d) ", i + 1);
        for(j = 0; j<pregraph->degree[i]; j++){
            if(ISMARKED_EDGES(i, j)) {
                fprintf(f, "%d ", colours[i][j]);
            } else {
                fprintf(f, "X ");
            }
        }
        fprintf(f, "\n");
    }
}

void initIsColourable(PREGRAPH *ppgraph) {
    RESETMARKS_EDGES;
    int i, j;
    for (i = 0; i < ppgraph->order; i++) {
        #ifdef _DEBUG
        //clean up mapping when we are debugging
        for(j = 0; j < ppgraph->order; j++) {
            neighbourToIndexMapping[i][j] = -1;
        }
        #endif
        coloursAroundVertex[i] = 0;
        for(j = 0; j < ppgraph->degree[i]; j++) {
            neighbourToIndexMapping[i][ppgraph->adjList[i][j]] = j;
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

void determineUncolouredVertex(int vertex, int *uncolouredVertex, int *missingColour, PREGRAPH *ppgraph) {
    DEBUGASSERT(coloursAroundVertex[vertex] == 2);

    int i;
    int sum_colours = 0;
    for(i = 0; i < ppgraph->degree[vertex]; i++) {
        if(!ISMARKED_EDGES(vertex, i)) {
            *uncolouredVertex = ppgraph->adjList[vertex][i];
        } else {
            sum_colours += colours[vertex][i];
        }
    }
    switch(sum_colours) {
        case 3:
            *missingColour = 3;
            break;
        case 4:
            *missingColour = 2;
            break;
        case 5:
            *missingColour = 1;
            break;
        default:
            fprintf(stderr, "Invalid sum_colours: %d\n", sum_colours);
            fprintf(stderr, "vertex %d\n", vertex+1);
            writePregraphTable(stderr, ppgraph);
            printColours(stderr, ppgraph);
            int k,l;
            fprintf(stderr, "    ");
            for(k=0; k<ppgraph->order; k++){
                fprintf(stderr, "%2d ", k+1);
            }
            fprintf(stderr, "\n");
            for(k=0; k<ppgraph->order; k++){
                fprintf(stderr, "%2d) ", k+1);
                for(l=0; l<ppgraph->order; l++){
                    if(neighbourToIndexMapping[k][l]>=0)
                        fprintf(stderr, "%2d ", neighbourToIndexMapping[k][l]+1);
                    else
                        fprintf(stderr, "   ");
                }
                fprintf(stderr, "\n");
            }
            exit(0);
    }

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
boolean propagateFixedColours(int currentVertex, int nonfree_labelled[][2], int *nonfree_labelled_size, PREGRAPH *ppgraph) {
    int uncolouredVertex, missingColour;
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
int tryExtendingColouring(int numberOfColouredEdges, int numberOfEdges, PREGRAPH *ppgraph) {
    if(numberOfColouredEdges != numberOfEdges) {
        int currentVertex;

        for(currentVertex = 1; currentVertex < ppgraph->order; currentVertex++) {
            if(coloursAroundVertex[currentVertex] == 1 && ppgraph->degree[currentVertex] != 1) {
                break;
            }
            DEBUGASSERT(coloursAroundVertex[currentVertex] != 2)
        }
        DEBUGASSERT(currentVertex < ppgraph->order);

        int usedColour; //the colour already used at this vertex
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
        availableVertices[0] = ppgraph->adjList[currentVertex][indexAvailableVertex0];
        availableVertices[1] = ppgraph->adjList[currentVertex][indexAvailableVertex1];

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
            boolean abort = 0;
            for(j = 0; j < 2;j++) {
                if(!propagateFixedColours(availableVertices[j], nonfree_labelled, &nonfree_labelled_size, ppgraph)) {
                    DEBUGASSERT(nonfree_labelled_size <= numberOfEdges - numberOfColouredEdges)
                    abort = TRUE;
                    break;
                }
                DEBUGASSERT(nonfree_labelled_size <= numberOfEdges - numberOfColouredEdges)
            }
            //fprintf(stderr, "nonfree_labelled_size: %d and abort == %d \n", nonfree_labelled_size, abort);
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
boolean isColourableGraph(PREGRAPH *ppgraph) {
    int degree1Count=0;
    int j;
    for(j=0; j<ppgraph->order; j++){
        if(ppgraph->degree[j]==1)
            degree1Count++;
    }
    if(degree1Count==1) return FALSE;

    DEBUGMSG("Initializing isColourable")
    initIsColourable(ppgraph);
    DEBUGMSG("IsColourable initialized")

    int currentVertex = 0;
    int i, neighbour, currentIndex;
    for(i = 0; i < ppgraph->degree[currentVertex]; i++) {
        colours[currentVertex][i] = i + 1;
        neighbour = ppgraph->adjList[currentVertex][i];
        currentIndex = neighbourToIndexMapping[neighbour][currentVertex];
        colours[neighbour][currentIndex] = i + 1;

        MARK_EDGES(currentVertex, i);
        MARK_EDGES(neighbour, currentIndex);
        coloursAroundVertex[neighbour] = 1;
    }
    coloursAroundVertex[currentVertex] = ppgraph->degree[currentVertex];
    int numberOfEdges = (3*(ppgraph->order-degree1Count) + degree1Count)/2;
    return tryExtendingColouring(ppgraph->degree[currentVertex], numberOfEdges, ppgraph);
}

//basically just strips the semi-edges from the graph and replaces multi-edges with diamonds
void translateToSimpleGraph(PREGRAPH *pregraph, PREGRAPH *simplegraph){
    simplegraph->order = pregraph->order;
    int i, j, v;
    int newVertex = pregraph->order;
    //first we replace semi-edges with vertices of degree 1
    for(i=0; i<pregraph->order; i++){
        simplegraph->degree[i] = 3;
        for(j=0; j<3; j++){
            if(pregraph->adjList[i][j]!=pregraph->order){
                simplegraph->adjList[i][j] = pregraph->adjList[i][j];
            } else {
                simplegraph->adjList[i][j] = newVertex;
                simplegraph->adjList[newVertex][0] = i;
                simplegraph->degree[newVertex] = 1;
                simplegraph->order++;
                newVertex++;
                if(newVertex>MAXN) {fprintf(stderr, "Graph too large.\n"); exit(1);}
            }
        }
    }
    //next we insert an edge between multi-edges
    for(v=0; v<simplegraph->order; v++){
        if(simplegraph->degree[v]==3){
            if(simplegraph->adjList[v][0]==simplegraph->adjList[v][1]){
                int neighbour = simplegraph->adjList[v][0];
                simplegraph->adjList[v][0] = newVertex;
                simplegraph->adjList[v][1] = newVertex + 1;
                simplegraph->adjList[newVertex][0] = newVertex + 1;
                simplegraph->adjList[newVertex][1] = v;
                simplegraph->adjList[newVertex][2] = neighbour;
                simplegraph->degree[newVertex] = 3;
                simplegraph->adjList[newVertex+1][0] = newVertex;
                simplegraph->adjList[newVertex+1][1] = v;
                simplegraph->adjList[newVertex+1][2] = neighbour;
                simplegraph->degree[newVertex+1] = 3;
                simplegraph->order+=2;
                if(simplegraph->adjList[neighbour][0]!=v){
                    simplegraph->adjList[neighbour][1] = newVertex;
                    simplegraph->adjList[neighbour][2] = newVertex + 1;
                } else if (simplegraph->adjList[neighbour][1]!=v){
                    simplegraph->adjList[neighbour][0] = newVertex;
                    simplegraph->adjList[neighbour][2] = newVertex + 1;
                } else {
                    simplegraph->adjList[neighbour][0] = newVertex;
                    simplegraph->adjList[neighbour][1] = newVertex + 1;
                }
                newVertex+=2;
                if(newVertex>MAXN) {fprintf(stderr, "Graph too large.\n"); exit(1);}
            } else if(simplegraph->adjList[v][0]==simplegraph->adjList[v][2]){
                int neighbour = simplegraph->adjList[v][0];
                simplegraph->adjList[v][0] = newVertex;
                simplegraph->adjList[v][2] = newVertex + 1;
                simplegraph->adjList[newVertex][0] = newVertex + 1;
                simplegraph->adjList[newVertex][1] = v;
                simplegraph->adjList[newVertex][2] = neighbour;
                simplegraph->degree[newVertex] = 3;
                simplegraph->adjList[newVertex+1][0] = newVertex;
                simplegraph->adjList[newVertex+1][1] = v;
                simplegraph->adjList[newVertex+1][2] = neighbour;
                simplegraph->degree[newVertex+1] = 3;
                simplegraph->order+=2;
                if(simplegraph->adjList[neighbour][0]!=v){
                    simplegraph->adjList[neighbour][1] = newVertex;
                    simplegraph->adjList[neighbour][2] = newVertex + 1;
                } else if (simplegraph->adjList[neighbour][1]!=v){
                    simplegraph->adjList[neighbour][0] = newVertex;
                    simplegraph->adjList[neighbour][2] = newVertex + 1;
                } else {
                    simplegraph->adjList[neighbour][0] = newVertex;
                    simplegraph->adjList[neighbour][1] = newVertex + 1;
                }
                newVertex+=2;
                if(newVertex>MAXN) {fprintf(stderr, "Graph too large.\n"); exit(1);}
            } else if(simplegraph->adjList[v][1]==simplegraph->adjList[v][2]){
                int neighbour = simplegraph->adjList[v][1];
                simplegraph->adjList[v][1] = newVertex;
                simplegraph->adjList[v][2] = newVertex + 1;
                simplegraph->adjList[newVertex][0] = newVertex + 1;
                simplegraph->adjList[newVertex][1] = v;
                simplegraph->adjList[newVertex][2] = neighbour;
                simplegraph->degree[newVertex] = 3;
                simplegraph->adjList[newVertex+1][0] = newVertex;
                simplegraph->adjList[newVertex+1][1] = v;
                simplegraph->adjList[newVertex+1][2] = neighbour;
                simplegraph->degree[newVertex+1] = 3;
                simplegraph->order+=2;
                if(simplegraph->adjList[neighbour][0]!=v){
                    simplegraph->adjList[neighbour][1] = newVertex;
                    simplegraph->adjList[neighbour][2] = newVertex + 1;
                } else if (simplegraph->adjList[neighbour][1]!=v){
                    simplegraph->adjList[neighbour][0] = newVertex;
                    simplegraph->adjList[neighbour][2] = newVertex + 1;
                } else {
                    simplegraph->adjList[neighbour][0] = newVertex;
                    simplegraph->adjList[neighbour][1] = newVertex + 1;
                }
                newVertex+=2;
                if(newVertex>MAXN) {fprintf(stderr, "Graph too large.\n"); exit(1);}
            }
        }
    }
}

boolean isColourable(PREGRAPH *pregraph){
    PREGRAPH simplegraph;

    translateToSimpleGraph(pregraph, &simplegraph);
    DEBUGMSG("Graph translated")

    return isColourableGraph(&simplegraph);
    
}

//-------------------End colouring methods---------------------------------

//======================================================================================================
//======================================================================================================

/*
print a usage message. name is the name of the current program.
 */
void usage(char *name) {
    fprintf(stderr, "Usage: %s\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

/*
print a help message. name is the name of the current program.
 */
void help(char *name) {
    fprintf(stderr, "The program %s .\n", name);
    fprintf(stderr, "Usage: %s [options] \n\n", name);
    fprintf(stderr, "Valid options:\n");
    fprintf(stderr, "  -c      : Only count the number of graphs.\n");
    fprintf(stderr, "  -f name : Prints statistics to file with geven name.\n");
    fprintf(stderr, "  -h      : Print this help and return.\n");
}

/*
 *
 */
int main(int argc, char** argv) {

    /*=========== commandline parsing ===========*/

    int c;
    char *name = argv[0];
    char *outputFileName = NULL;
    FILE *outputFile;
    int endian = defaultEndian;
    boolean onlyCount = FALSE;
    boolean negate = FALSE;

    while ((c = getopt(argc, argv, "hcf:n")) != -1) {
        switch (c) {
            case 'c':
                onlyCount = TRUE;
                break;
            case 'n':
                negate = TRUE;
                break;
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case 'f': //(defaults to stderr)
                outputFileName = optarg;
                break;
            default:
                usage(name);
                return EXIT_FAILURE;
        }
    }

    // check the non-option arguments
    if (argc - optind != 0) {
        usage(name);
        return EXIT_FAILURE;
    }

    /*=========== initialization ===========*/
    PREGRAPH pregraph;

    unsigned long count = 0;
    unsigned long allows3EdgeColouring = 0;
    unsigned long doesntAllow3EdgeColouring = 0;

    while(readPregraphCode(stdin, &pregraph, &endian, count)==(char)1){
        count++;
        DEBUGDUMP(count, "%ld")
        #ifdef _DEBUG
	fprintf(stderr, "Graph %2d\n========\n", count);
        #endif
        /**/

        if(negate){
            if(!isColourable(&pregraph)){
                #ifdef _DEBUG
                fprintf(stderr, "non-3-edge-colourable\n\n");
                #endif
                doesntAllow3EdgeColouring++;
                if(!onlyCount)
                    writePregraphCode(stdout, &pregraph, LITTLE_ENDIAN, doesntAllow3EdgeColouring);
            }
        } else {
            if(isColourable(&pregraph)){
                #ifdef _DEBUG
                fprintf(stderr, "3-edge-colourable\n\n");
                #endif
                allows3EdgeColouring++;
                if(!onlyCount)
                    writePregraphCode(stdout, &pregraph, LITTLE_ENDIAN, allows3EdgeColouring);
            }
        }
        /* */
    }

    if(outputFileName != NULL){
        outputFile = fopen(outputFileName, "a");
        if(outputFile == NULL){
            fprintf(stderr, "File %s couldn't be created: aborting!\n", outputFileName);
            return EXIT_FAILURE;
        }
    } else {
        outputFile = stderr;
    }

    fprintf(outputFile, "Read %ld graphs ", count);
    if(negate)
        fprintf(outputFile, "of which %ld graphs are non-3-edge-colourable.\n", doesntAllow3EdgeColouring);
    else
        fprintf(outputFile, "of which %ld graphs are 3-edge-colourable.\n", allows3EdgeColouring);
    return EXIT_SUCCESS;
}

