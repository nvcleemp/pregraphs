/* bipartite.c
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
 * File:   bipartite.c
 * Author: nvcleemp
 *
 * Created on November 3, 2010, 13:54
 *
 * Checks whether the given pregraph is bipartite.
 */
//#define _DEBUG

#include "bipartite.h"

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
        pregraph->adjMatrix[i]=0;
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
		pregraph->adjMatrix[i] = pregraph->adjMatrix[i] | (1<<n);
		pregraph->degree[i]++;
            } else {
                pregraph->adjList[i][pregraph->degree[i]] = number-1;
                pregraph->adjList[number-1][pregraph->degree[number-1]] = i;
		pregraph->adjMatrix[i] = pregraph->adjMatrix[i] | (1<<(number-1));
		pregraph->adjMatrix[number-1] = pregraph->adjMatrix[number-1] | (1<<i);
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

void printPregraph(FILE *f, PREGRAPH *pregraph){
    int i;
    for(i=0; i<pregraph->order;i++){
        fprintf(f, "%d) %d %d %d\n", i+1, pregraph->adjList[i][0]+1, pregraph->adjList[i][1]+1, pregraph->adjList[i][2]+1);
    }
}

//======================================================================================================
//======================================================================================================

boolean isBipartite(PREGRAPH *pregraph){

    int i, queueHead, queueTail;
    int queue[MAXN];
    boolean queued[MAXN];
    int colour[MAXN];
    for(i=0; i<pregraph->order; i++) {
        queued[i] = FALSE;
        colour[i] = -1;
    }

    queueHead = 0;
    queueTail = 1;
    queue[0] = 0;
    queued[0] = TRUE;
    
    while(queueTail > queueHead){
        int currentVertex = queue[queueHead++];
        int numberOfColouredNeighbours = 0, sumOfColours = 0;
        DEBUGDUMP(currentVertex, "%d")
        DEBUGDUMP(pregraph->degree[currentVertex], "%d")
        
        //put all neighbours in queue
        for(i=0; i<pregraph->degree[currentVertex]; i++){
            DEBUGDUMP(i, "%d")
            DEBUGDUMP(pregraph->adjList[currentVertex][i], "%d")
            DEBUGDUMP(colour[pregraph->adjList[currentVertex][i]], "%d")
            if(pregraph->adjList[currentVertex][i]!=pregraph->order){
                if(!queued[pregraph->adjList[currentVertex][i]]){
                    queued[pregraph->adjList[currentVertex][i]]=TRUE;
                    queue[queueTail++] = pregraph->adjList[currentVertex][i];
                } else if (colour[pregraph->adjList[currentVertex][i]]>-1){
                    numberOfColouredNeighbours++;
                    sumOfColours+=colour[pregraph->adjList[currentVertex][i]];
                }
            }
        }
        
        if(sumOfColours==0){
            //either all coloured neighbours have colour 0 or no coloured neighbours
            colour[currentVertex] = 1;
        } else if(sumOfColours==numberOfColouredNeighbours){
            //all coloured neighbours have colour 1
            colour[currentVertex] = 0;
        } else {
            DEBUGDUMP(sumOfColours, "%d")
            DEBUGDUMP(numberOfColouredNeighbours, "%d")
            DEBUGARRAYDUMP(queue, 6, "%d")
            DEBUGARRAYDUMP(queued, 6, "%d")
            DEBUGARRAYDUMP(colour, 6, "%d")
            //conflicting colours
            return FALSE;
        }
    }

    return TRUE;

}

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

    while ((c = getopt(argc, argv, "hcf:")) != -1) {
        switch (c) {
            case 'c':
                onlyCount = TRUE;
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
    unsigned long bipartiteGraps = 0;

    while(readPregraphCode(stdin, &pregraph, &endian, count)==(char)1){
        count++;
        //DEBUGDUMP(count, "%ld")
        #ifdef _DEBUG
	fprintf(stderr, "Graph %2d\n========\n", count);
        printPregraph(stderr, pregraph);
        #endif
	if(isBipartite(&pregraph)){
            #ifdef _DEBUG
	    fprintf(stderr, "Bipartite\n\n");
            #endif
	    bipartiteGraps++;
            if(!onlyCount)
                writePregraphCode(stdout, &pregraph, LITTLE_ENDIAN, bipartiteGraps);
	}
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
    fprintf(outputFile, "of which %ld graphs are bipartite.\n", bipartiteGraps);
    return EXIT_SUCCESS;
}

