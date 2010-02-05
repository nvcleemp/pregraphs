/*
 * File:   3regpregraphtranslator.c
 * Author: nvcleemp
 *
 * Created on February 4, 2010, 9:52 AM
 */
//#define _DEBUG

//TODO: currently can't handle fat K2

#include "3regpregraphtranslator.h"

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

char readPregraphCodeNoHeader(FILE *f, PREGRAPH *pregraph, int endian) {
    int i, j, n, dummyVertex;
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
    pregraph->ppgraph->order=n; //will be increased when we find semi-edges
    for(i = 0; i < MAXN; i++){
        pregraph->semiEdgeVertices[i]=FALSE;
        pregraph->ppgraph->degree[i]=0;
    }
    pregraph->ppgraph->degree1Count=0;
    pregraph->ppgraph->multiEdgeCount=0;

    i = 1;
    dummyVertex = n;
    while (i <= n) {
        if (read_old_or_new(f, signum == 0, endian, &number) == 2) {
            return (2);
        }
        //fprintf(stderr, "%d ", number);
        DEBUGDUMP(i, "%d")
        DEBUGDUMP(number, "%d")
        if (number != 0) {
            if(number == i){
                DEBUGMSG("loop")
                //we found a loop
                pregraph->ppgraph->degree1Count++;
            } else if (number == n+1){
                DEBUGMSG("semi-edge")
                //we found a semi-edge
                pregraph->ppgraph->degree1Count++;
                pregraph->ppgraph->order++;
                pregraph->ppgraph->adjList[(i-1)*3 + pregraph->ppgraph->degree[i-1]] = dummyVertex;
                pregraph->ppgraph->degree[i-1]++;
                pregraph->ppgraph->adjList[dummyVertex*3 + 0] = i-1;
                pregraph->ppgraph->degree[dummyVertex]=1;
                pregraph->semiEdgeVertices[dummyVertex]=TRUE;
                dummyVertex++;
            } else {
                //first check to see if we have a multi-edge
                j=0;
                while(j<pregraph->ppgraph->degree[i-1] && pregraph->ppgraph->adjList[(i-1)*3 + j] != (number-1)) j++;
                DEBUGDUMP(j, "%d")
                DEBUGDUMP(pregraph->ppgraph->degree[i-1], "%d")
                if(j != pregraph->ppgraph->degree[i-1]){
                    DEBUGMSG("multi-edge")
                    //we found a multi-edge (do not increase the degree)
                    pregraph->ppgraph->multiEdgeCount++;
                    pregraph->ppgraph->multiedge[i-1]=number-1;
                    pregraph->ppgraph->multiedge[number-1]=i-1;
                } else {
                    DEBUGMSG("simple edge")
                    pregraph->ppgraph->adjList[(i-1)*3 + pregraph->ppgraph->degree[i-1]] = number-1;
                    pregraph->ppgraph->adjList[(number-1)*3 + pregraph->ppgraph->degree[number-1]] = i-1;
                    pregraph->ppgraph->degree[i-1]++;
                    pregraph->ppgraph->degree[number-1]++;
                }
            }
        } else {
            DEBUGMSG("next vertex")
            i++;
            //fprintf(stderr, "%d \n", number);
        }
    }
    pregraph->ppgraph->multiEdgeCount/=2; //counted twice
    return (1);
}

char readPregraphCode(FILE *f, PREGRAPH *pregraph, int *endian, unsigned long count) {
    if (count == 0) {
        if (read_endian(f, endian) == 2) {
            return (2);
        }
    }
    return (readPregraphCodeNoHeader(f, pregraph, *endian));
}

char writePregraphTable(FILE *f, PREGRAPH *pregraph, unsigned long count) {
    fprintf(f, "==============================\n");
    fprintf(f, "|  Graph number: %10ld  |\n", count);
    fprintf(f, "|  Number of vertices: %4d  |\n", pregraph->order);
    fprintf(f, "==============================\n");

    unsigned short i, j;
    PRIMPREGRAPH *ppgraph = pregraph->ppgraph;
    int primPregraph2Pregraph[ppgraph->order];
    j = 1; //vertices are labeled starting from 1
    for (i = 0; i < ppgraph->order; i++) {
        if (pregraph->semiEdgeVertices[i]) {
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

char writePrimPregraphTable(FILE *f, PREGRAPH *pregraph, unsigned long count) {
    fprintf(f, "==============================\n");
    fprintf(f, "|  Graph number: %10ld  |\n", count);
    fprintf(f, "|  pregraph order: %8d  |\n", pregraph->order);
    fprintf(f, "|  primitive order: %7d  |\n", pregraph->ppgraph->order);
    fprintf(f, "==============================\n");

    unsigned short i, j;
    PRIMPREGRAPH *ppgraph = pregraph->ppgraph;
    for (i = 0; i < ppgraph->order; i++) {
            fprintf(f, "|%4d ||", i);
            for (j = 0; j < ppgraph->degree[i]; j++) {
                fprintf(f, " %4d |", ppgraph->adjList[i * 3 + j]);
            }
            if (j == 2) {
                //add multi-edge
                fprintf(f, "(%4d)|", ppgraph->multiedge[i]);
                j++;
            }
            for (; j < 3; j++) {
                fprintf(f, "      |");
            }
            fprintf(f,"|\n");
    }
    fprintf(f, "==============================\n");
    fprintf(f,"\n");
    return (ferror(f) ? 2 : 1);
}

char writeMulticode(FILE *f, PREGRAPH *pregraph, boolean header, int endian){
    //first we translate the pregraph to a simple graph
    PRIMPREGRAPH *ppgraph = pregraph->ppgraph;
    int loopCount = ppgraph->degree1Count - (ppgraph->order - pregraph->order);
    int newOrder = ppgraph->order + loopCount*2 + ppgraph->multiEdgeCount*2;
    DEBUGDUMP(loopCount, "%d")
    DEBUGDUMP(ppgraph->multiEdgeCount, "%d")
    DEBUGDUMP(ppgraph->degree1Count, "%d")
    DEBUGDUMP(ppgraph->order, "%d")
    DEBUGDUMP(pregraph->order, "%d")
    int adjacency[newOrder*3], degree[newOrder];
    int i;
    for (i = 0; i < newOrder; i++) {
        degree[i]=0;
    }
    int offset = ppgraph->order;

    for (i = 0; i < ppgraph->order; i++){
        if(ppgraph->degree[i]==1){
            adjacency[i*3+0]=ppgraph->adjList[i*3+0];
            if(pregraph->semiEdgeVertices[i]){
                degree[i]=1;
            } else {
                degree[i]=3;
                adjacency[i*3+1]=offset;
                adjacency[i*3+2]=offset+1;
                adjacency[offset*3+0]=i;
                adjacency[offset*3+1]=offset+1;
                adjacency[(offset+1)*3+0]=i;
                adjacency[(offset+1)*3+1]=offset;
                degree[offset]=degree[offset+1]=2;
                offset+=2;
            }
        } else if(ppgraph->degree[i]==2){
            degree[i]=3;
            if(ppgraph->adjList[i*3+0]==ppgraph->multiedge[i]){
                adjacency[i*3+0]=ppgraph->adjList[i*3+1];
            } else {
                adjacency[i*3+0]=ppgraph->adjList[i*3+0];
            }
            if(ppgraph->multiedge[i]>i){
                adjacency[i*3+1]=offset;
                adjacency[i*3+2]=offset+1;
                adjacency[(ppgraph->multiedge[i])*3+1]=offset;
                adjacency[(ppgraph->multiedge[i])*3+2]=offset+1;
                adjacency[offset*3+0]=i;
                adjacency[offset*3+1]=ppgraph->multiedge[i];
                adjacency[(offset+1)*3+0]=i;
                adjacency[(offset+1)*3+1]=ppgraph->multiedge[i];
                degree[offset]=degree[offset+1]=2;
                offset+=2;
            }
        } else {
            degree[i]=3;
            adjacency[i*3+0]=ppgraph->adjList[i*3+0];
            adjacency[i*3+1]=ppgraph->adjList[i*3+1];
            adjacency[i*3+2]=ppgraph->adjList[i*3+2];
        }
    }

    //next we write out the simple graph
    unsigned short j;
    if (header) { //if first graph
        fprintf(f, ">>multi_code %s<<", (endian == LITTLE_ENDIAN ? "le" : "be"));
    }
    if (debugOutput) {
        fprintf(f, "%d ", newOrder);
    } else if (newOrder <= UCHAR_MAX) {
        fprintf(f, "%c", (unsigned char) newOrder);
    } else {
        fprintf(f, "%c", 0);
        /* big graph */
        if (write_2byte_number(f, (unsigned short) newOrder, endian) == 2) {
            return (2);
        }
    }
    for (i = 0; i < newOrder-1; i++) {
        for (j = 0; j < degree[i]; j++) {
            if(adjacency[i * 3 + j]>i){
                if (debugOutput) {
                    fprintf(f, "%d ", adjacency[i * 3 + j] + 1);
                } else if (newOrder <= UCHAR_MAX) {
                    fprintf(f, "%c", (unsigned char) adjacency[i * 3 + j] + 1);
                } else if (write_2byte_number(f, adjacency[i * 3 + j] + 1, endian) == 2) {
                    return (2);
                }
            }
        }
        //closing 0
        if (debugOutput) {
            fprintf(f, "0 ");
        } else if (newOrder <= UCHAR_MAX) {
            fprintf(f, "%c", 0);
        } else {
            if (write_2byte_number(f, 0, endian) == 2) {
                return (2);
            }
        }
    }
    if (debugOutput) {
        fprintf(f, "\n");
    }
    return (ferror(f) ? 2 : 1);
}

/*
print a usage message. name is the name of the current program.
 */
void usage(char *name) {
    fprintf(stderr, "Usage: %s [t/m/h]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

/*
print a help message. name is the name of the current program.
 */
void help(char *name) {
    fprintf(stderr, "The program %s .\n", name);
    fprintf(stderr, "Usage: %s [t/m/h] \n\n", name);
    fprintf(stderr, "Valid options:\n");
    fprintf(stderr, "  -h          : Print this help and return.\n");
    fprintf(stderr, "  -t          : .\n");
    fprintf(stderr, "  -m          : .\n");
}

/*
 *
 */
int main(int argc, char** argv) {

    /*=========== commandline parsing ===========*/

    boolean table = FALSE;
    boolean multicode = FALSE;

    int c;
    char *name = argv[0];
    int endian = defaultEndian;

    while ((c = getopt(argc, argv, "htmd")) != -1) {
        switch (c) {
            case 't':
                table = TRUE;
                break;
            case 'm':
                multicode = TRUE;
                break;
            case 'd':
                debugOutput = TRUE;
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
    if (argc - optind != 0) {
        usage(name);
        return EXIT_FAILURE;
    }

    if(table == multicode) {
        usage(name);
        //use exactly one of t or m
        return EXIT_FAILURE;
    }

    /*=========== initialization ===========*/
    PREGRAPH pregraph;
    PRIMPREGRAPH ppgraph;

    pregraph.ppgraph = &ppgraph;

    unsigned long count = 0;

    /*
    outputFile = fopen("output.debug", "r");
    if(outputFile != NULL){
        fprintf(stderr, "File output.debug already exists: aborting!\n");
        return EXIT_FAILURE;
    }
    outputFile = fopen("output.debug", "w");
    */

    while(readPregraphCode(stdin, &pregraph, &endian, count)==(char)1){
        //fprintf(outputFile, "\n\n");
        count++;
        if(table)
            writePregraphTable(stdout, &pregraph, count);
        else if (multicode)
            writeMulticode(stdout, &pregraph, count==1, endian);
    }

    return EXIT_SUCCESS;
}

