/*
 * File:   pgfilter.c
 * Author: nvcleemp
 *
 * Created on January 5, 2011, 9:52 AM
 *
 * Utility program to filter some 3-regular pregraphs out of a file.
 */
//#define _DEBUG

#include "pgfilter.h"
#include <string.h>
#include <ctype.h>

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
		pregraph->adjMatrix[i] =
                        pregraph->adjMatrix[i] | (1<<(number-1));
		pregraph->adjMatrix[number-1] =
                        pregraph->adjMatrix[number-1] | (1<<i);
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

char readPregraphCode(FILE *f, PREGRAPH *pregraph, int *endian,
        unsigned long count) {
    if (count == 0) {
        if (read_endian(f, endian) == 2) {
            DEBUGMSG("Read error")
            return (2);
        }
    }
    return (readPregraphCodeNoHeader(f, pregraph, *endian));
}

int writePregraphCode(FILE *f, PREGRAPH *pregraph, int endian,
        unsigned long structureCount){
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
                    if (write_2byte_number(f,
                            (unsigned short) (pregraph->adjList[i][j]+1), endian) == 2) {
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
        for(j=0; j<3; j++){
            fprintf(f, "%d ", pregraph->adjList[i][j]+1);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
}

void outputPregraph(FILE *f, PREGRAPH *pregraph, int endian,
        unsigned long structureCount, boolean table){
    if(table){
        writePregraphTable(f, pregraph);
    } else {
        writePregraphCode(f, pregraph, endian, structureCount);
    }
}

//=============================================================================
//=============================================================================

/*
print a usage message. name is the name of the current program.
 */
void usage(char *name) {
    fprintf(stderr, "Usage: %s -n number [options]\n", name);
    fprintf(stderr, "       %s -L lower bound [options]\n", name);
    fprintf(stderr, "       %s -U upper bound [options]\n", name);
    fprintf(stderr, "       %s -L lower bound -U upper bound [options]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

/*
print a help message. name is the name of the current program.
 */
void help(char *name) {
    fprintf(stderr, "The program %s .\n", name);
    fprintf(stderr, "Usage: %s -n number [options]\n", name);
    fprintf(stderr, "       %s -L lower bound [options]\n", name);
    fprintf(stderr, "       %s -U upper bound [options]\n", name);
    fprintf(stderr, "       %s -L lower bound -U upper bound [options]\n", name);
    fprintf(stderr, "Valid options:\n");
    fprintf(stderr, "  -h      : Print this help and return.\n");
    fprintf(stderr, "  -t      : Output human readable tables instead of pregraph code.\n");
    fprintf(stderr, "  -n      : Output graph number n.\n");
    fprintf(stderr, "  -L      : Output only graphs starting from this graph.\n");
    fprintf(stderr, "  -U      : Output no graphs beyond this graph.\n");
}

/*
 *
 */
int main(int argc, char** argv) {

    /*=========== commandline parsing ===========*/

    int c;
    char *name = argv[0];
    int endian = defaultEndian;
    int graph = -1;
    int lower = -1;
    int upper = -1;
    boolean table = FALSE;

    while ((c = getopt(argc, argv, "hn:L:U:t")) != -1) {
        switch (c) {
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case 't':
                table = TRUE;
                break;
            case 'n':
                if(lower!=-1 || upper!=-1){
                    usage(name);
                    return EXIT_FAILURE;
                }
                graph = strtol(optarg, NULL, 10);
                break;
            case 'L':
                if(graph!=-1){
                    usage(name);
                    return EXIT_FAILURE;
                }
                lower = strtol(optarg, NULL, 10);
                break;
            case 'U':
                if(graph!=-1){
                    usage(name);
                    return EXIT_FAILURE;
                }
                upper = strtol(optarg, NULL, 10);
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

    if (graph==-1 && lower==-1 && upper==-1) {
        usage(name);
        return EXIT_FAILURE;
    }

    if(graph!=-1){
        fprintf(stderr, "filtering graph : %d\n", graph);
    } else {

        fprintf(stderr, "filtering from ");
        if(lower!=-1)
            fprintf(stderr, "%d", lower);
        else
            fprintf(stderr, "start");
        fprintf(stderr, " until ");
        if(upper!=-1)
            fprintf(stderr, "%d.\n", upper);
        else
            fprintf(stderr, "end.\n");
    }

    /*=========== initialization ===========*/
    PREGRAPH pregraph;

    unsigned long count = 0;
    unsigned long outputted = 0;

    while(readPregraphCode(stdin, &pregraph, &endian, count)==(char)1){
        count++;
        #ifdef _DEBUG
            fprintf(stderr, "Graph %2d\n========\n", count);
        #endif
        if(count==graph){
            outputPregraph(stdout, &pregraph, LITTLE_ENDIAN, ++outputted, table);
        } else if (lower!=-1){
            if(lower <= count && (upper==-1 || count <= upper)){
                outputPregraph(stdout, &pregraph, LITTLE_ENDIAN, ++outputted, table);
            }
        } else if (count <= upper){
            outputPregraph(stdout, &pregraph, LITTLE_ENDIAN, ++outputted, table);
        }
    }

    fprintf(stderr, "Read %ld graphs\n", count);
    fprintf(stderr, "Written %ld graphs\n", outputted);
    return EXIT_SUCCESS;
}

