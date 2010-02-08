/*
 * File:   multi2simple.c
 * Author: nvcleemp
 *
 * Created on February 4, 2010, 9:52 AM
 *
 * Program that reads 3 regular multigraphs in multicode and translates
 * them to simple graphs.
 */
//#define _DEBUG

//TODO: currently can't handle fat K2
//TODO: only files without header can be read at the moment

#include "multi2simple.h"

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

char readMultiCodeNoHeader(FILE *f, MULTIGRAPH *multigraph, int endian) {
    int i, j, n, dummyVertex;
    unsigned short signum, number;
    if (read_old_or_new(f, FALSE, endian, &signum) == 2) {
        return (feof(f) ? 0 : 2);
    }
    //if the code starts with a zero, all the entries are two bytes
    //else the number we just read was the order of the graph
    if (signum == 0) {
        if (read_old_or_new(f, TRUE, endian, &number) == 2) {
            DEBUGMSG("Error while reading")
            return (2);
        }
    } else {
        number = signum;
    }

    if ((n = (int) number) > MAXN) {
        return (3);
    }

    //initialize the pregraph
    multigraph->order=n;
    for(i = 0; i < MAXN; i++){
        multigraph->degree[i]=0;
    }
    multigraph->multiEdgeCount=0;

    i = 1;
    dummyVertex = n;
    while (i < n) {
        if (read_old_or_new(f, signum == 0, endian, &number) == 2) {
            DEBUGMSG("Error while reading")
            return (2);
        }
        //fprintf(stderr, "%d ", number);
        DEBUGDUMP(i, "%d")
        DEBUGDUMP(number, "%d")
        if (number != 0) {
            //first check to see if we have a multi-edge
            j=0;
            while(j<multigraph->degree[i-1] && multigraph->adjList[(i-1)*3 + j] != (number-1)) j++;
            DEBUGDUMP(j, "%d")
            if(j != multigraph->degree[i-1]){
                DEBUGMSG("multi-edge")
                //we found a multi-edge (do not increase the degree)
                multigraph->multiEdgeCount++;
                multigraph->multiedge[i-1]=number-1;
                multigraph->multiedge[number-1]=i-1;
            } else {
                DEBUGMSG("simple edge")
                multigraph->adjList[(i-1)*3 + multigraph->degree[i-1]] = number-1;
                multigraph->adjList[(number-1)*3 + multigraph->degree[number-1]] = i-1;
                multigraph->degree[i-1]++;
                multigraph->degree[number-1]++;
            }
        } else {
            DEBUGMSG("next vertex")
            i++;
            //fprintf(stderr, "%d \n", number);
        }
    }
    return (1);
}

char readMultiCode(FILE *f, MULTIGRAPH *multigraph, int *endian, unsigned long count) {
    if (count == 0) {
        if (read_endian(f, endian) == 2) {
            DEBUGMSG("Error while reading")
            return (2);
        }
    }
    return (readMultiCodeNoHeader(f, multigraph, *endian));
}


char writeMulticode(FILE *f, MULTIGRAPH *multigraph, boolean header, int endian){
    //first we translate the pregraph to a simple graph
    int newOrder = multigraph->order + 2*multigraph->multiEdgeCount;
    int adjacency[newOrder*3], degree[newOrder];
    int i;
    for (i = 0; i < newOrder; i++) {
        degree[i]=0;
    }
    int offset = multigraph->order;

    for (i = 0; i < multigraph->order; i++){
        if(multigraph->degree[i]==2){
            degree[i]=3;
            if(multigraph->adjList[i*3+0]==multigraph->multiedge[i]){
                adjacency[i*3+0]=multigraph->adjList[i*3+1];
            } else {
                adjacency[i*3+0]=multigraph->adjList[i*3+0];
            }
            if(multigraph->multiedge[i]>i){
                adjacency[i*3+1]=offset;
                adjacency[i*3+2]=offset+1;
                adjacency[(multigraph->multiedge[i])*3+1]=offset;
                adjacency[(multigraph->multiedge[i])*3+2]=offset+1;
                adjacency[offset*3+0]=i;
                adjacency[offset*3+1]=multigraph->multiedge[i];
                adjacency[(offset+1)*3+0]=i;
                adjacency[(offset+1)*3+1]=multigraph->multiedge[i];
                degree[offset]=degree[offset+1]=2;
                offset+=2;
            }
        } else {
            degree[i]=3;
            adjacency[i*3+0]=multigraph->adjList[i*3+0];
            adjacency[i*3+1]=multigraph->adjList[i*3+1];
            adjacency[i*3+2]=multigraph->adjList[i*3+2];
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
    fprintf(stderr, "Usage: %s [t/m/h][X]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

/*
print a help message. name is the name of the current program.
 */
void help(char *name) {
    fprintf(stderr, "The program %s .\n", name);
    fprintf(stderr, "Usage: %s [t/m/h][X] \n\n", name);
    fprintf(stderr, "Valid options:\n");
    fprintf(stderr, "  -h          : Print this help and return.\n");
    fprintf(stderr, "  -t          : .\n");
    fprintf(stderr, "  -m          : .\n");
    fprintf(stderr, "  -X          : don't include a header in case multicode is exported.\n");
}

/*
 *
 */
int main(int argc, char** argv) {

    /*=========== commandline parsing ===========*/

    boolean table = FALSE;
    boolean multicode = FALSE;
    boolean header = TRUE;

    int c;
    char *name = argv[0];
    int endian = defaultEndian;

    while ((c = getopt(argc, argv, "htmdX")) != -1) {
        switch (c) {
            case 't':
                fprintf(stderr, "-t currently not supported.\n");
                table = TRUE;
                break;
            case 'm':
                multicode = TRUE;
                break;
            case 'd':
                debugOutput = TRUE;
                break;
            case 'X':
                header = FALSE;
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
    MULTIGRAPH multigraph;

    unsigned long count = 0;

    while(readMultiCodeNoHeader(stdin, &multigraph, endian)==(char)1){
        count++;
        if(table){
            //writeGraphTable(stdout, &multigraph, count);
        } else if (multicode)
            writeMulticode(stdout, &multigraph, header && count==1, endian);
    }

    return EXIT_SUCCESS;
}
