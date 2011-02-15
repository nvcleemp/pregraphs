/*
 * File:   admissable_c.c
 * Author: nvcleemp
 *
 * Created on June 9, 2010, 9:52 AM
 *
 * Checks whether the given pregraph has an admissable colouring. The only 
 * accepted pregraphs are multigraphs with semi-edges. An admissable colouring
 * corresponds to a 2-factor where each component is a quotient of a 4-cycle.
 */
//#define _DEBUG

#include "admissable_c.h"

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

void writePregraphTableDebug(FILE *f, PREGRAPH *pregraph, boolean *removed){
    int i, j;
    for(i=0; i<pregraph->order; i++){
        if(!removed[i]){
            fprintf(f, "%d) ", i+1);
            for(j=0; j<3; j++){
                fprintf(f, "%d ", pregraph->adjList[i][j]+1);
            }
            fprintf(f, "\n");
        }
    }
    fprintf(f, "\n");
}

//=============================================================================
//=============================================================================

int getThirdNeighbour(PREGRAPH *pregraph, int v, int n1, int n2){
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
int getSquare(PREGRAPH *pregraph, int v, int n1, int n2){
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
void finishChain(PREGRAPH *pregraph, boolean *removed, int v1, int v2, int v3,
        int v4, boolean *admissable){
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

	if(pregraph->adjMatrix[nextUp] & (1<<nextDown)){
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
void finishChainWithout(PREGRAPH *pregraph, boolean *removed, int v1, int v2,
        int v3, int v4){
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

	if(pregraph->adjMatrix[nextUp] & (1<<nextDown)){
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
void detectAndRemoveChain(PREGRAPH *pregraph, boolean *removed, int v1, int v2,
        int v3, int v4, boolean *admissable){
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
        finishChain(pregraph, removed, v1, v2, v3, v4, admissable);
	return;
    } else if (n1==v3){
	finishChain(pregraph, removed, v1, v3, v2, v4, admissable);
	return;
    } else if(n4==v2){
        finishChain(pregraph, removed, v4, v2, v3, v1, admissable);
	return;
    } else if (n4==v3){
        finishChain(pregraph, removed, v4, v3, v2, v1, admissable);
	return;
    }

    //check for semi-edges
    if(n1 == semiEdge && n2==semiEdge){
        finishChain(pregraph, removed, v1, v2, v3, v4, admissable);
	return;
    } else if (n1==semiEdge && n3==semiEdge){
	finishChain(pregraph, removed, v1, v3, v2, v4, admissable);
	return;
    } else if(n4== semiEdge && n2==semiEdge){
        finishChain(pregraph, removed, v4, v2, v3, v1, admissable);
	return;
    } else if (n4==semiEdge && n3==semiEdge){
        finishChain(pregraph, removed, v4, v3, v2, v1, admissable);
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

    //determine sides of chains
    if(pregraph->adjMatrix[n1] & (1<<n2)){
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
    } else if(pregraph->adjMatrix[n1] & (1<<n3)){
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
    } else if(pregraph->adjMatrix[n4] & (1<<n2)){
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
                    admissable);
	    return;
	}
	if(nextUp==currentDown){
	    //end of chain is reached because multi-edge is found
	    //finish chain in opposite direction
	    finishChain(pregraph, removed, side1_a, side1_b, side2_a, side2_b,
                    admissable);
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
        } else if(pregraph->adjMatrix[nextUp] & (1<<nextDown)){
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
        } else if(pregraph->adjMatrix[nextUp] & (1<<nextDown)){
	    //still in the chain of square
	    prevUp = currentUp;
	    prevDown = currentDown;
	    currentUp = nextUp;
	    currentDown = nextDown;
	    nrOfSquares++;
	} else {
	    //end of chain reached because nextUp and nextDown aren't adjacent
	    *admissable = nrOfSquares % 2;
	    inChain = FALSE;
	}
    }
}

void removeChainsOfSquares(PREGRAPH *pregraph, boolean *removed,
        boolean *admissable){
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
                                oppositeCorner, admissable);
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
void checkAndRemoveMultiEdge(PREGRAPH *pregraph, boolean *removed, int v,
        boolean *admissable){
    if(pregraph->adjList[v][0]==pregraph->adjList[v][1]){
	removed[v]=TRUE;
	removed[pregraph->adjList[v][0]]=TRUE;
	*admissable = !(pregraph->adjMatrix[pregraph->adjList[v][0]] &
                (1<<pregraph->adjList[v][2]));
    } else if(pregraph->adjList[v][0]==pregraph->adjList[v][2]){
	removed[v]=TRUE;
	removed[pregraph->adjList[v][0]]=TRUE;
	*admissable = !(pregraph->adjMatrix[pregraph->adjList[v][0]] &
                (1<<pregraph->adjList[v][1]));
    } else if(pregraph->adjList[v][1]==pregraph->adjList[v][2]){
	removed[v]=TRUE;
	removed[pregraph->adjList[v][1]]=TRUE;
	*admissable = !(pregraph->adjMatrix[pregraph->adjList[v][1]] &
                (1<<pregraph->adjList[v][0]));
    }
}

void removeMultiEdges(PREGRAPH *pregraph, boolean *removed, boolean *admissable){
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
		checkAndRemoveMultiEdge(pregraph, removed, i, admissable);
		if(!(*admissable)) return;
	    }
	}
    }
}

//assumes that no multi-edges are left
boolean checkRemainder(PREGRAPH *pregraph, boolean *removed){
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
            int nextVertex;
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

boolean existsVertexWithMultipleSemiEdges(PREGRAPH *pregraph){
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

boolean isColourable(PREGRAPH *pregraph){
    int i;
    boolean removed[MAXN];
    for(i=0; i<pregraph->order; i++) removed[i] = FALSE;

    //if all vertices are adjacent with at least one semi-edge, then the graph
    //has an admissable colouring if there are an even number of vertices
    i=0;
    while(i<pregraph->order && (pregraph->adjMatrix[i] & (1<<(pregraph->order))))
        i++;
    if(i==pregraph->order) return !(pregraph->order%2) ||
            existsVertexWithMultipleSemiEdges(pregraph);

    boolean admissable = TRUE;

    removeChainsOfSquares(pregraph, removed, &admissable);

    if(!admissable) return FALSE;

    removeMultiEdges(pregraph, removed, &admissable);

    if(!admissable) return FALSE;

    return checkRemainder(pregraph, removed);
}

//=============================================================================
//=============================================================================

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
    fprintf(stderr, "The program %s filters pregraphs that have an admissable colouring.\n", name);
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
    unsigned long allowsAdmissableColouring = 0;

    while(readPregraphCode(stdin, &pregraph, &endian, count)==(char)1){
        count++;
        //DEBUGDUMP(count, "%ld")
	//determine admissable colourable
        #ifdef _DEBUG
	fprintf(stderr, "Graph %2d\n========\n", count);
        #endif
	if(isColourable(&pregraph)){
            #ifdef _DEBUG
	    fprintf(stderr, "Admissable\n\n");
            #endif
	    allowsAdmissableColouring++;
            if(!onlyCount)
                writePregraphCode(stdout, &pregraph, LITTLE_ENDIAN,
                        allowsAdmissableColouring);
        }
    }

    if(outputFileName != NULL){
        outputFile = fopen(outputFileName, "a");
        if(outputFile == NULL){
            fprintf(stderr, "File %s couldn't be created: aborting!\n",
                    outputFileName);
            return EXIT_FAILURE;
        }
    } else {
        outputFile = stderr;
    }

    fprintf(outputFile, "Read %ld graphs ", count);
    fprintf(outputFile, "of which %ld graphs allow an admissable colouring.\n",
            allowsAdmissableColouring);
    return EXIT_SUCCESS;
}

