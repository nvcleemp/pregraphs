/*
 * File:   multicode2dreadnaut.c
 * Author: nvcleemp
 *
 * Created on February 9, 2010, 9:52 AM
 *
 * Program that reads graphs in multicode and prints dreadnaut commands.
 */

#include "multicode2dreadnaut.h"

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
    //fprintf(outputFile, "%d ", *number);
    return (1);
}

char readAndWriteMultiCodeNoHeader(FILE *f, FILE *out, int endian) {
    int i, n;
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

    fprintf(out, "n=%d $=1 g\n", number);

    if ((n = (int) number) > MAXN) {
        return (3);
    }

    fprintf(out, "  1: ");
    i = 1;
    while (i < n) {
        if (read_old_or_new(f, signum == 0, endian, &number) == 2) {
            DEBUGMSG("Error while reading")
            return (2);
        }
        if (number != 0) {
            fprintf(out, "%d ", number);
        } else {
            DEBUGMSG("next vertex")
            i++;
            if(i==n)
                fprintf(out, ";\n");
            else
                fprintf(out, ";\n%3d: ", i);
        }
    }
    fprintf(out, "\n");
    return (1);
}

/*
 *
 */
int main(int argc, char** argv) {

    /*=========== commandline parsing ===========*/

    //int c;
    //char *name = argv[0];
    int endian = defaultEndian;

    /*
    while ((c = getopt(argc, argv, "h")) != -1) {
        switch (c) {
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            default:
                usage(name);
                return EXIT_FAILURE;
        }
    }
     */

    // check the non-option arguments
    //if (argc - optind != 0) {
    //    usage(name);
    //    return EXIT_FAILURE;
    //}

    unsigned long count = 0;

    while(readAndWriteMultiCodeNoHeader(stdin, stdout, endian)==(char)1){
        count++;
    }

    return EXIT_SUCCESS;
}
