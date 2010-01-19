/*
 *  util.h
 *
 *
 *  Created by Nico Van Cleemput on 25/05/09.
 *
 */

#ifndef _UTIL_H //if not defined
#define _UTIL_H

#define HALFFLOOR(n) (n%2==0 ? n/2 : (n-1)/2)

/******************Debugging macros**********************/

#define _DEBUG

#ifdef _DEBUG

#define DEBUGMSG(msg) { fprintf(stderr, "%s:%u %s\n", __FILE__, __LINE__, msg); fflush(stderr); }


#define DEBUGCONDITIONALMSG(condition, msg) if(condition){ fprintf(stderr, "%s:%u %s\n", __FILE__, __LINE__, msg); fflush(stderr);}

#define DEBUGDUMP(var, format) { fprintf(stderr, "%s:%u %s=" format "\n", __FILE__, __LINE__, #var, var); fflush(stderr); }

#define DEBUGARRAYDUMP(var, size, format) { \
                                            fprintf(stderr, "%s:%u %s= [" format, __FILE__, __LINE__, #var, var[0]);\
                                            int debugarraydumpcounter;\
                                            for(debugarraydumpcounter=1; debugarraydumpcounter<size-1; debugarraydumpcounter++){ \
                                                fprintf(stderr, ", " format, var[debugarraydumpcounter]);\
                                            }\
                                            fprintf(stderr, ", " format "]\n", var[size-1]); fflush(stderr);\
                                          }

#define DEBUG2DARRAYDUMP(var, size1, size2, format) { \
                                            fprintf(stderr, "%s:%u %s=\n", __FILE__, __LINE__, #var);\
                                            int debug2darraydumpcounter1,debug2darraydumpcounter2;\
                                            for(debug2darraydumpcounter1=0; debug2darraydumpcounter1<size1; debug2darraydumpcounter1++){ \
                                              fprintf(stderr, "[" format, var[debug2darraydumpcounter1][0]);\
                                              for(debug2darraydumpcounter2=1; debug2darraydumpcounter2<size2-1; debug2darraydumpcounter2++){ \
                                                 fprintf(stderr, ", " format, var[debug2darraydumpcounter1][debug2darraydumpcounter2]);\
                                              }\
                                              fprintf(stderr, ", " format "]\n", var[debug2darraydumpcounter1][size2-1]);\
                                            }\
                                            fflush(stderr);\
                                          }

#define DEBUGASSERT(assertion) if(!(assertion)) {fprintf(stderr, "%s:%u Assertion failed: %s\n", __FILE__, __LINE__, #assertion); fflush(stderr); exit(1);}

#else

#define DEBUGMSG(msg)

#define DEBUGCONDITIONALMSG(condition, msg)

#define DEBUGDUMP(var, format)

#define DEBUGARRAYDUMP(var, size, format)

#define DEBUG2DARRAYDUMP(var, size1, size2, format)

#define DEBUGASSERT(assertion)

#endif

#endif // end if not defined, and end the header file
