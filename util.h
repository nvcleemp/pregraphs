/*
 *  util.h
 *
 *
 *  Created by Nico Van Cleemput on 25/05/09.
 *
 */

#ifndef _UTIL_H //if not defined
#define _UTIL_H

#define HALFFLOOR(n) ((n)%2==0 ? (n)/2 : ((n)-1)/2)
#define MAX(a, b) (((a) >= (b)) ? (a) : (b))
#define MIN(a, b) (((a) <= (b)) ? (a) : (b))
#define BIT(i) (1 << i)

#define ERRORMSG(msg) { fprintf(stderr, "%s:%u %s\n", __FILE__, __LINE__, msg); fflush(stderr); exit(1); }

/******************Debugging macros**********************/

//#define _DEBUG

#ifdef _DEBUG

#define _DEBUGMESSAGES
#define _DEBUGDUMPS
#define _DEBUGASSERTS

#endif

//=============== MESSAGE MACRO'S ===============

#ifdef _DEBUGMESSAGES

#define DEBUGMSG(msg) { fprintf(stderr, "%s:%u %s\n", __FILE__, __LINE__, msg); fflush(stderr); }
#define DEBUGCONDITIONALMSG(condition, msg) if(condition){ fprintf(stderr, "%s:%u %s\n", __FILE__, __LINE__, msg); fflush(stderr);}

#else

#define DEBUGMSG(msg)
#define DEBUGCONDITIONALMSG(condition, msg)

#endif

//=============== DUMP MACRO'S ===============

#ifdef _DEBUGDUMPS

#define DEBUGDUMP(var, format) { fprintf(stderr, "%s:%u %s=" format "\n", __FILE__, __LINE__, #var, var); fflush(stderr); }

#define DEBUGARRAYDUMP(var, size, format) { \
                                            if(size > 0) {\
                                                fprintf(stderr, "%s:%u %s= [" format, __FILE__, __LINE__, #var, var[0]);\
                                                if(size > 0) {\
                                                    int debugarraydumpcounter;\
                                                    for(debugarraydumpcounter=1; debugarraydumpcounter<size-1; debugarraydumpcounter++){ \
                                                        fprintf(stderr, ", " format, var[debugarraydumpcounter]);\
                                                    }\
                                                    fprintf(stderr, ", " format "]", var[size-1]);\
                                                 }\
                                                 fprintf(stderr,"\n"); fflush(stderr);\
                                            } else {\
                                                fprintf(stderr, "%s:%u %s= []\n", __FILE__, __LINE__, #var);\
                                            }\
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

#else

#define DEBUGDUMP(var, format)
#define DEBUGARRAYDUMP(var, size, format)
#define DEBUG2DARRAYDUMP(var, size1, size2, format)

#endif

//=============== ASSERT MACRO'S ===============

#ifdef _DEBUGASSERTS

#define DEBUGASSERT(assertion) if(!(assertion)) {fprintf(stderr, "%s:%u Assertion failed: %s\n", __FILE__, __LINE__, #assertion); fflush(stderr); exit(1);}

#define DEBUGASSERTMSG(assertion, msg) if(!(assertion)) {fprintf(stderr, "%s:%u Assertion failed: %s\n", __FILE__, __LINE__, #assertion);\
                                                         fprintf(stderr, "%s:%u %s\n", __FILE__, __LINE__, msg); exit(1);}

#define DEBUGASSERTMSGDUMP(assertion, msg, var, format) if(!(assertion)) {fprintf(stderr, "%s:%u Assertion failed: %s\n", __FILE__, __LINE__, #assertion);\
                                                         fprintf(stderr, "%s:%u %s: " format "\n", __FILE__, __LINE__, msg, var); exit(1);}

#else

#define DEBUGASSERT(assertion)
#define DEBUGASSERTMSG(assertion, msg)
#define DEBUGASSERTMSGDUMP(assertion, msg, var, format)

#endif

#endif // end if not defined, and end the header file
