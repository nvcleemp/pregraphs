/* 
 * File:   multicode2dreadnaut.h
 * Author: nvcleemp
 *
 * Created on February 9, 2010, 9:23 AM
 */

#ifndef _MULTICODE2DREADNAUT_H
#define	_MULTICODE2DREADNAUT_H

#ifndef MAXN
#define MAXN 50
#endif

#include "util.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>

/******************Defines*******************************/

#define TRUE 1
#define FALSE 0

typedef int boolean;

/******************Global Variables**********************/

short defaultEndian = LITTLE_ENDIAN; // defines which endian should be used while exporting pregraph code

boolean debugOutput = FALSE;


#endif	/* _MULTICODE2DREADNAUT_H */

