/* 
 * File:   minibaumlib.h
 * Author: nvcleemp
 *
 * Created on February 23, 2011, 10:17 AM
 */

#ifndef MINIBAUMLIB_H
#define	MINIBAUMLIB_H

# define knoten 62   /* maximal value 63 in this version */
# define reg 3

int is_minibaum_available(int order);
void call_minibaum(int order, int isBipartite, int allGraphs);
void handle_minibaum_result(unsigned char minibaum_graph[knoten+1][reg], int order);

#endif	/* MINIBAUMLIB_H */

