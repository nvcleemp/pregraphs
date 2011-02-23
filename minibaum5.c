/* 2.12.2003: Speicherfehler korrigiert, der bei kleinen Knotenzahlen
(<10) und singleout auftrat. Fehler war nur bei den Parametern minibaum 4 3 s
festzustellen -- Effekt: Es wurde immer shortcode ausgegeben. Andere Fehler
sind aber fuer |V|<10 und Option s nicht auszuschliessen */


/* 26.10.2000: Problem with graphs being smaller than minlevel being
multiply outputted when using graph6code resolved (output line was at 
wrong place) */

/* From gunnar@Mathematik.Uni-Bielefeld.DE  Tue Aug 13 23:45:49 1996
   To: bdm@cs.anu.edu.au (Brendan McKay)
   Date: Tue, 13 Aug 1996 15:31:34 +0200 (MET DST)
*/
/*
Ausserdem sind noch die optionen "S" und "Zx" hinzugekommen. "S" fuer
snarks -- d.h. nicht kanten 3-faerbbare Graphen und "Zx" das fuer zyklisch
mindestens x-fach zusammenhaengende Graphen.
Diese beiden werden nur ueber filter implementiert und nicht besonders effektiv. */



/* dieses programm baut auf minibaum2.c auf. Es ist lediglich die option
   "i x" hinzugekommen, die es ermoeglicht, dass nur graphen erzeugt werden,
   die keine unabhaengige Menge der groesse x+1 haben.
   Der Testprozess wird aehnlich dem minimalitaetstest parallelisiert.

 */
/*****MINIBAUM.C*********/
/* Diese Version enthaelt eine subroutine, die schon relativ frueh
   verbotene Kanten markiert                                       */
/* Hier gibt es auch die Option v fuer die Erzeugung ausschliesslich
   vertex-transitiver Graphen -- funktioniert aber miserabel */


# include <stdio.h>
# include <stdlib.h>
//# include <malloc.h>
# include <time.h>
# include <sys/times.h>
# include <string.h>
# include <limits.h>
# include <unistd.h>

#include "minibaumlib.h"

/*******DEFINITIONEN***************************************/

/*******May*be*changed:********************************/

# define baumgrenze 30
# define indepspeichergrenze 22

#define listenlaenge 5000

/******Should*not*be*changed:*****************************/

#define Sicherungsintervall 1800

                /*At least "Sicherungsintervall" seconds between two savings*/
                /* 0 means: Do not save the lastgraph to recover */
/***********May*not*be*changed:*****************************/


#define unbelegt 98

#define leer 99

#define nil 0

#define codelaenge (knoten*reg)/2+1 /* +1 fuer shortcodes im worst case
				       und dort ja auch ohne -reg */

#define time_factor sysconf(_SC_CLK_TCK)


typedef unsigned char GRAPH[knoten+1][reg];


typedef unsigned char BILD[knoten+1];
typedef unsigned char ADJAZENZ[knoten+1];
typedef unsigned char BILDELEMENT;

typedef struct ind_set {
                unsigned long long int drin; /* Welche Knoten sind erlaubt ? */
                unsigned char groesse; /* wieviel Knoten sind schon in der
                                          unabhaengigen Menge */
	        } IND_SET;


typedef struct wurzel {
            struct zweig *(naechster[6]);
            unsigned char tot[6];
            unsigned char wurzeltot;
            unsigned char tiefe;   /* in welcher tiefe wurde diese wurzel
                                       initialisiert (d.h.: Wie gross ist die
                                       anzahl der kanten im graphen schon ) */
            } ROOT;

typedef struct zweig {
             signed char blatt;            /* 1, wenn es ein blatt ist, 0 sonst */
             BILD bild;
             BILD urbild;
             unsigned char totlinks;/* gibt an, ob es noch Sinn hat, in diese*/
             unsigned char totrechts; /* Richtung weiterzusuchen, oder ob da
                                       keine Moeglichkeit besteht */
             struct zweig *links;   /* Die Nachfolger */
             struct zweig *rechts;
             struct zweig *vorgaenger;
             ROOT *zurwurzel;
             unsigned char nextbild;/* die naechste zu vergebende nummer */
             unsigned char wo;      /* welchen "wahren" Namen hat der naechste
                                       Knoten, an dem noch wahlmoeglichkeiten
                                       fuer bild bestehen */
             unsigned char tiefe;   /* in welcher tiefe wurde dieser knoten
                                       eingefuegt (d.h.: Wie gross ist die
                                       anzahl der kanten im graphen schon ) */
             } BRANCH;

/* Ein neuer Zweig wird eingefuegt, wenn das bild sich aendert ( dann werden
   ein oder zwei nachfolger "gebaut" ) oder auch, wenn es einfach weitergeht --
   d.h.: Durch neu hinzugekommene kanten ist die adjazenz des naechsten
   Knotens jetzt reg und man kann den Baum ausbauen */


/**************************GLOBALE VARIABLEN************************/


ADJAZENZ adj;
int knotenzahl,taille, border, treeborder, indepborder, snarks, zyconnecti, zyminconn;
int noout_minibaum=0;
int indepstoreborder, maximalkanten;
int singleout_minibaum, lastout, recover, startvorgabe, bipartite; /* die Optionen */
int graph6_output_minibaum=0;
int partial, timerestricted, TIMERESTRICTED, observer, short_code;    /* die Optionen */
int vertex_trans, independent, modulo_minibaum, mod_minibaum, rest_minibaum;  /* die Optionen */
int transtest; /* eine variable, die beim Test auf vertextransitivitaet 
                   hilft */
int maxindep, splitlevel_minibaum, minwrite;
unsigned int splitzaehler;
unsigned char * globalliste[knoten+1];
int codelength[knoten+1];
unsigned long long int number_of_graphs;
unsigned char forbidden[knoten+1][knoten+1];
int startzahl; /* die kleinste ueberhaupt fuer die gesuchten graphen 
                  moegliche zahl */
char codename[45], lastgraphname[45], startgraphname[45], endgraphname[45];
unsigned long long int graphzahlen[knoten+1];
unsigned char * nextposition[knoten+1];
unsigned char oldcode1[knoten+1][codelaenge];
unsigned char oldcode2[knoten+1][codelaenge];
unsigned char whichcode[knoten+1];
unsigned long long int MASK[66], NMASK[66];
ROOT roots[knoten+1];
GRAPH lastgraph,endgraph;
signed char colour[knoten+1]; /* Fuer bipartit: Die Farben: 1 und 2;  0 bedeutet
                           unbelegt */
int anzahl_colours[3]; /* wieviel gibt es mit colour 1, 2 */
unsigned char mindestkantenzahl[knoten+2][11]; 
/* mindestkantenzahl[i][j] gibt an, wieviel kanten mindestens schon im Graphen
   sein muessen, wenn man bei knoten i ist, der graph taille j bekommen soll
   und option singleout benutzt wurde */


int blaetterzahl[knoten+1]; /* die jeweilige anzahl der in den baeumen
                               vorhandenen Blaetter; blaetterzahlen[0][0]
                               enthaelt den wert von g[0][1], wann die
                               baumanalyse durchgefuehrt wurde (default: 0) */
BRANCH **(blaetterkette[knoten+1]);
IND_SET *(indeps[knoten+1]); /* indeps[x] wird ein array aller unabhaengigen
                                Mengen sein, die maximal das element x 
                                enthalten und bei denen sich etwas gegenueber
                                Schritt x-1 geaendert hat -- d.h.: die x auch
                                wirklich enthalten. indeps[x] wird 
                                belegt, nachdem Knoten x fertig ist .
                                indeps[1]..indeps[x-1] sind auch noch "gueltig"
                                -- es sind halt die Mengen, die x nicht 
                                enthalten */

unsigned int numberofinds[knoten+2]; /* Bei der Belegung von indeps[x] wird
                                jeweils ueberprueft, ob x+1 verboten ist oder
                                nicht. Wenn nicht, so wird ein element fuer
                                indeps[x+1] gebraucht. */

unsigned char five14s[knoten+1];

unsigned char minimalwerte[knoten+1][knoten+1][knoten+1]; /* die wirklich
                             berechneten minimalwerte fuer independence --
                             eine optimierung gegenueber five14s */

int starttime, endtime; /* Die Zeit wann das programm gestartet, bzw */
                        /* abgeschossen wird */

unsigned int oldtime;
int S_intervall;


int cyc_test_knotenzahl, cyc_test_kantenzahl;

int forbiddencycles=0; /* contains the option and also the length of the
			  longest forbidden cycle */
char forbcyc[knoten+1]={0};


/********************ADDVERTEX***********************************/

void addvertex(GRAPH graph, int vertex, 
	       unsigned long long int neighbors, unsigned long long int *newneighbors,
	       unsigned int  long long elements, unsigned long long int *newelements, 
	       int *newnum_elements,
	       int *newlkantenzahl, int *newcutsize, 
	       unsigned  long long int noelements, char *fehler)
/* Fuegt einen Knoten in die ZK-Menge von Knoten 1 ein */

{
int j;
unsigned long long int pufferMASK;
int forced[2], forcedzaehler=0; /* fuer zwangslaeufig einzufuegende nachbarn */
int forb_neighb=0;

forced[0]=forced[1]=0;

*newelements |= MASK[vertex];
(*newnum_elements)++;

for (j=0; j<reg; j++) 
   { pufferMASK=MASK[graph[vertex][j]];
          if (pufferMASK & noelements) forb_neighb++;
	  if (!((pufferMASK & neighbors)||(pufferMASK & elements)))
	                 *newneighbors |= pufferMASK; 
	  if (pufferMASK & elements) { (*newlkantenzahl)++; (*newcutsize)--; }
	         else if (pufferMASK & neighbors)
		        { (*newcutsize)++;
			  /* dann muss aber fuer einen guten cut der nachbar
			     auch mit rein */
			  forced[forcedzaehler]=graph[vertex][j]; 
			  forcedzaehler++; }
	                 else /* d.h. ein neues cutelement -- weder
			         zu nachbar noch zu element */ 
			      (*newcutsize)++;
   }
/* jetzt sollten die Kanten abgearbeitet sein */

if (forb_neighb == 2) { *fehler=1; return; }
   /* es kann kein guter cut mehr werden */

for (j=0; j<forcedzaehler; j++)
    { elements |= *newelements; 
      neighbors = (neighbors | *newneighbors) & ~(*newelements);
      /* Das kann sich auch nach dem ersten Durchlauf der Schleife geaendert
	 haben. */
      if (MASK[forced[j]] & noelements) *fehler=1;
      /* wenn ein verbotener knoten gezwungen wird -- das kann nicht gut sein */
      if (*fehler) return; /* muss nicht von einer zeile drueber kommen */
      if (!(MASK[forced[j]] & elements)) /* fuer den 2. durchlauf wichtige
					    Abfrage */
      addvertex(graph, forced[j],neighbors,newneighbors,elements,newelements,
              newnum_elements,newlkantenzahl,newcutsize,noelements,fehler);
    }
}



void unsignedintegerschreib(unsigned int number, FILE *wohin)
{
  unsigned char puffer[4];
  int i, mark;;

      puffer[0]= (number >> 24) & 255U;
      puffer[1]= (number >> 16) & 255U;
      puffer[2]= (number >> 8) & 255U;
      puffer[3]= number & 255U;
      
      if (fwrite(puffer,1,4,wohin) != 4)
	{ for (i=mark=0; (i<10) && !mark; i++) 
	    if (fwrite(puffer,1,4,wohin) != 4) sleep(3);
	    else mark=1;
	  if (!mark) { 
	    fprintf(stderr,"Dangerous error -- could not write\n");
	    exit(0); }
	}
}

void unsigned_ll_integerschreib(unsigned long long int number, FILE *wohin)
{
  unsigned char puffer[8];
  int i, mark;;

  puffer[0]= (number >> 56) & 255ULL;
  puffer[1]= (number >> 48) & 255ULL;
  puffer[2]= (number >> 40) & 255ULL;
  puffer[3]= (number >> 32) & 255ULL;
  puffer[4]= (number >> 24) & 255ULL;
  puffer[5]= (number >> 16) & 255ULL;
  puffer[6]= (number >> 8) & 255ULL;
  puffer[7]= number & 255ULL;
      
      if (fwrite(puffer,1,8,wohin) != 8)
	{ for (i=mark=0; (i<10) && !mark; i++) 
	    if (fwrite(puffer,1,8,wohin) != 8) sleep(3);
	    else mark=1;
	  if (!mark) { 
	    fprintf(stderr,"Dangerous error -- could not write\n");
	    exit(0); }
	}
}

void unsignedintegerlies(unsigned int *number, FILE *woher)
{
  unsigned char puffer[4];
  int i, mark;

      if (fread(puffer,1,4,woher) != 4)
	{ for (i=mark=0; (i<10) && !mark; i++) 
	    if (fread(puffer,1,4,woher) != 4) sleep(3);
	    else mark=1;
	  if (!mark) { 
	    fprintf(stderr,"Dangerous error -- could not read\n");
	    exit(0); }
	}
  
      *number= (((unsigned int)puffer[0]) <<24) + (((unsigned int)puffer[1]) <<16) +
	(((unsigned int)puffer[2]) <<8) + ((unsigned int)puffer[3]); 
}

void unsigned_ll_integerlies(unsigned long long int *number, FILE *woher)
{
  unsigned char puffer[8];
  int i, mark;

  if (fread(puffer,1,8,woher) != 8)
    { for (i=mark=0; (i<10) && !mark; i++) 
	if (fread(puffer,1,8,woher) != 8) sleep(3);
	else mark=1;
      if (!mark) { 
	fprintf(stderr,"Dangerous error -- could not read\n");
	exit(0); }
    }
  
  *number= (((unsigned long long int)puffer[0]) <<56) + (((unsigned long long int)puffer[1]) <<48) +
    (((unsigned long long int)puffer[2]) <<40) + (((unsigned long long int)puffer[3]) <<32) +
    (((unsigned long long int)puffer[4]) <<24) + (((unsigned long long int)puffer[5]) <<16) +
    (((unsigned long long int)puffer[6]) <<8)  + ((unsigned long long int)puffer[7]); 


}

/**********************WEITERWURSTELN*****************************/

void weiterwursteln(GRAPH graph,unsigned  long long int neighbors, 
		    int num_elements, unsigned  long long int elements,
                    unsigned  long long int noelements, int cutsize,int lkantenzahl,
		    int *minimum, int *nurtriv)
{
int i,j,minsize,newlkantenzahl, defekt;
unsigned  long long int newneighbors, newelements;
int newnum_elements, newcutsize, puffer;
char fehler;

/* Erst testen, ob ueberhaupt noch die Moeglichkeit besteht, einen guten cut
   zu finden: */
if (num_elements> (cyc_test_knotenzahl/2)) return; /* konstruiert wird immer der kleinere
					     der beiden Teile */
/*if (restconnected(graph,elements)<cyc_test_knotenzahl-num_elements) 
  dann muessen alle bis auf eine der komponenten "einverleibt werden"*/

if (cyc_test_knotenzahl-num_elements > cyc_test_kantenzahl-lkantenzahl-cutsize) return;

defekt=cutsize-(*minimum)+1;
if (num_elements + defekt > cyc_test_knotenzahl/2) return;

minsize=0;
for (i=1; i<= cyc_test_knotenzahl; i++)
 { if (MASK[i] & noelements)
	{  for (j=0; j<reg; j++) 
	      if (MASK[graph[i][j]] & elements) minsize++; 
	      /* diese kanten liegen in jedem hieraus konstruierten cut */ }
    else
/* es darf nie ein element zu 2 verbotenen adjazent sein -- das waere nie ein
   guter cut */
     if (MASK[i] & elements)
	 { puffer=0;
	   for (j=0; j<reg; j++) 
	   if (MASK[graph[i][j]] & noelements) puffer++;
           if (puffer>=2) { return; }
       }
 }
      
if (minsize> *minimum) return; /* == minimum muss drinbleiben, um nichttriviale cuts 
				  zu finden, wenn schon ein trivialer gefunden wurde,
				  aber noch kein nicht-trivialer */
if ((*nurtriv==0) && (minsize == *minimum)) return;

/* Dann testen, ob ein besserer Cut gefunden wurde: */
if ((num_elements<cyc_test_knotenzahl) && (lkantenzahl>=num_elements) && 
           (cyc_test_kantenzahl-lkantenzahl-cutsize >= cyc_test_knotenzahl-num_elements))
/* d.h. beide Teile enthalten einen Zykel */
{ if (cutsize<*minimum)  { *minimum=cutsize; 
			   if (num_elements==cutsize) /* es ist die kleinere von beiden,
							 die andere muss also nicht ueberprueft
							 werden */
			     *nurtriv=1; else *nurtriv=0;
			 }
 else if (cutsize== *minimum)  { if (num_elements>cutsize) *nurtriv=0; }
}


/* Jetzt geht's los: */
for (i=2; i<=cyc_test_knotenzahl && (!(MASK[i] & neighbors) || (MASK[i] & noelements));
     i++); /* suche ersten noch nicht festgelegten nachbarn */

if (i>cyc_test_knotenzahl) {/*alles festgelegt:*/ return;}

     /* zunaechst Knoten rein: */
      fehler=0;
      newneighbors=newelements=0;
      newlkantenzahl=newcutsize=newnum_elements=0;
      addvertex(graph, i,neighbors,&newneighbors,elements,&newelements, 
	      &newnum_elements, &newlkantenzahl,&newcutsize,noelements,&fehler);
      if (!fehler)
      weiterwursteln(graph, (neighbors | newneighbors)& ~newelements, 
		     num_elements + newnum_elements,
		     elements | newelements, noelements, cutsize+newcutsize,
		     lkantenzahl+newlkantenzahl, minimum, nurtriv);
      if (minsize>=*minimum) return;
      /* ansonsten diesen knoten verbieten und weiter */
      weiterwursteln(graph,neighbors, num_elements, elements, 
		     noelements | MASK[i], cutsize, lkantenzahl, minimum, nurtriv);

}

/**************************CYC_CON********************************/

int cyc_con(GRAPH graph, int *nurtriv)
{
/* strategie: Ein minimaler cut teilt die vertexmenge in 2 zusammenhaengende
   Teile, von denen jeder mindestens einen Zykel enthaelt. Konstruiere den, 
   der am kleinsten ist. */


int i, minknoten;
unsigned long long int neighbors, elements, noelements;
int minimum=knoten, cutsize, lkantenzahl, num_elements, forb_n;

cyc_test_knotenzahl=graph[0][0];
cyc_test_kantenzahl=graph[0][1];


/* lkantenzahl ist die Anzahl der inneren Kanten in der konstruierten Komponente */

*nurtriv=1; /* wird gegebenenfalls auf 0 gesetzt */

for (minknoten=1; minknoten <= cyc_test_knotenzahl-3; minknoten++)
{
/* zuerst nur Knoten minknoten rein: */


elements=MASK[minknoten]; noelements=neighbors=0;
for (i=1; i<minknoten; i++) noelements |= MASK[i]; /* Elemente, die nicht drin sein DUERFEN */
forb_n=0;
for (i=0; i<reg; i++) { neighbors |= MASK[graph[minknoten][i]];
			if (MASK[graph[minknoten][i]] & noelements) forb_n++; }
if (forb_n<2) /* sonst braucht nichts gemacht zu werden, da es besser waere, den Knoten noch
	         mit zu dem anderen Cut zu nehmen */
{ cutsize=3; lkantenzahl=0; num_elements=1;

weiterwursteln(graph,neighbors,num_elements,elements,noelements,cutsize,lkantenzahl,&minimum,nurtriv); } /* ende if */
}/* ende for */


if (minimum<knoten) return(minimum); else return(0);
}




/*******************SIGNALHANDLER****************************/

#include <signal.h>
#include <sys/time.h>

/* SA_OLDSTYLE */

#ifndef	SA_RESTART
#define	SA_RESTART	0
#endif

#ifdef __STDC__
volatile	int	flag = 0;

void	sigcatch( int sig )
{
	flag = 1;
}

int		setup( void )

#else

int flag = 0;

void sigcatch (sig)
int sig;
{
	flag = 1;
}

int		setup()
#endif

{
	struct	sigaction	action, old;
	struct	itimerval	timer;

	action.sa_handler = sigcatch;
	action.sa_flags   = SA_RESTART;
	sigemptyset( &action.sa_mask );

	if( sigaction( SIGALRM, &action, &old ) == -1 )
	{
		perror( "sigaction" );
		return( -1 );
	}

	timer.it_interval.tv_sec  = S_intervall;
	timer.it_interval.tv_usec = 0;
	timer.it_value.tv_sec     = timer.it_interval.tv_sec;
	timer.it_value.tv_usec    = timer.it_interval.tv_usec;

	if( setitimer( ITIMER_REAL, &timer, NULL ) == -1 )
	{
		perror( "setitimer" );
		sigaction( SIGALRM, &old, NULL );
		return( -1 );
	}

	return( 0 );
}




/**************************DIE FUNKTIONEN*******************************/

void schreibemenge( menge )
unsigned  long long int menge;
{
int i;
fprintf(stderr,"\n");
for(i=1;i<=64;i++) if (menge & MASK[i]) fprintf(stderr," %d ",i);
fprintf(stderr,"\n");
}

void writeset( menge )
unsigned  long long int menge;
{
int i;
fprintf(stderr,"\n");
 for(i=1;i<=64;i++) if (menge & (1ULL<<i)) fprintf(stderr," %d ",i);
fprintf(stderr,"\n");
}

void schreibegraph(g)
GRAPH g;
{
 int x,y;
FILE *fil2;

fil2=stderr;   /*fopen("broadcastoutput","a");*/
fprintf(fil2," ");

fprintf(fil2,"|%2d",g[0][0]);

for(x=1; (x <= g[0][0])&&(x<=24); x++)  fprintf(fil2,"|%2d",x); fprintf(fil2,"|\n");

fprintf(fil2," ");

for(x=0; (x <= g[0][0])&&(x<=24); x++) fprintf(fil2,"|==");    fprintf(fil2,"|\n");

for(x=0; x < reg; x++)
  {
   fprintf(fil2," |  ");
   for(y=1; (y<=g[0][0])&&(y<=24); y++)
       fprintf(fil2,"|%2d",g[y][x]);
       fprintf(fil2,"|\n");
  }

if (g[0][0]>24) 
{
fprintf(fil2,"\n");

fprintf(fil2,"    ");

for(x=25; (x <= g[0][0])&&(x<=48); x++)  fprintf(fil2,"|%2d",x); fprintf(fil2,"|\n");

fprintf(fil2,"    ");

for(x=25; (x <= g[0][0])&&(x<=48); x++) fprintf(fil2,"|==");    fprintf(fil2,"|\n");

for(x=0; x < reg; x++)
  {
   fprintf(fil2,"    ");
   for(y=25; (y <= g[0][0])&&(y<=48); y++)
       fprintf(fil2,"|%2d",g[y][x]);
       fprintf(fil2,"|\n");
  }
}

if (g[0][0]>48) 
{
fprintf(fil2,"\n");

fprintf(fil2,"    ");

for(x=49; x<=g[0][0]; x++)  fprintf(fil2,"|%2d",x); fprintf(fil2,"|\n");

fprintf(fil2,"    ");

for(x=49; x <= g[0][0]; x++) fprintf(fil2,"|==");    fprintf(fil2,"|\n");

for(x=0; x < reg; x++)
  {
   fprintf(fil2,"    ");
   for(y=49; y<=g[0][0]; y++)
       fprintf(fil2,"|%2d",g[y][x]);
       fprintf(fil2,"|\n");
  }
}
/*fclose(fil2);*/
}

/****************************INDEPNORMAL**********************************/

void indepnormal( graph, set, wo, maxwo, abbruch)
/* wird nur aufgerufen fuer einen vertex wo, der die wahl laesst zwischen
   einfuegen und nicht einfuegen -- darf aber nicht fuer 1 aufgerufen werden */
GRAPH graph;
IND_SET set;
unsigned char wo, maxwo;
int *abbruch;

{
unsigned char i,j,k,rest,puffer;
signed char z;
unsigned char valf[4];

if (taille==3)
{
/* zunaechst nicht einfuegen: */
for (i=wo+1; (i<=maxwo) && !(set.drin & MASK[i]) ; i++);
/* sucht naechsten vertex, wo verzweigung moeglich ist */
if (i<=maxwo) indepnormal( graph, set, i, maxwo, abbruch );

/* dann doch einfuegen: */

if (!(*abbruch))
  {
 for (j=1; j<reg; j++)
     { if ((puffer=graph[wo][j])>wo) 
       { set.drin &= NMASK[puffer]; }
     }
  set.groesse++;
  rest=knotenzahl-graph[0][0];
            for (k=wo+1; k<=graph[0][0]; k++)
                    { if (set.drin & MASK[k]) rest++; }
  if (set.groesse > maxindep-five14s[rest]) *abbruch =1;
     else
       {
	 for (i=wo+1; (i<=maxwo) && !(set.drin &MASK[i]) ; i++);
	 /* sucht naechsten vertex, wo verzweigung moeglich ist */
	 if (i<=maxwo) indepnormal( graph, set, i, maxwo, abbruch );
       }
}/* ende if not abbruch */
} /* ende taille == 3 */
else /* d.h. taille >=4 -- hier koennen die optimierten werte benutzt werden */
{

/* zunaechst nicht einfuegen: */
for (i=wo+1; (i<=maxwo) && !(set.drin &MASK[i]); i++);
/* sucht naechsten vertex, wo verzweigung moeglich ist */
if (i<=maxwo) indepnormal( graph, set, i, maxwo, abbruch );

/* dann doch einfuegen: */

if (!(*abbruch))
  {
 for (j=1; j<reg; j++)
     { if ((puffer=graph[wo][j])>wo) 
       { set.drin &= NMASK[puffer]; }
     }
  set.groesse++;
  valf[3]=knotenzahl-graph[0][0]; valf[0]=valf[1]=valf[2]=0;
            for (k=wo+1; k<=graph[0][0]; k++)
                    { 
		      if (set.drin & MASK[k]) 
                             { for (z=reg-1;(z>=0) && (graph[k][z]>wo); z--); valf[reg-z-1]++;} 
		    }

  if (set.groesse > maxindep-minimalwerte[valf[1]][valf[2]][valf[3]]-valf[0])
         *abbruch =1;
     else
       {
for (i=wo+1; (i<=maxwo) && !(set.drin &MASK[i]) ; i++);
/* sucht naechsten vertex, wo verzweigung moeglich ist */
if (i<=maxwo) indepnormal( graph, set, i, maxwo, abbruch );
      }
}/* ende if not abbruch */
}/* ende taille >=4 */
}
/**************************INDEPENDENCE**********************************/

int independence(graph, vertex)

GRAPH graph;
unsigned char vertex;

/* Hier wird indeps[vertex] belegt. wird nur fuer vertex>1 aufgerufen */

{
unsigned char grenze, puffer, rest;
IND_SET *schreibezeiger, *lesezeiger, *ende;
int i,j,k, abbruch;
unsigned  long long int vertexmask, nextmask, stempel;
unsigned int *naechstezahl;
unsigned char valf[4];


if (vertex<=indepstoreborder)
{/*anfang <= indepstoreborder*/
if (taille==3)
{ /* anfang taille==3 */
vertexmask=MASK[vertex];
nextmask=MASK[vertex+1];
stempel= ~0ULL; /* die auszuschliessenden werte */
            for (j=1; j<reg; j++)
                { if ((puffer=graph[vertex][j])>vertex) 
                    { stempel ^= MASK[puffer];
                    /* Jeder kommt nur einmal vor -- deshalb tut es exor */
		    }
		}

grenze = knotenzahl-graph[0][0];
indeps[vertex] = (IND_SET *)malloc(numberofinds[vertex]*sizeof(IND_SET));
 if (indeps[vertex] ==NULL) { fprintf(stderr,"Can't get more space (1) -- exit\n"); exit(1); }
schreibezeiger=indeps[vertex];
naechstezahl = numberofinds + vertex + 1;
(*naechstezahl)=0;


for (i=1; i<vertex; i++)
  { lesezeiger=indeps[i];
    ende = lesezeiger + numberofinds[i];
    for ( ; lesezeiger < ende; lesezeiger++ )
      { 
        if ((lesezeiger->drin) & nextmask) (*naechstezahl)++;
        if ((lesezeiger->drin) & vertexmask) 
                                       /* d.h. kann dazugenommen werden */
	  { schreibezeiger->drin =lesezeiger->drin & stempel;
            schreibezeiger->groesse=lesezeiger->groesse +1;

            rest=grenze;
            for (k=vertex+1; k<=graph[0][0]; k++)
	      { if ((schreibezeiger->drin) & MASK[k]) rest++; }

            if (schreibezeiger->groesse > maxindep-five14s[rest]) return(0);

            if ((schreibezeiger->drin) & nextmask) 
                                                 (*naechstezahl)++;
            schreibezeiger++;
	  }
      } /* ende der inneren schleife -- d.h. fuer festes i alle indsets */
  } /* ende for alle i */



   }/* ende taille==3 */

else /* d.h. taille>=4 */
{
vertexmask=MASK[vertex];
nextmask=MASK[vertex+1];
stempel=~0ULL; /* die auszuschliessenden werte */
            for (j=1; j<reg; j++)
                { if ((puffer=graph[vertex][j])>vertex) 
                    { stempel ^= MASK[puffer];
                    /* Jeder kommt nur einmal vor -- deshalb tut es exor */
		    }
		}

grenze = knotenzahl-graph[0][0];
indeps[vertex] = (IND_SET *)malloc(numberofinds[vertex]*sizeof(IND_SET));
 if (indeps[vertex] ==NULL) { fprintf(stderr,"Can't get more space (2) -- exit\n"); exit(1); }
schreibezeiger=indeps[vertex];
naechstezahl = numberofinds + vertex + 1;
(*naechstezahl)=0;

for (i=1; i<vertex; i++)
  { lesezeiger=indeps[i];
    ende = lesezeiger + numberofinds[i];
    for ( ; lesezeiger < ende; lesezeiger++ )
      {
        if ((lesezeiger->drin) & nextmask) (*naechstezahl)++;
        if ((lesezeiger->drin) & vertexmask) 
                                       /* d.h. kann dazugenommen werden */
	  {
            schreibezeiger->drin = lesezeiger->drin & stempel;
            schreibezeiger->groesse=lesezeiger->groesse +1;

            valf[3]=grenze; valf[0]=valf[1]=valf[2]=0;
            for (k=vertex+1; k<=graph[0][0]; k++)
                    { if (schreibezeiger->drin  & MASK[k]) valf[reg-adj[k]]++; }
              if (schreibezeiger->groesse > 
                   maxindep-minimalwerte[valf[1]][valf[2]][valf[3]]-valf[0]) 
                                                                { if (!recover) return(0); }

            if ((schreibezeiger->drin) & nextmask) 
                                                 (*naechstezahl)++;
            schreibezeiger++;
	  }
      } /* ende der inneren schleife -- d.h. fuer festes i alle indsets */
  } /* ende for alle i */

   } /* ende taille>=4 */
   } /* ende vertex <=indepstoreborder */


else /* d.h. vertex > indepstoreborder */
{

if (recover) return(1); /* Falls die MINIMALWERTE datei vergroessert wurde,
                           wird eventuell schon bei recover abgebrochen, was
                           ueble Folgen haben kann. Deshalb bei recover nie
                           abbrechen */

abbruch=0;
for (i=1; (i<=indepstoreborder) && !abbruch; i++)
  { lesezeiger=indeps[i];
    ende = lesezeiger + numberofinds[i];
    for ( ; (lesezeiger < ende) && !abbruch; lesezeiger++ )
      { 
for (j=indepstoreborder+1; (j<=vertex) && !(lesezeiger->drin &MASK[j]); j++);
/* sucht naechsten vertex, wo verzweigung moeglich ist */
if (j<=vertex) { indepnormal( graph, *lesezeiger,(unsigned char)j, vertex, &abbruch );}
     }/* ende innere schleife -- d.h. fuer eine vorkonstruierte unabh. Menge */
  } /* ende aeussere schleife -- d.h. alle indepniveaus */


return(!abbruch);
} /* ende vertex > indepstoreborder */


return(1);
   } /* ende der funktion */


/******************************KONSTRMENGE******************************/

void konstrmenge(graph,vertex,max,tiefe,menge,vmg)

/* hier werden knoten mit abstand tiefe der Menge hinzugefuegt */

GRAPH graph;
int vertex, max, tiefe;
unsigned long long int *menge;
unsigned long long int *vmg; /* vergleichsmenge */

{  

*menge = *menge | MASK[graph[vertex][0]];
if ( (tiefe < max) && !(*menge & *vmg) )
  konstrmenge(graph,(int)graph[vertex][0],max,tiefe+1,menge,vmg);
	  

if (graph[vertex][1]!=leer)
{   
  *menge = *menge | MASK[graph[vertex][1]];
  if ( (tiefe < max) && !(*menge & *vmg)) 
    konstrmenge(graph,(int)graph[vertex][1],max,tiefe+1,menge,vmg);
}


if (graph[vertex][2]!=leer)
{
  *menge = *menge  | MASK[graph[vertex][2]];
  if ( (tiefe < max) && !(*menge & *vmg))
                 konstrmenge(graph,(int)graph[vertex][2],max,tiefe+1,menge,vmg);
}
}

/***************************SEARCHCYCLE****************************/

int searchcycle(GRAPH g,ADJAZENZ adj,ADJAZENZ mark, unsigned char vertex,
                int length, int max)
 
     /* length is the length of a cycle if there is an edge from
	vertex to ziel (that is the vertex marked 2). */

{

  int i,next;

  if (forbcyc[length])
    { for (i=0; i< adj[vertex]; i++)
      if (mark[g[vertex][i]]==2) return 1;
    }

  if (length<max)
      for (i=0; i< adj[vertex]; i++)
	{ 
	  next=g[vertex][i];
	  if (mark[next]==0)
	    { mark[next]=1;
	    if (searchcycle(g,adj,mark,next,length+1,max)) return 1;
	    mark[next]=0;
	    }
	}
  return 0;
}


/*******************************NOCYC******************************/

int nocyc(GRAPH g, ADJAZENZ adj, unsigned char start, unsigned char ziel)

/* checks whether an edge start to ziel would be part of a cycle of
forbidden length. The edge itself may not yet be in the graph. Returns
1 if it is OK to insert the edge (would not be on a forbidden cycle, 0
else. */

{
  int upto, i;
  ADJAZENZ mark;


  upto=forbiddencycles;
  if (upto>g[0][0])
    { for (upto=g[0][0]; (upto>=3) && (forbcyc[upto]==0); upto--);

    /* find largest possible size of a forbidden cycle */
    
    if (upto==2) return 1;
    }

  if (adj[start]>adj[ziel]) {i=start; start=ziel; ziel=i; }

  for (i=1; i<=g[0][0]; i++) mark[i]=0;
  mark[start]=1; mark[ziel]=2;

  for (i=0; i<adj[start]; i++)
    {
      mark[g[start][i]]=1;
    if (searchcycle(g,adj,mark,g[start][i],3,upto)) 
      { return 0; }
    mark[g[start][i]]=0;
    }

  return 1;
}




/*******************************ABSTANDOK******************************/

int abstandok(graph,v1,v2,taille)

GRAPH graph;
int v1, v2;
unsigned char taille;

/* ueberprueft, ob der Abstand von v1 und v2 in graph <= taille-1 ist */
/* das geschieht, indem die Menge m1 aller Knoten, die einen Abstand */
/* von (taille-2)/2 von v1 haben  ( und analog m2 bezueglich v2 ) konstruiert*/
/* werden und der schnitt daraufhin ueberprueft wird, ob er leer ist */
/* der wert von taille muss aber mindestens 3 sein -- wie das in einfachen */
/* Graphen natuerlich ohnehin der Fall ist */

{

unsigned long long int m1, m2;

if (taille==3) return(1);

m1=MASK[v1];

m2=MASK[v2];

konstrmenge(graph,v1,((taille-1)/2),1,&m1,&m2);


if (!(m1 & m2)) konstrmenge(graph,v2,((taille-2)/2),1,&m2,&m1);

 if (m1 & m2) return 0; else return 1;

}



/**********************************CODIERE********************************/

void codiere(graph,code)
GRAPH graph;
unsigned char code[codelaenge];

{ int i,j,k;

i=0;

for (j=2; j<= graph[0][0]; j++) for (k=1; k<reg; k++)
        if (graph[j][k]>j) { code[i]=graph[j][k]; i++; }
}

/*****************************SHORTCODIERE********************************/

/* codiert laenger -- aber fuer shortcode -- mit 1-2-3 */

void shortcodiere(graph,code)
GRAPH graph;
unsigned char code[codelaenge];

{ int i,j,k;

i=0;

for (j=2; j<= graph[0][0]; j++) for (k=1; k<reg; k++)
        if (graph[j][k]>j) { code[i]=graph[j][k]; i++; }
}


/********************************WEGSPEICHERN******************************/

void wegspeichern_minibaum(liste,graphzahlen,graph)


unsigned char * liste[knoten/2-1];
unsigned long long int graphzahlen[knoten];
GRAPH graph;
{

    return;
    //added return, so this method basically does nothing
    //this is because this is not needed for pregraphs

FILE *fil;
int zahl,i;
struct tm *zeitzeiger;
struct tms TMS;
time_t b;
unsigned int savetime;

/* Damit keine Unfaelle beim sichern entstehen: Ab 10 sec vor erwartetem kill*/
/* nicht mehr sichern  */
if (TIMERESTRICTED){ b=time(0);
                 zeitzeiger=localtime(&b);
                 i= (3600*zeitzeiger->tm_hour) + (60*zeitzeiger->tm_min)
                    + zeitzeiger->tm_sec;  }
         if ( !(TIMERESTRICTED) || (i<=endtime-10) ||
              (zeitzeiger->tm_wday == 0) || (zeitzeiger->tm_wday == 6) ||
              ( (endtime<starttime) && (i>=starttime) ) )
{
if (!graph6_output_minibaum && !noout_minibaum)
  {
    if (!singleout_minibaum)
      for (zahl=startzahl; zahl<knotenzahl; zahl+=2) 
	if (graphzahlen[zahl])
	  {
	    codename[6]=(zahl/10)+48; codename[7]=(zahl%10)+48;
	    fil = fopen(codename,"a");
	    if (fwrite(liste[zahl],sizeof(unsigned char),
		       nextposition[zahl]-liste[zahl],fil)
		< nextposition[zahl]-liste[zahl]) 
	      { fprintf(stderr,"Dangerous error ! Could not write !!\n\n");
		exit(1); }
	    fclose(fil);
	    nextposition[zahl]=liste[zahl];
	    graphzahlen[zahl]=0ULL;
	  }
    
    if (graphzahlen[knotenzahl])
      { if ( lastout ) fil=stdout; else 
	  { codename[6]=(knotenzahl/10)+48; codename[7]=(knotenzahl%10)+48; 
	    fil = fopen(codename,"a");  }
	if (fwrite(liste[knotenzahl],sizeof(unsigned char),
		   nextposition[knotenzahl]-liste[knotenzahl],fil)
	    < nextposition[knotenzahl]-liste[knotenzahl])
	  { fprintf(stderr,"Dangerous error ! Could not write !!\n\n");
	    exit(1); }
	if (!lastout) { fclose(fil); }
	number_of_graphs += graphzahlen[knotenzahl];
	graphzahlen[knotenzahl]=0ULL;
	nextposition[knotenzahl]=globalliste[knotenzahl];
      }
  } /* ende !graph6_output && !noout*/
 else if (noout_minibaum) { number_of_graphs += graphzahlen[knotenzahl]; graphzahlen[knotenzahl]=0ULL; }

if (S_intervall)
      {
      fil = fopen(lastgraphname,"w");
      if (fwrite(graph,1,sizeof(GRAPH),fil)<sizeof(GRAPH))
                  { fprintf(stderr,"Dangerous error ! Could not write !!\n\n");
		    exit(1); }

      times(&TMS);
      savetime=oldtime + (unsigned int) TMS.tms_utime;
      unsignedintegerschreib(savetime,fil);
      unsigned_ll_integerschreib(number_of_graphs,fil);
      if (modulo_minibaum) unsignedintegerschreib(splitzaehler,fil);
      fclose(fil);
      }
}
else { fprintf(stderr,"No Trying to save when killing expected.\n");
       exit(1);
     }
}

/*****************************WRITEG6************************************/

void writeg6(g)
GRAPH g;
{
int i,j,k,org,nv;
int nlen,bodylen;
static unsigned char g6bit[] = {32,16,8,4,2,1};
unsigned char code[20+knoten*(knoten-1)/12];
register unsigned char *body;

nv = g[0][0];
if (nv==knotenzahl) number_of_graphs++;

if (nv <= 62)
{
    code[0] = 63 + nv;
    nlen = 1;
}
else
{
    code[0] = 63 + 63;
    code[1] = 63 + 0;
    code[2] = 63 + (nv >> 6);
    code[3] = 63 + (nv & 0x3F);
    nlen = 4;
}

body = code + nlen;
bodylen = ((nv * (nv-1)) / 2 + 5) / 6;
for (i = 0; i < bodylen; ++i)
    body[i] = 63;
body[bodylen] = '\n';

for (i = org = 0; i < nv; org += i, ++i)
{
    j = g[i+1][0]-1;
    if (j < i)
    {
	k = org + j;
	body[k/6] += g6bit[k%6];
    }
    j = g[i+1][1]-1;
    if (j < i)  
    {
        k = org + j;
        body[k/6] += g6bit[k%6]; 
    }
    j = g[i+1][2]-1;
    if (j < i)  
    {
        k = org + j;
        body[k/6] += g6bit[k%6]; 
    }
}

j = nlen + bodylen + 1;
if (fwrite(code,(size_t)1,(size_t)j,stdout) != j)
{
    fprintf(stderr,"fwrite() failed in writeg6()\n");
    exit(1);
}

}


/****************************SONDERTEST************************/

int sondertest(GRAPH graph)

/* Ueberprueft, ob ein Graph zweifaerbbar ist -- wird nur aufgerufen, wenn maxvalenz
   zwei ist -- d.h. hier natuerlich zwei-regulaer */

{

int marks[knoten+1], knotenzahl;
int i, merke, laenge, dummy, lauf;

knotenzahl=graph[0][0];
for (i=1; i<=knotenzahl; i++) marks[i]=0;

for (i=1; i<=knotenzahl; i++)
  if (!marks[i])
    { merke=graph[i][0]; lauf=graph[i][1]; dummy=i; laenge=2;
      marks[merke]=marks[lauf]=1;
      while (lauf != merke) 
	{ if (graph[lauf][0]==dummy) /* Kante zurueck */
	    { dummy=lauf; lauf=graph[lauf][1]; }
	  else 
	    { dummy=lauf; lauf=graph[lauf][0]; }
	  laenge++;
	  marks[lauf]=1;
	}
      if (laenge%2) return(0); /* nicht (2-) faerbbar */
    }
return(1);
}



/*********************************WEITERMATCHEN***********************/

int weitermatchen(GRAPH graph, int coloured[],
		  int matching[][2],int matchingsize, int next)

/* versucht das matching auszubauen. Wenn es maximal ist und alle "forced"
   Knoten saturiert sind, wird die Entscheidung darauf reduziert, ob der
   Graph ohne dieses Matching faerbbar ist. Gibt 1 zurueck, wenn das bisherige
   Matching zu einer Faerbung fortsetzbar ist, 0 sonst */


{
int i,j,start,ziel, knotenzahl;
GRAPH neug;
unsigned char *run, *end;


knotenzahl=graph[0][0];
for ( ; (next<=knotenzahl) && coloured[next]; next++);

if (next>knotenzahl) /* Neuen Graphen fertigmachen und den auf Faerbbarkeit testen */
  { memcpy(neug,graph,reg*(knotenzahl+1)); 
    
    for (i=0; i<matchingsize; i++)
      { start=matching[i][0]; ziel=matching[i][1]; 
	for (j=0; neug[start][j] != ziel; j++); 
	for ( ; j<reg; j++) neug[start][j]=neug[start][j+1];
	for (j=0; neug[ziel][j] != start; j++); 
	for ( ; j<reg; j++) neug[ziel][j]=neug[ziel][j+1];
      }
      return(sondertest(neug));
  }


for (run=graph[next], end=run+reg; run != end; run++)
  if (!coloured[*run])
    { ziel=*run;
      coloured[next]=coloured[ziel]=1;
      matching[matchingsize][0]=next; matching[matchingsize][1]=ziel;
      if (weitermatchen(graph, coloured, matching, matchingsize+1, 
			next+1))
	{ coloured[next]=coloured[ziel]=0; return(1); }
      coloured[next]=coloured[ziel]=0;
    }

return(0);
}




/**************************FAERBBAR****************************/

int faerbbar(GRAPH graph)

/* stark geaendert fuer 3-regulaere Graphen, aber stammt aus multi_class2.c */

/* Ein Graph ist genau dann faerbbar, wenn man ein enthaltensmaximales Matching
   herausnehmen kann und der resultierende graph ist auch faerbbar (mit einer
   Farbe weniger). Hier wird versucht, ein solches Matching zu konstruieren.
   Dabei muessen natuerlich alle Knoten mit maximalem Grad im matching enthalten
   sein. Eine erste Kante kann sogar beliebig gewaehlt werden. Das wird hier eine 
   Kante zwischen zwei Knoten moeglichst kleiner Valenz um die Anforderung, dass
   alle Knoten maximaler Valenz saturiert werden muessen staerker zu gestalten. 
   gibt 1 zurueck, wenn er faerbbar ist, 0 sonst */


{
int i, knotenzahl;
int coloured[knoten+1]; /* ist der Knoten schon zu einer gefaerbten Kante adjazent ? */
int matching[knoten][2]; /* hier werden die Kanten gespeichert */
int matchingsize;


knotenzahl=graph[0][0];

for (i=1 ; i<=knotenzahl ; i++) coloured[i]=0;

matching[0][0]=knotenzahl; matching[0][1]=graph[knotenzahl][0];

/* jetzt hat man sich fuer die Kante entschieden, die man festlegen will. */


coloured[matching[0][0]]=coloured[matching[0][1]]=1;
matchingsize=1;

return (weitermatchen(graph,coloured,matching,matchingsize,1));

}






/**********************************AUFSCHREIBEN*****************************/


/*** Nimmt den Graphen und schreibt ihn vorerst in die globale Liste */
/* wenn die voll ist, wird sie weggesichert. */
/* bei graph6 output wird writeg6 direkt aufgerufen */
/* Overwritten for pregraphs: pass graph on to pregraphs */
void aufschreiben_minibaum(GRAPH g) {
    handle_minibaum_result(g, g[0][0]);
}

/****************************ENDSUCHE**********************************/

void endsuche(zweig,blaetterkette,zaehler)

BRANCH *zweig;
BRANCH **blaetterkette;
int *zaehler; /* das wievielte element der kette muss beschrieben werden */


{

if ( zweig->blatt ) 
   { blaetterkette[*zaehler]=zweig; (*zaehler)++; }

else {
     if ( (zweig->links != nil) && (!zweig->totlinks) )
                                  endsuche(zweig->links,blaetterkette,zaehler);
     if ( (zweig->rechts != nil) && (!zweig->totrechts) )
                                 endsuche(zweig->rechts,blaetterkette,zaehler);
     }
}



/**************************BAUMANALYSE**************************************/


void baumanalyse(bis_wo)

unsigned char bis_wo; /* gibt an, bis zu welcher zahl eine analyse noetig
                         ist. Im allgemeinen wird graph[0][0] uebergeben */

{
int i,j,zaehler;


for (i=1; i<=bis_wo; i++)
    { if ((roots[i].tiefe) && (!(roots[i].wurzeltot)))
	  { zaehler=0;

          blaetterkette[i]=(BRANCH **)malloc(sizeof(BRANCH *)*blaetterzahl[i]);
 if (blaetterkette[i] ==NULL) { fprintf(stderr,"Can't get more space (3) -- exit\n"); exit(1); }
            for (j=0; j<6; j++)
                if (!(roots[i].tot[j])) 
                     endsuche(roots[i].naechster[j],blaetterkette[i],&zaehler);

	}
     }
}

/**************************ZWEIGENTFERNEN*****************************/

void zweigentfernen(zweig)

BRANCH *zweig;

{ if (zweig->blatt) blaetterzahl[zweig->urbild[1]]--;
    else {
          if (zweig->links != nil) zweigentfernen(zweig->links);
          if (zweig->rechts != nil) zweigentfernen(zweig->rechts);
         }

  free(zweig);
}
 


/****************************ROOTENTFERNEN********************************/

void rootentfernen(i)

unsigned char i;

/* entfernt die i-te wurzel (genauer: entfernt alles, was daran haengt)
   und initialisiert die werte der wurzel */

{
int j;


roots[i].wurzeltot=roots[i].tiefe=0;

for (j=0; j<6; j++) 
      {
      zweigentfernen( roots[i].naechster[j] );
      roots[i].naechster[j]=nil;
      roots[i].tot[j]=0;
      }
blaetterzahl[i]=0;
}



/**********************ZWEIGAUFRAEUMEN**********************************/


BRANCH *zweigaufraeumen(zweig,k)

BRANCH *zweig;
unsigned char k;

/* ich komme hier nur hin, wenn diese abzweigung nicht schon in einem
   schritt vor k totgesetzt wurde */

{ 

if (zweig->tiefe == k) { 
                       if (zweig->blatt) blaetterzahl[(zweig->urbild)[1]]--;
                       else {
                       if (zweig->links != nil) zweigentfernen(zweig->links);
                       if (zweig->rechts != nil) zweigentfernen(zweig->rechts);
		            }
                       free(zweig); return(nil); }
else /* d.h. es ist ein aelterer zweig */
{

if (zweig->totlinks == k) zweig->totlinks=0;
if (zweig->totrechts == k) zweig->totrechts=0;

if ((zweig->links != nil) && (zweig->totlinks==0) )
zweig->links=zweigaufraeumen(zweig->links,k);

if ((zweig->rechts != nil) && (zweig->totrechts==0) )
zweig->rechts=zweigaufraeumen(zweig->rechts,k);



          if ( (zweig->links == nil) && (zweig->rechts == nil) 
               && (zweig->totlinks == 0) && (zweig->totrechts == 0)
                                     && (!(zweig->blatt)) )
              {  zweig->blatt = 1; blaetterzahl[zweig->urbild[1]]++; }
          return(zweig);
}/* ende else */
}




/****************************BAUMAUFRAEUMEN******************************/


void baumaufraeumen(k,vertexzahl)

unsigned char k;
unsigned char vertexzahl;

/* entfernt alle im schritt k gesetzten Blaetter und totmarken 
   ueberprueft dabei alle an den roots 1..vertexzahl haengenden Baeume .
   blaetter und totmarken mit einem index > k existieren zu diesem Zeitpunkt
   schon nicht mehr im baum 
   Hier wird auch das forbidden array wieder "aufgeraeumt" */

{

unsigned char i,j;


for (i=1; i<=vertexzahl; i++)
    { if (roots[i].tiefe)
         {
            if (roots[i].tiefe == k) rootentfernen(i);
              else /* d.h.: es ist eine schon aeltere wurzel */
                  {  if (roots[i].wurzeltot==k) roots[i].wurzeltot=0;
                     for (j=0; j<6; j++)
                        { if (roots[i].tot[j]==k) roots[i].tot[j]=0;
                          if (!roots[i].tot[j])
                              zweigaufraeumen((roots[i].naechster)[j],k);
/* wenn dieser ast schon vorher tot war, wurde er auch im k-ten schritt
   nicht mehr benutzt */
                        }
                  } /* ende aeltere wurzel */
          }
      }

/* jetzt noch forbidden aufraeumen: */

for (i=2; i<=vertexzahl; i++) for (j=i+1; j<=vertexzahl; j++)
  { if (forbidden[i][j]==k)
        forbidden[i][j]=0; }


} /* ende baumaufraeumen */

/*********************************ordnen******************************/

void ordnen(tripel)
/* ordnet einfach nur die naechsten 3 char nach dem Zeiger tripel */

unsigned char tripel[reg];

{

unsigned char puffer;

if (tripel[0]<tripel[1])
      { if (tripel[2]<tripel[0]) { puffer=tripel[0]; tripel[0]=tripel[2];
                                    tripel[2]=tripel[1]; tripel[1]=puffer; }
        else /* d.h. tripel[2]>=tripel[0]<tripel[1] */
              { if (tripel[2]<tripel[1]) { puffer=tripel[2]; 
                                   tripel[2]=tripel[1]; tripel[1]=puffer; }
               }
      }

else /* d.h. tripel[0]>=tripel[1] */
    { if (tripel[2]<= tripel[1]) { puffer=tripel[0]; tripel[0]=tripel[2];
                                   tripel[2]=puffer; }
       else /* d.h. tripel[0]>=tripel[1]<tripel[2] */ 
             { if (tripel[0]<=tripel[2]) { puffer=tripel[0]; 
                                       tripel[0]=tripel[1]; tripel[1]=puffer; }
                  else /* d.h. tripel[0]>tripel[2]>tripel[1] */
                      { puffer=tripel[0]; tripel[0]=tripel[1];
                        tripel[1]=tripel[2]; tripel[2]=puffer; }
             }
     }

}

/****************************COMPARE************************************/

int compare (bild,tr,br)

unsigned char bild[knoten+1];  /* die Bilder der Knoten */
unsigned char tr[reg];         /* die Reihe in der die Bilder eingesetzt */
                               /* werden muessen                         */
unsigned char br[reg];         /* die Reihe mit der die ersetzte verglichen */
                               /* wird */

/* vergleicht die Reihe tauschreihe, in der alle Elemente durch ihre Bilder */
/* ersetzt werden mit bestreihe. Gibt etwas negatives zur"uck wenn tr<br, 0 */
/* wenn tr=br und etwas positives sonst */

{

unsigned char test[reg];

test[0]=bild[tr[0]];
if (tr[1] != leer) test[1]=bild[tr[1]]; else test[1]=leer;
if (tr[2] != leer) test[2]=bild[tr[2]]; else test[2]=leer;
ordnen(test);
if (test[0] != br[0]) return(test[0]-br[0]);
   else if (test[1] != br[1]) return(test[1]-br[1]);
           else return(test[2]-br[2]);

}



/***********************NORMALWEITER********************************/


void normalweiter(graph,adj,bild,urbild,wo,naechstenummer,abbruch)

GRAPH graph;
ADJAZENZ adj;
BILD bild;
BILD urbild;
unsigned char wo; /* wahrer namen des hier zu behandelnden knotens */
unsigned char naechstenummer; /* naechste zu vergebende nummer */
int *abbruch;

{
int i,luecken,test;
int lueckenzeiger[2];

if (vertex_trans)
   if ((adj[wo]<reg) || (bild[wo]==graph[0][0])) transtest=1;

if (adj[wo]>0)
{

luecken = 0;
for (i=0; i<adj[wo]; i++) if (bild[graph[wo][i]]==unbelegt)
                          { lueckenzeiger[luecken]=graph[wo][i]; luecken++; }



switch(luecken)
   {
   case 0: { test = compare(bild,graph[wo],graph[bild[wo]]);
             if (test<0) *abbruch=1;
                 else if ((test==0) && (bild[wo]<graph[0][0]))
                   {
                   normalweiter(graph,adj,bild,urbild,
                                urbild[bild[wo]+1],naechstenummer,abbruch);
                   }
              break;
             } /* ende des falles "keine luecken" */

     case 1: { bild[lueckenzeiger[0]]=naechstenummer;
               urbild[naechstenummer]=lueckenzeiger[0];
               test = compare(bild,graph[wo],graph[bild[wo]]);
               if (test<0)
                   *abbruch=1;
                         else if ((test==0) && (bild[wo]<graph[0][0]))
                           {
                   normalweiter(graph,adj,bild,urbild,urbild[bild[wo]+1],
                                                 naechstenummer+1,abbruch);
                   }
              bild[lueckenzeiger[0]]=unbelegt;
              break;
             } /* ende des falles "genau eine luecke" */

                     
     case 2: { bild[lueckenzeiger[0]]=naechstenummer;
               bild[lueckenzeiger[1]]=naechstenummer+1;
               urbild[naechstenummer]=lueckenzeiger[0];
               urbild[naechstenummer+1]=lueckenzeiger[1];
               test = compare(bild,graph[wo],graph[bild[wo]]);
               if (test<0)
                  { *abbruch=1;
                  }
                      else if ((test==0) && (bild[wo]<graph[0][0]))
                           {
                   normalweiter(graph,adj,bild,urbild,urbild[bild[wo]+1],
                                                 naechstenummer+2,abbruch);

                   if (!(*abbruch)) /* dann muss das zweite blatt auch noch
                                       angefuegt werden */
                   {
                   bild[lueckenzeiger[0]]=naechstenummer+1;
                   bild[lueckenzeiger[1]]=naechstenummer;
                   urbild[naechstenummer]=lueckenzeiger[1];
                   urbild[naechstenummer+1]=lueckenzeiger[0];
                   normalweiter(graph,adj,bild,urbild,urbild[bild[wo]+1],
                                                   naechstenummer+2,abbruch);
                   }
                   }
              bild[lueckenzeiger[0]]=bild[lueckenzeiger[1]]=unbelegt;
              break;
             } /* ende des falles "zwei luecken" */
        } /* ende switch */
}/* ende adj[wo]>0 denn wenn die adjazenz 0 ist, muss nichts gemacht werden */
} /* ende normalweiter */


/******************************TESTBERECHNUNG***********************/

int testberechnung(graph,adj,wo,bild,urbild,naechster)

GRAPH graph;
ADJAZENZ adj;
unsigned char wo; /* der wirkliche name des zu behandelnden knotens */
BILD bild;
BILD urbild;
unsigned char *naechster; /* die adresse fuer den naechsten knoten, an dem
                             wahlmoeglichkeit besteht */

/* ist der letzte knoten an dem noch wahlmoeglichkeit besteht noch nicht
   voll saturiert, so wird nur ein negatives testergebnis mit einbezogen
   in die auswertung von test -- ein positives kann sich noch aendern    */


{
int test,vortest,ur,stufe,ende,j;


ur=graph[0][0]+1; /* das soll der defaultwert fuer *naechster sein, wenn der
                    graph schon komplett ist */
vortest=ende=0;
test = compare(bild,graph[wo],graph[bild[wo]]);

for (stufe=bild[wo]+1; (!ende) && (test==0) && (stufe<=graph[0][0]); )
 {
ur = urbild[stufe];

for (j=0; j<adj[ur]; j++) if (bild[graph[ur][j]]==unbelegt) { ende=1;
                                                              j=reg; }



if (!ende) vortest = compare(bild,graph[ur],graph[stufe]);
if (adj[ur]<reg) { ende=1; if (vortest<0) test=vortest; }
                else test=vortest;
stufe++; 
}
*naechster=ur;


return(test);
}

/*****************************WURZELTOTSETZEN**************************/

void wurzeltotsetzen(wurzel,zweig,tiefe)

ROOT *wurzel;
BRANCH *zweig;
unsigned char tiefe;

/* setzt die entsprechende Spalte in tot auf tiefe und ueberprueft, ob es
   an irgendeiner stelle noch weitergeht. Wenn nicht wird wurzeltot gesetzt.
   Achtung die reihenfolge in der schleife ist wichtig */

{
int test, i;

test=tiefe;
for (i=0; i<6; i++)
    { if ((wurzel->naechster)[i]==zweig) (wurzel->tot)[i]=tiefe;
      if ((wurzel->tot)[i]==0) test=0; }
wurzel->wurzeltot = test;
}



/**************************TOTSETZEN**********************************/

void totsetzen (zweig,woher,tiefe)

BRANCH *zweig;  /* welcher zweig soll bearbeitet werden */
BRANCH *woher;  /* woher kommt man ? */
unsigned char tiefe;  /* das zeichen, ab wann ein zweig tot ist */

/* diese funktion laeuft den baum rueckwaerts und setzt den totzeiger in
   die richtung aus der man kommt auf tiefe. Wenn der andere zeiger auch
   blockiert ist, geht es noch einen zweig zurueck, ...                */

{

if ( zweig->links==woher ) 
                   /* der Fall man kommt aus der ersten Richtung */
   { zweig->totlinks = tiefe; /* der war vorher zwangsweise noch 0 */
     if (zweig->totrechts) 
                 { if (zweig->vorgaenger != nil)
                                  totsetzen(zweig->vorgaenger,zweig,tiefe);
                      else wurzeltotsetzen(zweig->zurwurzel,zweig,tiefe);
                  }
   }

else               /* der Fall man kommt aus der zweiten Richtung */
   { zweig->totrechts = tiefe; /* der war vorher zwangsweise noch 0 */
     if (zweig->totlinks)
                 { if (zweig->vorgaenger != nil)
                                  totsetzen(zweig->vorgaenger,zweig,tiefe);
                      else wurzeltotsetzen(zweig->zurwurzel,zweig,tiefe);
                  }

   }

}


/*************************ZWEIGINITIALISIEREN***************************/

void zweiginitialisieren(graph,adj,zweig,sollvorgaenger,sollbild,sollurbild,
                                            sollnextbild,solltiefe,sollwo)

GRAPH graph;
ADJAZENZ adj;
BRANCH *zweig; /* der zu initialisierende zweig */
BRANCH *sollvorgaenger;
BILD sollbild;
BILD sollurbild;
unsigned char sollnextbild;
unsigned char solltiefe;
unsigned char sollwo; 

{

zweig->blatt = 1;
memcpy(zweig->bild,sollbild,sizeof(BILD));
memcpy(zweig->urbild,sollurbild,sizeof(BILD));
zweig->totlinks=zweig->totrechts=0;
zweig->links=zweig->rechts=nil;
zweig->vorgaenger=sollvorgaenger;
zweig->nextbild=sollnextbild;
zweig->tiefe = solltiefe;
/* achtung! zurwurzel enthaelt keinen sinnvollen wert ! */

zweig->wo = sollwo;

}/* ende zweiginitialisieren */

/**************************MEMOCOMP*************************************/


int memocomp(adr1,adr2, anzahl)
/* Der Ersatz fuer das auf einigen dec-rechnern nicht funktionierende memcmp */

unsigned char *adr1, *adr2;
int anzahl;

{
for ( ;(*adr1 == *adr2) && (anzahl>1); adr1++, adr2++, anzahl--);
return (int)(*adr1 - *adr2);
}


/**************************forbiddenbelegen**********************/

void forbiddenbelegen( zweig, graph, adj )

BRANCH *zweig;
GRAPH graph;
ADJAZENZ adj;


{
int i;
unsigned char testpuffer[3], wo;
BILDELEMENT *bild, *urbild;


wo=zweig->wo;
bild=zweig->bild;
urbild=zweig->urbild;

for (i=0; i<adj[wo]; i++) /* maximal eine (!) luecke moeglich */
   { if (bild[graph[wo][i]] == unbelegt) testpuffer[i]=zweig->nextbild;
                          else testpuffer[i]=bild[graph[wo][i]]; }
for( ; i<reg; i++) testpuffer[i]=leer;
ordnen(testpuffer);

if (memocomp(testpuffer,graph[bild[wo]],adj[wo])==0)
   { 
    for (i=1; i<graph[bild[wo]][adj[wo]]; i++)
       if ((urbild[i] <= graph[0][0]) && (bild[urbild[i]]==i) 
            && !(forbidden[wo][urbild[i]])
            && (adj[urbild[i]]<reg) )
/* && (wo != urbild[i]) && ( graph[wo][0] != urbild[i])
&& ( graph[wo][2] != urbild[i]) && ( graph[wo][1] != urbild[i])  */ 
         { 
           if (wo>urbild[i]) forbidden[urbild[i]][wo]=graph[0][1];
           else forbidden[wo][urbild[i]]=graph[0][1];
         }
   }
}



/***********************************BAUMBAU*****************************/


void baumbau(graph,adj,zweig,abbruch)

GRAPH graph;
ADJAZENZ adj;
BRANCH *zweig;
int *abbruch;

/* baumbau versucht, startend bei zweig den baum auszubauen -- 
   beziehungsweise in den normalmodus umzuschalten, falls noetig */

{
BILDELEMENT *bild, *urbild;
int i, luecken, test;
int lueckenzeiger[2];
unsigned char wo, naechster;


wo = zweig->wo;


if (adj[wo]<reg) 
      { 
       if (adj[wo]>0) 
       {  bild=zweig->bild;
          if (adj[bild[wo]]>adj[wo]) 
            forbiddenbelegen( zweig, graph, adj ); 
            normalweiter(graph,adj,zweig->bild,zweig->urbild,wo,
                                          zweig->nextbild,abbruch); }
       }
else
   {

   zweig->blatt=0; /* es ist auf jeden Fall kein blatt mehr */
   bild = zweig->bild;
   urbild = zweig->urbild;
   blaetterzahl[urbild[1]]--; /* wenn neue blaetter eingefuegt werden muss
                                 die zahl wieder erhoeht werden */

/* zuerst wird herausgefunden, wieviel zu diesem knoten adjacente knoten noch
   kein bild haben : */

   luecken =0;
   for (i=0; i<reg; i++)
       if (bild[graph[wo][i]]==unbelegt)
          { lueckenzeiger[luecken]=graph[wo][i]; luecken++; }


/* jetzt wird je nach der anzahl der luecken unterschiedlich verfahren: */

switch(luecken)
   {
   case 0: { test = testberechnung(graph,adj,wo,bild,urbild,&naechster);
             if (test<0)  *abbruch=1; 
                else
                    if (test>0) {
                               zweig->totlinks=zweig->totrechts=graph[0][1];
                               if (zweig->vorgaenger != nil)
                               totsetzen(zweig->vorgaenger,zweig,graph[0][1]);
                      else wurzeltotsetzen(zweig->zurwurzel,zweig,graph[0][1]);
                                }
                 else /*d.h. test==0 */ if (bild[wo]<graph[0][0])
                   {
                  zweig->totrechts=graph[0][1];/*es gibt nur einen nachfolger*/
                   zweig->links=(BRANCH *)malloc(sizeof(BRANCH));
 if (zweig->links ==NULL) { fprintf(stderr,"Can't get more space (4)-- exit\n"); exit(1); }
                   zweiginitialisieren(graph,adj,zweig->links,zweig,
                          bild,urbild,zweig->nextbild,graph[0][1],naechster);
                   blaetterzahl[urbild[1]]++;
                   baumbau(graph,adj,zweig->links,abbruch);
                   }
              break;
             } /* ende des falles "keine luecken" */

     case 1: { bild[lueckenzeiger[0]]=zweig->nextbild;
               urbild[zweig->nextbild]=lueckenzeiger[0];
                     /* achtung ! jetzt hat das bild auf diesem zweig einen
                        wert, der wieder geaendert werden muss */
               test = testberechnung(graph,adj,wo,bild,urbild,&naechster);
               if (test<0)
                  { *abbruch=1; 
                  }
                      else
                      if (test>0)
                         { 
                               zweig->totlinks=zweig->totrechts=graph[0][1];
                               if (zweig->vorgaenger != nil)
                               totsetzen(zweig->vorgaenger,zweig,graph[0][1]);
                      else wurzeltotsetzen(zweig->zurwurzel,zweig,graph[0][1]);

                         }
                         else /* d.h. test == 0 */  if (bild[wo]<graph[0][0])
                           {
                  zweig->totrechts=graph[0][1];/*es gibt nur einen nachfolger*/
                  zweig->links=(BRANCH *)malloc(sizeof(BRANCH));
 if (zweig->links ==NULL) { fprintf(stderr,"Can't get more space (5)-- exit\n"); exit(1); }
                  zweiginitialisieren(graph,adj,zweig->links,zweig,
                       bild,urbild,(zweig->nextbild)+1,graph[0][1],naechster);
                  blaetterzahl[urbild[1]]++;
                  baumbau(graph,adj,zweig->links,abbruch);
                   }
              bild[lueckenzeiger[0]]=unbelegt;
              break;
             } /* ende des falles "genau eine luecke" */

                     
     case 2: { bild[lueckenzeiger[0]]=zweig->nextbild;
               bild[lueckenzeiger[1]]=(zweig->nextbild)+1;
               urbild[zweig->nextbild]=lueckenzeiger[0];
               urbild[(zweig->nextbild)+1]=lueckenzeiger[1];
                     /* achtung ! jetzt hat das bild auf diesem zweig einen
                        wert, der wieder geaendert werden muss */
               test = testberechnung(graph,adj,wo,bild,urbild,&naechster);
               if (test<0) *abbruch=1;
                      else
                      if (test>0)
                         { zweig->totlinks=graph[0][1];
                         }
                         else /* d.h. test == 0 */  if (bild[wo]<graph[0][0])
                           {
                   zweig->links=(BRANCH *)malloc(sizeof(BRANCH));
 if (zweig->links ==NULL) { fprintf(stderr,"Can't get more space (6)-- exit\n"); exit(1); }
                   zweiginitialisieren(graph,adj,zweig->links,zweig,
                        bild,urbild,(zweig->nextbild)+2,graph[0][1],naechster);
                   blaetterzahl[urbild[1]]++;
                   baumbau(graph,adj,zweig->links,abbruch);
	                   }
                   if (!(*abbruch)) /* dann muss das zweite blatt auch noch
                                       angefuegt werden */
	      { 
              bild[lueckenzeiger[0]]=zweig->nextbild+1;
              bild[lueckenzeiger[1]]=(zweig->nextbild);
              urbild[zweig->nextbild]=lueckenzeiger[1];
              urbild[(zweig->nextbild)+1]=lueckenzeiger[0];
              test = testberechnung(graph,adj,wo,bild,urbild,&naechster);
/* da testberechnung hier nicht nur die eine stelle ueberprueft, muss auch
im zweiten fall eine erneute berechnung stattfinden -- anders als bei
 test=compare(), das nur die eine stelle ueberprueft */
               if (test<0) { *abbruch=1; 
			  }
                      else
                      if (test>0)
                         { zweig->totrechts=graph[0][1];
                            if (zweig->totlinks)
                              { if (zweig->vorgaenger != nil)
                                totsetzen(zweig->vorgaenger,zweig,graph[0][1]);
                      else wurzeltotsetzen(zweig->zurwurzel,zweig,graph[0][1]);
			      }
                         }
                   else /*d.h. test == 0 */
                   {
                   zweig->rechts=(BRANCH *)malloc(sizeof(BRANCH));
 if (zweig->rechts ==NULL) { fprintf(stderr,"Can't get more space (7)-- exit\n"); exit(1); }
                   zweiginitialisieren(graph,adj,zweig->rechts,zweig,
                       bild,urbild,(zweig->nextbild)+2,graph[0][1],naechster);
                   blaetterzahl[urbild[1]]++;
                   baumbau(graph,adj,zweig->rechts,abbruch);
                   }
	        }/* ende if not abbruch zweiter fall */
              bild[lueckenzeiger[0]]=bild[lueckenzeiger[1]]=unbelegt;
              break;
             } /* ende des falles "zwei luecken" */
        } /* ende switch */
   } /* ende else -- d.h. adjazenz == reg */


}/* ende baumbau */



/*****************************faerbe**********************************/


void faerbe(colour,graph,vertex,abbruchzeiger)


signed char colour[knoten+1];
GRAPH graph;
int vertex;
int *abbruchzeiger;

{
int i;

for (i=0; (graph[vertex][i]!=leer) && (i<reg); i++)
               if (!(*abbruchzeiger))
                  { if (colour[graph[vertex][i]]==0)
                         { if (colour[vertex]==1) colour[graph[vertex][i]]=2;
                                   else colour[graph[vertex][i]]=1;
                        faerbe(colour,graph,graph[vertex][i],abbruchzeiger); }
                       else { if (colour[graph[vertex][i]]==colour[vertex])
                                      {*abbruchzeiger=1; return; } }
                  }
               else return;
}

/*********************************BIPARTIT******************************/

int bipartit(graph)

GRAPH graph;

{
int i,abbruch;

abbruch=0;
for (i=2;i<=knoten;i++) colour[i]=0;
colour[1]=1;

faerbe(colour,graph,1,&abbruch);


return(!abbruch);
}



/***************************ENDMINI*****************************************/


int endmini ( graph, adj, ready )

GRAPH graph;
/* graph[0][0] ist die anzahl der knoten, graph[0][1] die anzahl der kanten
   ( also tiefe ) und graph[0][2] die taille, falls schon festgelegt */

ADJAZENZ adj; /* adj[i] ist der grad von knoten i */
int ready; /* is it a cubic graph (so no edges missing) */

{
int i,j,l,abbruch;
BILD image,preimage;

// is the vertex on a minimal cycle?
#define SET_MINCYC(i) _mincycs[i]=mincycvalue;
#define IS_MINCYC(i) (_mincycs[i]==mincycvalue)
#define NOT_MINCYC(i) (_mincycs[i]!=mincycvalue)
#define RESETMINCYCS { int _i; if (mincycvalue==INT_MAX) \
		    { for (_i=0;_i<=knoten;_i++) _mincycs[_i]=0; mincycvalue=1; } \
  else mincycvalue++; }
 static int _mincycs[knoten+1]={0}, mincycvalue=1;
 static unsigned long long int bitvek[knoten+1]={0ULL,28ULL,0ULL,0ULL};// 28=4+8+16

for (i=1; i<=graph[0][0]; i++) image[i]=unbelegt;

 if (!ready || (graph[0][2]>4)) 
   { for (i=1; i<=graph[0][0]; i++) SET_MINCYC(i); }
 else
   { RESETMINCYCS;
     //if (graph[0][2]<5)
       for (i=2; i<=graph[0][0];i++) // 1 always the same
	 { bitvek[i]= (1ULL<<graph[i][0]) | (1ULL<<graph[i][1]) | (1ULL<<graph[i][2]);
	 }
     if (graph[0][2]==3)
       {
	 for (i=1;i<=graph[0][0];i++) // sorted neighbours
	   //if NOT_MINCYC(i) marking from the smallest requires start from already marked one
	   { if (graph[i][1]>i)
	       { if (bitvek[graph[i][1]] & (1ULL<<graph[i][2]))
		   { SET_MINCYC(i); SET_MINCYC(graph[i][1]);
		     SET_MINCYC(graph[i][2]); 
		   }
	       }
	     if (graph[i][0]>i)
	       { if (bitvek[graph[i][0]] & (1ULL<<graph[i][1]))
		   { SET_MINCYC(i); SET_MINCYC(graph[i][0]);
		     SET_MINCYC(graph[i][1]); 
		   }
		 if (bitvek[graph[i][0]] & (1ULL<<graph[i][2]))
	       { SET_MINCYC(i); SET_MINCYC(graph[i][0]);
		 SET_MINCYC(graph[i][2]); 
	       }
	       }
	   }
       }// end girth 3
     else //if (graph[0][2]==4)
       {
	 for (i=1;i<=graph[0][0];i++) // sorted neighbours
	   //if NOT_MINCYC(i) marking from the smallest requires start from already marked one
	   { if (graph[i][1]>i)
	       { if (bitvek[graph[i][1]] & bitvek[graph[i][2]] & (~(1ULL<<i)))
		   { SET_MINCYC(i); SET_MINCYC(graph[i][1]);
		     SET_MINCYC(graph[i][2]);
		     l=graph[i][1];
		     // de eerste is zeker kleiner en is i of een andere kleinere die later
		     // de 4-hoek markeert
		     if (bitvek[graph[i][2]] & (1ULL<<graph[l][1])) SET_MINCYC(graph[l][1]);
		     if (bitvek[graph[i][2]] & (1ULL<<graph[l][2])) SET_MINCYC(graph[l][2]);
		   }
	       }
	     if (graph[i][0]>i)
	       { 
		 if (bitvek[graph[i][0]] & bitvek[graph[i][1]] & (~(1ULL<<i)))
		   { SET_MINCYC(i); SET_MINCYC(graph[i][0]);
		     SET_MINCYC(graph[i][1]);
 		     l=graph[i][1];
		     if (bitvek[graph[i][0]] & (1ULL<<graph[l][1])) SET_MINCYC(graph[l][1]);
		     if (bitvek[graph[i][0]] & (1ULL<<graph[l][2])) SET_MINCYC(graph[l][2]);

		   }
		 if (bitvek[graph[i][0]] & bitvek[graph[i][2]] & (~(1ULL<<i)))
		   { SET_MINCYC(i); SET_MINCYC(graph[i][0]);
		     SET_MINCYC(graph[i][2]);
  		     l=graph[i][2];
		     if (bitvek[graph[i][0]] & (1ULL<<graph[l][1])) SET_MINCYC(graph[l][1]);
		     if (bitvek[graph[i][0]] & (1ULL<<graph[l][2])) SET_MINCYC(graph[l][2]);

		   }
	       }
	   }
       }// end girth 4
     /*else // OK -- general BFS


doesn't help much -- not worth taking the risk of difficult testing

       {
#define MARK(i) _mark[i]=markvalue;
#define ISMARKED(i) (_mark[i]==markvalue)
#define UNMARKED(i) (_mark[i]!=markvalue)
#define RESETMARKS { int _i; if (markvalue==INT_MAX) \
		    { for (_i=0;_i<=knoten;_i++) _mark[_i]=0; markvalue=1; } \
  else markvalue++; }
 static int _mark[knoten+1]={0}, markvalue=1;
 int maxdistance, maxwritedistance, ancestor[knoten+1], dist[knoten+1], list[knoten+1];
 int run, length, d, top, m;


	 if (graph[0][2]%2)
	   {maxdistance=maxwritedistance=(graph[0][2]-1)/2;}
	 else
	   { maxdistance=(graph[0][2]-2)/2;
	     maxwritedistance=maxdistance+1;
	   }
	 for (i=1;i<=graph[0][0];i++)
	   if NOT_MINCYC(i)
	     { RESETMARKS;
	       MARK(i); ancestor[i]=0;
	       list[0]=graph[i][0];list[1]=graph[i][1];list[2]=graph[i][2];
	       MARK(graph[i][0]); MARK(graph[i][1]); MARK(graph[i][2]); 
	       dist[0]=dist[1]=dist[2]=1;
	       ancestor[graph[i][0]]=ancestor[graph[i][1]]=ancestor[graph[i][2]]=i;
	       length=3;
	       run=0;
	       while((run<length) && ((d=dist[run])<=maxdistance))
		 { 
		   top=list[run];
		   for (l=0;l<3;l++)
		     if (ISMARKED(graph[top][l]))
		       {
			 if (graph[top][l]!=ancestor[top])
			   { SET_MINCYC(i);
			     for (m=graph[top][l];m!=i;m=ancestor[m]) SET_MINCYC(i);
			     for (m=top;m!=i;m=ancestor[m]) SET_MINCYC(i);
			   }
		       }
		     else
		       if (d<maxwritedistance)
		       {
			 dist[length]=d+1;
			 list[length]=graph[top][l];
			 MARK(graph[top][l]);
			 ancestor[graph[top][l]]=top;
			 length++;
		       }
		   run++;
		 }
	     }
	 //for (i=1;i<=graph[0][0];i++) fprintf(stderr,"%d",IS_MINCYC(i));
	 //fprintf(stderr,"\n");
       }*/
 
   }
 
/*
 if (graph[0][0]==12 && graph[0][1]==18 && graph[2][2]==4 &&
     graph[9][1]==11 && graph[9][2]==12 && graph[7][2]==8 && graph[8][0]==7)
   { fprintf(stderr,"---------------------------------------\n");
     schreibegraph(graph);
     for (i=1; (i<=graph[0][0]); i++) fprintf(stderr,"%d: %d\n",i,IS_MINCYC(i));
     for (i=1; (i<=graph[0][0]); i++) writeset(bitvek[i]); exit(0);
   }
*/



abbruch = 0;

for (i=1; (i<=graph[0][0]) && (!abbruch); i++)
          if (adj[i]==reg)
	    {  
	       transtest=0;
                if (roots[i].tiefe) /*d.h.: schon initialisiert*/
                  { for (j=0; (j<blaetterzahl[i]) && (!abbruch); j++)
                     normalweiter(graph,adj,blaetterkette[i][j]->bild,
                           blaetterkette[i][j]->urbild,blaetterkette[i][j]->wo,
                           blaetterkette[i][j]->nextbild,&abbruch);
                     if (vertex_trans && (!transtest)) abbruch=1; }
                        else 
			  { 
 if (IS_MINCYC(i))
			     {

image[i]=1;
preimage[1]=i;

/* erster Fall */
if (IS_MINCYC(graph[i][0]) && IS_MINCYC(graph[i][1]))
  {
    image[graph[i][0]]=2;
    image[graph[i][1]]=3;
    image[graph[i][2]]=4;
    preimage[2]=graph[i][0];
    preimage[3]=graph[i][1];
    preimage[4]=graph[i][2];

    normalweiter(graph,adj,image,preimage,graph[i][0],5,&abbruch);
  }

/* zweiter Fall */
if (!abbruch && IS_MINCYC(graph[i][0]) && IS_MINCYC(graph[i][2]))
{
image[graph[i][0]]=2;
image[graph[i][1]]=4;
image[graph[i][2]]=3;
preimage[2]=graph[i][0];
preimage[3]=graph[i][2];
preimage[4]=graph[i][1];

normalweiter(graph,adj,image,preimage,graph[i][0],5,&abbruch);

}


/* dritter Fall */
if (!abbruch && IS_MINCYC(graph[i][0]) && IS_MINCYC(graph[i][2]))
{
image[graph[i][0]]=3;
image[graph[i][1]]=4;
image[graph[i][2]]=2;
preimage[2]=graph[i][2];
preimage[3]=graph[i][0];
preimage[4]=graph[i][1];
normalweiter(graph,adj,image,preimage,graph[i][2],5,&abbruch);

}


/* vierter Fall */
if (!abbruch && IS_MINCYC(graph[i][1]) && IS_MINCYC(graph[i][0]))
{
image[graph[i][0]]=3;
image[graph[i][1]]=2;
image[graph[i][2]]=4;
preimage[2]=graph[i][1];
preimage[3]=graph[i][0];
preimage[4]=graph[i][2];

normalweiter(graph,adj,image,preimage,graph[i][1],5,&abbruch);

}


/* funfter Fall */ 
if (!abbruch && IS_MINCYC(graph[i][1]) && IS_MINCYC(graph[i][2]))
{
image[graph[i][0]]=4;
image[graph[i][1]]=2;
image[graph[i][2]]=3;
preimage[2]=graph[i][1];
preimage[3]=graph[i][2];
preimage[4]=graph[i][0];

normalweiter(graph,adj,image,preimage,graph[i][1],5,&abbruch);

}



/* sechster Fall */
if (!abbruch && IS_MINCYC(graph[i][1]) && IS_MINCYC(graph[i][2]))
{
image[graph[i][0]]=4;
image[graph[i][1]]=3;
image[graph[i][2]]=2;
preimage[2]=graph[i][2];
preimage[3]=graph[i][1];
preimage[4]=graph[i][0];

normalweiter(graph,adj,image,preimage,graph[i][2],5,&abbruch);

} 

if (vertex_trans && (!transtest)) abbruch=1;
image[i]=unbelegt;
image[graph[i][0]]=unbelegt;
image[graph[i][1]]=unbelegt;
image[graph[i][2]]=unbelegt;

                         }
}
}
/*if (abbruch) notaccept++; else accept++;*/

return (!abbruch);
}

/**************************EINFUGEN************************************/

void einfugen (graph,adj,v,w)
GRAPH graph;
ADJAZENZ adj;
unsigned char v, w;
/* Fuegt die Kante (v,w) in den Graphen graph ein. Dabei wird aber davon */
/* ausgegangen, dass in adj die wirklich aktuellen werte fuer die */
/* Adjazenzen stehen. Die adjazenzen werden dann aktualisiert. */

{
graph[v][adj[v]]=w;
graph[w][adj[w]]=v;
adj[v]++;
adj[w]++;
graph[0][1]++;

}

/*************************ENTFERNEN*************************************/

void entfernen (graph,adj,v,w)
GRAPH graph;
ADJAZENZ adj;
unsigned char v, w;
/* Entfernt die Kante (v,w) aus dem Graphen graph. Dabei wird aber davon */
/* ausgegangen, dass in adj die wirklich aktuellen werte fuer die */
/* Adjazenzen stehen. Die adjazenzen werden dann aktualisiert. */

{
graph[v][adj[v]-1]=leer;
graph[w][adj[w]-1]=leer;
adj[v]--;
adj[w]--;
graph[0][1]--;
}




/******************************KETTENENTFERNEN***********************/


void kettenentfernen( bis_wo )

unsigned char bis_wo;

{
int i;
for (i=1; i<= bis_wo; i++)
      if (blaetterkette[i] != nil) { free(blaetterkette[i]);
                                     blaetterkette[i]=nil;
				   }
}

/**********************************SUCHE*****************************/

void suche (graph,adj,zweig,abbruch)

GRAPH graph;
ADJAZENZ adj;
BRANCH *zweig;
int *abbruch;

/* suche durchlaeuft einfach den Baum auf der suche nach blaettern -- bei einem
   blatt schaltet es dann in den konstruktionsmodus um */

{

if ( zweig->blatt ) baumbau(graph,adj,zweig,abbruch);

else {
     if ( (zweig->links != nil) && (!zweig->totlinks) )
{                   suche(graph,adj,zweig->links,abbruch);}
     if ( (zweig->rechts != nil) && (!zweig->totrechts) && (!(*abbruch)) )
{                   suche(graph,adj,zweig->rechts,abbruch);}
     }
}

/*****************************ROOTINITIALISIEREN************************/


void rootinitialisieren(graph,adj,which,number)

/* initialisiert eine wurzel und die nachfolgenden blaetter. Es wird
   vorausgesetzt, dass knoten number grad reg hat */

GRAPH graph;
ADJAZENZ adj;
ROOT *which;  /* zeiger auf die zu bearbeitende wurzel */
int number; /* welcher knoten soll die 1 als bild bekommen */

{
BILD image,preimage;
int i;


which->wurzeltot=0;
which->tiefe=graph[0][1];
for (i=0; i<6; i++) which->tot[i]=0;
for (i=0; i<=knoten; i++) image[i]=unbelegt;
image[number]=1;
preimage[1]=number;
blaetterzahl[number]=6;

/* erster Fall */
image[graph[number][0]]=2;
image[graph[number][1]]=3;
image[graph[number][2]]=4;
preimage[2]=graph[number][0];
preimage[3]=graph[number][1];
preimage[4]=graph[number][2];
which->naechster[0]=(BRANCH *)malloc(sizeof(BRANCH));
 if (which->naechster[0] ==NULL) { fprintf(stderr,"Can't get more space (8)-- exit\n"); exit(1); }
zweiginitialisieren(graph,adj,which->naechster[0],nil,image,preimage,5,
                                            graph[0][1],graph[number][0]);
which->naechster[0]->zurwurzel=which;


/* zweiter Fall */

image[graph[number][1]]=4;
image[graph[number][2]]=3;
preimage[3]=graph[number][2];
preimage[4]=graph[number][1];

which->naechster[1]=(BRANCH *)malloc(sizeof(BRANCH));
 if (which->naechster[1] ==NULL) { fprintf(stderr,"Can't get more space (9)-- exit\n"); exit(1); }
zweiginitialisieren(graph,adj,which->naechster[1],nil,image,preimage,5,
                                           graph[0][1],graph[number][0]);
which->naechster[1]->zurwurzel=which;



/* dritter Fall */

image[graph[number][0]]=3;
image[graph[number][2]]=2;
preimage[2]=graph[number][2];
preimage[3]=graph[number][0];

which->naechster[2]=(BRANCH *)malloc(sizeof(BRANCH));
 if (which->naechster[2] ==NULL) { fprintf(stderr,"Can't get more space (10)-- exit\n"); exit(1); }
zweiginitialisieren(graph,adj,which->naechster[2],nil,image,preimage,
                                            5,graph[0][1],graph[number][2]);
which->naechster[2]->zurwurzel=which;



/* vierter Fall */

image[graph[number][1]]=2;
image[graph[number][2]]=4;
preimage[2]=graph[number][1];
preimage[4]=graph[number][2];

which->naechster[3]=(BRANCH *)malloc(sizeof(BRANCH));
 if (which->naechster[3] ==NULL) { fprintf(stderr,"Can't get more space (11) -- exit\n"); exit(1); }
zweiginitialisieren(graph,adj,which->naechster[3],nil,image,preimage,
                                            5,graph[0][1],graph[number][1]);
which->naechster[3]->zurwurzel=which;


/* f"unfter Fall */

image[graph[number][0]]=4;
image[graph[number][2]]=3;
preimage[3]=graph[number][2];
preimage[4]=graph[number][0];

which->naechster[4]=(BRANCH *)malloc(sizeof(BRANCH));
 if (which->naechster[4] ==NULL) { fprintf(stderr,"Can't get more space (12)-- exit\n"); exit(1); }
zweiginitialisieren(graph,adj,which->naechster[4],nil,image,preimage,
                                           5,graph[0][1],graph[number][1]);
which->naechster[4]->zurwurzel=which;





/* sechster Fall */

image[graph[number][1]]=3;
image[graph[number][2]]=2;
preimage[2]=graph[number][2];
preimage[3]=graph[number][1];

which->naechster[5]=(BRANCH *)malloc(sizeof(BRANCH));
 if (which->naechster[5] ==NULL) { fprintf(stderr,"Can't get more space (13) -- exit\n"); exit(1); }
zweiginitialisieren(graph,adj,which->naechster[5],nil,image,preimage,
                                            5,graph[0][1],graph[number][2]);
which->naechster[5]->zurwurzel=which;


image[graph[number][0]]=unbelegt;
image[graph[number][1]]=unbelegt;
image[graph[number][2]]=unbelegt;


}





/*************************MINIDAR************************************/

int minidar ( graph, adj )

GRAPH graph;
/* graph[0][0] ist die anzahl der knoten, graph[0][1] die anzahl der kanten
   ( also tiefe ) und graph[0][2] die taille, falls schon festgelegt */

ADJAZENZ adj; /* adj[i] ist der grad von knoten i */


{
int i,j,abbruch;


abbruch = 0;


for (i=1; (i<=graph[0][0]) && (!abbruch); i++)
          if (adj[i]==reg)
             { if (roots[i].tiefe) /*d.h.: schon initialisiert*/
	       { if (!(roots[i].wurzeltot))
               for (j=0; (j<6) && (!abbruch); j++) 
                    if (!(roots[i].tot[j]))
                    suche(graph,adj,roots[i].naechster[j],&abbruch);
                 if (vertex_trans && roots[i].wurzeltot) abbruch=1; }

                        else 
                         { rootinitialisieren(graph,adj,roots+i,i);
                           for (j=0; j<6; j++) 
                           baumbau(graph,adj,roots[i].naechster[j],&abbruch);
                           if (vertex_trans && roots[i].wurzeltot) abbruch=1;
                         }
              }

/*if (abbruch) notaccept++; else accept++;*/

return (!abbruch);
}
    


/*************************TAILLENWEITE***********************************/


int taillenweite(graph)

/* berechnet die taillenweite eines Graphen, der in Minimaldarstellung ist */
/* und den im ganzen Programm benutzten Anforderungen an einen Teilgraphen */
/* eines regulaeren Graphen genuegt. Die taille wird in g[0][2] geschrieben */
/* die Funktion gibt 0 zurueck, wenn noch kein Zykel existiert, 1 sonst */

GRAPH graph;

{ unsigned char herkunft, tiefe, vertex;

tiefe=1; herkunft=1; graph[0][2]=0;
for (vertex=2; vertex!=1;) { if (graph[vertex][0] != herkunft)
			       { herkunft=vertex; vertex=graph[vertex][0]; }
                             else { if (graph[vertex][1]!=leer)
			       { herkunft=vertex; vertex=graph[vertex][1]; }
                                 else break; }
                              tiefe++;
                           }
if ((vertex==1) && (herkunft==3)) { graph[0][2]=tiefe; return(1); }

                                  else return(0);

}




/***************************KONSTRUKT*********************************/


void konstrukt(graph,adj,vertex)
/* konstruiert weiter am graphen und zwar am Knoten vertex. Es wird */
/* versucht, solange Knoten einzufuegen, bis der Grad von vertex reg ist */
GRAPH graph;
ADJAZENZ adj;
unsigned char vertex;

{
long int anfang, i, test,baumgebaut;

if (modulo_minibaum && adj[vertex]<3 && (graph[0][1]==splitlevel_minibaum) && !recover)
    { splitzaehler++; if (splitzaehler== mod_minibaum) splitzaehler=0;
      if (splitzaehler != rest_minibaum) return; }


if (graph[0][1] > maximalkanten) maximalkanten=graph[0][1];

if ( (!singleout_minibaum)  || recover
           || (graph[0][1]>=mindestkantenzahl[vertex][graph[0][2]]) )
if ( (!partial) || 
(memocomp(graph[2],endgraph[2],reg*(vertex-2)+adj[vertex])<=0) )

{
if (flag && (!recover))
  {  
      if ((vertex<=knotenzahl) && (adj[vertex]>0) && (adj[vertex]<reg))
        { flag=0; wegspeichern_minibaum(globalliste,graphzahlen,graph); }

/* Wenn der jetzige Graph ein fertiger graph waere und als lastgraph 
   weggeschrieben wuerde und dann sofort das Programm beendet wuerde,
   wuerde er beim restart nicht in die liste aufgenommen */
   
  }


if ((vertex==treeborder+1) && (blaetterzahl[0]==0))
    { blaetterzahl[0]=graph[0][1]; baumanalyse(graph[0][0]); }


if (vertex>knotenzahl)
                  { if (recover) { recover=0; }
                    else { if (endmini(graph,adj,1))
                             { if (independent)
                                { test=1;
                                for(i=indepborder+1; 
                                   (i<=indepstoreborder) && (i<vertex) &&
                                    test; i++)
                                        test= (independence(graph,i));
                                if (test && (indepstoreborder < knotenzahl))
                                        test= (independence(graph,knotenzahl));
                                if (test) aufschreiben_minibaum(graph);
                                for(i-- ; i>indepborder; i--) 
                                                  free(indeps[i]); }
                               else aufschreiben_minibaum(graph); }
			 }
		  }
else
   {
   if (adj[vertex]==0)
        { if (recover) { recover=0; }
             else   if (!singleout_minibaum)
                    { if (vertex<=border+2) aufschreiben_minibaum(graph);
		      else if (endmini(graph,adj,1)) aufschreiben_minibaum(graph); }
        }
   else
      if (adj[vertex]==reg) {
                             if (independent && (vertex<=indepborder))
                             { 
                                if (independence(graph,vertex))
                                konstrukt(graph,adj,vertex+1);
                                free(indeps[vertex]); 
                             }
                             else konstrukt(graph,adj,vertex+1);
                            }


         else                       /* d.h. 0 < adj[vertex] < reg */
            {
            if ((recover) && (lastgraph[vertex][adj[vertex]]==leer))
                     { recover=0; }
            anfang=graph[vertex][adj[vertex]-1]+1; 
            if (anfang<vertex+1) anfang=vertex+1;
                 if (graph[0][2]) /* d.h. taille ist schon festgelegt */
                     { for (i=anfang; i<=graph[0][0]; i++)
                        if ((adj[i]<reg) &&  (!forbidden[vertex][i]) && 
                                (abstandok(graph,vertex,i,graph[0][2])) &&
                            ( (!bipartite) || (colour[i]!=colour[vertex]) ) &&
			    (!forbiddencycles || nocyc(graph,adj,vertex,i)) )
                               { einfugen(graph,adj,vertex,i);
                                 baumgebaut=0;
                                 if (recover) test=
                                         (memocomp(graph[1],lastgraph[1],
                                           reg*(vertex-1)+adj[vertex])==0);
                                                 else test=1;
/* Achtung ! Die vorige Zeile benutzt explizit eine vorgegebene Form von
   GRAPH -- bei einer anderen Definition kommt es zu Fehlern ! */

                                   if (test)
				       { if (vertex<=treeborder)
					     { baumgebaut = 1;
                                               test = minidar(graph,adj);
					     }
                                         else 
                                            if ((vertex<=border) && (!recover))
					      test=endmini(graph,adj,0);
				        }

                                 if (test) konstrukt(graph,adj,vertex);
                                 if (baumgebaut)
                                 baumaufraeumen(graph[0][1],graph[0][0]);
                                 entfernen(graph,adj,vertex,i);
                               }
                       if (graph[0][0]+reg-adj[vertex]<=knotenzahl)
                       {  graph[0][0]++;
                          einfugen(graph,adj,vertex,graph[0][0]);
                          if (bipartite)
			    { if (colour[vertex]==1) 
				{ colour[graph[0][0]]=2; (anzahl_colours[2])++; }
			    else 
			      { colour[graph[0][0]]=1; (anzahl_colours[1])++; }
			    }
                          baumgebaut=0;
                          if (recover) test=(memocomp(graph[1],lastgraph[1],
                                             reg*(vertex-1)+adj[vertex])==0);
                                                 else test=1;
                                             /* siehe oben */
			  if (bipartite) 
			    { if ((anzahl_colours[1]>knotenzahl/2) || (anzahl_colours[2]>knotenzahl/2))
				test=0; }
                                   if (test)
				       { if (vertex<=treeborder)
					     { baumgebaut = 1;
                                               test = minidar(graph,adj);
					     }
                                         else 
                                            if ((vertex<=border) && (!recover))
					      test=endmini(graph,adj,0);
				        }

                          if (test) konstrukt(graph,adj,vertex);
                          if (baumgebaut)
                          baumaufraeumen(graph[0][1],graph[0][0]);
                          entfernen(graph,adj,vertex,graph[0][0]);
			  if (bipartite) (anzahl_colours[colour[graph[0][0]]])--;
                          graph[0][0]--;
                        }
                    }
                 else /* d.h. taille noch nicht belegt */
/* Der erste Zykel: */
/* hier kann besser ueberprueft werden,vertex->i kann naemlich nur die Werte
   2->3,3->5,5->7,7->11,11->15,15->23,23->31 und 31->47 annehmen */
                     {
                      switch(vertex)
			  {
			case 2: { i=3; break; }
		        case 3: { i=5; break; }
         	        case 5: { i=7; break; }
		        case 7: { i=11; break; }
		       case 11: { i=15; break; }
		       case 15: { i=23; break; }
		       case 23: { i=31; break; }
		       case 31: { i=47; break; }
		       default: { i=vertex; break; }
		         } 
                        if (graph[vertex][adj[vertex]-1]>i) i=vertex;

                        if ( (vertex != i) &&
                                (abstandok(graph,vertex,i,taille))&&
                                       (i-1+reg-adj[vertex]<=knotenzahl) &&
                            ( (!bipartite) || (colour[i]!=colour[vertex]) ) &&
			     (!forbiddencycles || nocyc(graph,adj,vertex,i)) )
                       /* Hier wird jetzt garantiert der erste Zykel erzeugt */
                               { einfugen(graph,adj,vertex,i);
                                 if (taillenweite(graph))
                                 {
                                 baumgebaut=0;
                                 if (recover) test=
                                              (memocomp(graph[1],lastgraph[1],
                                               reg*(vertex-1)+adj[vertex])==0);
                                                   else test=1;
                                               /* siehe oben */

                                   if (test)
				       { if (vertex<=treeborder)
					     { baumgebaut = 1;
                                               test = minidar(graph,adj);
					     }
                                         else 
                                            if ((vertex<=border) && (!recover))
					      test=endmini(graph,adj,0);
				       }
                                 
                                 if (test) konstrukt(graph,adj,vertex);
                                 graph[0][2]=0;
                                 if (baumgebaut)
                                 baumaufraeumen(graph[0][1],graph[0][0]);
			         }
                                 entfernen(graph,adj,vertex,i);
                               }
                       if (graph[0][0]+reg-adj[vertex]<=knotenzahl)
                       {  einfugen(graph,adj,vertex,graph[0][0]+1);
                          graph[0][0]++;
			  if (bipartite)
			    { if (colour[vertex]==1) 
				{ colour[graph[0][0]]=2; (anzahl_colours[2])++; }
			    else 
			      { colour[graph[0][0]]=1; (anzahl_colours[1])++; }
			    }
                          baumgebaut=0;
                          if (recover) test=(memocomp(graph[1],lastgraph[1],
                                             reg*(vertex-1)+adj[vertex])==0);
                                  else test=1;
                                               /* siehe oben */
			  if (bipartite) 
			    { if ((anzahl_colours[1]>knotenzahl/2) || (anzahl_colours[2]>knotenzahl/2))
				test=0; }
                                   if (test)
				       { if (vertex<=treeborder)
					     { baumgebaut = 1;
                                               test = minidar(graph,adj);
					     }
                                         else 
                                            if ((vertex<=border) && (!recover))
					      test=endmini(graph,adj,0);
				        }

                                 
                          if (test) konstrukt(graph,adj,vertex);
			  if (bipartite) (anzahl_colours[colour[graph[0][0]]])--;
                          if (baumgebaut)
                          baumaufraeumen(graph[0][1],graph[0][0]);
                          entfernen(graph,adj,vertex,graph[0][0]);
                          graph[0][0]--;
                        }
                    }

               } /* ende 0 < adj < reg */
         } /* ende vertex <= knotenzahl */

if ((vertex==treeborder+1) && (blaetterzahl[0]==graph[0][1])) 
                    {  kettenentfernen(graph[0][0]); blaetterzahl[0]=0; }

} /* ende if "not partial" oder "<=endgraph" */
} /* ende konstruct */












/**************************LESESTARTGRAPH********************************/

int lesestartgraph(graph,adj,startknoten)
GRAPH graph;
ADJAZENZ adj;
unsigned char *startknoten;

/* Belegt graph und adj der startvorgabe entsprechend */
/* fuer den Beginn der Konstruktion -- initialisierung von graph und adj*/
/* wird vorausgesetzt */
/* Achtung ! hier wird die genaue definition von GRAPH benutzt !!!!!!!*/

{
int lesetest, indtest;
FILE* fil;
int vertex1, vertex2, altvertex1, altvertex2, puf;

indtest=1;

fil = fopen(startgraphname,"r");
if (fil==0) { fprintf(stderr,"No startgraphfile there\n"); exit(0); }

altvertex1=altvertex2=graph[0][0]=1;

for (;(lesetest=fscanf(fil,"%d %d",&vertex1, &vertex2)) == 2; )
    { if ( (altvertex1>vertex1) || 
               ( (altvertex1==vertex1) && (altvertex2>=vertex2) ) )
       { fprintf(stderr,"Startgraph not given correctly (0)\n"); return(0); }
      if (vertex1>=vertex2) 
       { fprintf(stderr,"Startgraph not given correctly (1)\n"); return(0); }
      if (vertex1>1) if (adj[vertex1-1]!=3)
       { fprintf(stderr,"Startgraph not given correctly (2)\n"); return(0); }
      if ((adj[vertex1]>2) || (adj[vertex2]>2))
       { fprintf(stderr,"Startgraph not given correctly (3)\n"); return(0); }
      if ((vertex2>graph[0][0]+1) || (vertex2>knoten))
       { fprintf(stderr,"Startgraph not given correctly (4) (not maximal)\n");
             return(0); }

   if ((vertex1==treeborder+1) && (blaetterzahl[0]==0))
    {   if (!(minidar(graph,adj)))
        { fprintf(stderr,"Startgraph not given correctly (5) (not maximal)\n");
           return(0); }
        blaetterzahl[0]=graph[0][1]; baumanalyse(graph[0][0]);  }

      if (vertex2>graph[0][0]) graph[0][0]=vertex2;
      einfugen(graph,adj,(unsigned char)vertex1,(unsigned char)vertex2);
      altvertex1=vertex1; altvertex2=vertex2;
      if (graph[0][1] > maximalkanten) maximalkanten=graph[0][1];
      if ( independent && indtest && (vertex1 >=2) && (adj[vertex1]==3) )
                  for (puf=vertex1; (adj[puf]==3) && indtest; puf++)
                               indtest =  independence(graph,puf);

      }

if (lesetest==1) 
     { fprintf(stderr,"Startgraph not given correctly (6) -- odd number of entries\n");
           return(0); }

for (*startknoten=1; adj[*startknoten]==reg; (*startknoten)++);

if ((*startknoten)<2)
     { fprintf(stderr,"Startgraph not given correctly (7) -- too few edges\n");
           return(0); }


if ( (*startknoten)<=treeborder )
   {  if (!(minidar(graph,adj)))
        { fprintf(stderr,"Startgraph not given correctly (8) (not maximal)\n");
           return(0); }
   } else
if ((*startknoten==treeborder+1) && (blaetterzahl[0]==0))
    { if (!(minidar(graph,adj)))
        { fprintf(stderr,"Startgraph not given correctly (9) (not maximal)\n");
           return(0); }
      blaetterzahl[0]=graph[0][1]; baumanalyse(graph[0][0]); 
    } else
{
if (!recover)
  if (!(endmini(graph,adj,0)))
       { fprintf(stderr,"Startgraph not given correctly (10) (not maximal)\n");
           return(0); }
}

if (bipartite && !(bipartit(graph)))
   { fprintf(stderr,"Startgraph not given correctly: \n Option bipartite used, but startgraph not bipartite.\n"); return(0); }

taillenweite(graph);
if ((graph[0][2] < taille) && (graph[0][2]>0))
   { fprintf(stderr,"Startgraph not given correctly: Girth too small.\n"); 
     return(0); }

    if (!indtest)
   { fprintf(stderr,"Startgraph will always produce graphs with larger independent sets.\n"); 
     return(0); }
    

fclose(fil);

return(1);

}



/**************************LESEGRAPH********************************/

void lesegraph(graph,name)
GRAPH graph;
char *name;

/* Liest einen graphen aus der Datei "name" in graph                  */
/* Muss auf jeden Fall aufgerufen werden, bevor (!!!) der minibaum    */
/* aufgebaut wird */
/* Achtung ! hier wird die genaue definition von GRAPH benutzt !!!!!!!*/


{
int i,j,lesetest;
FILE* fil;
int vertex1, vertex2, altvertex1, altvertex2, x;
unsigned char adjaz[knoten+1];


for (i=0; i<=knoten; i++) { for (j=0; j<reg; j++) graph[i][j]=leer;
                            adjaz[i]=0; }
for (i=0; i<reg; i++) graph[0][i]=0;

fil = fopen(name,"r");
if (fil==0) { fprintf(stderr,"No startfile there\n"); exit(0); }

altvertex1=altvertex2=graph[0][0]=1;

for (;(lesetest=fscanf(fil,"%d %d",&vertex1, &vertex2)) == 2; )
    { if ( (altvertex1>vertex1) || 
               ( (altvertex1==vertex1) && (altvertex2>=vertex2) ) )
       { fprintf(stderr,"In %s : Graph not given correctly (0)\n", name); 
                                                                  exit(0); }
      if (vertex1>=vertex2) 
       { fprintf(stderr,"In %s :Graph not given correctly (1)\n",name); 
                                                                  exit(0); }
      if (vertex1>1) if (adjaz[vertex1-1]!=3)
       { fprintf(stderr,"In %s : Graph not given correctly (2)\n",name); 
                                                                  exit(0); }
      if ((adjaz[vertex1]>2) || (adjaz[vertex2]>2))
       { fprintf(stderr,"In %s : Graph not given correctly (3)\n",name); 
                                                                  exit(0); }
      if ((vertex2>graph[0][0]+1) || (vertex2>knoten))
{ fprintf(stderr,"In %s : Graph not given correctly (4) (not maximal)\n",name);
             exit(0); }

      if (vertex2>graph[0][0]) graph[0][0]=vertex2;
      einfugen(graph,adjaz,(unsigned char)vertex1,(unsigned char)vertex2);
      altvertex1=vertex1; altvertex2=vertex2;
      if (graph[0][1] > maximalkanten) maximalkanten=graph[0][1];
      }

if (lesetest==1) 
{ fprintf(stderr,"In %s : Graph not given correctly (5) -- odd number of entries\n",name);
           exit(0); }

if (strcmp(name,startgraphname)==0) 
{ for (x=1; adjaz[x]==reg; x++);
if (x<2)
     { fprintf(stderr,"Startgraph not given correctly (6) -- too few edges\n");
           exit(0); }
}

 if (!(endmini(graph,adjaz,0)))
{fprintf(stderr,"In %s : Graph not given correctly (7) (not maximal)\n",name);
           exit(0); }

if (strcmp(name,startgraphname)==0) 
{ if (bipartite && !(bipartit(graph)))
   { fprintf(stderr,"Startgraph not given correctly: \n Option bipartite used, but startgraph not bipartite.\n"); exit(0); }

taillenweite(graph);
if ((graph[0][2] < taille) && (graph[0][2]>0))
   { fprintf(stderr,"Startgraph not given correctly: Girth too small.\n"); 
     exit(0); }
}


fclose(fil);

}

/**************************INITIALIZE*********************************/

void initialize(graph,adj)
GRAPH graph;
ADJAZENZ adj;

/* initialisiert graph, adj und roots*/
/* Achtung ! hier wird die genaue definition von GRAPH benutzt !!!!!!!*/

{
int i,j;

blaetterzahl[0]=adj[0]=0;

for (i=1; i<=knoten; i++) { for (j=0; j<reg; j++) graph[i][j]=leer;
                            adj[i]=0;
                            roots[i].tiefe=0;
                            blaetterkette[i]=nil;
                            blaetterzahl[i]=0; }
for (i=0; i<reg; i++) graph[0][i]=0;

if (independent)
   { numberofinds[2]=1; numberofinds[1]=2;
     indeps[1] = (IND_SET *)malloc(numberofinds[1]*sizeof(IND_SET));
 if (indeps[1] ==NULL) { fprintf(stderr,"Can't get more space (14) -- exit\n"); exit(1); }
     (indeps[1][0]).drin = ~0ULL;
     (indeps[1][0]).groesse=0;
     (indeps[1][1]).drin = ~0ULL;
     (indeps[1][1]).drin ^= (MASK[2] | MASK[3] | MASK[4]);
     (indeps[1][1]).groesse=1; }


}

/*****************************DEFAULTWERTE*****************************/

void defaultwerte(graph,adj,startknoten)

GRAPH graph;
unsigned char *adj;
unsigned char *startknoten;

/* schreibt die defaultwerte in einen zuvor initialisierten graphen und
   aktualisiert auch adj */

{
int i;

*startknoten=2;

for (i=2; i<= reg+1; i++) { graph[1][i-2]=i;
                            graph[i][0]=1;
                            adj[i]=1;
                           }
adj[1]=reg;
graph[0][0]=reg+1;
graph[0][1]=reg;


i=bipartit(graph);

/* damit auch der minidar-baum vernuenftig initialisiert wird: */
if (treeborder>=1) minidar(graph,adj);
}

/******************************GRAPHERZEUGUNG**************************/

void grapherzeugung()
{
GRAPH graph;
unsigned char startknoten;
int test,i;

test=1;

initialize(graph,adj);

if (partial) { if (!recover) lesegraph(lastgraph,startgraphname);
               lesegraph(endgraph,endgraphname); 
               recover=1; }


if (startvorgabe) { test=lesestartgraph(graph,adj,&startknoten); 
		    if (modulo_minibaum && (splitlevel_minibaum < graph[0][1]))
			{ fprintf(stderr,"Error: Splitlevel too small !\n");
			  exit(0); } }
                  else { defaultwerte(graph,adj,&startknoten); 
                         bipartit(graph); }
/* Auch bei partial muessen die Defaultwerte gelesen werden */

if (bipartite) for (i=1; i<=graph[0][0]; i++) (anzahl_colours[colour[i]])++;

if (test) konstrukt(graph,adj,startknoten);

 wegspeichern_minibaum(globalliste,graphzahlen,graph);
}

/******************BELEGEMINIMALWERTE*********************************/

void belegeminimalwerte()
{
int i,j,k,gesamt,puffer, puffer2, kantenzahl2, dummy; 
                    /* in Kantenzahl2 ist immer 2*kantenzahl */
FILE *fil;
signed char x1,x2,x3,min;


 for (i=0; i<=knoten; i++) for (j=0; j<=knoten; j++) 
  for (k=0; k<=knoten; k++) 
     { gesamt = i+j+k;
       kantenzahl2= i + (2*j) + (3*k); /* das ist schon kantenzahl*2 */
/* Nun das Theorem von Kathy Jones: i>= (13v-2e)/28 */
puffer= (13*gesamt) - kantenzahl2;
if (!(puffer%28)) { puffer2= puffer/28; dummy=gesamt/8;
/* Nun das Resultat von B.Bajnok: */
		    if (gesamt && ((taille >= 5) || (gesamt % 8) || 
				      (kantenzahl2 != 20*dummy)))
			                       puffer2++;
/* Sicherstellen, dass es sich nicht um mehrere Kopien des H-Graphen handelt */
                    minimalwerte[i][j][k]=puffer2; }
else { puffer2=puffer/28;
       if (28*puffer2 > puffer) 
            { fprintf(stderr,"wrong truncation direction\n");
                    minimalwerte[i][j][k]=puffer2; }
       else { puffer2++;
                    minimalwerte[i][j][k]=puffer2; }
	  } 
     } /* ende 3-fache for-schleife */


/* nun die genau berechneten werte: */

fil= fopen("MINIMALWERTE","r");
if (fil==0) { fprintf(stderr,"Die Datei MINIMALWERTE fehlt !!\n");
              /*exit(0);*/ }
else
  { /* it is assumed that the file is ok and contains n complete sets of data 
       that is 4*n bytes for some n */
  while ((x1=getc(fil)) != EOF)
    { x2=getc(fil); x3=getc(fil); min=getc(fil);
      puffer=minimalwerte[x1][x2][x3];
      if ((puffer>min) && (taille==4)) { fprintf(stderr,"Dangerous error !\n");
                        exit(80); }
      if (puffer < min)
      minimalwerte[x1][x2][x3]=min;
/*fprintf(stderr,"%d %d %d : %d \n",x1,x2,x3,min);*/
    }

  }
fclose(fil);

/*fprintf(stderr,"Ende !!\n");*/

}

void call_minibaum(int order, int isBipartite, int allGraphs){
    int i, j;

    anzahl_colours[0] = anzahl_colours[1] = anzahl_colours[2] = 0;

    for (i = 0; i <= knoten; i++) {
        colour[i] = 0;
        for (j = 0; j <= knoten; j++)
            forbidden[i][j] = 0;
    }

    MASK[1] = 1ULL;
    for (i = 2; i <= 64; i++) MASK[i] = MASK[i - 1] << 1;
    for (i = 1; i <= 64; i++) NMASK[i] = ~MASK[i];

    timerestricted = bipartite = recover = startvorgabe = singleout_minibaum = lastout = partial = 0;
    observer = TIMERESTRICTED = vertex_trans = independent = maximalkanten = modulo_minibaum = 0;
    splitzaehler = oldtime = short_code = 0;
    snarks = zyconnecti = zyminconn = 0;
    number_of_graphs = 0ULL;

    knotenzahl = order;
    taille = (isBipartite ? 4 : 3);
    bipartite = (isBipartite ? 1 : 0);
    if(allGraphs){
       noout_minibaum = 1;
       lastout = singleout_minibaum = 0;
    } else {
        noout_minibaum = singleout_minibaum = 1;
        lastout = 0;
    }


    S_intervall = 0; //cen be removed since wegspeichern isn't called anyway

    for (i = 4; i <= knoten; i += 2) {
        graphzahlen[i] = 0ULL;
    }

    /* das erste +1 fuer shortcodes worst case */

    for (i = 0; i <= knotenzahl + 1; i++) for (j = 0; j <= 10; j++)
            mindestkantenzahl[i][j] = 0;

    /* Die folgenden Zahlen beruhen einfach auf den Ueberlegungen, wann ein Zykel
       zu geringer Laenge schon im Rest entstehen muesste -- unabhaengig von dem
       schon erzeugten Graphen zum Zeitpunkt der Abfrage */

    if (singleout_minibaum) {
        int puffer = knotenzahl / 2;
        puffer = puffer * 3;
        for (i = 3; i <= 10; i++) /* i ist die taille */ {
            for (j = 0; j < i - 1; j++)
                if ((knotenzahl - j) >= 0) mindestkantenzahl[knotenzahl - j][i] = puffer - j;
            if ((knotenzahl - i + 1) >= 0) mindestkantenzahl[knotenzahl - i + 1][i] = puffer - i;
        }
        /* Jetzt die Faelle einzeln: */
        /* "knotenzahl-k" steht immer bei "k+1" restknoten */
        /* Taillenweite 3:  Leider keine Ideen */
        if (knotenzahl >= 6) {
            /* Taillenweite 4: */ mindestkantenzahl[knotenzahl - 4][4] = puffer - 6;
        }
        if (knotenzahl >= 10) {
            /* Taillenweite 5: */ mindestkantenzahl[knotenzahl - 5][5] = puffer - 6;
            mindestkantenzahl[knotenzahl - 6][5] = puffer - 8;
            mindestkantenzahl[knotenzahl - 7][5] = puffer - 10;
        }
        if (knotenzahl >= 14) {
            /* Taillenweite 6: */ mindestkantenzahl[knotenzahl - 6][6] = puffer - 7;
            mindestkantenzahl[knotenzahl - 7][6] = puffer - 9;
            mindestkantenzahl[knotenzahl - 8][6] = puffer - 10;
            mindestkantenzahl[knotenzahl - 9][6] = puffer - 12;
        }
        if (knotenzahl >= 22) {
            /* Taillenweite 7: */ mindestkantenzahl[knotenzahl - 7][7] = puffer - 8;
            mindestkantenzahl[knotenzahl - 8][7] = puffer - 9;
            mindestkantenzahl[knotenzahl - 9][7] = puffer - 11;
            mindestkantenzahl[knotenzahl - 10][7] = puffer - 12;
            mindestkantenzahl[knotenzahl - 11][7] = puffer - 14;
        }
        if (knotenzahl >= 30) {
            /* Taillenweite 8: */ mindestkantenzahl[knotenzahl - 8][8] = puffer - 9;
            mindestkantenzahl[knotenzahl - 9][8] = puffer - 10;
            mindestkantenzahl[knotenzahl - 10][8] = puffer - 12;
            mindestkantenzahl[knotenzahl - 11][8] = puffer - 13;
            mindestkantenzahl[knotenzahl - 12][8] = puffer - 14;
            mindestkantenzahl[knotenzahl - 13][8] = puffer - 16;
        }
        if (knotenzahl >= 46) {
            /* Taillenweite 9: */ mindestkantenzahl[knotenzahl - 9][9] = puffer - 10;
            mindestkantenzahl[knotenzahl - 10][9] = puffer - 11;
            mindestkantenzahl[knotenzahl - 11][9] = puffer - 12;
            mindestkantenzahl[knotenzahl - 12][9] = puffer - 14;
            mindestkantenzahl[knotenzahl - 13][9] = puffer - 15;
            mindestkantenzahl[knotenzahl - 14][9] = puffer - 16;
            mindestkantenzahl[knotenzahl - 15][9] = puffer - 18;
        }
        if (knotenzahl >= 60) {
            /* Taillenweite 10: */ mindestkantenzahl[knotenzahl - 10][10] = puffer - 11;
            mindestkantenzahl[knotenzahl - 11][10] = puffer - 12;
            mindestkantenzahl[knotenzahl - 12][10] = puffer - 13;
            mindestkantenzahl[knotenzahl - 13][10] = puffer - 15;
            mindestkantenzahl[knotenzahl - 14][10] = puffer - 16;
            mindestkantenzahl[knotenzahl - 15][10] = puffer - 17;
            mindestkantenzahl[knotenzahl - 16][10] = puffer - 18;
            mindestkantenzahl[knotenzahl - 17][10] = puffer - 20;
        }

    }

    switch (taille) {
        case 3: startzahl = 4;
            break;
        case 4: startzahl = 6;
            break;
        case 5: startzahl = 10;
            break;
        case 6: startzahl = 14;
            break;
        case 7: startzahl = 24;
            break;
        case 8: startzahl = 30;
            break;
        case 9: startzahl = 58;
            break; /* information by Brendan McKay */
        case 10: startzahl = 70;
    }

    if (knotenzahl < startzahl) {
        //no graphs
        return;
    }

    switch (taille) {
        case 3: if (knotenzahl < 16) border = knotenzahl - 6;
            else border = knotenzahl - 7;
            break;
        case 4: border = knotenzahl - 7;
            break;
        case 5: border = knotenzahl - 8;
            break;
        case 6: border = knotenzahl - 9;
            break;
        case 7: if (knotenzahl < 28) border = knotenzahl - 12;
            else border = knotenzahl - 13;
            break;
        case 8: border = knotenzahl - 17;
            break;
        case 9: border = knotenzahl - 22;
            break; /* just a guess */
        case 10: border = knotenzahl - 29;
            break; /* just a guess */
    }

    if (bipartite) border--;


    if (vertex_trans) border--; /* wundert mich selbst -- ich haette eher erwartet,
                               dass man border erhoehen muss */

    treeborder = (9 * border) / 10;

    /* PLAY with treeborder and border (always border>=treeborder) between here */

    //border= border-2;
    //treeborder= treeborder-2;

    /* and here */

    if (treeborder > baumgrenze) treeborder = baumgrenze;

    //fprintf(stderr,"%d %d\n",treeborder,border);

    //if (border<=baumgrenze) treeborder=border;
    //                      else treeborder=baumgrenze;


    if (treeborder < 0) treeborder = 0;

    indepborder = border + 2;
    if (indepborder <= indepspeichergrenze) indepstoreborder = indepborder;
    else indepstoreborder = indepspeichergrenze;


    grapherzeugung();
}

int is_minibaum_available(int order){
    if (sizeof(unsigned long long int)!=8 || sizeof(unsigned int)!=4){
        return 0;
    }
    if(order>knoten){
        return 0;
    }
    return 1;
}

/*********LAST*BUT*NOT*LEAST:****************MAIN***********************/

#ifdef MINIBAUM_NO_MAIN
    #define MINIBAUM_MAIN_FUNCTION minibaumnomain
#else
    #define MINIBAUM_MAIN_FUNCTION main
#endif
int MINIBAUM_MAIN_FUNCTION(argc,argv)

int argc;
char *argv[];


{
        fprintf(stderr, "This version has been modified for pregraphs and can't be used as standalone generator -- sorry.\n");
        exit(0);

    int i, j, puffer, system(), getpid(), dummy;
    FILE *fil;
    int init = 1, remarks = 0;
    char strpuf1[60], startname[45], unterbrechenname[50];
    char kommandoname[60], observername[50];
    struct tm *zeitzeiger;
    struct tms TMS;
    time_t b;
    unsigned int savetime;
    signed char test;
    int numforbcycles = 0;

    if (sizeof (unsigned long long int) != 8) {
        fprintf(stderr, "This version relies on 64 bit unsigned long longs -- sorry.\n");
        exit(0);
    }
    if (sizeof (unsigned int) != 4) {
        fprintf(stderr, "This version relies on 32 bit unsigned ints -- sorry.\n");
        exit(0);
    }


    anzahl_colours[0] = anzahl_colours[1] = anzahl_colours[2] = 0;

    for (i = 0; i <= knoten; i++) {
        colour[i] = 0;
        for (j = 0; j <= knoten; j++)
            forbidden[i][j] = 0;
    }

    for (j = 4; j <= knotenzahl; j += 2)
        for (i = 0; i < codelaenge; i++) oldcode1[j][i] = oldcode2[j][i] = whichcode[j] = 0;

    MASK[1] = 1ULL;
    for (i = 2; i <= 64; i++) MASK[i] = MASK[i - 1] << 1;
    for (i = 1; i <= 64; i++) NMASK[i] = ~MASK[i];



    timerestricted = bipartite = recover = startvorgabe = singleout_minibaum = lastout = partial = 0;
    observer = TIMERESTRICTED = vertex_trans = independent = maximalkanten = modulo_minibaum = 0;
    splitzaehler = oldtime = short_code = 0;
    snarks = zyconnecti = zyminconn = 0;
    number_of_graphs = 0ULL;

    if (argc < 3) {
        fprintf(stderr, "Wrong input data. At least the number of vertices and the minimum\n");
        fprintf(stderr, "girth must be specified.\n");
        exit(0);
    }
    else {
        sscanf(argv[1], "%d", &knotenzahl);
        sscanf(argv[2], "%d", &taille);
        for (i = 3; i < argc; i++) {
            switch (argv[i][0]) {
                case 'b':
                {
                    bipartite = 1;
                    break;
                }
                case 'i':
                {
                    independent = i;
                    singleout_minibaum = 1; /* fuer die kleineren werte wuerden
                                         durch die Belegung von five14s
                                         falsche Zahlen erzeugt */
                    maxindep = atoi(argv[i] + 1);
                    break;
                }
                case 'g':
                {
                    graph6_output_minibaum = 1;
                    break;
                }
                case 'm':
                {
                    modulo_minibaum = 1;
                    i++;
                    rest_minibaum = atoi(argv[i]);
                    i++;
                    mod_minibaum = atoi(argv[i]);
                    if (rest_minibaum >= mod_minibaum) {
                        fprintf(stderr, "rest must be smaller than mod !\n");
                        exit(1);
                    }
                    init = 0;
                    break;
                }
                case 'R':
                {
                    remarks = 1;
                    break;
                }
                case 'r':
                {
                    recover = 1;
                    break;
                }
                case 'c':
                {
                    startvorgabe = i; /* gleich zum merken statt 1 */
                    break;
                }
                case 'o':
                {
                    lastout = 1;
                    break;
                }
                case 's':
                {
                    if (argv[i][1] == 'h') short_code = 1;
                    else if (argv[i][1] == 0) singleout_minibaum = 1;
                    else {
                        fprintf(stderr, "wrong option\n");
                        exit(0);
                    }
                    break;
                }
                case 'S':
                {
                    snarks = 1;
                    break;
                }
                case 'v':
                {
                    vertex_trans = 1;
                    break;
                }
                case 'p':
                {
                    partial = i;
                    break;
                }
                case 'P':
                {
                    sprintf(strpuf1, "/etc/renice %s %d > yqvj ; rm yqvj ", argv[i] + 1, getpid());
                    dummy = system(strpuf1);
                    break;
                }
                case 't':
                {
                    timerestricted = i;
                    /*so merkt man sich gleich den index*/
                    i += 2;
                    break;
                }
                case 'T':
                {
                    TIMERESTRICTED = i;
                    /*so merkt man sich gleich den index*/
                    i += 2;
                    break;
                }

                case 'O':
                {
                    observer = 1;
                    break;
                }

                case 'n':
                {
                    if (argv[i][1] == 'i') init = 0;
                    else
                        if (strncmp(argv[i], (char *) "nocyc", (size_t) 5) == 0) {
                        numforbcycles++;
                        puffer = atoi(argv[i] + 5);
                        if ((puffer <= knotenzahl) && (puffer >= 3)) {
                            if (puffer > forbiddencycles) forbiddencycles = puffer;
                            lastout = singleout_minibaum = 1;
                            forbcyc[puffer] = 1;
                        }
                    } else
                        if (strncmp(argv[i], (char *) "noout", (size_t) 5) == 0) {
                        //noout=lastout=singleout=1;
                        noout_minibaum = 1;
                        lastout = singleout_minibaum = 0;
                    } else {
                        fprintf(stderr, "%s nonidentified option\n", argv[i]);
                        exit(0);
                    }
                    break;
                }

                case 'Z':
                {
                    zyconnecti = i;
                    zyminconn = atoi(argv[i] + 1);
                    if (taille < zyminconn) taille = zyminconn;
                    if (knotenzahl > 8 * sizeof (unsigned long long int)) {
                        fprintf(stderr, "Sorry -- testing cyclic connectivity is only possible\n");
                        fprintf(stderr, "up to 8*sizeof(unsigned long long int)=%ld vertices.\n",
                                8 * sizeof (unsigned long long int));
                        exit(0);
                    }
                    break;
                }


                default:
                {
                    fprintf(stderr, "%s nonidentified option\n", argv[i]);
                    exit(0);
                }
            }
        }
    }

    if (independent && partial) {
        fprintf(stderr, "Independent and partial can not be chosen together in this version.");
        fprintf(stderr, "\nSorry\n");
        exit(0);
    }

    if (lastout || graph6_output_minibaum) S_intervall = 0;
    else S_intervall = Sicherungsintervall;
    if (S_intervall) setup();

    if (independent) {
        if (taille > 3) {
            belegeminimalwerte();
            for (i = 0; i <= knoten; i++) {
                if (i % 14) {
                    five14s[i] = 5 * i;
                    five14s[i] = (five14s[i] / 14) + 1;
                } else five14s[i] = (5 * i) / 14;
            }
        } else if (maxindep > 1) {
            for (i = 0; i <= knoten; i++) if (i % 3) five14s[i] = (i / 3) + 1;
                else five14s[i] = i / 3;
        } /* Der name passt dann
                                            natuerlich nicht gerade */

        if (maxindep == 1) {
            for (i = 0; i <= knoten; i++) five14s[i] = 0;
        }
    }


    if (modulo_minibaum) {
        switch (taille) {
            case 3:
            {
                splitlevel_minibaum = 20;
                break;
            }
            case 4:
            {
                splitlevel_minibaum = 22;
                break;
            }
            case 5:
            {
                splitlevel_minibaum = 24;
                break;
            }
            case 6:
            {
                splitlevel_minibaum = 29;
                break;
            }
            case 7:
            {
                splitlevel_minibaum = 34;
                break;
            }
            case 8:
            {
                splitlevel_minibaum = 40; /* just a guess */
                fprintf(stderr, "contact gunnar@mathematik.uni-bielefeld.de ");
                fprintf(stderr,"whether the faster program for this task is ready!\n");
                break;
            }
            case 9:
            {
                splitlevel_minibaum = 54; /* just a guess */
                fprintf(stderr, "contact gunnar@mathematik.uni-bielefeld.de ");
                fprintf(stderr,"whether the faster program for this task is ready or better forget it!\n");
                break;
            }
            case 10:
            {
                splitlevel_minibaum = 76; /* just a guess */
                fprintf(stderr, "contact gunnar@mathematik.uni-bielefeld.de ");
                fprintf(stderr,"whether the faster program for this task is ready or better forget it!\n");
                break;
            }
            default:
            {
                fprintf(stderr, "Girth too large. \n");
                exit(1);
            }
        }

        if (bipartite && (taille < 6)) splitlevel_minibaum += 5;
        if (independent) splitlevel_minibaum += 3;
        if (forbiddencycles) splitlevel_minibaum += 8;
        if (snarks) splitlevel_minibaum += 3;

        if ((3 * knotenzahl) / 2 <= splitlevel_minibaum) splitlevel_minibaum = (3 * knotenzahl) / 2 - 1;

        minwrite = splitlevel_minibaum * 2 / 3 + 1;

    } /* ende if modulo */





    if ((TIMERESTRICTED && timerestricted) || (TIMERESTRICTED && observer)) {
        fprintf(stderr, "The option T may not be used together with t or O.\n");
        exit(0);
    }



    if (partial && startvorgabe) {
        fprintf(stderr, "The options p and c can not be used simultaneously\n");
        exit(0);
    }

    if (knotenzahl > knoten) {
        fprintf(stderr, "Knotenzahl zu gross !\n\n");
        exit(0);
    }
    if ((taille < 3) || (taille > 10)) {
        fprintf(stderr, "The girth must be between 3 and 10 !\n\n");
        exit(0);
    }

    for (i = 4; i <= knoten; i += 2) {
        codelength[i] = (reg * i) / 2 - reg;
        graphzahlen[i] = 0ULL;
        nextposition[i] = globalliste[i] =
                (unsigned char*) malloc((codelength[i] + 1) * listenlaenge + 4);
        if (nextposition[i] == NULL) {
            fprintf(stderr, "Can't get more space (15)-- exit\n");
            exit(1);
        }
    }

    /* das erste +1 fuer shortcodes worst case */

    for (i = 0; i <= knotenzahl + 1; i++) for (j = 0; j <= 10; j++)
            mindestkantenzahl[i][j] = 0;

    /* Die folgenden Zahlen beruhen einfach auf den Ueberlegungen, wann ein Zykel
       zu geringer Laenge schon im Rest entstehen muesste -- unabhaengig von dem
       schon erzeugten Graphen zum Zeitpunkt der Abfrage */

    if (singleout_minibaum) {
        puffer = knotenzahl / 2;
        puffer = puffer * 3;
        for (i = 3; i <= 10; i++) /* i ist die taille */ {
            for (j = 0; j < i - 1; j++)
                if ((knotenzahl - j) >= 0) mindestkantenzahl[knotenzahl - j][i] = puffer - j;
            if ((knotenzahl - i + 1) >= 0) mindestkantenzahl[knotenzahl - i + 1][i] = puffer - i;
        }
        /* Jetzt die Faelle einzeln: */
        /* "knotenzahl-k" steht immer bei "k+1" restknoten */
        /* Taillenweite 3:  Leider keine Ideen */
        if (knotenzahl >= 6) {
            /* Taillenweite 4: */ mindestkantenzahl[knotenzahl - 4][4] = puffer - 6;
        }
        if (knotenzahl >= 10) {
            /* Taillenweite 5: */ mindestkantenzahl[knotenzahl - 5][5] = puffer - 6;
            mindestkantenzahl[knotenzahl - 6][5] = puffer - 8;
            mindestkantenzahl[knotenzahl - 7][5] = puffer - 10;
        }
        if (knotenzahl >= 14) {
            /* Taillenweite 6: */ mindestkantenzahl[knotenzahl - 6][6] = puffer - 7;
            mindestkantenzahl[knotenzahl - 7][6] = puffer - 9;
            mindestkantenzahl[knotenzahl - 8][6] = puffer - 10;
            mindestkantenzahl[knotenzahl - 9][6] = puffer - 12;
        }
        if (knotenzahl >= 22) {
            /* Taillenweite 7: */ mindestkantenzahl[knotenzahl - 7][7] = puffer - 8;
            mindestkantenzahl[knotenzahl - 8][7] = puffer - 9;
            mindestkantenzahl[knotenzahl - 9][7] = puffer - 11;
            mindestkantenzahl[knotenzahl - 10][7] = puffer - 12;
            mindestkantenzahl[knotenzahl - 11][7] = puffer - 14;
        }
        if (knotenzahl >= 30) {
            /* Taillenweite 8: */ mindestkantenzahl[knotenzahl - 8][8] = puffer - 9;
            mindestkantenzahl[knotenzahl - 9][8] = puffer - 10;
            mindestkantenzahl[knotenzahl - 10][8] = puffer - 12;
            mindestkantenzahl[knotenzahl - 11][8] = puffer - 13;
            mindestkantenzahl[knotenzahl - 12][8] = puffer - 14;
            mindestkantenzahl[knotenzahl - 13][8] = puffer - 16;
        }
        if (knotenzahl >= 46) {
            /* Taillenweite 9: */ mindestkantenzahl[knotenzahl - 9][9] = puffer - 10;
            mindestkantenzahl[knotenzahl - 10][9] = puffer - 11;
            mindestkantenzahl[knotenzahl - 11][9] = puffer - 12;
            mindestkantenzahl[knotenzahl - 12][9] = puffer - 14;
            mindestkantenzahl[knotenzahl - 13][9] = puffer - 15;
            mindestkantenzahl[knotenzahl - 14][9] = puffer - 16;
            mindestkantenzahl[knotenzahl - 15][9] = puffer - 18;
        }
        if (knotenzahl >= 60) {
            /* Taillenweite 10: */ mindestkantenzahl[knotenzahl - 10][10] = puffer - 11;
            mindestkantenzahl[knotenzahl - 11][10] = puffer - 12;
            mindestkantenzahl[knotenzahl - 12][10] = puffer - 13;
            mindestkantenzahl[knotenzahl - 13][10] = puffer - 15;
            mindestkantenzahl[knotenzahl - 14][10] = puffer - 16;
            mindestkantenzahl[knotenzahl - 15][10] = puffer - 17;
            mindestkantenzahl[knotenzahl - 16][10] = puffer - 18;
            mindestkantenzahl[knotenzahl - 17][10] = puffer - 20;
        }




    }


    switch (taille) {
        case 3: startzahl = 4;
            break;
        case 4: startzahl = 6;
            break;
        case 5: startzahl = 10;
            break;
        case 6: startzahl = 14;
            break;
        case 7: startzahl = 24;
            break;
        case 8: startzahl = 30;
            break;
        case 9: startzahl = 58;
            break; /* information by Brendan McKay */
        case 10: startzahl = 70;
    }

    if (knotenzahl < startzahl) {
        fprintf(stderr, "Cubic graphs of the given size and girth do not exist.\n");
        exit(0);
    }

    if (short_code) sprintf(codename, "Sodes_00.03.00");
    else sprintf(codename, "Codes_00.03.00");
    codename[6] = (knotenzahl / 10) + 48;
    codename[7] = (knotenzahl % 10) + 48;
    codename[12] = (taille / 10) + 48;
    codename[13] = (taille % 10) + 48;
    if (bipartite) {
        i = strlen(codename);
        codename[i] = '.';
        codename[i + 1] = 'b';
        codename[i + 2] = 0;
    }
    if (vertex_trans) {
        i = strlen(codename);
        codename[i] = '.';
        codename[i + 1] = 'v';
        codename[i + 2] = 0;
    }


    if (independent) {
        sprintf(strpuf1, ".");
        strcat(codename, strpuf1);
        strcat(codename, argv[independent]);
    }

    if (snarks) {
        sprintf(strpuf1, ".sn");
        strcat(codename, strpuf1);
    }

    if (zyconnecti) {
        sprintf(strpuf1, ".cyc");
        strcat(codename, strpuf1);
        strcat(codename, argv[zyconnecti] + 1);
    }



    if (startvorgabe) {
        sprintf(strpuf1, ".");
        strcat(codename, strpuf1);
        strcat(codename, argv[startvorgabe]);
    }
    if (partial) {
        sprintf(strpuf1, ".");
        strcat(codename, strpuf1);
        strcat(codename, argv[partial]);
    }

    if (modulo_minibaum) { /* must be the last item added to the name */
        sprintf(strpuf1, ".m_%d_%d", rest_minibaum, mod_minibaum);
        strcat(codename, strpuf1);
    }



    sprintf(lastgraphname, "lastgraph");
    strcat(lastgraphname, codename + 5);

    if (S_intervall) {
        if ((fil = fopen(lastgraphname, "r")) != NULL)
            if (fread(lastgraph, 1, sizeof (GRAPH), fil) == sizeof (GRAPH)) {
                recover = 1;
                fprintf(stderr, "Using recover-file !\n");
                fclose(fil);
            }
    }

    if (startvorgabe || partial) {
        sprintf(startgraphname, "startgraph");
        strcat(startgraphname, codename + 5);
        if (modulo_minibaum) {
            for (i = 12; (startgraphname[i] != 'm') ||
                    (startgraphname[i - 1] != '.'); i++);
            startgraphname[i - 1] = 0;
            /*abschneiden*/
        }
        fil = fopen(startgraphname, "r");
        if (fil == 0) {
            fprintf(stderr, "No startfile %s there\n", startgraphname);
            exit(0);
        }
    }



    if (recover) {
        fil = fopen(lastgraphname, "r");
        if (fil == 0) {
            fprintf(stderr, "No lastgraph there ! \n");
            exit(0);
        } else fclose(fil);
    } else
        if (S_intervall) {
        fil = fopen(lastgraphname, "w");
        fclose(fil);
    }


    if (partial) {
        sprintf(endgraphname, "endgraph");
        strcat(endgraphname, codename + 5);
    }


    if (observer) {
        sprintf(strpuf1, " ");
        kommandoname[0] = 0;
        for (i = 0; i < argc; i++) {
            strcat(kommandoname, argv[i]);
            strcat(kommandoname, strpuf1);
        }
        if (!recover) {
            sprintf(strpuf1, "r ");
            strcat(kommandoname, strpuf1);
        }

        for (i = 1; i <= 6; i++) {
            sprintf(observername, "observer%1d", i);
            strcat(observername, codename + 5);

            sprintf(strpuf1, "observern%1d", i);
            strcat(strpuf1, codename + 5);
            fil = fopen(strpuf1, "w");
            fprintf(fil, "sh -c %s\n", observername);
            fclose(fil);

            fil = fopen(observername, "w");

            sprintf(strpuf1, "observern%1d", (i % 6) + 1);
            strcat(strpuf1, codename + 5);

            fprintf(fil, " \n exec > at.out 2>at.errout\n");
            fprintf(fil, " \n if test -f %s \n", lastgraphname);
            fprintf(fil, "    then \n   if ps -%d \n", getpid());
            fprintf(fil, "         then at %4d %s \n", ((i * 400) + 100) % 2400,
                    strpuf1);
            fprintf(fil, "         else %s \n         fi\n", kommandoname);
            fprintf(fil, "fi \n \n");
            fclose(fil);

            sprintf(strpuf1, "chmod 755 %s", observername);
            dummy = system(strpuf1);

        }

        b = time(0);
        zeitzeiger = localtime(&b);
        sprintf(observername, "sh -c observer%1d",
                ((zeitzeiger->tm_hour / 4) + 1));
        strcat(observername, codename + 5);

        dummy = system(observername);

    }



    if (timerestricted) {
        i = timerestricted;
        i++;
        starttime = (argv[i][0] - 48)*36000 +
                (argv[i][1] - 48)*3600 + (argv[i][2] - 48)*600
                + (argv[i][3] - 48)*60;
        i++;
        endtime = (argv[i][0] - 48)*36000 +
                (argv[i][1] - 48)*3600 + (argv[i][2] - 48)*600
                + (argv[i][3] - 48)*60;

        sprintf(startname, "start");
        strcat(startname, codename + 5);

        sprintf(unterbrechenname, "interrupt");
        strcat(unterbrechenname, codename + 5);

        sprintf(strpuf1, "startn");
        strcat(strpuf1, codename + 5);
        fil = fopen(strpuf1, "w");
        fprintf(fil, "sh -c %s\n", startname);
        fclose(fil);


        fil = fopen(unterbrechenname, "w");
        fprintf(fil, "exec > at.t.out 2>at.t.errout \n");
        fprintf(fil, "kill -STOP %d && at %s %s \n"
                , getpid(), argv[i - 1], strpuf1);
        fclose(fil);
        sprintf(strpuf1, "chmod 755 %s", unterbrechenname);
        dummy = system(strpuf1);

        sprintf(strpuf1, "interruptn");
        strcat(strpuf1, codename + 5);
        fil = fopen(strpuf1, "w");
        fprintf(fil, "sh -c %s\n", unterbrechenname);
        fclose(fil);


        fil = fopen(startname, "w");
        fprintf(fil, "exec > at.t.out 2>at.t.errout\n");
        fprintf(fil, "kill -CONT %d && at %s %s"
                , getpid(), argv[i], strpuf1);
        fclose(fil);
        sprintf(strpuf1, "chmod 755 %s", startname);
        dummy = system(strpuf1);

        b = time(0);
        zeitzeiger = localtime(&b);
        j = (3600 * zeitzeiger->tm_hour) + (60 * zeitzeiger->tm_min)
                + zeitzeiger->tm_sec;

        if (starttime < endtime) {
            if ((j > starttime) && (j < endtime))
                sprintf(strpuf1, "at %s interruptn%s", argv[i], codename + 5);
            else sprintf(strpuf1, "sh -c %s", unterbrechenname);
        } else /* d.h. starttime>endtime, z.B. ueber nacht */ {
            if ((j > starttime) || (j < endtime))
                sprintf(strpuf1, "at %s interruptn%s", argv[i], codename + 5);
            else sprintf(strpuf1, "sh -c %s", unterbrechenname);
        }
        dummy = system(strpuf1);

    }

    if (TIMERESTRICTED) {
        i = TIMERESTRICTED;
        i++;
        starttime = (argv[i][0] - 48)*36000 +
                (argv[i][1] - 48)*3600 + (argv[i][2] - 48)*600
                + (argv[i][3] - 48)*60;
        i++;
        endtime = (argv[i][0] - 48)*36000 +
                (argv[i][1] - 48)*3600 + (argv[i][2] - 48)*600
                + (argv[i][3] - 48)*60;

        sprintf(startname, "start");
        strcat(startname, codename + 5);

        sprintf(unterbrechenname, "interrupt");
        strcat(unterbrechenname, codename + 5);

        sprintf(strpuf1, "startn");
        strcat(strpuf1, codename + 5);
        fil = fopen(strpuf1, "w");
        fprintf(fil, "sh -c %s\n", startname);
        fclose(fil);

        fil = fopen(unterbrechenname, "w");
        fprintf(fil, "exec > at.t.out 2>at.t.errout \n");
        fprintf(fil, "kill -9 %d \n", getpid());
        fprintf(fil, " if test -f %s \n ", lastgraphname);
        fprintf(fil, " then at %s %s \n fi \n", argv[i - 1], strpuf1);
        fclose(fil);
        sprintf(strpuf1, "chmod 755 %s", unterbrechenname);
        dummy = system(strpuf1);

        sprintf(strpuf1, "interruptn");
        strcat(strpuf1, codename + 5);
        fil = fopen(strpuf1, "w");
        fprintf(fil, "sh -c %s\n", unterbrechenname);
        fclose(fil);

        sprintf(strpuf1, " ");
        kommandoname[0] = 0;
        for (j = 0; j < argc; j++) {
            strcat(kommandoname, argv[j]);
            strcat(kommandoname, strpuf1);
        }
        if (!recover) {
            sprintf(strpuf1, "r ");
            strcat(kommandoname, strpuf1);
        }

        fil = fopen(startname, "w");
        fprintf(fil, "exec > at.t.out 2>at.t.errout\n");
        fprintf(fil, "%s \n", kommandoname);
        fclose(fil);
        sprintf(strpuf1, "chmod 755 %s", startname);
        dummy = system(strpuf1);

        b = time(0);
        zeitzeiger = localtime(&b);
        j = (3600 * zeitzeiger->tm_hour) + (60 * zeitzeiger->tm_min)
                + zeitzeiger->tm_sec;

        if (starttime < endtime) {
            if (((j > starttime) && (j < endtime)) ||
                    ((zeitzeiger->tm_wday == 0) || (zeitzeiger->tm_wday == 6))) {
                if ((zeitzeiger->tm_wday == 0) ||
                        (zeitzeiger->tm_wday == 6))
                    sprintf(strpuf1, "at 0000 mon interruptn%s", codename + 5);
                else sprintf(strpuf1, "at %s interruptn%s", argv[i], codename + 5);
            } else sprintf(strpuf1, "sh -c %s", unterbrechenname);
        } else /* d.h. starttime>endtime, z.B. ueber nacht */ {
            if (((j > starttime) || (j < endtime)) ||
                    ((zeitzeiger->tm_wday == 0) || (zeitzeiger->tm_wday == 6))) {
                if ((zeitzeiger->tm_wday == 5) ||
                        (zeitzeiger->tm_wday == 6))
                    sprintf(strpuf1, "at %s mon interruptn%s", argv[i], codename + 5);
                else sprintf(strpuf1, "at %s interruptn%s", argv[i], codename + 5);
            } else sprintf(strpuf1, "sh -c %s", unterbrechenname);
        }
        fprintf(stderr, "%s\n", strpuf1);
        dummy = system(strpuf1);
        dummy = system("at -l");

    }




    if (!recover) {
        if (!singleout_minibaum && !noout_minibaum)
            for (i = startzahl; i < knotenzahl; i += 2) {
                codename[6] = (i / 10) + 48;
                codename[7] = (i % 10) + 48;
                if (!graph6_output_minibaum) {
                    if (init) {
                        fil = fopen(codename, "w");
                        fclose(fil);
                    } else if ((fil = fopen(codename, "r")) != NULL) {
                        fprintf(stderr, "There is a file %s. Is this OK ? (y/n)\n", codename);
                        test = 0;
                        while ((test != 'y') && (test != 'n')) test = getchar();
                        if (test == 'n') exit(0);
                    }
                }
            }

        if (!lastout) {
            codename[6] = (knotenzahl / 10) + 48;
            codename[7] = (knotenzahl % 10) + 48;
            if (!graph6_output_minibaum) {
                if (init) {
                    fil = fopen(codename, "w");
                    fclose(fil);
                } else if ((fil = fopen(codename, "r")) != NULL) {
                    fprintf(stderr, "There is a file %s. Is this OK ? (y/n)\n", codename);
                    test = 0;
                    while ((test != 'y') && (test != 'n')) test = getchar();
                    if (test == 'n') exit(0);
                }
            }
        }
    }





    if (recover) {
        fil = fopen(lastgraphname, "r");
        if (fread(lastgraph, 1, sizeof (GRAPH), fil)<sizeof (GRAPH)) {
            fprintf(stderr, "No lastgraph in file ! Setting recover=0.\n");
            recover = 0;
        }/* d.h.: Da stand noch gar kein graph drin */
        else {
            unsignedintegerlies(&oldtime, fil);
            unsigned_ll_integerlies(&number_of_graphs, fil);
            if (modulo_minibaum) unsignedintegerlies(&splitzaehler, fil);
        }
        fclose(fil);
    }

    switch (taille) {
        case 3: if (knotenzahl < 16) border = knotenzahl - 6;
            else border = knotenzahl - 7;
            break;
        case 4: border = knotenzahl - 7;
            break;
        case 5: border = knotenzahl - 8;
            break;
        case 6: border = knotenzahl - 9;
            break;
        case 7: if (knotenzahl < 28) border = knotenzahl - 12;
            else border = knotenzahl - 13;
            break;
        case 8: border = knotenzahl - 17;
            break;
        case 9: border = knotenzahl - 22;
            break; /* just a guess */
        case 10: border = knotenzahl - 29;
            break; /* just a guess */
    }

    if (bipartite) border--;

    if (numforbcycles >= 3) border -= 5;
    else if (numforbcycles) border -= 2;

    if (vertex_trans) border--; /* wundert mich selbst -- ich haette eher erwartet,
                               dass man border erhoehen muss */

    treeborder = (9 * border) / 10;

    if (numforbcycles >= 3) treeborder = (8 * border) / 10;

    /* PLAY with treeborder and border (always border>=treeborder) between here */

    //border= border-2;
    //treeborder= treeborder-2;

    /* and here */

    if (treeborder > baumgrenze) treeborder = baumgrenze;

    //fprintf(stderr,"%d %d\n",treeborder,border);

    //if (border<=baumgrenze) treeborder=border;
    //                      else treeborder=baumgrenze;



    if (treeborder < 0) treeborder = 0;

    indepborder = border + 2;
    if (indepborder <= indepspeichergrenze) indepstoreborder = indepborder;
    else indepstoreborder = indepspeichergrenze;


    grapherzeugung();

    sprintf(strpuf1, "rm "); /* um nicht noch extra hierfuer
                                            Speicher zu missbrauchen, wird
                                            strpuf1 zweckentfremdet */
    strcat(strpuf1, lastgraphname);

    if ((fil = fopen(lastgraphname, "r")) != NULL) {
        fclose(fil);
        dummy = system(strpuf1);
    }



    if (startvorgabe && (!modulo_minibaum)) {
        sprintf(strpuf1, "rm "); /* um nicht noch extra hierfuer
                                            Speicher zu missbrauchen, wird
                                            strpuf1 zweckentfremdet */
        strcat(strpuf1, startgraphname);

        dummy = system(strpuf1);
    }







    if (timerestricted || TIMERESTRICTED) {
        sprintf(strpuf1, "rm %s %s", startname, unterbrechenname);
        dummy = system(strpuf1);
        sprintf(strpuf1, "rm startn");
        strcat(strpuf1, codename + 5);
        dummy = system(strpuf1);
        sprintf(strpuf1, "rm interruptn");
        strcat(strpuf1, codename + 5);
        dummy = system(strpuf1);
    }

    if (observer) {
        for (i = 1; i <= 6; i++) {
            sprintf(observername, "observer%1d", i);
            strcat(observername, codename + 5);
            sprintf(strpuf1, "rm %s", observername);
            dummy = system(strpuf1);
            sprintf(observername, "observern%1d", i);
            strcat(observername, codename + 5);
            sprintf(strpuf1, "rm %s", observername);
            dummy = system(strpuf1);

        }
    }


    times(&TMS);
    savetime = oldtime + (unsigned int) TMS.tms_utime;

    if (remarks) {
        codename[6] = (knotenzahl / 10) + 48;
        codename[7] = (knotenzahl % 10) + 48;
        /* falls zuletzt eine andere Zahl geschrieben wurde */
        sprintf(strpuf1, "REMARKS_%s", codename + 5);
        fil = fopen(strpuf1, "w");
        for (i = 0; i < argc; i++) fprintf(fil, "%s ", argv[i]);
        fprintf(fil, "\n");
        fprintf(fil, "Erzeugung beendet. Maximal %d Kanten erreicht.\n", maximalkanten);
#ifdef time_factor
        fprintf(fil, "User CPU-Zeit (gesamt): %.1f Sekunden \n", (double) savetime / time_factor);
#endif
        fprintf(fil, "Generated %llu graphs with %d vertices. \n", number_of_graphs, knotenzahl);
        fclose(fil);
    }

    for (i = 0; i < argc; i++) fprintf(stderr, "%s ", argv[i]);
    fprintf(stderr, ": ");
#ifdef time_factor
    fprintf(stderr, "Erzeugung beendet. Maximal %d Kanten erreicht -- %.1f Sekunden.\n"
            , maximalkanten, (double) savetime / time_factor);
#else
    fprintf(stderr, "Erzeugung beendet. Maximal %d Kanten erreicht.\n", maximalkanten);
#endif
    fprintf(stderr, "Generated %llu graphs with %d vertices. \n", number_of_graphs, knotenzahl);


    return (0);

}

