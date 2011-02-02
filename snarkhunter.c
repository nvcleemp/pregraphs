/**
 * snarkhunter.c
 *
 * Snarkhunter: a generator for cubic graphs and snarks.
 *
 * Author: Jan Goedgebeur (Jan.Goedgebeur@UGent.be)
 * In collaboration with Gunnar Brinkmann
 *
 */

/**
 * Compile with:
 *
 * cc -DWORDSIZE=32 -DMAXN=WORDSIZE -O4 snarkhunter.c nautil.c nausparse.c nauty.c -o snarkhunter
 *
 * Or:
 *
 * cc -DWORDSIZE=64 -DMAXN=WORDSIZE -O4 snarkhunter.c nautil.c nausparse.c nauty.c -o snarkhunter
 *
 * -DWORDSIZE=32 is slightly faster, but limits the order of the graphs to 32.
 *
 */

//If -DMAXN=WORDSIZE is omitted MAXN will be 0 in nauty.c (which would also work, but it is slower)

//The following is not sufficient (since MAXN will still be 0 in nauty):
//#ifndef MAXN  /* maximum allowed n value; use 0 for dynamic sizing. */
//#define MAXN WORDSIZE
//#endif

#include <limits.h>
#include <sys/times.h>
#include "nauty/nausparse.h" /* which includes nauty.h */
#include "snarkhunter.h"
#include "util.h"


/****************************Some useful macros********************************/

/**
 * Returns 1 if v2 is a neighbour of v1, else 0.
 */
#define is_neighbour_old(v1, v2) (current_graph[v1][0] == v2 ? 1 : current_graph[v1][1] == v2 ? 1 : current_graph[v1][2] == v2)
#define is_neighbour(v1, v2) ((vertex_neighbourhood[v1] & BIT(v2)) > 0)

/**
 * Increases nauty_calls and calls nauty.
 */
//No difference in timing if usin counter or not
#define nauty_sh(g_arg, lab, ptn, active_arg, orbits_arg, options, stats_arg,\
        ws_arg, worksize, m_arg, n_arg, canong_arg) { nauty_calls++;\
        nauty(g_arg, lab, ptn, active_arg, orbits_arg, options, stats_arg, ws_arg, worksize, m_arg, n_arg, canong_arg);}

/**
 * Replaces neighbour v1 of vertex by v2.
 */
//Using the macro is slower
//#define replace_neighbour(vertex, v1, v2) (current_graph[vertex][0] == v1 ? (current_graph[vertex][0] = v2) : current_graph[vertex][1] == v1 ? (current_graph[vertex][1] = v2) : (current_graph[vertex][2] = v2))
void replace_neighbour(unsigned char vertex, unsigned char v1, unsigned char v2) {
    if(current_graph[vertex][0] == v1)
        current_graph[vertex][0] = v2;
    else if(current_graph[vertex][1] == v1)
        current_graph[vertex][1] = v2;
    else
        current_graph[vertex][2] = v2;
}

/**************Methods for the generation of irreducible graphs****************/

void init_generate_irreducible_graphs() {
    //Startgraph is the K4
    int i, j;
    for(i = 0; i < 4; i++) {
        for(j = 0; j < REG; j++) {
            current_graph[i][j] = (i + j + 1) % 4;
        }
        edge_diamonds[0][i] = i;
    }
    current_number_of_vertices = 4;

    current_number_of_edges = 0;
    edge_labels[0][1] = current_number_of_edges;
    edge_labels[1][0] = current_number_of_edges++;
    edge_labels[0][2] = current_number_of_edges;
    edge_labels[2][0] = current_number_of_edges++;
    edge_labels[0][3] = current_number_of_edges;
    edge_labels[3][0] = current_number_of_edges++;
    edge_labels[1][2] = current_number_of_edges;
    edge_labels[2][1] = current_number_of_edges++;
    edge_labels[1][3] = current_number_of_edges;
    edge_labels[3][1] = current_number_of_edges++;
    edge_labels[2][3] = current_number_of_edges;
    edge_labels[3][2] = current_number_of_edges++;

    number_of_edge_diamonds = 1;

    eligible_edges[0][0] = 0;
    eligible_edges[0][1] = 3;
    eligible_edges_size = 1;

    number_of_nonadj_edge_diamonds = 0;
    number_of_lollipop_diamonds = 0;
}

/**
 * Main method of the generation algorithm.
 * Order is the number of vertices the generated graphs should have.
 * Min_girth can at the moment only be 3, 4 or 5.
 * The last parameter is a pointer to a function which is called when a new graph is generated.
 * This pointer may also be NULL.
 */
void generate_irreducible_graphs(int order, int min_girth, void (*userproc) (unsigned char (*)[REG + 1], int)) {
    DEBUGASSERT(min_girth >= 3 && min_girth <= 5);

    number_of_vertices = order;
    if(girth == 3)
        max_number_of_vertices_irred_graph = order;
    else
        max_number_of_vertices_irred_graph = order - 2;

    if(!noout && output_to_file && output_irred_graphs) {
        sprintf(outputfilename_irred, "Irreducible_graphs.%d", max_number_of_vertices_irred_graph);
        if(graph6_output) {
            strcat(outputfilename_irred, ".g6");
        }

        //To make sure the output appears in a new file instead of being appended to an existing file
        //FILE *file = fopen(outputfilename_irred, "w");
        //fclose(file);
        outputfile_irred = fopen(outputfilename_irred, "w");
        if(outputfile_irred == NULL) {
            fprintf(stderr, "Error: could not create outputfile for irreducible graphs\n");
            exit(1);
        }
    } else {
        outputfile_irred = stdout;
    }

    init_nauty_options();
    calculate_binom_coefficients(number_of_vertices - 1);

    vertexset_index = malloc(sizeof(unsigned int) * MAXN * MAXN * MAXN * MAXN);
    if(vertexset_index == NULL) {
        fprintf(stderr, "Error: Can't get enough memory while creating vertexset_index\n");
        exit(1);
    }

    int i = 0;
    int number_of_edges;
    for(i = 4; i <= number_of_vertices - 2; i += 2) {
        number_of_edges = 3 * i / 2;
        max_edgepairlist_size[i] = number_of_edges * (number_of_edges - 5) / 2;
    }

    max_edgepairlist_size_triangle = (3 * max_number_of_vertices_irred_graph / 2 - 5) * 3;

    if(userproc != NULL) {
        userproc_sh = userproc;
    }

    init_generate_irreducible_graphs();

    extend_irreducible_graph(INIT);

    //Write any remaining graphs to the file
    if(!noout) {
        i = min_order;
        if(singleout)
            i = order;
        while(i <= order) {
            wegspeichern(graphlist, i);
	    graphlist_number_of_graphs[i] = 0;
            i += 2;
        }
    }

    free(vertexset_index);

    SG_FREE(sg);
    SG_FREE(sg_canon);

    if(!noout && output_to_file && output_irred_graphs) {
        fclose(outputfile_irred);
    }

}

//Returns the last edge diamond. This is the diamond which has the border vertex with the largest label the canonical graph
int determine_last_edge_diamond(int lab[]) {
    setword diamonds_border_vertices_bitvector = (setword) 0;
    int i;
    unsigned char neighbour0, neighbour1;
    for(i = 0; i < number_of_edge_diamonds; i++) {
        neighbour0 = determine_external_diamond_neighbour_index(edge_diamonds[i], 0);
        neighbour1 = determine_external_diamond_neighbour_index(edge_diamonds[i], 3);


        //Will never be neighbours since the graph won't contain any reducible edges
        DEBUGASSERT(!is_neighbour_old(neighbour0, neighbour1));

        //Needed because destroying a lollipop can generate an "irreducible" edge diamond
        //reducing that diamond, would yield a loop
        if(neighbour0 != neighbour1)
            diamonds_border_vertices_bitvector |= BIT(edge_diamonds[i][0]) | BIT(edge_diamonds[i][3]);
    }
    for(i = current_number_of_vertices - 1; i >= 0; i--) {
        if((BIT(lab[i]) & diamonds_border_vertices_bitvector) > 0)
            return lab[i];
    }
    fprintf(stderr, "Error: no border vertex found -- edge diamond (can never happen)\n");
    exit(1);

}

int is_major_edge_diamond() {
    if(number_of_edge_diamonds == 1 && current_number_of_vertices == max_number_of_vertices_irred_graph) {
        return 1;
    } else {
        number_of_generators = 0;

        options.getcanon = TRUE;
        options.defaultptn = TRUE;
        copy_sparse_graph();

        nauty_sh((graph*) & sg, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, current_number_of_vertices, (graph*) & sg_canon);

        DEBUGASSERT(number_of_nonadj_edge_diamonds == 0 && number_of_lollipop_diamonds == 0);

        /**
         * TODO: could immediately return 1 if there is only 1 diamond edge
         * or if all diamond edges are in the same orbit. But this won't be
         * much faster since this method certainly isn't the bottleneck.
         *
         * The same goes for lollipops and nonadj edge diamonds.
         */

        int last_edge_diamond_border_vertex = determine_last_edge_diamond(lab);
        int last_diamond_orbit = orbits[last_edge_diamond_border_vertex];

        return orbits[current_number_of_vertices - 4] == last_diamond_orbit || orbits[current_number_of_vertices - 1] == last_diamond_orbit;
    }
}

int determine_last_lollipop(int lab[]) {
    setword lollipop_central_vertices_bitvector = (setword) 0;
    int i;
    unsigned char central_vertex;
    for(i = 0; i < number_of_lollipop_diamonds; i++) {
        central_vertex = determine_external_diamond_neighbour(lollipop_diamonds[i]);
        lollipop_central_vertices_bitvector |= BIT(central_vertex);
    }
    for(i = current_number_of_vertices - 1; i >= 0; i--) {
        if((BIT(lab[i]) & lollipop_central_vertices_bitvector) > 0)
            return lab[i];
    }
    fprintf(stderr, "Error: no central vertex found (can never happen)\n");
    exit(1);

}

int is_major_lollipop_diamond() {
    if(number_of_lollipop_diamonds == 1 && current_number_of_vertices == max_number_of_vertices_irred_graph) {
        return 1;
    } else {
        number_of_generators = 0;

        options.getcanon = TRUE;
        options.defaultptn = TRUE;
        copy_sparse_graph();

        nauty_sh((graph*) & sg, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, current_number_of_vertices, (graph*) & sg_canon);

        int last_lollipop = determine_last_lollipop(lab);
        int last_lollipop_orbit = orbits[last_lollipop];

        return orbits[current_number_of_vertices - 5] == last_lollipop_orbit;
    }
}

int determine_last_nonadj_edge_diamond(int lab[]) {
    setword diamonds_border_vertices_bitvector = (setword) 0;
    /**
     * Remark: would be more efficient and more elegant to go through all nonadjacent
     * edge diamonds and the labels of their extremal vertices and chose the
     * one with the biggest label.
     * But didn't change the code since it hardly uses any cputime (>>0.01%).
     */
    int i;
    for(i = 0; i < number_of_nonadj_edge_diamonds; i++) {
        diamonds_border_vertices_bitvector |= BIT(nonadj_edge_diamonds[i][0]) | BIT(nonadj_edge_diamonds[i][3]);
    }
    for(i = current_number_of_vertices - 1; i >= 0; i--) {
        if((BIT(lab[i]) & diamonds_border_vertices_bitvector) > 0)
            return lab[i];
    }
    fprintf(stderr, "Error no border vertex found -- nonadjacent edge diamond (can never happen)\n");
    exit(1);
}

int is_major_nonadj_edge_diamond() {
    if(number_of_nonadj_edge_diamonds == 1 && current_number_of_vertices == max_number_of_vertices_irred_graph) {
        return 1;
    } else {
        number_of_generators = 0;

        options.getcanon = TRUE;
        options.defaultptn = TRUE;
        copy_sparse_graph();

        nauty_sh((graph*) & sg, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, current_number_of_vertices, (graph*) & sg_canon);

        DEBUGASSERT(number_of_lollipop_diamonds == 0);

        int last_edge_diamond_border_vertex = determine_last_nonadj_edge_diamond(lab);
        int last_diamond_orbit = orbits[last_edge_diamond_border_vertex];

        return orbits[current_number_of_vertices - 5] == last_diamond_orbit || orbits[current_number_of_vertices - 2] == last_diamond_orbit;
    }
}

/**
 * Determines the neighbourhood of each vertex (i.e. the vertices on distance 1).
 */
void calculate_vertex_neighbourhood() {
    unsigned char i, j;
    for(i = 0; i < current_number_of_vertices; i++) {
        vertex_neighbourhood[i] = (setword) 0;
        for(j = 0; j < REG; j++) {
            vertex_neighbourhood[i] |= BIT(current_graph[i][j]);
        }
    }
}

void irreducible_graph_generated() {
    if(output_irred_graphs && !noout && output_to_file)
        aufschreiben_irred();

    //Take backups
    permutation generators_local[number_of_generators][MAXN];
    memcpy(generators_local, generators, sizeof(permutation) * number_of_generators * MAXN);
    int number_of_generators_local = number_of_generators;

    number_of_reducible_triangles = 0;

    //irreducible_triangles_bitvector is already up to date
    fill_list_of_irreducible_triangles();

    //list of bridges are up to date as well

    //Set is_bridge
    //This is not too expensive, since there are hardly any prime graphs
    int j, k;
    for(j = 0; j < number_of_vertices; j++) {
        for(k = j + 1; k < number_of_vertices; k++) {
            is_bridge[j][k] = 0;
        }
    }
    for(j = 0; j < number_of_bridges; j++)
        is_bridge[bridges[j][0]][bridges[j][1]] = 1;

    //Set variables
    previous_min_edge[0] = 0;
    previous_min_edge[1] = 0;
    min_edge[0] = 0;
    min_edge[1] = 0;

    min_colour_one = MAX_EDGE_COLOUR_TWO;
    //Is actually not necessary
    min_edge_is_part_of_square = -1;

    calculate_vertex_neighbourhood();

    extend(INIT, 0);

    //Restore backups
    number_of_generators = number_of_generators_local;
    memcpy(generators, generators_local, sizeof(permutation) * number_of_generators * MAXN);

    //Bridges are already restored by extend
/*
    number_of_bridges = number_of_bridges_local;
    if(number_of_bridges_local > 0)
        memcpy(bridges, bridges_local, sizeof(EDGE) * number_of_bridges_local);
*/
}

void extend_irreducible_graph(int edge_inserted) {
    DEBUGASSERT(current_number_of_vertices <= max_number_of_vertices_irred_graph);

    update_irreducible_triangles_bitvector();

    //Check if graph is irreducible
    int i;
    for(i = 0; i < eligible_edges_size; i++) {
        if(!is_part_of_irreducible_triangle_bitvector(eligible_edges[i][0]) &&
                !is_part_of_irreducible_triangle_bitvector(eligible_edges[i][1]) &&
                !is_a_bridge_list(eligible_edges[i][0], eligible_edges[i][1])) {
            /**
             * I.e. graph contains reducible edges.
             * Only irreducible graphs allowed because the canonical parent of an
             * irreducible graph will always be an irreducible graph. This is because
             * the nonadj diamond operation and the lollipop operation have a higher priority
             * than the edge diamond operation.
             * Only the nonadj diamond operation and the lollipop operation can
             * yield reducible edges
             */

            return;
        }
    }

    //Check if last operation was canonical
    if(edge_inserted == EDGE_DIAMOND_INSERTED) {
        if(!is_major_edge_diamond()) {
            return;
        }
    } else if(edge_inserted == LOLLIPOP_DIAMOND_INSERTED) {
        if(!is_major_lollipop_diamond()) {
            return;
        }
    } else if(edge_inserted == NONADJ_EDGE_DIAMOND_INSERTED) {
        if(!is_major_nonadj_edge_diamond()) {
            return;
        }
    } else if(edge_inserted == INIT) {
        number_of_generators = 0;

        options.getcanon = FALSE;
        options.defaultptn = TRUE;
        copy_sparse_graph();

        nauty_sh((graph*) & sg, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, current_number_of_vertices, (graph*) & sg_canon);
    }

    irreducible_graph_generated();

    if(current_number_of_vertices <= max_number_of_vertices_irred_graph - 4) {
        //Make copy of generators because they will be modified by edge_extend and triangle_extend
        permutation generators_local[number_of_generators][MAXN];
        memcpy(generators_local, generators, sizeof(permutation) * number_of_generators * MAXN);
        int number_of_generators_local = number_of_generators;

        //Diamond edges
        EDGE eligible_edges_list[eligible_edges_size];
        int eligible_edges_list_size = 0;
        generate_eligible_diamond_edges(eligible_edges_list, &eligible_edges_list_size);

        if(eligible_edges_list_size > 0) {
            edge_diamond_extend(eligible_edges_list, eligible_edges_list_size);

            number_of_generators = number_of_generators_local;
            memcpy(generators, generators_local, sizeof(permutation) * number_of_generators * MAXN);
        }

        if(current_number_of_vertices <= max_number_of_vertices_irred_graph - 6) {
            //Nonadj edges
            //List is actually a lot smaller, but is no problem since there are hardly any prime graphs
            EDGEPAIR edge_pairs_list[max_edgepairlist_size[current_number_of_vertices]];
            int edge_pairs_list_size;
            generate_non_adjacent_diamond_edge_pairs(edge_pairs_list, &edge_pairs_list_size);

            if(edge_pairs_list_size > 0) {
                nonadj_edge_diamond_extend(edge_pairs_list, edge_pairs_list_size);

                number_of_generators = number_of_generators_local;
                memcpy(generators, generators_local, sizeof(permutation) * number_of_generators * MAXN);
            }

            //Lollipops
            eligible_edges_list_size = 0;
            generate_eligible_lollipop_edges(eligible_edges_list, &eligible_edges_list_size);

            if(eligible_edges_list_size > 0) {
                edge_lollipop_extend(eligible_edges_list, eligible_edges_list_size);
            }
        }
    }
}

void add_edge_diamond_to_list(unsigned char v0, unsigned char v1, unsigned char v2, unsigned char v3) {
    edge_diamonds[number_of_edge_diamonds][0] = v0;
    edge_diamonds[number_of_edge_diamonds][1] = v1;
    edge_diamonds[number_of_edge_diamonds][2] = v2;
    edge_diamonds[number_of_edge_diamonds][3] = v3;
    number_of_edge_diamonds++;
}

void remove_edge_diamond_from_list(int index) {
    DEBUGASSERT(number_of_edge_diamonds > 0);
    number_of_edge_diamonds--;
    if(number_of_edge_diamonds > 0 && index < number_of_edge_diamonds) {
        edge_diamonds[index][0] = edge_diamonds[number_of_edge_diamonds][0];
        edge_diamonds[index][1] = edge_diamonds[number_of_edge_diamonds][1];
        edge_diamonds[index][2] = edge_diamonds[number_of_edge_diamonds][2];
        edge_diamonds[index][3] = edge_diamonds[number_of_edge_diamonds][3];
    }
}

void add_eligible_edge(unsigned char from, unsigned char to) {
    eligible_edges[eligible_edges_size][0] = from;
    eligible_edges[eligible_edges_size][1] = to;
    eligible_edges_size++;
}

void replace_eligible_edge(int old_from, int old_to, int from, int to) {
    int i;
    for(i = 0; i < eligible_edges_size; i++) {
        if(eligible_edges[i][0] == old_from && eligible_edges[i][1] == old_to) {
            eligible_edges[i][0] = from;
            eligible_edges[i][1] = to;
            break;
        }
    }
    if(i >= eligible_edges_size) {
        fprintf(stderr, "Error: eligible edge not found\n");
        exit(1);
    }
}

/**
 * Adds an edge diamonds. After adding this diamond, all lollipops or nonadj edges
 * should be destroyed.
 */
void add_edge_diamond(EDGE edge) {
    current_graph[current_number_of_vertices][0] = edge[0];
    current_graph[current_number_of_vertices][1] = current_number_of_vertices + 1;
    current_graph[current_number_of_vertices][2] = current_number_of_vertices + 2;

    current_graph[current_number_of_vertices + 1][0] = current_number_of_vertices;
    current_graph[current_number_of_vertices + 1][1] = current_number_of_vertices + 2;
    current_graph[current_number_of_vertices + 1][2] = current_number_of_vertices + 3;

    current_graph[current_number_of_vertices + 2][0] = current_number_of_vertices;
    current_graph[current_number_of_vertices + 2][1] = current_number_of_vertices + 1;
    current_graph[current_number_of_vertices + 2][2] = current_number_of_vertices + 3;

    current_graph[current_number_of_vertices + 3][0] = edge[1];
    current_graph[current_number_of_vertices + 3][1] = current_number_of_vertices + 1;
    current_graph[current_number_of_vertices + 3][2] = current_number_of_vertices + 2;


    replace_neighbour(edge[0], edge[1], current_number_of_vertices);
    replace_neighbour(edge[1], edge[0], current_number_of_vertices + 3);

    //Update labels (needed to calculate edgepair orbits)
    //Recycle old label
    edge_labels[edge[0]][current_number_of_vertices] = edge_labels[edge[0]][edge[1]];
    edge_labels[current_number_of_vertices][edge[0]] = edge_labels[edge[0]][edge[1]];

    //New labels
    edge_labels[edge[1]][current_number_of_vertices + 3] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 3][edge[1]] = current_number_of_edges++;

    edge_labels[current_number_of_vertices][current_number_of_vertices + 1] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 1][current_number_of_vertices] = current_number_of_edges++;

    edge_labels[current_number_of_vertices][current_number_of_vertices + 2] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 2][current_number_of_vertices] = current_number_of_edges++;

    edge_labels[current_number_of_vertices + 1][current_number_of_vertices + 2] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 2][current_number_of_vertices + 1] = current_number_of_edges++;

    edge_labels[current_number_of_vertices + 1][current_number_of_vertices + 3] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 3][current_number_of_vertices + 1] = current_number_of_edges++;

    edge_labels[current_number_of_vertices + 2][current_number_of_vertices + 3] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 3][current_number_of_vertices + 2] = current_number_of_edges++;

    //Update bridges
    if(is_a_bridge_list(edge[0], edge[1])) {
        replace_bridge(edge[0], edge[1], edge[0], current_number_of_vertices);
        add_bridge(edge[1], current_number_of_vertices + 3);
    }

    add_edge_diamond_to_list(current_number_of_vertices, current_number_of_vertices + 1, current_number_of_vertices + 2, current_number_of_vertices + 3);

    //Update eligible edges
    replace_eligible_edge(edge[0], edge[1], edge[0], current_number_of_vertices);
    add_eligible_edge(edge[1], current_number_of_vertices + 3);

    //Update lollipop en nonadj diamond edges (zou nu allemaal 0 moeten zijn)
    DEBUGASSERT(number_of_nonadj_edge_diamonds <= 1 && number_of_lollipop_diamonds <= 2);
    if(number_of_nonadj_edge_diamonds == 1) {
        add_edge_diamond_to_list(nonadj_edge_diamonds[0][0], nonadj_edge_diamonds[0][1], nonadj_edge_diamonds[0][2], nonadj_edge_diamonds[0][3]);
        number_of_nonadj_edge_diamonds = 0;
    } else if(number_of_lollipop_diamonds <= 2) {
        int i;
        for(i = 0; i < number_of_lollipop_diamonds; i++)
            add_edge_diamond_to_list(lollipop_diamonds[i][0], lollipop_diamonds[i][1], lollipop_diamonds[i][2], lollipop_diamonds[i][3]);
        number_of_lollipop_diamonds = 0;
    }

    //Adding a an edge diamond can never yield a lollipop or nonadj edge diamond
    //update_edge_diamonds();

    current_number_of_vertices += 4;

}

void remove_diamond_edge(EDGE edge) {
    replace_neighbour(edge[0], current_number_of_vertices - 4, edge[1]);
    replace_neighbour(edge[1], current_number_of_vertices - 1, edge[0]);

    /**
     * Important!: the edge_labels of the original edges do not need to be restored, because
     * they still have their original value. This is because there will never be a
     * new edge added which is identical to the orginal edge.
     */

    current_number_of_vertices -= 4;
    current_number_of_edges -= 6;
}

void add_lollipop_to_list(unsigned char v0, unsigned char v1, unsigned char v2, unsigned char v3) {
    lollipop_diamonds[number_of_lollipop_diamonds][0] = v0;
    lollipop_diamonds[number_of_lollipop_diamonds][1] = v1;
    lollipop_diamonds[number_of_lollipop_diamonds][2] = v2;
    lollipop_diamonds[number_of_lollipop_diamonds][3] = v3;
    number_of_lollipop_diamonds++;
}

void remove_lollipop_from_list(int index) {
    number_of_lollipop_diamonds--;
    if(number_of_lollipop_diamonds > 0 && index < number_of_lollipop_diamonds) {
        int i;
        for(i = 0; i < 4; i++)
            lollipop_diamonds[index][i] = lollipop_diamonds[number_of_lollipop_diamonds][i];
    }
}

/**
 * The insertion of a lollipop may destroy an existing lollipop. This method will
 * detected if an existing lollipop will be destroyed by inserting a new lollipop in "edge".
 * And if this is the case, this existing lollipop will also be removed from the list
 * of lollipops and will be added to the list of edge diamonds.
 */
void update_existing_lollipops(EDGE edge) {
    int i;
    for(i = 0; i < number_of_lollipop_diamonds; i++) {
        if(edge[0] == lollipop_diamonds[i][0] || edge[0] == lollipop_diamonds[i][3] ||
                edge[1] == lollipop_diamonds[i][0] || edge[1] == lollipop_diamonds[i][3]) {
            add_edge_diamond_to_list(lollipop_diamonds[i][0], lollipop_diamonds[i][1], lollipop_diamonds[i][2], lollipop_diamonds[i][3]);
            remove_lollipop_from_list(i);
            break;
        }
    }
}

void add_lollipop_edge(EDGE edge) {
    replace_neighbour(edge[0], edge[1], current_number_of_vertices);
    replace_neighbour(edge[1], edge[0], current_number_of_vertices);

    current_graph[current_number_of_vertices][0] = edge[0];
    current_graph[current_number_of_vertices][1] = edge[1];
    current_graph[current_number_of_vertices][2] = current_number_of_vertices + 1;

    current_graph[current_number_of_vertices + 1][0] = current_number_of_vertices;
    current_graph[current_number_of_vertices + 1][1] = current_number_of_vertices + 2;
    current_graph[current_number_of_vertices + 1][2] = current_number_of_vertices + 5;

    current_graph[current_number_of_vertices + 2][0] = current_number_of_vertices + 1;
    current_graph[current_number_of_vertices + 2][1] = current_number_of_vertices + 3;
    current_graph[current_number_of_vertices + 2][2] = current_number_of_vertices + 4;

    current_graph[current_number_of_vertices + 3][0] = current_number_of_vertices + 2;
    current_graph[current_number_of_vertices + 3][1] = current_number_of_vertices + 4;
    current_graph[current_number_of_vertices + 3][2] = current_number_of_vertices + 5;

    current_graph[current_number_of_vertices + 4][0] = current_number_of_vertices + 2;
    current_graph[current_number_of_vertices + 4][1] = current_number_of_vertices + 3;
    current_graph[current_number_of_vertices + 4][2] = current_number_of_vertices + 5;

    current_graph[current_number_of_vertices + 5][0] = current_number_of_vertices + 1;
    current_graph[current_number_of_vertices + 5][1] = current_number_of_vertices + 3;
    current_graph[current_number_of_vertices + 5][2] = current_number_of_vertices + 4;

    //Update labels (needed to calculate edgepair orbits)
    //One old label and 9 new labels
    //Recycle old label
    edge_labels[edge[0]][current_number_of_vertices] = edge_labels[edge[0]][edge[1]];
    edge_labels[current_number_of_vertices][edge[0]] = edge_labels[edge[0]][edge[1]];

    //New labels
    edge_labels[edge[1]][current_number_of_vertices] = current_number_of_edges;
    edge_labels[current_number_of_vertices][edge[1]] = current_number_of_edges++;

    edge_labels[current_number_of_vertices][current_number_of_vertices + 1] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 1][current_number_of_vertices] = current_number_of_edges++;

    edge_labels[current_number_of_vertices + 1][current_number_of_vertices + 2] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 2][current_number_of_vertices + 1] = current_number_of_edges++;

    edge_labels[current_number_of_vertices + 1][current_number_of_vertices + 5] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 5][current_number_of_vertices + 1] = current_number_of_edges++;

    edge_labels[current_number_of_vertices + 2][current_number_of_vertices + 3] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 3][current_number_of_vertices + 2] = current_number_of_edges++;

    edge_labels[current_number_of_vertices + 2][current_number_of_vertices + 4] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 4][current_number_of_vertices + 2] = current_number_of_edges++;

    edge_labels[current_number_of_vertices + 3][current_number_of_vertices + 4] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 4][current_number_of_vertices + 3] = current_number_of_edges++;

    edge_labels[current_number_of_vertices + 3][current_number_of_vertices + 5] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 5][current_number_of_vertices + 3] = current_number_of_edges++;

    edge_labels[current_number_of_vertices + 4][current_number_of_vertices + 5] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 5][current_number_of_vertices + 4] = current_number_of_edges++;

    //Update bridges
    if(is_a_bridge_list(edge[0], edge[1])) {
        replace_bridge(edge[0], edge[1], edge[0], current_number_of_vertices);
        add_bridge(edge[1], current_number_of_vertices);
    }
    add_bridge(current_number_of_vertices, current_number_of_vertices + 1);

    add_lollipop_to_list(current_number_of_vertices + 2, current_number_of_vertices + 3, current_number_of_vertices + 4, current_number_of_vertices + 5);

    update_existing_lollipops(edge);

    update_edge_diamonds();

    //Update eligible edges
    replace_eligible_edge(edge[0], edge[1], edge[0], current_number_of_vertices);
    add_eligible_edge(edge[1], current_number_of_vertices);
    add_eligible_edge(current_number_of_vertices, current_number_of_vertices + 1);
    add_eligible_edge(current_number_of_vertices + 1, current_number_of_vertices + 2);
    add_eligible_edge(current_number_of_vertices + 1, current_number_of_vertices + 5);

    current_number_of_vertices += 6;
    //current_number_of_edges += 9;
}

void remove_lollipop_edge(EDGE edge) {
    replace_neighbour(edge[0], current_number_of_vertices - 6, edge[1]);
    replace_neighbour(edge[1], current_number_of_vertices - 6, edge[0]);

    current_number_of_vertices -= 6;
    current_number_of_edges -= 9;
}

unsigned char determine_external_diamond_neighbour(IRRED_TRIANGLE diamond) {
    int i;
    for(i = 0; i < REG; i++) {
        if(current_graph[diamond[0]][i] != diamond[1] && current_graph[diamond[0]][i] != diamond[2])
            return current_graph[diamond[0]][i];
    }
    fprintf(stderr, "Error: no external neighbour found (should never happen)\n");
    exit(1);
}

unsigned char determine_external_diamond_neighbour_index(IRRED_TRIANGLE diamond, int index) {
    DEBUGASSERT(index == 0 || index == 3);
    int i;
    for(i = 0; i < REG; i++) {
        if(current_graph[diamond[index]][i] != diamond[1] && current_graph[diamond[index]][i] != diamond[2])
            return current_graph[diamond[index]][i];
    }
    fprintf(stderr, "Error: no external neighbour found -- index (should never happen)\n");
    exit(1);
}

/**
 * Returns the neighbour of a vertex which is not an extremal vertex of a diamond.
 */
unsigned char determine_nondiamond_neighbour(IRRED_TRIANGLE diamond, unsigned char vertex) {
    int i;
    for(i = 0; i < REG; i++) {
        if(current_graph[vertex][i] != diamond[0] && current_graph[vertex][i] != diamond[3])
            return current_graph[vertex][i];
    }
    fprintf(stderr, "Error: nondiamond neighbour found\n");
    exit(1);
}

void generate_eligible_diamond_edges(EDGE eligible_diamond_edges[], int *eligible_diamond_edges_size) {
    *eligible_diamond_edges_size = 0;
    if(number_of_nonadj_edge_diamonds == 0 && number_of_lollipop_diamonds == 0) {
        *eligible_diamond_edges_size = eligible_edges_size;
        memcpy(eligible_diamond_edges, eligible_edges, sizeof(EDGE) * eligible_edges_size);
    } else if(number_of_nonadj_edge_diamonds == 1 && number_of_lollipop_diamonds == 0) {
        int i;
        for(i = 0; i < 2; i++) {
            eligible_diamond_edges[*eligible_diamond_edges_size][0] = determine_external_diamond_neighbour_index(nonadj_edge_diamonds[0], 3*i);
            eligible_diamond_edges[*eligible_diamond_edges_size][1] = nonadj_edge_diamonds[0][3*i];
            transform_edge_into_canonical_form(eligible_diamond_edges[*eligible_diamond_edges_size]);
            (*eligible_diamond_edges_size)++;
        }
    } else if(number_of_lollipop_diamonds == 1 && number_of_nonadj_edge_diamonds == 0) {
        unsigned char central_vertex = determine_external_diamond_neighbour(lollipop_diamonds[0]);
        int i;
        for(i = 0; i < 2; i++) {
            eligible_diamond_edges[*eligible_diamond_edges_size][0] = central_vertex;
            eligible_diamond_edges[*eligible_diamond_edges_size][1] = lollipop_diamonds[0][3*i];
            DEBUGASSERT(eligible_diamond_edges[*eligible_diamond_edges_size][0] < eligible_diamond_edges[*eligible_diamond_edges_size][1]);
            (*eligible_diamond_edges_size)++;
        }
        DEBUGASSERT(determine_external_diamond_neighbour_index(lollipop_diamonds[0], 0) == determine_external_diamond_neighbour_index(lollipop_diamonds[0], 3));

        //Different edge (which will always be in a different orbit)
        eligible_diamond_edges[*eligible_diamond_edges_size][0] = central_vertex;
        eligible_diamond_edges[*eligible_diamond_edges_size][1] = determine_nondiamond_neighbour(lollipop_diamonds[0], central_vertex);
        transform_edge_into_canonical_form(eligible_diamond_edges[*eligible_diamond_edges_size]);
        (*eligible_diamond_edges_size)++;

        DEBUGASSERT(*eligible_diamond_edges_size == 3);

    } else if(number_of_lollipop_diamonds == 2 && number_of_nonadj_edge_diamonds == 0) {
        //Only possibility to destroy 2 lollipops is if their central vertices are neighbours
        //(Only happens in one case where there are 10 vertices)
        unsigned char central_vertex0 = determine_external_diamond_neighbour(lollipop_diamonds[0]);
        unsigned char central_vertex1 = determine_external_diamond_neighbour(lollipop_diamonds[1]);
        if(is_neighbour_old(central_vertex0, central_vertex1)) {
            DEBUGASSERT(current_number_of_vertices == 10);
            eligible_diamond_edges[0][0] = central_vertex0;
            eligible_diamond_edges[0][1] = central_vertex1;
            transform_edge_into_canonical_form(eligible_diamond_edges[0]);
            *eligible_diamond_edges_size = 1;
        }
    }
    //Else no eligible edge diamonds

    //Set the edge_index, needed to calculate the orbits!
    if(*eligible_diamond_edges_size > 0) {
        int i;
        for(i = 0; i < *eligible_diamond_edges_size; i++) {
            DEBUGASSERT(eligible_diamond_edges[i][0] < eligible_diamond_edges[i][1]);
            edge_index[eligible_diamond_edges[i][0]][eligible_diamond_edges[i][1]] = i;
        }
    }

}

void generate_eligible_lollipop_edges(EDGE eligible_lollipop_edges[], int *eligible_lollipop_edges_size) {
    *eligible_lollipop_edges_size = eligible_edges_size;
    memcpy(eligible_lollipop_edges, eligible_edges, sizeof(EDGE) * eligible_edges_size);

    if(*eligible_lollipop_edges_size > 0) {
        int i;
        for(i = 0; i < *eligible_lollipop_edges_size; i++) {
            DEBUGASSERT(eligible_lollipop_edges[i][0] < eligible_lollipop_edges[i][1]);
            edge_index[eligible_lollipop_edges[i][0]][eligible_lollipop_edges[i][1]] = i;
        }
    }
}

void edge_diamond_extend(EDGE eligible_diamond_edges[], int eligible_diamond_edges_size) {
    int edge_diamond_orbits[eligible_diamond_edges_size+1];
    int number_of_edge_diamond_orbits = 0;

    int is_trivial_group = 0;
    if(number_of_generators > 0 && eligible_diamond_edges_size > 1) {
        //The generators are still ok, because edge_extension is done before triangle extension
        determine_edge_orbits(eligible_diamond_edges, eligible_diamond_edges_size, edge_diamond_orbits, &number_of_edge_diamond_orbits);
    } else {
        is_trivial_group = 1;
        number_of_edge_diamond_orbits = eligible_diamond_edges_size;
    }

    //Take backups
    int number_of_bridges_local = number_of_bridges;
    EDGE bridges_local[number_of_bridges+1];
    if(number_of_bridges_local > 0)
        memcpy(bridges_local, bridges, sizeof(EDGE) * number_of_bridges_local);

    EDGE eligible_edges_local[eligible_edges_size+1];
    int eligible_edges_size_local = eligible_edges_size;
    if(eligible_edges_size_local > 0)
        memcpy(eligible_edges_local, eligible_edges, sizeof(EDGE) * eligible_edges_size_local);

    IRRED_TRIANGLE edge_diamonds_local[number_of_edge_diamonds+1];
    int number_of_edge_diamonds_local = number_of_edge_diamonds;
    if(number_of_edge_diamonds_local > 0)
        memcpy(edge_diamonds_local, edge_diamonds, sizeof(IRRED_TRIANGLE) * number_of_edge_diamonds_local);

    IRRED_TRIANGLE lollipop_diamonds_local[number_of_lollipop_diamonds+1];
    int number_of_lollipop_diamonds_local = number_of_lollipop_diamonds;
    if(number_of_lollipop_diamonds_local > 0)
        memcpy(lollipop_diamonds_local, lollipop_diamonds, sizeof(IRRED_TRIANGLE) * number_of_lollipop_diamonds);

    IRRED_TRIANGLE nonadj_edge_diamonds_local[number_of_nonadj_edge_diamonds+1];
    int number_of_nonadj_edge_diamonds_local = number_of_nonadj_edge_diamonds;
    if(number_of_nonadj_edge_diamonds_local > 0)
        memcpy(nonadj_edge_diamonds_local, nonadj_edge_diamonds, sizeof(IRRED_TRIANGLE) * number_of_nonadj_edge_diamonds);

    int i;
    int num_orbits = 0;
    for(i = 0; i < eligible_diamond_edges_size; i++) {
        if(is_trivial_group || edge_diamond_orbits[i] == i) {
            add_edge_diamond(eligible_diamond_edges[i]);
            DEBUGASSERT(number_of_nonadj_edge_diamonds == 0 && number_of_lollipop_diamonds == 0);

            extend_irreducible_graph(EDGE_DIAMOND_INSERTED);

            remove_diamond_edge(eligible_diamond_edges[i]);

            //Restore lists
            number_of_bridges = number_of_bridges_local;
            if(number_of_bridges_local > 0)
                memcpy(bridges, bridges_local, sizeof(EDGE) * number_of_bridges_local);

            eligible_edges_size = eligible_edges_size_local;
            if(eligible_edges_size_local > 0)
                memcpy(eligible_edges, eligible_edges_local, sizeof(EDGE) * eligible_edges_size_local);

            number_of_edge_diamonds = number_of_edge_diamonds_local;
            if(number_of_edge_diamonds > 0) {
                memcpy(edge_diamonds, edge_diamonds_local, sizeof (IRRED_TRIANGLE) * number_of_edge_diamonds);
            }

            number_of_lollipop_diamonds = number_of_lollipop_diamonds_local;
            if(number_of_lollipop_diamonds > 0) {
                memcpy(lollipop_diamonds, lollipop_diamonds_local, sizeof (IRRED_TRIANGLE) * number_of_lollipop_diamonds);
            }

            number_of_nonadj_edge_diamonds = number_of_nonadj_edge_diamonds_local;
            if(number_of_nonadj_edge_diamonds > 0) {
                memcpy(nonadj_edge_diamonds, nonadj_edge_diamonds_local, sizeof (IRRED_TRIANGLE) * number_of_nonadj_edge_diamonds);
            }

            num_orbits++;
            if(num_orbits == number_of_edge_diamond_orbits)
                break;
        }
    }
    DEBUGASSERT(num_orbits == number_of_edge_diamond_orbits);
}

void edge_lollipop_extend(EDGE eligible_lollipop_edges[], int eligible_lollipop_edges_size) {
    int edge_diamond_orbits[eligible_lollipop_edges_size+1];
    int number_of_edge_diamond_orbits = 0;

    int is_trivial_group = 0;
    if(number_of_generators > 0 && eligible_lollipop_edges_size > 1) {
        //The generators are still ok, because edge_extension is done before triangle extension
        determine_edge_orbits(eligible_lollipop_edges, eligible_lollipop_edges_size, edge_diamond_orbits, &number_of_edge_diamond_orbits);
    } else {
        is_trivial_group = 1;
        number_of_edge_diamond_orbits = eligible_lollipop_edges_size;
    }

    int number_of_bridges_local = number_of_bridges;
    EDGE bridges_local[number_of_bridges+1];
    if(number_of_bridges_local > 0)
        memcpy(bridges_local, bridges, sizeof(EDGE) * number_of_bridges_local);

    EDGE eligible_edges_local[eligible_edges_size+1];
    int eligible_edges_size_local = eligible_edges_size;
    if(eligible_edges_size_local > 0)
        memcpy(eligible_edges_local, eligible_edges, sizeof(EDGE) * eligible_edges_size_local);

    IRRED_TRIANGLE edge_diamonds_local[number_of_edge_diamonds+1];
    int number_of_edge_diamonds_local = number_of_edge_diamonds;
    if(number_of_edge_diamonds_local > 0)
        memcpy(edge_diamonds_local, edge_diamonds, sizeof(IRRED_TRIANGLE) * number_of_edge_diamonds_local);

    IRRED_TRIANGLE lollipop_diamonds_local[number_of_lollipop_diamonds+1];
    int number_of_lollipop_diamonds_local = number_of_lollipop_diamonds;
    if(number_of_lollipop_diamonds_local > 0)
        memcpy(lollipop_diamonds_local, lollipop_diamonds, sizeof(IRRED_TRIANGLE) * number_of_lollipop_diamonds);

    IRRED_TRIANGLE nonadj_edge_diamonds_local[number_of_nonadj_edge_diamonds+1];
    int number_of_nonadj_edge_diamonds_local = number_of_nonadj_edge_diamonds;
    if(number_of_nonadj_edge_diamonds_local > 0)
        memcpy(nonadj_edge_diamonds_local, nonadj_edge_diamonds, sizeof(IRRED_TRIANGLE) * number_of_nonadj_edge_diamonds);

    int i;
    int num_orbits = 0;
    for(i = 0; i < eligible_lollipop_edges_size; i++) {
        if(is_trivial_group || edge_diamond_orbits[i] == i) {
            add_lollipop_edge(eligible_lollipop_edges[i]);

            extend_irreducible_graph(LOLLIPOP_DIAMOND_INSERTED);

            remove_lollipop_edge(eligible_lollipop_edges[i]);

            //Restore lists
            number_of_bridges = number_of_bridges_local;
            if(number_of_bridges_local > 0)
                memcpy(bridges, bridges_local, sizeof(EDGE) * number_of_bridges_local);

            eligible_edges_size = eligible_edges_size_local;
            if(eligible_edges_size_local > 0)
                memcpy(eligible_edges, eligible_edges_local, sizeof(EDGE) * eligible_edges_size_local);

            number_of_edge_diamonds = number_of_edge_diamonds_local;
            if(number_of_edge_diamonds > 0) {
                memcpy(edge_diamonds, edge_diamonds_local, sizeof (IRRED_TRIANGLE) * number_of_edge_diamonds);
            }

            number_of_lollipop_diamonds = number_of_lollipop_diamonds_local;
            if(number_of_lollipop_diamonds > 0) {
                memcpy(lollipop_diamonds, lollipop_diamonds_local, sizeof (IRRED_TRIANGLE) * number_of_lollipop_diamonds);
            }

            number_of_nonadj_edge_diamonds = number_of_nonadj_edge_diamonds_local;
            if(number_of_nonadj_edge_diamonds > 0) {
                memcpy(nonadj_edge_diamonds, nonadj_edge_diamonds_local, sizeof (IRRED_TRIANGLE) * number_of_nonadj_edge_diamonds);
            }

            num_orbits++;
            if(num_orbits == number_of_edge_diamond_orbits)
                break;
        }
    }
    DEBUGASSERT(num_orbits == number_of_edge_diamond_orbits);
}

/**
 * Fully sets irreducible_triangles_bitvector.
 */
void update_irreducible_triangles_bitvector() {
    irreducible_triangles_bitvector = 0;
    int i, j;
    for(i = 0; i < number_of_edge_diamonds; i++) {
        for(j = 0; j < 4; j++) {
            irreducible_triangles_bitvector |= BIT(edge_diamonds[i][j]);
        }
    }
    for(i = 0; i < number_of_lollipop_diamonds; i++) {
        for(j = 0; j < 4; j++) {
            irreducible_triangles_bitvector |= BIT(lollipop_diamonds[i][j]);
        }
    }
    for(i = 0; i < number_of_nonadj_edge_diamonds; i++) {
        for(j = 0; j < 4; j++) {
            irreducible_triangles_bitvector |= BIT(nonadj_edge_diamonds[i][j]);
        }
    }
}

/**
 * Update the edge diamonds.
 * It's possible that an edge diamond becomes a nonadj edge diamond after the insertion
 * of a lollipop or a nonadj edge.
 */
void update_edge_diamonds() {
    update_irreducible_triangles_bitvector();

    int previous_i = 0;
    int abort = 0;
    int i = 0;
    int j, k;
    EDGE external_vertices;
    EDGE temp_edge;
    setword external_neighbours[2];
    while(i < number_of_edge_diamonds) {
        for(i = previous_i; i < number_of_edge_diamonds; i++) {
            for(j = 0; j < 2; j++) {
                external_vertices[j] = determine_external_diamond_neighbour_index(edge_diamonds[i], 3 * j);

                temp_edge[0] = edge_diamonds[i][3 * j];
                temp_edge[1] = external_vertices[j];
                transform_edge_into_canonical_form(temp_edge);
                if((BIT(external_vertices[j]) & irreducible_triangles_bitvector) > 0 ||
                        is_a_bridge_list(temp_edge[0], temp_edge[1])) { //Not part of a diamond and not a bridge
                    abort = 1;
                    break;
                }
                external_neighbours[j] = BIT(external_vertices[j]);
                for(k = 0; k < REG; k++) {
                    external_neighbours[j] |= BIT(current_graph[external_vertices[j]][k]);
                }
            }
            if(!abort) {
                //Check if edgepairs are nonadjacent
                if((external_neighbours[0] & external_neighbours[1]) == 0) {
                    //Add nonadj edge diamond
                    add_nonadj_edge_diamond_to_list(edge_diamonds[i][0], edge_diamonds[i][1], edge_diamonds[i][2], edge_diamonds[i][3]);

                    //Remove edge diamond
                    remove_edge_diamond_from_list(i);

                    previous_i = i;
                    break;
                } else if(external_vertices[0] == external_vertices[1]) {
                    //Might be a lollipop now
                    unsigned char central_vertex_neighbour = determine_nondiamond_neighbour(edge_diamonds[i], external_vertices[0]);
                    //Neighbour of central vertex isn't part of a diamond
                    if((BIT(central_vertex_neighbour) & irreducible_triangles_bitvector) == 0) {
                        //Add lollipop
                        add_lollipop_to_list(edge_diamonds[i][0], edge_diamonds[i][1], edge_diamonds[i][2], edge_diamonds[i][3]);

                        //Remove edge diamond
                        remove_edge_diamond_from_list(i);

                        previous_i = i;
                        break;
                    }
                }
            } else
                abort = 0;
        }
    }
}

void add_nonadj_edge_diamond_to_list(unsigned char v0, unsigned char v1, unsigned char v2, unsigned char v3) {
    nonadj_edge_diamonds[number_of_nonadj_edge_diamonds][0] = v0;
    nonadj_edge_diamonds[number_of_nonadj_edge_diamonds][1] = v1;
    nonadj_edge_diamonds[number_of_nonadj_edge_diamonds][2] = v2;
    nonadj_edge_diamonds[number_of_nonadj_edge_diamonds][3] = v3;
    number_of_nonadj_edge_diamonds++;
}

/**
 * Returns 0 if there were any lollipops generated by adding this nonadj diamond edge
 * (in this case the nonadj diamond edge can't be canonical), else returns 1.
 */
int add_nonadj_diamond_edge(EDGEPAIR edge_pair) {
    DEBUGASSERT(edge_pair[2] < current_number_of_vertices && edge_pair[3] < current_number_of_vertices);

    current_graph[current_number_of_vertices][0] = edge_pair[0];
    current_graph[current_number_of_vertices][1] = edge_pair[1];
    current_graph[current_number_of_vertices][2] = current_number_of_vertices + 1;

    current_graph[current_number_of_vertices + 1][0] = current_number_of_vertices;
    current_graph[current_number_of_vertices + 1][1] = current_number_of_vertices + 2;
    current_graph[current_number_of_vertices + 1][2] = current_number_of_vertices + 3;

    current_graph[current_number_of_vertices + 2][0] = current_number_of_vertices + 1;
    current_graph[current_number_of_vertices + 2][1] = current_number_of_vertices + 3;
    current_graph[current_number_of_vertices + 2][2] = current_number_of_vertices + 4;

    current_graph[current_number_of_vertices + 3][0] = current_number_of_vertices + 1;
    current_graph[current_number_of_vertices + 3][1] = current_number_of_vertices + 2;
    current_graph[current_number_of_vertices + 3][2] = current_number_of_vertices + 4;

    current_graph[current_number_of_vertices + 4][0] = current_number_of_vertices + 2;
    current_graph[current_number_of_vertices + 4][1] = current_number_of_vertices + 3;
    current_graph[current_number_of_vertices + 4][2] = current_number_of_vertices + 5;

    current_graph[current_number_of_vertices + 5][0] = edge_pair[2];
    current_graph[current_number_of_vertices + 5][1] = edge_pair[3];
    current_graph[current_number_of_vertices + 5][2] = current_number_of_vertices + 4;

    replace_neighbour(edge_pair[0], edge_pair[1], current_number_of_vertices);
    replace_neighbour(edge_pair[1], edge_pair[0], current_number_of_vertices);

    replace_neighbour(edge_pair[2], edge_pair[3], current_number_of_vertices + 5);
    replace_neighbour(edge_pair[3], edge_pair[2], current_number_of_vertices + 5);

    //Update labels (needed to calculate edgepair orbits) (2 recycled and 9 new)
    //Recycle old labels
    edge_labels[edge_pair[0]][current_number_of_vertices] = edge_labels[edge_pair[0]][edge_pair[1]];
    edge_labels[current_number_of_vertices][edge_pair[0]] = edge_labels[edge_pair[0]][edge_pair[1]];

    edge_labels[edge_pair[1]][current_number_of_vertices] = edge_labels[edge_pair[2]][edge_pair[3]];
    edge_labels[current_number_of_vertices][edge_pair[1]] = edge_labels[edge_pair[2]][edge_pair[3]];

    //New labels
    edge_labels[current_number_of_vertices][current_number_of_vertices + 1] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 1][current_number_of_vertices] = current_number_of_edges++;

    edge_labels[current_number_of_vertices + 1][current_number_of_vertices + 2] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 2][current_number_of_vertices + 1] = current_number_of_edges++;

    edge_labels[current_number_of_vertices + 1][current_number_of_vertices + 3] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 3][current_number_of_vertices + 1] = current_number_of_edges++;

    edge_labels[current_number_of_vertices + 2][current_number_of_vertices + 3] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 3][current_number_of_vertices + 2] = current_number_of_edges++;

    edge_labels[current_number_of_vertices + 2][current_number_of_vertices + 4] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 4][current_number_of_vertices + 2] = current_number_of_edges++;

    edge_labels[current_number_of_vertices + 3][current_number_of_vertices + 4] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 4][current_number_of_vertices + 3] = current_number_of_edges++;

    edge_labels[current_number_of_vertices + 4][current_number_of_vertices + 5] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 5][current_number_of_vertices + 4] = current_number_of_edges++;

    edge_labels[edge_pair[2]][current_number_of_vertices + 5] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 5][edge_pair[2]] = current_number_of_edges++;

    edge_labels[edge_pair[3]][current_number_of_vertices + 5] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 5][edge_pair[3]] = current_number_of_edges++;

    /**
     * Could return 0 if a reducible edge was generated.
     * But this is checked at extend_irreducible_graph() and the generation
     * of prime graphs is certainly not the bottleneck.
     */
    /* Update bridges */
    int res = contains_bridge(edge_pair);
    if(res == 0) {
        //update_bridges_add_edge();
    } else if(res == 1 || res == 2) {
        res--;
        if(is_a_bridge(edge_pair[2*res], current_number_of_vertices + 5*res))
            replace_bridge(edge_pair[2*res], edge_pair[2*res + 1], edge_pair[2*res], current_number_of_vertices + 5*res);
        else
            replace_bridge(edge_pair[2*res], edge_pair[2*res + 1], edge_pair[2*res + 1], current_number_of_vertices + 5*res);
    } else {
        //Either (edge_pair[0], new) or (edge_pair[1], new) is now a bridge
        if(is_a_bridge(edge_pair[0], current_number_of_vertices))
            replace_bridge(edge_pair[0], edge_pair[1], edge_pair[0], current_number_of_vertices);
        else
            replace_bridge(edge_pair[0], edge_pair[1], edge_pair[1], current_number_of_vertices);

        //Either (edge_pair[2], new) or (edge_pair[3], new) is now a bridge
        if(is_a_bridge(edge_pair[2], current_number_of_vertices + 5))
            replace_bridge(edge_pair[2], edge_pair[3], edge_pair[2], current_number_of_vertices + 5);
        else
            replace_bridge(edge_pair[2], edge_pair[3], edge_pair[3], current_number_of_vertices + 5);
    }

    update_bridges_add_edge();

    //Adding a nonadj edge diamond can destroy at most 2 lollipops
    if(number_of_lollipop_diamonds > 0) {
        DEBUGASSERT(number_of_lollipop_diamonds <= 2);
        int i;
        for(i = 0; i < number_of_lollipop_diamonds; i++)
            add_edge_diamond_to_list(lollipop_diamonds[i][0], lollipop_diamonds[i][1], lollipop_diamonds[i][2], lollipop_diamonds[i][3]);
        number_of_lollipop_diamonds = 0;
    }

    add_nonadj_edge_diamond_to_list(current_number_of_vertices + 1, current_number_of_vertices + 2, current_number_of_vertices + 3, current_number_of_vertices + 4);

    //Of beter voor de lollipop doen, want die zullen toch nooit nonadj zijn?
    update_edge_diamonds();

    //Update eligible edges
    replace_eligible_edge(edge_pair[0], edge_pair[1], edge_pair[0], current_number_of_vertices);
    add_eligible_edge(edge_pair[1], current_number_of_vertices);
    add_eligible_edge(current_number_of_vertices, current_number_of_vertices + 1);

    replace_eligible_edge(edge_pair[2], edge_pair[3], edge_pair[2], current_number_of_vertices + 5);
    add_eligible_edge(edge_pair[3], current_number_of_vertices + 5);
    add_eligible_edge(current_number_of_vertices + 4, current_number_of_vertices + 5);


    current_number_of_vertices += 6;

    return number_of_lollipop_diamonds == 0;
}

void remove_nonadj_diamond_edge(EDGEPAIR edge_pair) {
    replace_neighbour(edge_pair[0], current_number_of_vertices - 6, edge_pair[1]);
    replace_neighbour(edge_pair[1], current_number_of_vertices - 6, edge_pair[0]);
    replace_neighbour(edge_pair[2], current_number_of_vertices - 1, edge_pair[3]);
    replace_neighbour(edge_pair[3], current_number_of_vertices - 1, edge_pair[2]);

    current_number_of_vertices -= 6;
    current_number_of_edges -= 9;
}

/**
 * Fills the list of irreducible_triangles with the edge_diamonds, lollipop_diamonds
 * and nonadj_edge_diamonds.
 */
void fill_list_of_irreducible_triangles() {
    number_of_irreducible_triangles = 0;
    int i, j;
    for(i = 0; i < number_of_edge_diamonds; i++) {
        for(j = 0; j < 4; j++)
            irreducible_triangles[number_of_irreducible_triangles][j] = edge_diamonds[i][j];
        number_of_irreducible_triangles++;
    }
    for(i = 0; i < number_of_lollipop_diamonds; i++) {
        for(j = 0; j < 4; j++)
            irreducible_triangles[number_of_irreducible_triangles][j] = lollipop_diamonds[i][j];
        number_of_irreducible_triangles++;
    }
    for(i = 0; i < number_of_nonadj_edge_diamonds; i++) {
        for(j = 0; j < 4; j++)
            irreducible_triangles[number_of_irreducible_triangles][j] = nonadj_edge_diamonds[i][j];
        number_of_irreducible_triangles++;
    }
    DEBUGASSERT(number_of_irreducible_triangles == number_of_edge_diamonds + number_of_lollipop_diamonds + number_of_nonadj_edge_diamonds);
}

void generate_non_adjacent_diamond_edge_pairs(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size) {
    *edge_pair_list_size = 0;
    /**
     * If number_of_lollipop_diamonds == 2, the lollipops can be destroyed,
     * but the extended graph will contain reducible edges.
     */
    if(number_of_lollipop_diamonds <= 1) {
        fill_list_of_irreducible_triangles();
        if(number_of_lollipop_diamonds == 0) {
            generate_all_nonadj_diamond_edge_pairs(edge_pairs_list, edge_pair_list_size);
        } else { //number_of_lollipop_diamonds == 1
            generate_nonadj_diamond_edge_pairs_one_lollipop(edge_pairs_list, edge_pair_list_size);
        }
    }
}

void generate_nonadj_diamond_edge_pairs_one_lollipop(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size) {
    DEBUGASSERT(number_of_lollipop_diamonds == 1);
    unsigned char central_vertex = determine_external_diamond_neighbour(lollipop_diamonds[0]);

    int i, j, k, next, diamond;
    for(i = 0; i < current_number_of_vertices - 1; i++) {
        if(i != central_vertex) {
            for(j = 0; j < REG; j++) {
                next = current_graph[i][j];
                if(i < next && next != central_vertex) {
                    if(!(is_part_of_irreducible_triangle_diamond(i, &diamond) && is_part_of_same_irreducible_triangle(next, diamond))) {
                        for(k = 0; k < 2; k++) {
                            edge_pairs_list[*edge_pair_list_size][0] = central_vertex;
                            edge_pairs_list[*edge_pair_list_size][1] = lollipop_diamonds[0][3*k];
                            edge_pairs_list[*edge_pair_list_size][2] = i;
                            edge_pairs_list[*edge_pair_list_size][3] = next;

                            transform_edgepair_into_canonical_form(edge_pairs_list[*edge_pair_list_size]);
                            edgepair_index[edge_labels[edge_pairs_list[*edge_pair_list_size][0]][edge_pairs_list[*edge_pair_list_size][1]]][edge_labels[edge_pairs_list[*edge_pair_list_size][2]][edge_pairs_list[*edge_pair_list_size][3]]] = *edge_pair_list_size;
                            (*edge_pair_list_size)++;
                        }
                    }
                }
            }
        }
    }
}

void generate_all_nonadj_diamond_edge_pairs(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size) {
    DEBUGASSERT(number_of_lollipop_diamonds == 0);

    //Generate the edgepairs that contain no edge which is fully in a diamond
    int diamond;
    int i, j, next, k, l, next0;
    for(i = 0; i < current_number_of_vertices - 3; i++) {
        for(j = 0; j < REG; j++) {
            next = current_graph[i][j];
            if(i < next) {
                if(!(is_part_of_irreducible_triangle_diamond(i, &diamond) && is_part_of_same_irreducible_triangle(next, diamond))) {
                    for(k = i + 1; k < current_number_of_vertices - 1; k++) {
                        if(k != next) { // k != i is automatically implied because k >= i+1
                            for(l = 0; l < REG; l++) {
                                next0 = current_graph[k][l];
                                if(k < next0 && next0 != next) {
                                    if(!(is_part_of_irreducible_triangle_diamond(k, &diamond) && is_part_of_same_irreducible_triangle(next0, diamond))) {
                                        //Add edgepair i next k next0
                                        edge_pairs_list[*edge_pair_list_size][0] = i;
                                        edge_pairs_list[*edge_pair_list_size][1] = next;
                                        edge_pairs_list[*edge_pair_list_size][2] = k;
                                        edge_pairs_list[*edge_pair_list_size][3] = next0;

                                        int index0 = edge_labels[i][next];
                                        int index1 = edge_labels[k][next0];

                                        edgepair_index[index0][index1] = *edge_pair_list_size;
                                        (*edge_pair_list_size)++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void nonadj_edge_diamond_extend(EDGEPAIR edge_pairs_list[], int edge_pair_list_size) {
    int edgepair_orbits[edge_pair_list_size+1];
    int number_of_edgepair_orbits = 0;

    int is_trivial_group = 0;
    if(number_of_generators > 0 && edge_pair_list_size > 1) {
        determine_edgepair_orbits(edge_pairs_list, edge_pair_list_size, edgepair_orbits, &number_of_edgepair_orbits);
    } else {
        is_trivial_group = 1;
        number_of_edgepair_orbits = edge_pair_list_size;
    }

    int number_of_bridges_local = number_of_bridges;
    EDGE bridges_local[number_of_bridges+1];
    if(number_of_bridges_local > 0)
        memcpy(bridges_local, bridges, sizeof(EDGE) * number_of_bridges_local);

    EDGE eligible_edges_local[eligible_edges_size+1];
    int eligible_edges_size_local = eligible_edges_size;
    if(eligible_edges_size_local > 0)
        memcpy(eligible_edges_local, eligible_edges, sizeof(EDGE) * eligible_edges_size_local);

    IRRED_TRIANGLE edge_diamonds_local[number_of_edge_diamonds+1];
    int number_of_edge_diamonds_local = number_of_edge_diamonds;
    if(number_of_edge_diamonds_local > 0)
        memcpy(edge_diamonds_local, edge_diamonds, sizeof(IRRED_TRIANGLE) * number_of_edge_diamonds_local);

    IRRED_TRIANGLE lollipop_diamonds_local[number_of_lollipop_diamonds+1];
    int number_of_lollipop_diamonds_local = number_of_lollipop_diamonds;
    if(number_of_lollipop_diamonds_local > 0)
        memcpy(lollipop_diamonds_local, lollipop_diamonds, sizeof(IRRED_TRIANGLE) * number_of_lollipop_diamonds);

    IRRED_TRIANGLE nonadj_edge_diamonds_local[number_of_nonadj_edge_diamonds+1];
    int number_of_nonadj_edge_diamonds_local = number_of_nonadj_edge_diamonds;
    if(number_of_nonadj_edge_diamonds_local > 0)
        memcpy(nonadj_edge_diamonds_local, nonadj_edge_diamonds, sizeof(IRRED_TRIANGLE) * number_of_nonadj_edge_diamonds);

    int i;
    int num_orbits = 0;
    for(i = 0; i < edge_pair_list_size; i++) {
        if(is_trivial_group || edgepair_orbits[i] == i) {
            if(add_nonadj_diamond_edge(edge_pairs_list[i])) {
                DEBUGASSERT(number_of_lollipop_diamonds == 0);
                extend_irreducible_graph(NONADJ_EDGE_DIAMOND_INSERTED);
            }

            remove_nonadj_diamond_edge(edge_pairs_list[i]);

            //Restore lists
            number_of_bridges = number_of_bridges_local;
            if(number_of_bridges_local > 0)
                memcpy(bridges, bridges_local, sizeof(EDGE) * number_of_bridges_local);

            eligible_edges_size = eligible_edges_size_local;
            if(eligible_edges_size_local > 0)
                memcpy(eligible_edges, eligible_edges_local, sizeof(EDGE) * eligible_edges_size_local);

            number_of_edge_diamonds = number_of_edge_diamonds_local;
            if(number_of_edge_diamonds > 0) {
                memcpy(edge_diamonds, edge_diamonds_local, sizeof (IRRED_TRIANGLE) * number_of_edge_diamonds);
            }

            number_of_lollipop_diamonds = number_of_lollipop_diamonds_local;
            if(number_of_lollipop_diamonds > 0) {
                memcpy(lollipop_diamonds, lollipop_diamonds_local, sizeof (IRRED_TRIANGLE) * number_of_lollipop_diamonds);
            }

            number_of_nonadj_edge_diamonds = number_of_nonadj_edge_diamonds_local;
            if(number_of_nonadj_edge_diamonds > 0) {
                memcpy(nonadj_edge_diamonds, nonadj_edge_diamonds_local, sizeof (IRRED_TRIANGLE) * number_of_nonadj_edge_diamonds);
            }

            num_orbits++;
            if(num_orbits == number_of_edgepair_orbits)
                break;
        }
    }
    DEBUGASSERT(num_orbits == number_of_edgepair_orbits);
}

/************End methods for the generation of irreducible graphs**************/

/**
 * Determines partitions for nauty in case of EDGE_INSERTED.
 * Vertices with the same colour are put in the same partition.
 */
void determine_vertex_partitions() {
    int i, j;

    int min_colour = MAX_VERTEX_COLOUR_DISTANCE_TWO;
    int max_colour = 0;

    int min_colour_marked = MAX_VERTEX_COLOUR_DISTANCE_THREE;
    int max_colour_marked = 0;

    int number_not_marked = 0;
    unsigned char vertex_colours_local[current_number_of_vertices];
    for(i = 0; i < current_number_of_vertices; i++) {
        vertex_colours_local[i] = 0;
        if(!ISMARKED(i)) {
            for(j = 0; j < REG; j++) {
                if(i < current_graph[i][j]) {
                    vertex_colours_local[i] += colours_one[i][current_graph[i][j]];
                } else {
                    vertex_colours_local[i] += colours_one[current_graph[i][j]][i];
                }
            }
            if(vertex_colours_local[i] < min_colour)
                min_colour = vertex_colours_local[i];
            if(vertex_colours_local[i] > max_colour)
                max_colour = vertex_colours_local[i];

            number_not_marked++;
        } else {
            for(j = 0; j < REG; j++) {
                /**
                 * Remark: if current_graph[i][j] is not marked,
                 * colours_two[i][current_graph[i][j]] will always be 0.
                 */
                if(i < current_graph[i][j]) {
                    vertex_colours_local[i] += colours_two[i][current_graph[i][j]];
                } else {
                    vertex_colours_local[i] += colours_two[current_graph[i][j]][i];
                }
            }
            //vertex_colours_local[i] = POPCOUNT(vertex_colours_long[i]);
            if(vertex_colours_local[i] < min_colour_marked)
                min_colour_marked = vertex_colours_local[i];
            if(vertex_colours_local[i] > max_colour_marked)
                max_colour_marked = vertex_colours_local[i];
        }
    }

    //Always at least one marked, namely the inserted edge
    DEBUGASSERT(max_colour_marked >= min_colour_marked);
    DEBUGASSERT(number_not_marked < current_number_of_vertices);

    DEBUGASSERT(max_colour <= MAX_VERTEX_COLOUR_DISTANCE_TWO);
    DEBUGASSERT(max_colour_marked <= MAX_VERTEX_COLOUR_DISTANCE_THREE);

    RESETMARKS_VERTEX_COLOUR_2;
    RESETMARKS_VERTEX_COLOUR_3;
    unsigned char colour;
    for(i = 0; i < current_number_of_vertices; i++) {
        colour = vertex_colours_local[i];
        if(!ISMARKED(i)) {
            if(!ISMARKED_VERTEX_COLOUR_2(colour)) {
                vertex_colours_one_size[colour] = 0;
                MARK_VERTEX_COLOUR_2(colour);
            }
            vertex_colours_one[colour][vertex_colours_one_size[colour]++] = i;
        } else {
            if(!ISMARKED_VERTEX_COLOUR_3(colour)) {
                vertex_colours_marked_size[colour] = 0;
                MARK_VERTEX_COLOUR_3(colour);
            }
            vertex_colours_marked[colour][vertex_colours_marked_size[colour]++] = i;
        }
    }

    //Important: the partitions must be "in order", otherwise there may be isomorphic copies
    //But the order of the vertices in the same partition doesn't matter

    //First determine the partitions of the "nonminimal" vertices
    int lab_index = 0;
    for(i = min_colour; i <= max_colour && lab_index < number_not_marked; i++) {
        if(ISMARKED_VERTEX_COLOUR_2(i)) { //i.e. vertex_colours_one_size[i] > 0
            DEBUGASSERT(vertex_colours_one_size[i] > 0);
            for(j = 0; j < vertex_colours_one_size[i]; j++) {
                lab[lab_index] = vertex_colours_one[i][j];
                ptn[lab_index] = 1;
                lab_index++;
            }
            ptn[lab_index - 1] = 0;
        }
    }

    //Now determine the partitions of the marked vertices
    for(i = min_colour_marked; i <= max_colour_marked && lab_index < current_number_of_vertices; i++) {
        if(ISMARKED_VERTEX_COLOUR_3(i)) { //i.e. vertex_colours_marked_size[i] > 0
            DEBUGASSERT(vertex_colours_marked_size[i] > 0);
            for(j = 0; j < vertex_colours_marked_size[i]; j++) {
                lab[lab_index] = vertex_colours_marked[i][j];
                ptn[lab_index] = 1;
                lab_index++;
            }
            ptn[lab_index - 1] = 0;
        }
    }

    DEBUGASSERT(lab_index == current_number_of_vertices);

    options.defaultptn = FALSE;

}

/**
 * Partitions the vertices of the graph using colours_one.
 *
 * Important: This should only be used in case of MAJOR_EDGE_INSERTED.
 * Because in this method the resulting partitions are not sorted.
 * But this is no problem since MAJOR_EDGE_INSERTED is always canonical.
 */
//Max edge colour is 14 = 1110, so shift should be at least 4
//#define PARTITION_SHIFT 4
void determine_major_edge_partitions() {
    //DEBUGASSERT(options.getcanon == FALSE);

    int i, j;
    int neighbour;
    int colour_list[current_number_of_vertices];
    for(i = 0; i < current_number_of_vertices; i++) {
        colour_list[i] = 0;
        for(j = 0; j < REG; j++) {
            neighbour = current_graph[i][j];
            if(i < neighbour) //Is slightly faster than UNION
                colour_list[i] += colours_one[i][neighbour];
            else
                colour_list[i] += colours_one[neighbour][i];
        }
    }

    //Is not any faster (or slower):
/*
    int i, j;
    int neighbour;
    int colour_list[current_number_of_vertices];
    TRIANGLE current_colour;
    for(i = 0; i < current_number_of_vertices; i++) {
        for(j = 0; j < REG; j++) {
            neighbour = current_graph[i][j];
            if(i < neighbour)
                current_colour[j] = colours_one[i][neighbour];
            else
                current_colour[j] = colours_one[neighbour][i];
        }
        transform_triangle_into_canonical_form_full(current_colour);

        //Now determine the vertex colour
        //Max edge colour is 14 = 1110, so shift should be at least 4
        colour_list[i] = current_colour[0];
        for(j = 1; j < 3; j++) {
            UNION(colour_list[i], current_colour[j] << (j * PARTITION_SHIFT));
        }
    }
*/

    int lab_index = 0;
    RESETMARKS;
    for(i = 0; i < current_number_of_vertices && lab_index < current_number_of_vertices; i++) {
        if(!ISMARKED(i)) {
            MARK(i);
            lab[lab_index] = i;
            ptn[lab_index] = 1;
            lab_index++;
            for(j = i + 1; j < current_number_of_vertices; j++)
                if(colour_list[j] == colour_list[i]) {
                    DEBUGASSERT(!ISMARKED(j));
                    MARK(j);
                    lab[lab_index] = j;
                    ptn[lab_index] = 1;
                    lab_index++;
                }
            ptn[lab_index - 1] = 0;
        }
    }

    DEBUGASSERT(lab_index == current_number_of_vertices);

    options.defaultptn = FALSE;
}

/**
 * Calculates the partitions in case of (triangle inserted && !trivial_group)
 * or in case of major_edge_inserted and calls nauty after the partitions have
 * been calculated.
 */
void call_nauty_major_edge_or_triangle_inserted(int edge_inserted, int trivial_group) {
    if(edge_inserted == TRIANGLE_INSERTED && !trivial_group) {
        number_of_generators = 0;
        options.getcanon = FALSE;

        //ptn should already have been set by triangle_extend
        DEBUGASSERT(options.defaultptn == FALSE);

        copy_sparse_graph();

        nauty_sh((graph*) & sg, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, current_number_of_vertices, NULL);
    } else if(edge_inserted == MAJOR_EDGE_INSERTED) {
        options.getcanon = FALSE;
        determine_major_edge_partitions();
        //Calling determine_vertex_partitions instead is not faster
        //Calling determine_vertex_colours_shiher is a lot slower!
        //determine_vertex_partitions();

        number_of_generators = 0;

        copy_sparse_graph();

        nauty_sh((graph*) & sg, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, current_number_of_vertices, NULL);
    }
}

/**
 * In case of snarks, nauty isn't called on level n-2 and n.
 * So afterwards one has to check the generated snark children to make sure
 * no isomorphic snarks are output.
 */
void output_nonisomorphic_children(int children_can_be_isomorphic) {
    DEBUGASSERT(current_number_of_vertices == number_of_vertices - 2);

    if(children_can_be_isomorphic && graphlist_snarks_size > 1) {
        /**
         * TODO: Could use a resizable array of marks for this, but this method is
         * only very rarely called and doesn't cost much.
         */
        int disabled[graphlist_snarks_size];
        int i;
        for(i = 0; i < graphlist_snarks_size; i++) {
            disabled[i] = 0;
        }

        //Check graphs pairwise to see if their canonical form is the same
        int j;
        for(i = 0; i < graphlist_snarks_size; i++) {
            if(!disabled[i]) {
                for(j = i + 1; j < graphlist_snarks_size; j++) {
                    if(!disabled[j]) {
                        //compare canon i en j
                        //necessary to call aresame_sg, because the order of the neighbours might be different
                        if(aresame_sg(&canon_graphs[i], &canon_graphs[j])) {
                            disabled[j] = 1;
                        }
                    }
                }
            }
        }
        for(i = 0; i < graphlist_snarks_size; i++) {
            if(!disabled[i]) {
                aufschreiben_already_filtered(graphlist_snarks[i], number_of_vertices);
            }
        }
    } else {
        int i;
        for(i = 0; i < graphlist_snarks_size; i++)
            aufschreiben_already_filtered(graphlist_snarks[i], number_of_vertices);
    }
}

/**
 * The main method of the generation algorithm.
 * This methods checks if the current graph is canonical and determines the
 * possible extensions for it.
 */
void extend(int edge_inserted, int trivial_group) {
    if(modulo && current_number_of_vertices == splitlevel) {
        //Is cheaper than just doing mod
        splitlevel_counter++;
        if(splitlevel_counter == mod)
            splitlevel_counter = 0;
        if(splitlevel_counter != rest)
            return;
    }

    if(girth > 3) {
        //Some extra bounding criterion
        if(current_number_of_vertices == number_of_vertices - 2) {
            if(number_of_reducible_triangles == 1 && number_of_irreducible_triangles == 1)
                return;
            if(number_of_irreducible_triangles == 2)
                return;
        }

        int vertices_required = ((number_of_irreducible_triangles + number_of_reducible_triangles + 1) / 2) * 2;
        if(girth > 4)
            vertices_required *= 2;
        if(current_number_of_vertices + vertices_required > number_of_vertices) {
            return;
        }
    }

    /* First check if the current graph is canonical */
    still_has_to_check_if_graph_is_canonical = 0;
    if(edge_inserted == EDGE_INSERTED) {
        //In case of snarks, nauty isnt called on level >= n - 2
        if(!snarks || current_number_of_vertices < number_of_vertices - 2) {
            //Could check if can generate valid g5 graph here if g5 and cv == n - 2
            //But this isn't any faster (thanks to number_of_edges_part_of_square in has_min_colour_cycle)

            if(!is_major_edge(0, &sg_canon)) {
                return;
            } else {
                //The newly inserted edge  was (a) minimal edge
                if(current_number_of_vertices < number_of_vertices) {
                    min_edge[0] = current_number_of_vertices - 2;
                    min_edge[1] = current_number_of_vertices - 1;
                }
            }
        } else { //snarks && current_number_of_vertices >= n - 2
            still_has_to_check_if_graph_is_canonical = 1;
        }
    } else if(edge_inserted == INIT) {

/*
        if(modulo && rest != 0 && current_number_of_vertices > splitlevel) {
            return;
        }
*/

        /**
         * To make sure no duplicate graphs are generated when using modulo,
         * irred startgraphs with >= splitlevel vertices are only accepted
         * if rest == splitlevel_counter (this will spread the graphs more evenly
         * than when using rest == 0)
         */
        if(modulo && current_number_of_vertices > splitlevel && splitlevel_counter != rest) {
            return;
        }

        //Only needs to check in case of INIT, in case of edge or triangle extend it was
        //already checked during extension
        if((check_cyclic_connectivity && min_cyclic_connectivity > 1)) {
            int vertices_required = (min_cyclic_connectivity - 1) * 2 * MIN(number_of_bridges, 1);
            if(current_number_of_vertices + vertices_required > number_of_vertices) {
                return;
            }
        }

        if(current_number_of_vertices < number_of_vertices) {
            //Necessary because generators of irreducible graph aren't always calculated because of optimization:
            //e.g. if on last irred level and just 1 lollipop, then nauty isn't called

            number_of_generators = 0;

            options.getcanon = FALSE;
            options.defaultptn = TRUE;
            copy_sparse_graph();

            nauty_sh((graph*) & sg, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, current_number_of_vertices, (graph*) & sg_canon);
        }

    }

    DEBUGASSERT(current_number_of_vertices <= number_of_vertices);

    /* Graph is canonical */
    if(current_number_of_vertices == number_of_vertices) {
        aufschreiben();
    } else {
        if(!singleout) {
            if(!still_has_to_check_if_graph_is_canonical) {
                aufschreiben();
            } else {
                DEBUGASSERT(edge_inserted == EDGE_INSERTED && snarks && current_number_of_vertices == number_of_vertices - 2);
                if(!is_major_edge(0, &sg_canon)) {
                    return;
                } else {
                    //Could first check if it can yield a valid g5 graph
                    //This helps significantly, but the code becomes ugly
                    //and the singleout case then also becomes slightly slower (due to optimizer?)

                    aufschreiben();
                    //Could reuse generates but this doesn't help much
                }
            }
        }

        /* Determine the possible extension */

        EDGEPAIR edge_pairs_list[max_edgepairlist_size[current_number_of_vertices]];
        int edge_pair_list_size;

        DEBUGASSERT(current_number_of_vertices < number_of_vertices);

        if(girth > 3 && current_number_of_vertices == number_of_vertices - 2) {
            if(girth == 4) {
                if(!generate_edgepairs_penultimate_level_girth4(edge_pairs_list, &edge_pair_list_size))
                    return;
            } else if(girth == 5) {
                //Already taken care of by lookahead
                DEBUGASSERT(number_of_reducible_triangles == 0 && number_of_irreducible_triangles == 0);
                if(!generate_edgepairs_penultimate_level_girth5_no_triangles(edge_pairs_list, &edge_pair_list_size)) {
                    return;
                } else {
                    edgepair_list_index_squares = edge_pair_list_size;
                }
            }

            //edge_pair_list_size > 0 if generate_edgepairs_penultimate_level returns 1.
            DEBUGASSERT(edge_pair_list_size > 0 && edge_pair_list_size <= max_edgepairlist_size[current_number_of_vertices]);

            if(!snarks) {
                call_nauty_major_edge_or_triangle_inserted(edge_inserted, trivial_group);

                if(trivial_group) {
                    number_of_generators = 0;
                }
            } else {
                //Expand all edgepairs (act as if it's a trivial group)
                //In most cases the group will be trivial anyways
                number_of_generators = 0;
            }

            if(snarks && edge_inserted == EDGE_INSERTED) {
                //Calling determine_vertex_partitions here isn't too expensive, since it doesn't happen very often
                determine_vertex_partitions();
                memcpy(lab_parent, lab, sizeof(int) * current_number_of_vertices);
                memcpy(ptn_parent, ptn, sizeof(int) * current_number_of_vertices);

                memcpy(colours_three_parent, colours_three, sizeof(unsigned char) * current_number_of_vertices * MAXN);
                min_colour_three_parent = min_colour_three;

                memcpy(edgelist_parent, edgelist, sizeof(EDGE) * edgelist_size);
                edgelist_parent_size = edgelist_size;
                memcpy(edge_index_parent, edge_index, sizeof(unsigned char) * current_number_of_vertices * MAXN);
            }

            graphlist_snarks_size = 0;

            //No need to copy generators if only edge_extend is done
            edge_extend(edge_pairs_list, edge_pair_list_size);

            //Will only be > 0 in case of snarks
            if(graphlist_snarks_size > 0) {
                //Check if parent was canonical
                if(edge_inserted == EDGE_INSERTED) {
                    memcpy(lab, lab_parent, sizeof(int) * current_number_of_vertices);
                    memcpy(ptn, ptn_parent, sizeof(int) * current_number_of_vertices);
                    options.defaultptn = FALSE;

                    memcpy(colours_three, colours_three_parent, sizeof(unsigned char) * current_number_of_vertices * MAXN);
                    min_colour_three = min_colour_three_parent;

                    edgelist_size = edgelist_parent_size;
                    memcpy(edgelist, edgelist_parent, sizeof(EDGE) * edgelist_size);

                    memcpy(edge_index, edge_index_parent, sizeof(unsigned char) * current_number_of_vertices * MAXN);
                    if(!is_major_edge(1, &sg_canon)) {
                        return;
                    }
                    trivial_group = number_of_generators == 0;
                } // else parent was certainly canonical

                //Could also calculate the group in case of MAJOR_EDGE_INSERTED (will often be trivial)
                //But this isn't any faster, since output_nonisomorphic_children is very cheap (only rarely called)

                //Output the nonisomorph children.
                //If the group was trivial, the children will certainly be nonisomorphic
                output_nonisomorphic_children(!trivial_group);
            }

        } else if(girth == 3 || current_number_of_vertices < number_of_vertices - 2) {
            //Always needs to call nauty even if edge_pair_list is empty, because triangle_extend is always called
            call_nauty_major_edge_or_triangle_inserted(edge_inserted, trivial_group);

            find_edge_pairs(edge_pairs_list, &edge_pair_list_size);

            DEBUGASSERT(edge_pair_list_size >= 0 && edge_pair_list_size <= max_edgepairlist_size[current_number_of_vertices]);

            //Make copy of generators because they will be modified by edge_extend and triangle_extend
            permutation generators_local[number_of_generators+1][MAXN];
            int vertex_orbits_local[current_number_of_vertices];
            int number_of_vertex_orbits;
            int groupsize;
            if(!trivial_group) {
                memcpy(generators_local, generators, sizeof(permutation) * number_of_generators * MAXN);

                memcpy(vertex_orbits_local, orbits, sizeof(int) * current_number_of_vertices);
                number_of_vertex_orbits = stats.numorbits;

                //stats.grpsize2 should normally always be 0, but just checking to be sure
                if(stats.grpsize2 == 0) {
                    if(stats.grpsize1 > MAXVAL) {
                        fprintf(stderr, "Error: groupsize doesn't fit in an integer (should never happen)\n");
                        exit(1);
                    }
                    //int groupsize_no_round = stats.grpsize1;
                    groupsize = stats.grpsize1 + 0.1;
                    //Should always be the same
/*
                    if(groupsize != groupsize_no_round) {
                        fprintf(stderr, "Error: groupsize != groupsize_no_round: %d vs %d\n", groupsize, groupsize_no_round);
                        exit(1);
                    }
*/
                } else {
                    fprintf(stderr, "Warning: groupsize2 != 0: %d\n", stats.grpsize2);
                    groupsize = stats.grpsize1 * power(10, stats.grpsize2) + 0.1;
                }
            } else {
                number_of_generators = 0;
                DEBUGASSERT(edge_inserted == TRIANGLE_INSERTED);
            }
            int number_of_generators_local = number_of_generators;

            if(edge_pair_list_size != 0)
                edge_extend(edge_pairs_list, edge_pair_list_size);
            triangle_extend_all(generators_local, number_of_generators_local, vertex_orbits_local, number_of_vertex_orbits, groupsize);
        }
    }

}

/**
 * Last_edge will be set to the reducible edge with min_colour_three which has
 * the biggest label in the canonical graph.
 *
 * Remark: lastedge will be in canonical form.
 */
void determine_last_edge(sparsegraph sparse_graph_canon, int lab[], EDGE lastedge) {
    /**
     * Remark: first the biggest label is determined and only afterwards one tests
     * if it's in the list of minimal edges. This approach is very inefficient.
     * It would be a lot more efficient only to check the labels of the edges
     * with min colour (i.e. the edges in "edgelist") and chose the one with the
     * largest label from it.
     *
     * But didn't change this approach since it only consumes very few cpu
     * according to the profiler (0.2% for 24 3 and descending).
     */

    int i, j, neighbour;
    TRIANGLE neighbours;

    for(i = current_number_of_vertices - 1; i >= 1; i--) {
        if(!is_part_of_irreducible_triangle_bitvector(lab[i])) {
            for(j = 0; j < REG; j++) {
                neighbours[j] = sparse_graph_canon.e[i * REG + j];
            }

            /**
             * Remark: would be better just to determine max neighbour, instead of
             * first sorting the array (of 3 elements). Only in very few cases the
             * edge with the biggest label won't be reducible.
             * But left it like this since the transform_triangle_into_canonical_form_full
             * hardly consumes any cputime according to the profiler (0.15% of the time
             * for 22 3 and descending).
             */
            //Is needed, otherwise the "last" edge might not always be the same (the order of the edges doesnt matter in the e-list)
            transform_triangle_into_canonical_form_full(neighbours);

            for(j = REG - 1; j >= 0; j--) {
                neighbour = neighbours[j];

                if(!is_part_of_irreducible_triangle_bitvector(lab[neighbour])) {
                    /**
                     * lab is the order in which the vertices of the orginal graph should be relabelled
                     * in order to obtain the canonical graph.
                     * So lab[k] is the original label of vertex k
                     */

                    /* Transform into canonical form */
                    if(lab[i] < lab[neighbour]) {
                        lastedge[0] = lab[i];
                        lastedge[1] = lab[neighbour];
                    } else {
                        lastedge[0] = lab[neighbour];
                        lastedge[1] = lab[i];
                    }
                    if(colours_three[lastedge[0]][lastedge[1]] == min_colour_three) {
                        //fprintf(stderr, "Corresponding lastedge in org graph: %d %d\n", lastedge[0], lastedge[1]);
                        return;
                    }
                }
            }
        }
    }
    fprintf(stderr, "Error: last edge not found\n");
    exit(1);
}

/**
 * Returns 1 if the last inserted edge is the major edge, else returns 0.
 * (Using McKay's canonical construction path method).
 */
int is_major_edge(int partitions_already_set, sparsegraph *sparsegraph_canon) {
    int i;

    if(!partitions_already_set)
        determine_vertex_partitions();

    number_of_generators = 0;

    options.getcanon = TRUE;
    copy_sparse_graph();

    //nauty_sh((graph*) & sg, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, current_number_of_vertices, (graph*) & sg_canon);
    nauty_sh((graph*) & sg, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, current_number_of_vertices, (graph*) sparsegraph_canon);

    /* Check if the inserted edge was a major edge */

    DEBUGASSERT(number_of_reducible_triangles == 0);

    DEBUGASSERT(edgelist_size > 1 && edgelist_size <= 3 * current_number_of_vertices / 2);

    /* Last edge is edge = {current_number_of_vertices - 2, current_number_of_vertices - 1} */

    /* Determine the edge canonically labelled last by nauty  */
    //Actually we dont need to know the last edge, we are only interested in it's index.
    EDGE lastedge;
    determine_last_edge(*sparsegraph_canon, lab, lastedge);

    if(number_of_generators == 0) {
        return lastedge[0] == current_number_of_vertices - 2 && lastedge[1] == current_number_of_vertices - 1;
    } else {
        int edge_orbits[edgelist_size+1];
        int number_of_edge_orbits = 0;
        determine_edge_orbits(edgelist, edgelist_size, edge_orbits, &number_of_edge_orbits);

        if(number_of_edge_orbits == 1) {
            return 1;
        }

        i = edge_index[lastedge[0]][lastedge[1]];

        /**
         * Important: the inserted edge is always at the last position of the edgelist,
         * because in has_min_colour() it's the edge which was investigated last.
         */
        int new_edge = edgelist_size - 1;

        DEBUGASSERT(edge_index[current_number_of_vertices - 2][current_number_of_vertices - 1] == edgelist_size - 1);

        DEBUGASSERT(i < edgelist_size);

        if(edge_orbits[new_edge] != edge_orbits[i]) {
            return 0; //Last edge wasn't a major edge, so reject
        } else {
            return 1;
        }
    }
}

/**
 * Generates the edgepairs on level n-2 in case of girth = 4.
 *
 * The edgepairs must destroy all triangles.
 * Returns 1 if such edgepairs were found, else returns 0.
 */
int generate_edgepairs_penultimate_level_girth4(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size) {
    edgepair_list_index_squares = 0;
    if(number_of_reducible_triangles > 0 || number_of_irreducible_triangles > 0) {

        //Doesn't happen anymore due to bounding criterion
/*
        if((number_of_reducible_triangles == 1 && number_of_irreducible_triangles == 1) ||
                (number_of_reducible_triangles == 0 && number_of_irreducible_triangles == 2)) {
            return 0;
        }
*/

        //i.e. number of diamonds > 0. Cannot destroy a diamond without generating a square
        if(snarks && number_of_reducible_triangles == 0 && is_colourable()) {

            return 0;
        }
        *edge_pair_list_size = 0;

        DEBUGASSERT(number_of_reducible_triangles <= 2);

        //if(number_of_reducible_triangles >= 1 && snarks)
        //    DEBUGASSERT(!is_colourable());

        if(number_of_reducible_triangles == 1 && number_of_irreducible_triangles == 0) {
            generate_edgepairs_one_triangle(edge_pairs_list, edge_pair_list_size);

            //Can be < 15 in case there is a square
            DEBUGASSERT(*edge_pair_list_size <= 15);
        }
/*
        else if(number_of_reducible_triangles == 1 && number_of_irreducible_triangles == 1) {
            //Colour inserted edge is 11, but is never accepted because there will be an edge with colour 10
            return 0;
        }
*/
        else if(number_of_reducible_triangles == 0 && number_of_irreducible_triangles == 1) {
            //There is only one possible edgepair in this  case, otherwise the inserted edge won't have the min colour
            generate_edgepairs_triangle_free_one_diamond(edge_pairs_list, edge_pair_list_size);
        }
/*
        else if(number_of_reducible_triangles == 0 && number_of_irreducible_triangles == 2) {
            //Colour inserted edge is 10, but there will be an other edge which has colour 10 AND is part of a square
            //generate_edgepairs_triangle_free_two_diamonds(edge_pairs_list, edge_pair_list_size);

            return 0;
        }
*/
        else if(number_of_reducible_triangles == 2) {
            DEBUGASSERT(number_of_irreducible_triangles == 0);

            //Are always adjacent triangles if girth > 3
            //if(are_adjacent_triangles()) {
                //Edgepair must be part of square or pentagon, otherwise it wont have min colour:
                //there will be an edge with colour 10 which is part of 2 squares!
                generate_edgepairs_two_triangles(edge_pairs_list, edge_pair_list_size);
            //}

        } else {
            fprintf(stderr, "Error: too much triangles (should never happen because of bounding criterion)\n");
            exit(1);
        }
        
    } else {
        if(snarks && is_colourable()) {
            //If there are any squares left, all squares must be destroyed and the inserted edge may not be part of a square
            generate_edgepairs_penultimate_level_girth5_no_triangles(edge_pairs_list, edge_pair_list_size);
            edgepair_list_index_squares = *edge_pair_list_size;
        } else {
            generate_edgepairs_no_triangles(edge_pairs_list, edge_pair_list_size);
        }
    }

    return *edge_pair_list_size > 0;
}

/**
 * Generates the edgepairs on level n-2 in case of girth = 5.
 *
 * It is assumed that current_graph contains no triangles. The edgepairs must
 * destroy all squares.
 * Returns 1 if such edgepairs were found, else returns 0.
 */
int generate_edgepairs_penultimate_level_girth5_no_triangles(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size) {
    *edge_pair_list_size = 0;

    if(!find_squares(squares_global, squares_global_bitvectors, &squares_global_size, adjacent_squares_global, &adjacent_squares_global_size)) {
        return 0;
    }

    DEBUGASSERT(squares_global_size <= 4);
    DEBUGASSERT(current_number_of_vertices == number_of_vertices - 2);

    if(squares_global_size == 0) {
        int check_remove = 0;
        //This method is only called in case of girth 4 if the graph is colourable
        if(snarks && (girth == 4 || is_colourable())) {
            init_search_cycles();
            search_cycles();

            check_remove = 1;
            generate_non_adjacent_edge_pairs_square_free(edge_pairs_list, edge_pair_list_size, check_remove);

            if(*edge_pair_list_size > 0 && modify_existing_colouring()) {
                int i;
                int zero_indices[*edge_pair_list_size];
                int one_indices[*edge_pair_list_size];
                for(i = 0; i < *edge_pair_list_size; i++) {
                    zero_indices[i] = edge_labels[edge_pairs_list[i][0]][edge_pairs_list[i][1]];
                    one_indices[i] = edge_labels[edge_pairs_list[i][2]][edge_pairs_list[i][3]];
                }

                search_cycles();

                int current_size = 0;
                for(i = 0; i < *edge_pair_list_size; i++) {
                    if((edge_in_cycles[zero_indices[i]] & edge_in_cycles[one_indices[i]]) == 0) {
                        edge_pairs_list[current_size][0] = edge_pairs_list[i][0];
                        edge_pairs_list[current_size][1] = edge_pairs_list[i][1];
                        edge_pairs_list[current_size][2] = edge_pairs_list[i][2];
                        edge_pairs_list[current_size][3] = edge_pairs_list[i][3];
                        edgepair_index[zero_indices[i]][one_indices[i]] = current_size;

                        zero_indices[current_size] = zero_indices[i];
                        one_indices[current_size] = one_indices[i];

                        current_size++;
                    }
                }
                *edge_pair_list_size = current_size;
                
                //5 is more or less optimal
                if(*edge_pair_list_size > 5 && is_colourable_other_colouring()) {
                    search_cycles();
                    current_size = 0;
                    for(i = 0; i < *edge_pair_list_size; i++) {
                        if((edge_in_cycles[zero_indices[i]] & edge_in_cycles[one_indices[i]]) == 0) {
                            edge_pairs_list[current_size][0] = edge_pairs_list[i][0];
                            edge_pairs_list[current_size][1] = edge_pairs_list[i][1];
                            edge_pairs_list[current_size][2] = edge_pairs_list[i][2];
                            edge_pairs_list[current_size][3] = edge_pairs_list[i][3];
                            edgepair_index[zero_indices[i]][one_indices[i]] = current_size;

                            zero_indices[current_size] = zero_indices[i];
                            one_indices[current_size] = one_indices[i];

                            current_size++;
                        }
                    }
                    *edge_pair_list_size = current_size;
                }

            }
        } else
            generate_non_adjacent_edge_pairs_square_free(edge_pairs_list, edge_pair_list_size, check_remove);
    } else if(squares_global_size == 1) {
        DEBUGASSERT(adjacent_squares_global_size == 0);
        generate_non_adjacent_edge_pairs_one_square(edge_pairs_list, edge_pair_list_size, squares_global[0]);
    } else if(squares_global_size == 2) {
        DEBUGASSERT(adjacent_squares_global_size <= 1);
        if(adjacent_squares_global_size == 0) {
            if(squares_have_multiple_common_neighbours(squares_global)) {
                generate_non_adjacent_edge_pairs_two_squares(edge_pairs_list, edge_pair_list_size, squares_global[0], squares_global[1]);
            } else {
                *edge_pair_list_size = 0;
            }
        } else { //adjacent_squares_sizes == 1
            EDGE common_edge;
            determine_common_edge(squares_global_bitvectors[adjacent_squares_global[0][0]] & squares_global_bitvectors[adjacent_squares_global[0][1]], common_edge);
            generate_non_adjacent_edge_pairs_two_adjacent_squares(edge_pairs_list, edge_pair_list_size, squares_global, common_edge);
        }
    } else if(squares_global_size == 4) {
        DEBUGASSERT(adjacent_squares_global_size == 2);
        //Else can never yield an edge with minimal colour
        if(current_number_of_vertices == 8)
            generate_non_adjacent_edge_pairs_four_squares(edge_pairs_list, edge_pair_list_size, squares_global_bitvectors, adjacent_squares_global, adjacent_squares_global_size);
        else
            *edge_pair_list_size = 0;
    } else {
        //squares_size == 3
        fprintf(stderr, "Error: too much squares: %d\n", squares_global_size);
        exit(1);
    }

    /**
     * In case of squares_global_size > 0 && snarks it's faster to check first if there are
     * any possible edgepairs. And only if that's the case removing the edgepairs
     * which are on the same colour cycle.
     */
    /**
     * In case of girth 5, it's best to call is_colourable() only if edge_pairs_list
     * contains quite some elements.
     */
    if(squares_global_size > 0 && *edge_pair_list_size > min_edgepairlist_size && snarks && (girth == 4 || is_colourable())) {
        int i;
        int zero_indices[*edge_pair_list_size];
        int one_indices[*edge_pair_list_size];
        for(i = 0; i < *edge_pair_list_size; i++) {
            zero_indices[i] = edge_labels[edge_pairs_list[i][0]][edge_pairs_list[i][1]];
            one_indices[i] = edge_labels[edge_pairs_list[i][2]][edge_pairs_list[i][3]];
        }

        init_search_cycles();
        search_cycles_square(squares_global, squares_global_size);

        int current_size = 0;
        for(i = 0; i < *edge_pair_list_size; i++) {
            if((edge_in_cycles[zero_indices[i]] & edge_in_cycles[one_indices[i]]) == 0) {
                edge_pairs_list[current_size][0] = edge_pairs_list[i][0];
                edge_pairs_list[current_size][1] = edge_pairs_list[i][1];
                edge_pairs_list[current_size][2] = edge_pairs_list[i][2];
                edge_pairs_list[current_size][3] = edge_pairs_list[i][3];
                edgepair_index[zero_indices[i]][one_indices[i]] = current_size;

                current_size++;
            }
        }
        *edge_pair_list_size = current_size;

        //Remark: also trying a modified colouring is slightly slower
    }

    return *edge_pair_list_size > 0;
}

/**
 * Generates the edgepairs on any level where the graph after extension
 * should have girth 5. So all squares must be destroyed.
 * It is assumed that current_graph contains no reducible triangles, but
 * it MAY contain diamonds.
 *
 * Returns 1 if such edgepairs were found, else returns 0.
 */
int generate_edgepairs_penultimate_level_girth5_no_triangles_but_diamonds(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size) {
    SQUARE squares[4];
    int squares_size = 0;

    *edge_pair_list_size = 0;

    setword squares_bitvectors[4];

    EDGE adjacent_squares[2];
    int adjacent_squares_size = 0;

    if(!find_squares_contains_diamonds(squares, squares_bitvectors, &squares_size, adjacent_squares, &adjacent_squares_size)) {
        return 0;
    }

    DEBUGASSERT(squares_size <= 4);

    if(squares_size == 0) {
        generate_non_adjacent_edge_pairs_square_free_contains_diamonds(edge_pairs_list, edge_pair_list_size);
    } else if(squares_size == 1) {
        DEBUGASSERT(adjacent_squares_size == 0);
        generate_non_adjacent_edge_pairs_one_square_contains_diamonds(edge_pairs_list, edge_pair_list_size, squares[0]);
    } else if(squares_size == 2) {
        DEBUGASSERT(adjacent_squares_size <= 1);
        if(adjacent_squares_size == 0) {
            if(squares_have_multiple_common_neighbours(squares)) {
                generate_non_adjacent_edge_pairs_two_squares(edge_pairs_list, edge_pair_list_size, squares[0], squares[1]);
            } else {
                *edge_pair_list_size = 0;
            }
        } else { //adjacent_squares_sizes == 1
            EDGE common_edge;
            determine_common_edge(squares_bitvectors[adjacent_squares[0][0]] & squares_bitvectors[adjacent_squares[0][1]], common_edge);
            generate_non_adjacent_edge_pairs_two_adjacent_squares(edge_pairs_list, edge_pair_list_size, squares, common_edge);
        }
    } else if(squares_size == 4) {
        DEBUGASSERT(adjacent_squares_size == 2);
        //Can never yield an edge with minimal colour
        if(current_number_of_vertices == 8)
            generate_non_adjacent_edge_pairs_four_squares(edge_pairs_list, edge_pair_list_size, squares_bitvectors, adjacent_squares, adjacent_squares_size);
        else
            *edge_pair_list_size = 0;
    } else {
        //squares_size == 3
        fprintf(stderr, "Error: too much squares: %d (case diamonds)\n", squares_global_size);
        exit(1);
    }

    return *edge_pair_list_size > 0;
}

/**
 * Returns 1 if the edge was already marked by a certain colour,
 * else returns 0.
 */
int ismarked_colour(int from, int index_other, int colour) {
    if(colour == 1) {
        return ISMARKED_CYCLE_COLOUR1(from, index_other);
    } else if(colour == 2) {
        return ISMARKED_CYCLE_COLOUR2(from, index_other);
    } else if(colour == 3) {
        return ISMARKED_CYCLE_COLOUR3(from, index_other);
    } else {
        fprintf(stderr, "Error: invalid colour: %d\n", colour);
        exit(1);
    }

}

/**
 * Marks the colour of an edge.
 * Remark: it doesnt matter if from is > or < to.
 */
void mark_colour(int from, int to, int colour) {
    int to_index = neighbour_index[from][to];
    int from_index = neighbour_index[to][from];
    if(colour == 1) {
        MARK_CYCLE_COLOUR1(from, to_index);
        MARK_CYCLE_COLOUR1(to, from_index);
    } else if(colour == 2) {
        MARK_CYCLE_COLOUR2(from, to_index);
        MARK_CYCLE_COLOUR2(to, from_index);
    } else if(colour == 3) {
        MARK_CYCLE_COLOUR3(from, to_index);
        MARK_CYCLE_COLOUR3(to, from_index);
    } else {
        fprintf(stderr, "Error: invalid colour: %d\n", colour);
        exit(1);
    }
}

/**
 * Determines an even colour cycle of colours[0] and colours[1] which contains
 * start_edge.
 * Important: it is assumed that start_edge has colour colours[0].
 */
void determine_even_cycle(EDGE start_edge, EDGE colours, unsigned char cycle[], int *cycle_size) {
    DEBUGASSERT(colours_snarks[start_edge[0]][neighbour_index[start_edge[0]][start_edge[1]]] == colours[0]);

    cycle[0] = start_edge[0];
    //cycle[1] = start_edge[1];
    *cycle_size = 1;
    int current_vertex = start_edge[1];
    int current_colour_index = 1;
    while(current_vertex != start_edge[0]) {
        cycle[*cycle_size] = current_vertex;
        (*cycle_size)++;
        current_vertex = find_next_vertex(current_vertex, colours[current_colour_index]);
        current_colour_index = (current_colour_index + 1) % 2;
    }

    //Cycle will always be even because of alternating colours
    DEBUGASSERT(*cycle_size % 2 == 0);
}

void init_search_cycles() {
    int i;
    for(i = 0; i < 3 * current_number_of_vertices / 2; i++)
        edge_in_cycles[i] = (setword) 0;
    current_cycle_number = 0;
}

/**
 * Forms all colour cycles of current_graph.
 * For each edge a bitvector is kept with the cycles where it is part of.
 * Pairs of edges which are part of the same colour cycle do not have to be extended,
 * because the extended graph will also be colourable.
 *
 * Important: current_graph must already be coloured before calling this method
 */
void search_cycles() {
    unsigned char cycle[current_number_of_vertices];
    int neighbour, colour, cycle_size;
    EDGE available_colours, start_edge, temp_edge, colours_even_cycle;
    int i, j, k;

    RESETMARKS_CYCLE_COLOUR1;
    RESETMARKS_CYCLE_COLOUR2;
    RESETMARKS_CYCLE_COLOUR3;

    for(i = 0; i < current_number_of_vertices - 1; i++) {
        for(j = 0; j < REG; j++) {
            neighbour = current_graph[i][j];
            if(i < neighbour) {
                start_edge[0] = i;
                start_edge[1] = neighbour;
                colour = colours_snarks[i][j];
                colours_even_cycle[0] = colour;
                determine_available_colours(colour, available_colours);
                for(k = 0; k < 2; k++) {
                    if(!ismarked_colour(i, j, available_colours[k])) {
                        cycle_size = 0;
                        colours_even_cycle[1] = available_colours[k];
                        determine_even_cycle(start_edge, colours_even_cycle, cycle, &cycle_size);
                        if(current_cycle_number >= MAX_NUMBER_OF_CYCLES) {
                            fprintf(stderr, "Error: too much colour cycles found\n");
                            exit(1);
                        }

                        int l;
                        for(l = 0; l < cycle_size; l++) {
                            temp_edge[0] = cycle[l];
                            temp_edge[1] = cycle[(l + 1) % cycle_size];
                            //transform_edge_into_canonical_form(temp_edge);
                            mark_colour(temp_edge[0], temp_edge[1], colours_even_cycle[(l + 1) % 2]);

                            edge_in_cycles[edge_labels[temp_edge[0]][temp_edge[1]]] |= BIT(current_cycle_number);
                        }
                        current_cycle_number++;
                    }
                }
            }
        }
    }
}


/**
 * Forms all colour cycles which contain an edge which is part of a square.
 * For each edge a bitvector is kept with the cycles where it is part of.
 * Pairs of edges which are part of the same colour cycle do not have to be extended,
 * because the extended graph will also be colourable.
 *
 * Important: current_graph must already be coloured before calling this method
 */
void search_cycles_square(SQUARE squares[], int number_of_squares) {
    unsigned char cycle[current_number_of_vertices];
    int colour, cycle_size;
    EDGE available_colours, start_edge, temp_edge, colours_even_cycle;
    int i, j, k, l;

    RESETMARKS_CYCLE_COLOUR1;
    RESETMARKS_CYCLE_COLOUR2;
    RESETMARKS_CYCLE_COLOUR3;

    for(l = 0; l < number_of_squares; l++)
        for(i = 0; i < 4; i++) {
            start_edge[0] = squares[l][i];
            start_edge[1] = squares[l][(i + 1) % 4];

            colour = colours_snarks[start_edge[0]][neighbour_index[start_edge[0]][start_edge[1]]];
            colours_even_cycle[0] = colour;
            determine_available_colours(colour, available_colours);
            for(j = 0; j < 2; j++) {
                if(!ismarked_colour(start_edge[0], neighbour_index[start_edge[0]][start_edge[1]], available_colours[j])) {
                    cycle_size = 0;
                    colours_even_cycle[1] = available_colours[j];
                    determine_even_cycle(start_edge, colours_even_cycle, cycle, &cycle_size);

                    if(current_cycle_number >= MAX_NUMBER_OF_CYCLES) {
                        fprintf(stderr, "Error: too much colour cycles found\n");
                        exit(1);
                    }

                    for(k = 0; k < cycle_size; k++) {
                        temp_edge[0] = cycle[k];
                        temp_edge[1] = cycle[(k + 1) % cycle_size];
                        //transform_edge_into_canonical_form(temp_edge);
                        mark_colour(temp_edge[0], temp_edge[1], colours_even_cycle[(k + 1) % 2]);

                        edge_in_cycles[edge_labels[temp_edge[0]][temp_edge[1]]] |= BIT(current_cycle_number);
                    }
                    current_cycle_number++;
                } 
            }
        }
    
}

/**
 * Extends a list of edgepairs.
 */
void edge_extend(EDGEPAIR edge_pairs_list[], int edge_pair_list_size) {
    DEBUGASSERT(edge_pair_list_size > 0);

    int edgepair_orbits[edge_pair_list_size+1];
    int number_of_edgepair_orbits = 0;

    int is_trivial_group = 0;
    if(number_of_generators > 0) {
        //The generators are still valid, because edge_extension is done before triangle extension
        determine_edgepair_orbits(edge_pairs_list, edge_pair_list_size, edgepair_orbits, &number_of_edgepair_orbits);
    } else {
        is_trivial_group = 1;
        number_of_edgepair_orbits = edge_pair_list_size;
    }

    //Backup triangle lists
    //Somehow is slightly faster with +1
    int number_of_reducible_triangles_local = number_of_reducible_triangles;
    TRIANGLE reducible_triangles_local[number_of_reducible_triangles+1];
    if(number_of_reducible_triangles_local > 0)
        memcpy(reducible_triangles_local, reducible_triangles, sizeof (TRIANGLE) * number_of_reducible_triangles_local);

    int number_of_irreducible_triangles_local = number_of_irreducible_triangles;
    setword irreducible_triangles_bitvector_local;
    IRRED_TRIANGLE irreducible_triangles_local[number_of_irreducible_triangles+1];
    if(number_of_irreducible_triangles_local > 0) {
        memcpy(irreducible_triangles_local, irreducible_triangles, sizeof (IRRED_TRIANGLE) * number_of_irreducible_triangles_local);
        irreducible_triangles_bitvector_local = irreducible_triangles_bitvector;
    }

    int number_of_bridges_local = number_of_bridges;
    EDGE bridges_local[number_of_bridges+1];
    if(number_of_bridges_local > 0)
        memcpy(bridges_local, bridges, sizeof(EDGE) * number_of_bridges_local);

    EDGE min_edge_local;
    min_edge_local[0] = min_edge[0];
    min_edge_local[1] = min_edge[1];

    int edgepair_list_index_squares_local = edgepair_list_index_squares;

    int i, j;
    int num_orbits = 0;
    int edge_inserted;

    EDGE inserted_edge;

    EDGE previous_min_edge_local;

    previous_min_edge_local[0] = 0;
    previous_min_edge_local[1] = 0;

    for(i = 0; i < edge_pair_list_size; i++) {
        if(is_trivial_group || edgepair_orbits[i] == i) {
            add_edge(edge_pairs_list[i]);

            //Bounding criterion
            int abort = 0;
            if(girth > 3 && current_number_of_vertices < number_of_vertices) {
                int vertices_required = ((number_of_irreducible_triangles + number_of_reducible_triangles + 1) / 2) * 2;
                if(girth > 4)
                    vertices_required *= 2;
                abort = current_number_of_vertices + vertices_required > number_of_vertices;
            }

            /**
             * This only really helps in case of snarks, otherwise there are just too
             * many graphs which have connectivity > 2.
             */
            if((check_cyclic_connectivity && min_cyclic_connectivity > 1) && !abort)  {
                int vertices_required = (min_cyclic_connectivity - 1) * 2 * MIN(number_of_bridges, 1);
                if(current_number_of_vertices + vertices_required > number_of_vertices) {
                    abort = 1;
                }
            }

            if(!abort) {
                all_edges_are_reducible = number_of_bridges == 0 && number_of_irreducible_triangles == 0;

                inserted_edge[0] = current_number_of_vertices - 2;
                inserted_edge[1] = current_number_of_vertices - 1;

                previous_min_edge[0] = previous_min_edge_local[0];
                previous_min_edge[1] = previous_min_edge_local[1];

                if(i >= edgepair_list_index_squares) {
                    //i.e. the inserted edge will be part of a square
                    if(has_min_colour_cycle(inserted_edge, &edge_inserted)) {
                        extend(edge_inserted, 0);

                        min_edge[0] = min_edge_local[0];
                        min_edge[1] = min_edge_local[1];
                        edgepair_list_index_squares = edgepair_list_index_squares_local;
                    } else {
                        previous_min_edge_local[0] = previous_min_edge[0];
                        previous_min_edge_local[1] = previous_min_edge[1];
                    }
                } else {
                    //i.e. the extended graph doesn't contain any squares

                    EDGE neighbours[4];
                    for(j = 0; j < 4; j++) {
                        neighbours[j][0] = edge_pairs_list[i][j];
                        neighbours[j][1] = inserted_edge[j >= 2];
                    }

                    //Extend is only called if the inserted edge has minimal colour
                    if(has_min_colour(inserted_edge, &edge_inserted, neighbours)) {
                        extend(edge_inserted, 0);

                        min_edge[0] = min_edge_local[0];
                        min_edge[1] = min_edge_local[1];
                        edgepair_list_index_squares = edgepair_list_index_squares_local;
                    } else {
                        previous_min_edge_local[0] = previous_min_edge[0];
                        previous_min_edge_local[1] = previous_min_edge[1];
                    }
                }
            }

            remove_edge(edge_pairs_list[i]);

            //Restore triangle lists
            //Red triangles only have to be restored after the last recursive call (because number of red triangles will always be 0 after add_edge)
            if(num_orbits == number_of_edgepair_orbits - 1) {
                number_of_reducible_triangles = number_of_reducible_triangles_local;
                if(number_of_reducible_triangles_local > 0)
                    memcpy(reducible_triangles, reducible_triangles_local, sizeof (TRIANGLE) * number_of_reducible_triangles_local);
            }

            number_of_irreducible_triangles = number_of_irreducible_triangles_local;
            if(number_of_irreducible_triangles_local > 0) {
                memcpy(irreducible_triangles, irreducible_triangles_local, sizeof (IRRED_TRIANGLE) * number_of_irreducible_triangles_local);
                irreducible_triangles_bitvector = irreducible_triangles_bitvector_local;
            }

            //Restore is_bridge
            for(j = 0; j < number_of_bridges; j++)
                is_bridge[bridges[j][0]][bridges[j][1]] = 0;

            number_of_bridges = number_of_bridges_local;
            if(number_of_bridges_local > 0)
                memcpy(bridges, bridges_local, sizeof(EDGE) * number_of_bridges_local);

            //Restore is_bridge
            for(j = 0; j < number_of_bridges; j++)
                is_bridge[bridges[j][0]][bridges[j][1]] = 1;

            num_orbits++;
            if(num_orbits == number_of_edgepair_orbits)
                break;
        }
    }
    DEBUGASSERT(num_orbits == number_of_edgepair_orbits);
}

/**
 * Returns 1 if vertexset[vertexset_index] is the only extremal vertex of a diamond of
 * vertexset, else returns 0.
 * If 1 is returned, irred_triangle points to the irreducible triangle which contains
 * vertexset[vertexset_index] as extremal vertex.
 */
int is_only_vertex_from_irred_triangle(unsigned char vertexset[], int num_vertices_in_set, int vertexset_index, int *irred_triangle) {
    int i;
    for(i = 0; i < number_of_irreducible_triangles; i++) {
        if(vertexset[vertexset_index] == irreducible_triangles[i][0]) {
            if(is_only_vertex_from_diamond(vertexset, num_vertices_in_set, i, 0)) {
                *irred_triangle = i;
                return 1;
            } else
                return 0;
        } else if(vertexset[vertexset_index] == irreducible_triangles[i][3]) {
            if(is_only_vertex_from_diamond(vertexset, num_vertices_in_set, i, 3)) {
                *irred_triangle = i;
                return 1;
            } else
                return 0;
        }
    }
    return 0;
}

/**
 * Determines partitions for nauty for a current graph which contains triangles.
 *
 * Remark: order of partitions in lab doesn't matter because triangle_extend is always canonical
 */
void determine_nauty_partitions_triangle(unsigned char vertexset[], int vertexset_size, int vertex_orbits_local[], int number_of_vertex_orbits) {
    int i, j;
    unsigned char already_visited[vertexset_size+1];
    for(i = 0; i < vertexset_size; i++)
        already_visited[i] = 0;

    //Could use marks for this
    unsigned char is_part_of_vertexset[current_number_of_vertices];
    for(i = 0; i < current_number_of_vertices; i++) {
        is_part_of_vertexset[i] = 0;
    }
    for(i = 0; i < vertexset_size; i++) {
        is_part_of_vertexset[vertexset[i]] = 1;
    }

    int lab_index = 0;
    if(number_of_irreducible_triangles > 0) {
        //Put all vertices of the diamonds of which exactly one extremal vertex is expanded in the same partition
        //The partition could be finer, but this case only rarely occurs and this is a cheap way to handle this rather rare case
        int irred_triangle;
        for(i = 0; i < vertexset_size; i++) {
            if(is_only_vertex_from_irred_triangle(vertexset, vertexset_size, i, &irred_triangle)) {
                already_visited[i] = 1;
                for(j = 0; j < 4; j++) {
                    lab[lab_index] = irreducible_triangles[irred_triangle][j];
                    ptn[lab_index] = 1;
                    lab_index++;

                    //To make sure this vertex won't be put twice in the list
                    is_part_of_vertexset[irreducible_triangles[irred_triangle][j]] = 1;
                }
                lab[lab_index] = current_number_of_vertices + 2 * i;
                ptn[lab_index] = 1;
                lab_index++;

                lab[lab_index] = current_number_of_vertices + 2 * i + 1;
                ptn[lab_index] = 1;
                lab_index++;
            }
        }
        if(lab_index > 0)
            ptn[lab_index - 1] = 0;
    }

    /**
     * Put the remaining vertices from the vertexset (and their generated vertices)
     * in the right partition.
     * Vertices of the vertexset which are in different orbits before expansion,
     * will also be in different orbits after expansion.
     */
    for(i = 0; i < vertexset_size; i++) {
        if(!already_visited[i]) {
            lab[lab_index] = vertexset[i];
            ptn[lab_index] = 1;
            lab_index++;
            for(j = 0; j < 2; j++) {
                lab[lab_index] = current_number_of_vertices + 2 * i + j;
                ptn[lab_index] = 1;
                lab_index++;
            }
            for(j = i + 1; j < vertexset_size; j++) {
                //already_visited is needed when vertexset size > 2
                if(!already_visited[j] && vertex_orbits_local[vertexset[j]] == vertex_orbits_local[vertexset[i]]) {
                    already_visited[j] = 1;
                    lab[lab_index] = vertexset[j];
                    ptn[lab_index] = 1;
                    lab_index++;
                    int k;
                    for(k = 0; k < 2; k++) {
                        lab[lab_index] = current_number_of_vertices + 2 * j + k;
                        ptn[lab_index] = 1;
                        lab_index++;
                    }
                }
            }
            ptn[lab_index - 1] = 0;
        }
    }

    /* Put vertices which aren't part of the vertexset in their partitions (according to their old orbits) */
    //Could do this more efficiently, but this method hardly requires any cpu according to the profiler
    for(i = 0; i < current_number_of_vertices; i++) {
        for(j = 0; j < current_number_of_vertices; j++) {
            if(!is_part_of_vertexset[j] && vertex_orbits_local[j] == i) {
                lab[lab_index] = j;
                ptn[lab_index] = 1;
                lab_index++;
            }
        }
        DEBUGASSERT(lab_index > 0);
        ptn[lab_index - 1] = 0;
    }

    DEBUGASSERT(lab_index == current_number_of_vertices + 2 * vertexset_size);

    options.defaultptn = FALSE;
}

/**
 * Returns 1 if at least two elements of the vertexset are in the same orbit.
 */
int contains_vertices_which_are_in_same_orbit(unsigned char vertexset[], int vertexset_size, int vertex_orbits[]) {
    if(vertexset_size == 1)
        return 0;
    else {
        RESETMARKS;
        int i;
        for(i = 0; i < vertexset_size; i++) {
            if(ISMARKED(vertex_orbits[vertexset[i]]))
                return 1;
            else
                MARK(vertex_orbits[vertexset[i]]);
        }
        return 0;
    }
}

/**
 * Determines possible vertexsets and transforms the vertices of those vertexsets
 * into triangles. All reducible triangles of the old graph have to be destroyed
 * (so each reducible triangle must contain at least 1 vertex which will be blown up).
 * In this case the resulting graph is always canonical (so nauty doesn't have to be called).
 */
void triangle_extend_all(permutation generators_local[][MAXN], int number_of_generators_local,
        int vertex_orbits_local[], int number_of_vertex_orbits, int groupsize) {
    int max_num_of_vertices = number_of_vertices;
    if(girth > 3) {
        if(number_of_reducible_triangles > 1)
            return;
        max_num_of_vertices -= 2;
    }
    //Splitlevel is <= number_of_vertices - 2
    if(modulo && splitlevel_counter != rest && current_number_of_vertices < splitlevel) {
        //if(modulo && rest != 0 && current_number_of_vertices < splitlevel) {
        max_num_of_vertices = splitlevel;
    }
    int max_vertexset_size = (max_num_of_vertices - current_number_of_vertices) / 2;
    max_vertexset_size = MIN(max_vertexset_size, current_number_of_vertices);

    if(girth > 3)
        max_vertexset_size = MIN(max_vertexset_size, 2); //A new edge can at most reduce 2 triangles

    int min_i = MAX(number_of_reducible_triangles, 1);
    if(min_i > max_vertexset_size) {

        return;
    }

    if (girth > 3) {
        int vertices_required = ((number_of_irreducible_triangles + number_of_reducible_triangles + 1) / 2) * 2;
        if (girth > 4)
            vertices_required *= 2;
        if (current_number_of_vertices + vertices_required + 2 * min_i > number_of_vertices) {
            return;
        }
        vertices_required = ((min_i + 1) / 2) * 2;
        if (girth > 4)
            vertices_required *= 2;
        if (current_number_of_vertices + vertices_required + 2 * min_i > number_of_vertices) {
            return;
        }
    }

    //Backup triangle lists
    int number_of_reducible_triangles_local;
    TRIANGLE reducible_triangles_local[number_of_reducible_triangles+1];

    int number_of_irreducible_triangles_local;
    setword irreducible_triangles_bitvector_local;
    IRRED_TRIANGLE irreducible_triangles_local[number_of_irreducible_triangles+1];


    int number_of_bridges_local;
    EDGE bridges_local[number_of_bridges+1];

    if(current_number_of_vertices + 2 * min_i < number_of_vertices) {
        //Backup not needed on last level

        number_of_reducible_triangles_local = number_of_reducible_triangles;
        if(number_of_reducible_triangles_local > 0)
            memcpy(reducible_triangles_local, reducible_triangles, sizeof (TRIANGLE) * number_of_reducible_triangles_local);

        number_of_irreducible_triangles_local = number_of_irreducible_triangles;
        if(number_of_irreducible_triangles_local > 0) {
            memcpy(irreducible_triangles_local, irreducible_triangles, sizeof (IRRED_TRIANGLE) * number_of_irreducible_triangles_local);
            irreducible_triangles_bitvector_local = irreducible_triangles_bitvector;
        }

        number_of_bridges_local = number_of_bridges;
        if(number_of_bridges_local > 0)
            memcpy(bridges_local, bridges, sizeof (EDGE) * number_of_bridges_local);
    }

    int num_orbits;
    int number_of_generators_new;

    int i, j;
    unsigned char orbit_frequencies[current_number_of_vertices];
    //Determine the frequency of each vertex orbit
    if(number_of_generators_local > 0) {
        num_orbits = 0;
        int k;
        for(k = 0; k < current_number_of_vertices; k++) {
            if(vertex_orbits_local[k] == k) {
                int frequency = 1;
                for(j = k + 1; j < current_number_of_vertices; j++) {
                    if(vertex_orbits_local[j] == k)
                        frequency++;
                }
                orbit_frequencies[k] = frequency;
                num_orbits++;
                if(num_orbits == number_of_vertex_orbits)
                    break;
            }
        }
    }

    for(i = min_i; i <= max_vertexset_size; i++) {
        if(girth > 3 && i > min_i) {
            int vertices_required = ((number_of_irreducible_triangles + number_of_reducible_triangles + 1) / 2) * 2;
            if(girth > 4)
                vertices_required *= 2;
            if(current_number_of_vertices + vertices_required + 2 * i > number_of_vertices) {
                return;
            }
            vertices_required = ((i + 1) / 2) * 2;
            if(girth > 4)
                vertices_required *= 2;
            if(current_number_of_vertices + vertices_required + 2 * i > number_of_vertices) {
                return;
            }
        }

        if(check_cyclic_connectivity && min_cyclic_connectivity > 1) {
            int vertices_required = (min_cyclic_connectivity - 1) * 2 * MIN(number_of_bridges, 1);
            //Triangle extend never destroys any bridges
            if(current_number_of_vertices + 2 * i + vertices_required > number_of_vertices) {
                return;
            }
        }

        //An upperbound for the max number of vertexsets
        int max_number_of_vertexsets = binom_coefficients[current_number_of_vertices][i];
        //+1 needed because of memcpy in generate_all_vertexsets
        unsigned char vertexset[max_number_of_vertexsets + 1][i];
        int num_of_vertexsets;

        if(girth == 3 || i == 1) {
            if(snarks && current_number_of_vertices == number_of_vertices - 4 && is_colourable()) {
                return; //A single triangle can only be destroyed if the inserted edge will be part of a square, but in that case the graph is still colourable
            } else {
                //Special cases hardly yield any speedup
                find_vertexsets(i, vertexset, &num_of_vertexsets);
            }
        } else { //girth > 3 && i == 2
            DEBUGASSERT(girth > 3 && i == 2);
            if(snarks && current_number_of_vertices == number_of_vertices - 6 && is_colourable()) {
                //If a graph is colourable, it will remain colourable after expanding.
                //The triangles must be destroyed, so there will be edges which are part of a square
                //This means that the inserted edge also should be part of a square (otherwise it won't be minimal), but then the graph will still be colourable

                return;
            } else {
                find_vertexsets_girth4(i, vertexset, &num_of_vertexsets);
            }
        }

        DEBUGASSERT(num_of_vertexsets <= max_number_of_vertexsets);

        int vertexset_orbits[num_of_vertexsets+1];
        int number_of_vertexset_orbits = 0;

        if(number_of_generators_local > 0) {
            //Could reuse the orbits if i == 1 && num_red_tria == 0 && num_irred_tria == 0
            //But won't yield a decent speedup, since the percentage of nontrivial groups goes to 0% as the order of the graphs increases
            determine_vertexset_orbits(generators_local, number_of_generators_local, i, vertexset, num_of_vertexsets, vertexset_orbits, &number_of_vertexset_orbits);
        } else {
            number_of_vertexset_orbits = num_of_vertexsets;
        }

        num_orbits = 0;
        for(j = 0; j < num_of_vertexsets; j++) {
            if(number_of_generators_local == 0 || vertexset_orbits[j] == j) {
                number_of_generators_new = number_of_generators_local;
                if(number_of_generators_local > 0 && current_number_of_vertices + i * 2 < number_of_vertices) {
                    /**
                     * If the frequency of the orbit of a vertex v which will be
                     * blown up is equal to the groupsize, this means there is
                     * no automorphism g \in Aut(G) for which g(v) = v.
                     * So after the triangle extension the group will be trivial.
                     *
                     * Of course this doesn't hold if the vertexset contains
                     * vertices which are in the same orbit.
                     */
                    /**
                     * This optimization doesn't yield much speedup.
                     * Could also first check if orbitsize == grpsize and afterwards
                     * also check if there isn't another vertex in the vertexset
                     * which is in the same orbit. But this won't be much faster.
                     */
                    if(!contains_vertices_which_are_in_same_orbit(vertexset[j], i, vertex_orbits_local)) {
                        int k, orbit;
                        for(k = 0; k < i; k++) {
                            orbit = vertex_orbits_local[vertexset[j][k]];

                            if(orbit_frequencies[orbit] == groupsize) {
                                number_of_generators_new = 0;
                                break;
                            }
                        }
                    }

                    if(number_of_generators_new > 0) {
                        determine_nauty_partitions_triangle(vertexset[j], i, vertex_orbits_local, number_of_vertex_orbits);
                    }
                }

                transform_vertexset_into_triangles(vertexset[j], i);

                //To make sure no lookahead is applied
                min_colour_one = MAX_EDGE_COLOUR_TWO;
                min_edge_is_part_of_square = -1;
                min_edge[0] = 0;
                min_edge[1] = 0;

                extend(TRIANGLE_INSERTED, number_of_generators_new == 0);

                undo_vertexset_triangle_transformation(vertexset[j], i);

                if(current_number_of_vertices + i * 2 < number_of_vertices) {
                    //Restore triangle lists
                    number_of_reducible_triangles = number_of_reducible_triangles_local;
                    if(number_of_reducible_triangles_local > 0)
                        memcpy(reducible_triangles, reducible_triangles_local, sizeof(TRIANGLE) * number_of_reducible_triangles_local);

                    number_of_irreducible_triangles = number_of_irreducible_triangles_local;
                    if(number_of_irreducible_triangles_local > 0) {
                        memcpy(irreducible_triangles, irreducible_triangles_local, sizeof(IRRED_TRIANGLE) * number_of_irreducible_triangles_local);
                        irreducible_triangles_bitvector = irreducible_triangles_bitvector_local;
                    }

                    //Restore is_bridge
                    int k;
                    for(k = 0; k < number_of_bridges; k++)
                        is_bridge[bridges[k][0]][bridges[k][1]] = 0;

                    number_of_bridges = number_of_bridges_local;
                    if(number_of_bridges_local > 0)
                        memcpy(bridges, bridges_local, sizeof(EDGE) * number_of_bridges_local);

                    //Restore is_bridge
                    for(k = 0; k < number_of_bridges; k++)
                        is_bridge[bridges[k][0]][bridges[k][1]] = 1;
                }
                num_orbits++;
                if(num_orbits == number_of_vertexset_orbits)
                    break;
            }
        }
        DEBUGASSERT(num_orbits == number_of_vertexset_orbits);

    }
}

/**
 * Generate all sets which consist of num_vertices_in_set vertices which eliminate all reducible triangles.
 */
void find_vertexsets(int num_vertices_in_set, unsigned char vertexset[][num_vertices_in_set], int *vertexset_size) {
    *vertexset_size = 0;
    DEBUGASSERT(number_of_reducible_triangles <= num_vertices_in_set);
    if(number_of_reducible_triangles > 0)
        generate_triangle_vertexsets(0, num_vertices_in_set, vertexset, vertexset_size);
    else
        generate_all_vertexsets(0, num_vertices_in_set, vertexset, vertexset_size);
}

/**
 * Generate all vertexsets which will destroy all reducible triangles.
 * It is assumed that number_of_reducible_triangles <= num_vertices_in_set.
 */
void generate_triangle_vertexsets(int current_index, int num_vertices_in_set, unsigned char vertexset[][num_vertices_in_set], int *vertexset_size) {
    if(current_index < number_of_reducible_triangles) { //First destroy the reducible triangles
        int i;
        for(i = 0; i < 3; i++) {
            vertexset[*vertexset_size][current_index] = reducible_triangles[current_index][i];
            generate_triangle_vertexsets(current_index + 1, num_vertices_in_set, vertexset, vertexset_size);
        }
    } else {
        DEBUGASSERT(current_index == number_of_reducible_triangles);

        //Once the reducible triangles are destroyed, the remaining vertices can be chosen freely
        generate_triangle_vertexsets_remaining(current_index, num_vertices_in_set, vertexset, vertexset_size);
    }
}

/**
 * Adds all possible remaining vertices to the vertexsets till it has the right size.
 * The remaining vertices may be part of a reducible triangle, but can't be identical to
 * a vertex which was already present in the vertexset.
 *
 * Remark: in the first call current_index should be number_of_reducible_triangles.
 */
void generate_triangle_vertexsets_remaining(int current_index, int num_vertices_in_set, unsigned char vertexset[][num_vertices_in_set], int *vertexset_size) {
    int i;
    if(current_index < num_vertices_in_set) {
        int previous = -1;
        if(current_index > number_of_reducible_triangles)
            previous = vertexset[*vertexset_size][current_index - 1];
        for(i = previous + 1; i < current_number_of_vertices; i++) {
            int triangle;
            //Valid if it's not part of a reducible triangle or it is and it can be added (its in canonical form)
            if(!is_part_of_reducible_triangle(i, &triangle) || can_be_added_to_reducible_triangle_vertexset(vertexset[*vertexset_size], current_index, i, triangle)) {
                vertexset[*vertexset_size][current_index] = i;
                generate_triangle_vertexsets_remaining(current_index + 1, num_vertices_in_set, vertexset, vertexset_size);
            }
        }
    } else {
        (*vertexset_size)++;

        //Copy elements of the array
        //Warning: make sure array doesn't go out of bounds! (should have size *vertexset_size + 1)
        memcpy(vertexset[*vertexset_size], vertexset[*vertexset_size - 1], sizeof(unsigned char) * (num_vertices_in_set - 1));

        transform_vertexset_into_canonical_form(vertexset[*vertexset_size - 1], num_vertices_in_set);
        if(num_vertices_in_set <= 4) {
            //Save vertexset_index in an array for faster access in determine_vertexset_orbits
            SQUARE vertexset_elements;
            int j;
            for(j = 0; j < 4; j++)
                if(j < num_vertices_in_set)
                    vertexset_elements[j] = vertexset[*vertexset_size - 1][j];
                else
                    vertexset_elements[j] = 0;

            (*vertexset_index)[vertexset_elements[3]][vertexset_elements[2]][vertexset_elements[1]][vertexset_elements[0]] = *vertexset_size - 1;
        }
    }
}

/**
 * Generate all sets which consist of num_vertices_in_set vertices. It is
 * assumed that the graph contains no reducible triangles.
 *
 * Important: *vertexset_size should be initialized to 0 before the first call.
 * And current_index should also be 0 in the first call.
 */
//Could often stop earlier (eg if one knows there aren't enough elements left), but hardly uses any cpu
void generate_all_vertexsets(int current_index, int num_vertices_in_set,
        unsigned char vertexset[][num_vertices_in_set], int *vertexset_size) {
    int i;
    if(current_index < num_vertices_in_set) {
        int previous = -1;
        if(current_index > 0)
            previous = vertexset[*vertexset_size][current_index - 1];
        for(i = previous + 1; i < current_number_of_vertices; i++) {
                vertexset[*vertexset_size][current_index] = i;
                generate_all_vertexsets(current_index + 1, num_vertices_in_set, vertexset, vertexset_size);
        }
    } else {
        //DEBUGARRAYDUMP(vertexset[*vertexset_size], num_vertices_in_set, "%d");
        (*vertexset_size)++;

        //Copy elements of the array
        //Warning: make sure array doesn't go out of bounds! (should have size *vertexset_size + 1)
        memcpy(vertexset[*vertexset_size], vertexset[*vertexset_size - 1], sizeof(unsigned char) * (num_vertices_in_set - 1));

        if(num_vertices_in_set <= 4) {
            SQUARE vertexset_elements;
            int j;
            for(j = 0; j < 4; j++)
                if(j < num_vertices_in_set)
                    vertexset_elements[j] = vertexset[*vertexset_size - 1][j];
                else
                    vertexset_elements[j] = 0;

            (*vertexset_index)[vertexset_elements[3]][vertexset_elements[2]][vertexset_elements[1]][vertexset_elements[0]] = *vertexset_size - 1;
        }
    }
}

/**
 * Generate all sets of vertices which consist of num_vertices_in_set vertices which eliminate all reducible triangles.
 * Only generating vertexsets which can later still become girth 4 graphs.
 *
 * Remark it is assumed that num_vertices_in_set = 2
 */
void find_vertexsets_girth4(int num_vertices_in_set, unsigned char vertexset[][num_vertices_in_set], int *vertexset_size) {
    *vertexset_size = 0;
    DEBUGASSERT(num_vertices_in_set == 2);

    DEBUGASSERT(number_of_reducible_triangles <= num_vertices_in_set);
    if(number_of_reducible_triangles > 0)
        generate_triangle_vertexsets_girth4(num_vertices_in_set, vertexset, vertexset_size);
    else
        generate_all_vertexsets_girth4(num_vertices_in_set, vertexset, vertexset_size);
}

/**
 * Generate all sets of vertices which consist of 2 vertices which eliminate all reducible triangles.
 * Only generating vertexsets which can later still become girth 4 graphs.
 *
 * Remark: it is assumed that number_of_reducible_triangles <= num_vertices_in_set.
 */
void generate_triangle_vertexsets_girth4(int num_vertices_in_set, unsigned char vertexset[][num_vertices_in_set], int *vertexset_size) {
    DEBUGASSERT(number_of_reducible_triangles <= 2);

    //Remark: must also be called if there are diamonds, otherwise not all graphs will be generated
    //DEBUGASSERT(number_of_irreducible_triangles == 0);
    if(number_of_reducible_triangles == 1) {
        //Generate all 3 edges from the triangle
        int i;
        for(i = 0; i < 2; i++) {
            vertexset[*vertexset_size][0] = reducible_triangles[0][i];
            vertexset[*vertexset_size][1] = reducible_triangles[0][i + 1];
            (*vertexset_index)[0][0][reducible_triangles[0][i + 1]][reducible_triangles[0][i]] = *vertexset_size;
            (*vertexset_size)++;
        }
        vertexset[*vertexset_size][0] = reducible_triangles[0][0];
        vertexset[*vertexset_size][1] = reducible_triangles[0][2];
        (*vertexset_index)[0][0][reducible_triangles[0][2]][reducible_triangles[0][0]] = *vertexset_size;
        (*vertexset_size)++;

        DEBUGASSERT(*vertexset_size == 3);
    }
    //No edgepairs generated if number_of_reducible_triangles > 1, because then it cant be part of 2 squares
}

/**
 * Generate all sets of vertices which consist of 2 vertices.
 * Only generating vertexsets which can later still become girth 4 graphs.
 *
 * Remark: it is assumed that current_graph contains no reducible triangles
 */
void generate_all_vertexsets_girth4(int num_vertices_in_set, unsigned char vertexset[][num_vertices_in_set], int *vertexset_size) {
    DEBUGASSERT(num_vertices_in_set == 2);
    int i, j, k;
    //Blowing up the vertices of the 5 edges which are in the diamond.
    if(current_number_of_vertices > 4) {
        for(i = 0; i < number_of_irreducible_triangles; i++) {
            for(j = 0; j < 3; j++) {
                int max;
                if(j != 0)
                    max = 4;
                else
                    max = 3;
                for(k = j + 1; k < max; k++) {
                    if(irreducible_triangles[i][j] < irreducible_triangles[i][k]) {
                        vertexset[*vertexset_size][0] = irreducible_triangles[i][j];
                        vertexset[*vertexset_size][1] = irreducible_triangles[i][k];
                    } else {
                        vertexset[*vertexset_size][0] = irreducible_triangles[i][k];
                        vertexset[*vertexset_size][1] = irreducible_triangles[i][j];
                    }
                    (*vertexset_index)[0][0][vertexset[*vertexset_size][1]][vertexset[*vertexset_size][0]] = *vertexset_size;
                    (*vertexset_size)++;
                }
            }
        }
        DEBUGASSERT(*vertexset_size == number_of_irreducible_triangles * 5);
    } else { //0 3 also needed in case of K4
        for(j = 0; j < 3; j++) {
            for(k = j + 1; k < 4; k++) {
                if(irreducible_triangles[0][j] < irreducible_triangles[0][k]) {
                    vertexset[*vertexset_size][0] = irreducible_triangles[0][j];
                    vertexset[*vertexset_size][1] = irreducible_triangles[0][k];
                } else {
                    vertexset[*vertexset_size][0] = irreducible_triangles[0][k];
                    vertexset[*vertexset_size][1] = irreducible_triangles[0][j];
                }
                (*vertexset_index)[0][0][vertexset[*vertexset_size][1]][vertexset[*vertexset_size][0]] = *vertexset_size;
                (*vertexset_size)++;
            }
        }
        DEBUGASSERT(*vertexset_size == 6);
    }
}

/**
 * Returns 1 if the vertex is part of a diamond, else 0.
 */
int is_part_of_irreducible_triangle_bitvector(int vertex) {
    return (irreducible_triangles_bitvector & BIT(vertex)) > 0;
}

/**
 * Returns 1 if vertex is part of the diamond at index diamond, else returns 0.
 */
int is_part_of_same_irreducible_triangle(unsigned char vertex, int diamond) {
    DEBUGASSERT(diamond >= 0 && diamond < number_of_irreducible_triangles);
    return irreducible_triangles[diamond][0] == vertex || irreducible_triangles[diamond][1] == vertex ||
                    irreducible_triangles[diamond][2] == vertex || irreducible_triangles[diamond][3] == vertex;
}

/**
 * Returns 1 if vertex is part of a diamond, else returns 0.
 * If 1 is returned, diamond points to the index of the diamond in the list of
 * irreducible triangles of which vertex is a part.
 */
int is_part_of_irreducible_triangle_diamond(int vertex, int *diamond) {
    if(number_of_irreducible_triangles == 0)
        return 0;
    else {
        int i;
        for(i = 0; i < number_of_irreducible_triangles; i++) {
            if(irreducible_triangles[i][0] == vertex || irreducible_triangles[i][1] == vertex ||
                    irreducible_triangles[i][2] == vertex || irreducible_triangles[i][3] == vertex) {
                *diamond = i;
                return 1;
            }
        }
        return 0;
    }
}

/**
 * Returns 1 if vertex belongs to a reducible triangle, else 0.
 * If 1 is returned, triangle points to the triangle of which vertex is a part.
 */
int is_part_of_reducible_triangle(int vertex, int *triangle) {
    if(number_of_reducible_triangles == 0)
        return 0;
    else {
        int i;
        for(i = 0; i < number_of_reducible_triangles; i++) {
            if(reducible_triangles[i][0] == vertex || reducible_triangles[i][1] == vertex ||
                    reducible_triangles[i][2] == vertex) {
                *triangle = i;
                return 1;
            }
        }
        return 0;
    }
}

/**
 * Returns 1 if there is no vertex in current_vertexset which is part of triangle
 * and which has a bigger label than vertex, else returns 0.
 *
 * Remark: it is assumed that vertex is part of "triangle"
 */
int can_be_added_to_reducible_triangle_vertexset(unsigned char current_vertexset[], int current_size, int vertex, int triangle) {
    DEBUGASSERT(triangle >= 0 && triangle < number_of_reducible_triangles);
    int i;
    for(i = 0; i < current_size;i++) {
        if(is_part_of_same_reducible_triangle(current_vertexset[i], triangle) && current_vertexset[i] >= vertex)
            return 0;
    }
    return 1;
}

int is_part_of_same_reducible_triangle(int vertex, int triangle) {
    return reducible_triangles[triangle][0] == vertex || reducible_triangles[triangle][1] == vertex || reducible_triangles[triangle][2] == vertex;
}

/**
 * Returns 1 if the edge (from, to) is a bridge, else 0.
 *
 * This is determined using a BFS algorithm.
 */
int is_a_bridge(unsigned char from, unsigned char to) {
    //Could define globally (with MAXN), but is faster as local local variable
    //number_of_vertices instead of current_number_of_vertices because current_number_of_vertices may not be incremented yet
    unsigned char queue[number_of_vertices];
    int i, j, next;

    RESETMARKS;

    int queue_size = 0;
    queue[queue_size++] = from;
    MARK(from);

    //i points to the current element in the queue
    i = 0;
    while(i < queue_size) {
        for(j = 0; j < REG; j++) {
            next = current_graph[queue[i]][j];
            if(!ISMARKED(next)) {
                if(next != to) {
                    queue[queue_size++] = next;
                    MARK(next);
                } else if(i > 0)
                    return 0;
            }

        }
        i++;
    }
    DEBUGASSERT(queue_size < number_of_vertices);

    return 1;
}

/**
 * Update the bridges after the triangle insertion. It is possible that a triangle insertion
 * modifies the label of one vertex of a bridge (but a triangle insertion can never
 * destroy a bridge).
 *
 * org_vertex is the vertex which is transformed into a triangle
 * new_vertex is the new label of this vertex in the edge (org_vertex, fixed_vertex)
 * fixed_vertex is a neighbour of org_vertex
 *
 * Returns 1 if a bridge was modified, else 0.
 */
int update_bridges_triangle_insert(unsigned char org_vertex, unsigned char new_vertex, unsigned char fixed_vertex) {
    if(number_of_bridges == 0)
        return 0;
    int i;

    for(i = 0; i < number_of_bridges; i++)
        if(bridges[i][0] == org_vertex && bridges[i][1] == fixed_vertex) {
            is_bridge[org_vertex][fixed_vertex] = 0;
            is_bridge[fixed_vertex][new_vertex] = 1;
            //Always swapped order here, because new_vertex will always be > fixed_vertex
            //bridges[i][0] = new_vertex;
            //transform_edge_into_canonical_form(bridges[i]);
            bridges[i][0] = fixed_vertex;
            bridges[i][1] = new_vertex;

            return 1;
        } else if(bridges[i][0] == fixed_vertex && bridges[i][1] == org_vertex) {
            is_bridge[fixed_vertex][org_vertex] = 0;
            is_bridge[fixed_vertex][new_vertex] = 1;
            //Never swapped order here, because new_vertex will always be > fixed_vertex
            bridges[i][1] = new_vertex;
            //transform_edge_into_canonical_form(bridges[i]);

            return 1;
        }

    return 0;
}

/**
 * Checks if the bridges of the old graph are still bridges in the modified graph.
 * The edges which are no longer bridges are removed from the list of bridges
 */
void update_bridges_add_edge() {
    int i = 0, previous_i = 0;
    while(i < number_of_bridges) {
        for(i = previous_i; i < number_of_bridges; i++)
            if(!is_a_bridge(bridges[i][0], bridges[i][1])) {
                remove_bridge(i);

                previous_i = i;
                break;
            }
    }
}

/**
 * Adds a bridge to the list of bridges.
 */
void add_bridge(unsigned char from, unsigned char to) {
    bridges[number_of_bridges][0] = from;
    bridges[number_of_bridges][1] = to;
    number_of_bridges++;
}

/**
 * Removes the bridge at that index.
 */
//TODO: could just update pointers? But doesn't use much cpu
void remove_bridge(int index) {
    DEBUGASSERT(number_of_bridges > 0);
    is_bridge[bridges[index][0]][bridges[index][1]] = 0;
    number_of_bridges--;
    if(number_of_bridges > 0 && index < number_of_bridges) {
        bridges[index][0] = bridges[number_of_bridges][0];
        bridges[index][1] = bridges[number_of_bridges][1];
    }
}
/**
 * If old_from, old_to was a bridge in the parent graph, it is
 * updated to from, to.
 */
void replace_bridge(int old_from, int old_to, int from, int to) {
    int i;
    is_bridge[old_from][old_to] = 0;
    is_bridge[from][to] = 1;
    for(i = 0; i < number_of_bridges; i++) {
        if(bridges[i][0] == old_from && bridges[i][1] == old_to) {
            bridges[i][0] = from;
            bridges[i][1] = to;
            break;
        }
    }
    DEBUGASSERT(i < number_of_bridges);
}

/**
 * Returns 1 if edge is in the list of bridges, else 0 is returned.
 *
 * From is assumed to be < than to.
 */
int is_a_bridge_list(unsigned char from, unsigned char to) {
    //Remark: this method is not very efficient, but is only used for the generation of prime graphs,
    //so it hardly consumes any cpu-time
    if(number_of_bridges == 0)
        return 0;
    int i;
    for(i = 0; i < number_of_bridges; i++)
        if(bridges[i][0] == from && bridges[i][1] == to)
            return 1;
    return 0;
}

/**
 * Returns 1 is the edge is reducible, else 0.
 *
 * Warning: doesn't check if edge is part of a reducible triangle, so shouldn't be
 * called when there are reducible triangles present.
 */
//Using a macro is not faster
int is_reducible_edge(EDGE edge) {
    /**
     * It's slightly faster to check diamonds first, since the percentage of graphs which
     * contain diamonds is bigger than the percentage of graphs which contain bridges:
     * 13 vs 1.5% (for 26 vertices)
     */
    //return !is_bridge[edge[0]][edge[1]] && (irreducible_triangles_bitvector & (BIT(edge[0]) | BIT(edge[1]))) == 0;
    return (irreducible_triangles_bitvector & (BIT(edge[0]) | BIT(edge[1]))) == 0 && !is_bridge[edge[0]][edge[1]];
}


/**
 * Returns the neighbours on distance <= 2 of "vertex"
 */
setword determine_vertex_neighbours_distance_two(unsigned char vertex) {
    DEBUGASSERT(vertex < current_number_of_vertices);
    setword neighbours = (setword) 0;
    unsigned char i;
    for(i = 0; i < REG; i++) {
        neighbours |= vertex_neighbourhood[current_graph[vertex][i]];
    }

    vertex_colours_long_two[vertex] = neighbours;
    MARK(vertex);

    return neighbours;
}

/**
 * Colour is the number of vertices on distance <= 2 of that edge.
 */
int determine_edge_colour(EDGE edge) {
    setword neighbours0, neighbours1;
    if(ISMARKED(edge[0]))
        neighbours0 = vertex_colours_long_two[edge[0]];
    else
        neighbours0 = determine_vertex_neighbours_distance_two(edge[0]);

    if(ISMARKED(edge[1]))
        neighbours1 = vertex_colours_long_two[edge[1]];
    else
        neighbours1 = determine_vertex_neighbours_distance_two(edge[1]);

    int colour = POPCOUNT(neighbours0 | neighbours1);

    DEBUGASSERT(colour <= MAX_EDGE_COLOUR_TWO);

    colours_one[edge[0]][edge[1]] = colour;

    return colour;
}

/**
 * Returns 1 if edge is part of a square, else returns 0.
 */
int edge_is_part_of_square(EDGE edge) {
    setword neighbours = (setword) 0;
    unsigned char i, neighbour;
    int number_visited = 0;
    for(i = 0; i < REG; i++) {
        neighbour = current_graph[edge[0]][i];
        if(neighbour != edge[1]) {
            neighbours |= vertex_neighbourhood[neighbour];
            if(number_visited == 1)
                break; //Is slightly faster like this
            number_visited++;
        }
    }

/*
    //Is slightly slower
    neighbours &= ~BIT(edge[0]);

    return (neighbours & vertex_neighbourhood[edge[1]]) > 0;
*/

    return (neighbours & vertex_neighbourhood[edge[1]]) != BIT(edge[0]);
}

/**
 * Returns of how much squares edge is a part.
 */
int edge_is_part_of_how_much_squares(EDGE edge) {
    setword neighbours[2];
    unsigned char neighbour, current_neighbour_index = 0;
    int i;
    for(i = 0; i < REG; i++) {
        neighbour = current_graph[edge[0]][i];
        if(neighbour != edge[1]) {
            neighbours[current_neighbour_index] = vertex_neighbourhood[neighbour];
            if(current_neighbour_index == 1)
                break;
            current_neighbour_index++;
        }
    }

    setword neighbours_edge1 = vertex_neighbourhood[edge[1]];

    neighbours_edge1 &= ~BIT(edge[0]);
    return POPCOUNT(neighbours[0] & neighbours_edge1) + POPCOUNT(neighbours[1] & neighbours_edge1);
}

/**
 * Determines the number of neighbours on distance two of the vertices which weren't
 * marked yet.
 * This is necessary for determine_edge_colour_expensive, because this method
 * reuses the information of determine_vertex_neighbours_distance_two.
 */
void determine_colours_of_remaining_vertices() {
    int i;
    for(i = 0; i < current_number_of_vertices; i++) {
        if(!ISMARKED(i))
            determine_vertex_neighbours_distance_two(i);
    }
}

/**
 * Tests if the inserted edge has minimal colour.
 * Returns 1 if it has minimal colour, else 0.
 * If 1 is returned, edge_inserted is also set to EDGE_INSERTED or MAJOR_EDGE_INSERTED
 * depending on how much edges have minimal colour.
 *
 * Important: this method assumes that current_graph has girth 5.
 *
 * The first colour is the number of vertices on distance <= 2 of an edge
 */
int has_min_colour(EDGE inserted_edge, int *edge_inserted, EDGE neighbours[4]) {
    //It is assumed that the graph has girth 5
    //DEBUGASSERT(!contains_squares());

    RESETMARKS;
    int colour_inserted_edge = determine_edge_colour(inserted_edge);

    //Checking min_edge first is slightly faster
    //is_reducible_edge(min_edge) not needed because if it's reducible in the parent graph, it will remain reducible (or will no longer be an edge)
    if(is_neighbour(min_edge[0], min_edge[1])) {
        unsigned char min_colour_new = determine_edge_colour(min_edge);
        if(min_colour_new < colour_inserted_edge)
            return 0;
        //Is slower to check number of squares here
    }

    if(is_neighbour(previous_min_edge[0], previous_min_edge[1]) && (all_edges_are_reducible || is_reducible_edge(previous_min_edge))) {
        unsigned char prev_colour = determine_edge_colour(previous_min_edge);
        if(prev_colour < colour_inserted_edge)
            return 0;
    }

    //Helps a little
    int k;
    for(k = 0; k < 4; k++) {
        if(all_edges_are_reducible || is_reducible_edge(neighbours[k])) {
            if(determine_edge_colour(neighbours[k]) < colour_inserted_edge)
                return 0;
        }
    }

    int colour;
    EDGE edge;

    min_edges_size = 0;

    unsigned char i, j;
    unsigned char neighbour;
    //current_number_of_vertices - 2 to make sure inserted edge won't be put twice in the list
    for(i = 0; i < current_number_of_vertices - 2; i++) {
        for(j = 0; j < REG; j++) {
            if(i < current_graph[i][j]) {
                neighbour = current_graph[i][j];
                edge[0] = i;
                edge[1] = neighbour;

                if(all_edges_are_reducible || is_reducible_edge(edge)) {
                    colour = determine_edge_colour(edge);
                    if(colour < colour_inserted_edge) {
                        previous_min_edge[0] = i;
                        previous_min_edge[1] = neighbour;

                        return 0;
                    } else if(colour == colour_inserted_edge) {
                        min_edges[min_edges_size][0] = i;
                        min_edges[min_edges_size][1] = neighbour;
                        min_edges_size++;
                    } else {
                        colours_two[i][neighbour] = 0;
                        colours_three[i][neighbour] = 0;
                    }
                } else { //Irreducible edges have colour 0
                    //Needed because otherwise the "third" colour will be incorrect
                    colours_three[i][neighbour] = 0;
                    colours_two[i][neighbour] = 0;
                    colours_one[i][neighbour] = 0;
                }
            }
        }
    }

    /**
     * Imporant remark, the inserted edge (i.e. {n-2, n-1}) will always be
     * the last edge in min_edges.
     */
    min_edges[min_edges_size][0] = inserted_edge[0];
    min_edges[min_edges_size][1] = inserted_edge[1];
    min_edges_size++;

    min_colour_one = colour_inserted_edge;
    min_edge_is_part_of_square = 0;
    if(min_edges_size > 1) {
        //Needed because vertex_colours_long_two must be up to date
        determine_colours_of_remaining_vertices();
        return has_min_colour_distance_three(inserted_edge, edge_inserted, min_edges, min_edges_size);
    } else {
        //In other case min_edge will be set by the extend method, but not here because of the INIT
        if(current_number_of_vertices < number_of_vertices) {
            min_edge[0] = inserted_edge[0];
            min_edge[1] = inserted_edge[1];
        }
        *edge_inserted = MAJOR_EDGE_INSERTED;
    }
    return 1;
}

/**
 * Tests if the inserted edge has minimal colour.
 * Returns 1 if it has minimal colour, else 0.
 *
 * The colour here is the number of vertices on distance <= 2 of an edge.
 * The colour with the minimal frequency is the canonical colour.
 *
 * Important: it is assumed that all edges in min_colour_edges are part of a square.
 * And the inserted edge is always the last edge of min_colour_edges.
 */
int has_min_colour_distance_two(EDGE inserted_edge, int *edge_inserted, EDGE min_colour_edges[], int number_of_min_colour_edges) {
    //Already done by has_min_colour_cycle
    //RESETMARKS;
    int colour_inserted_edge = determine_edge_colour(inserted_edge);
    int is_part_of_num_squares = -1;

    //Previous_min_edge also already done by has_min_colour_cycle

    int colour;
    min_edges_size = 0;

    int i;
    //number_of_min_colour_edges - 1 because the inserted edge is always the last one in the list
    for(i = 0; i < number_of_min_colour_edges - 1; i++) {
        colour = determine_edge_colour(min_colour_edges[i]);
        if(colour < colour_inserted_edge) {
            previous_min_edge[0] = min_colour_edges[i][0];
            previous_min_edge[1] = min_colour_edges[i][1];

            return 0;
        } else if(colour == colour_inserted_edge) {
            if(is_part_of_num_squares == -1) {
                is_part_of_num_squares = edge_is_part_of_how_much_squares(inserted_edge);
            }
            int is_part_of_num_squares_other = edge_is_part_of_how_much_squares(min_colour_edges[i]);
            if(is_part_of_num_squares == is_part_of_num_squares_other) {
                min_edges[min_edges_size][0] = min_colour_edges[i][0];
                min_edges[min_edges_size][1] = min_colour_edges[i][1];
                min_edges_size++;
            } else if(is_part_of_num_squares > is_part_of_num_squares_other) {
                //Needed because otherwise the "third" colour will be incorrect
                colours_two[min_colour_edges[i][0]][min_colour_edges[i][1]] = 0;
                colours_three[min_colour_edges[i][0]][min_colour_edges[i][1]] = 0;
            } else { //part_of_square_inserted < part_of_square_other
                DEBUGASSERT(is_part_of_num_squares < is_part_of_num_squares_other);
/*
                previous_min_edge[0] = min_colour_edges[i][0];
                previous_min_edge[1] = min_colour_edges[i][1];
*/

                return 0;
            }
        } else {
            colours_two[min_colour_edges[i][0]][min_colour_edges[i][1]] = 0;
            colours_three[min_colour_edges[i][0]][min_colour_edges[i][1]] = 0;
        }
    }

    /**
     * Imporant remark, the inserted edge (i.e. {n-2, n-1}) will always be
     * the last edge in min_edges.
     */
    min_edges[min_edges_size][0] = inserted_edge[0];
    min_edges[min_edges_size][1] = inserted_edge[1];
    min_edges_size++;

    min_colour_one = colour_inserted_edge;
    if(min_edges_size > 1) {

        DEBUGASSERT(is_part_of_num_squares > -1);
        //Needed because vertex_colours_long_two must be up to date
        determine_colours_of_remaining_vertices();
        return has_min_colour_distance_three(inserted_edge, edge_inserted, min_edges, min_edges_size);
    } else {
        //In other case min_edge will be set by the extend method, but not here because of the INIT
        if(current_number_of_vertices < number_of_vertices) {
            min_edge[0] = inserted_edge[0];
            min_edge[1] = inserted_edge[1];
        }
        *edge_inserted = MAJOR_EDGE_INSERTED;
    }
    return 1;
}

/**
 * Tests if the inserted edge has minimal colour.
 * Returns 1 if it has minimal colour, else 0.
 * If 1 is returned, edge_inserted is also set to EDGE_INSERTED or MAJOR_EDGE_INSERTED
 * depending on how much edges have minimal colour.
 *
 * Important: this method that the inserted edge is part of a square.
 *
 * Minimal edges must be part of a square and have as few as possible vertices on
 * distance <= 2.
 */
int has_min_colour_cycle(EDGE inserted_edge, int *edge_inserted) {
    DEBUGASSERT(girth < 5 || current_number_of_vertices < number_of_vertices);

    //For determine_edge_colour
    RESETMARKS;

    int part_of_num_squares_inserted_edge = edge_is_part_of_how_much_squares(inserted_edge);
    DEBUGASSERT(part_of_num_squares_inserted_edge > 0);
    //DEBUGASSERT(edge_is_part_of_square(inserted_edge));

    int colour_inserted_edge = determine_edge_colour(inserted_edge);

    //is_reducible_edge(min_edge) not needed because if it's reducible in the parent graph, it will remain reducible (or will no longer be an edge)
    if(is_neighbour(min_edge[0], min_edge[1])) {
        if(edge_is_part_of_square(min_edge)) {
            if(determine_edge_colour(min_edge) < colour_inserted_edge)
                return 0;
            //Is slower to check number of squares here
        }
    }

    //Fastest is checking min_edge first
    if(is_neighbour(previous_min_edge[0], previous_min_edge[1]) && (all_edges_are_reducible || is_reducible_edge(previous_min_edge))) {
        if(edge_is_part_of_square(previous_min_edge)) {
            if(determine_edge_colour(previous_min_edge) < colour_inserted_edge) {
                return 0;
            }
        }
    }

    //Is slower when checking neighbours here!

    if(girth == 5 && current_number_of_vertices == number_of_vertices - 2 && current_number_of_vertices != 8) {
        return has_min_colour_cycle_girth_5(inserted_edge, edge_inserted);
    } else {
        EDGE edge;

        //For the lookahead
        min_edges_size = 0;

        //EDGE min_colour_edges[current_number_of_edges];
        //int number_of_min_colour_edges = 0;

        unsigned char i, j;
        unsigned char neighbour;
        int part_of_num_squares, colour;
        //current_number_of_vertices - 2 to make sure inserted edge won't be put twice in the list
        for(i = 0; i < current_number_of_vertices - 2; i++) {
            for(j = 0; j < REG; j++) {
                if(i < current_graph[i][j]) {
                    //if(i < neighbour && !(i == inserted_edge[0] && neighbour == inserted_edge[1])) {
                    neighbour = current_graph[i][j];
                    edge[0] = i;
                    edge[1] = neighbour;

                    if(all_edges_are_reducible || is_reducible_edge(edge)) {
                        //part_of_num_squares = edge_is_part_of_how_much_squares(edge);
                        if(edge_is_part_of_square(edge)) {
                            colour = determine_edge_colour(edge);
                            if(colour < colour_inserted_edge) {
                                previous_min_edge[0] = i;
                                previous_min_edge[1] = neighbour;

                                return 0;
                            } else if(colour == colour_inserted_edge) {
                                //Ok, is better to test here, than to postpone this till the colours of all edges have been calculated.
                                //(And then filtering the ones which are part of too few or too much squares)
                                part_of_num_squares = edge_is_part_of_how_much_squares(edge);
                                if(part_of_num_squares == part_of_num_squares_inserted_edge) {
                                    min_edges[min_edges_size][0] = i;
                                    min_edges[min_edges_size][1] = neighbour;
                                    min_edges_size++;
                                } else if(part_of_num_squares < part_of_num_squares_inserted_edge) {
                                    //and !part_of_square_other

                                    //Needed because otherwise the "third" colour will be incorrect
                                    colours_two[i][neighbour] = 0;
                                    colours_three[i][neighbour] = 0;
                                } else {
                                    //Is slightly faster when setting this, so ok!
                                    previous_min_edge[0] = i;
                                    previous_min_edge[1] = neighbour;

                                    return 0;
                                }
                            } else {
                                colours_two[i][neighbour] = 0;
                                colours_three[i][neighbour] = 0;
                            }
                        } else {
                            colours_one[i][neighbour] = 0;
                            colours_two[i][neighbour] = 0;
                            colours_three[i][neighbour] = 0;
                        }
                    } else { //Irreducible edges have colour 0
                        //Needed because otherwise the "third" colour will be incorrect
                        colours_three[i][neighbour] = 0;
                        colours_two[i][neighbour] = 0;
                        colours_one[i][neighbour] = 0;
                    }
                }
            }
        }

        /**
         * Imporant remark, the inserted edge (i.e. {n-2, n-1}) will always be
         * the last edge in min_edges.
         */
        min_edges[min_edges_size][0] = inserted_edge[0];
        min_edges[min_edges_size][1] = inserted_edge[1];
        min_edges_size++;

        min_colour_one = colour_inserted_edge;
        min_edge_is_part_of_square = 1;
        if(min_edges_size > 1) {
            //Needed because vertex_colours_long_two must be up to date
            determine_colours_of_remaining_vertices();
            return has_min_colour_distance_three(inserted_edge, edge_inserted, min_edges, min_edges_size);
        } else {
            //In other case min_edge will be set by the extend method, but not here because of the INIT
            if(current_number_of_vertices < number_of_vertices) {
                min_edge[0] = inserted_edge[0];
                min_edge[1] = inserted_edge[1];
            }
            *edge_inserted = MAJOR_EDGE_INSERTED;
            return 1;
        }
    }
}

/**
 * Tests if the inserted edge has minimal colour.
 * Returns 1 if it has minimal colour, else 0.
 * If 1 is returned, edge_inserted is also set to EDGE_INSERTED or MAJOR_EDGE_INSERTED
 * depending on how much edges have minimal colour.
 *
 * Important: this method assumes that current_graph contains squares and is
 * on level n-2 of girth == 5 and that the inserted edge is part of a square.
 *
 * Minimal edges must be part of a square and have as few as possible vertices on
 * distance <= 2.
 */
int has_min_colour_cycle_girth_5(EDGE inserted_edge, int *edge_inserted) {
    DEBUGASSERT(girth == 5 && current_number_of_vertices == number_of_vertices - 2 && current_number_of_vertices != 8);

    //Previous min edge was already checked

    /**
     * Not using has_min_colour_cycle here, because first checking if there aren't
     * too much squares and only afterwards calculating colour_one is more efficient
     * in case of girth == 5.
     */

    EDGE edge;

    //For the lookahead
    min_edges_size = 0;

    EDGE min_colour_edges[current_number_of_edges];
    int number_of_min_colour_edges = 0;
    int number_of_edges_part_of_square = 0;

    unsigned char i, j;
    unsigned char neighbour;
    //current_number_of_vertices - 2 to make sure inserted edge won't be put twice in the list
    for(i = 0; i < current_number_of_vertices - 2; i++) {
        for(j = 0; j < REG; j++) {
            if(i < current_graph[i][j]) {
                //if(i < neighbour && !(i == inserted_edge[0] && neighbour == inserted_edge[1])) {
                neighbour = current_graph[i][j];
                edge[0] = i;
                edge[1] = neighbour;

                if(all_edges_are_reducible || is_reducible_edge(edge)) {
                    if(edge_is_part_of_square(edge)) {
                        if(number_of_edges_part_of_square > 7) {
                            //previous min edge not modified!
                            return 0;
                        }
                        number_of_edges_part_of_square++;
                        min_colour_edges[number_of_min_colour_edges][0] = i;
                        min_colour_edges[number_of_min_colour_edges][1] = neighbour;
                        number_of_min_colour_edges++;
                    } else {
                        colours_one[i][neighbour] = 0;
                        colours_two[i][neighbour] = 0;
                        colours_three[i][neighbour] = 0;
                    }
                } else { //Irreducible edges have colour 0
                    //Needed because otherwise the "third" colour will be incorrect
                    colours_three[i][neighbour] = 0;
                    colours_two[i][neighbour] = 0;
                    colours_one[i][neighbour] = 0;
                }
            }
        }
    }

    /**
     * Imporant remark, the inserted edge (i.e. {n-2, n-1}) will always be
     * the last edge in min_edges.
     */
    min_colour_edges[number_of_min_colour_edges][0] = inserted_edge[0];
    min_colour_edges[number_of_min_colour_edges][1] = inserted_edge[1];
    number_of_min_colour_edges++;

    //This method only called if the graph contains squares!
    DEBUGASSERT(number_of_edges_part_of_square > 0);

    //Also checking g4 && snarks here is more expensive because in that case is_colourable has to be called
    //while edge is mostly rejected by more expensive colour

    min_edge_is_part_of_square = 1;
    //min_colour_one = colour_inserted_edge;
    if(number_of_min_colour_edges > 1) {
        return has_min_colour_distance_two(inserted_edge, edge_inserted, min_colour_edges, number_of_min_colour_edges);
    } else {
        //Never happens because a reducible edge is part of a square, there will be at least 4 number_of_min_colour_edges
        fprintf(stderr, "Error: there is only one reducible edge which is part of a square\n");
        exit(1);

/*
        //In other case min_edge will be set by the extend method, but not here because of the INIT
        if(current_number_of_vertices < number_of_vertices) {
            min_edge[0] = inserted_edge[0];
            min_edge[1] = inserted_edge[1];
        }

        min_colour_one = MAX_EDGE_COLOUR_TWO;


        //To make sure the partitions will be ok
        //The partitions won't be that precise, but this case only rarely happens
        colours_one[inserted_edge[0]][inserted_edge[1]] = UCHAR_MAX;

        *edge_inserted = MAJOR_EDGE_INSERTED;

        return 1;
*/
    }
}

/**
 * Returns the neighbours on distance <= 3 of "vertex".
 *
 * Remark: information of vertex_colours_long_two is reused, so the content of that
 * array should be valid.
 */
setword determine_vertex_neighbours_distance_three(unsigned char vertex) {
    DEBUGASSERT(vertex < current_number_of_vertices);
    unsigned char i;

    setword neighbours = vertex_colours_long_two[vertex];
    for(i = 0; i < REG; i++) {
        neighbours |= vertex_colours_long_two[current_graph[vertex][i]];
    }

    vertex_colours_long_three[vertex] = neighbours;
    MARK(vertex);

    return neighbours;
}

/**
 * Determines the number of vertices on distance <= 3 of the edge.
 */
int determine_edge_colour_expensive(EDGE edge) {
    setword neighbours0, neighbours1;
    if(ISMARKED(edge[0]))
        neighbours0 = vertex_colours_long_three[edge[0]];
    else
        neighbours0 = determine_vertex_neighbours_distance_three(edge[0]);

    if(ISMARKED(edge[1]))
        neighbours1 = vertex_colours_long_three[edge[1]];
    else
        neighbours1 = determine_vertex_neighbours_distance_three(edge[1]);

    int colour = POPCOUNT(neighbours0 | neighbours1);

    DEBUGASSERT(colour <= MAX_EDGE_COLOUR_DISTANCE_THREE);

    colours_two[edge[0]][edge[1]] = colour;

    return colour;
}

/**
 * Tests if the inserted edge has minimal colour.
 * Returns 1 if it has minimal colour, else 0.
 *
 * The colour here is the number of vertices on distance <= 3 of an edge.
 * The colour with the minimal frequency is the canonical colour.
 *
 * Important: the inserted edge is always the last edge of min_colour_edges.
 */
//Important: inserted_edge is included in min_colour_edges
int has_min_colour_distance_three(EDGE inserted_edge, int *edge_inserted, EDGE min_colour_edges[], int number_of_min_colour_edges) {
    /**
     * First calculating the expensive colour of all edges. The colour with the minimal frequency,
     * is the canonical colour.
     * This is much faster (5%) than just taking the minimal colour as canonical colour.
     */
    RESETMARKS;
    int i;
    for(i = MIN_EDGE_COLOUR_DISTANCE_THREE; i < MAX_EDGE_COLOUR_DISTANCE_THREE + 1; i++) {
        frequencies[i] = 0;
    }

    int colour;
    for(i = 0; i < number_of_min_colour_edges; i++) {
        colour = determine_edge_colour_expensive(min_colour_edges[i]);
        frequencies[colour]++;
    }

    //int expensive_colour_inserted_edge = determine_edge_colour_expensive(inserted_edge);
    int expensive_colour_inserted_edge = colour; //Inserted edge is last one in the list

    unsigned char min_frequency = UCHAR_MAX;
    //If two sets have the same frequency, the set with the smallest colour will be selected
    for(i = MIN_EDGE_COLOUR_DISTANCE_THREE; i < MAX_EDGE_COLOUR_DISTANCE_THREE + 1; i++) {
        if(frequencies[i] > 0 && frequencies[i] < min_frequency) {
            min_frequency = frequencies[i];
            min_colour_two = i;
        }
    }

    if(expensive_colour_inserted_edge != min_colour_two) {
        return 0;
    } else if(frequencies[expensive_colour_inserted_edge] > 1) {
        EDGE min_colour_edges_new[number_of_min_colour_edges];
        int number_of_min_colour_edges_new = 0;
        unsigned char v1, v2;
        for(i = 0; i < number_of_min_colour_edges;i++) {
            v1 = min_colour_edges[i][0];
            v2 = min_colour_edges[i][1];
            if(colours_two[v1][v2] == min_colour_two) {
                min_colour_edges_new[number_of_min_colour_edges_new][0] = v1;
                min_colour_edges_new[number_of_min_colour_edges_new][1] = v2;

                number_of_min_colour_edges_new++;
            } else {
                colours_three[v1][v2] = 0;
            }
        }

        return has_min_colour_combination(inserted_edge, edge_inserted, min_colour_edges_new, number_of_min_colour_edges_new);
    } else {
        //In other case min_edge will be set by the extend method, but not here because of the INIT
        if(current_number_of_vertices < number_of_vertices) {
            min_edge[0] = inserted_edge[0];
            min_edge[1] = inserted_edge[1];
        }
        *edge_inserted = MAJOR_EDGE_INSERTED; //Inserted edge is only one in the smallest set

        return 1;
    }
}

/**
 * This colour is the sum of colours_one and colours_two of the 4 outgoing edges of "edge" + 1.
 */
int determine_edge_colour_combination(EDGE edge) {
    //Otherwise colour might be 0 and then it isn't possible to destinguish this from an irreducible edge
    int colour = 1;
    int colour2 = 1;

    EDGE temp_edge;
    int i, j;
    for(i = 0; i < 2; i++) {
        for(j = 0; j < REG; j++) {
            if(current_graph[edge[i]][j] != edge[(i + 1) % 2]) {
                temp_edge[0] = edge[i];
                temp_edge[1] = current_graph[edge[i]][j];

                //Could avoid this, but transform_edge_into_canonical_form is not often called
                transform_edge_into_canonical_form(temp_edge);

                /**
                 * Colour1 is actually useful, because 2 edges may have the same
                 * value for colours_two, but a different sum of colours_one.
                 */
                colour += colours_one[temp_edge[0]][temp_edge[1]];
                colour2 += colours_two[temp_edge[0]][temp_edge[1]];
            }
        }
    }

    DEBUGASSERT(colour <= 57); //4*14 + 1
    DEBUGASSERT(colour2 <= 121);

    int colour_final = colour + colour2;

    DEBUGASSERT(colour_final <= MAX_EDGE_COLOUR_COMBINATION);

    colour = colour_final;

    colours_three[edge[0]][edge[1]] = colour;

    return colour;
}

/**
 * Tests if the inserted edge has minimal colour.
 * Returns 1 if it has minimal colour, else 0.
 *
 * Important: the inserted edge is always the last edge of min_colour_edges.
 */
int has_min_colour_combination(EDGE inserted_edge, int *edge_inserted, EDGE min_colour_edges[], int number_of_min_colour_edges) {
    int expensive_colour_inserted_edge = determine_edge_colour_combination(inserted_edge);

    int i;
    int colour;
    edgelist_size = 0;
    for(i = 0; i < number_of_min_colour_edges - 1; i++) {
        colour = determine_edge_colour_combination(min_colour_edges[i]);
        if(colour < expensive_colour_inserted_edge) {
            return 0;
        } else if(colour == expensive_colour_inserted_edge) {
            edgelist[edgelist_size][0] = min_colour_edges[i][0];
            edgelist[edgelist_size][1] = min_colour_edges[i][1];

            //Could wait to put it in list till it is sure that it isn't rejected
            //But this doesn't make any difference in the timing
            edge_index[min_colour_edges[i][0]][min_colour_edges[i][1]] = edgelist_size;
            edgelist_size++;
        }
    }

    /**
     * Imporant remark, the inserted edge (i.e. {n-2, n-1}) must always be
     * the last edge in edgelist.
     */
    edgelist[edgelist_size][0] = min_colour_edges[i][0];
    edgelist[edgelist_size][1] = min_colour_edges[i][1];
    edge_index[min_colour_edges[i][0]][min_colour_edges[i][1]] = edgelist_size;
    edgelist_size++;

    min_colour_three = expensive_colour_inserted_edge;
    if(edgelist_size > 1) {
        *edge_inserted = EDGE_INSERTED; //More than one in smallest set
    } else {
        //In other case min_edge will be set by the extend method, but not here because of the INIT
        if(current_number_of_vertices < number_of_vertices) {
            min_edge[0] = inserted_edge[0];
            min_edge[1] = inserted_edge[1];
        }
        *edge_inserted = MAJOR_EDGE_INSERTED; //Inserted edge is only one in the smallest set
    }

    return 1;
}

/**
 * Determines orbits of edges.
 *
 * Uses union by size with path compression for this.
 */
void determine_edge_orbits(EDGE edgelist[], int edgelist_size, int edge_orbits[], int *number_of_orbits) {
    int i, j, k, temp;
    EDGE edge;

    int root_orbits_size[edgelist_size+1];
    for(i = 0; i < edgelist_size; i++) {
        edge_orbits[i] = i;
        root_orbits_size[i] = 1;
    }
    *number_of_orbits = edgelist_size;

    permutation *perm;
    for(i = 0; i < number_of_generators; i++) {
        if(*number_of_orbits == 1)
            break;
        perm = generators[i];
        for(j = 0; j < edgelist_size; j++) {
            if(*number_of_orbits == 1)
                break;

            int changed = 0;
            for(k = 0; k < 2; k++) {
                edge[k] = perm[edgelist[j][k]];
                if(!changed && perm[edgelist[j][k]] != edgelist[j][k])
                    changed = 1;
            }
            if(!changed) //Go to the next iteration if the edge wasn't modified by the generator
                continue;

            // Transform edgepair into canonical form
            if(edge[0] > edge[1]) {
                temp = edge[1];
                edge[1] = edge[0];
                edge[0] = temp;
            }

            //Is always the case except when generating prime graphs
            //DEBUGASSERT(colours_three[edge[0]][edge[1]] == min_colour_three);

            union_elements(edge_orbits, root_orbits_size, number_of_orbits, j, edge_index[edge[0]][edge[1]]);
        }
    }

    //Make sure edge_orbits[] really points to the roots
    //Important: this is necessary for is_major_edge
    if(*number_of_orbits > 1) {
        for(i = 0; i < edgelist_size; i++)
            edge_orbits[i] = find_root(edge_orbits, i);
    } else {
        for(i = 0; i < edgelist_size; i++)
            edge_orbits[i] = 0;
    }
}

/**
 * Determines the orbit of edge pairs.
 *
 * Uses union by size with path compression for this.
 */
void determine_edgepair_orbits(EDGEPAIR edge_pairs_list[], int edge_pair_list_size, int edgepair_orbits[], int *number_of_orbits) {
    int i, j, k;

    EDGEPAIR edgepair;

    int root_orbits_size[edge_pair_list_size+1];
    for(i = 0; i < edge_pair_list_size; i++) {
        edgepair_orbits[i] = i;
        root_orbits_size[i] = 1;
    }
    *number_of_orbits = edge_pair_list_size;

    //Is needed e.g. in case of 1 diamond on last level
    if(edge_pair_list_size == 1)
        return;

    permutation *perm;
    int index0, index1;
    for(i = 0; i < number_of_generators; i++) {
        if(*number_of_orbits == 1)
            break;
        perm = generators[i];
        for(j = 0; j < edge_pair_list_size; j++) {
            if(*number_of_orbits == 0)
                break;
            for(k = 0; k < 4; k++) {
                edgepair[k] = perm[edge_pairs_list[j][k]];
            }
            transform_edgepair_into_canonical_form(edgepair);
            if(edgepair[0] == edge_pairs_list[j][0] && edgepair[1] == edge_pairs_list[j][1] &&
                    edgepair[2] == edge_pairs_list[j][2] && edgepair[3] == edge_pairs_list[j][3])
                continue;


            index0 = edge_labels[edgepair[0]][edgepair[1]];
            //DEBUGASSERT(edge_labels[edgepair[1]][edgepair[0]] == index0);

            index1 = edge_labels[edgepair[2]][edgepair[3]];
            //DEBUGASSERT(edge_labels[edgepair[3]][edgepair[2]] == index1);

            /**
             * Remark: it never happens that an edgepair is in the same orbit
             * of a removed edgepair (because of colourability or new_edge_has_min_colour).
             * Because if that would be the case, this means that the edgepair
             * is part of a colour cycle, so then it was already removed before!
             * (because all edgepairs part of even cycles are removed in
             * remove_edgepairs_from_even_cycle.
             */
            k = edgepair_index[index0][index1];

            if(k < edge_pair_list_size) {
                union_elements(edgepair_orbits, root_orbits_size, number_of_orbits, j, k);
            } else {
                fprintf(stderr, "Error: transformed edgepair not found in edgepair list\n");
                printgraph();
                fprintf(stderr, "Original edgepair: %d %d %d %d\n", edge_pairs_list[j][0], edge_pairs_list[j][1], edge_pairs_list[j][2], edge_pairs_list[j][3]);
                fprintf(stderr, "Transformed edgepair: %d %d %d %d\n", edgepair[0], edgepair[1], edgepair[2], edgepair[3]);
                exit(1);
            }

        }
    }

    //Orbits don't have to point to actual roots since only the elements with orbit[i] = i are expanded

}

/**
 * Determines the "orbits" of vertexset. Those are not necessarily all "real" orbits
 * (in the sense of automorphism groups):
 * A graph with an extended diamond can be obtained in 2 ways. And the vertexsets
 * that would yield that graph are not necessarily in the same orbit.
 *
 * Uses union by size with path compression for this.
 */
void determine_vertexset_orbits(permutation generators_local[][MAXN], int number_of_generators_local,
        int num_vertices_in_set, unsigned char vertexset[][num_vertices_in_set], int vertexset_size, int vertexset_orbits[], int *number_of_orbits) {
    int i, j, k, l;
    unsigned char single_vertexset[num_vertices_in_set];
    SQUARE vertexset_elements;

    int root_orbits_size[vertexset_size+1];
    for(i = 0; i < vertexset_size; i++) {
        vertexset_orbits[i] = i;
        root_orbits_size[i] = 1;
    }
    *number_of_orbits = vertexset_size;

    permutation *perm;
    for(i = 0; i < number_of_generators_local; i++) {
        if(*number_of_orbits == 1)
            break;

        perm = generators_local[i];
        for(j = 0; j < vertexset_size; j++) {
            if(*number_of_orbits == 1)
                break;

            for(k = 0; k < num_vertices_in_set; k++) {
                single_vertexset[k] = perm[vertexset[j][k]];
            }

            transform_vertexset_into_canonical_form(single_vertexset, num_vertices_in_set);

            for(k = 0; k < num_vertices_in_set; k++) {
                if(single_vertexset[k] != vertexset[j][k]) {
                    break;
                }
            }
            if(k == num_vertices_in_set) //Go to the next iteration if the vertexset wasn't modified by the generator
                continue;

            if(num_vertices_in_set <= 4) {
                /**
                 * Cheap way to determine the index in the list of vertexsets of the transformed vertex.
                 * Can't do this for vertexsets with more vertices, because this would require more memory.
                 * But usually the size of a verteset is <= 4, so this shouldn't be a problem.
                 */
                for(l = 0; l < 4; l++) {
                    if(l < num_vertices_in_set)
                        vertexset_elements[l] = single_vertexset[l];
                    else
                        vertexset_elements[l] = 0;
                }
                //Leading zeros in case of num_vertices_in_set < 4. In that case the calculation of the addresses is cheaper
                k = (*vertexset_index)[vertexset_elements[3]][vertexset_elements[2]][vertexset_elements[1]][vertexset_elements[0]];
            } else {
                for(k = 0; k < vertexset_size; k++) {
                    for(l = 0; l < num_vertices_in_set; l++) {
                        if(single_vertexset[l] != vertexset[k][l])
                            break;
                    }
                    if(l == num_vertices_in_set) //Vertexset found
                        break;
                }
            }
            if(k < vertexset_size) {
                union_elements(vertexset_orbits, root_orbits_size, number_of_orbits, j, k);
                //break; //break not allowed!
            } else {
                fprintf(stderr, "Error: transformed vertexset not found in vertexset list\n");
                printgraph();
                exit(1);
            }
        }
    }

/*
    //Not needed: if the vertexset doesn't point to itself, it will never point to itself
    if(*number_of_orbits > 1) {
        for(i = 0; i < vertexset_size; i++)
             vertexset_orbits[i] = find_root(vertexset_orbits, i);
    }
*/

    /**
     * If girth4 and num_vertices_in_set = 2, then either both vertices of the set are part of the same diamond
     * or none of them are part of a diamond. This is because of the colour requirements.
     * So no need to apply the special "generators" in that case.
     */
    if(!(girth > 3 && num_vertices_in_set == 2)) {
        //Apply special "generators". This is to make sure that graphs with "expanded diamonds" won't be generated multiple times
        for(i = 0; i < number_of_irreducible_triangles; i++) {
            if(*number_of_orbits == 1)
                break;

            for(j = 0; j < vertexset_size; j++) {
                if(*number_of_orbits == 1)
                    break;
                //if(vertexset_orbits[j] == j) { //Don't uncomment this, otherwise too few orbits found (only occurs from 32 3)
                    if(apply_special_generator(vertexset[j], num_vertices_in_set, i, single_vertexset)) {
                        transform_vertexset_into_canonical_form(single_vertexset, num_vertices_in_set);

                        if(num_vertices_in_set <= 4) {
                            for(l = 0; l < 4; l++) {
                                if(l < num_vertices_in_set)
                                    vertexset_elements[l] = single_vertexset[l];
                                else
                                    vertexset_elements[l] = 0;
                            }
                            //Leading zeros in case of num_vertices_in_set < 4. In that case the calculation of the addresses is cheaper
                            k = (*vertexset_index)[vertexset_elements[3]][vertexset_elements[2]][vertexset_elements[1]][vertexset_elements[0]];
                        } else {
                            for(k = 0; k < vertexset_size; k++) {
                                for(l = 0; l < num_vertices_in_set; l++) {
                                    if(single_vertexset[l] != vertexset[k][l])
                                        break;
                                }
                                if(l == num_vertices_in_set) //Vertexset found
                                    break;
                            }
                        }
                        if(k < vertexset_size) {
                            union_elements(vertexset_orbits, root_orbits_size, number_of_orbits, j, k);
                        } else {
                            fprintf(stderr, "Error: transformed vertexset not found in vertexset list (2)\n");
                            printgraph();
                            exit(1);
                        }
                    }
                //}
            }
        }
    }


    /**
     * Is only used once (for K4). Otherwise the orbits would point to 2 instead of 0.
     * (And this would cause a problem in transform vertexset into triangles).
     */
    if(*number_of_orbits <= 1 && current_number_of_vertices == 4 && num_vertices_in_set == 1) {
        for(i = 0; i < vertexset_size; i++)
            vertexset_orbits[i] = 0;
    }
}

/**
 * Applies the special "generator" for extended diamonds to the vertexset.
 *
 * Returns 1 if the org_vertexset was modified, else 0. The org_vertexset can only
 * be modified if it contains a vertex at index 0 or 3 of the diamond (i.e. is an
 * extremal vertex) and that vertex must be the only vertex which is part of the diamond.
 */
int apply_special_generator(unsigned char org_vertexset[], int num_vertices_in_set, int diamond, unsigned char permutated_vertexset[]) {
    int i, j;
    for(i = 0; i < num_vertices_in_set; i++) {
        if(org_vertexset[i] == irreducible_triangles[diamond][0]) {
            if(is_only_vertex_from_diamond(org_vertexset, num_vertices_in_set, diamond, 0)) {
                //Could use memcpy for this
                for(j = 0; j < num_vertices_in_set; j++)
                    permutated_vertexset[j] = org_vertexset[j];
                permutated_vertexset[i] = irreducible_triangles[diamond][3];
                return 1;
            } else
                return 0;
        } else if(org_vertexset[i] == irreducible_triangles[diamond][3]) {
            if(is_only_vertex_from_diamond(org_vertexset, num_vertices_in_set, diamond, 3)) {
                for(j = 0; j < num_vertices_in_set; j++)
                    permutated_vertexset[j] = org_vertexset[j];
                permutated_vertexset[i] = irreducible_triangles[diamond][0];
                return 1;
            } else
                return 0;
        }
    }
    return 0;
}

/**
 * Returns 1 if the vertexset contains other vertices of the diamond (other than the vertex at index index of the diamond),
 * else returns 0.
 *
 * Remark: index can be 0 or 3.
 */
int is_only_vertex_from_diamond(unsigned char vertexset[], int num_vertices_in_set, int irred_triangle, int index) {
    DEBUGASSERT(irred_triangle >= 0 && irred_triangle < number_of_irreducible_triangles);
    DEBUGASSERT(index >= 0 && index <= 3);
    int i;
    for(i = 0; i < num_vertices_in_set;i++)
        if(vertexset[i] == irreducible_triangles[irred_triangle][(index + 1) % 4] || vertexset[i] == irreducible_triangles[irred_triangle][(index + 2) % 4] ||
                vertexset[i] == irreducible_triangles[irred_triangle][(index + 3) % 4])
            return 0;
    return 1;
}

/**
 * Unions the roots of a and b (if they weren't already equal) using union by size.
 */
void union_elements(int orbits[], int root_orbits_size[], int *number_of_orbits, int a, int b) {
    //DEBUGASSERT(a != b);
    int root_a = find_root(orbits, a);
    int root_b = find_root(orbits, b);
    //fprintf(stderr, "Union %d and %d (roots: %d and %d)\n", a, b, root_a, root_b);
    if(root_a != root_b) {
        if(root_orbits_size[root_a] < root_orbits_size[root_b]) {
            orbits[root_a] = root_b;
            root_orbits_size[root_b] += root_orbits_size[root_a];
        } else {
            orbits[root_b] = root_a;
            root_orbits_size[root_a] += root_orbits_size[root_b];
        }
        (*number_of_orbits)--;
    }
/*
    else {
        fprintf(stderr, "root of %d and %d are equal: %d\n", a, b, root_a);
    }
*/
}

/**
 * Find with simple path compression.
 */
int find_root(int orbits[], int i) {
    while(i != orbits[i]) {
        orbits[i] = orbits[orbits[i]]; //Simple variant of path compression
        i = orbits[i];
    }
    return i;

}

/**
 * Transforms an edge into its canonical form: edge[0] < edge[1]
 */
void transform_edge_into_canonical_form(EDGE edge) {
    if(edge[0] > edge[1]) {
        unsigned char temp;
        temp = edge[0];
        edge[0] = edge[1];
        edge[1] = temp;
    }
}

/**
 * Transforms an edgepair into its canonical form (i.e. lexiconographically smallest):
 * a b c d -> a < b and c < d and a < c.
 */
void transform_edgepair_into_canonical_form(EDGEPAIR edgepair) {
    int temp, temp2;
    if(edgepair[0] > edgepair[1]) {
        temp = edgepair[1];
        edgepair[1] = edgepair[0];
        edgepair[0] = temp;
    }
    if(edgepair[2] > edgepair[3]) {
        temp = edgepair[3];
        edgepair[3] = edgepair[2];
        edgepair[2] = temp;
    }
    if(edgepair[0] > edgepair[2]) {
        temp = edgepair[2];
        temp2 = edgepair[3];
        edgepair[2] = edgepair[0];
        edgepair[3] = edgepair[1];
        edgepair[0] = temp;
        edgepair[1] = temp2;
    }
}

/**
 * Transforms a triangle into its canonical form (i.e. lexiconographically smallest).
 *
 * Important: it is assumed that triangle[1] < triangle[2]
 */
void transform_triangle_into_canonical_form(TRIANGLE triangle) {
    //Using variant of bubblesort
    int temp;
    if(triangle[0] > triangle[1]) {
        temp = triangle[1];
        triangle[1] = triangle[0];
        triangle[0] = temp;
    }
    if(triangle[1] > triangle[2]) {
        temp = triangle[2];
        triangle[2] = triangle[1];
        triangle[1] = temp;
    }
}

/**
 * Transforms a diamond into its canonical form (i.e. lexiconographically smallest).
 * No assuptions are made about the triangle.
 */
void transform_triangle_into_canonical_form_full(TRIANGLE triangle) {
    //Using variant of bubblesort
    int temp;
    if(triangle[0] > triangle[1]) {
        temp = triangle[1];
        triangle[1] = triangle[0];
        triangle[0] = temp;
    }
    if(triangle[1] > triangle[2]) {
        temp = triangle[2];
        triangle[2] = triangle[1];
        triangle[1] = temp;
    }
    if(triangle[0] > triangle[1]) {
        temp = triangle[1];
        triangle[1] = triangle[0];
        triangle[0] = temp;
    }
}

/**
 * Transforms the vertexset into canonical form (using bubblesort).
 */
//Could use mergesort or quicksort, but won't make a big difference since the arrays are very small
void transform_vertexset_into_canonical_form(unsigned char vertexset[], int size) {
    int i, j, temp;
    for(j = 0; j < size; j++) {
        for(i = 1; i < size - j; i++) {
            if(vertexset[i - 1] > vertexset[i]) {
                temp = vertexset[i];
                vertexset[i] = vertexset[i - 1];
                vertexset[i - 1] = temp;
            }
        }
    }
}

/**
 * Transforms the vertices in vertexset into triangles (in one operation).
 * This operation must destroy all reducible triangles, so the resulting is always
 * canonical.
 */
void transform_vertexset_into_triangles(unsigned char vertexset[], int vertexset_size) {
    number_of_reducible_triangles = 0; //All existing reducible triangles will be destroyed

    int i, j, neighbour;
    if(current_number_of_vertices + (2 * vertexset_size) < number_of_vertices) {

        for(j = 0; j < vertexset_size; j++) {
            current_graph[current_number_of_vertices][0] = vertexset[j];
            current_graph[current_number_of_vertices][1] = current_number_of_vertices + 1;

            current_graph[current_number_of_vertices + 1][0] = vertexset[j];
            current_graph[current_number_of_vertices + 1][1] = current_number_of_vertices;

            //New labels
            edge_labels[vertexset[j]][current_number_of_vertices] = current_number_of_edges;
            edge_labels[current_number_of_vertices][vertexset[j]] = current_number_of_edges++;

            edge_labels[vertexset[j]][current_number_of_vertices + 1] = current_number_of_edges;
            edge_labels[current_number_of_vertices + 1][vertexset[j]] = current_number_of_edges++;

            edge_labels[current_number_of_vertices][current_number_of_vertices + 1] = current_number_of_edges;
            edge_labels[current_number_of_vertices + 1][current_number_of_vertices] = current_number_of_edges++;

            vertex_neighbourhood[current_number_of_vertices] = BIT(vertexset[j]) | BIT(current_graph[vertexset[j]][0]) | BIT(current_number_of_vertices + 1);
            vertex_neighbourhood[current_number_of_vertices + 1] = BIT(vertexset[j]) | BIT(current_graph[vertexset[j]][1]) | BIT(current_number_of_vertices);
            vertex_neighbourhood[vertexset[j]] = BIT(current_graph[vertexset[j]][2]) | BIT(current_number_of_vertices) | BIT(current_number_of_vertices + 1);

            vertex_neighbourhood[current_graph[vertexset[j]][0]] &= ~BIT(vertexset[j]);
            vertex_neighbourhood[current_graph[vertexset[j]][0]] |= BIT(current_number_of_vertices);

            vertex_neighbourhood[current_graph[vertexset[j]][1]] &= ~BIT(vertexset[j]);
            vertex_neighbourhood[current_graph[vertexset[j]][1]] |= BIT(current_number_of_vertices + 1);

            //The third neighbour remains the same
            for(i = 0; i < 2; i++) {
                neighbour = current_graph[vertexset[j]][i];

                current_graph[current_number_of_vertices + i][2] = neighbour;
                replace_neighbour(neighbour, vertexset[j], current_number_of_vertices + i);

                edge_labels[neighbour][current_number_of_vertices + i] = edge_labels[vertexset[j]][neighbour];
                edge_labels[current_number_of_vertices + i][neighbour] = edge_labels[vertexset[j]][neighbour];

                update_bridges_triangle_insert(vertexset[j], current_number_of_vertices + i, neighbour);
            }

            current_graph[vertexset[j]][0] = current_number_of_vertices;
            current_graph[vertexset[j]][1] = current_number_of_vertices + 1;

            //Add triangle to list of reducible triangles
            reducible_triangles[number_of_reducible_triangles][0] = vertexset[j];
            reducible_triangles[number_of_reducible_triangles][1] = current_number_of_vertices;
            reducible_triangles[number_of_reducible_triangles][2] = current_number_of_vertices + 1;
            number_of_reducible_triangles++;

            /*
             * Check if the vertex wasn't part of a diamond.
             * If it was, that irreducible triangle is no longer irreducible.
             */
            for(i = 0; i < number_of_irreducible_triangles; i++) {
                if(irreducible_triangles[i][0] == vertexset[j]) {
                    if(!other_element_of_vertexset_is_part_of_diamond(irreducible_triangles[i], vertexset, vertexset_size, 0)) {
                        //Diamond is transformed into an extended diamond
                        reducible_triangles[number_of_reducible_triangles][0] = irreducible_triangles[i][3];
                        reducible_triangles[number_of_reducible_triangles][1] = irreducible_triangles[i][1];
                        reducible_triangles[number_of_reducible_triangles][2] = irreducible_triangles[i][2];
                        transform_triangle_into_canonical_form(reducible_triangles[number_of_reducible_triangles]);
                        number_of_reducible_triangles++;
                    }
                    remove_irreducible_triangle(i);
                    break;
                } else if(irreducible_triangles[i][3] == vertexset[j]) {
                    if(!other_element_of_vertexset_is_part_of_diamond(irreducible_triangles[i], vertexset, vertexset_size, 3)) {
                        //Diamond is transformed into an extended diamond
                        reducible_triangles[number_of_reducible_triangles][0] = irreducible_triangles[i][0];
                        reducible_triangles[number_of_reducible_triangles][1] = irreducible_triangles[i][1];
                        reducible_triangles[number_of_reducible_triangles][2] = irreducible_triangles[i][2];
                        transform_triangle_into_canonical_form(reducible_triangles[number_of_reducible_triangles]);
                        number_of_reducible_triangles++;
                    }
                    remove_irreducible_triangle(i);
                    break;
                } else if(irreducible_triangles[i][1] == vertexset[j]
                        || irreducible_triangles[i][2] == vertexset[j]) {
                    //Somehow irreducible_triangles[i][2] == vertexset[j] never happens

                    //Diamond is transformed into two squares
                    remove_irreducible_triangle(i);
                    break;
                }
            }
            current_number_of_vertices += 2;
            //DEBUGASSERT(current_number_of_edges == 3 * current_number_of_vertices / 2);
        }
    } else { //All this extra stuff (such as updating bridges or nauty_graph) doesn't have to be done for the last step
        DEBUGASSERT(girth == 3);

        for(j = 0; j < vertexset_size; j++) {
            current_graph[current_number_of_vertices][0] = vertexset[j];
            current_graph[current_number_of_vertices][1] = current_number_of_vertices + 1;

            current_graph[current_number_of_vertices + 1][0] = vertexset[j];
            current_graph[current_number_of_vertices + 1][1] = current_number_of_vertices;

            //The third neighbour remains the same
            for(i = 0; i < 2; i++) {
                neighbour = current_graph[vertexset[j]][i];

                current_graph[current_number_of_vertices + i][2] = neighbour;
                replace_neighbour(neighbour, vertexset[j], current_number_of_vertices + i);
            }
            current_graph[vertexset[j]][0] = current_number_of_vertices;
            current_graph[vertexset[j]][1] = current_number_of_vertices + 1;

            current_number_of_vertices += 2;
        }
    }
}

/**
 * Returns 1 if an other element of the vertexset is part of the diamond, else returns 0.
 * Index can be 0 or 3.
 */
int other_element_of_vertexset_is_part_of_diamond(IRRED_TRIANGLE diamond, unsigned char vertexset[], int vertexset_size, int index) {
    int i;
    for(i = 0; i < vertexset_size;i++)
        if(vertexset[i] == diamond[(index + 1) % 4] || vertexset[i] == diamond[(index + 2) % 4] || vertexset[i] == diamond[(index + 3) % 4])
            return 1;
    return 0;
}

/**
 * Removes the irreducible triangle at index index in the list of irreducible triangles..
 */
void remove_irreducible_triangle(int index) {
    number_of_irreducible_triangles--;
    if(number_of_irreducible_triangles == 0)
        irreducible_triangles_bitvector = 0;
    else {
        irreducible_triangles_bitvector &= ~(BIT(irreducible_triangles[index][0]) | BIT(irreducible_triangles[index][1]) | BIT(irreducible_triangles[index][2]) | BIT(irreducible_triangles[index][3]));
    }
    if(number_of_irreducible_triangles > 0 && index < number_of_irreducible_triangles) {
        //TODO: could just update the pointers, but this doesn't cost much
        int i;
        for(i = 0; i < 4; i++)
            irreducible_triangles[index][i] = irreducible_triangles[number_of_irreducible_triangles][i];
    }
}

/**
 * Undoes the vertexset into triangles operation.
 */
void undo_vertexset_triangle_transformation(unsigned char vertexset[], int vertexset_size) {
    if(current_number_of_vertices < number_of_vertices)
        current_number_of_edges -= (3 * vertexset_size);

    int restore_neighbourhoods = (current_number_of_vertices < number_of_vertices);

    int i, j;
    int new_neighbour, org_neighbour;
    //Has to be undone in reverse order (i.e. LIFO)
    for(j = vertexset_size - 1; j >= 0; j--) {
        for(i = 0; i < 2; i++) { //The modified neighbours are at index 0 and 1
            new_neighbour = current_graph[vertexset[j]][i];

            //The neighbour at index 2 is the original neighbour
            org_neighbour = current_graph[new_neighbour][2];
            current_graph[vertexset[j]][i] = org_neighbour;

            replace_neighbour(org_neighbour, new_neighbour, vertexset[j]);
        }
        current_number_of_vertices -= 2;

        if(restore_neighbourhoods) {
            vertex_neighbourhood[vertexset[j]] = BIT(current_graph[vertexset[j]][0]) | BIT(current_graph[vertexset[j]][1]) | BIT(current_graph[vertexset[j]][2]);


            vertex_neighbourhood[current_graph[vertexset[j]][0]] &= ~BIT(current_number_of_vertices);
            vertex_neighbourhood[current_graph[vertexset[j]][0]] |= BIT(vertexset[j]);

            vertex_neighbourhood[current_graph[vertexset[j]][1]] &= ~BIT(current_number_of_vertices + 1);
            vertex_neighbourhood[current_graph[vertexset[j]][1]] |= BIT(vertexset[j]);
        }
    }

    /* Important: no need to restore the colours of the edge_labels of the orginial edges (for the same reason as with add_edge) */
    //current_number_of_vertices -= (2 * vertexset_size);
    //DEBUGASSERT(current_number_of_edges == 3 * current_number_of_vertices / 2);
}

/**
 * Adds an edge between 2 edges.
 * It is assumed that the edges are not adjacent.
 */
void add_edge(EDGEPAIR edge_pair) {
    int i;

/*
    for(i = 0; i < 2; i++) {
        //adj[current_number_of_vertices + i] = 0; //Not needed, because adj is only used for init
        current_graph[current_number_of_vertices + i][0] = edge_pair[2 * i];
        current_graph[current_number_of_vertices + i][1] = edge_pair[2 * i + 1];
        current_graph[current_number_of_vertices + i][2] = current_number_of_vertices + (1 + i) % 2;

        //Recycle old label
        edge_labels[edge_pair[2 * i]][current_number_of_vertices + i] = edge_labels[edge_pair[2 * i]][edge_pair[2 * i + 1]];
        edge_labels[current_number_of_vertices + i][edge_pair[2 * i]] = edge_labels[edge_pair[2 * i]][edge_pair[2 * i + 1]];

        //New label
        edge_labels[edge_pair[2 * i + 1]][current_number_of_vertices + i] = current_number_of_edges;
        edge_labels[current_number_of_vertices + i][edge_pair[2 * i + 1]] = current_number_of_edges++;


        replace_neighbour(edge_pair[2 * i], edge_pair[2 * i + 1], current_number_of_vertices + i);
        replace_neighbour(edge_pair[2 * i + 1], edge_pair[2 * i], current_number_of_vertices + i);
    }
*/

    current_graph[current_number_of_vertices][0] = edge_pair[0];
    current_graph[current_number_of_vertices][1] = edge_pair[1];
    current_graph[current_number_of_vertices][2] = current_number_of_vertices + 1;

    //Recycle old label
    edge_labels[edge_pair[0]][current_number_of_vertices] = edge_labels[edge_pair[0]][edge_pair[1]];
    edge_labels[current_number_of_vertices][edge_pair[0]] = edge_labels[edge_pair[0]][edge_pair[1]];

    //New label
    edge_labels[edge_pair[1]][current_number_of_vertices] = current_number_of_edges;
    edge_labels[current_number_of_vertices][edge_pair[1]] = current_number_of_edges++;

    replace_neighbour(edge_pair[0], edge_pair[1], current_number_of_vertices);
    replace_neighbour(edge_pair[1], edge_pair[0], current_number_of_vertices);


    current_graph[current_number_of_vertices + 1][0] = edge_pair[2];
    current_graph[current_number_of_vertices + 1][1] = edge_pair[3];
    current_graph[current_number_of_vertices + 1][2] = current_number_of_vertices;

    //Recycle old label
    edge_labels[edge_pair[2]][current_number_of_vertices + 1] = edge_labels[edge_pair[2]][edge_pair[3]];
    edge_labels[current_number_of_vertices + 1][edge_pair[2]] = edge_labels[edge_pair[2]][edge_pair[3]];

    //New label
    edge_labels[edge_pair[3]][current_number_of_vertices + 1] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 1][edge_pair[3]] = current_number_of_edges++;

    replace_neighbour(edge_pair[2], edge_pair[3], current_number_of_vertices + 1);
    replace_neighbour(edge_pair[3], edge_pair[2], current_number_of_vertices + 1);


    vertex_neighbourhood[current_number_of_vertices] = BIT(edge_pair[0]) | BIT(edge_pair[1]) | BIT(current_number_of_vertices + 1);
    vertex_neighbourhood[current_number_of_vertices + 1] = BIT(edge_pair[2]) | BIT(edge_pair[3]) | BIT(current_number_of_vertices);

    vertex_neighbourhood[edge_pair[0]] &= ~BIT(edge_pair[1]);
    vertex_neighbourhood[edge_pair[0]] |= BIT(current_number_of_vertices);

    vertex_neighbourhood[edge_pair[1]] &= ~BIT(edge_pair[0]);
    vertex_neighbourhood[edge_pair[1]] |= BIT(current_number_of_vertices);

    vertex_neighbourhood[edge_pair[2]] &= ~BIT(edge_pair[3]);
    vertex_neighbourhood[edge_pair[2]] |= BIT(current_number_of_vertices + 1);

    vertex_neighbourhood[edge_pair[3]] &= ~BIT(edge_pair[2]);
    vertex_neighbourhood[edge_pair[3]] |= BIT(current_number_of_vertices + 1);


    /* In total 5 edges need to be relabelled: 2 recycle an old label and 3 get a new label */
    edge_labels[current_number_of_vertices][current_number_of_vertices + 1] = current_number_of_edges;
    edge_labels[current_number_of_vertices + 1][current_number_of_vertices] = current_number_of_edges++;


    //Could do this more efficiently, but this isn't the bottleneck since there are usually no bridges
    /* Update bridges */
    int res = contains_bridge(edge_pair);
    if(res == 0) {
        //update_bridges_add_edge();
    } else if(res == 1 || res == 2) {
        res--;
        if(is_a_bridge(edge_pair[2*res], current_number_of_vertices + res))
            replace_bridge(edge_pair[2*res], edge_pair[2*res + 1], edge_pair[2*res], current_number_of_vertices + res);
        else
            replace_bridge(edge_pair[2*res], edge_pair[2*res + 1], edge_pair[2*res + 1], current_number_of_vertices + res);
    } else {
        //Either (edge_pair[0], new) or (edge_pair[1], new) is now a bridge
        if(is_a_bridge(edge_pair[0], current_number_of_vertices))
            replace_bridge(edge_pair[0], edge_pair[1], edge_pair[0], current_number_of_vertices);
        else
            replace_bridge(edge_pair[0], edge_pair[1], edge_pair[1], current_number_of_vertices);

        //Either (edge_pair[2], new) or (edge_pair[3], new) is now a bridge
        if(is_a_bridge(edge_pair[2], current_number_of_vertices + 1))
            replace_bridge(edge_pair[2], edge_pair[3], edge_pair[2], current_number_of_vertices + 1);
        else
            replace_bridge(edge_pair[2], edge_pair[3], edge_pair[3], current_number_of_vertices + 1);
    }

    //Even when one edge of the edgepair is a bridge, an other bridge can still be undone by inserting the new edge
    update_bridges_add_edge();

    if(girth == 3 || current_number_of_vertices < number_of_vertices - 2) {
        /**
         * Can only remove at most 1 diamond at a time.
         * Only possible edgepairs are 01 23 and 02 13.
         *
         * Since the element at index 0 is always < the elements at index 1, 2 or 3
         * it should only be checked if edgepair[0, 1] = 01 or 02
         */
        for(i = 0; i < number_of_irreducible_triangles; i++) {
            if(edge_pair[0] == irreducible_triangles[i][0] && (edge_pair[1] == irreducible_triangles[i][1] || edge_pair[1] == irreducible_triangles[i][2])) {
                remove_irreducible_triangle(i);
                break;
            }
        }
    } else { //If girth4 and number_of_vertices - 2, then all diamonds must be removed
        number_of_irreducible_triangles = 0;
        irreducible_triangles_bitvector = 0;
    }
    //All reducible triangles (at most 2) must be destroyed
    number_of_reducible_triangles = 0;

    current_number_of_vertices += 2;
}

/**
 * Removes the edge (and its vertices) which was added between edge_pair by the add_edge()
 */
void remove_edge(EDGEPAIR edge_pair) {
/*
    int i, j;

    for(i = 0; i < 2; i++) {
        for(j = 0; j < 2; j++) {
            replace_neighbour(edge_pair[2 * i + j], current_number_of_vertices - 2 + i, edge_pair[2 * i + (1 + j) % 2]);
        }
    }
*/
    replace_neighbour(edge_pair[0], current_number_of_vertices - 2, edge_pair[1]);
    replace_neighbour(edge_pair[1], current_number_of_vertices - 2, edge_pair[0]);

    replace_neighbour(edge_pair[2], current_number_of_vertices - 1, edge_pair[3]);
    replace_neighbour(edge_pair[3], current_number_of_vertices - 1, edge_pair[2]);


    /**
     * Important!: the edge_labels of the original edges do not need to be restored, because
     * they still have their original value. This is because there will never be a
     * new edge added which is identical to the orginal edge.
     */

    current_number_of_vertices -= 2;
    current_number_of_edges -= 3;

    vertex_neighbourhood[edge_pair[0]] &= ~BIT(current_number_of_vertices);
    vertex_neighbourhood[edge_pair[0]] |= BIT(edge_pair[1]);

    vertex_neighbourhood[edge_pair[1]] &= ~BIT(current_number_of_vertices);
    vertex_neighbourhood[edge_pair[1]] |= BIT(edge_pair[0]);

    vertex_neighbourhood[edge_pair[2]] &= ~BIT(current_number_of_vertices + 1);
    vertex_neighbourhood[edge_pair[2]] |= BIT(edge_pair[3]);

    vertex_neighbourhood[edge_pair[3]] &= ~BIT(current_number_of_vertices + 1);
    vertex_neighbourhood[edge_pair[3]] |= BIT(edge_pair[2]);

    //DEBUGASSERT(current_number_of_edges == 3 * current_number_of_vertices / 2);
}

/**
 * Returns 0 if the edgepair doesn't consist out of any bridges.
 * Returns 1 if the first edge of the edgepair is a bridge,
 * returns 2 if the second edge of the edgepair is a bridge and
 * returns 3 if both edges of the edgepair are bridges.
 */
int contains_bridge(EDGEPAIR edge_pair) {
    int res = 0;
    int i;
    for(i = 0; i < number_of_bridges;i++) {
        if(bridges[i][0] == edge_pair[0] && bridges[i][1] == edge_pair[1])
            res++;
        else if(bridges[i][0] == edge_pair[2] && bridges[i][1] == edge_pair[3])
            res += 2;
    }
    if(res > 3) {
        fprintf(stderr, "Error: list of bridges contains duplicates\n");
        exit(1);
    }
    return res;
}

/**
 * Returns 1 if the edge between this edgepair will be part of a square, else
 * returns 0.
 *
 * Nothing is assumed about the edgepair (it doesn't have to be in canonical form).
 */
int inserted_edge_will_be_part_of_square(EDGEPAIR edgepair) {
    //Old way
/*
    //Will be part of a square if at least one of the vertices edgepair[2] or edgepair[3] is a neighbour of edgepair[0] or edgepair[1]
    int i, j;

    setword neighbours = 0;
    for(i = 0; i < 2; i++) {
        for(j = 0; j < REG; j++) {
            neighbours |= BIT(current_graph[edgepair[i]][j]);
        }
    }
    //setword other_bitvector = BIT(edgepair[2]) | BIT(edgepair[3]);
    return (neighbours & (BIT(edgepair[2]) | BIT(edgepair[3]))) > 0;
*/
    return ((vertex_neighbourhood[edgepair[0]] | vertex_neighbourhood[edgepair[1]]) & (BIT(edgepair[2]) | BIT(edgepair[3]))) > 0;
}

/**
 * Returns 1 if the inserted edge will be part of at least 2 pentagons but won't
 * be part of a square, else returns 0.
 */
int inserted_edge_will_be_part_of_at_least_2_pentagons_but_not_of_square(EDGEPAIR edgepair) {
    //Will be part of a square if at least one of the vertices edgepair[2] or edgepair[3] is a neighbour of edgepair[0] or edgepair[1]
    setword neighbours0 = vertex_neighbourhood[edgepair[0]] | vertex_neighbourhood[edgepair[1]];

    if((neighbours0 & (BIT(edgepair[2]) | BIT(edgepair[3]))) > 0)
        return 0; //inserted edge will be part of square

    setword neighbours1 = vertex_neighbourhood[edgepair[2]] | vertex_neighbourhood[edgepair[3]];
    //return (neighbours0 & neighbours1) > 0;

    return POPCOUNT(neighbours0 & neighbours1) > 1;
}

/**
 * Generates the edgepairs in case there is 1 triangle.
 */
void generate_edgepairs_one_triangle(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size) {
    DEBUGASSERT(number_of_reducible_triangles == 1);

    int i, j;
    unsigned char third_vertex;

    //Neighbours which aren't part of the triangle
    TRIANGLE external_neighbours;
    for(i = 0; i < 3; i++) {
        edge_pairs_list[*edge_pair_list_size][0] = reducible_triangles[0][i];
        edge_pairs_list[*edge_pair_list_size][1] = reducible_triangles[0][(i + 1) % 3];

        third_vertex = reducible_triangles[0][(i + 2) % 3];
        edge_pairs_list[*edge_pair_list_size][2] = third_vertex;
        for(j = 0; j < REG; j++) {
            if(current_graph[third_vertex][j] != reducible_triangles[0][i] && current_graph[third_vertex][j] != reducible_triangles[0][(i + 1) % 3])
            //if(!is_part_of_same_reducible_triangle(current_graph[third_vertex][j], 0))
                break;
        }
        edge_pairs_list[*edge_pair_list_size][3] = current_graph[third_vertex][j];
        external_neighbours[(i + 2) % 3] = current_graph[third_vertex][j];

        transform_edgepair_into_canonical_form(edge_pairs_list[*edge_pair_list_size]);

        int index0 = edge_labels[edge_pairs_list[*edge_pair_list_size][0]][edge_pairs_list[*edge_pair_list_size][1]];
        int index1 = edge_labels[edge_pairs_list[*edge_pair_list_size][2]][edge_pairs_list[*edge_pair_list_size][3]];

        edgepair_index[index0][index1] = *edge_pair_list_size;

        (*edge_pair_list_size)++;
    }
    for(i = 0; i < 3; i++) {
        //Edgepair is in a square
        if(is_neighbour(external_neighbours[i], external_neighbours[(i + 1) % 3])) {
            edge_pairs_list[*edge_pair_list_size][0] = reducible_triangles[0][i];
            edge_pairs_list[*edge_pair_list_size][1] = reducible_triangles[0][(i + 1) % 3];
            edge_pairs_list[*edge_pair_list_size][2] = external_neighbours[i];
            edge_pairs_list[*edge_pair_list_size][3] = external_neighbours[(i + 1) % 3];

            transform_edgepair_into_canonical_form(edge_pairs_list[*edge_pair_list_size]);

            int index0 = edge_labels[edge_pairs_list[*edge_pair_list_size][0]][edge_pairs_list[*edge_pair_list_size][1]];
            int index1 = edge_labels[edge_pairs_list[*edge_pair_list_size][2]][edge_pairs_list[*edge_pair_list_size][3]];

            edgepair_index[index0][index1] = *edge_pair_list_size;

            (*edge_pair_list_size)++;
        }
    }

}


/**
 * Generates the edgepairs in case there is one diamond and no reducible triangles.
 * The edgepairs consist of the common edge and any other edge which isnt adjacent to it.
 */
//Only one edgepair is generated, all other edgepairs will yield an edge that isn't minimal
void generate_edgepairs_triangle_free_one_diamond(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size) {
    DEBUGASSERT(number_of_reducible_triangles == 0 && number_of_irreducible_triangles == 1);
    /**
     * Adding one edgepair is sufficient, because 1 and 2 will always be in the same (vertex)orbit.
     * So edpairs 01 23 and 02 13 will always be in the same orbit.
     *
     * If the edgelist only contains 1 element, determine_edge_orbits won't
     * apply the generators to the edgepair.
     */
    //Edgepair 01 23
    edge_pairs_list[*edge_pair_list_size][0] = irreducible_triangles[0][0];
    edge_pairs_list[*edge_pair_list_size][1] = irreducible_triangles[0][1];
    edge_pairs_list[*edge_pair_list_size][2] = irreducible_triangles[0][2];
    edge_pairs_list[*edge_pair_list_size][3] = irreducible_triangles[0][3];
    transform_edgepair_into_canonical_form(edge_pairs_list[*edge_pair_list_size]);

    //edgepair_index not needed because *edge_pair_list_size == 1, so no need to determine edge orbits
    (*edge_pair_list_size)++;

    DEBUGASSERT(*edge_pair_list_size == 1);
}

/**
 * Returns 1 if the edgepair is part of a square or a pentagon.
 */
int edgepair_is_part_of_square_or_pentagon(EDGEPAIR edgepair) {
    //The neighbours (on distance 1) of the edge (excluding the 2 vertices of the edge)
    setword neighbours[2];
    int i;
    for(i = 0; i < 2; i++) {
        neighbours[i] = vertex_neighbourhood[edgepair[2 * i]] | vertex_neighbourhood[edgepair[2 * i + 1]];
        neighbours[i] &= ~(BIT(edgepair[2 * i]) | BIT(edgepair[2 * i + 1]));
    }
    setword edgepair_bitvector = BIT(edgepair[0]) | BIT(edgepair[1]) | BIT(edgepair[2]) | BIT(edgepair[3]);
    int number_common = POPCOUNT((neighbours[0] | neighbours[1]) & edgepair_bitvector);

    //DEBUGASSERT(number_common % 2 == 0);

    if(number_common == 4) {
        //Is part of square
        return 1;
    } else if(number_common == 2) {
        //Might be part of pentagon

        return (neighbours[0] & neighbours[1]) > 0;
    } else
        return 0;
}

/**
 * Generates edgepairs in case there are 2 reducible triangles.
 */
void generate_edgepairs_two_triangles(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size) {
    DEBUGASSERT(number_of_reducible_triangles == 2);

    int i, j;
    EDGEPAIR edgepair;
    for(i = 0; i < REG; i++) {
        edgepair[0] = reducible_triangles[0][i];
        edgepair[1] = reducible_triangles[0][(i + 1) % 3];

        for(j = 0; j < REG; j++) {
            edgepair[2] = reducible_triangles[1][j];
            edgepair[3] = reducible_triangles[1][(j + 1) % 3];
            //If edgepair is not part of square or pentagon, there will be an edge
            //with colour 10 which is part of 2 squares
            if(edgepair_is_part_of_square_or_pentagon(edgepair)) {
                edge_pairs_list[*edge_pair_list_size][0] = edgepair[0];
                edge_pairs_list[*edge_pair_list_size][1] = edgepair[1];
                edge_pairs_list[*edge_pair_list_size][2] = edgepair[2];
                edge_pairs_list[*edge_pair_list_size][3] = edgepair[3];

                transform_edgepair_into_canonical_form(edge_pairs_list[*edge_pair_list_size]);

                edgepair_index[edge_labels[edge_pairs_list[*edge_pair_list_size][0]][edge_pairs_list[*edge_pair_list_size][1]]][edge_labels[edge_pairs_list[*edge_pair_list_size][2]][edge_pairs_list[*edge_pair_list_size][3]]] = *edge_pair_list_size;

                (*edge_pair_list_size)++;
            }
        }
    }
    DEBUGASSERT(*edge_pair_list_size <= 9);
}

/**
 * Simple filter for squares.
 * Returns 1 if current_graphs contains squares, else returns 0.
 */
//Is faster than using the method from find_squares
int contains_squares() {
    EDGE edge;
    unsigned char i, j;
    unsigned char neighbour;
    for(i = 0; i < current_number_of_vertices - 1; i++) {
        for(j = 0; j < REG; j++) {
            if(i < current_graph[i][j]) {
                neighbour = current_graph[i][j];
                edge[0] = i;
                edge[1] = neighbour;
                if(edge_is_part_of_square(edge)) {
                    return 1;
                }
            }
        }
    }

    return 0;
}

int square_was_already_adjacent_to_other_square(int square_index, EDGE adjacent_squares[], int adjacent_squares_sizes, int *index_adj_square) {
    int i;
    for(i = 0; i < adjacent_squares_sizes; i++) {
        if(adjacent_squares[i][0] == square_index) {
            *index_adj_square = adjacent_squares[i][1];
            return 1;
        } else if(adjacent_squares[i][1] == square_index) {
            *index_adj_square = adjacent_squares[i][0];
            return 1;
        }
/*
        if(adjacent_squares[i][0] == square_index || adjacent_squares[i][1] == square_index) {
            return 1;
        }
*/
    }
    return 0;
}

/**
 * Updates the list of squares with a new square v1, v2, v3, v4.
 * Returns 1 if the squares can still be destroyed by adding an edge, else returns 0.
 */
int update_adjacent_squares(SQUARE squares[], setword squares_bitvectors[], int *squares_size, EDGE adjacent_squares[], int *adjacent_squares_sizes, unsigned char v1, unsigned char v2, unsigned char v3, unsigned char v4) {
    squares_bitvectors[*squares_size] = BIT(v1) | BIT(v2) | BIT(v3) | BIT(v4);
    int i, common_vertices, adj_square_index;
    for(i = 0; i < *squares_size; i++) {
        common_vertices = POPCOUNT((squares_bitvectors[i] & squares_bitvectors[*squares_size]));
        DEBUGASSERT(common_vertices == 0 || common_vertices == 2 || common_vertices == 3);
        if(common_vertices == 2) {
            if(!square_was_already_adjacent_to_other_square(i, adjacent_squares, *adjacent_squares_sizes, &adj_square_index)) {
                adjacent_squares[*adjacent_squares_sizes][0] = i;
                adjacent_squares[*adjacent_squares_sizes][1] = *squares_size;
                (*adjacent_squares_sizes)++;
                break;
            } else {
                //Special case 3-chain which cannot be made square free by adding a single edge
                //(but this case only rarely happens), so is not really needed to filter here
                if((squares_bitvectors[*squares_size] & squares_bitvectors[adj_square_index]) > 0)
                    return 0;
            }
        } else if(common_vertices == 3) {
            //In this case the squares can't be destroyed by adding a single edge
            return 0;
        }
    }
    squares[*squares_size][0] = v1;
    squares[*squares_size][1] = v2;
    squares[*squares_size][2] = v3;
    squares[*squares_size][3] = v4;
    (*squares_size)++;

    return 1;
}

/**
 * Finds all squares and determines which squares are adjacent.
 * The squares are in "canonical form": if the square is "a b c d",
 * a = min(a, b, c, d) and b < d and both b and d are neighbours of a.
 *
 * Returns 1 if the squares can still be destroyed by adding an edge, else returns 0.
 *
 * Remark: it is assumed that current_graph contains no diamonds.
 */
int find_squares(SQUARE squares[], setword squares_bitvectors[], int *squares_size, EDGE adjacent_squares[], int *adjacent_squares_sizes) {
    int i, j, k, l;
    unsigned char neighbour, neighbour2, neighbour3;
    *squares_size = 0;
    *adjacent_squares_sizes = 0;
    for(i = 0; i < current_number_of_vertices - 1; i++) {
        for(j = 0; j < REG; j++) {
            neighbour = current_graph[i][j];
            if(i < neighbour) {
               for(k = 0; k < REG; k++) {
                   neighbour2 = current_graph[neighbour][k];
                   if(neighbour2 > i) {
                       for(l = 0; l < REG; l++) {
                           neighbour3 = current_graph[neighbour2][l];
                            if(neighbour3 > neighbour) { //Which implies neighbour3 > i
                                if(is_neighbour(neighbour3, i)) {
                                    /**
                                     * If there are more than 2 squares, it's never possible to
                                     * find an edgepair which will destroy all squares AND
                                     * which will have the minimal colour (except in the case
                                     * where there are 8 vertices and 4 squares).
                                     */
                                    if(*squares_size == 2 && current_number_of_vertices > 8)
                                        return 0;
                                    else {
                                        //No special checks needed in case of current_number_of_vertices == 8, update_adjacent_squares suffices
                                        if(!update_adjacent_squares(squares, squares_bitvectors, squares_size, adjacent_squares, adjacent_squares_sizes, i, neighbour, neighbour2, neighbour3)) {
                                            return 0;
                                        }
                                    }
                                }
                            }
                       }
                   }
               }
            }
        }
    }
    return 1;
}

/**
 * Same as find_squares(), but here current_graph can still contain diamonds.
 */
int find_squares_contains_diamonds(SQUARE squares[], setword squares_bitvectors[], int *squares_size, EDGE adjacent_squares[], int *adjacent_squares_sizes) {
    int i, j, k, l;
    unsigned char neighbour, neighbour2, neighbour3;
    *squares_size = 0;
    for(i = 0; i < current_number_of_vertices - 1; i++) {
        if((BIT(i) & irreducible_triangles_bitvector) == 0) {
            for(j = 0; j < REG; j++) {
                neighbour = current_graph[i][j];
                if(i < neighbour && (BIT(neighbour) & irreducible_triangles_bitvector) == 0) {
                    for(k = 0; k < REG; k++) {
                        neighbour2 = current_graph[neighbour][k];
                        if(neighbour2 > i && (BIT(neighbour2) & irreducible_triangles_bitvector) == 0) {
                            for(l = 0; l < REG; l++) {
                                neighbour3 = current_graph[neighbour2][l];
                                if(neighbour3 > neighbour && (BIT(neighbour3) & irreducible_triangles_bitvector) == 0) { //Which implies neighbour3 > i
                                    if(is_neighbour(neighbour3, i)) {
                                        /**
                                         * If there are more than 2 squares, it's never possible to
                                         * find an edgepair which will destroy all squares AND
                                         * which will have the minimal colour (except in the case
                                         * where there are 8 vertices and 4 squares).
                                         */
                                        if(*squares_size == 2 && current_number_of_vertices > 8)
                                            return 0;
                                        else {
                                            //No special checks needed in case of current_number_of_vertices == 8, update_adjacent_squares suffices
                                            if(!update_adjacent_squares(squares, squares_bitvectors, squares_size, adjacent_squares, adjacent_squares_sizes, i, neighbour, neighbour2, neighbour3)) {
                                                return 0;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return 1;
}

/**
 * Returns 1 if the 2 reducible triangles are adjacent, else returns 0.
 */
int are_adjacent_triangles() {
    DEBUGASSERT(number_of_reducible_triangles == 2);

    setword neighbours = vertex_neighbourhood[reducible_triangles[0][0]] | vertex_neighbourhood[reducible_triangles[0][1]] | vertex_neighbourhood[reducible_triangles[0][2]];
    setword other_triangle = BIT(reducible_triangles[1][0]) | BIT(reducible_triangles[1][1]) | BIT(reducible_triangles[1][2]);

    return (neighbours & other_triangle) > 0;
}

/**
 * Returns edgepairs which can remove all reducible triangles by adding an edge between
 * an edgepair.
 *
 * The edgepairs are in canonical form (i.e. lexiconographically smallest)
 */
void find_edge_pairs(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size) {
    *edge_pair_list_size = 0;
    edgepair_list_index_squares = 0;
    if(number_of_reducible_triangles == 2) {
        //Are always adjacent triangles if girth > 3
        if(girth > 3 || are_adjacent_triangles()) {
            generate_edgepairs_two_triangles(edge_pairs_list, edge_pair_list_size);
        }
    } else if(number_of_reducible_triangles == 1) {
        generate_edgepairs_one_triangle(edge_pairs_list, edge_pair_list_size);

        //Can be < 15 in case there is a square
        DEBUGASSERT(*edge_pair_list_size <= 15);

    } else if(number_of_reducible_triangles == 0) { /* number_of_reducible_triangles == 0 */
        generate_edgepairs_no_triangles(edge_pairs_list, edge_pair_list_size);
    }
}

/**
 * Fills the list of min_edge_bitvectors (for new_edge_has_min_colour).
 */
void update_min_edges() {
    DEBUGASSERT(number_of_reducible_triangles == 0);

    if(min_colour_one < MAX_EDGE_COLOUR_TWO) {
        int i;
        for(i = 0; i < min_edges_size; i++) {
            min_edge_bitvectors[i] = BIT(min_edges[i][0]) | BIT(min_edges[i][1]);
        }
    } else {
        min_edges_size = 0;
    }

}

/**
 * Generates the edgepairs which destroy all squares.
 * In this case it is assumed that current_graph contains 1 square.
 */
void generate_non_adjacent_edge_pairs_one_square(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size, SQUARE square) {
    //Could do this using BFS

    *edge_pair_list_size = 0;
    int i, j, k;
    int ep0, ep1, ep3;
    for(i = 0; i < 4; i++) {
        ep0 = square[i];
        ep1 = square[(i + 1) % 4];
        for(j = 0; j < current_number_of_vertices - 1; j++) {
            if(j != ep0 && j != ep1) {
                for(k = 0; k < REG; k++) {
                    ep3 = current_graph[j][k];
                    if(j < ep3 && ep3 != ep0 && ep3 != ep1) {
                        edge_pairs_list[*edge_pair_list_size][0] = ep0;
                        edge_pairs_list[*edge_pair_list_size][1] = ep1;
                        edge_pairs_list[*edge_pair_list_size][2] = j;
                        edge_pairs_list[*edge_pair_list_size][3] = ep3;

                        //Can't be part of square and must be part of pentagon, otherwise colour won't be minimal
                        //Don't perform new_edge_has_min_colour_no_squares here, is slower
                        if(inserted_edge_will_be_part_of_at_least_2_pentagons_but_not_of_square(edge_pairs_list[*edge_pair_list_size])) {
                            transform_edgepair_into_canonical_form(edge_pairs_list[*edge_pair_list_size]);
                            int index0 = edge_labels[edge_pairs_list[*edge_pair_list_size][0]][edge_pairs_list[*edge_pair_list_size][1]];
                            int index1 = edge_labels[edge_pairs_list[*edge_pair_list_size][2]][edge_pairs_list[*edge_pair_list_size][3]];

                            edgepair_index[index0][index1] = *edge_pair_list_size;

                            (*edge_pair_list_size)++;
                        }
                    }
                }
            }
        }
    }
}

/**
 * Generates the edgepairs which destroy all squares.
 * In this case it is assumed that current_graph contains 1 square and it
 * may also contain diamonds.
 */
void generate_non_adjacent_edge_pairs_one_square_contains_diamonds(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size, SQUARE square) {
    *edge_pair_list_size = 0;
    int i, j, k;
    int ep0, ep1, ep3;
    int diamond;
    for(i = 0; i < 4; i++) {
        ep0 = square[i];
        ep1 = square[(i + 1) % 4];
        for(j = 0; j < current_number_of_vertices - 1; j++) {
            if(j != ep0 && j != ep1) {
                for(k = 0; k < REG; k++) {
                    ep3 = current_graph[j][k];
                    if(j < ep3 && ep3 != ep0 && ep3 != ep1
                            && !(is_part_of_irreducible_triangle_diamond(j, &diamond) && is_part_of_same_irreducible_triangle(ep3, diamond))) {
                        edge_pairs_list[*edge_pair_list_size][0] = ep0;
                        edge_pairs_list[*edge_pair_list_size][1] = ep1;
                        edge_pairs_list[*edge_pair_list_size][2] = j;
                        edge_pairs_list[*edge_pair_list_size][3] = ep3;

                        //Can't be part of square and must be part of pentagon, otherwise colour won't be minimal
                        if(inserted_edge_will_be_part_of_at_least_2_pentagons_but_not_of_square(edge_pairs_list[*edge_pair_list_size])) {
                            transform_edgepair_into_canonical_form(edge_pairs_list[*edge_pair_list_size]);

                            int index0 = edge_labels[edge_pairs_list[*edge_pair_list_size][0]][edge_pairs_list[*edge_pair_list_size][1]];
                            int index1 = edge_labels[edge_pairs_list[*edge_pair_list_size][2]][edge_pairs_list[*edge_pair_list_size][3]];
                            edgepair_index[index0][index1] = *edge_pair_list_size;

                            (*edge_pair_list_size)++;
                        }
                    }
                }
            }
        }
    }
}

/**
 * Returns 1 if the 2 squares have multiple common neighbours, else returns 0.
 */
int squares_have_multiple_common_neighbours(SQUARE squares[2]) {
    setword neighbours[2];
    int i, j;
    for(i = 0; i < 2; i++) {
        neighbours[i] = (setword) 0;
        for(j = 0; j < 4; j++) {
            neighbours[i] |= vertex_neighbourhood[squares[i][j]];
        }
    }
    return POPCOUNT(neighbours[0] & neighbours[1]) >= 2;
}

/**
 * Generates the edgepairs which destroy all squares without generating any new squares.
 * In this case it is assumed that current_graph contains 2 squares which
 * are not adjacent.
 */
void generate_non_adjacent_edge_pairs_two_squares(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size, SQUARE square1, SQUARE square2) {
    *edge_pair_list_size = 0;
    int i, j;
    int ep0, ep1;
    for(i = 0; i < 4; i++) {
        ep0 = square1[i];
        ep1 = square1[(i + 1) % 4];
        for(j = 0; j < 4; j++) {
            edge_pairs_list[*edge_pair_list_size][0] = ep0;
            edge_pairs_list[*edge_pair_list_size][1] = ep1;
            edge_pairs_list[*edge_pair_list_size][2] = square2[j];
            edge_pairs_list[*edge_pair_list_size][3] = square2[(j + 1) % 4];

            //Can't be part of square and must be part of pentagon, otherwise colour won't be minimal
            if(inserted_edge_will_be_part_of_at_least_2_pentagons_but_not_of_square(edge_pairs_list[*edge_pair_list_size])) {
                transform_edgepair_into_canonical_form(edge_pairs_list[*edge_pair_list_size]);
                int index0 = edge_labels[edge_pairs_list[*edge_pair_list_size][0]][edge_pairs_list[*edge_pair_list_size][1]];
                int index1 = edge_labels[edge_pairs_list[*edge_pair_list_size][2]][edge_pairs_list[*edge_pair_list_size][3]];

                edgepair_index[index0][index1] = *edge_pair_list_size;

                (*edge_pair_list_size)++;
            }
        }
    }
}

/**
 * Generates the edgepairs which destroy all squares without generating any new squares.
 * In this case it is assumed that current_graph contains 2 squares which
 * are adjacent.
 */
void generate_non_adjacent_edge_pairs_two_adjacent_squares(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size, SQUARE squares[], EDGE common_edge) {
    *edge_pair_list_size = 0;
    //Connecting the two outer edges of the nonadjacent squares
    int i, j;
    for(i = 0; i < 2; i++) {
        for(j = 0; j < 4; j++)
            if(squares[i][j] != common_edge[0] && squares[i][j] != common_edge[1] &&
                    squares[i][(j + 1) % 4] != common_edge[0] && squares[i][(j + 1) % 4] != common_edge[1]) {
                edge_pairs_list[0][2 * i] = squares[i][j];
                edge_pairs_list[0][2 * i + 1] = squares[i][(j + 1) % 4];
                break;
            }
    }

    //Can't be part of square and must be part of pentagon, otherwise colour won't be minimal
    if(!inserted_edge_will_be_part_of_square(edge_pairs_list[*edge_pair_list_size])) {
        transform_edgepair_into_canonical_form(edge_pairs_list[*edge_pair_list_size]);
        int index0 = edge_labels[edge_pairs_list[*edge_pair_list_size][0]][edge_pairs_list[*edge_pair_list_size][1]];
        int index1 = edge_labels[edge_pairs_list[*edge_pair_list_size][2]][edge_pairs_list[*edge_pair_list_size][3]];

        edgepair_index[index0][index1] = *edge_pair_list_size;

        (*edge_pair_list_size)++;
    }

    //Other edgepairs won't have minimal colour

}

/**
 * Returns the next element in the bitvector at position >= pos.
 */
int nextelement_long(setword bitvector, setword pos) {
    setword bitvector_local;

    if (pos < 0)
        bitvector_local = bitvector;
    else
        bitvector_local = bitvector & BITMASK(pos);

    if (bitvector_local == 0)
        return -1;
    else
        return FIRSTBIT(bitvector_local);
}

#define MAXPOS (WORDSIZE - 1)

/**
 * bitvector is assumed to be a bitvector with 2 nonzero elements.
 * When this method is finished, common_edge will contain those 2 elements.
 * Common_edge will be in canonical form.
 */
void determine_common_edge(setword bitvector, EDGE common_edge) {
    DEBUGASSERT(POPCOUNT(bitvector) == 2);

/*
    int pos = -1;
    //int pos = MAXPOS - current_number_of_vertices -1;
    int i = 0;
    while((pos = nextelement_long(bitvector, pos)) >= 0) {
        common_edge[i] = MAXPOS - pos;
        i++;
    }
    DEBUGASSERT(i == 2);
*/

    //Slightly faster:
    int pos = MAXPOS - current_number_of_vertices - 1;
    pos = nextelement_long(bitvector, pos);
    common_edge[1] = MAXPOS - pos;
    common_edge[0] = MAXPOS - nextelement_long(bitvector, pos);

}

/**
 * Generates the edgepairs which destroy all squares without generating any new squares.
 * In this case it is assumed that current_graph contains a chain of 4 adjacent
 * squares. These can only be destroyed if current_number_of_vertices == 8.
 */
void generate_non_adjacent_edge_pairs_four_squares(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size, setword squares_bitvectors[], EDGE adjacent_squares[], int adjacent_squares_sizes) {
    //DEBUGASSERT(snarks || girth == 5);
    DEBUGASSERT(adjacent_squares_sizes == 2);

    //Remark: this method is not very efficient, but it is only called in 1 special case
    *edge_pair_list_size = 0;
    int i, j, k, l, next, next0, index0, index1;
    for(i = 0; i < current_number_of_vertices - 1; i++) {
        for(j = 0; j < REG; j++) {
            next = current_graph[i][j];
            if(i < next) {
                for(k = i + 1; k < current_number_of_vertices - 1; k++) {
                    if(k != next) { // k != i is automatically implied because k >= i+1
                        for(l = 0; l < REG; l++) {
                            next0 = current_graph[k][l];
                            if(k < next0 && next0 != next) {
                                //Add edgepair i next j next0
                                edge_pairs_list[*edge_pair_list_size][0] = i;
                                edge_pairs_list[*edge_pair_list_size][1] = next;
                                edge_pairs_list[*edge_pair_list_size][2] = k;
                                edge_pairs_list[*edge_pair_list_size][3] = next0;


                                if(inserted_edge_will_be_part_of_at_least_2_pentagons_but_not_of_square(edge_pairs_list[*edge_pair_list_size])) {
                                    transform_edgepair_into_canonical_form(edge_pairs_list[*edge_pair_list_size]);
                                    int index0 = edge_labels[edge_pairs_list[*edge_pair_list_size][0]][edge_pairs_list[*edge_pair_list_size][1]];
                                    int index1 = edge_labels[edge_pairs_list[*edge_pair_list_size][2]][edge_pairs_list[*edge_pair_list_size][3]];
                                    edgepair_index[index0][index1] = *edge_pair_list_size;

                                    (*edge_pair_list_size)++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //And both pairs of edges are in the same orbit
    DEBUGASSERT(*edge_pair_list_size == 2);
}

/**
 * Generates the edgepairs which destroy all squares without generating any new squares.
 * In this case it is assumed that current_graph contains no squares.
 *
 * Remark: it is assumed that current_graph contains no diamonds.
 */
void generate_non_adjacent_edge_pairs_square_free(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size, int check_remove) {
    update_min_edges();

    //DEBUGASSERT(snarks || girth == 5);

    *edge_pair_list_size = 0;
    int i, j, k, l, next, next0, index0, index1;
    for(i = 0; i < current_number_of_vertices - 1; i++) {
        for(j = 0; j < REG; j++) {
            next = current_graph[i][j];
            if(i < next) {
                for(k = i + 1; k < current_number_of_vertices - 1; k++) {
                    if(k != next) { // k != i is automatically implied because k >= i+1
                        for(l = 0; l < REG; l++) {
                            next0 = current_graph[k][l];
                            if(k < next0 && next0 != next) {
                                //Add edgepair i next j next0
                                edge_pairs_list[*edge_pair_list_size][0] = i;
                                edge_pairs_list[*edge_pair_list_size][1] = next;
                                edge_pairs_list[*edge_pair_list_size][2] = k;
                                edge_pairs_list[*edge_pair_list_size][3] = next0;

                                if(check_remove) {
                                    index0 = edge_labels[i][next];
                                    index1 = edge_labels[k][next0];
                                }

                                //Als inserted edge deel van square zal zijn, zal de graaf nog altijd kleurbaar zijn
                                //if((!check_remove || !ISREMOVED_EDGEPAIR(index0, index1)) && !inserted_edge_will_be_part_of_square(edge_pairs_list[*edge_pair_list_size])) {
                                if((!check_remove || ((edge_in_cycles[index0] & edge_in_cycles[index1]) == 0)) && !inserted_edge_will_be_part_of_square(edge_pairs_list[*edge_pair_list_size])) {
                                    if(new_edge_has_min_colour_no_squares(edge_pairs_list[*edge_pair_list_size])) {
                                        if(!check_remove) {
                                            index0 = edge_labels[i][next];
                                            index1 = edge_labels[k][next0];
                                        }

                                        edgepair_index[index0][index1] = *edge_pair_list_size;

                                        (*edge_pair_list_size)++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 * Generates the edgepairs which destroy all squares without generating any new squares.
 * In this case it is assumed that current_graph contains no squares.
 *
 * Remark: current_graph may contain diamonds.
 */
void generate_non_adjacent_edge_pairs_square_free_contains_diamonds(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size) {
    update_min_edges();

    *edge_pair_list_size = 0;
    int i, j, k, l, next, next0;
    int diamond;
    for(i = 0; i < current_number_of_vertices - 3; i++) {
        for(j = 0; j < REG; j++) {
            next = current_graph[i][j];
            if(i < next) {
                if(!(is_part_of_irreducible_triangle_diamond(i, &diamond) && is_part_of_same_irreducible_triangle(next, diamond))) {
                    for(k = i + 1; k < current_number_of_vertices - 1; k++) {
                        if(k != next) { // k != i is automatically implied because k >= i+1
                            for(l = 0; l < REG; l++) {
                                next0 = current_graph[k][l];
                                if(k < next0 && next0 != next) {
                                    if(!(is_part_of_irreducible_triangle_diamond(k, &diamond) && is_part_of_same_irreducible_triangle(next0, diamond))) {
                                        //Add edgepair i next k next0
                                        edge_pairs_list[*edge_pair_list_size][0] = i;
                                        edge_pairs_list[*edge_pair_list_size][1] = next;
                                        edge_pairs_list[*edge_pair_list_size][2] = k;
                                        edge_pairs_list[*edge_pair_list_size][3] = next0;
                                        if(!inserted_edge_will_be_part_of_square(edge_pairs_list[*edge_pair_list_size])
                                                && new_edge_has_min_colour_no_squares(edge_pairs_list[*edge_pair_list_size])) {
                                            int index0 = edge_labels[i][next];
                                            int index1 = edge_labels[k][next0];
                                            edgepair_index[index0][index1] = *edge_pair_list_size;

                                            (*edge_pair_list_size)++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 * Generate all possible edgepairs which may have minimal colour.
 * Remark: it is assumed that current_graph contains no triangles, but it may contain diamonds.
 */
void generate_edgepairs_no_triangles(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size) {
    DEBUGASSERT(number_of_reducible_triangles == 0);

    generate_edgepairs_penultimate_level_girth5_no_triangles_but_diamonds(edge_pairs_list, edge_pair_list_size);
    edgepair_list_index_squares = *edge_pair_list_size;

    update_min_edges();

    //Generate edgepairs which are part of a diamond and which can yield a minimal edge
    //The edge connecting 2 diamonds will never be minimal, only possibilities are ep's 01 23 and 02 13
    int i, j, k, l, next, next0;
    for(i = 0; i < number_of_irreducible_triangles; i++) {
        edge_pairs_list[*edge_pair_list_size][0] = irreducible_triangles[i][0];
        edge_pairs_list[*edge_pair_list_size][1] = irreducible_triangles[i][1];
        edge_pairs_list[*edge_pair_list_size][2] = irreducible_triangles[i][2];
        edge_pairs_list[*edge_pair_list_size][3] = irreducible_triangles[i][3];
        transform_edgepair_into_canonical_form(edge_pairs_list[*edge_pair_list_size]);

        edgepair_index[edge_labels[edge_pairs_list[*edge_pair_list_size][0]][edge_pairs_list[*edge_pair_list_size][1]]][edge_labels[edge_pairs_list[*edge_pair_list_size][2]][edge_pairs_list[*edge_pair_list_size][3]]] = *edge_pair_list_size;
        (*edge_pair_list_size)++;

        //Will always be in the same orbit as 01 23
        //But must be generated else there will be problems when determining the orbits
        edge_pairs_list[*edge_pair_list_size][0] = irreducible_triangles[i][0];
        edge_pairs_list[*edge_pair_list_size][1] = irreducible_triangles[i][2];
        edge_pairs_list[*edge_pair_list_size][2] = irreducible_triangles[i][1];
        edge_pairs_list[*edge_pair_list_size][3] = irreducible_triangles[i][3];
        transform_edgepair_into_canonical_form(edge_pairs_list[*edge_pair_list_size]);

        edgepair_index[edge_labels[edge_pairs_list[*edge_pair_list_size][0]][edge_pairs_list[*edge_pair_list_size][1]]][edge_labels[edge_pairs_list[*edge_pair_list_size][2]][edge_pairs_list[*edge_pair_list_size][3]]] = *edge_pair_list_size;
        (*edge_pair_list_size)++;

        if(current_number_of_vertices == 4) {
            edge_pairs_list[*edge_pair_list_size][0] = irreducible_triangles[i][0];
            edge_pairs_list[*edge_pair_list_size][1] = irreducible_triangles[i][3];
            edge_pairs_list[*edge_pair_list_size][2] = irreducible_triangles[i][1];
            edge_pairs_list[*edge_pair_list_size][3] = irreducible_triangles[i][2];
            transform_edgepair_into_canonical_form(edge_pairs_list[*edge_pair_list_size]);

            edgepair_index[edge_labels[edge_pairs_list[*edge_pair_list_size][0]][edge_pairs_list[*edge_pair_list_size][1]]][edge_labels[edge_pairs_list[*edge_pair_list_size][2]][edge_pairs_list[*edge_pair_list_size][3]]] = *edge_pair_list_size;
            (*edge_pair_list_size)++;
        }
    }

    //Generate the edgepairs that contain no edge which is fully in a diamond
    int diamond;
    for(i = 0; i < current_number_of_vertices - 3; i++) {
        for(j = 0; j < REG; j++) {
            next = current_graph[i][j];
            if(i < next) {
                if(!(is_part_of_irreducible_triangle_diamond(i, &diamond) && is_part_of_same_irreducible_triangle(next, diamond))) {
                    for(k = i + 1; k < current_number_of_vertices - 1; k++) {
                        if(k != next) { // k != i is automatically implied because k >= i+1
                            for(l = 0; l < REG; l++) {
                                next0 = current_graph[k][l];
                                if(k < next0 && next0 != next) {
                                    if(!(is_part_of_irreducible_triangle_diamond(k, &diamond) && is_part_of_same_irreducible_triangle(next0, diamond))) {
                                        //Add edgepair i next k next0
                                        edge_pairs_list[*edge_pair_list_size][0] = i;
                                        edge_pairs_list[*edge_pair_list_size][1] = next;
                                        edge_pairs_list[*edge_pair_list_size][2] = k;
                                        edge_pairs_list[*edge_pair_list_size][3] = next0;

                                        if(inserted_edge_will_be_part_of_square(edge_pairs_list[*edge_pair_list_size])
                                                && new_edge_has_min_colour(edge_pairs_list[*edge_pair_list_size])) {
                                            int index0 = edge_labels[i][next];
                                            int index1 = edge_labels[k][next0];

                                            edgepair_index[index0][index1] = *edge_pair_list_size;

                                            (*edge_pair_list_size)++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 * Returns the neighbours on distance one of the edge.
 */
setword get_neighbours_distance_one(EDGE edge) {
    return vertex_neighbourhood[edge[0]] | vertex_neighbourhood[edge[1]];
}

/**
 * Lookahead for the cheap colour (i.e. the number of vertices on distance <= 2 the edge).
 * Remark: it is assumed that the inserted edge will be part of a square.
 *
 * Returns 0 if the new edge will certainly not have the min colour.
 * Returns 1 if it MIGHT have the min colour.
 */
int new_edge_has_min_colour(EDGEPAIR edge_pair) {
    //If the current min colour == the max colour, then the inserted edge will also have the min colour
    if(min_edges_size > 0) {
        DEBUGASSERT(min_edge_is_part_of_square > -1);

        if(min_edge_is_part_of_square) { //i.e. inserted_edge_part_of_square == min_edge_is_part_of_square
            EDGE edge;
            edge[0] = edge_pair[0];
            edge[1] = edge_pair[1];
            setword neighbours0 = get_neighbours_distance_one(edge);

            edge[0] = edge_pair[2];
            edge[1] = edge_pair[3];
            setword neighbours1 = get_neighbours_distance_one(edge);

            int colour_inserted = MAX_EDGE_COLOUR_TWO;

            int number_overlapping = POPCOUNT(neighbours0 & neighbours1);

            DEBUGASSERT(number_overlapping >= 0 && number_overlapping <= 6);

            colour_inserted -= number_overlapping;

            if(colour_inserted > min_colour_one) {
                int i;
                for(i = 0; i < min_edges_size; i++) {
                    if((neighbours0 & min_edge_bitvectors[i]) == 0 && (neighbours1 & min_edge_bitvectors[i]) == 0) {
                        return 0;
                    }
                }

            } else
                return 1;

            //Also checking number of squares isn't any faster

        }

        return 1;
    } else {
        return 1;
    }
}

/**
 * Lookahead for the cheap colour (i.e. the number of vertices on distance <= 2 the edge).
 * Remark: it is assumed that the graph after extension won't contain any squares.
 *
 * Returns 0 if the new edge will certainly not have the min colour.
 * Returns 1 if it MIGHT have the min colour.
 */
int new_edge_has_min_colour_no_squares(EDGEPAIR edge_pair) {
    if(min_edges_size > 0) {
        DEBUGASSERT(min_edge[0] < current_number_of_vertices && min_edge[1] < current_number_of_vertices);

        EDGE edge;
        edge[0] = edge_pair[0];
        edge[1] = edge_pair[1];
        setword neighbours0 = get_neighbours_distance_one(edge);

        edge[0] = edge_pair[2];
        edge[1] = edge_pair[3];
        setword neighbours1 = get_neighbours_distance_one(edge);

        int colour_inserted = MAX_EDGE_COLOUR_TWO;

        int number_overlapping = POPCOUNT(neighbours0 & neighbours1);

        DEBUGASSERT(number_overlapping >= 0 && number_overlapping <= 6);

        colour_inserted -= number_overlapping;

        if(colour_inserted > min_colour_one) {
            int i;
            for(i = 0; i < min_edges_size; i++) {
                if((neighbours0 & min_edge_bitvectors[i]) == 0 && (neighbours1 & min_edge_bitvectors[i]) == 0) {
                    return 0;
                }
            }
        } else
            return 1;
        return 1;
    } else {
        return 1;
    }
}

/*****************Methods to test if the graph is colourable*******************/

void init_is_colourable(unsigned char number_of_colours[]) {
    RESETMARKS_SNARKS;
    int i, j;
    for (i = 0; i < current_number_of_vertices; i++) {
        number_of_colours[i] = 0;
        for(j = 0; j < REG; j++) {
            neighbour_index[i][current_graph[i][j]] = j;
        }
    }

}

/**
 * Returns 1 if current_graph is 3-edge colourable, else returns 0.
 */
int is_colourable() {
    if(number_of_bridges > 0)
        return 0;

    init_is_colourable(number_of_colours_snarks);

    int current_vertex = 0;
    int i, neighbour, current_index;
    for(i = 0; i < REG; i++) {
        colours_snarks[current_vertex][i] = i + 1;
        neighbour = current_graph[current_vertex][i];
        current_index = neighbour_index[neighbour][current_vertex];
        colours_snarks[neighbour][current_index] = i + 1;

        MARK_SNARKS(current_vertex, i);
        MARK_SNARKS(neighbour, current_index);
        number_of_colours_snarks[neighbour] = 1;
    }
    number_of_colours_snarks[current_vertex] = 3;
    return colour_next_free_choice(3);
}

/**
 * Returns 1 if the colours were swapped and returns 0 if the colours weren't swapped
 * because the cycle is hamiltonian.
 */
int swap_colours(unsigned char start_vertex) {
    unsigned char cycle[current_number_of_vertices];
    int cycle_size;

    EDGE start_edge;
    start_edge[0] = start_vertex;
    start_edge[1] = current_graph[start_vertex][0];

    EDGE colours, available_colours;

    unsigned char colour = colours_snarks[start_edge[0]][neighbour_index[start_edge[0]][start_edge[1]]];
    colours[0] = colour;
    determine_available_colours(colour, available_colours);
    colours[1] = available_colours[0];

    determine_even_cycle(start_edge, colours, cycle, &cycle_size);

    if(cycle_size == current_number_of_vertices) {
        //colour and available_colours[1]
        colours[1] = available_colours[1];
        determine_even_cycle(start_edge, colours, cycle, &cycle_size);

        if(cycle_size < current_number_of_vertices) {
            int i;
            for(i = 0; i < cycle_size; i++) {
                colours_snarks[cycle[i]][neighbour_index[cycle[i]][cycle[(i + 1) % cycle_size]]] = colours[(i + 1) % 2];
                colours_snarks[cycle[(i + 1) % cycle_size]][neighbour_index[cycle[(i + 1) % cycle_size]][cycle[i]]] = colours[(i + 1) % 2];
            }
            return 1;
        } else
            return 0;
    } else {
        //Swap
        int i;
        for(i = 0; i < cycle_size; i++) {
            colours_snarks[cycle[i]][neighbour_index[cycle[i]][cycle[(i + 1) % cycle_size]]] = colours[(i + 1) % 2];
            colours_snarks[cycle[(i + 1) % cycle_size]][neighbour_index[cycle[(i + 1) % cycle_size]][cycle[i]]] = colours[(i + 1) % 2];
        }

        colours[0] = colours[1];
        colours[1] = available_colours[1];
        determine_even_cycle(start_edge, colours, cycle, &cycle_size);

        //Another swap
        if(cycle_size < current_number_of_vertices)
            for(i = 0; i < cycle_size; i++) {
                colours_snarks[cycle[i]][neighbour_index[cycle[i]][cycle[(i + 1) % cycle_size]]] = colours[(i + 1) % 2];
                colours_snarks[cycle[(i + 1) % cycle_size]][neighbour_index[cycle[(i + 1) % cycle_size]][cycle[i]]] = colours[(i + 1) % 2];
            }

        return 1;
    }
}

/**
 * Determines a new colouring in a cheap way by swapping the colours of some colour chains.
 */
int modify_existing_colouring() {
    //modify_colouring(0);

    //return swap_colours(0);

    int i;
    for(i = 0; i < current_number_of_vertices; i++) {
        if(swap_colours(i))
            return 1;
    }
    //fprintf(stderr, "only very rarely happens\n");

    return 0;
}

void init_is_colourable_other(unsigned char number_of_colours[]) {
    RESETMARKS_SNARKS;
    int i;
    for (i = 0; i < current_number_of_vertices; i++) {
        number_of_colours[i] = 0;
        //neighbour_index was already set before
    }
}

/**
 * Returns 1 if current_graph is 3-edge colourable using another
 * colouring then is_colourable, else returns 0.
 *
 * Important: it is assumed that current_graph contains no triangles
 */
int is_colourable_other_colouring() {
    init_is_colourable_other(number_of_colours_snarks);

    int current_vertex = 0;
    int i, neighbour, current_index;
    for(i = 0; i < REG; i++) {
        colours_snarks[current_vertex][i] = i + 1;
        neighbour = current_graph[current_vertex][i];
        current_index = neighbour_index[neighbour][current_vertex];
        colours_snarks[neighbour][current_index] = i + 1;

        MARK_SNARKS(current_vertex, i);
        MARK_SNARKS(neighbour, current_index);
        number_of_colours_snarks[neighbour] = 1;
    }
    number_of_colours_snarks[current_vertex] = 3;


    //Guarantee that other colouring is different from first colouring
    for(current_vertex = 1; current_vertex < current_number_of_vertices; current_vertex++) {
        if(number_of_colours_snarks[current_vertex] == 1) {
            break;
        }
    }

    int used_colour;
    for(i = 0; i < REG; i++) {
        if(ISMARKED_SNARKS(current_vertex, i)) {
            used_colour = colours_snarks[current_vertex][i];
            break;
        }
    }

    EDGE available_vertices;
    int index_available_vertex0 = (i + 1) % 3;
    int index_available_vertex1 = (i + 2) % 3;
    available_vertices[0] = current_graph[current_vertex][index_available_vertex0];
    available_vertices[1] = current_graph[current_vertex][index_available_vertex1];

    EDGE available_colours;
    determine_available_colours(used_colour, available_colours);

    //Never a conflict since current_graph contains no triangles
    if(is_conflicting_colouring(colours_snarks, available_vertices[0], available_colours[1])
            || is_conflicting_colouring(colours_snarks, available_vertices[1], available_colours[0])) {
        fprintf(stderr, "Error: confliciting colouring (should never happen)\n");
        exit(1);
    }

    int index_current_vertex0 = neighbour_index[available_vertices[0]][current_vertex];
    int index_current_vertex1 = neighbour_index[available_vertices[1]][current_vertex];

    colours_snarks[current_vertex][index_available_vertex0] = available_colours[1];
    colours_snarks[current_vertex][index_available_vertex1] = available_colours[0];
    colours_snarks[available_vertices[0]][index_current_vertex0] = available_colours[1];
    colours_snarks[available_vertices[1]][index_current_vertex1] = available_colours[0];

    MARK_SNARKS(current_vertex, index_available_vertex0);
    MARK_SNARKS(current_vertex, index_available_vertex1);
    MARK_SNARKS(available_vertices[0], index_current_vertex0);
    MARK_SNARKS(available_vertices[1], index_current_vertex1);

    number_of_colours_snarks[current_vertex] = 3;
    number_of_colours_snarks[available_vertices[0]]++;
    number_of_colours_snarks[available_vertices[1]]++;


    for(i = 0; i < 2; i++) {
        if(number_of_colours_snarks[available_vertices[i]] != 1) {
            fprintf(stderr, "Error: Available vertex has unexpected number of coloured edges\n");
            exit(1);
        }
    }

    return colour_next_free_choice(5);
}

/**
 * Returns 1 if current_vertex already has a coloured edge with colour "colour",
 * else returns 0.
 */
int is_conflicting_colouring(unsigned char colours[][REG], int current_vertex, int colour) {
    int i;
    for(i = 0; i < REG; i++) {
        if(ISMARKED_SNARKS(current_vertex, i) && colours[current_vertex][i] == colour)
            return 1;
    }
    return 0;

}

/**
 * Determines the 2 available colours which are different from used_colour.
 */
void determine_available_colours(int used_colour, EDGE available_colours) {
    switch(used_colour) {
        case 1:
            available_colours[0] = 2;
            available_colours[1] = 3;
            break;
        case 2:
            available_colours[0] = 1;
            available_colours[1] = 3;
            break;
        case 3:
            available_colours[0] = 1;
            available_colours[1] = 2;
            break;
        default:
            fprintf(stderr, "Error: invalid previous colour\n");
            exit(1);
    }
}

/**
 * It is assumed that vertex already has 2 coloured incident edges. This method
 * determines the uncoloured neighbour of vertex and the missing colour.
 */
void determine_uncoloured_vertex(int vertex, int *uncoloured_vertex, int *missing_colour) {
    DEBUGASSERT(number_of_colours_snarks[vertex] == 2);

    int i;
    int sum_colours = 0;
    for(i = 0; i < REG; i++) {
        if(!ISMARKED_SNARKS(vertex, i)) {
            *uncoloured_vertex = current_graph[vertex][i];
        } else {
            sum_colours += colours_snarks[vertex][i];
        }
    }
    switch(sum_colours) {
        case 3:
            *missing_colour = 3;
            break;
        case 4:
            *missing_colour = 2;
            break;
        case 5:
            *missing_colour = 1;
            break;
        default:
            fprintf(stdout, "Error: invalid sum_colours\n");
            exit(1);
    }

}

/**
 * Resets the colour of the edges which were coloured by label_nonfree_choices.
 */
void unmark_colours(EDGE nonfree_labelled[], int nonfree_labelled_size) {
    int i;
    int vertex0, vertex1;
    for(i = 0; i < nonfree_labelled_size; i++) {
        vertex0 = nonfree_labelled[i][0];
        vertex1 = nonfree_labelled[i][1];
        UNMARK_SNARKS(vertex0, neighbour_index[vertex0][vertex1]);
        UNMARK_SNARKS(vertex1, neighbour_index[vertex1][vertex0]);
        number_of_colours_snarks[vertex0]--;
        number_of_colours_snarks[vertex1]--;
    }
}

/**
 * Labels all edges with the colour that was imposed by the current colouring
 * (i.e. if a vertex has 2 colours, it's third neighbour must be coloured with the
 * missing colour).
 *
 * The edges which were coloured like this, are put in nonfree_labelled.
 */
int label_nonfree_choices(int current_vertex, EDGE nonfree_labelled[], int *nonfree_labelled_size) {
    int uncoloured_vertex;
    int missing_colour;
    while(number_of_colours_snarks[current_vertex] == 2) {
        determine_uncoloured_vertex(current_vertex, &uncoloured_vertex, &missing_colour);
        if(!is_conflicting_colouring(colours_snarks, uncoloured_vertex, missing_colour)) {
            int index_uncoloured_vertex = neighbour_index[current_vertex][uncoloured_vertex];
            int index_current_vertex = neighbour_index[uncoloured_vertex][current_vertex];
            colours_snarks[current_vertex][index_uncoloured_vertex] = missing_colour;
            colours_snarks[uncoloured_vertex][index_current_vertex] = missing_colour;

            DEBUGASSERT(missing_colour > 0 && missing_colour < 4);

            MARK_SNARKS(current_vertex, index_uncoloured_vertex);
            MARK_SNARKS(uncoloured_vertex, index_current_vertex);
            number_of_colours_snarks[current_vertex] = 3;
            number_of_colours_snarks[uncoloured_vertex]++;

            nonfree_labelled[*nonfree_labelled_size][0] = current_vertex;
            nonfree_labelled[*nonfree_labelled_size][1] = uncoloured_vertex;
            (*nonfree_labelled_size)++;

            current_vertex = uncoloured_vertex;
        } else {
            unmark_colours(nonfree_labelled, *nonfree_labelled_size);
            return 0;
        }
    }
    return 1;

}

/**
 * Main (recursive) method of the colour routine.
 */
int colour_next_free_choice(int number_of_coloured_edges) {
    if(number_of_coloured_edges != current_number_of_edges) {
        int current_vertex;
        for(current_vertex = 1; current_vertex < current_number_of_vertices; current_vertex++) {
            if(number_of_colours_snarks[current_vertex] == 1) {
                break;
            }
            DEBUGASSERT(number_of_colours_snarks[current_vertex] != 2);
        }
        DEBUGASSERT(current_vertex < current_number_of_vertices);

        int used_colour;
        int i;
        for(i = 0; i < REG; i++) {
            if(ISMARKED_SNARKS(current_vertex, i)) {
                used_colour = colours_snarks[current_vertex][i];
                break;
            }
        }
        DEBUGASSERT(i < REG);

        EDGE available_vertices;
        int index_available_vertex0 = (i + 1) % 3;
        int index_available_vertex1 = (i + 2) % 3;
        available_vertices[0] = current_graph[current_vertex][index_available_vertex0];
        available_vertices[1] = current_graph[current_vertex][index_available_vertex1];

        EDGE available_colours;
        determine_available_colours(used_colour, available_colours);

        DEBUGASSERT((available_colours[0] > 0 && available_colours[0] < 4) && (available_colours[1] > 0 && available_colours[1] < 4));

        EDGE nonfree_labelled[current_number_of_edges - number_of_coloured_edges];
        int j;
        for(i = 0; i < 2; i++) {
            if(is_conflicting_colouring(colours_snarks, available_vertices[0], available_colours[i]) ||
                    is_conflicting_colouring(colours_snarks, available_vertices[1], available_colours[(i + 1) % 2])) {
                continue;
            }

            int index_current_vertex0 = neighbour_index[available_vertices[0]][current_vertex];
            int index_current_vertex1 = neighbour_index[available_vertices[1]][current_vertex];

            colours_snarks[available_vertices[0]][index_current_vertex0] = available_colours[i];
            colours_snarks[available_vertices[1]][index_current_vertex1] = available_colours[(i + 1) % 2];

            MARK_SNARKS(available_vertices[0], index_current_vertex0);
            MARK_SNARKS(available_vertices[1], index_current_vertex1);

            //number_of_colours_snarks[current_vertex] = 3;
            number_of_colours_snarks[available_vertices[0]]++;
            number_of_colours_snarks[available_vertices[1]]++;

            int nonfree_labelled_size = 0;
            int abort = 0;
            for(j = 0; j < 2;j++) {
                if(!label_nonfree_choices(available_vertices[j], nonfree_labelled, &nonfree_labelled_size)) {
                    DEBUGASSERT(nonfree_labelled_size <= current_number_of_edges - number_of_coloured_edges);
                    abort = 1;
                    break;
                }
                DEBUGASSERT(nonfree_labelled_size <= current_number_of_edges - number_of_coloured_edges);
            }

            if(!abort) {
                colours_snarks[current_vertex][index_available_vertex0] = available_colours[i];
                colours_snarks[current_vertex][index_available_vertex1] = available_colours[(i + 1) % 2];
                MARK_SNARKS(current_vertex, index_available_vertex0);
                MARK_SNARKS(current_vertex, index_available_vertex1);
                number_of_colours_snarks[current_vertex] = 3;

                if(colour_next_free_choice(number_of_coloured_edges + nonfree_labelled_size + 2)) {
                    return 1;
                } else {
                    unmark_colours(nonfree_labelled, nonfree_labelled_size);
                }
                UNMARK_SNARKS(current_vertex, index_available_vertex0);
                UNMARK_SNARKS(current_vertex, index_available_vertex1);
                number_of_colours_snarks[current_vertex] = 1;
            }
            UNMARK_SNARKS(available_vertices[0], index_current_vertex0);
            UNMARK_SNARKS(available_vertices[1], index_current_vertex1);

            //number_of_colours_snarks[current_vertex] = 1;
            number_of_colours_snarks[available_vertices[0]]--;
            number_of_colours_snarks[available_vertices[1]]--;
        }
        return 0;
    } else {
        return 1;
    }
}


/**********************Methods to test connectivity****************************/

void label_dfs_cut(GRAPH graph, int vertex, int previous_vertex, int number[], int *new_number, int *cut, int *lowpt) {
    int i, neighbour, number_neighbour;
    int locallowpt_vertex, lowpt_neighbour;

    number[vertex] = locallowpt_vertex = *new_number;
    (*new_number)++;
    for(i = 0; i < 3; i++) {
        neighbour = graph[vertex][i];
        if((number_neighbour = number[neighbour])) {
            if((neighbour != previous_vertex) && (number_neighbour < locallowpt_vertex))
                locallowpt_vertex = number_neighbour;
        } else {
            label_dfs_cut(graph, neighbour, vertex, number, new_number, cut, &lowpt_neighbour);
            if(lowpt_neighbour < locallowpt_vertex)
                locallowpt_vertex = lowpt_neighbour;
        }
    }

    if(locallowpt_vertex == number[vertex])
        *cut = 1;
    *lowpt = locallowpt_vertex;

}

/**
 * Checks if the graph without vertices forbidden_vertex1 and forbidden_vertex2 has a bridge.
 * Returns 0 if it has a bridge, else returns 1.
 *
 * If both forbidden_vertex1 >= 0 and forbidden_vertex2 >= 0, forbidden_vertex1
 * will always be < forbidden_vertex2.
 * If forbidden_vertex1 and forbidden_vertex2 are both < 0, no vertices will be removed.
 * If only forbidden_vertex2 < 0 (and forbidden_vertex1 >= 0), only forbidden_vertex1
 * will be removed.
 * If forbidden_vertex1 < 0, forbidden_vertex2 must also be < 0.
 *
 * Remark: because we're dealing with cubic graphs, a graph has a bridge iff
 * it has a cutvertex.
 */
int is_twoconnected(GRAPH graph, int forbidden_vertex1, int forbidden_vertex2) {
    int i, number[number_of_vertices], start;
    int cut, dummy, nextnumber;

    cut = 0;
    for(i = 1; i < number_of_vertices; i++) number[i] = 0;
    if(forbidden_vertex1 >= 0) {
        number[forbidden_vertex1] = INT_MAX;
        if(forbidden_vertex2 >= 0) {
            if(forbidden_vertex1 >= forbidden_vertex2) {
                fprintf(stderr, "Error: forbidden_vertex1 must always be smaller than forbidden_vertex2\n");
                exit(1);
            }
            number[forbidden_vertex2] = INT_MAX;
        }
        if(forbidden_vertex1 == 0) {
            for(start = 1; start < number_of_vertices; start++)
                if(start != forbidden_vertex2)
                    break;
        } else
            start = 0;
        number[start] = 1;
        nextnumber = 2;
        if(graph[start][0] != forbidden_vertex1 && graph[start][0] != forbidden_vertex2)
            label_dfs_cut(graph, graph[start][0], start, number, &nextnumber, &cut, &dummy);
        else if(graph[start][1] != forbidden_vertex1 && graph[start][1] != forbidden_vertex2)
            label_dfs_cut(graph, graph[start][1], start, number, &nextnumber, &cut, &dummy);
        else
            label_dfs_cut(graph, graph[start][2], start, number, &nextnumber, &cut, &dummy);
        if(cut || (number[graph[start][1]] == 0) || (number[graph[start][2]] == 0))
            return 0;
        /* ook een goede test als graaf[start][0]==top */
    } else {
        number[0] = 1;
        nextnumber = 2;
        label_dfs_cut(graph, graph[0][0], 0, number, &nextnumber, &cut, &dummy);
        if(cut || (number[graph[0][1]] == 0) || (number[graph[0][2]] == 0))
            return 0;
    }

    return 1;

}

/**
 * Returns 1 if graph has a twocut, else returns 0.
 *
 * Algorithm: removes 1 vertex and tests if the resulting graph has a bridge.
 */
int has_twocut(GRAPH graph) {
    int i;

    for(i = 0; i < number_of_vertices; i++)
        if(!is_twoconnected(graph, i, -1))
            return 1;

    return 0;
}

/**
 * Returns 1 if current_graph has a nontrivial threecut, else returns 0.
 *
 * Algorithm: removes 2 vertices and tests if the resulting graph has a bridge.
 *
 * This algorithm is not very efficient but it's simple.
 */
int has_nontrivial_threecut() {
    //If the graph contains triangles, it certainly has a nontrivial threecut
    if(number_of_reducible_triangles > 0 || number_of_irreducible_triangles > 0)
        return number_of_vertices > 4; //No nontrivial threecut in case of the K4

    //vertex_neighbourhoods are not updated on last level of triangle_extend, but no problem since if the graph contains
    //triangles, it certainly has a nontrivial threecut

    int i, j;
    for(i = 0; i < number_of_vertices - 1; i++) {
        for(j = i + 1; j < number_of_vertices; j++) {
            if((vertex_neighbourhood[i] & vertex_neighbourhood[j]) == 0) { //Don't continue if there is a common neighbour: if common neighbour, there will be a trivial threecut
                //Don't test connectivity if i and j are neighbours, because if it's not twoconnected after
                //the removal of i and j, there will also be other nonadjacent cutvertices which will yield a bridge
                if((vertex_neighbourhood[i] & BIT(j)) == 0 && !is_twoconnected(current_graph, i, j)) {
                    return 1;
                }
            }
        }
    }

    return 0;
}


/***************************Auxiliary methods**********************************/


/**
 * Code the cubic graph using the multicode format.
 */
void code_multicode(unsigned char *code, GRAPH g, int num_of_vertices) {
    int j;
    unsigned char *p;

    *(code++) = num_of_vertices;
    for(j = 0, p = g[0]; j < num_of_vertices - 1; j++) {
        if(j < *p)
            *(code++) = (*p) + 1;
        p++;
        if(j < *p)
            *(code++) = (*p) + 1;
        p++;
        if(j < *p)
            *(code++) = (*p) + 1;
        p += 2; //Because graph is REG + 1: element at last position is empty
        *(code++) = 0;
    }
}

/**
 * Code the graph using the graph6 format.
 */
void code_graph6(unsigned char *code, GRAPH g, int num_of_vertices) {
    int i, j, k, org;
    int nlen, bodylen;
    static unsigned char g6bit[] = {32, 16, 8, 4, 2, 1};
    register unsigned char *body;

    if(num_of_vertices <= 62) {
        code[0] = 63 + num_of_vertices;
        nlen = 1;
    } else {
        code[0] = 63 + 63;
        code[1] = 63 + 0; //Will always be zero if num_of_vertices < 4096
        code[2] = 63 + (num_of_vertices >> 6);
        code[3] = 63 + (num_of_vertices & 0x3F);
        nlen = 4;
    }

    body = code + nlen;
    bodylen = ((num_of_vertices * (num_of_vertices - 1)) / 2 + 5) / 6;
    for(i = 0; i < bodylen; ++i)
        body[i] = 63;
    body[bodylen] = '\n';

    for(i = org = 0; i < num_of_vertices; org += i, ++i) {
        j = g[i][0];
        if(j < i) {
            k = org + j;
            body[k / 6] += g6bit[k % 6];
        }
        j = g[i][1];
        if(j < i) {
            k = org + j;
            body[k / 6] += g6bit[k % 6];
        }
        j = g[i][2];
        if(j < i) {
            k = org + j;
            body[k / 6] += g6bit[k % 6];
        }
    }
}

/**
 * Writes the codes from the graphlist to the file.
 * All graphs with the same number of vertices are written to the same file.
 */
void wegspeichern(unsigned char *liste[MAXN], int num_of_vertices) {
    if(fwrite(liste[num_of_vertices], codelength[num_of_vertices], graphlist_number_of_graphs[num_of_vertices], outputfile[num_of_vertices])
            < graphlist_number_of_graphs[num_of_vertices]) {
        fprintf(stderr, "Error: couldn't write all graphs! \n\n");
        exit(1);
    }
}

/**
 * Writes the code of a prime graph to it's outputfile.
 */
void wegspeichern_irred(unsigned char *liste, int codelength, int num_of_vertices) {
    if(fwrite(liste, codelength, 1, outputfile_irred) < 1) {
        fprintf(stderr, "Error: couldn't write irreducible graph! \n\n");
        exit(1);
    }
}

/**
 * Calculates the code of g and writes it to the right buffer of graphs to be output.
 * If the buffer is full, the codes of the graphs are actually output.
 *
 * If userproc_sh != NULL, this method is also called with the generated graph
 * as argument.
 *
 * It is assumed that the graphs are already filtered here (i.e. that g already
 * fulfills the necessary requirements, eg is a snark, has the right girth, ...).
 */
void aufschreiben_already_filtered(GRAPH g, int num_vertices) {
    if(num_vertices == number_of_vertices) {
        number_of_graphs_found++;
    }

    if(userproc_sh != NULL) {
        (*userproc_sh) (g, num_vertices);
    }

    if(!noout) {
        if(!graph6_output) //I.e. multicode
            code_multicode(graphlist[num_vertices] + (graphlist_number_of_graphs[num_vertices] * codelength[num_vertices]), g, num_vertices);
        else
            code_graph6(graphlist[num_vertices] + (graphlist_number_of_graphs[num_vertices] * codelength[num_vertices]), g, num_vertices);

        graphlist_number_of_graphs[num_vertices]++;
        if(graphlist_number_of_graphs[num_vertices] == listlength) {
            wegspeichern(graphlist, num_vertices);
            graphlist_number_of_graphs[num_vertices] = 0;
        }
    }
}

/**
 * First checks if current_graph fullfills the necessary requirements
 * (has the right girth and connectivity, is a snark, ...).
 * If it does, aufschreiben_already_filtered is called.
 */
void aufschreiben() {
    //In case of modulo, graphs below splitlevel are only written if rest == 0
    if(modulo && rest != 0 && current_number_of_vertices < splitlevel) {
        return;
    }

    //Only needed if current_number_of_vertices < number_of_vertices
    if(current_number_of_vertices < number_of_vertices && girth > 3) {
        if((number_of_reducible_triangles > 0 || number_of_irreducible_triangles > 0)
                || (girth == 5 && contains_squares()))
            return;
    }

    if(snarks) {
        if(is_colourable())
            return;
    }

    /**
     * Make sure snarks are checked before checking cyclic connectivity, because
     * C3 and expecially C4 are very expensive to check (much more expensive then
     * testing colourability).
     */
    if(check_cyclic_connectivity) {
        if((min_cyclic_connectivity >= 2 && number_of_bridges > 0)
                || (min_cyclic_connectivity >= 3 && has_twocut(current_graph))
                || (min_cyclic_connectivity >= 4 && has_nontrivial_threecut()))
            return;
    }

    if(!snarks || current_number_of_vertices < number_of_vertices) {
        aufschreiben_already_filtered(current_graph, current_number_of_vertices);
    } else { //snarks && current_number_of_vertices == number_of_vertices
        //
        /**
         * In case of snarks, the canonicity check of the graphs is delayed.
         * So the generated snarks may contain duplicates. Therefore they are
         * written in a seperate list and later in output_nonisomorphic_children
         * that list is filtered before actually outputting the snarks.
         */

        if(still_has_to_check_if_graph_is_canonical && !is_major_edge(0, &canon_graphs[graphlist_snarks_size])) {
            return;
        }

        DEBUGASSERT(graphlist_snarks_size < max_graphlist_snarks_size);
        //Already has to expand here, not at == max_graphlist_snarks_size
        //Because canon_graphs[graphlist_snarks_size] may already be set by is_major_edge_snarks
        if(graphlist_snarks_size == max_graphlist_snarks_size - 1) {
            GRAPH graphlist_snarks_temp[graphlist_snarks_size];
            memcpy(graphlist_snarks_temp, graphlist_snarks, graphlist_snarks_size * sizeof(GRAPH));

            free(graphlist_snarks);

            //Could also expand with a certain percentage
            max_graphlist_snarks_size += DEFAULT_MAX_GRAPHLIST_SNARKS_SIZE;

            graphlist_snarks = (GRAPH *) malloc(sizeof(GRAPH) * max_graphlist_snarks_size);
            if(graphlist_snarks == NULL) {
                fprintf(stderr, "Error: out of memory while expanding graphlist_snarks\n");
                exit(1);
            }

            memcpy(graphlist_snarks, graphlist_snarks_temp, graphlist_snarks_size * sizeof(GRAPH));


            //Also expand canon_graphs
            sparsegraph canon_graphs_temp[graphlist_snarks_size + 1];
            memcpy(canon_graphs_temp, canon_graphs, (graphlist_snarks_size + 1) * sizeof(sparsegraph));

            free(canon_graphs);

            canon_graphs = (sparsegraph *) malloc(sizeof(sparsegraph) * max_graphlist_snarks_size);
            if(canon_graphs == NULL) {
                fprintf(stderr, "Error: out of memory while expanding canon_graphs\n");
                exit(1);
            }

            memcpy(canon_graphs, canon_graphs_temp, (graphlist_snarks_size + 1) * sizeof(sparsegraph));

            //Init the new sparse graphs
            int i;
            for(i = graphlist_snarks_size + 1; i < max_graphlist_snarks_size; i++) {
                SG_INIT(canon_graphs[i]);
                SG_ALLOC(canon_graphs[i], number_of_vertices, 3 * number_of_vertices, "malloc");
            }
        }
        if(!still_has_to_check_if_graph_is_canonical) { //else canon_graphs[graphlist_snarks_size] was already filled
            //Can only be MAJOR_EDGE_INSERTED

            options.getcanon = TRUE;

            //Important: calling determine_major_edge_partitions() is no longer valid, since there may be isomorphic graphs!
            //determine_major_edge_partitions();

            determine_vertex_partitions();

            number_of_generators = 0;

            copy_sparse_graph();

            nauty_sh((graph*) & sg, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, current_number_of_vertices, (graph*) & canon_graphs[graphlist_snarks_size]);
        }

        /**
         * Could check if parent is actually canonical after 1 snarkchild has been found,
         * but since the percentage of cases where there are > 0 snarkchildren is small and even decreases
         * this is no bottleneck. So checking canonicity after 1 child only yield an extremely small speedup.
         */

        memcpy(graphlist_snarks[graphlist_snarks_size], current_graph, sizeof(GRAPH));
        //memcpy(graphlist_snarks + graphlist_snarks_size * sizeof(GRAPH), current_graph, sizeof(GRAPH));
        graphlist_snarks_size++;
    }
}

/**
 * Encodes and outputs a prime graph.
 */
void aufschreiben_irred() {
    DEBUGASSERT(!noout && output_irred_graphs && output_to_file);

    if(!graph6_output) {
        unsigned char code[(5 * current_number_of_vertices) / 2];
        code_multicode(code, current_graph, current_number_of_vertices);
        wegspeichern_irred(code, (5 * current_number_of_vertices) / 2, current_number_of_vertices);
    } else {
        int nlen = (current_number_of_vertices <= 62) ? 1 : 4;
        int codelength = nlen + ((current_number_of_vertices * (current_number_of_vertices - 1)) / 2 + 5) / 6 + 1;
        unsigned char code[codelength];
        code_graph6(code, current_graph, current_number_of_vertices);
        wegspeichern_irred(code, codelength, current_number_of_vertices);
    }

}

/**
 * Calculates all binomial coefficients with n <= max_n and saves them into
 * binom_coefficients.
 * Is used to determine the max number of vertexsets.
 */
void calculate_binom_coefficients(int max_n) {
  binom_coefficients[0][0] = 1;
  int i,j;
    for(i = 1; i <= max_n; i++) {
        binom_coefficients[i][i] = 1;
        for(j = i - 1; j > 0; j--)
            binom_coefficients[i][j] = binom_coefficients[i - 1][j] + binom_coefficients[i - 1][j - 1];
        binom_coefficients[i][0] = 1;
    }
}

int power(int base, int exponent) {
    int result = 1;
    int i;
    for(i = 0; i < exponent; i++) {
        result *= base;
    }
    return result;
}

/* Methods to print various types of graphs */

void printgraph_nauty(graph *g, int current_number_of_vertices) {
    int i;
    fprintf(stderr, "Printing graph nauty:\n");
    for(i = 0; i < current_number_of_vertices; i++) {
        fprintf(stderr, "%d :", i);
        set *gv = GRAPHROW(g, i, MAXM);
        int neighbour = -1;
        while((neighbour = nextelement(gv, MAXM, neighbour)) >= 0) {
            fprintf(stderr, " %d", neighbour);
        }
        fprintf(stderr, "\n");
    }
}

void print_sparse_graph_nauty(sparsegraph sparse_graph) {
    int i, j;
    fprintf(stderr, "Printing sparse graph nauty:\n");
    for(i = 0; i < sparse_graph.nv; i++) {
        fprintf(stderr, "%d :", i);
        for(j = 0; j < sparse_graph.d[i];j++) {
            //fprintf(stderr, " %d", sparse_graph.e[i * REG + j]);
            fprintf(stderr, " %d", sparse_graph.e[sparse_graph.v[i] + j]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "Number of directed edges: %lu\n", (unsigned long) sparse_graph.nde);
}

void printgraph() {
    int i, j;
    fprintf(stderr, "Printing graph:\n");
    for(i = 0; i < current_number_of_vertices; i++) {
        fprintf(stderr, "%d :", i);
        for(j = 0; j < REG; j++)
            fprintf(stderr, " %d", current_graph[i][j]);
        fprintf(stderr, "\n");
    }
    for(i = 0; i < number_of_irreducible_triangles; i++)
        fprintf(stderr, "%d %d %d %d\n", irreducible_triangles[i][0], irreducible_triangles[i][1], irreducible_triangles[i][2], irreducible_triangles[i][3]);
    for(i = 0; i < number_of_reducible_triangles; i++)
        fprintf(stderr, "%d %d %d\n", reducible_triangles[i][0], reducible_triangles[i][1], reducible_triangles[i][2]);
    for(i = 0; i < number_of_bridges;i++)
        fprintf(stderr, "%d %d\n", bridges[i][0], bridges[i][1]);
}

void printgraph_irred() {
    int i, j;
    fprintf(stderr, "Printing graph:\n");
    for(i = 0; i < current_number_of_vertices; i++) {
        fprintf(stderr, "%d :", i);
        for(j = 0; j < REG; j++)
            fprintf(stderr, " %d", current_graph[i][j]);
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "edge diamonds:\n");
    for(i = 0; i < number_of_edge_diamonds; i++)
        fprintf(stderr, "%d %d %d %d\n", edge_diamonds[i][0], edge_diamonds[i][1], edge_diamonds[i][2], edge_diamonds[i][3]);
    fprintf(stderr, "lollipops:\n");
    for(i = 0; i < number_of_lollipop_diamonds; i++)
        fprintf(stderr, "%d %d %d %d\n", lollipop_diamonds[i][0], lollipop_diamonds[i][1], lollipop_diamonds[i][2], lollipop_diamonds[i][3]);
    fprintf(stderr, "nonadj edge diamonds:\n");
    for(i = 0; i < number_of_nonadj_edge_diamonds;i++)
        fprintf(stderr, "%d %d %d %d\n", nonadj_edge_diamonds[i][0], nonadj_edge_diamonds[i][1], nonadj_edge_diamonds[i][2], nonadj_edge_diamonds[i][3]);
}

/* Methods for nauty */

/**
 * Initializes the nauty options and the nauty datastructures.
 */
void init_nauty_options() {
    //options.getcanon = TRUE;
    options.userautomproc = save_generators;

    /* Init the nauty datastructures */
    SG_INIT(sg);
    SG_ALLOC(sg, number_of_vertices, 3 * number_of_vertices, "malloc");

    //sg.v and sg.d only have to be set once
    int i;
    for(i = 0; i < number_of_vertices; i++) {
        sg.v[i] = i * REG;
        sg.d[i] = 3;
    }

    SG_INIT(sg_canon);
    SG_ALLOC(sg_canon, number_of_vertices, 3 * number_of_vertices, "malloc");
}

/**
 * Method which is called each time nauty finds a generator.
 */
void save_generators(int count, permutation perm[], nvector orbits[],
        int numorbits, int stabvertex, int n) {
    memcpy(generators + number_of_generators, perm, sizeof (permutation) * n);

    number_of_generators++;
}

//Just for testing
/*
void handle_3_regular_result(unsigned char snarkhunter_graph[MAXN][REG + 1], int order) {
    fprintf(stderr, "Graph generated with order %d\n", order);
}
*/

/**
 * Copies current_graph into sg.
 */
void copy_sparse_graph() {
    sg.nv = current_number_of_vertices;
    sg.nde = 3 * current_number_of_vertices;

    int i, j;
    for(i = 0; i < current_number_of_vertices;i++) {
        //These values were already set in init_nauty_options()
        //sg.v[i] = i * REG;
        //sg.d[i] = 3;
        for(j = 0; j < REG; j++) {
            sg.e[i * REG + j] = current_graph[i][j];
        }
    }
}

void print_help(char * argv0) {
    fprintf(stderr, "Usage: %s <number_of_vertices> <min_girth> [options]\n", argv0);
    fprintf(stderr, "At least the number of vertices and the minimum girth must be specified.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Possible options are:\n");
    fprintf(stderr, "  S: only output class 2 graphs (i.e. graphs with chromatic index 4)\n");
    fprintf(stderr, "  C<x>: only output graphs with cyclic edge connectivity >= x\n");
    fprintf(stderr, "  s: only output the graphs with 'n' vertices (instead of all graphs with <= n vertices)\n");
    fprintf(stderr, "  o: the graphs are written to stdout instead of to files\n");
    fprintf(stderr, "  g: the graphs will be output in graph6 format instead of multicode format\n");
    fprintf(stderr, "  n: no graphs will be output (i.e. the graphs are only counted)\n");
    fprintf(stderr, "  i: prime graphs will also be output (to a separate file)\n");
    fprintf(stderr, "  m <rest> <modulo>: splits the generation in <modulo> (more or less equally big) parts. Here part <rest> will be executed.\n");
}

void call_snarkhunter(int n, int g, void (*userproc) (unsigned char (*)[REG + 1], int)){
    /* Checks to test if the WORDSIZE and sizeof setwords is valid */
    if(WORDSIZE != 32 && WORDSIZE != 64) {
        fprintf(stderr, "Error: invalid value for wordsize: %d\n", WORDSIZE);
        fprintf(stderr, "Valid values are 32 and 64\n");
        exit(1);
    }

   if(MAXN != WORDSIZE) {
       fprintf(stderr, "Error: MAXN is not equal to WORDSIZE: %d vs %d\n", MAXN, WORDSIZE);
       fprintf(stderr, "Please compile with option -DMAXN=WORDSIZE\n");
       exit(1);
   }

    if(WORDSIZE == 32) {
        //fprintf(stderr, "32 bit mode\n");
        if(sizeof(unsigned int) < 4) {
            fprintf(stderr, "Error: unsigned ints should be at least 32 bit -- sorry.\n");
            exit(1);
        }
        if(sizeof(setword) < 4) {
            fprintf(stderr, "Error: this version relies on 32 bit setwords -- sorry.\n");
            exit(1);
        }
    } else {
        //fprintf(stderr, "64 bit mode\n");
        if(sizeof(setword) != 8) {
            fprintf(stderr, "Error: this version relies on 64 bit setwords -- sorry.\n");
            exit(1);
        }
        if(sizeof(unsigned long long int) != 8) {
            fprintf(stderr, "Error: this version relies on 64 bit unsigned long long ints -- sorry.\n");
            exit(1);
        }
    }
	
	number_of_vertices = n;
	girth = g;
    noout = singleout = 1;

	switch(girth) {
		case 3:
		{
			min_order = 4;
			break;
		}
		case 4:
		{
			min_order = 6;
			break;
		}
		case 5:
		{
			min_order = 10;
			break;
		}
		default:
		{
			fprintf(stderr, "Error: invalid minimum girth: %d (valid girths are 3, 4 or 5).\n", girth);
			exit(1);
		}
	}
	
    nauty_check(WORDSIZE, MAXM, MAXN, NAUTYVERSIONID);
	
    generate_irreducible_graphs(number_of_vertices, girth, userproc);
}


#ifndef SNARKHUNTERMAIN
  #define SNARKHUNTERMAIN main
#endif
int SNARKHUNTERMAIN(int argc, char *argv[]) {
    int i;
    char strbuffer[50];

    /* Checks to test if the WORDSIZE and sizeof setwords is valid */
    if(WORDSIZE != 32 && WORDSIZE != 64) {
        fprintf(stderr, "Error: invalid value for wordsize: %d\n", WORDSIZE);
        fprintf(stderr, "Valid values are 32 and 64\n");
        exit(1);
    }

   if(MAXN != WORDSIZE) {
       fprintf(stderr, "Error: MAXN is not equal to WORDSIZE: %d vs %d\n", MAXN, WORDSIZE);
       fprintf(stderr, "Please compile with option -DMAXN=WORDSIZE\n");
       exit(1);
   }

    if(WORDSIZE == 32) {
        fprintf(stderr, "32 bit mode\n");
        if(sizeof(unsigned int) < 4) {
            fprintf(stderr, "Error: unsigned ints should be at least 32 bit -- sorry.\n");
            exit(1);
        }
        if(sizeof(setword) < 4) {
            fprintf(stderr, "Error: this version relies on 32 bit setwords -- sorry.\n");
            exit(1);
        }
    } else {
        fprintf(stderr, "64 bit mode\n");
        if(sizeof(setword) != 8) {
            fprintf(stderr, "Error: this version relies on 64 bit setwords -- sorry.\n");
            exit(1);
        }
        if(sizeof(unsigned long long int) != 8) {
            fprintf(stderr, "Error: this version relies on 64 bit unsigned long long ints -- sorry.\n");
            exit(1);
        }
    }

    DEBUGMSG("Debug is on");


    if(argc < 3) {
        if(argc == 2 && argv[1][0] == '-' && argv[1][1] == 'h') {
            print_help(argv[0]);
        } else {
            fprintf(stderr, "Error: invalid number of arguments. At least the number of vertices and the minimum girth must be specified.\n");
            fprintf(stderr, "Usage: %s <number_of_vertices> <min_girth> [options]\n", argv[0]);
            fprintf(stderr, "Execute '%s -h' for more extensive help.\n", argv[0]);
        }
        exit(1);
    } else {
        if(sscanf(argv[1], "%d", &number_of_vertices) && (number_of_vertices % 2 == 0)) {
            if(number_of_vertices > MAXN) {
                if(MAXN == 32) {
                    fprintf(stderr, "Error: maximum number of vertices is 32 in this version.\n");
                    fprintf(stderr, "Compile with option -DWORDSIZE=64 in order to be able to generate graphs with a higher order\n");
                    fprintf(stderr, "Or alternatively use the command 'make 64bit'\n");
                    exit(1);
                } else {
                    fprintf(stderr, "Error: number of vertices is too big (limit is %d).\n", MAXN);
                    exit(1);
                }
            }
            if(MAXN == 64 && number_of_vertices <= 32) {
                fprintf(stderr, "Info: it is recommended to use the 32 bit version for generating graphs with order <= 32\n");
            }

            sscanf(argv[2], "%d", &girth);
        } else {
            fprintf(stderr, "Error: number of vertices should be an (even) integer.\n");
            exit(EXIT_FAILURE);
        }

        for(i = 3; i < argc; i++) {
            switch(argv[i][0]) {
                case 's':
                {
                    singleout = 1;
                    break;
                }
                case 'S':
                {
                    if(girth == 3) {
                        fprintf(stderr, "Snarks with girth 3 can be easily obtained from other snarks.\n");
                        fprintf(stderr, "Therefore this generator can only generate snarks with girth > 3.\n");
                        exit(1);
                    }
                    snarks = 1;
                    break;
                }
                case 'i':
                {
                    output_irred_graphs = 1;
                    break;
                }
                case 'o':
                {
                    output_to_file = 0;
                    break;
                }
                case 'n':
                {
                    //Overrides output to file
                    noout = singleout = 1;
                    break;
                }
                case 'm':
                {
                    modulo = 1;
                    i++;
                    rest = atoi(argv[i]);
                    i++;
                    mod = atoi(argv[i]);
                    if (rest >= mod) {
                        fprintf(stderr, "Error: rest (%d) must be smaller than mod (%d).\n", rest, mod);
                        exit(1);
                    }
                    break;
                }
                case 'C':
                {
                    check_cyclic_connectivity = 1;
                    min_cyclic_connectivity = atoi(argv[i] + 1);
                    if(min_cyclic_connectivity < 2 || min_cyclic_connectivity > 4) {
                        fprintf(stderr, "Error: cyclic edge connectivity should be at least 2 and at most 4.\n");
                        fprintf(stderr, "If you want the graphs with cyclic edge connectivity at least 1, you don't need to add the 'C'-option.\n");
                        exit(1);
                    }
                    if(girth == 3 && min_cyclic_connectivity == 4) {
                        fprintf(stderr, "Remark: graphs with girth 3 are never cyclically 4-edge connected, so setting girth to 4.\n");
                        girth = 4;
                    }
                    break;
                }
                case 'g':
                {
                    graph6_output = 1;
                    break;
                }
                default:
                {
                    fprintf(stderr, "Error: invalid option: %s\n", argv[i]);
                    fprintf(stderr, "Execute '%s -h' for a list of possible options.\n", argv[0]);
                    exit(1);
                }
            }
        }

        switch(girth) {
            case 3:
            {
                min_order = 4;
                break;
            }
            case 4:
            {
                min_order = 6;
                break;
            }
            case 5:
            {
                min_order = 10;
                break;
            }
            default:
            {
                fprintf(stderr, "Error: invalid minimum girth: %d (valid girths are 3, 4 or 5).\n", girth);
                exit(1);
            }
        }

        if(min_order > number_of_vertices) {
            fprintf(stderr, "Cubic graphs of the given order and girth do not exist.\n");
            exit(1);
        }

        if(output_irred_graphs && !output_to_file) {
            fprintf(stderr, "Warning: option 'i' and 'o' cannot be used at the same time because irreducible graphs cannot be output to stdout. Therefore parameter 'i' was ignored.\n");
            output_irred_graphs = 0;
        }
        if(noout && !output_to_file) {
            fprintf(stderr, "Warning: option 'o' was ignored because of option 'n'.\n");
        }
        if(noout && graph6_output) {
            fprintf(stderr, "Warning: option 'g' was ignored because of option 'n'.\n");
        }
    }

    if(modulo) {
        /**
         * The level (i.e. number of vertices) where the recursion tree is split for the different cases.
         * A higher level means a more equal distribution of the cases, but the common part will take longer.
         */
        splitlevel = number_of_vertices - 6;
        splitlevel = MAX(splitlevel, 6);
        //24 seems to be optimal for all girths (takes less than 5 minutes)
        splitlevel = MIN(splitlevel, 24);
        if(splitlevel % 2 != 0) {
            fprintf(stderr, "Error: splitlevel has to be even.\n");
            exit(1);
        }
        fprintf(stderr, "Splitlevel is %d.\n", splitlevel);
        //Splitlevel should be at least 2 smaller than the number of vertices
        if(splitlevel > number_of_vertices - 2) {
            fprintf(stderr, "Error: splitlevel is too big\n");
            exit(1);
        }
    }

    nauty_check(WORDSIZE, MAXM, MAXN, NAUTYVERSIONID);

    //For the timing
    struct tms TMS;
    unsigned int oldtime = 0;

    if(snarks) {
        graphlist_snarks = (GRAPH *) malloc(sizeof(GRAPH) * max_graphlist_snarks_size);
        if(graphlist_snarks == NULL) {
            fprintf(stderr, "Error: out of memory while creating graphlist_snarks\n");
            exit(1);
        }
        canon_graphs = (sparsegraph *) malloc(sizeof(sparsegraph) * max_graphlist_snarks_size);
        if(canon_graphs == NULL) {
            fprintf(stderr, "Error: out of memory while creating canon_graphs\n");
            exit(1);
        }
        for(i = 0; i < max_graphlist_snarks_size; i++) {
            SG_INIT(canon_graphs[i]);
            SG_ALLOC(canon_graphs[i], number_of_vertices, 3 * number_of_vertices, "malloc");
        }
    }

    if(!noout) {
        if(snarks) {
            listlength = LISTLENGTH_SNARKS;
            if(!singleout) {
                fprintf(stderr, "Warning: Generating snarks ('S') without the singleout ('s') option is significantly slower.\n");
            }

            /**
             * These values (which were empirically determined) seem to be optimal.
             */
            if(girth == 4) {
                min_edgepairlist_size = 0;
            } else if(girth >= 5) {
                min_edgepairlist_size = 9;
            }
        } else
            listlength = LISTLENGTH;

        i = min_order;
        if(singleout)
            i = number_of_vertices;
        while(i <= number_of_vertices) {
            if(graph6_output) {
                int nlen = (i <= 62) ? 1 : 4;
                codelength[i] = nlen + ((i * (i - 1)) / 2 + 5) / 6 + 1;
            } else //Codelength of the multigraph format (i.e. number of vertices + number of edges)
                codelength[i] = i + (REG * i) / 2;
            graphlist_number_of_graphs[i] = 0;
            graphlist[i] = (unsigned char *) malloc(codelength[i] * listlength);
            if(graphlist[i] == NULL) {
                fprintf(stderr, "Error: out of memory while creating graphlist\n");
                exit(1);
            }

            if(output_to_file) {
                sprintf(outputfilename[i], "Generated_graphs.00.00");

                if(snarks) {
                    sprintf(strbuffer, ".sn");
                    strcat(outputfilename[i], strbuffer);
                }

                if(check_cyclic_connectivity) {
                    sprintf(strbuffer, ".cyc%d", min_cyclic_connectivity);
                    strcat(outputfilename[i], strbuffer);
                }

                if(modulo) { /* must be the last item added to the name */
                    sprintf(strbuffer, ".m_%d_%d", rest, mod);
                    strcat(outputfilename[i], strbuffer);
                }

                if(graph6_output) {
                    sprintf(strbuffer, ".g6");
                    strcat(outputfilename[i], strbuffer);
                }
                outputfilename[i][17] = (i / 10) + 48;
                outputfilename[i][18] = (i % 10) + 48;

                outputfilename[i][21] = girth + 48;

                //To make sure the output appears in a new file instead of being appended to an existing file
                outputfile[i] = fopen(outputfilename[i], "w");
                if(outputfile == NULL) {
                    fprintf(stderr, "Error: could not create outputfile\n");
                    exit(1);
                }
            } else {
                outputfile[i] = stdout;
            }

            i += 2;
        }
    }


    for(i = 0; i < argc; i++)
        fprintf(stderr, "%s ", argv[i]);
    fprintf(stderr, "\n");

    generate_irreducible_graphs(number_of_vertices, girth, NULL);

    if(!noout && output_to_file) {
        i = min_order;
        if(singleout)
            i = number_of_vertices;
        while(i <= number_of_vertices) {
            fclose(outputfile[i]);
            i += 2;
        }
    }

    times(&TMS);
    unsigned int savetime = oldtime + (unsigned int) TMS.tms_utime;

    fprintf(stderr, "All done, %llu graphs generated with %d vertices.\n", number_of_graphs_found, number_of_vertices);
    fprintf(stderr, "CPU time: %.1f seconds.\n", (double) savetime / time_factor);
    //fprintf(stderr, "Number of nauty calls: %llu\n", nauty_calls);


    return (EXIT_SUCCESS);
}
