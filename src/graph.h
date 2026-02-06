#ifndef READ_GRAPH_H
#define READ_GRAPH_H

#include <stdlib.h>
#include <stdio.h>


/* 
 * Efficiently store a pair of 16-bit digits in a 32-bit integer
 * 
 * WARNING: this won't work if n > 65535
 *
 */
typedef int pair_t;
#define make_pair(X,Y) (X << 16) | (Y & 0xFFFF)
#define src(P) (P >> 16) & 0xFFFF
#define snk(P) (P&0xFFFF)
/* Does pair P = (X,Y) or (Y,X) ?*/
#define pair_eq(P,X,Y) ((((P>>16)&0xFFFF)==X)&&((P&0xFFFF)==Y))||((((P>>16)&0xFFFF)==Y)&&((P&0xFFFF)==X))

/* 
 * Data structure for graph information
 * 
 * n - number of vertices
 * m - number of edges
 * k - size of target set
 * edge_list - a list of m pairs containing the edges
 * adj_list[v] adjacency list for vertex v
 * adj_list_len[v] length of adjacency list for vertex v
 * 
 * p3_vlist[v] list of P3s (pairs u,w) for vertex v
 * p3_vlist_len[v] length of p3_vlist for v
 * 
 * p3_elist[e] list of P3 partners (edge f for P3 e,f) for edge e
 * p3_elist_len[e] length of p3_elist for edge e
 *
 * triangle_elist[e] a list of pairs (f,g) where edges e-f-g form a triangle
 * triangle_elist_len[e] length of triangle_elist for edge e
 *
 */
typedef struct {
  int n, m, k;

  pair_t* edge_list;
  
  int** adj_list;
  int*  adj_list_len;
  
  pair_t** p3_vlist;
  int*     p3_vlist_len;

  int** p3_elist;
  int*  p3_elist_len;

  pair_t** triangle_elist;
  int*     triangle_elist_len;
} graph_data;


int read_graph(graph_data* G, FILE* file);

int read_edges_from_plaintext(int* edgebuf, const int bufsize, FILE* file);

int read_edges_from_compressed(int* edgebuf, const int bufsize, FILE* file);

void free_graph(graph_data* G);


#endif
