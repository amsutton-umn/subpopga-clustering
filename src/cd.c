
#include <stdbool.h>
#include <string.h>
#include <assert.h>

#include "chromosome.h"
#include "packed_set.h"
#include "graph.h"
#include "cd.h"

/* Is edge in G[V] - S ? 
 * 
 * not in S
 * both source and sink in V
 */
bool edge_in_graph(int edge, const chromosome* offspr, const graph_data* G)
{
  return !ps_read(&offspr->S,edge) &&	\
    ps_read(&offspr->V,src(G->edge_list[edge])) &&	\
    ps_read(&offspr->V,snk(G->edge_list[edge])); 
}

/* Determine if G[V] - S is a cluster graph */
bool cd_feasible(const chromosome* offspr, const graph_data* G)
{
  int i, j;
  int count=0;
  for (i=0; i<G->m; i++){
    // Is this edge in S and G[V]?
    if (ps_read(&offspr->S,i) && ps_read(&offspr->V,src(G->edge_list[i])) && ps_read(&offspr->V,snk(G->edge_list[i]))) count++;
    if (count > G->k) return false; // More than k edges in S \cap G[V]
    
    if (edge_in_graph(i,offspr,G)) {

      // go through p3 pairs of edge i and see if any of them are also in the graph
      for (j=0; j<G->p3_elist_len[i]; j++){
	if (edge_in_graph(G->p3_elist[i][j],offspr,G)){
	  return false;
	}	  
      }

      // go through the triangles of edge i and see if only one of them is in the graph
      for (j=0; j<G->triangle_elist_len[i]; j++){
	int e1 = src(G->triangle_elist[i][j]);
	int e2 = snk(G->triangle_elist[i][j]);
	if (edge_in_graph(e1,offspr,G) != edge_in_graph(e2,offspr,G)){
	  return false;
	}
      }
    }
  }
  return true;
}

/* Repair operator for CD */
bool cd_repair(chromosome* offspr, const packed_set* x, const packed_set* y, const packed_set* t, const graph_data* G)
{

  int i, j, len;
  
  /* A is the set of edges that cannot be deleted */
  static packed_set A;

  static int* Alist;
  
  /* D is the set of edges to delete */
  static packed_set D;
  
  /* Initialization flag */
  static bool initialized = false;

  /* Initialize static data structures the first time repair is called */
  if (!initialized){
    ps_init(&A,G->m);
    ps_init(&D,G->m);
    Alist = malloc(ps_capacity(&A)*sizeof(int));
    initialized = true;
  }


  /* If offspr is already feasible, success */
  if (cd_feasible(offspr, G)) return true;
  
  /* Determine set A = (x | y | template) - self.S */
  ps_zero(&A);
  ps_union(&A,x,y);
  ps_union(&A,&A,t);
  ps_subtract(&A,&offspr->S);

  
  /* for each e in A, if there is an f in E(H)-A such that (e,f) is a P3 add f to z */
  ps_zero(&D);
  ps_contents(Alist,&len,&A);
  for (i=0; i<len; i++){
    if (!edge_in_graph(Alist[i],offspr,G)) continue;
    for (j=0; j<G->p3_elist_len[Alist[i]]; j++){
      int f = G->p3_elist[Alist[i]][j];
      if (edge_in_graph(f,offspr,G) && !ps_read(&A,f)) {
	ps_store(&D,f);
      }
    }
  }

  ps_union(&offspr->S,&offspr->S,&D);
    
  return true;
}

/*
 * Create a template parent z for CD
 */
void cd_template(packed_set* z, const chromosome* offspr, const chromosome* p1, const chromosome* p2, const graph_data* G)
{
  int e,f,i;
  assert(p1->S.capacity == p2->S.capacity && p1->S.capacity == z->capacity);
  ps_zero(z);
  /* while there are P3s, take a P3 out */
  for (e=0; e<G->m; e++){
    if (edge_in_graph(e,offspr,G) && !ps_read(&p1->S,e) && !ps_read(&p2->S,e) && !ps_read(z,e)){
      for (i=0; i<G->p3_elist_len[e]; i++){
	f = G->p3_elist[e][i];
	if (edge_in_graph(f,offspr,G) && !ps_read(&p1->S,f) && !ps_read(&p2->S,f) && !ps_read(z,f)){
	  ps_store(z,e);
	  ps_store(z,f);
	  break;
	}
      }      
    }
  }
}

