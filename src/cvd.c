
#include <stdbool.h>
#include <string.h>
#include <assert.h>

/* Uses GNU Linear Programming Kit */
#include <glpk.h> 

#include "chromosome.h"
#include "packed_set.h"
#include "graph.h"
#include "cvd.h"

/* 
 * Determine if G[offspr->V & A] is a cluster graph
 */
bool cvd_cluster_graph(const chromosome* offspr, const packed_set* A, const graph_data* G)
{
  int i, j;
  for (i=0; i<offspr->cached_Vlist_len; i++){
    /* for each v in V*/
    if (ps_read(A,offspr->cached_Vlist[i])){ /* if v is in A */
      int v = offspr->cached_Vlist[i];
      /* check all p3s of v if they are also in V & A */
      for (j=0; j<G->p3_vlist_len[v]; j++){
	int u = src(G->p3_vlist[v][j]);
	int w = snk(G->p3_vlist[v][j]);
	//int u = G->p3_vlist[v][2*j];
	//int w = G->p3_vlist[v][2*j+1];
	if (ps_read(A,u) && ps_read(&offspr->V,u) && ps_read(A,w) && ps_read(&offspr->V,w)) {
	  return false;
	}
      }
    }
  }
  /* no P3s found */
  return true;
}

/* Determine if G[V \ S] is a cluster graph */
bool cvd_feasible(const chromosome* offspr, const graph_data* G)
{
  /* TODO: probably can improve memory management here... */
  packed_set tmp;
  bool is_cluster_graph;
  if (ps_popcount_and(&offspr->S,&offspr->V) > (size_t)G->k) return false;
  ps_copyinit(&tmp,&offspr->V);
  ps_subtract(&tmp,&offspr->S);
  is_cluster_graph = cvd_cluster_graph(offspr,&tmp,G);  
  ps_free(&tmp);
  return is_cluster_graph;
}

bool cvd_repair(chromosome* offspr, const packed_set* x, const packed_set* y, const packed_set* t, const graph_data* G)
{
  int i,j,c,m,l;

  /* LP problem structure */
  glp_prob *lp;

  /* Parameters for simplex solver */
  glp_smcp solver_params;

  /* A is the set of vertices that cannot be deleted */
  static packed_set A;
  
  /* D is the set of vertices to delete */
  static packed_set D;

  /* A set structure for storing vertices */
  static packed_set vertex_store;

  /* Map V -> cluster number*/
  static int* AClusterMap;

  /* Apartition[i] is the set of A vertices in A-cluster i */
  static packed_set* Apartition;

  /* Bpartition[i] is the set of B = V \ A vertices in B-cluster i */
  static packed_set* Bpartition;

  /* Cpartition[i] is the set of B = V \ A vertices in C-partition i */
  static packed_set* Cpartition;

  /* Initialization flag */
  static bool initialized = false;

  /* Initialize static data structures the first time repair is called */
  if (!initialized){
    ps_init(&A,G->n);
    ps_init(&D,G->n);
    ps_init(&vertex_store,G->n);
    AClusterMap = malloc(sizeof(int)*G->n);
    Apartition = malloc(sizeof(packed_set)*G->n);
    Bpartition = malloc(sizeof(packed_set)*G->n);
    Cpartition = malloc(sizeof(packed_set)*G->n);
    for (i=0; i<G->n; i++){
      ps_init(&Apartition[i],G->n);
      ps_init(&Bpartition[i],G->n);
      ps_init(&Cpartition[i],G->n);
    }
    initialized = true;
  }

  /* Determine set A = (x | y | template) - self.S */
  ps_zero(&A);
  ps_union(&A,x,y);
  ps_union(&A,&A,t);
  ps_subtract(&A,&offspr->S);
  ps_intersect(&A,&A,&offspr->V); /* keep only vertices in V */

  /* If G[A] is not a cluster graph, then fail */
  if (!cvd_cluster_graph(offspr,&A,G)) return false;


  /* ********************** */
  /* Find the clusters of A */
  /* ********************** */
  
  /* Cluster 0 is reserved for "not in A" */
  memset(AClusterMap,0,sizeof(int)*G->n);
  for (i=0,c=1; i<offspr->cached_Vlist_len; i++){
    int v = offspr->cached_Vlist[i];
    if (ps_read(&A,v) && AClusterMap[v] == 0){
      AClusterMap[v] = c;
      for (j=0; j<G->adj_list_len[v]; j++){
	int u = G->adj_list[v][j];
	if (ps_read(&A,u)) {
	  /* XXX BUG: this seems to intermittently fail */
	  assert(AClusterMap[u] == 0 || AClusterMap[u] == c);
	  AClusterMap[u] = c;
	}
      }
      c++;
    }
  }
  l = c-1;


  /* Zero out structures */
  for (i=0; i<G->n; i++){
    ps_zero(&Apartition[i]);
    ps_zero(&Bpartition[i]);
    ps_zero(&Cpartition[i]);
  }
  ps_zero(&D);
  
  /* Compute the A partitions */
  for (i=0; i<offspr->cached_Vlist_len; i++){
    int v = offspr->cached_Vlist[i];
    if (ps_read(&A,v)) ps_store(&Apartition[AClusterMap[v]],v);
  }

  /* ******************************* */
  /* Classify the remaining vertices */
  /* ******************************* */
  
  /* For each u in V \ A */
  for (i=0; i<offspr->cached_Vlist_len; i++){
    int u = offspr->cached_Vlist[i];        
    if (!ps_read(&A,u)){
      int adjacentClusters =0;
      int c=-1;
      
      /* Store the neighbors of u that are in V */
      ps_zero(&vertex_store);
      for (j=0; j<G->adj_list_len[u]; j++){
	int w = G->adj_list[u][j];
	if (ps_read(&offspr->V,w)){
	  ps_store(&vertex_store,w);
	}
      }
      /* See how many clusters of A u is adjacent to */
      for (j=1; j<=l; j++){
	if (ps_popcount_and(&vertex_store,&Apartition[j])){
	  adjacentClusters++;
	  c=j;
	}
      }
      if (adjacentClusters == 0){
	/* u is not adjacent to any clusters of A, store in C0 */
	ps_store(&Cpartition[0],u);
      }
      else if (adjacentClusters == 1){
	assert(c > 0);	
	/* check if Apartition[c] is a subset of stored neighbors of u */
	if (ps_compare(&Apartition[c],&vertex_store) >= 0){
	  /* if so, store u into Cpartition[c] */
	  ps_store(&Cpartition[c],u);
	}
	else {
	  /* u is not adjacent to all members of the cluster: delete */
	  ps_store(&D,u);
	}
      }
      else {
	/* u straddles more than one cluster: delete */
	ps_store(&D,u);
      }
    }
  }
                        
  /* Now find the clusters of B */
  ps_zero(&vertex_store);	  
  for (i=0,c=1; i<offspr->cached_Vlist_len; i++){
    int v = offspr->cached_Vlist[i];
    if (!ps_read(&A,v) && !ps_read(&D,v) && !ps_read(&vertex_store,v)){
      ps_store(&Bpartition[c],v);
      ps_store(&vertex_store,v);
      for (j=0; j<G->adj_list_len[v]; j++){
	int u = G->adj_list[v][j];
	if (!ps_read(&A,u) && !ps_read(&D,u)  && !ps_read(&vertex_store,u) && ps_read(&offspr->V,u)) {
	  ps_store(&Bpartition[c],u);
	  ps_store(&vertex_store,u);
	}
      }
      c++;
    }
  }
  m = c-1;

  /* ******************************** */
  /* If l,m > 0 solve the LP instance */
  /* ******************************** */
  if (l*m > 0){
    packed_set tmp;
    int* row_idx, *col_idx;
    double* constr_mat;
    int el;
    
    ps_init(&tmp,G->n);
    for (i=1; i<=m; i++){
      for (j=0; j<l+1; j++){
	ps_copy(&tmp,&Bpartition[i]);
	ps_subtract(&tmp,&Cpartition[j]);
      }
    }

    /* ************************** */
    /* Create LP problem instance */
    /* ************************** */
    glp_init_smcp(&solver_params);
    solver_params.msg_lev = GLP_MSG_ERR;
    solver_params.meth = GLP_PRIMAL;
    /* solver_params.it_lim: iteration limit (default INT_MAX) */
    /* solver_params.tm_lim: time limit (ms) (default INT_MAX) */
    /* solver_params.presolve:  LP presolver (default GLP_OFF) */

    lp = glp_create_prob();
    glp_set_obj_dir(lp,GLP_MIN);
    glp_add_cols(lp,m*(l+1));
    glp_add_rows(lp,m+(l+1)); 

    /* there are exactly 2lm+m nonzero entries in constraint matrix */
    row_idx = calloc(1+2*l*m+m,sizeof(int));
    col_idx = calloc(1+2*l*m+m,sizeof(int));
    constr_mat = calloc(1+2*l*m+m,sizeof(double));


    /* ***************************** */
    /* Compute the constraint matrix */
    /* ***************************** */    
    el=1;
    /* First m constraints: each cluster Bi is assigned exactly one class */
    for (i=1; i<=m; i++){
      for (j=0; j<=l; j++){      
	int idx = (i-1)*(l+1) + (j+1);
	ps_copy(&tmp,&Bpartition[i]);
	ps_subtract(&tmp,&Cpartition[j]);
	glp_set_obj_coef(lp,idx,(double)ps_popcount(&tmp));	
	glp_set_col_bnds(lp,idx,GLP_DB,0.0,1.0);
	row_idx[el] = i;
	col_idx[el] = idx;
	constr_mat[el] = 1.0;
	el++;
      }
      glp_set_row_bnds(lp,i,GLP_FX,1.0,1.0);
    }
    /* Next l constraints: class Cj (j > 0) has at most one cluster assigned to it */
    for (j=1; j<=l; j++){
      for (i=1; i<=m; i++){
	int idx = (i-1)*(l+1) + (j+1);
	row_idx[el] = m+j;
	col_idx[el] = idx;
	constr_mat[el] = 1.0;
	el++;
      }
      glp_set_row_bnds(lp,m+j,GLP_UP,0,1.0);
    }
    glp_load_matrix(lp,el-1,&row_idx[0],&col_idx[0],&constr_mat[0]);

    /* ************ */
    /* Solve the LP */
    /* ************ */    
    glp_simplex(lp,&solver_params);

    /* ******************************************* */
    /* Retrieve the variables and translate to set */
    /* ******************************************* */
    for (i=1; i<=m; i++){
      for (j=0; j<=l; j++){
	//fprintf(stderr,"x%d%d=%f\n",i,j,glp_get_col_prim(lp,(i-1)*(l+1) + (j+1)));
	if (glp_get_col_prim(lp,(i-1)*(l+1) + (j+1)) > 0.99){
	  ps_copy(&tmp,&Bpartition[i]);
	  ps_subtract(&tmp,&Cpartition[j]);
	  ps_union(&D,&D,&tmp);
	}
      }
    }   
    ps_union(&offspr->S,&offspr->S,&D);


    /* ************************* */
    /* Clean up after the solver */
    /* ************************* */
    ps_free(&tmp);
    glp_delete_prob(lp);
    free(row_idx);
    free(col_idx);
    free(constr_mat);
  }
    
  return true;
}

/*
 * Create a template parent z for CVD
 */
void cvd_template(packed_set* z, const chromosome* offspr, const chromosome* p1, const chromosome* p2, const graph_data* G)
{
  int i,j;
  assert(p1->S.capacity == p2->S.capacity && p1->S.capacity == z->capacity);
  ps_zero(z);

  /* while there are P3s, take a P3 out */
  for (i=0; i<offspr->cached_Vlist_len; i++){
    int v = offspr->cached_Vlist[i];
    if (!ps_read(&p1->S,v) && !ps_read(&p2->S,v) && !ps_read(z,v)){
      for (j=0; j<G->p3_vlist_len[v]; j++){
	int u = src(G->p3_vlist[v][j]);
	int w = snk(G->p3_vlist[v][j]);

	if (ps_read(&offspr->V,u) && ps_read(&offspr->V,w)){
	  if (!ps_read(&p1->S,u) && !ps_read(&p2->S,u) && !ps_read(z,u) &&
	      !ps_read(&p1->S,w) && !ps_read(&p2->S,w) && !ps_read(z,w) ){
	    ps_store(z,u);
	    ps_store(z,v);
	    ps_store(z,w);
	    break;
	  }
	}
      }
    }
  }
}
