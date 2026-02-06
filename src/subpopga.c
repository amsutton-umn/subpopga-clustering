
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>

#include "pcg64_rng.h"
#include "graph.h"
#include "chromosome.h"
#include "cvd.h"
#include "cd.h"
#include "params.h"

static pcg64_random_t rng;

/* Standard uniform mutation */
void mutate(packed_set* x,double rate)
{
  size_t i = pcg64_random_geom(&rng,rate) - 1;
  while (i < x->capacity){
    ps_flip(x,i);
    i += pcg64_random_geom(&rng,rate);
  }    
}

/* Uniform 3-way crossover */
void crossover(packed_set* x, packed_set* p1, packed_set* p2, packed_set* p3)
{
  size_t i;
  packed_set* parents[3] = {p1,p2,p3};
  assert(x->capacity == p1->capacity && p1->capacity == p2->capacity);
  ps_zero(x);
  for (i =0; i<x->capacity; i++){    
    ps_copy_bit(x,parents[pcg64_random_bounded(&rng,3)],i);
  }
}

/* Select two distinct elements of [0,n) */
void select2(int* a, int* b, int n)
{
  assert(n > 1);
  if (n == 2){
    *a = pcg64_random_fast(&rng) % 2;
    *b = 1 - *a;
  }
  else {
    *a = pcg64_random_bounded(&rng,n);
    do {
      *b = pcg64_random_bounded(&rng,n);
    } while (*a == *b);
  }
}
  
/* Return how many elements in set, otherwise -1 if not feasible */
int calculate(chromosome* chr, graph_data*G, bool (*feasible)(const chromosome*, const graph_data*))
{
  if (!feasible(chr,G)) return -1;
  return ps_popcount(&chr->S);
}


int main(int argc, char** argv)
{
  params prm;
  char* typestr;
  int i;
  size_t t, popsize, setlen;
  graph_data G;
  chromosome** P;
  chromosome offspr;
  packed_set tau;
  FILE* file;
  bool solved=false;
  size_t cutoff = 0;
  
  /* function pointers */
  bool (*feasible)(const chromosome*, const graph_data*);
  bool (*repair)(chromosome*, const packed_set*, const packed_set*, const packed_set*, const graph_data*);
  void (*template)(packed_set*, const chromosome*, const chromosome*, const chromosome*, const graph_data*);


  /* set parameters from command line arguments */
  set_params_from_args(&prm,argc,argv);


  switch(prm.type){
  case CVD:
    feasible = &cvd_feasible;
    repair = &cvd_repair;
    template = &cvd_template;
    typestr="cvd";
    break;
  case CD:
    feasible = &cd_feasible;
    repair = &cd_repair;
    template = &cd_template;
    typestr="cd";
    break;
  default:
    perror(strerror(ENOSYS));
    exit(EXIT_FAILURE);
  }
  
  /* open graph file and read graph */
  file = fopen(prm.input_filename,"r");  
  if (file == NULL) {
    perror(strerror(ENOENT));
    fprintf(stderr, "Unable to open graph file '%s'\n",prm.input_filename);
    exit(EXIT_FAILURE);
  }
  read_graph(&G,file);
  fclose(file);
 
  fprintf(stderr,"Loaded a graph with %d vertices and %d edges\n",G.n,G.m);
  popsize = G.n;
  G.k = prm.k;
  cutoff = prm.cutoff;
  if (prm.type == CVD) setlen = G.n;
  if (prm.type == CD) setlen = G.m;
    
  /* seed random number generator */
  fprintf(stderr,"Seeding random number generator...\n");
  pcg64_getentropy(&rng);

  /* initialize population */
  fprintf(stderr,"Initializing population...\n");
  P = malloc(popsize*sizeof(chromosome*));
  for (i=0; i<G.n; i++){
    P[i] = malloc(sizeof(chromosome));
    chromosome_init(P[i],setlen,G.n);
    chromosome_seed(P[i],i);
    ps_randomize(&P[i]->S,&rng);
  }
  chromosome_init(&offspr,setlen,G.n);
  ps_init(&tau,setlen);

  fprintf(stderr,"Starting run with n=%d, k=%d, popsize=%lu, cutoff=%lu\n",G.n,G.k,popsize,cutoff);

  /* main loop */
  t = 0;
  while( popsize > 1){
    uint64_t parent[2];
    int r;

    /* choose parents */
    pcg64_random_choose2(&rng,parent,popsize);

    /* crossover */
    if (pcg64_random_unif(&rng) < 0.8){

      /* offspring vertex set is union of parent vertex sets */
      chromosome_vmerge(&offspr,P[parent[0]],P[parent[1]]);

      /* compute template parent */
      template(&tau,&offspr,P[parent[0]],P[parent[1]],&G);

      /* 3-way uniform crossover */
      crossover(&offspr.S,&P[parent[0]]->S,&P[parent[1]]->S,&tau);

      /* repair operator */
      repair(&offspr,&P[parent[0]]->S,&P[parent[1]]->S,&tau,&G);

      /* determine feasibility */
      r = calculate(&offspr,&G,feasible);

      if (r >= 0) {
	/* offspring was feasible, it must dominate both parents */
	chromosome_copy(P[parent[0]],&offspr);
	chromosome* tmp = P[popsize - 1];
	P[popsize-1] = P[parent[1]];
	P[parent[1]] = tmp;
	popsize--;
      }
    }
    /* mutation */
    else {
      /* copy and flip each bit of parent 0 to create offspring */
      chromosome_copy(&offspr,P[parent[0]]);     
      mutate(&offspr.S,1.0/setlen);

      /* determine feasibility */
      r = calculate(&offspr,&G,feasible);
      
      if (r >= 0 && r <= calculate(P[parent[0]],&G,feasible)){
	/* offspring was feasible and dominates parent */
	chromosome_copy(P[parent[0]],&offspr);
      }      
    }
    if (++t >= cutoff) break;
    
    if (popsize == 1) {
      solved = true;
      break;
    }
  }

  /* output results */
  fprintf(stdout,"%s,%d,%d,%d,%s,%lu,%d,%lu,%lu\n",prm.input_filename,G.n,G.m,G.k,typestr,t,solved,popsize,cutoff);

  /* save solution if requested */
  if (solved && prm.save_solution){
    FILE* save_chr = fopen(prm.solution_filename,"w");
    if (!save_chr) {
      fprintf(stderr, "ERROR: fopen failed for '%s' (%s)\n", prm.solution_filename, strerror(errno));
    }
    else{
      int* A = malloc(ps_capacity(&P[0]->S)*sizeof(int));
      int len;
      ps_contents(A,&len,&P[0]->S);
      if (prm.type == CVD){
	for (i=0; i<len; i++) fprintf(save_chr,"%d\n",A[i]);
      }
      else if (prm.type == CD){
	for (i=0; i<len; i++) fprintf(save_chr,"%d %d\n",src(G.edge_list[A[i]]),snk(G.edge_list[A[i]]));
      }
      else {
	perror(strerror(ENOSYS));
	exit(EXIT_FAILURE);
      }
      fclose(save_chr);
      fprintf(stderr,"wrote length %d solution to '%s'\n",len, prm.solution_filename);
    }
  }
  
  //fprintf(stderr, "++++++++++ FINAL POPULATION ++++++++++\n");
  if (!solved){
    fprintf(stderr,"Unsolved; final population\n");
    for (i=0; i<(int)popsize; i++){
      chromosome_debug(P[i]);
    }
  }

  ps_free(&tau);
  free_graph(&G);
  

}

