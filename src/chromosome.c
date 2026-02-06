#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "chromosome.h"
#include "packed_set.h"

/* Initialize a chromosome with solution length len and n vertices */
void chromosome_init(chromosome* chr, size_t len, size_t n)
{
  ps_init(&chr->S,len);
  ps_init(&chr->V,n);
  chr->cached_Vlist = malloc(n*sizeof(int));
  chr->cached_Vlist_len = 0;
}

/* Set this chromosome's subgraph to contain a single vertex */
void chromosome_seed(chromosome* chr, int vertex)
{
  ps_zero(&chr->V);
  ps_store(&chr->V,vertex);
  chr->cached_Vlist[0] = vertex;
  chr->cached_Vlist_len = 1;  
}

/* Free the memory associated with this chromosome */
void chromosome_free(chromosome*  chr)
{
  ps_free(&chr->S);
  ps_free(&chr->V);
  free(chr->cached_Vlist);
}

/* Copy a chromosome, destroying contents of chr */
void chromosome_copy(chromosome* chr, const chromosome* copy)
{
  /* Copy sets */
  ps_copy(&chr->S,&copy->S);
  ps_copy(&chr->V,&copy->V);
  
  /* Update Vlist cache */
  chromosome_update_cache(chr);
}

/* Merge the vertex sets of chr1 and chr2 into dest, destroying contents of dest */
void chromosome_vmerge(chromosome* dest, chromosome* chr1, chromosome* chr2)
{
  /* Merge vertex sets */
  ps_union(&dest->V,&chr1->V,&chr2->V);

  /* Update Vlist cache */
  chromosome_update_cache(dest);
}

void chromosome_update_cache(chromosome* chr)
{
  ps_contents(chr->cached_Vlist,&chr->cached_Vlist_len,&chr->V);  
}

void chromosome_debug(chromosome* chr)
{
  int i, len;
  int* A = malloc(ps_capacity(&chr->S)*sizeof(int));

  fprintf(stderr,"--- CHROMOSOME ---\nS: ");
  ps_contents(A,&len,&chr->S);
  for (i=0; i<len; i++) fprintf(stderr,"%d ",A[i]);
  fprintf(stderr,"\nV: ");
  ps_contents(A,&len,&chr->V);
  for (i=0; i<len; i++) fprintf(stderr,"%d, ",A[i]);
  fprintf(stderr,"\n");

  fprintf(stderr,"(cached V)\nC: ");
  for (i=0; i<chr->cached_Vlist_len; i++) fprintf(stderr,"%d, ",chr->cached_Vlist[i]);
  fprintf(stderr,"\n");

  free(A);
}
