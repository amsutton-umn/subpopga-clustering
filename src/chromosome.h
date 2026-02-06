/*
 * Definitions for chromosome structure
 */

#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include <stdlib.h>
#include <stdbool.h>

#include "packed_set.h"

typedef struct {
  packed_set S;
  packed_set V;
  int* cached_Vlist;
  int  cached_Vlist_len;
  int* cached_Elist;
  int  cached_Elist_len;
} chromosome;

/* Initialize chromosome chr with solution size len and vertex size n */
void chromosome_init(chromosome* chr, size_t len, size_t n);

/* "Seed" a chromosome chr by inserting a single vertex v into its graph */
void chromosome_seed(chromosome* chr, int v);

/* Free memory associated with chromosome */
void chromosome_free(chromosome* chr);

/* Copy a chromosome, destroying contents of dest */
void chromosome_copy(chromosome* dest, const chromosome* src);

/* Merge the vertex sets of p1 and p2 into dest */
void chromosome_vmerge(chromosome* dest, chromosome* p1, chromosome*p2);

/* Update the cached vertex/edge lists (this must be done any time chr->V is changed, which isn't often */
void chromosome_update_cache(chromosome* chr);

void chromosome_debug(chromosome*);

#endif
