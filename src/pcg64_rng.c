#include <unistd.h>
#include "pcg64_rng.h"
#include <stdio.h>
#include <math.h>

/* SPDX-License-Identifier: 0BSD OR MIT-0 */

/*
 * pcg64 dxsm by Melissa O'Neill
 */
uint64_t pcg64_random_fast(pcg64_random_t *rng) {
	/* cheap (half-width) multiplier */
	const uint64_t mul = 15750249268501108917ULL;
	/* linear congruential generator */
	pcg_ulong_t state = rng->state;
	rng->state = state * mul + rng->inc;
	/* DXSM (double xor shift multiply) permuted output */
	uint64_t hi = (uint64_t)(state >> 64);
	uint64_t lo = (uint64_t)(state | 1);
	hi ^= hi >> 32;
	hi *= mul;
	hi ^= hi >> 48;
	hi *= lo;
	return (hi);
}

/* sample uniformly from (0,1) */
double pcg64_random_unif(pcg64_random_t *rng)
{
  return pcg64_random_fast(rng)/(double)UINT64_MAX;
}

/* sample from a geometric distribution */
size_t pcg64_random_geom(pcg64_random_t *rng, double p)
{
  return 1 + (int)(log1p(-pcg64_random_unif(rng))/log1p(-p));
}


/* sample uniformly from {0,1,...,bound-1} */
uint64_t pcg64_random_bounded(pcg64_random_t *rng, uint64_t bound)
{
  uint64_t threshold = -bound % bound;
  for (;;) {
    uint64_t r = pcg64_random_fast(rng);
    if (r >= threshold)
      return r % bound;
  }
}

/* Select two distinct elements from {0,1,...,bound-1} */
/* (has no effect if bound < 2) */
void pcg64_random_choose2(pcg64_random_t *rng, uint64_t elements[], uint64_t bound)
{
  if (bound == 2){
    elements[0] = pcg64_random_fast(rng) % 2;
    elements[1] = 1 - elements[0];
  }
  else {
    elements[0] = pcg64_random_bounded(rng,bound);
    do {
      elements[1] = pcg64_random_bounded(rng,bound);
    } while (elements[0] == elements[1]);
  }
}


/* seed rng */
void pcg64_getentropy(pcg64_random_t* rng)
{
  const pcg_ulong_t inc = PCG_INCREMENT;
  getentropy(rng,sizeof(rng));
  /* must ensure rng->inc is odd */
  rng->inc = (rng->inc > 0) ? (rng->inc << 1) | 1 : inc;
  rng->state += rng->inc;
  pcg64_random_fast(rng);
}


