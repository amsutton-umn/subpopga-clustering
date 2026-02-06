#ifndef PCG64_RNG_H
#define PCG64_RNG_H

#include <stdint.h>

#define PCG_INCREMENT    ((pcg_ulong_t)6364136223846793005ULL << 64 |\
                          (pcg_ulong_t)1442695040888963407ULL)

#define pcg_ulong_t      __uint128_t

typedef struct { pcg_ulong_t state, inc; } pcg64_random_t;

uint64_t pcg64_random_fast(pcg64_random_t *rng);
double pcg64_random_unif(pcg64_random_t *rng);
uint64_t pcg64_random_geom(pcg64_random_t *rng, double p);
uint64_t pcg64_random_bounded(pcg64_random_t *rng, uint64_t bound);
void pcg64_random_choose2(pcg64_random_t *rng, uint64_t elements[], uint64_t bound);
void pcg64_getentropy(pcg64_random_t* rng);

#endif
