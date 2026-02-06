#ifndef PACKED_SET_H
#define PACKED_SET_H

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "pcg64_rng.h"

typedef uint64_t word;
#define NBYTES 6lu
#define MODMASK 63lu

typedef struct {
  word* data;
  size_t capacity;
  size_t word_cnt;
} packed_set;

void ps_init(packed_set* s, size_t size);
void ps_free(packed_set* s);
void ps_copyinit(packed_set* dest, const packed_set* src);
void ps_copy(packed_set* dest, const packed_set* src);
void ps_copy_bit(packed_set* dest, const packed_set* src, int i);
void ps_zero(packed_set* s);
size_t ps_capacity(const packed_set* s);
bool ps_read(const packed_set* s, int i);
void ps_store(packed_set* s, int i);
void ps_clear(packed_set* s, int i);
void ps_flip(packed_set* s, int i);
void ps_union(packed_set* s, const packed_set* s1, const packed_set* s2);
void ps_intersect(packed_set* s, const packed_set* s1, const packed_set* s2);
void ps_subtract(packed_set* s, const packed_set* r);
void ps_contents(int* list, int* len, const packed_set* s);
size_t ps_popcount(const packed_set* s);
size_t ps_popcount_and(const packed_set* s1, const packed_set* s2);
int ps_compare(const packed_set* s1, const packed_set* s2);
void ps_randomize(packed_set* s, pcg64_random_t* rng);
void ps_debug(const packed_set* s);
#endif
