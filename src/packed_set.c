#include <stdio.h>

#include <stddef.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

#include "packed_set.h"
#include "pcg64_rng.h"

/* Allocate memory to contain n elements */
void ps_init(packed_set* s, size_t n)
{
  assert(n > 0);
  s->word_cnt = ((n-1) >> NBYTES) + 1; 
  s->data = calloc(s->word_cnt,sizeof(word));
  s->capacity = n;
}

/* Free memory */
void ps_free(packed_set* s)
{
  free(s->data);
  s->word_cnt = s->capacity = 0;  
}

/* Copy src to test, initializing src */
void ps_copyinit(packed_set* dest, const packed_set* src)
{
  ps_init(dest,src->capacity);
  memcpy(dest->data,src->data,sizeof(word)*src->word_cnt);
}

/* Copy src to dest, assuming dest is already initialized */
void ps_copy(packed_set* dest, const packed_set* src)
{
  assert(src->capacity == dest->capacity);
  memcpy(dest->data,src->data,sizeof(word)*src->word_cnt);
}

/* Copy one bit in src to dest, assuming dest is initialized */
void ps_copy_bit(packed_set* dest, const packed_set* src, int i)
{
  assert(src->capacity == dest->capacity);
  word mask = (1lu << (i&MODMASK));
  dest->data[i>>NBYTES] |= (src->data[i>>NBYTES]&mask);
}

/* Remove all elements */
void ps_zero(packed_set* s)
{
  size_t i;
  for (i=0; i<s->word_cnt; i++) s->data[i] = 0;
}

/* Return the capacity */
size_t ps_capacity(const packed_set* s)
{
  return s->capacity;
}

/* Read the i-th bit / check if i is in the set */
bool ps_read(const packed_set* s, int i)
{
  return (s->data[i >> NBYTES] >> (i&MODMASK)) & 1lu;
}

/* Set the i-th bit / store i in the set: precondition i<s->capacity */
void ps_store(packed_set* s, int i)
{
  s->data[i >> NBYTES] |= 1lu << (i&MODMASK);
}

/* Clear the i-th bit / remove i from the set */
void ps_clear(packed_set* s, int i)
{
  s->data[i >> NBYTES] &= ~(1lu<<(i&MODMASK));
}

/* Flip the i-th bit */
void ps_flip(packed_set* s, int i)
{
  s->data[i >> NBYTES] ^= 1lu<<(i&MODMASK);
}

/* Union the sets s1 and s2, storing in s (s == s1 is allowed) */
void ps_union(packed_set* s, const packed_set* s1, const packed_set* s2)
{
  size_t i;
  assert(s->capacity == s1->capacity && s1->capacity == s2->capacity);
  for (i=0; i<s->word_cnt; i++) s->data[i] = s1->data[i] | s2->data[i];
}

/* Intersect the sets s1 and s2, storing in s (s == s1 is allowed) */
void ps_intersect(packed_set* s, const packed_set* s1, const packed_set* s2)
{
  size_t i;
  assert(s->capacity == s1->capacity && s1->capacity == s2->capacity);
  for (i=0; i<s->word_cnt; i++) s->data[i] = s1->data[i] & s2->data[i]; 
}

/* Remove elements of r from s */
void ps_subtract(packed_set* s, const packed_set* r)
{
  size_t i;
  assert(s->capacity == r->capacity);
  for (i=0; i<s->word_cnt; i++) s->data[i] &= ~r->data[i]; 
}

/* number of bits set */
size_t ps_popcount(const packed_set* s)
{
  size_t i, sum;
  for (i=0,sum=0; i<s->word_cnt; i++) sum+=__builtin_popcountl(s->data[i]);
  return sum;
}

/* number of bits set in the logical an of s1 and s2 */
size_t ps_popcount_and(const packed_set* s1, const packed_set* s2)
{
  size_t i, sum;
  assert(s1->capacity == s2->capacity);
  for (i=0,sum=0; i<s1->word_cnt; i++) sum+=__builtin_popcountl(s1->data[i]&s2->data[i]);
  return sum;
}


/* save {i : s_i == 1} to list and store corresponding list length len */
void ps_contents(int* list, int* len, const packed_set* s)
{
  size_t i,j;
  uint64_t x;
  *len=0;
  for (i=0; i<s->word_cnt; i++){
    x = s->data[i];
    j=0;
    while(x){
      if (x&1) list[(*len)++] = (i << NBYTES) + j;
      x >>= 1;
      j++;
    }
  }  
}

/* write the contents of s to standard out */
void ps_debug(const packed_set* s)
{
  size_t i;
  for (i=0; i<s->capacity; i++) if (ps_read(s,i)) fprintf(stderr,"%lu ",i);
  fprintf(stderr,"\n");
}

/*
 * compare sets s1 and s2
 * return value:
 *   0 if s1 == s2
 *   1 if s1 is a proper subset of s2
 *   -1 otherwise
 */
int ps_compare(const packed_set* s1, const packed_set* s2)
{
  size_t i;
  bool eq=true;
  assert(s1->capacity == s2->capacity);  
  for (i=0; i<s1->word_cnt; i++){
    if (s1->data[i] != s2->data[i]){
      eq=false;
      if ((s1->data[i] & s2->data[i]) != s1->data[i]){
	return -1;
      }
    }
  }
  if (eq) return 0;
  return 1;
}

/* randomize s */
void ps_randomize(packed_set* s, pcg64_random_t* rng)
{
  size_t i;
  for (i=0; i<s->word_cnt; i++) s->data[i] = pcg64_random_fast(rng) % UINT64_MAX;
  /* clear bits beyond capacity */
  s->data[s->word_cnt-1] &= ~(UINT64_MAX << (s->capacity&MODMASK));
}
