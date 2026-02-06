#ifndef PARAMS_H
#define PARAMS_H

#include <stdlib.h>
#include <stdbool.h>
#include <argp.h>


typedef enum {
  NONE, CVD, CD, CEP
} prob_type;

typedef struct {
  char* input_filename;
  char* solution_filename;
  size_t k;
  size_t cutoff;
  prob_type type;
  bool save_solution;
} params;

void init_params(params*);
bool set_params_from_args(params*,int, char**);

#endif
