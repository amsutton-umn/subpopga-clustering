#include <stdio.h>
#include <string.h>
#include <argp.h>

#include "params.h"


static char doc[] = "\nSubgraph-Population Genetic Algorithm (subpopga)\n"\
  "------------------------------------------------";

static struct argp_option options[] = {
  {
    .name  = NULL,
    .key   = 'i',
    .arg   = "<file>",
    .flags = 0,
    .doc   = "input file",
    .group = 0
  },
  {
    .name  = NULL,
    .key   = 'k',
    .arg   = "<k>",
    .flags = 0,
    .doc   = "solution size",
    .group = 0
  },
  {
    .name  = NULL,
    .key   = 'c',
    .arg   = "<cutoff>",
    .flags = 0,
    .doc   = "maximum number of generations",
    .group = 0
  },
  {
    .name  = NULL,
    .key   = 't',
    .arg   = "<type>",
    .flags = 0,
    .doc   = "specify the problem type: cvd | cd | cep",
    .group = 0
  },
  {0,0,0,0,"Optional arguments:",1},   
  {
    .name  = 0,
    .key   = 's',
    .arg   = "FILE",
    .flags = 0,
    .doc   = "save solution in file (if found)",
    .group = 1
  },
  { NULL, 'h', 0, OPTION_HIDDEN, NULL, -1 },  
};

void init_params(params* prm)
{
  prm->input_filename = NULL;
  prm->solution_filename = NULL;
  prm->k = 0;
  prm->cutoff = 0;
  prm->type = NONE;
  prm->save_solution = false;
}

/* Parse a single option. */
static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
  params *prm = state->input;

  switch(key) { 
  case 'i':
    prm->input_filename = arg;
    break;
  case 'k':
    prm->k = atoi(arg);
    break;
  case 's':
    prm->solution_filename = arg;
    prm->save_solution = true;
    break;
  case 'c':
    prm->cutoff = atoi(arg);
    break;
  case 't':
    if (strcmp(arg,"cvd") == 0){
      prm->type = CVD;
    }
    else if (strcmp(arg,"cd") == 0) {
      prm->type = CD;      
    }
    else if (strcmp(arg,"cep") == 0) {
      prm->type = CEP;
    }
    else {
      fprintf(state->err_stream,"\n");
      argp_failure(state,0,EINVAL,"ERROR: problem type '%s'",arg);
      argp_state_help(state, state->out_stream, ARGP_HELP_STD_HELP);
      return EINVAL;
    }

    break;

  case 'h':
    argp_state_help(state, state->out_stream, ARGP_HELP_STD_HELP);
    break;
  case ARGP_KEY_ARG:
    argp_usage (state);
    break;
  case ARGP_KEY_END:
    /* Check for missing arguments */
    if (!prm->input_filename) {
      fprintf(state->err_stream,"\n");
      argp_failure(state,0,0,"ERROR: missing input file");
      argp_state_help(state, state->out_stream, ARGP_HELP_STD_HELP);
      return EINVAL;
    }
    if (!prm->k){
      fprintf(state->err_stream,"\n");
      argp_failure(state,0,0,"ERROR: missing k value");
      argp_state_help(state, state->out_stream, ARGP_HELP_STD_HELP);
      return EINVAL;
    }
    if (!prm->cutoff) {
      fprintf(state->err_stream,"\n");
      argp_failure(state,0,0,"ERROR: missing cutoff value");
      argp_state_help(state, state->out_stream, ARGP_HELP_STD_HELP);
      return EINVAL;
    }
    if (!prm->type) {
      fprintf(state->err_stream,"\n");
      argp_failure(state,0,0,"ERROR: missing problem instance type");
      argp_state_help(state, state->out_stream, ARGP_HELP_STD_HELP);
      return EINVAL;
    }
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}


bool set_params_from_args(params* prm, int argc, char** argv)
{
  init_params(prm);
  struct argp argp = { options, parse_opt, 0, doc, 0, 0, 0 };
  argp_parse (&argp, argc, argv, ARGP_NO_ARGS, 0, prm);
  
  return true;
}
