#ifndef CD_H
#define CD_H

#include <stdbool.h>
#include "chromosome.h"
#include "packed_set.h"
#include "graph.h"

bool cd_feasible(const chromosome* offspr, const graph_data* G);
bool cd_cluster_graph(const chromosome* offspr, const packed_set* A, const graph_data* G);
bool cd_repair(chromosome* offspr, const packed_set* x, const packed_set* y, const packed_set* t, const graph_data* G);
void cd_template(packed_set* z, const chromosome* offspr, const chromosome* p1, const chromosome* p2, const graph_data* G);
#endif
