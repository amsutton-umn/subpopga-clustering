#ifndef CVD_H
#define CVD_H

#include <stdbool.h>
#include "chromosome.h"
#include "packed_set.h"
#include "graph.h"

bool cvd_feasible(const chromosome* offspr, const graph_data* G);
bool cvd_cluster_graph(const chromosome* offspr, const packed_set* A, const graph_data* G);
bool cvd_repair(chromosome* offspr, const packed_set* x, const packed_set* y, const packed_set* t, const graph_data* G);
void cvd_template(packed_set* z, const chromosome* offspr, const chromosome* p1, const chromosome* p2, const graph_data* G);

#endif
