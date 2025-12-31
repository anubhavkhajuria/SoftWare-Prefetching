
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/time.h>
#include <vector>

#include "../pass/reorder.h"

using namespace std;

int *IA = nullptr;
int *JA = nullptr;
int *out_edges = nullptr;
float *page_rank_prev = nullptr;
float *page_rank = nullptr;

int nodes = 0;
int max_nodes = 400000;
float threshold = 1e-9, damping_factor = 0.85;
int start_node = 0;
int iterations = 0;

// Reorder mode: 0=none, 1=rcm, 2=gorder, 3=hub, 4=dbg, 5=hubgorder
int reorder_mode = 3; // default to hub (best for power-law graphs)

void compute_page_rank() {
  float err;
  const int maxIter = 200;

  do {
    err = 0.0f;

    double dangling_mass = 0.0;
    for (int i = start_node; i < nodes; i++)
      if (out_edges[i] == 0)
        dangling_mass += page_rank_prev[i];
    dangling_mass /= nodes;

    for (int i = start_node; i < nodes; i++) {
      float new_rank = (1.0f - damping_factor) / nodes;
      new_rank += damping_factor * dangling_mass;

      for (int j = IA[i]; j < IA[i + 1]; j++) {
        int src = JA[j];
        // if (out_edges[src] > 0)
        new_rank += (page_rank_prev[src] / out_edges[src]) * damping_factor;
      }

      err = std::max(err, fabs(new_rank - page_rank_prev[i]));
      page_rank[i] = new_rank;
    }

    // copy new ranks
    for (int i = start_node; i < nodes; i++)
      page_rank_prev[i] = page_rank[i];

    iterations++;
  } while (err > threshold && iterations < maxIter);
}

int main(int argc, char **argv) {
  if (argc < 4) {
    cerr << "Usage: ./pagerank <graphfile> <start_node> <max_nodes_guess> "
            "[rcm|gorder|hub|dbg|hubgorder|none]\n";
    return 1;
  }

  string filename = argv[1];
  start_node = atoi(argv[2]);
  max_nodes = atoi(argv[3]);

  // Parse reorder mode (optional, default hub)
  if (argc >= 5) {
    if (strcmp(argv[4], "rcm") == 0) {
      reorder_mode = 1;
    } else if (strcmp(argv[4], "gorder") == 0) {
      reorder_mode = 2;
    } else if (strcmp(argv[4], "hub") == 0) {
      reorder_mode = 3;
    } else if (strcmp(argv[4], "dbg") == 0) {
      reorder_mode = 4;
    } else if (strcmp(argv[4], "hubgorder") == 0) {
      reorder_mode = 5;
    } else if (strcmp(argv[4], "none") == 0) {
      reorder_mode = 0;
    }
  }

  ifstream file(filename);
  if (!file.is_open()) {
    cerr << "Error opening file: " << filename << endl;
    return 1;
  }

  // --- read edges (first pass) ---
  vector<pair<int, int>> edges;
  edges.reserve(1000000);
  string line;
  int s, d;
  int max_node_id = -1;
  while (getline(file, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    istringstream iss(line);
    if (!(iss >> s >> d))
      continue;
    // store as (dst, src) for pull-based PageRank later
    edges.push_back({d, s});
    if (s > max_node_id)
      max_node_id = s;
    if (d > max_node_id)
      max_node_id = d;
  }
  file.close();

  nodes = max_node_id + 1;
  int edge_count = (int)edges.size();
  cout << "[read] nodes=" << nodes << " edges=" << edge_count << endl;

  // count incoming per node and out-degree per node
  vector<int> counts(nodes, 0);
  out_edges = new int[nodes](); // initialize to 0

  for (auto &e : edges) {
    int dst = e.first, src = e.second;
    counts[dst]++;
    out_edges[src]++;
  }

  // allocate IA with CSR convention IA[0] = 0, IA[i+1] = IA[i] + counts[i]
  IA = new int[nodes + 1];
  IA[0] = 0;
  for (int i = 0; i < nodes; ++i)
    IA[i + 1] = IA[i] + counts[i];

  // allocate JA and fill using fillPos = IA (copy)
  JA = new int[edge_count];
  int *fillPos = new int[nodes];
  for (int i = 0; i < nodes; ++i)
    fillPos[i] = IA[i];

  for (auto &e : edges) {
    int dst = e.first, src = e.second;
    JA[fillPos[dst]++] = src;
  }
  delete[] fillPos;

  // --- initialize PageRank arrays ---
  page_rank_prev = new float[nodes];
  page_rank = new float[nodes];
  for (int i = 0; i < nodes; i++) {
    page_rank_prev[i] = 1.0f / nodes;
    page_rank[i] = 0.0f;
  }

  // --- Graph reordering ---
  if (reorder_mode == 1) {
    // RCM reordering
    cerr << "Applying RCM reordering...\n";
    int *perm = new int[nodes];
    reorder::rcm(nodes, IA, JA, perm);

    int *IA_new = new int[nodes + 1];
    int *JA_new = new int[edge_count];
    int *out_new = new int[nodes];
    reorder::apply_permutation(nodes, perm, IA, JA, out_edges, IA_new, JA_new,
                               out_new);

    delete[] IA;
    delete[] JA;
    delete[] out_edges;
    IA = IA_new;
    JA = JA_new;
    out_edges = out_new;

    float *page_rank_prev_new = new float[nodes];
    for (int new_i = 0; new_i < nodes; new_i++)
      page_rank_prev_new[new_i] = page_rank_prev[perm[new_i]];
    delete[] page_rank_prev;
    page_rank_prev = page_rank_prev_new;
    delete[] perm;
  } else if (reorder_mode >= 2 && reorder_mode <= 5) {
    const char *method_names[] = {"",    "",          "Gorder", "Hub Sorting",
                                  "DBG", "Hub+Gorder"};
    cerr << "Applying " << method_names[reorder_mode] << " reordering...\n";

    int *perm = new int[nodes];

    if (reorder_mode == 2) {
      reorder::gorder(nodes, IA, JA, perm, 5);
    } else if (reorder_mode == 3) {
      reorder::hub_sort(nodes, IA, JA, perm);
    } else if (reorder_mode == 4) {
      reorder::dbg(nodes, IA, JA, perm, 8);
    } else {
      reorder::gorder(nodes, IA, JA, perm, 8);
    }

    int *IA_new = new int[nodes + 1];
    int *JA_new = new int[edge_count];
    int *out_new = new int[nodes];

    reorder::apply_permutation(nodes, perm, IA, JA, out_edges, IA_new, JA_new,
                               out_new);

    delete[] IA;
    delete[] JA;
    delete[] out_edges;
    IA = IA_new;
    JA = JA_new;
    out_edges = out_new;

    float *page_rank_prev_new = new float[nodes];
    for (int new_i = 0; new_i < nodes; new_i++)
      page_rank_prev_new[new_i] = page_rank_prev[perm[new_i]];
    delete[] page_rank_prev;
    page_rank_prev = page_rank_prev_new;
    delete[] perm;
  } else {
    cerr << "No reordering applied.\n";
  }

  // --- run PageRank ---
  struct timeval startwtime, endwtime;
  gettimeofday(&startwtime, NULL);
  compute_page_rank();
  gettimeofday(&endwtime, NULL);

  double time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 +
                         endwtime.tv_sec - startwtime.tv_sec);

  cout << "Time : " << time << ", Iterations : " << iterations << "\n";

  double total_rank = 0;
  for (int i = 0; i < nodes; i++)
    total_rank += page_rank[i];
  cout << "Sum of ranks = " << total_rank << endl;

  delete[] IA;
  delete[] JA;
  delete[] out_edges;
  delete[] page_rank;
  delete[] page_rank_prev;

  return 0;
}
