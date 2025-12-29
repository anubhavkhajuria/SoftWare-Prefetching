
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <sys/time.h>
#include <vector>

#include "gorder.h"

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

void reverse_cuthill_mckee(int n, int *IA_in, int *JA_in, int *perm) {
  vector<int> degree(n);
  vector<vector<int>> adj(n);

  for (int i = 0; i < n; i++) {
    for (int j = IA_in[i]; j < IA_in[i + 1]; j++) {
      int src = JA_in[j];
      // build undirected adjacency (avoid duplicates if graph had both
      // directions)
      adj[i].push_back(src);
      adj[src].push_back(i);
    }
  }

  for (int i = 0; i < n; i++)
    degree[i] = (int)adj[i].size();

  vector<char> visited(n, 0);
  vector<int> order;
  order.reserve(n);

  vector<int> nodes_sorted(n);
  iota(nodes_sorted.begin(), nodes_sorted.end(), 0);
  sort(nodes_sorted.begin(), nodes_sorted.end(),
       [&](int a, int b) { return degree[a] < degree[b]; });

  vector<int> queue(n);

  for (int s_i = 0; s_i < n; ++s_i) {
    int start = nodes_sorted[s_i];
    if (visited[start])
      continue;

    int front = 0, back = 0;
    queue[back++] = start;
    visited[start] = 1;

    while (front < back) {
      int u = queue[front++];
      order.push_back(u);

      sort(adj[u].begin(), adj[u].end(),
           [&](int a, int b) { return degree[a] < degree[b]; });

      for (int v : adj[u]) {
        if (!visited[v]) {
          visited[v] = 1;
          queue[back++] = v;
        }
      }
    }
  }

  // reverse order into perm
  for (int i = 0; i < n; i++)
    perm[i] = order[n - 1 - i];
}

void apply_permutation_rcm(int n, int *perm, int *IA_old, int *JA_old,
                           int *out_old, int *IA_new, int *JA_new,
                           int *out_new) {
  int *inv_perm = new int[n];
  for (int new_id = 0; new_id < n; new_id++)
    inv_perm[perm[new_id]] = new_id;

  for (int old = 0; old < n; old++) {
    int new_idx = inv_perm[old];
    out_new[new_idx] = out_old[old];
  }

  IA_new[0] = 0;
  for (int new_i = 0; new_i < n; new_i++) {
    int old_i = perm[new_i];
    IA_new[new_i + 1] = IA_new[new_i] + (IA_old[old_i + 1] - IA_old[old_i]);
  }

  for (int new_i = 0; new_i < n; new_i++) {
    int old_i = perm[new_i];
    int start_new = IA_new[new_i];
    for (int j = IA_old[old_i]; j < IA_old[old_i + 1]; j++) {
      int old_src = JA_old[j];
      int new_src = inv_perm[old_src];
      JA_new[start_new++] = new_src;
    }
  }

  delete[] inv_perm;
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
    reverse_cuthill_mckee(nodes, IA, JA, perm);

    int *IA_new = new int[nodes + 1];
    int *JA_new = new int[edge_count];
    int *out_new = new int[nodes];
    apply_permutation_rcm(nodes, perm, IA, JA, out_edges, IA_new, JA_new,
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
  } else if (reorder_mode == 2) {
    // Gorder reordering
    cerr << "Applying Gorder reordering...\n";

    // Convert arrays to vectors for gorder
    vector<int> IA_vec(IA, IA + nodes + 1);
    vector<int> JA_vec(JA, JA + edge_count);
    vector<int> out_vec(out_edges, out_edges + nodes);

    vector<int> perm = gorder::compute_gorder(nodes, IA_vec, JA_vec, 5);

    vector<int> IA_new, JA_new, out_new;
    gorder::apply_gorder_permutation(perm, IA_vec, JA_vec, out_vec, IA_new,
                                     JA_new, out_new);

    delete[] IA;
    delete[] JA;
    delete[] out_edges;

    IA = new int[nodes + 1];
    JA = new int[edge_count];
    out_edges = new int[nodes];
    copy(IA_new.begin(), IA_new.end(), IA);
    copy(JA_new.begin(), JA_new.end(), JA);
    copy(out_new.begin(), out_new.end(), out_edges);

    float *page_rank_prev_new = new float[nodes];
    for (int new_i = 0; new_i < nodes; new_i++)
      page_rank_prev_new[new_i] = page_rank_prev[perm[new_i]];
    delete[] page_rank_prev;
    page_rank_prev = page_rank_prev_new;
  } else if (reorder_mode >= 3 && reorder_mode <= 5) {
    // Hub-based reordering methods
    const char *method_names[] = {"",    "",          "", "Hub Sorting",
                                  "DBG", "Hub+Gorder"};
    cerr << "Applying " << method_names[reorder_mode] << " reordering...\n";

    // Convert arrays to vectors
    vector<int> IA_vec(IA, IA + nodes + 1);
    vector<int> JA_vec(JA, JA + edge_count);
    vector<int> out_vec(out_edges, out_edges + nodes);

    vector<int> perm;
    if (reorder_mode == 3) {
      perm = gorder::compute_hub_sorting(nodes, IA_vec, JA_vec);
    } else if (reorder_mode == 4) {
      perm = gorder::compute_dbg(nodes, IA_vec, JA_vec, 8);
    } else {
      perm = gorder::compute_hub_gorder(nodes, IA_vec, JA_vec, 8);
    }

    vector<int> IA_new, JA_new, out_new;
    gorder::apply_permutation(perm, IA_vec, JA_vec, out_vec, IA_new, JA_new,
                              out_new);

    delete[] IA;
    delete[] JA;
    delete[] out_edges;

    IA = new int[nodes + 1];
    JA = new int[edge_count];
    out_edges = new int[nodes];
    copy(IA_new.begin(), IA_new.end(), IA);
    copy(JA_new.begin(), JA_new.end(), JA);
    copy(out_new.begin(), out_new.end(), out_edges);

    float *page_rank_prev_new = new float[nodes];
    for (int new_i = 0; new_i < nodes; new_i++)
      page_rank_prev_new[new_i] = page_rank_prev[perm[new_i]];
    delete[] page_rank_prev;
    page_rank_prev = page_rank_prev_new;
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
