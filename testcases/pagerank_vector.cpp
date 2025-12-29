// Compile:  g++ pagerank_vector.cpp
// Run:      ./a.out graph.txt 0 400000
// ============================

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <vector>
using namespace std;

// ------------------- Globals -------------------
vector<int> IA, JA, out_edges;
vector<float> page_rank_prev, page_rank;
int nodes = 0;
int max_nodes = 400000;
float threshold = 1e-9, damping_factor = 0.85;
int start_node = 0;
int iterations = 0;

// ------------------- PageRank -------------------
void compute_page_rank() {
  float err;
  const int maxIter = 200;

  do {
    err = 0.0f;

    // --- dangling mass ---
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
        new_rank += (page_rank_prev[src] / out_edges[src]) * damping_factor;
      }

      err = max(err, fabs(new_rank - page_rank_prev[i]));
      page_rank[i] = new_rank;
    }

    for (int i = start_node; i < nodes; i++)
      page_rank_prev[i] = page_rank[i];

    iterations++;
  } while (err > threshold && iterations < maxIter);
}

// ------------------- RCM -------------------
vector<int> reverse_cuthill_mckee(int n, const vector<int> &IA,
                                  const vector<int> &JA) {
  vector<vector<int>> adj(n);
  for (int i = 0; i < n; i++)
    for (int j = IA[i]; j < IA[i + 1]; j++) {
      int src = JA[j];
      adj[i].push_back(src);
      adj[src].push_back(i);
    }

  vector<int> degree(n);
  for (int i = 0; i < n; i++)
    degree[i] = (int)adj[i].size();

  vector<char> visited(n, 0);
  vector<int> order;
  order.reserve(n);

  vector<int> nodes_sorted(n);
  iota(nodes_sorted.begin(), nodes_sorted.end(), 0);
  sort(nodes_sorted.begin(), nodes_sorted.end(),
       [&](int a, int b) { return degree[a] < degree[b]; });

  deque<int> q;
  for (int v : nodes_sorted) {
    if (visited[v])
      continue;
    q.push_back(v);
    visited[v] = 1;
    while (!q.empty()) {
      int u = q.front();
      q.pop_front();
      order.push_back(u);

      vector<int> neighs = adj[u];
      sort(neighs.begin(), neighs.end(),
           [&](int a, int b) { return degree[a] < degree[b]; });

      for (int vv : neighs)
        if (!visited[vv]) {
          visited[vv] = 1;
          q.push_back(vv);
        }
    }
  }

  reverse(order.begin(), order.end());
  return order;
}

// ------------------- Apply permutation -------------------
void apply_permutation_rcm(const vector<int> &perm, const vector<int> &IA_old,
                           const vector<int> &JA_old,
                           const vector<int> &out_old, vector<int> &IA_new,
                           vector<int> &JA_new, vector<int> &out_new) {
  int n = perm.size();
  vector<int> inv_perm(n);
  for (int new_id = 0; new_id < n; new_id++)
    inv_perm[perm[new_id]] = new_id;

  out_new.resize(n);
  for (int old = 0; old < n; old++)
    out_new[inv_perm[old]] = out_old[old];

  IA_new.assign(n + 1, 0);
  for (int new_i = 0; new_i < n; new_i++) {
    int old_i = perm[new_i];
    IA_new[new_i + 1] = IA_new[new_i] + (IA_old[old_i + 1] - IA_old[old_i]);
  }

  JA_new.resize(JA_old.size());
  for (int new_i = 0; new_i < n; new_i++) {
    int old_i = perm[new_i];
    int pos = IA_new[new_i];
    for (int j = IA_old[old_i]; j < IA_old[old_i + 1]; j++) {
      int old_src = JA_old[j];
      JA_new[pos++] = inv_perm[old_src];
    }
  }
}

// ------------------- Main -------------------
int main(int argc, char **argv) {
  if (argc != 4) {
    cerr << "Usage: ./pagerank <graphfile> <start_node> <max_nodes>\n";
    return 1;
  }

  string filename = argv[1];
  start_node = atoi(argv[2]);
  max_nodes = atoi(argv[3]);

  ifstream file(filename);
  if (!file.is_open()) {
    cerr << "Error opening file: " << filename << endl;
    return 1;
  }

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
    edges.emplace_back(d, s); // store as (dst, src)
    max_node_id = max({max_node_id, s, d});
  }
  file.close();

  nodes = max_node_id + 1;
  int edge_count = (int)edges.size();
  cerr << "[read] nodes=" << nodes << " edges=" << edge_count << endl;

  vector<int> counts(nodes, 0);
  out_edges.assign(nodes, 0);

  for (auto &e : edges) {
    int dst = e.first, src = e.second;
    counts[dst]++;
    out_edges[src]++;
  }

  IA.assign(nodes + 1, 0);
  for (int i = 0; i < nodes; i++)
    IA[i + 1] = IA[i] + counts[i];

  JA.resize(edge_count);
  vector<int> fillPos = IA;
  for (auto &e : edges) {
    int dst = e.first, src = e.second;
    JA[fillPos[dst]++] = src;
  }

  page_rank_prev.assign(nodes, 1.0f / nodes);
  page_rank.assign(nodes, 0.0f);

  // --- RCM reorder ---
  cerr << "[info] Applying RCM reordering...\n";
  vector<int> perm = reverse_cuthill_mckee(nodes, IA, JA);

  vector<int> IA_new, JA_new, out_new;
  apply_permutation_rcm(perm, IA, JA, out_edges, IA_new, JA_new, out_new);

  IA.swap(IA_new);
  JA.swap(JA_new);
  out_edges.swap(out_new);

  vector<float> page_rank_prev_new(nodes);
  for (int new_i = 0; new_i < nodes; new_i++)
    page_rank_prev_new[new_i] = page_rank_prev[perm[new_i]];
  page_rank_prev.swap(page_rank_prev_new);

  struct timeval startwtime, endwtime;
  gettimeofday(&startwtime, NULL);
  compute_page_rank();
  gettimeofday(&endwtime, NULL);

  double time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 +
                         endwtime.tv_sec - startwtime.tv_sec);

  cout << "Time : " << time << ", Iterations : " << iterations << "\n";

  double total_rank = 0;
  for (float v : page_rank)
    total_rank += v;
  cout << "Sum of ranks = " << total_rank << endl;
  return 0;
}
