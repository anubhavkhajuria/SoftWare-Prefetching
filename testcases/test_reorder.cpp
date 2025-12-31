#include "../pass/reorder.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/time.h>
#include <vector>

using namespace std;

int main(int argc, char **argv) {
  if (argc < 2) {
    cerr << "Usage: ./test_reorder <graphfile> [algorithm]\n";
    cerr << "Algorithms: rcm, hub, dbg, gorder, all\n";
    return 1;
  }

  string filename = argv[1];
  string algo = (argc >= 3) ? argv[2] : "all";

  ifstream infile(filename);
  if (!infile) {
    cerr << "Cannot open file: " << filename << "\n";
    return 1;
  }

  vector<pair<int, int>> edges;
  int max_node = 0;
  string line;

  while (getline(infile, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    istringstream iss(line);
    int src, dst;
    if (iss >> src >> dst) {
      edges.push_back({src, dst});
      max_node = max(max_node, max(src, dst));
    }
  }
  infile.close();

  int n = max_node + 1;
  int m = (int)edges.size();
  cout << "Nodes: " << n << ", Edges: " << m << "\n";

  vector<int> out_degree(n, 0);
  vector<vector<int>> adj(n);
  for (auto &e : edges) {
    out_degree[e.first]++;
    adj[e.second].push_back(e.first);
  }

  vector<int> IA(n + 1), JA, out_edges(n);
  IA[0] = 0;
  for (int i = 0; i < n; i++) {
    IA[i + 1] = IA[i] + (int)adj[i].size();
    for (int src : adj[i]) {
      JA.push_back(src);
    }
    out_edges[i] = out_degree[i];
  }

  vector<int> perm(n);
  struct timeval start, end;

  auto test_algo = [&](const string &name, auto func) {
    cout << "\nTesting " << name << "...\n";
    gettimeofday(&start, nullptr);
    func(n, IA.data(), JA.data(), perm.data());
    gettimeofday(&end, nullptr);
    double elapsed =
        (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6;
    cout << name << " completed in " << elapsed << "s\n";

    bool valid = true;
    vector<bool> seen(n, false);
    for (int i = 0; i < n; i++) {
      if (perm[i] < 0 || perm[i] >= n || seen[perm[i]]) {
        valid = false;
        break;
      }
      seen[perm[i]] = true;
    }
    cout << "Permutation valid: " << (valid ? "YES" : "NO") << "\n";

    if (valid) {
      vector<int> IA_new(n + 1), JA_new(JA.size()), out_new(n);
      reorder::apply_permutation(n, perm.data(), IA.data(), JA.data(),
                                 out_edges.data(), IA_new.data(), JA_new.data(),
                                 out_new.data());
      cout << "Apply permutation: OK\n";
    }
  };

  if (algo == "all" || algo == "rcm") {
    test_algo("RCM", [](int n, const int *IA, const int *JA, int *p) {
      reorder::rcm(n, IA, JA, p);
    });
  }

  if (algo == "all" || algo == "hub") {
    test_algo("Hub Sort", [](int n, const int *IA, const int *JA, int *p) {
      reorder::hub_sort(n, IA, JA, p);
    });
  }

  if (algo == "all" || algo == "dbg") {
    test_algo("DBG", [](int n, const int *IA, const int *JA, int *p) {
      reorder::dbg(n, IA, JA, p);
    });
  }

  if (algo == "all" || algo == "gorder") {
    test_algo("Gorder", [](int n, const int *IA, const int *JA, int *p) {
      reorder::gorder(n, IA, JA, p);
    });
  }

  cout << "\nAll tests completed.\n";
  return 0;
}
