#ifndef REORDER_H
#define REORDER_H

#include <algorithm>
#include <utility>
#include <vector>

namespace reorder {

inline void rcm(int n, const int *IA, const int *JA, int *perm) {
  std::vector<int> degree(n);
  for (int i = 0; i < n; i++) {
    degree[i] = IA[i + 1] - IA[i];
  }

  std::vector<bool> visited(n, false);
  std::vector<int> order;
  order.reserve(n);

  std::vector<std::pair<int, int>> nodes_by_degree(n);
  for (int i = 0; i < n; i++) {
    nodes_by_degree[i] = {degree[i], i};
  }
  std::sort(nodes_by_degree.begin(), nodes_by_degree.end());

  std::vector<int> queue(n);

  for (auto &p : nodes_by_degree) {
    int start = p.second;
    if (visited[start])
      continue;

    int front = 0, back = 0;
    queue[back++] = start;
    visited[start] = true;

    while (front < back) {
      int u = queue[front++];
      order.push_back(u);

      std::vector<std::pair<int, int>> neighbors;
      for (int j = IA[u]; j < IA[u + 1]; j++) {
        int v = JA[j];
        if (v >= 0 && v < n && !visited[v]) {
          neighbors.push_back({degree[v], v});
        }
      }
      std::sort(neighbors.begin(), neighbors.end());

      for (auto &nb : neighbors) {
        int v = nb.second;
        if (!visited[v]) {
          visited[v] = true;
          queue[back++] = v;
        }
      }
    }
  }

  int sz = (int)order.size();
  for (int i = 0; i < sz; i++) {
    perm[i] = order[sz - 1 - i];
  }
  for (int i = sz; i < n; i++) {
    perm[i] = i;
  }
}

inline void hub_sort(int n, const int *IA, const int *JA, int *perm) {
  std::vector<std::pair<int, int>> nodes(n);
  long long total = 0;

  for (int i = 0; i < n; i++) {
    int deg = IA[i + 1] - IA[i];
    nodes[i] = {deg, i};
    total += deg;
  }

  int avg = (n > 0) ? (int)(total / n) : 0;
  int threshold = std::max(avg, 1);

  std::vector<int> hubs, cold;
  for (int i = 0; i < n; i++) {
    if (nodes[i].first > threshold) {
      hubs.push_back(i);
    } else {
      cold.push_back(i);
    }
  }

  std::sort(hubs.begin(), hubs.end(), [&](int a, int b) {
    return (IA[a + 1] - IA[a]) > (IA[b + 1] - IA[b]);
  });

  std::sort(cold.begin(), cold.end(), [&](int a, int b) {
    return (IA[a + 1] - IA[a]) < (IA[b + 1] - IA[b]);
  });

  int idx = 0;
  for (int v : hubs)
    perm[idx++] = v;
  for (int v : cold)
    perm[idx++] = v;
}

inline void dbg(int n, const int *IA, const int *JA, int *perm,
                int num_buckets = 8) {
  int max_deg = 0;
  for (int i = 0; i < n; i++) {
    int deg = IA[i + 1] - IA[i];
    max_deg = std::max(max_deg, deg);
  }

  std::vector<std::vector<int>> buckets(num_buckets);
  for (int i = 0; i < n; i++) {
    int deg = IA[i + 1] - IA[i];
    int bucket = (max_deg > 0) ? std::min(num_buckets - 1,
                                          (deg * num_buckets) / (max_deg + 1))
                               : 0;
    buckets[bucket].push_back(i);
  }

  for (int b = 0; b < num_buckets; b++) {
    std::sort(buckets[b].begin(), buckets[b].end(), [&](int a, int bb) {
      return (IA[a + 1] - IA[a]) > (IA[bb + 1] - IA[bb]);
    });
  }

  int idx = 0;
  for (int b = num_buckets - 1; b >= 0; b--) {
    for (int v : buckets[b]) {
      perm[idx++] = v;
    }
  }
}

inline void gorder(int n, const int *IA, const int *JA, int *perm,
                   int window = 5) {
  std::vector<int> degree(n);
  for (int i = 0; i < n; i++) {
    degree[i] = IA[i + 1] - IA[i];
  }

  std::vector<std::pair<int, int>> nodes(n);
  for (int i = 0; i < n; i++) {
    nodes[i] = {degree[i], i};
  }
  std::sort(nodes.begin(), nodes.end(),
            [](auto &a, auto &b) { return a.first > b.first; });

  std::vector<bool> placed(n, false);
  std::vector<int> order;
  order.reserve(n);

  for (auto &nd : nodes) {
    int v = nd.second;
    if (placed[v])
      continue;

    order.push_back(v);
    placed[v] = true;

    std::vector<std::pair<int, int>> neighbors;
    for (int j = IA[v]; j < IA[v + 1]; j++) {
      int u = JA[j];
      if (u >= 0 && u < n && !placed[u]) {
        neighbors.push_back({degree[u], u});
      }
    }
    std::sort(neighbors.begin(), neighbors.end(),
              [](auto &a, auto &b) { return a.first > b.first; });

    int added = 0;
    for (auto &nb : neighbors) {
      if (added >= window)
        break;
      int u = nb.second;
      if (!placed[u]) {
        order.push_back(u);
        placed[u] = true;
        added++;
      }
    }
  }

  for (int i = 0; i < n; i++) {
    if (!placed[i]) {
      order.push_back(i);
    }
  }

  for (int i = 0; i < n; i++) {
    perm[i] = order[i];
  }
}

inline void apply_permutation(int n, const int *perm, const int *IA_old,
                              const int *JA_old, const int *out_old,
                              int *IA_new, int *JA_new, int *out_new) {
  std::vector<int> inv_perm(n);
  for (int i = 0; i < n; i++) {
    inv_perm[perm[i]] = i;
  }

  for (int old_id = 0; old_id < n; old_id++) {
    out_new[inv_perm[old_id]] = out_old[old_id];
  }

  IA_new[0] = 0;
  for (int new_i = 0; new_i < n; new_i++) {
    int old_i = perm[new_i];
    int valid_count = 0;
    for (int j = IA_old[old_i]; j < IA_old[old_i + 1]; j++) {
      int old_src = JA_old[j];
      if (old_src >= 0 && old_src < n)
        valid_count++;
    }
    IA_new[new_i + 1] = IA_new[new_i] + valid_count;
  }

  for (int new_i = 0; new_i < n; new_i++) {
    int old_i = perm[new_i];
    int pos = IA_new[new_i];
    for (int j = IA_old[old_i]; j < IA_old[old_i + 1]; j++) {
      int old_src = JA_old[j];
      if (old_src >= 0 && old_src < n) {
        JA_new[pos++] = inv_perm[old_src];
      }
    }
  }
}

} // namespace reorder

#endif
