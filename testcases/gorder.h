
#ifndef GORDER_H
#define GORDER_H

#include <algorithm>
#include <numeric>
#include <queue>
#include <unordered_set>
#include <vector>

namespace gorder {

inline std::vector<int> compute_hub_sorting(int n, const std::vector<int> &IA,
                                            const std::vector<int> &JA) {
  // Compute in-degrees (for pull-based PageRank, this is what matters)
  std::vector<int> in_degree(n, 0);
  for (int i = 0; i < n; i++) {
    in_degree[i] = IA[i + 1] - IA[i];
  }

  // Compute average degree
  long long total_degree = 0;
  for (int i = 0; i < n; i++) {
    total_degree += in_degree[i];
  }
  int avg_degree = static_cast<int>(total_degree / n);
  int hub_threshold = std::max(avg_degree, 1); // At least 1

  // Separate hubs from cold vertices
  std::vector<int> hubs, cold;
  for (int i = 0; i < n; i++) {
    if (in_degree[i] > hub_threshold) {
      hubs.push_back(i);
    } else {
      cold.push_back(i);
    }
  }

  // Sort hubs by descending degree (most frequently accessed first)
  std::sort(hubs.begin(), hubs.end(),
            [&](int a, int b) { return in_degree[a] > in_degree[b]; });

  // Sort cold vertices by ascending degree
  std::sort(cold.begin(), cold.end(),
            [&](int a, int b) { return in_degree[a] < in_degree[b]; });

  // Build permutation: hubs first, then cold
  std::vector<int> perm;
  perm.reserve(n);
  for (int v : hubs)
    perm.push_back(v);
  for (int v : cold)
    perm.push_back(v);

  return perm;
}

inline std::vector<int> compute_dbg(int n, const std::vector<int> &IA,
                                    const std::vector<int> &JA,
                                    int num_groups = 8) {
  std::vector<int> in_degree(n, 0);
  int max_degree = 0;
  for (int i = 0; i < n; i++) {
    in_degree[i] = IA[i + 1] - IA[i];
    max_degree = std::max(max_degree, in_degree[i]);
  }

  // Create buckets based on degree ranges
  std::vector<std::vector<int>> buckets(num_groups);
  for (int i = 0; i < n; i++) {
    int bucket_idx =
        (max_degree > 0)
            ? std::min(num_groups - 1,
                       (in_degree[i] * num_groups) / (max_degree + 1))
            : 0;
    buckets[bucket_idx].push_back(i);
  }

  // Sort within each bucket by degree (descending for high buckets)
  for (int b = 0; b < num_groups; b++) {
    if (b >= num_groups / 2) {
      // High-degree buckets: sort descending
      std::sort(buckets[b].begin(), buckets[b].end(),
                [&](int a, int bb) { return in_degree[a] > in_degree[bb]; });
    } else {
      // Low-degree buckets: sort ascending
      std::sort(buckets[b].begin(), buckets[b].end(),
                [&](int a, int bb) { return in_degree[a] < in_degree[bb]; });
    }
  }

  // Build permutation: high-degree buckets first
  std::vector<int> perm;
  perm.reserve(n);
  for (int b = num_groups - 1; b >= 0; b--) {
    for (int v : buckets[b]) {
      perm.push_back(v);
    }
  }

  return perm;
}

inline std::vector<int> compute_gorder(int n, const std::vector<int> &IA,
                                       const std::vector<int> &JA,
                                       int window_size = 8) {
  // Build undirected adjacency for scoring
  std::vector<std::unordered_set<int>> adj(n);
  for (int i = 0; i < n; i++) {
    for (int j = IA[i]; j < IA[i + 1]; j++) {
      int neighbor = JA[j];
      adj[i].insert(neighbor);
      adj[neighbor].insert(i);
    }
  }

  // Compute degrees
  std::vector<int> degree(n);
  for (int i = 0; i < n; i++) {
    degree[i] = static_cast<int>(adj[i].size());
  }

  std::vector<bool> ordered(n, false);
  std::vector<int> perm;
  perm.reserve(n);

  std::deque<int> window;

  using Entry = std::tuple<int, int, int>;
  std::priority_queue<Entry> pq;

  // Start with minimum degree vertex
  int start = 0;
  for (int i = 1; i < n; i++) {
    if (degree[i] < degree[start])
      start = i;
  }

  pq.push({0, -degree[start], start});

  while (perm.size() < static_cast<size_t>(n)) {
    while (!pq.empty()) {
      auto [score, neg_deg, v] = pq.top();
      pq.pop();
      if (!ordered[v]) {
        ordered[v] = true;
        perm.push_back(v);

        window.push_back(v);
        if (window.size() > static_cast<size_t>(window_size)) {
          window.pop_front();
        }

        for (int neighbor : adj[v]) {
          if (!ordered[neighbor]) {
            int new_score = 0;
            for (int w : window) {
              if (adj[w].count(neighbor) > 0) {
                new_score++;
              }
            }
            pq.push({new_score, -degree[neighbor], neighbor});
          }
        }
        break;
      }
    }

    if (pq.empty() && perm.size() < static_cast<size_t>(n)) {
      for (int i = 0; i < n; i++) {
        if (!ordered[i]) {
          pq.push({0, -degree[i], i});
          break;
        }
      }
    }
  }

  return perm;
}

inline std::vector<int> compute_hub_gorder(int n, const std::vector<int> &IA,
                                           const std::vector<int> &JA,
                                           int window_size = 8) {
  std::vector<int> in_degree(n, 0);
  for (int i = 0; i < n; i++) {
    in_degree[i] = IA[i + 1] - IA[i];
  }

  long long total = 0;
  for (int i = 0; i < n; i++)
    total += in_degree[i];
  int avg = static_cast<int>(total / n);
  int threshold = std::max(avg * 2, 1); // Higher threshold for hubs

  std::vector<int> hubs, rest;
  for (int i = 0; i < n; i++) {
    if (in_degree[i] > threshold) {
      hubs.push_back(i);
    } else {
      rest.push_back(i);
    }
  }

  // Sort hubs by descending degree
  std::sort(hubs.begin(), hubs.end(),
            [&](int a, int b) { return in_degree[a] > in_degree[b]; });

  // Apply Gorder to the rest (need to remap indices)
  // For simplicity, just sort rest by degree ascending
  std::sort(rest.begin(), rest.end(),
            [&](int a, int b) { return in_degree[a] < in_degree[b]; });

  std::vector<int> perm;
  perm.reserve(n);
  for (int v : hubs)
    perm.push_back(v);
  for (int v : rest)
    perm.push_back(v);

  return perm;
}

inline void
apply_permutation(const std::vector<int> &perm, const std::vector<int> &IA_old,
                  const std::vector<int> &JA_old,
                  const std::vector<int> &out_old, std::vector<int> &IA_new,
                  std::vector<int> &JA_new, std::vector<int> &out_new) {
  int n = static_cast<int>(perm.size());

  std::vector<int> inv_perm(n);
  for (int new_id = 0; new_id < n; new_id++) {
    inv_perm[perm[new_id]] = new_id;
  }

  out_new.resize(n);
  for (int old_id = 0; old_id < n; old_id++) {
    out_new[inv_perm[old_id]] = out_old[old_id];
  }

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
    // Sort adjacency list for better spatial locality during traversal
    std::sort(JA_new.begin() + IA_new[new_i],
              JA_new.begin() + IA_new[new_i + 1]);
  }
}

// Backwards compatibility alias
inline void apply_gorder_permutation(const std::vector<int> &perm,
                                     const std::vector<int> &IA_old,
                                     const std::vector<int> &JA_old,
                                     const std::vector<int> &out_old,
                                     std::vector<int> &IA_new,
                                     std::vector<int> &JA_new,
                                     std::vector<int> &out_new) {
  apply_permutation(perm, IA_old, JA_old, out_old, IA_new, JA_new, out_new);
}

} // namespace gorder

#endif // GORDER_H
