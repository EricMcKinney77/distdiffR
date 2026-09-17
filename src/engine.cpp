#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

// Fenwick Tree for O(log N) prefix sums
class FenwickTree {
  std::vector<int> tree;

public:
  FenwickTree(int n) : tree(n + 1, 0) {}
  void add(int i, int val) {
    for (; i < (int)tree.size(); i += i & -i)
      tree[i] += val;
  }
  int query(int i) {
    int sum = 0;
    for (; i > 0; i -= i & -i)
      sum += tree[i];
    return sum;
  }
  void clear() { std::fill(tree.begin(), tree.end(), 0); }
};

struct Event {
  double x;
  double y;
  int type; // 0 for data point, 1 for eval point
  int id;
  bool operator<(const Event &other) const {
    if (x != other.x)
      return x < other.x;
    return type < other.type;
  }
};

std::vector<double>
compute_bcdf_offline(const std::vector<std::pair<double, double>> &data,
                     const std::vector<std::pair<double, double>> &eval) {
  int n = data.size();
  int m = eval.size();
  if (n == 0)
    return std::vector<double>(m, 0.0);

  // Coordinate compression: collect all unique y-coordinates to map them to ranks 1..N
  std::vector<double> y_coords;
  y_coords.reserve(n + m);
  for (auto &p : data)
    y_coords.push_back(p.second);
  for (auto &p : eval)
    y_coords.push_back(p.second);
  std::sort(y_coords.begin(), y_coords.end());
  y_coords.erase(std::unique(y_coords.begin(), y_coords.end()), y_coords.end());

  // Sweep-line algorithm: sort events by x-coordinate. 
  // Process data points (type 0) and evaluation points (type 1) in order.
  std::vector<Event> events;
  events.reserve(n + m);
  for (int i = 0; i < n; ++i)
    events.push_back({data[i].first, data[i].second, 0, i});
  for (int i = 0; i < m; ++i)
    events.push_back({eval[i].first, eval[i].second, 1, i});
  std::sort(events.begin(), events.end());

  FenwickTree ft(y_coords.size());
  std::vector<double> results(m, 0.0);

  for (const auto &e : events) {
    if (e.type == 0) {
      // For data points, increment the count at the rank of its y-coordinate
      int rank = std::lower_bound(y_coords.begin(), y_coords.end(), e.y) -
                 y_coords.begin() + 1;
      ft.add(rank, 1);
    } else {
      // For evaluation points, query the Fenwick Tree for the number of data points 
      // with y-coordinate <= current evaluation y
      int rank = std::lower_bound(y_coords.begin(), y_coords.end(), e.y) -
                 y_coords.begin() + 1;
      results[e.id] = (double)ft.query(rank) / n;
    }
  }
  return results;
}

double calculate_psi(const std::vector<std::pair<double, double>> &s1,
                     const std::vector<std::pair<double, double>> &s2,
                     int stat_type) {
  int n1 = s1.size();
  int n2 = s2.size();

  std::vector<double> b11 = compute_bcdf_offline(s1, s1);
  std::vector<double> b21 = compute_bcdf_offline(s2, s1);
  std::vector<double> b22 = compute_bcdf_offline(s2, s2);
  std::vector<double> b12 = compute_bcdf_offline(s1, s2);

  double sum1 = 0, sum2 = 0;
  // Determine if we use absolute difference (L1) or squared difference (L2)
  bool is_abs = (stat_type < 3);

  for (int i = 0; i < n1; ++i) {
    double diff = b11[i] - b21[i];
    sum1 += is_abs ? std::abs(diff) : diff * diff;
  }
  for (int i = 0; i < n2; ++i) {
    double diff = b22[i] - b12[i];
    sum2 += is_abs ? std::abs(diff) : diff * diff;
  }

  // Select weighting mode:
  // Mode 0: Weighted by n2, n1 respectively (Complementary)
  // Mode 1: Weighted by n1, n2 respectively (Direct)
  // Mode 2: Unweighted sum
  int mode = stat_type % 3;
  if (mode == 0)
    return (n2 * sum1 + n1 * sum2) / (n1 + n2);
  if (mode == 1)
    return (n1 * sum1 + n2 * sum2) / (n1 + n2);
  return sum1 + sum2;
}

double calculate_psi_grouped(const std::vector<std::pair<double, double>> &s1,
                             const std::vector<std::pair<double, double>> &s2,
                             const std::vector<int> &subjects_s1,
                             const std::vector<int> &subjects_s2) {
  int n1 = s1.size();
  int n2 = s2.size();

  int max_subj = 0;
  for (int s : subjects_s1)
    max_subj = std::max(max_subj, s);
  for (int s : subjects_s2)
    max_subj = std::max(max_subj, s);
  if (max_subj == 0)
    return 0.0;

  std::vector<double> b11_avg(n1, 0.0), b21_avg(n1, 0.0), b22_avg(n2, 0.0),
      b12_avg(n2, 0.0);

  // Treat each subject equally by averaging their respective BECDF contributions
  for (int j = 1; j <= max_subj; ++j) {
    std::vector<std::pair<double, double>> subj_s1, subj_s2;
    for (size_t i = 0; i < s1.size(); ++i)
      if (subjects_s1[i] == j)
        subj_s1.push_back(s1[i]);
    for (size_t i = 0; i < s2.size(); ++i)
      if (subjects_s2[i] == j)
        subj_s2.push_back(s2[i]);

    if (!subj_s1.empty()) {
      // Average BECDF of sample 1 vs sample 1 and sample 2 vs sample 1
      std::vector<double> b = compute_bcdf_offline(subj_s1, s1);
      for (int i = 0; i < n1; ++i)
        b11_avg[i] += b[i];
      std::vector<double> b2 = compute_bcdf_offline(subj_s1, s2);
      for (int i = 0; i < n2; ++i)
        b12_avg[i] += b2[i];
    }
    if (!subj_s2.empty()) {
      // Average BECDF of sample 2 vs sample 1 and sample 2 vs sample 2
      std::vector<double> b = compute_bcdf_offline(subj_s2, s1);
      for (int i = 0; i < n1; ++i)
        b21_avg[i] += b[i];
      std::vector<double> b2 = compute_bcdf_offline(subj_s2, s2);
      for (int i = 0; i < n2; ++i)
        b22_avg[i] += b2[i];
    }
  }

  double inv_subj = 1.0 / max_subj;
  double sum1 = 0, sum2 = 0;
  for (int i = 0; i < n1; ++i) {
    double diff = (b11_avg[i] * inv_subj) - (b21_avg[i] * inv_subj);
    sum1 += diff * diff;
  }
  for (int i = 0; i < n2; ++i) {
    double diff = (b22_avg[i] * inv_subj) - (b12_avg[i] * inv_subj);
    sum2 += diff * diff;
  }

  return (n2 * sum1 + n1 * sum2) / (n1 + n2);
}

struct PermutationWorker : public Worker {
  const std::vector<std::pair<double, double>> &pooled_data;
  const std::vector<int> &original_subjects;
  int n1, n2;
  int num_rot, num_shifts;
  int stat_type;
  int seedNum;
  std::vector<std::pair<double, double>> pivots;

  RcppParallel::RVector<double> perm_psi;

  PermutationWorker(const std::vector<std::pair<double, double>> &data,
                    const std::vector<int> &subjects, int n1, int n2,
                    int num_rot, int num_shifts, int stat_type, int seedNum,
                    const std::vector<std::pair<double, double>> &pivots,
                    NumericVector out_psi)
      : pooled_data(data), original_subjects(subjects), n1(n1), n2(n2),
        num_rot(num_rot), num_shifts(num_shifts), stat_type(stat_type),
        seedNum(seedNum), pivots(pivots), perm_psi(out_psi) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<int> subjects = original_subjects;

    for (std::size_t i = begin; i < end; ++i) {
      std::mt19937 g(seedNum + i);
      std::shuffle(subjects.begin(), subjects.end(), g);

      double trial_psi_sum = 0;
      int total_configs = 0;

      for (int r = 0; r < num_rot; ++r) {
        // Apply rotation to data to capture distributional differences across angles
        double angle = 2.0 * M_PI * r / num_rot;
        double cos_r = std::cos(angle);
        double sin_r = std::sin(angle);

        std::vector<std::pair<double, double>> rotated_data;
        rotated_data.reserve(pooled_data.size());
        for (auto &p : pooled_data) {
          rotated_data.push_back({p.first * cos_r - p.second * sin_r,
                                   p.first * sin_r + p.second * cos_r});
        }

        // Define current bounding box for toroidal wrapping
        double min_x = rotated_data[0].first, max_x = rotated_data[0].first;
        double min_y = rotated_data[0].second, max_y = rotated_data[0].second;
        for (auto &p : rotated_data) {
          min_x = std::min(min_x, p.first);
          max_x = std::max(max_x, p.first);
          min_y = std::min(min_y, p.second);
          max_y = std::max(max_y, p.second);
        }
        double rot_width = max_x - min_x;
        double rot_height = max_y - min_y;

        for (int s = 0; s < (int)pivots.size(); ++s) {
          // Toroidal shift: wrap data points relative to a pivot point
          auto pivot = pivots[s];
          double px = pivot.first * cos_r - pivot.second * sin_r;
          double py = pivot.first * sin_r + pivot.second * cos_r;

          std::vector<std::pair<double, double>> shifted_data;
          shifted_data.reserve(rotated_data.size());
          for (auto &p : rotated_data) {
            double x = p.first;
            double y = p.second;
            if (x < px)
              x += rot_width;
            if (y < py)
              y += rot_height;
            shifted_data.push_back({x, y});
          }

          std::vector<std::pair<double, double>> s1, s2;
          for (size_t j = 0; j < shifted_data.size(); ++j) {
            if (subjects[j] == 1)
              s1.push_back(shifted_data[j]);
            else
              s2.push_back(shifted_data[j]);
          }

          trial_psi_sum += calculate_psi(s1, s2, stat_type);
          total_configs++;
        }
      }

      perm_psi[i] = total_configs > 0 ? trial_psi_sum / total_configs : 0.0;
    }
  }
};

// [[Rcpp::export]]
List distdiffR_engine(NumericMatrix data, IntegerVector subjects, int num_rot,
                      int num_shifts, int stat_type, int num_perms,
                      int seedNum) {

  // Convert NumericMatrix to a vector of pairs for faster C++ processing
  int n_pooled = data.nrow();
  std::vector<std::pair<double, double>> pooled_data;
  pooled_data.reserve(n_pooled);
  for (int i = 0; i < n_pooled; ++i) {
    pooled_data.push_back({data(i, 0), data(i, 1)});
  }

  std::vector<int> subj_vec = as<std::vector<int>>(subjects);

  // Randomly select points from the pooled data to serve as toroidal shift pivots
  std::vector<std::pair<double, double>> pivots;
  std::vector<int> indices(n_pooled);
  std::iota(indices.begin(), indices.end(), 0);
  std::mt19937 g(seedNum);
  std::shuffle(indices.begin(), indices.end(), g);
  for (int i = 0; i < std::min(num_shifts, n_pooled); ++i) {
    pivots.push_back(pooled_data[indices[i]]);
  }

  // Calculate the Psi statistic for the original (non-permuted) data
  double true_psi_sum = 0;
  int total_configs = 0;
  for (int r = 0; r < num_rot; ++r) {
    double angle = 2.0 * M_PI * r / num_rot;
    double cos_r = std::cos(angle);
    double sin_r = std::sin(angle);

    std::vector<std::pair<double, double>> rotated_data;
    rotated_data.reserve(n_pooled);
    for (auto &p : pooled_data) {
      rotated_data.push_back({p.first * cos_r - p.second * sin_r,
                              p.first * sin_r + p.second * cos_r});
    }

    double min_x = rotated_data[0].first, max_x = rotated_data[0].first;
    double min_y = rotated_data[0].second, max_y = rotated_data[0].second;
    for (auto &p : rotated_data) {
      min_x = std::min(min_x, p.first);
      max_x = std::max(max_x, p.first);
      min_y = std::min(min_y, p.second);
      max_y = std::max(max_y, p.second);
    }
    double rot_width = max_x - min_x;
    double rot_height = max_y - min_y;

    for (int s = 0; s < (int)pivots.size(); ++s) {
      auto pivot = pivots[s];
      double px = pivot.first * cos_r - pivot.second * sin_r;
      double py = pivot.first * sin_r + pivot.second * cos_r;

      std::vector<std::pair<double, double>> shifted_data;
      shifted_data.reserve(n_pooled);
      for (auto &p : rotated_data) {
        double x = p.first;
        double y = p.second;
        if (x < px)
          x += rot_width;
        if (y < py)
          y += rot_height;
        shifted_data.push_back({x, y});
      }

      std::vector<std::pair<double, double>> s1, s2;
      for (int j = 0; j < n_pooled; ++j) {
        if (subj_vec[j] == 1)
          s1.push_back(shifted_data[j]);
        else
          s2.push_back(shifted_data[j]);
      }
      true_psi_sum += calculate_psi(s1, s2, stat_type);
      total_configs++;
    }
  }
  double true_psi = total_configs > 0 ? true_psi_sum / total_configs : 0.0;

  NumericVector perm_psi_res(num_perms);
  PermutationWorker worker(pooled_data, subj_vec, 0, 0, num_rot, num_shifts,
                           stat_type, seedNum, pivots, perm_psi_res);
  RcppParallel::parallelFor(static_cast<size_t>(0),
                            static_cast<size_t>(num_perms), worker);

  return List::create(Named("psiStat") = true_psi,
                      Named("permPsi") = perm_psi_res);
}

// [[Rcpp::export]]
List grouped_distdiffR_engine(NumericMatrix data, IntegerVector subjects,
                              IntegerVector subjNums, int num_rot,
                              int num_shifts, int stat_type, int num_perms,
                              int seedNum) {

  int n_pooled = data.nrow();
  std::vector<std::pair<double, double>> pooled_data;
  pooled_data.reserve(n_pooled);
  for (int i = 0; i < n_pooled; ++i) {
    pooled_data.push_back({data(i, 0), data(i, 1)});
  }

  std::vector<int> subj_vec = as<std::vector<int>>(subjects);
  std::vector<int> snum_vec = as<std::vector<int>>(subjNums);

  std::vector<std::pair<double, double>> pivots;
  std::vector<int> indices(n_pooled);
  std::iota(indices.begin(), indices.end(), 0);
  std::mt19937 g(seedNum);
  std::shuffle(indices.begin(), indices.end(), g);
  for (int i = 0; i < std::min(num_shifts, n_pooled); ++i) {
    pivots.push_back(pooled_data[indices[i]]);
  }

  double true_psi_sum = 0;
  int total_configs = 0;
  for (int r = 0; r < num_rot; ++r) {
    double angle = 2.0 * M_PI * r / num_rot;
    double cos_r = std::cos(angle);
    double sin_r = std::sin(angle);

    std::vector<std::pair<double, double>> rotated_data;
    rotated_data.reserve(n_pooled);
    for (auto &p : pooled_data) {
      rotated_data.push_back({p.first * cos_r - p.second * sin_r,
                              p.first * sin_r + p.second * cos_r});
    }

    double min_x = rotated_data[0].first, max_x = rotated_data[0].first;
    double min_y = rotated_data[0].second, max_y = rotated_data[0].second;
    for (auto &p : rotated_data) {
      min_x = std::min(min_x, p.first);
      max_x = std::max(max_x, p.first);
      min_y = std::min(min_y, p.second);
      max_y = std::max(max_y, p.second);
    }
    double rot_width = max_x - min_x;
    double rot_height = max_y - min_y;

    for (int s = 0; s < (int)pivots.size(); ++s) {
      auto pivot = pivots[s];
      double px = pivot.first * cos_r - pivot.second * sin_r;
      double py = pivot.first * sin_r + pivot.second * cos_r;

      std::vector<std::pair<double, double>> shifted_data;
      shifted_data.reserve(n_pooled);
      for (auto &p : rotated_data) {
        double x = p.first;
        double y = p.second;
        if (x < px)
          x += rot_width;
        if (y < py)
          y += rot_height;
        shifted_data.push_back({x, y});
      }

      std::vector<std::pair<double, double>> s1, s2;
      std::vector<int> subs1, subs2;
      for (int j = 0; j < n_pooled; ++j) {
        if (subj_vec[j] == 1) {
          s1.push_back(shifted_data[j]);
          subs1.push_back(snum_vec[j]);
        } else {
          s2.push_back(shifted_data[j]);
          subs2.push_back(snum_vec[j]);
        }
      }
      true_psi_sum += calculate_psi_grouped(s1, s2, subs1, subs2);
      total_configs++;
    }
  }
  double true_psi = total_configs > 0 ? true_psi_sum / total_configs : 0.0;

  NumericVector perm_psi_res(num_perms);
  std::mt19937 g_perm(seedNum);
  std::vector<int> shuffled_subjects = subj_vec;
  for (int i = 0; i < num_perms; ++i) {
    std::shuffle(shuffled_subjects.begin(), shuffled_subjects.end(), g_perm);
    double trial_psi_sum = 0;
    int t_configs = 0;
    for (int r = 0; r < num_rot; ++r) {
      double angle = 2.0 * M_PI * r / num_rot;
      double cos_r = std::cos(angle);
      double sin_r = std::sin(angle);
      std::vector<std::pair<double, double>> rotated_data;
      rotated_data.reserve(n_pooled);
      for (auto &p : pooled_data)
        rotated_data.push_back({p.first * cos_r - p.second * sin_r,
                                p.first * sin_r + p.second * cos_r});
      double min_x = rotated_data[0].first, max_x = rotated_data[0].first;
      double min_y = rotated_data[0].second, max_y = rotated_data[0].second;
      for (auto &p : rotated_data) {
        min_x = std::min(min_x, p.first);
        max_x = std::max(max_x, p.first);
        min_y = std::min(min_y, p.second);
        max_y = std::max(max_y, p.second);
      }
      double rot_width = max_x - min_x;
      double rot_height = max_y - min_y;
      for (int s = 0; s < (int)pivots.size(); ++s) {
        auto pivot = pivots[s];
        double px = pivot.first * cos_r - pivot.second * sin_r;
        double py = pivot.first * sin_r + pivot.second * cos_r;
        std::vector<std::pair<double, double>> shifted_data;
        shifted_data.reserve(n_pooled);
        for (auto &p : rotated_data) {
          double x = p.first;
          double y = p.second;
          if (x < px)
            x += rot_width;
          if (y < py)
            y += rot_height;
          shifted_data.push_back({x, y});
        }
        std::vector<std::pair<double, double>> s1, s2;
        std::vector<int> subs1, subs2;
        for (int j = 0; j < n_pooled; ++j) {
          if (shuffled_subjects[j] == 1) {
            s1.push_back(shifted_data[j]);
            subs1.push_back(snum_vec[j]);
          } else {
            s2.push_back(shifted_data[j]);
            subs2.push_back(snum_vec[j]);
          }
        }
        trial_psi_sum += calculate_psi_grouped(s1, s2, subs1, subs2);
        t_configs++;
      }
    }
    perm_psi_res[i] = t_configs > 0 ? trial_psi_sum / t_configs : 0.0;
  }

  return List::create(Named("psiStat") = true_psi,
                      Named("permPsi") = perm_psi_res);
}
