library(distdiffR)

# Benchmark settings
sample_sizes <- seq(100, 2000, by = 200) # Test from 100 to 2000
num_perms_benchmark <- 10 # Use a small number of permutations for speed
num_perms_target <- 999 # Actual target number of permutations
num_runs <- 5 # Number of runs per sample size for CI
seedNum <- 123

# Matrix to store times: rows = runs, cols = sample sizes
times_matrix <- matrix(NA, nrow = num_runs, ncol = length(sample_sizes))

cat("Starting benchmark with", num_runs, "runs per sample size...\n")
cat(sprintf("%-15s | %-15s | %-15s\n", "Sample Size (n)", "Mean Time (s)", "Std Error (s)"))
cat(paste0(rep("-", 48), collapse = ""), "\n")

for (j in seq_along(sample_sizes)) {
  n <- sample_sizes[j]

  for (i in 1:num_runs) {
    # Generate fresh random bivariate data for each run to avoid caching effects
    data1 <- matrix(rnorm(n * 2), ncol = 2)
    data2 <- matrix(rnorm(n * 2), ncol = 2)

    start_time <- Sys.time()
    distdiffr(data1,
      data2,
      numPerms = num_perms_benchmark,
      seedNum = seedNum
    )
    end_time <- Sys.time()

    times_matrix[i, j] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  }

  m <- mean(times_matrix[, j])
  se <- sd(times_matrix[, j]) / sqrt(num_runs)
  cat(sprintf("%-15d | %-15.4f | %-15.4f\n", n, m, se))
}

means <- colMeans(times_matrix)
sds <- apply(times_matrix, 2, sd)
ses <- sds / sqrt(num_runs)

# Plotting with error bars
plot(sample_sizes, means,
  type = "b", pch = 19, col = "blue",
  xlab = "Sample Size (n per group)",
  ylab = "Time for 10 permutations (sec)",
  main = "Computational Scaling of distdiffR with 95% CI"
)
# Add 95% CI error bars (approx 1.96 * SE)
arrows(sample_sizes, means - 1.96 * ses, sample_sizes, means + 1.96 * ses,
  code = 3, angle = 90, length = 0.05, col = "blue"
)
grid()

# Estimation for N = 7000
n_last <- tail(sample_sizes, 1)
m_last <- tail(means, 1)
se_last <- tail(ses, 1)
n_target <- 7000

# Complexity factor for extrapolation
# Time(target) = Time(last) * (n_target/n_last) * (log(n_target)/log(n_last))
scale_factor <- (n_target / n_last) * (log(n_target) / log(n_last))
perm_factor <- (num_perms_target / num_perms_benchmark)

est_mean_10_perms <- m_last * scale_factor
est_se_10_perms <- se_last * scale_factor

est_total_mean_secs <- est_mean_10_perms * perm_factor
est_total_se_secs <- est_se_10_perms * perm_factor

# 95% Confidence Interval
lower_ci <- est_total_mean_secs - 1.96 * est_total_se_secs
upper_ci <- est_total_mean_secs + 1.96 * est_total_se_secs

cat("\n--- Estimation for N = 7000 (95% Confidence Interval) ---\n")
cat(sprintf(
  "Estimated Mean Time: %.2f seconds (%.2f minutes)\n",
  est_total_mean_secs, est_total_mean_secs / 60
))
cat(sprintf("95%% CI: [%.2f, %.2f] seconds\n", lower_ci, upper_ci))
