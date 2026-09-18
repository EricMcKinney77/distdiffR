library(distdiffR)

# Benchmark settings
sample_sizes <- c(25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600, 1000, 2000)
num_perms_benchmark <- 999 # Use a small number of permutations for speed
sample_size_target <- 7000 # Target sample size for extrapolation
num_runs <- 10 # Number of runs per sample size for CI
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

# Save results to CSV
results <- data.frame(
  "Sample Size (n)" = sample_sizes,
  "Mean Time (s)" = means,
  "Std Error (s)" = ses
)
write.csv(results, "benchmark_results.csv", row.names = FALSE)
cat("\nResults saved to benchmark_results.csv\n")
