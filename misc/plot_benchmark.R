library(ggplot2)

# Load benchmark results
if (!file.exists("misc/benchmark_results.csv")) {
  stop("benchmark_results.csv not found. Please run benchmark_distdiffR.R first.")
}
results <- read.csv("misc/benchmark_results.csv", check.names = FALSE)
sample_sizes <- results$"Sample Size (n)"
means <- results$"Mean Time (s)"
ses <- results$"Std Error (s)"

# Settings for extrapolation
sample_size_target <- 7000

# Create ggplot object
p <- ggplot(results, aes(x = `Sample Size (n)`, y = `Mean Time (s)`)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  geom_errorbar(aes(ymin = `Mean Time (s)` - 1.96 * `Std Error (s)`, 
                    ymax = `Mean Time (s)` + 1.96 * `Std Error (s)`), 
                width = 0.2, color = "blue") +
  scale_x_continuous(breaks = sample_sizes, limits = c(10, 2015), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 62), breaks = seq(0, 60, by = 5), expand = expansion(mult = c(0, 0))) +
  labs(x = "Sample Size (n per sample)", 
       y = "Time (sec)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.grid.major = element_line(color = "grey", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  )



# Save plot to PDF
ggsave("misc/benchmark_plot.pdf", plot = p, width = 12, height = 6)
cat("\nPlot saved to benchmark_plot.pdf\n")

# Estimation for Target Sample Size
n_last <- tail(sample_sizes, 1)
m_last <- tail(means, 1)
se_last <- tail(ses, 1)
n_target <- sample_size_target

# Complexity factor for extrapolation
# Time(target) = Time(last) * (n_target/n_last) * (log(n_target)/log(n_last))
scale_factor <- (n_target / n_last) * (log(n_target) / log(n_last))

est_total_mean_secs <- m_last * scale_factor
est_total_se_secs <- se_last * scale_factor

# 95% Confidence Interval
lower_ci <- est_total_mean_secs - 1.96 * est_total_se_secs
upper_ci <- est_total_mean_secs + 1.96 * est_total_se_secs

cat(sprintf("\n--- Estimation for N = %d (95%% Confidence Interval) ---\n", n_target))
cat(sprintf(
  "Estimated Mean Time: %.2f seconds (%.2f minutes)\n",
  est_total_mean_secs, est_total_mean_secs / 60
))
cat(sprintf("95%% CI: [%.2f, %.2f] seconds\n", lower_ci, upper_ci))
