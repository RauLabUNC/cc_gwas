#!/usr/bin/env Rscript
# generate_report.R - Generate HTML resource usage report from SLURM data

args <- commandArgs(trailingOnly = TRUE)
in_tsv <- if (length(args) >= 1) args[[1]] else "results/resource_reports/slurm_usage.tsv"
out_html <- if (length(args) >= 2) args[[2]] else "results/resource_reports/usage_report.html"

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  
  # Try to load gt, but don't fail if it's missing
  gt_available <- requireNamespace("gt", quietly = TRUE)
  if (gt_available) {
    library(gt)
  } else {
    warning("Package 'gt' not available. Will create simpler HTML tables.")
  }
})

# Helper functions
to_seconds <- function(x) {
  x <- ifelse(is.na(x) | x == "" | x == "Unknown", NA_character_, x)
  vapply(x, function(s) {
    if (is.na(s)) return(NA_real_)
    if (grepl("-", s)) {
      parts <- strsplit(s, "-", fixed = TRUE)[[1]]
      days <- as.numeric(parts[1])
      hms <- parts[2]
    } else {
      days <- 0
      hms <- s
    }
    hhmmss <- as.numeric(strsplit(hms, ":", fixed = TRUE)[[1]])
    if (length(hhmmss) != 3) return(NA_real_)
    days * 86400 + hhmmss[1] * 3600 + hhmmss[2] * 60 + hhmmss[3]
  }, numeric(1))
}

parse_mem_mb <- function(x) {
  ifelse(is.na(x) | x == "" | x == "Unknown", NA_real_,
    vapply(x, function(s) {
      s <- trimws(s)
      num <- as.numeric(sub("^([0-9.]+).*", "\\1", s))
      unit <- sub("^[0-9.]+\\s*([KkMmGgTt]).*$", "\\1", s)
      if (is.na(unit) || nchar(unit) == 0) unit <- "M"
      mult <- switch(toupper(unit),
        "K" = 1/1024,
        "M" = 1,
        "G" = 1024,
        "T" = 1024 * 1024,
        1
      )
      num * mult
    }, numeric(1))
  )
}

fmt_secs <- function(x) {
  h <- floor(x / 3600)
  m <- floor((x - h * 3600) / 60)
  s <- round(x - h * 3600 - m * 60)
  sprintf("%02d:%02d:%02d", h, m, s)
}

# Read and process data
dat <- read_delim(in_tsv, delim = "|", col_types = cols(.default = "c"))

# If no SLURM data, try to use benchmark files
if (nrow(dat) == 0) {
  cat("No SLURM data found. Checking for benchmark files...\n")
  
  # Look for benchmark files
  bench_dir <- dirname(in_tsv)
  if (dir.exists(file.path(bench_dir, "bench"))) {
    bench_files <- list.files(file.path(bench_dir, "bench"), 
                              pattern = "\\.txt$", 
                              recursive = TRUE, 
                              full.names = TRUE)
    
    if (length(bench_files) > 0) {
      cat("Found", length(bench_files), "benchmark files. Creating report from benchmarks.\n")
      
      # Create simple HTML report from benchmark data
      bench_data <- lapply(bench_files, function(f) {
        df <- read.table(f, header = TRUE, sep = "\t")
        df$file <- basename(dirname(f))
        df$rule <- basename(dirname(dirname(f)))
        df
      })
      bench_df <- do.call(rbind, bench_data)
      
      # Create basic HTML report
      html <- paste0(
        "<!doctype html><html><head><meta charset='utf-8'>",
        "<title>Resource Report (Benchmark Data)</title>",
        "<style>body{font-family:sans-serif;margin:24px;}table{border-collapse:collapse;}",
        "th,td{border:1px solid #ddd;padding:8px;text-align:left;}th{background:#f2f2f2;}</style>",
        "</head><body>",
        "<h1>Resource Usage Report (from Snakemake Benchmarks)</h1>",
        "<p>Note: SLURM accounting data was not available. Showing Snakemake benchmark data instead.</p>",
        "<table><tr><th>Rule</th><th>Job</th><th>Runtime (s)</th><th>Max RSS (MB)</th><th>CPU Time (s)</th></tr>"
      )
      
      for (i in 1:nrow(bench_df)) {
        html <- paste0(html, 
          "<tr><td>", bench_df$rule[i], "</td>",
          "<td>", bench_df$file[i], "</td>",
          "<td>", round(bench_df$s[i], 1), "</td>",
          "<td>", round(bench_df$max_rss[i], 0), "</td>",
          "<td>", round(bench_df$cpu_time[i], 1), "</td></tr>"
        )
      }
      
      html <- paste0(html, "</table></body></html>")
      writeLines(html, out_html)
      message(sprintf("Created fallback report from benchmarks: %s", out_html))
      quit(status = 0)
    }
  }
  
  cat("No data available for report generation.\n")
  quit(status = 0)
}

dat2 <- dat %>%
  mutate(
    AllocCPUS = suppressWarnings(as.numeric(AllocCPUS)),
    CPUTimeRAW = suppressWarnings(as.numeric(CPUTimeRAW)),
    elapsed_sec = to_seconds(Elapsed),
    user_sec = to_seconds(UserCPU),
    sys_sec = to_seconds(SystemCPU),
    timelimit_sec = to_seconds(Timelimit),
    req_mem_mb = parse_mem_mb(ReqMem),
    maxrss_mb = parse_mem_mb(MaxRSS),
    averss_mb = parse_mem_mb(AveRSS),
    avevm_mb = parse_mem_mb(AveVMSize),
    cpu_eff = ifelse(!is.na(CPUTimeRAW) & !is.na(AllocCPUS) & !is.na(elapsed_sec) & 
                     AllocCPUS > 0 & elapsed_sec > 0,
                     CPUTimeRAW / (AllocCPUS * elapsed_sec), NA_real_),
    cpu_util_pct = 100 * cpu_eff,
    mem_eff = ifelse(!is.na(maxrss_mb) & !is.na(req_mem_mb) & req_mem_mb > 0,
                     maxrss_mb / req_mem_mb, NA_real_),
    walltime_used = ifelse(!is.na(elapsed_sec) & !is.na(timelimit_sec) & timelimit_sec > 0,
                           elapsed_sec / timelimit_sec, NA_real_),
    rule = sub("^smk-([^ -]+).*", "\\1", JobName)
  )

# Per-job summary
per_job <- dat2 %>%
  filter(!is.na(elapsed_sec), !is.na(AllocCPUS)) %>%
  transmute(
    JobID, JobName, State, AllocCPUS,
    Elapsed = elapsed_sec,
    CPU_utilization = cpu_util_pct,
    CPUTimeRAW, UserCPU = user_sec, SystemCPU = sys_sec,
    ReqMem_MB = req_mem_mb, MaxRSS_MB = maxrss_mb,
    Mem_efficiency = mem_eff, Walltime_used = walltime_used, ExitCode
  )

# Per-rule summary
per_rule <- per_job %>%
  mutate(rule = sub("^smk-([^ -]+).*", "\\1", JobName)) %>%
  group_by(rule) %>%
  summarise(
    jobs = n(),
    alloc_cpus_mean = mean(AllocCPUS, na.rm = TRUE),
    elapsed_p50 = median(Elapsed, na.rm = TRUE),
    elapsed_p95 = quantile(Elapsed, 0.95, na.rm = TRUE),
    cpu_util_p50 = median(CPU_utilization, na.rm = TRUE),
    cpu_util_p95 = quantile(CPU_utilization, 0.95, na.rm = TRUE),
    mem_eff_p50 = median(Mem_efficiency, na.rm = TRUE),
    total_elapsed = sum(Elapsed, na.rm = TRUE),
    total_cpu_time = sum(CPUTimeRAW, na.rm = TRUE)
  ) %>%
  arrange(desc(total_elapsed))

# Create plots
dir.create("results/resource_reports/plots", showWarnings = FALSE, recursive = TRUE)

# Plot 1: Total wallclock by rule
if (nrow(per_rule) > 0) {
  p1 <- ggplot(per_rule, aes(x = reorder(rule, total_elapsed), y = total_elapsed / 3600)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "Total Cluster Wallclock by Rule",
         x = "Rule", y = "Hours") +
    theme_minimal()
  ggsave("results/resource_reports/plots/wallclock_by_rule.png", p1,
         width = 8, height = max(4, nrow(per_rule) * 0.3), dpi = 120)
}

# Plot 2: CPU utilization histogram
if (sum(!is.na(per_job$CPU_utilization)) > 0) {
  p2 <- per_job %>%
    filter(!is.na(CPU_utilization)) %>%
    ggplot(aes(x = CPU_utilization)) +
    geom_histogram(bins = 40, fill = "coral", color = "white") +
    labs(title = "Per-job CPU Utilization Distribution",
         x = "CPU Utilization (%)", y = "Count") +
    theme_minimal() +
    geom_vline(xintercept = 70, linetype = "dashed", color = "green", alpha = 0.5) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "orange", alpha = 0.5)
  ggsave("results/resource_reports/plots/cpu_util_hist.png", p2,
         width = 7, height = 4, dpi = 120)
}

# Plot 3: Memory efficiency scatter
if (sum(!is.na(per_job$Mem_efficiency)) > 0) {
  p3 <- per_job %>%
    filter(!is.na(Mem_efficiency)) %>%
    mutate(rule = sub("^smk-([^ -]+).*", "\\1", JobName)) %>%
    ggplot(aes(x = ReqMem_MB, y = MaxRSS_MB, color = rule)) +
    geom_point(alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_x_log10() + scale_y_log10() +
    labs(title = "Memory Usage: Requested vs Actual",
         x = "Requested Memory (MB)", y = "Max RSS (MB)") +
    theme_minimal()
  ggsave("results/resource_reports/plots/mem_efficiency.png", p3,
         width = 8, height = 6, dpi = 120)
}

# Format tables for HTML
per_job_tbl <- per_job %>%
  mutate(
    Elapsed = fmt_secs(Elapsed),
    CPU_utilization = ifelse(is.na(CPU_utilization), NA,
                            sprintf("%.1f%%", CPU_utilization)),
    Mem_efficiency = ifelse(is.na(Mem_efficiency), NA,
                           sprintf("%.0f%%", 100 * Mem_efficiency)),
    Walltime_used = ifelse(is.na(Walltime_used), NA,
                          sprintf("%.0f%%", 100 * Walltime_used)),
    ReqMem_MB = round(ReqMem_MB, 0),
    MaxRSS_MB = round(MaxRSS_MB, 0)
  ) %>%
  gt() %>%
  tab_header(title = md("**Per-job Resource Usage (SLURM)**")) %>%
  cols_label(
    JobID = "Job ID",
    JobName = "Job Name",
    State = "State",
    AllocCPUS = "CPUs",
    Elapsed = "Elapsed",
    CPU_utilization = "CPU Util",
    CPUTimeRAW = "CPU Time (s)",
    UserCPU = "User (s)",
    SystemCPU = "System (s)",
    ReqMem_MB = "Req Mem (MB)",
    MaxRSS_MB = "Max RSS (MB)",
    Mem_efficiency = "Mem Eff",
    Walltime_used = "Time Used",
    ExitCode = "Exit"
  ) %>%
  tab_options(table.font.size = px(12), data_row.padding = px(4))

per_rule_tbl <- per_rule %>%
  mutate(
    elapsed_p50 = fmt_secs(elapsed_p50),
    elapsed_p95 = fmt_secs(elapsed_p95),
    cpu_util_p50 = sprintf("%.1f%%", cpu_util_p50),
    cpu_util_p95 = sprintf("%.1f%%", cpu_util_p95),
    mem_eff_p50 = sprintf("%.0f%%", 100 * mem_eff_p50),
    total_elapsed = fmt_secs(total_elapsed),
    total_cpu_time = fmt_secs(total_cpu_time),
    alloc_cpus_mean = round(alloc_cpus_mean, 1)
  ) %>%
  gt() %>%
  tab_header(title = md("**Per-rule Summary (Bottlenecks & Efficiency)**")) %>%
  cols_label(
    rule = "Rule",
    jobs = "Jobs",
    alloc_cpus_mean = "Avg CPUs",
    elapsed_p50 = "Median Time",
    elapsed_p95 = "P95 Time",
    cpu_util_p50 = "Median CPU",
    cpu_util_p95 = "P95 CPU",
    mem_eff_p50 = "Median Mem",
    total_elapsed = "Total Wall",
    total_cpu_time = "Total CPU"
  ) %>%
  tab_options(table.font.size = px(12), data_row.padding = px(4))

# Generate HTML report
html <- paste0(
  "<!doctype html><html><head><meta charset='utf-8'>",
  "<title>Snakemake SLURM Resource Report</title>",
  "<style>",
  "body{font-family:system-ui,-apple-system,sans-serif;margin:24px;background:#f8f9fa;}",
  "h1{color:#2c3e50;border-bottom:3px solid #3498db;padding-bottom:10px;}",
  "h2{color:#34495e;margin-top:30px;}",
  ".note{color:#7f8c8d;margin:10px 0;font-size:14px;}",
  ".container{background:white;padding:20px;border-radius:8px;box-shadow:0 2px 4px rgba(0,0,0,0.1);margin:20px 0;}",
  "img{max-width:100%;height:auto;border:1px solid #ddd;margin:12px 0;padding:8px;border-radius:6px;background:white;}",
  ".legend{background:#ecf0f1;padding:10px;border-radius:4px;margin:10px 0;}",
  ".legend h3{margin:0 0 10px 0;color:#2c3e50;}",
  ".legend ul{margin:5px 0;padding-left:20px;}",
  "</style>",
  "</head><body>",
  "<h1>🔬 Snakemake SLURM Resource Usage Report</h1>",
  sprintf("<div class='note'>📊 Source: %s • Generated: %s</div>",
          in_tsv, format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  
  "<div class='container'>",
  "<h2>📈 Per-rule Totals & Bottlenecks</h2>",
  as_raw_html(per_rule_tbl),
  "<img src='plots/wallclock_by_rule.png' alt='Total wallclock by rule'>",
  "</div>",
  
  "<div class='container'>",
  "<h2>💻 CPU Utilization Analysis</h2>",
  "<div class='legend'>",
  "<h3>Interpretation Guide:</h3>",
  "<ul>",
  "<li><strong>&lt; 50%</strong> (orange line): Over-threaded or I/O bound → reduce threads</li>",
  "<li><strong>70-90%</strong> (green line): Healthy utilization</li>",
  "<li><strong>&gt; 95%</strong>: CPU-bound → consider increasing threads if scalable</li>",
  "</ul>",
  "</div>",
  "<img src='plots/cpu_util_hist.png' alt='CPU utilization histogram'>",
  "</div>",
  
  "<div class='container'>",
  "<h2>🧠 Memory Efficiency</h2>",
  "<div class='legend'>",
  "<h3>How to read:</h3>",
  "<ul>",
  "<li>Points <strong>below the line</strong>: Over-provisioned memory → reduce allocation</li>",
  "<li>Points <strong>on the line</strong>: Perfect estimation</li>",
  "<li>Points <strong>above the line</strong>: Under-provisioned (should not happen)</li>",
  "</ul>",
  "</div>",
  "<img src='plots/mem_efficiency.png' alt='Memory efficiency scatter'>",
  "</div>",
  
  "<div class='container'>",
  "<h2>📋 Per-job Details</h2>",
  as_raw_html(per_job_tbl),
  "</div>",
  
  "<div class='container'>",
  "<h2>📖 Metrics Glossary</h2>",
  "<ul>",
  "<li><strong>CPU Utilization</strong> = CPUTimeRAW / (AllocCPUS × Elapsed)</li>",
  "<li><strong>Memory Efficiency</strong> = MaxRSS / ReqMem</li>",
  "<li><strong>Walltime Used</strong> = Elapsed / Timelimit</li>",
  "<li><strong>P50/P95</strong> = 50th/95th percentile values</li>",
  "</ul>",
  "</div>",
  
  "</body></html>"
)

writeLines(html, out_html)
message(sprintf("✅ Report generated: %s", out_html))