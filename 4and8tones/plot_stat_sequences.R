# %% Load Libraries and Setup
# HDI + ROPE Analysis for Markov Sequences in R

# Load required libraries
library(tidyverse)
library(ggplot2)
library(HDInterval)
library(gridExtra)

# Clear workspace
rm(list = ls())

# %% Set Directories and Load Data
in_dir <- "/home/gdemarchi/Work/Data/markov_tng/sequences/"
seqs4_file <- "4_tone_sequences.txt"
seqs8_file <- "8_tone_sequences.txt"
out_dir <- "/home/gdemarchi/Work/Data/markov_tng/figures/"

# Load data
seqs4_data <- read.table(paste0(in_dir, seqs4_file), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
seqs8_data <- read.table(paste0(in_dir, seqs8_file), sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# %% Check Subjects
subjList <- c(
  '19920916slsa', '19940511krln', '20010127urdr', '20010618ttbl',
  '20010629tnli', '20011116mrsr', '20020521bidv', '20020827ktsh',
  '20021005asap', '20021230urre', '20030208crsh', '20030225aehf',
  '20030328ptbn', '20030502rnhn', '20030516ktba', '20030717begl',
  '20030812argn', '20030829mroc', '20030901adjg', '20030919mrgs',
  '20031221mlht', '20040210gidn', '20040308pgwl', '20040330brhr',
  '20040416hlln', '20031013aglv', '20011126sbgi'
)

# Check for missing subjects
missing_in_4 <- setdiff(subjList, seqs4_data$V1)
missing_in_8 <- setdiff(subjList, seqs8_data$V1)
extra_in_4 <- setdiff(seqs4_data$V1, subjList)
extra_in_8 <- setdiff(seqs8_data$V1, subjList)

if (length(missing_in_4) == 0 && length(missing_in_8) == 0 && 
    length(extra_in_4) == 0 && length(extra_in_8) == 0) {
  cat('No missing data!\n')
} else {
  if (length(missing_in_4) > 0) { cat('Missing in 4-tone:\n'); print(missing_in_4) }
  if (length(missing_in_8) > 0) { cat('Missing in 8-tone:\n'); print(missing_in_8) }
  if (length(extra_in_4) > 0) { cat('Extra in 4-tone:\n'); print(extra_in_4) }
  if (length(extra_in_8) > 0) { cat('Extra in 8-tone:\n'); print(extra_in_8) }
}

# %% Setup Tone Mappings and Extract Sequences
trigs4 <- c(129, 130, 131, 132)
freqs4 <- c(200, 431, 928, 2000)
trigs8 <- c(193, 194, 195, 196, 197, 198, 199, 200)
freqs8 <- round(10^(seq(log10(200), log10(2000), length.out = 8)))

n_subj4 <- nrow(seqs4_data)
n_subj8 <- nrow(seqs8_data)
x4 <- 1:4
x8 <- 1:8

# Extract sequences for 4-tone step by step
seq_list4 <- matrix(0, nrow = n_subj4, ncol = 4)
for (i in 1:n_subj4) {
  seq_string <- seqs4_data$V2[i]
  seq_parts <- strsplit(seq_string, ",")[[1]]
  seq_numeric <- as.numeric(seq_parts)
  seq_list4[i, ] <- seq_numeric
}

# Extract sequences for 8-tone step by step
seq_list8 <- matrix(0, nrow = n_subj8, ncol = 8)
for (i in 1:n_subj8) {
  seq_string <- seqs8_data$V2[i]
  seq_parts <- strsplit(seq_string, ",")[[1]]
  seq_numeric <- as.numeric(seq_parts)
  seq_list8[i, ] <- seq_numeric
}

# %% Perform Regressions
# 4-tone data
y_freq4_all <- matrix(0, nrow = n_subj4, ncol = 4)
y_log4_all <- matrix(0, nrow = n_subj4, ncol = 4)
slopes4_log <- numeric(n_subj4)

for (i in 1:n_subj4) {
  seq <- seq_list4[i, ]
  for (j in 1:4) {
    y_freq4_all[i, j] <- freqs4[which(trigs4 == seq[j])]
  }
  y_log4_all[i, ] <- log10(y_freq4_all[i, ])
  
  # Linear regression
  fit <- lm(y_log4_all[i, ] ~ x4)
  slopes4_log[i] <- coef(fit)[2]  # slope 
}

# 8-tone data
y_freq8_all <- matrix(0, nrow = n_subj8, ncol = 8)
y_log8_all <- matrix(0, nrow = n_subj8, ncol = 8)
slopes8_log <- numeric(n_subj8)

for (i in 1:n_subj8) {
  seq <- seq_list8[i, ]
  for (j in 1:8) {
    y_freq8_all[i, j] <- freqs8[which(trigs8 == seq[j])]
  }
  y_log8_all[i, ] <- log10(y_freq8_all[i, ])
  
  # Linear regression 
  fit <- lm(y_log8_all[i, ] ~ x8)
  slopes8_log[i] <- coef(fit)[2]  # slope
}

# Classical ordered sequences (like in the original paper)
y_ord_log4 <- log10(freqs4)
y_ord_log8 <- log10(freqs8)
fit_ord4 <- lm(y_ord_log4 ~ x4)
fit_ord8 <- lm(y_ord_log8 ~ x8)
p_ord4_slope <- coef(fit_ord4)[2]
p_ord8_slope <- coef(fit_ord8)[2]

# Print results
cat(sprintf("Classical ordered slope 4-tone: %.3f\n", p_ord4_slope))
cat(sprintf("Classical ordered slope 8-tone: %.3f\n", p_ord8_slope))


# %% Contact Sheet Plots for 4-tone and 8-tone Regressions

# Create contact sheet plots for regression visualization - 4 tones
n_plots_4 <- n_subj4  # Show all subjects
plot_list_4 <- list()

# Individual subject plots for 4-tone
for (i in 1:n_plots_4) {
  df <- data.frame(x = x4, y = y_log4_all[i, ])
  
  plot_list_4[[i]] <- ggplot(df, aes(x = x, y = y)) +
    geom_point(size = 2, color = "blue") +
    geom_line(linetype = "dashed", color = "blue", size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "black", size = 1) +
    labs(title = seqs4_data$V1[i], x = "Tone Position", y = "Log10(Frequency)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, hjust = 0.5)) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = sprintf("r = %.2f", slopes4_log[i]), size = 3, bg = "white")
}

# Classical ordered sequence plot for 4-tone (last plot)
df_classical_4 <- data.frame(x = x4, y = y_ord_log4)
plot_list_4[[n_plots_4 + 1]] <- ggplot(df_classical_4, aes(x = x, y = y)) +
  geom_point(size = 2, color = "red") +
  geom_line(linetype = "dashed", color = "red", size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 1) +
  labs(title = "Classical", x = "Tone Position", y = "Log10(Frequency)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, hjust = 0.5)) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = sprintf("r = %.2f", p_ord4_slope), size = 3, bg = "white")

# Arrange 4-tone plots in a grid
contact_sheet_4 <- grid.arrange(grobs = plot_list_4, ncol = 4,
                                top = "4-tone Regression Contact Sheet")


# 8-tone contact sheet
n_plots_8 <- n_subj8  # Show all subjects
plot_list_8 <- list()

# Individual subject plots for 8-tone
for (i in 1:n_plots_8) {
  df <- data.frame(x = x8, y = y_log8_all[i, ])
  
  plot_list_8[[i]] <- ggplot(df, aes(x = x, y = y)) +
    geom_point(size = 2, color = "blue") +
    geom_line(linetype = "dashed", color = "blue", size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "black", size = 1) +
    labs(title = seqs8_data$V1[i], x = "Tone Position", y = "Log10(Frequency)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, hjust = 0.5)) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = sprintf("r = %.2f", slopes8_log[i]), size = 3, bg = "white")
}

# Classical ordered sequence plot for 8-tone (last plot)
df_classical_8 <- data.frame(x = x8, y = y_ord_log8)
plot_list_8[[n_plots_8 + 1]] <- ggplot(df_classical_8, aes(x = x, y = y)) +
  geom_point(size = 2, color = "red") +
  geom_line(linetype = "dashed", color = "red", size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 1) +
  labs(title = "Classical", x = "Tone Position", y = "Log10(Frequency)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, hjust = 0.5)) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = sprintf("r = %.2f", p_ord8_slope), size = 3, bg = "white")

# Arrange 8-tone plots in a grid
contact_sheet_8 <- grid.arrange(grobs = plot_list_8, ncol = 4,
                                top = "8-tone Regression Contact Sheet")

# # %% Save contact sheets
# ggsave(paste0(out_dir, "contact_sheet_4tone_regressions.png"), contact_sheet_4, 
#        width = 16, height = 12, dpi = 300)
# ggsave(paste0(out_dir, "contact_sheet_8tone_regressions.png"), contact_sheet_8, 
#        width = 16, height = 12, dpi = 300)

# %% HDI + ROPE  Analysis 
rope_width <- 0.10  # width of practical equivalence zone around zero -> s. Kruschke, 2018
rope_low <- -rope_width
rope_high <- rope_width
hdi_prob <- 0.89


# Reuse slopes already computed above
data4 <- slopes4_log[is.finite(slopes4_log)]
data8 <- slopes8_log[is.finite(slopes8_log)]

# Compute 89% HDI using HDInterval package
# Following Kruschke (2018) Bayesian approach for posterior credible intervals
hdi4 <- hdi(data4, credMass = hdi_prob)
hdi8 <- hdi(data8, credMass = hdi_prob)

# Probability calculations (within ROPE)
p4_in_rope <- mean((data4 >= rope_low) & (data4 <= rope_high))
p8_in_rope <- mean((data8 >= rope_low) & (data8 <= rope_high))


cat(sprintf("4-tone: %.1f%% subjects have slopes practically equivalent to ZERO\n", p4_in_rope*100))
cat(sprintf("8-tone: %.1f%% subjects have slopes practically equivalent to ZERO\n", p8_in_rope*100))

# %% Create Visualization - ROPE centered on ZERO
# Create HDI + ROPE plots showing subjects are centered around zero

# Calculate the maximum density for both datasets to set consistent y-axis
max_density_4 <- max(density(data4)$y)
max_density_8 <- max(density(data8)$y)
y_max <- max(max_density_4, max_density_8) * 1.1  # Add 10% padding

# Plot for 4-tone
p1 <- ggplot(data.frame(slopes = data4), aes(x = slopes)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", alpha = 0.7, color = "black") +
  geom_density(color = "blue", size = 1.2) +
  
  # ROPE 
  annotate("rect", xmin = rope_low, xmax = rope_high, ymin = 0, ymax = Inf, 
           fill = "green", alpha = 0.3) +
  
  # HDI
  annotate("rect", xmin = hdi4[1], xmax = hdi4[2], ymin = 0, ymax = Inf, 
           fill = "blue", alpha = 0.1) +
  
  # Zero line 
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", size = 0.5) +
  
  # Classical ordered slope line
  geom_vline(xintercept = p_ord4_slope, color = "red", linetype = "dashed", size = 1.5) +
  
  labs(title = "4-tone Slopes",
       x = "Slope", y = "Density") +
  xlim(-0.5, 0.5) +
  ylim(0, y_max) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold")) +
  
  # Add ROPE probability text in top left
  annotate("text", x = -0.45, y = y_max * 0.95, vjust = 1, hjust = 0,
           label = sprintf("Pr(ROPE) = %.2f", p4_in_rope),
           size = 4, bg = "white") +
  
  # Add classical slope value near red line
  annotate("text", x = p_ord4_slope, y = y_max * 0.8, vjust = 1, hjust = -0.1,
           label = bquote(r[c] == .(sprintf("%.2f", p_ord4_slope))),
           size = 4, color = "red", bg = "white")

# Plot for 8-tone  
p2 <- ggplot(data.frame(slopes = data8), aes(x = slopes)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", alpha = 0.7, color = "black") +
  geom_density(color = "blue", size = 1.2) +
  
  # ROPE (centered on ZERO)
  annotate("rect", xmin = rope_low, xmax = rope_high, ymin = 0, ymax = Inf, 
           fill = "green", alpha = 0.3) +
  
  # HDI
  annotate("rect", xmin = hdi8[1], xmax = hdi8[2], ymin = 0, ymax = Inf, 
           fill = "blue", alpha = 0.1) +
  
  # Zero line (where subjects actually cluster)
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", size = 0.5) +
  
  # Theoretical ordered slope line (far from where subjects are!)
  geom_vline(xintercept = p_ord8_slope, color = "red", linetype = "dashed", size = 1.5) +
  
  labs(title = "8-tone Slopes",
       x = "Slope", y = "Density") +
  xlim(-0.5, 0.5) +
  ylim(0, y_max) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold")) +
  
  # Add ROPE probability text in top left
  annotate("text", x = -0.45, y = y_max * 0.95, vjust = 1, hjust = 0,
           label = sprintf("Pr(ROPE) = %.2f", p8_in_rope),
           size = 4, bg = "white") +
  
  # Add classical slope value near red line
  annotate("text", x = p_ord8_slope, y = y_max * 0.8, vjust = 1, hjust = -0.1,
           label = bquote(r[c] == .(sprintf("%.2f", p_ord8_slope))),
           size = 4, color = "red", bg = "white")

# Combine plots
combined_plot <- grid.arrange(p1, p2, ncol = 2, 
                             top = "Subjects slopes vs 'classical ordered'")

# %% Save Results
# first save the two plot p1,p2 separated as svg,pdf,eps
ggsave(paste0(out_dir, "slopes_4tone.png"), p1, 
       width = 6, height = 6, dpi = 300)
ggsave(paste0(out_dir, "slopes_4tone.pdf"), p1, 
       width = 6, height = 6)
ggsave(paste0(out_dir, "slopes_4tone.svg"), p1, 
       width = 6, height = 6, dpi = 300)
ggsave(paste0(out_dir, "slopes_4tone.eps"), p1, 
       width = 6, height = 6)

ggsave(paste0(out_dir, "slopes_8tone.png"), p2, 
       width = 6, height = 6, dpi = 300)
ggsave(paste0(out_dir, "slopes_8tone.pdf"), p2, 
       width = 6, height = 6)
ggsave(paste0(out_dir, "slopes_8tone.svg"), p2, 
       width = 6, height = 6, dpi = 300)
ggsave(paste0(out_dir, "slopes_8tone.eps"), p2, 
       width = 6, height = 6)

# Save the plot showing subjects cluster around zero, not theory
ggsave(paste0(out_dir, "slopes_ropehdi.png"), combined_plot, 
       width = 12, height = 6, dpi = 300)
ggsave(paste0(out_dir, "slopes_ropehdi.pdf"), combined_plot, 
       width = 12, height = 6)
ggsave(paste0(out_dir, "slopes_ropehdi.svg"), combined_plot, 
       width = 12, height = 6, dpi = 300)
ggsave(paste0(out_dir, "slopes_ropehdi.eps"), combined_plot, 
       width = 12, height = 6)
 
