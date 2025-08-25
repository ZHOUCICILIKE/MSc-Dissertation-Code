# --- MCMC最终修正版 v3 (调用 ookk.stan) ---

# --- 0. 加载必要的库 ---
library(rstan)
library(dplyr)
if (!require("loo")) {
  install.packages("loo")
  library(loo)
}

# 设定选项
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# --- 1. 加载并预处理你的 gait.csv 数据 ---
cat("Step 1: Loading and preprocessing gait.csv...\n")
gait_data <- read.csv("gait.csv")

# (这部分代码和之前完全一样)
gait_data$subject <- factor(gait_data$subject)
gait_data$condition <- factor(gait_data$condition, levels = 1:3, labels = c("unbraced", "knee.brace", "ankle.brace"))
gait_data$leg <- factor(gait_data$leg, levels = 1:2, labels = c("left", "right"))
gait_data$joint <- factor(gait_data$joint, levels = 1:3, labels = c("ankle", "knee", "hip"))


# --- 2. 为特定条件离散化数据 ---
cat("Step 2: Discretizing data for the selected condition...\n")
# (这部分代码和之前完全一样)
cond_to_test <- "knee.brace"
leg_to_test <- "right"
joint_to_test <- "hip"
num_states <- 10

subset_data <- gait_data %>%
  filter(condition == cond_to_test, leg == leg_to_test, joint == joint_to_test)

hist_obj <- hist(subset_data$angle, breaks = num_states, plot = FALSE)
bin_breaks <- hist_obj$breaks

all_trajectories <- list()
for (s in unique(subset_data$subject)) {
  for (r in 1:10) {
    traj <- subset_data %>%
      filter(subject == s, replication == r) %>%
      arrange(time) %>%
      pull(angle)
    
    if (length(traj) > 1) { 
      state_sequence_raw <- findInterval(traj, bin_breaks, all.inside = TRUE)
      state_sequence_corrected <- pmax(1, pmin(state_sequence_raw, num_states))
      all_trajectories[[length(all_trajectories) + 1]] <- state_sequence_corrected
    }
  }
}


# --- 3. 将数据格式化为 Stan 需要的输入 ---
cat("Step 3: Formatting data for Stan...\n")
# (这部分代码和之前完全一样)
traj_lengths <- sapply(all_trajectories, length)
fixed_length <- as.integer(median(traj_lengths))
num_trajectories <- length(all_trajectories)
y_matrix <- matrix(1, nrow = num_trajectories, ncol = fixed_length) 
for (i in 1:num_trajectories) {
  traj <- all_trajectories[[i]]
  len <- length(traj)
  if (len >= fixed_length) {
    y_matrix[i, ] <- traj[1:fixed_length]
  } else {
    y_matrix[i, ] <- c(traj, rep(traj[len], fixed_length - len))
  }
}


# --- 4. 使用 MCMC 运行 Stan 模型 ---
cat("Step 4: Running MCMC analysis...\n")

# **为 k=1 运行MCMC模型**
cat("Running MCMC analysis for K=1 (this will take some time)...\n")
stan_data_k1 <- list(K = 1, S = num_states, N = num_trajectories, T = fixed_length, y = y_matrix)
# 【修正处】使用新的文件名
fit_mcmc_k1 <- stan(
  file = "ookk.stan", 
  data = stan_data_k1, 
  chains = 4, 
  iter = 1000, 
  warmup = 500, 
  seed = 123
)

# **为 k=2 运行MCMC模型**
cat("Running MCMC analysis for K=2 (this will also take some time)...\n")
stan_data_k2 <- list(K = 2, S = num_states, N = num_trajectories, T = fixed_length, y = y_matrix)
# 【修正处】使用新的文件名
fit_mcmc_k2 <- stan(
  file = "ookk.stan", 
  data = stan_data_k2, 
  chains = 4,
  iter = 1000,
  warmup = 500,
  seed = 123
)


# --- 5. 使用 loo 包比较 MCMC 模型 ---
cat("Step 5: Comparing model fits using LOOIC and concluding...\n\n")
loo_k1 <- loo(fit_mcmc_k1, save_psis = TRUE)
loo_k2 <- loo(fit_mcmc_k2, save_psis = TRUE)

loo_comparison <- loo_compare(loo_k1, loo_k2)

cat("--- Stan MCMC Model Comparison using LOOIC ---\n")
cat("Condition tested:", cond_to_test, "-", leg_to_test, "-", joint_to_test, "\n\n")
cat("LOOIC Comparison Table:\n")
print(loo_comparison)

cat("\n--- Interpretation ---\n")
cat("The table ranks models from best to worst. The model at the top is the preferred model.\n")
cat("'elpd_diff' shows the difference in expected log predictive density. A negative value means the model is worse than the model above it.\n")
cat("'se_diff' is the standard error of that difference. If elpd_diff is more than a few times its standard error, the difference is considered meaningful.\n")