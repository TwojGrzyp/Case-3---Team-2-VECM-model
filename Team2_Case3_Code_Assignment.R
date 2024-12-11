setwd("C:/Data/Szkoła/Magisterka/Semestr_1/Time_series_econometrics/Assignment")
rm(list=ls()) # removing everything in your working memory

# Install and load required packages
install.packages("statespacer")
install.packages("data.table")
install.packages("urca")
library(statespacer)
library(data.table)
library(urca)
library(ggplot2)
library(bootUR)

#####
# 1. Data Exploration and plotting a graph
# Load and Inspect Data
data("FedYieldCurve")
summary(FedYieldCurve)

# Rename columns to standard format
colnames(FedYieldCurve) <- c("Date", "M3", "M6", "Y1", "Y2", "Y3", "Y5", "Y7", "Y10")

# Ensure numeric data for plotting
FedYieldCurve$M3 <- as.numeric(FedYieldCurve$M3)
FedYieldCurve$M6 <- as.numeric(FedYieldCurve$M6)
FedYieldCurve$Y1 <- as.numeric(FedYieldCurve$Y1)
FedYieldCurve$Y2 <- as.numeric(FedYieldCurve$Y2)
FedYieldCurve$Y2 <- as.numeric(FedYieldCurve$Y3)
FedYieldCurve$Y5 <- as.numeric(FedYieldCurve$Y5)
FedYieldCurve$Y2 <- as.numeric(FedYieldCurve$Y7)
FedYieldCurve$Y10 <- as.numeric(FedYieldCurve$Y10)

# Visualize Data
ggplot(FedYieldCurve, aes(x = Date)) +
  geom_line(aes(y = M3, color = "M3")) +
  geom_line(aes(y = M6, color = "M6")) +
  geom_line(aes(y = Y1, color = "Y1")) +
  geom_line(aes(y = Y2, color = "Y2")) +
  geom_line(aes(y = Y3, color = "Y3")) +
  geom_line(aes(y = Y5, color = "Y5")) +
  geom_line(aes(y = Y7, color = "Y7")) +
  geom_line(aes(y = Y10, color = "Y10")) +
  labs(title = "Yield Curves by Maturity", x = "Time", y = "Interest Rate") +
  theme_minimal()


####
# Code based on Smeekes and Wilms, 2023
# Function for the Pantula principle with clear hypothesis and order determination
pantula_principle <- function(data, min_lag = 0, max_lag = NULL, criterion = "AIC") {
  if (NCOL(data) > 1) {
    stop("Multiple time series not allowed. Provide a univariate time series.")
  }
  
  if (is.null(max_lag)) {
    n <- length(data)
    max_lag <- round(12 * (n / 100)^(1 / 4))
  }
  
  trend_options <- c("none", "intercept", "trend")
  results <- list()
  
  for (trend in trend_options) {
    test_result <- adf(
      data = data,
      deterministics = trend,
      min_lag = min_lag,
      max_lag = max_lag,
      criterion = criterion
    )
    results[[trend]] <- list(
      statistic = test_result$statistic,
      p_value = test_result$p.value
    )
  }
  
  # Results structures for reporting
  results_df <- data.frame(
    Deterministic = trend_options,
    Statistic = sapply(results, function(x) x$statistic),
    P_Value = sapply(results, function(x) x$p_value)
  )
  
  # Sequential testing logic based on Pantula principle
  rejected_trends <- which(results_df$P_Value < 0.05)
  order_of_integration <- NA
  
  if (length(rejected_trends) == 0) {
    decision <- "Fail to reject the null hypothesis of unit root under all trend assumptions."
    order_of_integration <- "I(1) or higher"
  } else {
    first_rejected <- rejected_trends[1]
    decision <- paste0(
      "Reject the null hypothesis of unit root under '",
      trend_options[first_rejected], "' deterministic assumption."
    )
    order_of_integration <- ifelse(first_rejected == 1, "I(0)", "I(1)")
  }
  
  return(list(
    Results = results_df,
    Decision = decision,
    Order = order_of_integration
  ))
}

# FedYieldCurve time series
series_names <- c("M3", "M6", "Y1", "Y2", "Y3", "Y5", "Y7", "Y10")
final_results <- data.frame(
  Series = character(),
  Deterministic = character(),
  Statistic = numeric(),
  P_Value = numeric(),
  Order = character()
)

for (series in series_names) {
  example_data <- FedYieldCurve[[series]]
  pantula_results <- pantula_principle(example_data)
  
  series_results <- pantula_results$Results
  series_results$Series <- series
  series_results$Order <- pantula_results$Order
  
  final_results <- rbind(final_results, series_results)
}

# Final table with results and order of integration
print(final_results)

###Bootstrap part
# Pantula Principle using Bootstrap ADF with AWB
apply_pantula_bootstrap <- function(series, bootstrap = "AWB", B = 999) {
  # Perform bootstrap ADF tests with different deterministic components
  adf_boot_none <- boot_adf(data = series, deterministics = "none", bootstrap = bootstrap, B = B)
  adf_boot_drift <- boot_adf(data = series, deterministics = "intercept", bootstrap = bootstrap, B = B)
  adf_boot_trend <- boot_adf(data = series, deterministics = "trend", bootstrap = bootstrap, B = B)
  
  #Test statistics, p-values, and critical values
  bootstrap_results <- data.table(
    Model = c("None", "Drift", "Trend"),
    TestStatistic = c(adf_boot_none$statistic, adf_boot_drift$statistic, adf_boot_trend$statistic),
    PValue = c(adf_boot_none$p.value, adf_boot_drift$p.value, adf_boot_trend$p.value),
    CriticalValue1pct = c(adf_boot_none$crit.val[1], adf_boot_drift$crit.val[1], adf_boot_trend$crit.val[1]),
    CriticalValue5pct = c(adf_boot_none$crit.val[2], adf_boot_drift$crit.val[2], adf_boot_trend$crit.val[2]),
    CriticalValue10pct = c(adf_boot_none$crit.val[3], adf_boot_drift$crit.val[3], adf_boot_trend$crit.val[3])
  )
  
  # Identify the first model where the null is rejected (p-value < 0.05)
  bootstrap_results[, RejectNull := PValue < 0.05]
  first_reject_bootstrap <- bootstrap_results[RejectNull == TRUE, .(Model)][1, Model]
  
  # Add the selected model to the results
  return(data.table(SelectedModel = ifelse(is.na(first_reject_bootstrap), "No Rejection", first_reject_bootstrap), bootstrap_results))
}

# Apply the bootstrap function to each series and store results
bootstrap_results_list <- lapply(colnames(FedYieldCurve), function(col) {
  series_bootstrap_results <- apply_pantula_bootstrap(FedYieldCurve[, col]) # Use the bootstrap function
  series_bootstrap_results[, Series := col]  # Add series name to results
  return(series_bootstrap_results)
})

# Combine all bootstrap results into one table
final_bootstrap_results <- rbindlist(bootstrap_results_list)
final_bootstrap_results <- final_bootstrap_results[-c(1:3), ]


# Display bootstrap results
print(final_bootstrap_results)



# Apply the Bootstrap Union Test to each series in the FedYieldCurve dataset
apply_union_bootstrap <- function(series, bootstrap = "MBB", B = 999) {
  # Perform the Union of Rejections Test
  union_results <- boot_union(data = series, bootstrap = bootstrap, B = B)
  
  # Extract test statistics, p-value, and critical values
  union_test_results <- data.table(
    Model = "Union of Rejections",
    TestStatistic = union_results$statistic,
    PValue = union_results$p.value,
    CriticalValue1pct = union_results$crit.val[1],
    CriticalValue5pct = union_results$crit.val[2],
    CriticalValue10pct = union_results$crit.val[3],
    RejectNull = union_results$p.value < 0.05
  )
  
  return(union_test_results)
}

###BOOTSTRAP UNION OF REJECTION TEST

# Bootstrap Union Test for each series in the FedYieldCurve dataset
apply_union_bootstrap <- function(series, bootstrap = "MBB", B = 999) {
  # Perform the Union of Rejections Test
  union_results <- boot_union(data = series, bootstrap = bootstrap, B = B)
  
  # Extract test statistics, p-value, and critical values
  union_test_results <- data.table(
    Model = "Union of Rejections",
    TestStatistic = union_results$statistic,
    PValue = union_results$p.value,
    CriticalValue1pct = union_results$crit.val[1],
    CriticalValue5pct = union_results$crit.val[2],
    CriticalValue10pct = union_results$crit.val[3],
    RejectNull = union_results$p.value < 0.05
  )
  
  return(union_test_results)
}

# Apply the function to all series
union_results_list <- lapply(colnames(FedYieldCurve), function(col) {
  series_union_results <- apply_union_bootstrap(FedYieldCurve[, col])
  series_union_results[, Series := col]  # Add series name
  return(series_union_results)
})

# Combining results into a single table
final_union_results <- rbindlist(union_results_list)

# Display results
print(final_union_results)

###DETERMINING THE COINTEGRATION RANK

##TRACE TEST
##Short-term Yields

FedYieldCurve <- FedYieldCurve[, sapply(FedYieldCurve, is.numeric)]  
FedYieldCurve <- na.omit(FedYieldCurve)  
ts_data <- ts(FedYieldCurve, start = c(2000, 1), frequency = 12) 

# Assume you have already selected the appropriate I(1) series under the "Trend" assumption
selected_series <- c("M3", "M6", "Y1")  # Adjust as needed
ts_data_filtered <- ts_data[, selected_series]  # Subset the selected series
ts_data_filtered <- na.omit(as.matrix(ts_data_filtered))  # Ensure clean numeric data

# Step 2: Determine optimal lag length
library(vars)
lag_selection <- VARselect(ts_data_filtered, lag.max = 10, type = "const")
optimal_lag <- lag_selection$selection["AIC(n)"]  

# Step 3: Perform the Johansen Test with the "Trend" assumption
library(urca)
johansen_result_trace_short <- ca.jo(ts_data_filtered, type = "trace", ecdet = "trend", K = optimal_lag)

# Step 4: Summarize results
summary(johansen_result_trace_short)

# Step 5: Extract results into a table
# Extract test statistics and critical values
trace_stats <- johansen_result_trace_short@teststat
critical_values <- johansen_result_trace_short@cval

# Create a summary table
hypotheses <- paste("r <=", 0:(length(trace_stats) - 1))  # Hypotheses tested
results_table_trace_short <- data.frame(
  Hypothesis = hypotheses,
  TraceStatistic = trace_stats,
  CriticalValue1pct = critical_values[, 1],  # 1% critical value
  CriticalValue5pct = critical_values[, 2],  # 5% critical value
  CriticalValue10pct = critical_values[, 3], # 10% critical value
  Decision = ifelse(trace_stats > critical_values[, 2], "Reject H0", "Fail to Reject H0")
)

# Add a column for the cointegration rank decision
results_table_trace_short$CointegrationRank <- NA
for (j in 1:nrow(results_table_trace_short)) {
  if (results_table_trace_short$Decision[j] == "Fail to Reject H0") {
    results_table_trace_short$CointegrationRank <- j - 1  # Cointegration rank is the first 'Fail to Reject H0'
    break
  }
}

# If no 'Fail to Reject H0', set the rank to max possible
if (is.na(results_table_trace_short$CointegrationRank[nrow(results_table_trace_short)])) {
  results_table_trace_short$CointegrationRank[nrow(results_table_trace_short)] <- nrow(results_table_trace_short) - 1
}

# Print the summary table
print(results_table_trace_short)

#####MAXIMUM EIGENVALUE TEST
##Short-term Yields

# Step 1: Prepare the data
# Assume you have already selected the appropriate I(1) series under the "Trend" assumption
selected_series <- c("M3", "M6", "Y1")  # Adjust as needed
ts_data_filtered <- ts_data[, selected_series]  # Subset the selected series
ts_data_filtered <- na.omit(as.matrix(ts_data_filtered))  # Ensure clean numeric data

# Step 2: Determine optimal lag length
library(vars)
lag_selection <- VARselect(ts_data_filtered, lag.max = 10, type = "const")
optimal_lag <- lag_selection$selection["AIC(n)"]  # Use AIC to select lag length

# Step 3: Perform the Johansen Test with the "Trend" assumption
library(urca)
johansen_result <- ca.jo(ts_data_filtered, type = "eigen", ecdet = "trend", K = optimal_lag)

# Step 4: Summarize results
summary(johansen_result)

# Step 5: Extract results into a table
# Extract test statistics and critical values
eigen_stats <- johansen_result@teststat
critical_values <- johansen_result@cval

# Create a summary table
hypotheses <- paste("r <=", 0:(length(eigen_stats) - 1))  # Hypotheses tested
results_table_eigen_short <- data.frame(
  Hypothesis = hypotheses,
  MaxEigenStatistic = eigen_stats,
  CriticalValue1pct = critical_values[, 1],  # 1% critical value
  CriticalValue5pct = critical_values[, 2],  # 5% critical value
  CriticalValue10pct = critical_values[, 3], # 10% critical value
  Decision = ifelse(eigen_stats > critical_values[, 2], "Reject H0", "Fail to Reject H0")
)

# Add a column for the cointegration rank decision
results_table_eigen_short$CointegrationRank <- NA
for (j in 1:nrow(results_table_eigen_short)) {
  if (results_table_eigen_short$Decision[j] == "Fail to Reject H0") {
    results_table_eigen_short$CointegrationRank <- j - 1  # Cointegration rank is the first 'Fail to Reject H0'
    break
  }
}

# If no 'Fail to Reject H0', set the rank to max possible
if (is.na(results_table_eigen_short$CointegrationRank[nrow(results_table_eigen_short)])) {
  results_table_eigen_short$CointegrationRank[nrow(results_table_eigen_short)] <- nrow(results_table_eigen_short) - 1
}

# Print the summary table
print(results_table_eigen_short)


####TRACE TEST
##Long-term Yields

# Step 1: Prepare the data
# Assume you have already selected the appropriate I(1) series under the "Trend" assumption
selected_series <- c("Y2", "Y3", "Y5")  # Adjust as needed
ts_data_filtered <- ts_data[, selected_series]  # Subset the selected series
ts_data_filtered <- na.omit(as.matrix(ts_data_filtered))  # Ensure clean numeric data

# Step 2: Determine optimal lag length
library(vars)
lag_selection <- VARselect(ts_data_filtered, lag.max = 10, type = "const")
optimal_lag <- lag_selection$selection["AIC(n)"]  

# Step 3: Perform the Johansen Test with the "Trend" assumption
library(urca)
johansen_result_trace_long <- ca.jo(ts_data_filtered, type = "trace", ecdet = "trend", K = optimal_lag)

# Step 4: Summarize results
summary(johansen_result_trace_long)

# Step 5: Extract results into a table
# Extract test statistics and critical values
trace_stats <- johansen_result_trace_long@teststat
critical_values <- johansen_result_trace_long@cval

# Create a summary table
hypotheses <- paste("r <=", 0:(length(trace_stats) - 1))  
results_table_trace_long <- data.frame(
  Hypothesis = hypotheses,
  TraceStatistic = trace_stats,
  CriticalValue1pct = critical_values[, 1],  # 1% critical value
  CriticalValue5pct = critical_values[, 2],  # 5% critical value
  CriticalValue10pct = critical_values[, 3], # 10% critical value
  Decision = ifelse(trace_stats > critical_values[, 2], "Reject H0", "Fail to Reject H0")
)

# Add a column for the cointegration rank decision
results_table_trace_long$CointegrationRank <- NA
for (j in 1:nrow(results_table_trace_long)) {
  if (results_table_trace_long$Decision[j] == "Fail to Reject H0") {
    results_table_trace_long$CointegrationRank <- j - 1  # Cointegration rank is the first 'Fail to Reject H0'
    break
  }
}

# If no 'Fail to Reject H0', set the rank to max possible
if (is.na(results_table_trace_long$CointegrationRank[nrow(results_table_trace_long)])) {
  results_table_trace_long$CointegrationRank[nrow(results_table_trace_long)] <- nrow(results_table_trace_long) - 1
}

# Print the summary table
print(results_table_trace_long)

#####MAXIMUM EIGENVALUE TEST
##Long-term Yields

# Step 1: Prepare the data
# Assume you have already selected the appropriate I(1) series under the "Trend" assumption
selected_series <- c("Y2", "Y3", "Y5")  # Adjust as needed
ts_data_filtered <- ts_data[, selected_series]  # Subset the selected series
ts_data_filtered <- na.omit(as.matrix(ts_data_filtered))  # Ensure clean numeric data

# Step 2: Determine optimal lag length
library(vars)
lag_selection <- VARselect(ts_data_filtered, lag.max = 10, type = "const")
optimal_lag <- lag_selection$selection["AIC(n)"]  # Use AIC to select lag length

# Step 3: Perform the Johansen Test with the "Trend" assumption
library(urca)
johansen_result <- ca.jo(ts_data_filtered, type = "eigen", ecdet = "trend", K = optimal_lag)

# Step 4: Summarize results
summary(johansen_result)

# Step 5: Extract results into a table
# Extract test statistics and critical values
eigen_stats <- johansen_result@teststat
critical_values <- johansen_result@cval

# Create a summary table
hypotheses <- paste("r <=", 0:(length(eigen_stats) - 1))  # Hypotheses tested
results_table_eigen <- data.frame(
  Hypothesis = hypotheses,
  MaxEigenStatistic = eigen_stats,
  CriticalValue1pct = critical_values[, 1],  # 1% critical value
  CriticalValue5pct = critical_values[, 2],  # 5% critical value
  CriticalValue10pct = critical_values[, 3], # 10% critical value
  Decision = ifelse(eigen_stats > critical_values[, 2], "Reject H0", "Fail to Reject H0")
)

# Add a column for the cointegration rank decision
results_table_eigen$CointegrationRank <- NA
for (j in 1:nrow(results_table_eigen)) {
  if (results_table_eigen$Decision[j] == "Fail to Reject H0") {
    results_table_eigen$CointegrationRank <- j - 1  # Cointegration rank is the first 'Fail to Reject H0'
    break
  }
}

# If no 'Fail to Reject H0', set the rank to max possible
if (is.na(results_table_eigen$CointegrationRank[nrow(results_table_eigen)])) {
  results_table_eigen$CointegrationRank[nrow(results_table_eigen)] <- nrow(results_table_eigen) - 1
}

# Print the summary table
print(results_table_eigen)

# Perform Johansen Cointegration Test
library(urca)

cointegration_data <- FedYieldCurve[, sapply(FedYieldCurve, is.numeric)]

# Remove rows with missing values
cointegration_data <- na.omit(cointegration_data)

# Convert to a numeric matrix
cointegration_data <- as.matrix(cointegration_data)

# Print the resulting matrix
print(cointegration_data)

#Ensure the data is clean and numeric
cointegration_data <- na.omit(cointegration_data)  # Remove missing values
cointegration_data <- as.matrix(cointegration_data)

#Check for multicollinearity and reduce dimensions if necessary
cor_matrix <- cor(cointegration_data)
if (any(abs(cor_matrix[upper.tri(cor_matrix)]) > 0.95)) {
  cat("High multicollinearity detected. Applying PCA...\n")
  pca_result <- prcomp(cointegration_data, scale. = TRUE)
  cointegration_data <- pca_result$x[, 1:5]  # Use first 5 PCs (adjust as needed)
}

#Check dimensions and reduce variables/lag length if necessary
n <- ncol(cointegration_data)
T <- nrow(cointegration_data)
if (T <= n * 2 + 1) {
  cat("Not enough observations for lag length. Reducing variables...\n")
  cointegration_data <- cointegration_data[, 1:5]  # Use fewer variables
}

#Johansen Test
library(urca)
johansen_test <- ca.jo(cointegration_data, type = "trace", ecdet = "trend", K = 2)
summary(johansen_test)

# Determine the rank (number of cointegrating relationships)
cointegration_rank <- sum(johansen_test@teststat > johansen_test@cval[, 2])  # Compare with 5% critical values
cat("\nIdentified Cointegration Rank:", cointegration_rank, "\n")

# Perform VECM if cointegration is found
if (cointegration_rank > 0) {
  library(vars)
  
  # Fit the VECM using the number of cointegrating vectors identified in Johansen test
  vecm_model <- cajorls(johansen_test, r = cointegration_rank)  # r is the cointegration rank
  cat("\nVECM Model Results:\n")
  print(vecm_model)
  
  # Analyze the VECM results
  cat("\nLong-Run Relationships (Beta):\n")
  print(vecm_model$beta)  # Long-run relationships (cointegration vectors)
  
  cat("\nAdjustment Coefficients (Alpha):\n")
  print(vecm_model$alpha)  # Speed of adjustment to equilibrium
} else {
  cat("\nNo cointegration relationships were identified.\n")
}

##VECM model

# Use the results of the Johansen test
vecm_model <- cajorls(johansen_result_trace_short)

# Summarize the results
summary(vecm_model)

# Extract the Alpha (Adjustment Coefficients) Matrix
alpha_matrix <- vecm_model$alpha
cat("\nAdjustment Coefficients (Alpha):\n")
print(alpha_matrix)

# Extract the Long-Run Cointegration Coefficients (Beta)
beta_matrix <- vecm_model$beta
cat("Long-Run Cointegration Coefficients (Beta):\n")
print(beta_matrix)

# Extract the Short-Run Dynamics Coefficients (Gamma)
gamma_matrix <- vecm_model$rlm$coefficients
cat("\nShort-Run Dynamics Coefficients (Gamma):\n")
print(gamma_matrix)

# Extract all short-run dynamics coefficients (Gamma)
gamma_matrix <- vecm_model$rlm$coefficients[-c(1, 2), ]  # Remove ECT1 and constant rows

# Filter only dl1 coefficients (first lag of differenced terms)
gamma_coefficients <- gamma_matrix[grep("\\.dl1$", rownames(gamma_matrix)), ]
cat("\nFirst Lag (dl1) Coefficients:\n")
print(gamma_coefficients)


#IRF graphs
# Convert the VECM into a VAR for IRF analysis
library(vars)
vecm_as_var <- vec2var(johansen_result_trace_short)

# Generate impulse response functions
irf_results_M3 <- irf(vecm_as_var, impulse = "M3", response = "M6", n.ahead = 10, boot = TRUE)

# Plot the IRF
plot(irf_results_M3)

# Generate impulse response functions
irf_results_M6 <- irf(vecm_as_var, impulse = "M6", response = "M3", n.ahead = 10, boot = TRUE)

# Plot the IRF
plot(irf_results_M6)

# Generate impulse response functions
irf_results_M3_2 <- irf(vecm_as_var, impulse = "M3", response = "Y1", n.ahead = 10, boot = TRUE)

# Plot the IRF
plot(irf_results_M3_2)

# Generate impulse response functions
irf_results_Y1 <- irf(vecm_as_var, impulse = "M6", response = "Y1", n.ahead = 10, boot = TRUE)

# Plot the IRF
plot(irf_results_Y1)

# Generate impulse response functions
irf_results_Y1_2 <- irf(vecm_as_var, impulse = "Y1", response = "M6", n.ahead = 10, boot = TRUE)

# Plot the IRF
plot(irf_results_Y1_2)

# Generate impulse response functions
irf_results_M3_3 <- irf(vecm_as_var, impulse = "Y1", response = "M3", n.ahead = 10, boot = TRUE)

# Plot the IRF
plot(irf_results_M3_3)

