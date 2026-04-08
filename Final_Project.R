## Libraries
install.packages("changepoint")
library(changepoint)
library(tseries)
library(forecast)
library(stats)
library(readr)
library(tidyverse)
library(lmtest)
library(knitr)
library(rmarkdown)
library(strucchange)
library(segmented)
library(ggplot2)
library(dplyr)
library(tidyr)
library(uroot)
library(depmixS4)
library(gridExtra)
library(dlm)


##############################################     
## TASK #1 Data acquisition and exploration ##
##############################################

### Load Data
gistemp <- read_csv("Data/gistemp.txt")
ghcn <- read_csv("Data/ghcn.txt")

### Create time series objects
#gistemp
long_data_gistemp <- gistemp[,1:13] %>%
  pivot_longer(-Year, names_to = "Month", values_to = "Temperature") %>%
  mutate(Month = match(Month, month.abb)) %>% 
  arrange(Year, Month)

gistemp_ts <- ts(long_data_gistemp$Temperature, start = c(min(long_data_gistemp$Year), 1), frequency = 12)
plot(gistemp_ts, ylab = "Temperature", xlab = "",
     col = "darkgrey", main = "GISTEMP Data: Average Global Surface Temperature")
grid()

#ghcn
data_ghcn <- ghcn[, c(2, 6:8)] %>% 
  mutate(TMIN = TMIN / 10, TMAX = TMAX / 10) %>%
  filter(station == "SAN FRANCISCO DWTN") 

ghcn_max_ts <- ts(data_ghcn$TMAX, start = c(as.numeric(format(min(data_ghcn$date), "%Y")), 
                                            as.numeric(format(min(data_ghcn$date), "%j"))), 
                  frequency = 365)
plot(ghcn_max_ts)
ghcn_min_ts <- ts(data_ghcn$TMIN, start = c(as.numeric(format(min(data_ghcn$date), "%Y")), 
                                            as.numeric(format(min(data_ghcn$date), "%j"))), 
                  frequency = 365)
plot(ghcn_min_ts)

# plot for ghcn
plot(ghcn_max_ts, type = "l", col = "orange", ylab = "Temperature", xlab = '',
     main = "GHCN Data: Daily Min and Max Temperatures",
     ylim = c(-5, 45))
lines(ghcn_min_ts, col = "navy")
legend("topright", legend = c("Minimum Temperature", "Maximum Temperature"), col = c("navy", "orange"), lty = 1)
grid(col = "lightgrey", lty = "dotted", lwd = 0.8)

### Time Series Decomposition: use STL to allow for changes in seasonality
# for GISTEMP

# seasonality changes over 2 years
decomp_stl_gistemp <- stl(gistemp_ts, s.window = 27) 
plot(decomp_stl_gistemp, col = "darkgrey",
     main = "STL Decomposition of GISTEMP Data")

# for GHCN
# seasonality changes over 1 year
decomp_stl_ghcn_max <- stl(ghcn_max_ts, s.window = 365) 
plot(decomp_stl_ghcn_max, main = "STL Decomposition of Max Temperatures", col = 'orange')

decomp_stl_ghcn_min <- stl(ghcn_min_ts, s.window = 365) 
plot(decomp_stl_ghcn_min, main = "STL Decomposition of Min Temperatures", col = 'navy')


###########################################################     
## TASK #2 GISTEMP: Do the data document global warming? ##
###########################################################

### Obtain the seasonally adjusted time series

adj_data <- seasadj(decomp_stl_gistemp)



##################################     
## Step 1. Hidden Markov Models ##
##################################

# Clean data
data_clean <- na.omit(adj_data)
time_idx <- as.numeric(time(gistemp_ts))
hmm_data <- data.frame(temp = as.numeric(data_clean), time = time_idx)

##################################################################################################
##### a state-dependent level: Yt | St = j ∼ N (μj , σ2j ) 

mod_level <- list()
post_states <- list()

for (i in 2:5) {
  mod_level[[i]] <- depmix(temp ~ 1, data = hmm_data, nstates = i, family = gaussian())
  fit_level <- fit(mod_level[[i]])
  post_states[[i]] <- posterior(fit_level) 
}


# Plot with colors
par(mfrow = c(2, 2),        
    oma = c(1, 1, 1, 1),  
    mar = c(2.5, 4, 2, 2))  

for (i in 2:5) {
  
  col_index <- (i - 2) %% 2 + 1 
  
  plot(hmm_data$time, hmm_data$temp, type = 'l',
       main = paste("Level-only HMM with", i, "States"),
       ylab = if (col_index == 1) "Temperature" else "",
       xlab = "", xaxt = "n", yaxt = if (col_index == 1) "s" else "n")
  grid()
  
  x_ticks <- seq(floor(min(hmm_data$time)), floor(max(hmm_data$time)), by = 20)
  axis(1, at = x_ticks, labels = FALSE)  # hide default labels
  text(x = x_ticks, y = par("usr")[3] - 0.15, labels = x_ticks,
       srt = 45, xpd = TRUE, adj = 1, cex = 0.9)
  
  # Add points colored by state
  points(hmm_data$time, hmm_data$temp,
         col = rainbow(i)[post_states[[i]]$state],
         pch = 20, cex = 0.6)
}

par(mfrow = c(1, 1))


plots_level <- list()

for (i in 2:4) {
  print(summary(fit(mod_level[[i]])))  # print model summary
  
  decoded_states_level_i <- post_states[[i]]$state
  time_index <- hmm_data$time
  state_df_level_i <- data.frame(time_index = time_index, state = decoded_states_level_i)
  
  plots_level[[i]] <- ggplot(state_df_level_i, aes(x = time_index, y = state)) +
    geom_step() +
    scale_y_continuous(breaks = 1:max(decoded_states_level_i)) +
    theme_minimal() +
    labs(title = paste(i, "States of Level-only HMM"),
         x = "Year", y = "State")
  
  print(plots_level[[i]])  # display the plot in each iteration
}
grid.arrange(plot21, plot31, plot41, ncol = 1)




##################################################################################################
# or with a state-dependent linear trend  Yt | St = j ∼ N (αj + βj t, σ2

# Model with linear trend (time-dependent means)
mod_trend <- list()
post_states_trend <- list()

for (i in 2:5) {
  mod_trend[[i]] <- depmix(temp ~ time, data = hmm_data, nstates = i, family = gaussian())
  fit_trend <- fit(mod_trend[[i]])
  mod_trend[[i]] <- fit_trend 
  post_states_trend[[i]] <- posterior(fit_trend)
}


par(mfrow = c(2, 2),        
    oma = c(1, 1, 1, 1),    
    mar = c(2.5, 4, 2, 2))  

for (i in 2:5) {
  col_index <- (i - 2) %% 2 + 1 
  
  plot(hmm_data$time, hmm_data$temp, type = 'l',
       main = paste("Linear Trend HMM with", i, "States"),
       ylab = if (col_index == 1) "Temperature" else "",
       xlab = "", xaxt = "n", yaxt = if (col_index == 1) "s" else "n")
  grid()
  
  x_ticks <- seq(floor(min(hmm_data$time)), floor(max(hmm_data$time)), by = 20)
  axis(1, at = x_ticks, labels = FALSE)
  text(x = x_ticks,
       y = par("usr")[3] - 0.15,
       labels = x_ticks,
       srt = 45, xpd = TRUE, adj = 1, cex = 0.9)
  
  points(hmm_data$time, hmm_data$temp,
         col = rainbow(i)[post_states_trend[[i]]$state],
         pch = 20, cex = 0.6)
}
par(mfrow = c(1, 1))

# Summary 
summary(fit(mod_trend[[2]])) 

# Plot for 2-state TREND model
decoded_states_trend_2 <- post_states_trend[[2]]$state
state_df_trend_2 <- data.frame(time_index = time_index, state = decoded_states_trend_2)

ggplot(state_df_trend_2, aes(x = time_index, y = state)) +
  geom_step() +
  scale_y_continuous(breaks = 1:max(decoded_states_trend_2)) +
  theme_minimal() +
  labs(title = "2 States of Linear Trend HMM",
       x = "Year", y = "State")


### BIC/AIC criteria

# Level-Only Models:
n_states <- c()
aic_vals <- c()
bic_vals <- c()

for (i in 2:5) {
  fit_model <- mod_level[[i]]  
  n_states <- c(n_states, i)
  aic_vals <- c(aic_vals, AIC(fit_model))
  bic_vals <- c(bic_vals, BIC(fit_model))
}

model_comparison <- data.frame(
  States = n_states,
  AIC = round(aic_vals, 2),
  BIC = round(bic_vals, 2)
)
print(model_comparison)


# Trend Models: 
n_states_trend <- c()
aic_vals_trend <- c()
bic_vals_trend <- c()

for (i in 2:5) {
  fit_model <- mod_trend[[i]]  
  n_states_trend <- c(n_states_trend, i)
  aic_vals_trend <- c(aic_vals_trend, AIC(fit_model))
  bic_vals_trend <- c(bic_vals_trend, BIC(fit_model))
}

trend_model_comparison <- data.frame(
  States = n_states_trend,
  AIC = round(aic_vals_trend, 2),
  BIC = round(bic_vals_trend, 2)
)

print(trend_model_comparison)

## Comparison of the Models:
combined_summary <- merge(
  model_comparison, trend_model_comparison,
  by = "States", suffixes = c("_Level", "_Trend")
)
print(combined_summary)






##################################     
## Step 2. Dynamic Linear Models #
##################################

##################################     
## Randow Walk plus Noise 
################################## 

dlm_data <- adj_data

build_rw <- function(param){
  dlmModPoly(order=1, dV=param[1], dW=param[2], m0= dlm_data[1])}
MLE_rw <- dlmMLE(dlm_data, parm = rep(1740, 2), build_rw, lower=c(0.00001, 0))

rw_model <- build_rw(MLE_rw$par)

#filtering
filt_rw <- dlmFilter(dlm_data,rw_model) 

# smoothing 
smoothed_rw <- dlmSmooth(filt_rw)

smoothed_mean_rw <- drop(smoothed_rw$s[-1])
smoothed_sd_rw <- sqrt(sapply(dlmSvd2var(smoothed_rw$U.S, smoothed_rw$D.S)[-1], function(x) x[1, 1]))

ts.plot(dlm_data,
        smoothed_mean_rw,
        smoothed_mean_rw + 1.96 * smoothed_sd_rw,
        smoothed_mean_rw - 1.96 * smoothed_sd_rw,
        col = c("lightgrey", "blue", "red", "red"),
        lty = c(1, 1, 2, 2),
        lwd = c(1, 2, 1, 1),
        ylab = "Temperature",
        main = "Smoothed State Estimates with 95% Confidence Interval")
abline(v = 1965, col = "black", lty = 3)
axis(1, at = 1965, labels = "1965")  # add custom x-axis tick at 1970
legend("topleft",
       legend = c("Observed Temperature", "Smoothed Estimate", "95% CI"),
       col = c("lightgrey", "blue", "red"),
       lty = c(1, 1, 2),
       lwd = c(1, 2, 1),
       bg = "white",
       box.lty = 0)

### Kalman Smoother for structural breaks

smoothed_level <- dropFirst(smoothed_rw$s)
smoother_innov <- diff(smoothed_level) 
smoother_innov_ts <- ts(smoother_innov, start = 2, frequency = frequency(dlm_data))
#plot(smoother_innov_ts, type = "p", col = "blue")
#abline(h = 0, col = "lightgrey", lty = 2)
  # basically, no visible change points
### Change points - formally - not found
cpt <- cpt.mean(smoother_innov, method = "PELT")
change_points <- cpts(cpt)
print(change_points)
time_points <- time(smoother_innov_ts)[change_points]
print(time_points) 

### Model Checking
# normality of residuals

res_rw <- residuals(filt_rw, sd = FALSE) 

#QQ-plot
qqnorm(res_rw)
qqline(res_rw)

# residuals distribution
hist(res_rw) 

# Other diagnostic plots
tsdiag(filt_rw)

#Shapiro test
shapiro.test(res_rw)

# lack of serial correlation
# Box-Ljung test
Box.test(res_rw, lag = 20, type = "Ljung")

##################################     
## Locally Linear Trend Model 
##################################

build_llt <- function(param){
  dlmModPoly(order=2, dV=param[1], dW=param[2:3], m0= c(dlm_data[1], 0))} 
#assuming trend starts at 0

MLE_llt <- dlmMLE(dlm_data, parm = rep(1740, 3), build_llt, lower=c(0.00001, 0, 0))

llt_model <- build_llt(MLE_llt$par)

#filtering
filt_llt <- dlmFilter(dlm_data,llt_model) 

filterEstimates_llt= window(filt_llt$m,start=start(dlm_data)[1]) 
filterEstimates_llt_theta <- filterEstimates_llt[,1]
filterEstimates_llt_beta <- filterEstimates_llt[,2]

### beta
# if we throw out first 100 "large" estimates, 
# it is easy to see that the level grows from year to year (by tiny percentage)


par(mfrow = c(2, 1), mar = c(3, 4, 3, 2))

start_index <- 1
start_year <- 1880 + (start_index - 1) %/% 12
start_month <- ((start_index - 1) %% 12) + 1

filterEstimates_llt_beta_ts <- ts(filterEstimates_llt_beta[1:150],
                                  start = c(start_year, start_month),
                                  frequency = 12)

plot(filterEstimates_llt_beta_ts, type = "l", col = "blue", lwd = 2, xlab = "",
     ylab = "", main = "Estimates of β for 1880–1892")
abline(h = 0, lty = 3, col = "black")

start_index <- 150
start_year <- 1880 + (start_index - 1) %/% 12
start_month <- ((start_index - 1) %% 12) + 1

filterEstimates_llt_beta_ts <- ts(filterEstimates_llt_beta[150:1740],
                                  start = c(start_year, start_month),
                                  frequency = 12)

plot(filterEstimates_llt_beta_ts, type = "l", col = "blue", lwd = 2, xlab = "",
     ylab = "", main = "Estimates of β for 1892–2025")
abline(h = 0, lty = 3, col = "black")

par(mfrow = c(1, 1))


# smoothing
smoothed_llt <- dlmSmooth(filt_llt)
smoothedEstimates_llt= window(smoothed_llt$s,start=start(dlm_data)[1]) 

smoothedEstimates_llt_theta <- smoothedEstimates_llt[,1]
smoothedEstimates_llt_beta <- smoothedEstimates_llt[,2] # beta is constant, its sd is 0

smoothed_mean_llt <- drop(smoothed_llt$s[-1,1])
smoothed_sd_llt <- sqrt(sapply(dlmSvd2var(smoothed_llt$U.S, smoothed_llt$D.S)[-1], function(x) x[1, 1]))

ts.plot(dlm_data,
        smoothed_mean_llt,
        smoothed_mean_llt + 1.96 * smoothed_sd_llt,
        smoothed_mean_llt - 1.96 * smoothed_sd_llt,
        lty = c(1, 1, 3, 3),
        col = c("lightgrey", "blue", "red", "red"),
        lwd = c(1, 2, 1, 1),
        ylab = "Temperature",
        main = "Smoothed State Estimates with 95% CI")

### Change point detection - not found
smoothed_level <- dropFirst(smoothed_llt$s)[, 1]
smoother_innov <- diff(smoothed_level)
smoother_innov_ts <- ts(smoother_innov, start = 2, frequency = frequency(dlm_data))
cpt <- cpt.mean(smoother_innov, method = "PELT")
change_points <- cpts(cpt)
time_points <- time(smoother_innov_ts)[change_points]
print(time_points)

### model checking
# residuals distribution
res_llt <- residuals(filt_llt, sd = FALSE) 

hist(res_llt) 

#QQ-plot
qqnorm(res_llt)
qqline(res_llt)

# Other diagnostic plots
tsdiag(filt_llt)

#Shapiro test

shapiro.test(res_llt)

# lack of serial correlation
# Box-Ljung test
Box.test(res_llt, lag = 20, type = "Ljung")








###########################################################     
## TASK #3 GHCN: Weather preditcion.
###########################################################

###### DATA PREPARATION ######
decomp_stl_ghcn_max <- stl(ghcn_max_ts, s.window = 365) 
adjusted_max <- seasadj(decomp_stl_ghcn_max)
decomp_stl_ghcn_min <- stl(ghcn_min_ts, s.window = 365) 
adjusted_min <- seasadj(decomp_stl_ghcn_min)

ghcn_ts <- cbind(adjusted_min, adjusted_max)

##################################     
## Step 2. Model A 
##################################


# Step 1: Use empirical variances of TS for initial guesses
# we would initially set relativiely small signal-to-noise ratio
obs_var <- apply(ghcn_ts, 2, var)
init_dV <- obs_var / 10                   
init_dW <- obs_var / 100                 

# Step 2: Build the model
build_rw_bivariate <- function(param) {
  mod <- dlmModPoly(1)
  
  mod$FF <- mod$FF %x% diag(2) 
  mod$GG <- mod$GG %x% diag(2)
  
  V <- diag(c(param[1], param[2]))   
  W <- diag(c(param[3], param[4]))
  mod$V <- V         
  mod$W <- W         
  
  mod$m0 <- rep(ghcn_ts[1, ], 1)             
  mod$C0 <- diag(c(init_dV[1], init_dV[2])) 
  return(mod)
}

# Step 3: MLE Estimation
start_vals <- c(init_dV, init_dW)
lower_bounds <- rep(0.0001, 4)

MLE_biv_rw <- dlmMLE(ghcn_ts, parm = start_vals, build = build_rw_bivariate, lower = lower_bounds)
rw_bivariate_model <- build_rw_bivariate(MLE_biv_rw$par)


# Step 4: Filtering
filt_rw_bivariate <- dlmFilter(ghcn_ts,rw_bivariate_model) 

# Step 5: Smoothing 
smooth_rw_bivariate <- dlmSmooth(filt_rw_bivariate)


# Step 6: Forecasting 
# one-step forecast (for the last 100 observations)
f_rw_bivariate <- dropFirst(filt_rw_bivariate$f)
q_rw_bivariate <- filt_rw_bivariate$D.R[-1, ]

ghcn_ts1 <- ghcn_ts[, 1]
ghcn_ts2 <- ghcn_ts[, 2]

f1_rw_bivariate <- ts(f_rw_bivariate[, 1], start = start(ghcn_ts1), frequency = frequency(ghcn_ts1))
q1_rw_bivariate <- ts(q_rw_bivariate[, 1], start = start(ghcn_ts1), frequency = frequency(ghcn_ts1))

f2_rw_bivariate <- ts(f_rw_bivariate[, 2], start = start(ghcn_ts2), frequency = frequency(ghcn_ts2))
q2_rw_bivariate <- ts(q_rw_bivariate[, 2], start = start(ghcn_ts2), frequency = frequency(ghcn_ts2))

n_total <- length(f1_rw_bivariate)

start_index <- n_total - 99  

ghcn_ts1_last100 <- window(ghcn_ts1, start = time(ghcn_ts1)[start_index])
f1_rw_bivariate_last100 <- window(f1_rw_bivariate, start = time(f1_rw_bivariate)[start_index])
q1_rw_bivariate_last100 <- window(q1_rw_bivariate, start = time(q1_rw_bivariate)[start_index])

ghcn_ts2_last100 <- window(ghcn_ts2, start = time(ghcn_ts2)[start_index])
f2_rw_bivariate_last100 <- window(f2_rw_bivariate, start = time(f2_rw_bivariate)[start_index])
q2_rw_bivariate_last100 <- window(q2_rw_bivariate, start = time(q2_rw_bivariate)[start_index])

par(mfrow = c(2, 1))

ts.plot(ghcn_ts1_last100,
        f1_rw_bivariate_last100,
        f1_rw_bivariate_last100 + 1.96 * q1_rw_bivariate_last100,
        f1_rw_bivariate_last100 - 1.96 * q1_rw_bivariate_last100,
        lty = c(1, 1, 3, 3),
        col = c("lightgrey", "blue", "red", "red"),
        ylab = "Temperature",
        main = "Forecast for T min (Last 100 Observations)")

legend("topleft", legend = c("Observed", "Forecast", "95% CI"),
       lty = c(1, 1, 3), col = c("lightgrey", "blue", "red"), cex=0.5)

ts.plot(ghcn_ts2_last100,
        f2_rw_bivariate_last100,
        f2_rw_bivariate_last100 + 1.96 * q2_rw_bivariate_last100,
        f2_rw_bivariate_last100 - 1.96 * q2_rw_bivariate_last100,
        lty = c(1, 1, 3, 3),
        col = c("lightgrey", "blue", "red", "red"),
        ylab = "Temperature",
        main = "Forecast for T max (Last 100 Observations)")

legend("topleft", legend = c("Observed", "Forecast", "95% CI"),
       lty = c(1, 1, 3), col = c("lightgrey", "blue", "red"), cex=0.5)

par(mfrow = c(1, 1))


# Step 7: Checking 

# normality of residuals

res_rw_bivariate <- residuals(filt_rw_bivariate, sd = FALSE) 

#QQ-plot
qqnorm(res_rw_bivariate)
qqline(res_rw_bivariate)

# residuals distribution
hist(res_rw_bivariate) 


# Other diagnostic plots
tsdiag(filt_rw_bivariate)


#Shapiro test with a random subset
set.seed(123)  # for reproducibility
sample_res_rw_bivariate <- sample(res_rw_bivariate, size = 5000)

shapiro.test(sample_res_rw_bivariate)

# lack of serial correlation
# Box-Ljung test
Box.test(res_rw_bivariate[,1], lag = 20, type = "Ljung")

# Step 8: Forecast accuracy measures 

# MAE
mae_rw       <- mean(abs(filt_rw_bivariate$f - ghcn_ts))
mae_rw_min   <- mean(abs(filt_rw_bivariate$f[,1] - ghcn_ts[,1]))
mae_rw_max   <- mean(abs(filt_rw_bivariate$f[,2] - ghcn_ts[,2]))

# MSE
mse_rw       <- mean((filt_rw_bivariate$f - ghcn_ts)^2)
mse_rw_min   <- mean((filt_rw_bivariate$f[,1] - ghcn_ts[,1])^2)
mse_rw_max   <- mean((filt_rw_bivariate$f[,2] - ghcn_ts[,2])^2)

# MAPE
mape_rw      <- mean(abs(filt_rw_bivariate$f - ghcn_ts) / abs(ghcn_ts), na.rm = TRUE)
mape_rw_min  <- mean(abs(filt_rw_bivariate$f[,1] - ghcn_ts[,1]) / abs(ghcn_ts[,1]), na.rm = TRUE)
mape_rw_max  <- mean(abs(filt_rw_bivariate$f[,2] - ghcn_ts[,2]) / abs(ghcn_ts[,2]), na.rm = TRUE)

# Table
metrics_table_rw <- data.frame(
  Metric  = c("MAE", "MSE", "MAPE"),
  Overall = c(mae_rw, mse_rw, mape_rw),
  Min     = c(mae_rw_min, mse_rw_min, mape_rw_min),
  Max     = c(mae_rw_max, mse_rw_max, mape_rw_max)
)

print(metrics_table_rw)

####################################################################

# Hold-out validation

# Step 0: Split the data
n <- nrow(ghcn_ts)
train_index <- floor(0.8 * n)

train_data <- ghcn_ts[1:train_index, ]
train_data <- ts(train_data, start = c(as.numeric(format(min(ghcn$date), "%Y")), 
                                       as.numeric(format(min(ghcn$date), "%j"))), 
                 frequency = 365)
test_data <- ghcn_ts[(train_index + 1):n, ]
test_data <- ts(test_data, start = c(as.numeric(format(min(ghcn$date), "%Y")), 
                                     as.numeric(format(min(ghcn$date), "%j"))), 
                frequency = 365)
# if I do plot with test data, I need to adjust years on X-axis

n_forecast <- nrow(test_data)

# Step 1: Fit model on training set
MLE_holdout <- dlmMLE(train_data, parm = start_vals, build = build_rw_bivariate, lower = lower_bounds)
model_holdout <- build_rw_bivariate(MLE_holdout$par)
filt_holdout <- dlmFilter(test_data, model_holdout)

# Step 2: Forecast the test set 
# one-step forecast (comes from filtering)
forecast_holdout <- filt_holdout$f
# One-step-ahead forecast standard deviations from SVD of forecast error variances
q_rw <- filt_rw_bivariate$D.R[-1, ]


# Step 4: Accuracy metrics for hold-out forecast
mae_cv <- mean(abs(forecast_holdout - test_data))
mse_cv <- mean((forecast_holdout - test_data)^2)
mape_cv <- mean(abs(forecast_holdout - test_data) / test_data)

# Step 5: Accuracy metrics for in-sample (full model)
mae_rw <- mean(abs(filt_rw_bivariate$f - ghcn_ts))
mse_rw <- mean((filt_rw_bivariate$f - ghcn_ts)^2)
mape_rw <- mean(abs(filt_rw_bivariate$f - ghcn_ts) / ghcn_ts)

# Step 6: Summary table
comparison <- data.frame(
  Model = c("Without HOV", "HOV"),
  MAE = c(mae_rw, mae_cv),
  MSE = c(mse_rw, mse_cv),
  MAPE = c(mape_rw, mape_cv)
)

print(comparison)



##################################     
## Step 2. Model B 
##################################

# Step 1: Initial parameter guesses
obs_var <- apply(ghcn_ts, 2, var)
init_dV <- obs_var / 10
init_dW <- obs_var / 100

# Start values: V (2) + W (3: var1, var2, cov)
start_vals <- c(init_dV[1], init_dV[2], init_dW[1], init_dW[2], 0)  # last 0 is initial guess for covariance
lower_bounds <- c(0.0001, 0.0001, 0.0001, 0.0001, -Inf)

# Step 2: Build model with full W matrix
build_surw_bivariate <- function(param) {
  mod <- dlmModPoly(1)
  
  mod$FF <- mod$FF %x% diag(2)
  mod$GG <- mod$GG %x% diag(2)
  
  # V remains diagonal
  V <- diag(c(param[1], param[2]))
  
  # W is full symmetric matrix
  w11 <- param[3]
  w22 <- param[4]
  w12 <- param[5]
  W <- matrix(c(w11, w12, w12, w22), nrow = 2)
  
  mod$V <- V
  mod$W <- W
  mod$m0 <- rep(ghcn_ts[1, ], 1)
  mod$C0 <- diag(c(init_dV[1], init_dV[2]))
  
  return(mod)
}

# Step 3: MLE Estimation
MLE_surw_bivariate <- dlmMLE(ghcn_ts, parm = start_vals, build = build_surw_bivariate, lower = lower_bounds)
surw_bivariate_model <- build_surw_bivariate(MLE_surw_bivariate$par)

# Step 4: Filtering
filt_surw_bivariate <- dlmFilter(ghcn_ts,surw_bivariate_model) 

# Step 5: Smoothing 
smooth_surw_bivariate <- dlmSmooth(filt_surw_bivariate)


# Step 6: Forecasting 

f_surw_bivariate <- dropFirst(filt_surw_bivariate$f)
q_surw_bivariate <- filt_surw_bivariate$D.R[-1, ]

ghcn_ts1 <- ghcn_ts[, 1]
ghcn_ts2 <- ghcn_ts[, 2]

f1_surw_bivariate <- ts(f_surw_bivariate[, 1], start = start(ghcn_ts1), frequency = frequency(ghcn_ts1))
q1_surw_bivariate <- ts(q_surw_bivariate[, 1], start = start(ghcn_ts1), frequency = frequency(ghcn_ts1))

f2_surw_bivariate <- ts(f_surw_bivariate[, 2], start = start(ghcn_ts2), frequency = frequency(ghcn_ts2))
q2_surw_bivariate <- ts(q_surw_bivariate[, 2], start = start(ghcn_ts2), frequency = frequency(ghcn_ts2))

n_total <- length(f1_surw_bivariate)

start_index <- n_total - 99  

ghcn_ts1_last100 <- window(ghcn_ts1, start = time(ghcn_ts1)[start_index])
f1_surw_bivariate_last100 <- window(f1_surw_bivariate, start = time(f1_surw_bivariate)[start_index])
q1_surw_bivariate_last100 <- window(q1_surw_bivariate, start = time(q1_surw_bivariate)[start_index])

ghcn_ts2_last100 <- window(ghcn_ts2, start = time(ghcn_ts2)[start_index])
f2_surw_bivariate_last100 <- window(f2_surw_bivariate, start = time(f2_surw_bivariate)[start_index])
q2_surw_bivariate_last100 <- window(q2_surw_bivariate, start = time(q2_surw_bivariate)[start_index])

par(mfrow = c(2, 1))

ts.plot(ghcn_ts1_last100,
        f1_surw_bivariate_last100,
        f1_surw_bivariate_last100 + 1.96 * q1_surw_bivariate_last100,
        f1_surw_bivariate_last100 - 1.96 * q1_surw_bivariate_last100,
        lty = c(1, 1, 3, 3),
        col = c("lightgrey", "blue", "red", "red"),
        ylab = "Temperature",
        main = "Forecast for T min (Last 100 Observations)")

legend("topleft", legend = c("Observed", "Forecast", "95% CI"),
       lty = c(1, 1, 3), col = c("lightgrey", "blue", "red"), cex = 0.5)

ts.plot(ghcn_ts2_last100,
        f2_surw_bivariate_last100,
        f2_surw_bivariate_last100 + 1.96 * q2_surw_bivariate_last100,
        f2_surw_bivariate_last100 - 1.96 * q2_surw_bivariate_last100,
        lty = c(1, 1, 3, 3),
        col = c("lightgrey", "blue", "red", "red"),
        ylab = "Temperature",
        main = "Forecast for T max (Last 100 Observations)")

legend("topleft", legend = c("Observed", "Forecast", "95% CI"),
       lty = c(1, 1, 3), col = c("lightgrey", "blue", "red"), cex = 0.5)

par(mfrow = c(1, 1))



# Step 7: Checking 

# normality of residuals

res_surw_bivariate <- residuals(filt_surw_bivariate, sd = FALSE) 

#QQ-plot
qqnorm(res_surw_bivariate)
qqline(res_surw_bivariate)

# residuals distribution
hist(res_surw_bivariate) 


# Other diagnostic plots
tsdiag(filt_surw_bivariate)


#Shapiro test with a random subset
set.seed(123)  # for reproducibility
sample_res_surw_bivariate <- sample(res_surw_bivariate, size = 5000)

shapiro.test(sample_res_surw_bivariate)

# lack of serial correlation
# Box-Ljung test
Box.test(res_surw_bivariate[,1], lag = 20, type = "Ljung")


# Step 8: Forecast accuracy measures 

# MAE
mae_surw       <- mean(abs(filt_surw_bivariate$f - ghcn_ts))
mae_surw_min   <- mean(abs(filt_surw_bivariate$f[,1] - ghcn_ts[,1]))
mae_surw_max   <- mean(abs(filt_surw_bivariate$f[,2] - ghcn_ts[,2]))

# MSE
mse_surw       <- mean((filt_surw_bivariate$f - ghcn_ts)^2)
mse_surw_min   <- mean((filt_surw_bivariate$f[,1] - ghcn_ts[,1])^2)
mse_surw_max   <- mean((filt_surw_bivariate$f[,2] - ghcn_ts[,2])^2)

# MAPE
mape_surw      <- mean(abs(filt_surw_bivariate$f - ghcn_ts) / abs(ghcn_ts), na.rm = TRUE)
mape_surw_min  <- mean(abs(filt_surw_bivariate$f[,1] - ghcn_ts[,1]) / abs(ghcn_ts[,1]), na.rm = TRUE)
mape_surw_max  <- mean(abs(filt_surw_bivariate$f[,2] - ghcn_ts[,2]) / abs(ghcn_ts[,2]), na.rm = TRUE)

# Table
metrics_table_surw <- data.frame(
  Metric  = c("MAE", "MSE", "MAPE"),
  Overall = c(mae_surw, mse_surw, mape_surw),
  Min     = c(mae_surw_min, mse_surw_min, mape_surw_min),
  Max     = c(mae_surw_max, mse_surw_max, mape_surw_max)
)

print(metrics_table_surw)



##################################     
## Step 2. Model C 
##################################


# Step 1: Empirical variances as starting point
obs_var <- apply(ghcn_ts, 2, var)
init_dV <- obs_var / 10
init_dW <- mean(obs_var) / 100

# Step 2: Model builder
build_latent_rw <- function(param) {
  alpha1 <- param[1]
  alpha2 <- param[2]
  beta   <- param[3]
  sigmaV1 <- param[4]
  sigmaV2 <- param[5]
  sigmaW  <- param[6]
  
  # Observation matrix: F %*% theta_t
  FF <- matrix(c(alpha1, beta,
                 alpha2, 1 / beta), nrow = 2, byrow = TRUE)
  
  # Transition matrix: only latent state evolves
  GG <- diag(2)
  
  # State noise: only second component (latent factor) has noise
  W <- matrix(0, 2, 2)
  W[2, 2] <- sigmaW^2
  
  V <- diag(c(sigmaV1^2, sigmaV2^2))  # observation noise
  
  # DLM object
  mod <- dlm(
    FF = FF,
    GG = GG,
    V = V,
    W = W,
    m0 = c(1, 0),                      # initial state: fixed intercept, latent = 0
    C0 = diag(c(0, init_dW))       # no uncertainty on fixed 1, more on latent
  )
  
  return(mod)
}

# Step 3: MLE Estimation
start_vals <- c(0, 0, 1, sqrt(init_dV[1]), sqrt(init_dV[2]), sqrt(init_dW))  # alpha1, alpha2, beta, sigmaV1, sigmaV2, sigmaW
lower_bounds <- rep(0.0001, length(start_vals))  # all parameters must be > 0

MLE_latent_rw <- dlmMLE(ghcn_ts, parm = start_vals, build = build_latent_rw, lower = lower_bounds)
latent_model <- build_latent_rw(MLE_latent_rw$par)

# Extract raw parameter vector
estimated_params <- MLE_latent_rw$par
names(estimated_params) <- c("alpha1", "alpha2", "beta", "sigmaV1", "sigmaV2", "sigmaW")
print(estimated_params)

# Step 4: Filtering 
filt_latent_rw_bivariate <- dlmFilter(ghcn_ts, latent_model)

# Step 5: Smoothing
smooth_latent_rw_bivariate <- dlmSmooth(filt_latent_rw_bivariate)

latent_component <- smooth_latent_rw_bivariate$s[-1, 2]
latent_component_ts <- ts(latent_component, start = c(1921, 1), frequency = 365)

plot(latent_component_ts, type = "l",
     col = 'grey', lty = 1,
     ylab = expression(xi[t]), xlab = "Year",
     main = "Estimated Latent Factor (ξ)")

# Fit linear trend
time_index <- time(latent_component_ts)
trend_model <- lm(latent_component ~ time_index)

abline(trend_model, col = "blue", lty = 2, lwd = 2)
abline(h = 0, col = "black", lty = 1)

legend("topleft", legend = c("Estimates", "Trend"),
       lty = c(1, 2), col = c("grey", "blue"), cex = 0.8)

summary(trend_model)

slope <- coef(trend_model)[2]
pval <- summary(trend_model)$coefficients[2, 4]
cat("Estimated trend slope:", slope, "\n")
cat("P-value for slope:", pval, "\n")
if (pval < 0.05 && slope > 0) {
  cat("→ The latent component shows a statistically significant upward trend.\n")
} else {
  cat("→ No significant upward trend detected.\n")
}

# Step 6: Forecasting 

f_latent_rw_bivariate <- dropFirst(filt_latent_rw_bivariate$f)
q_latent_rw_bivariate <- filt_latent_rw_bivariate$D.R[-1, ]

ghcn_ts1 <- ghcn_ts[, 1]
ghcn_ts2 <- ghcn_ts[, 2]

f1_latent_rw_bivariate <- ts(f_latent_rw_bivariate[, 1], start = start(ghcn_ts1), frequency = frequency(ghcn_ts1))
q1_latent_rw_bivariate <- ts(q_latent_rw_bivariate[, 1], start = start(ghcn_ts1), frequency = frequency(ghcn_ts1))

f2_latent_rw_bivariate <- ts(f_latent_rw_bivariate[, 2], start = start(ghcn_ts2), frequency = frequency(ghcn_ts2))
q2_latent_rw_bivariate <- ts(q_latent_rw_bivariate[, 2], start = start(ghcn_ts2), frequency = frequency(ghcn_ts2))

n_total <- length(f1_latent_rw_bivariate)

start_index <- n_total - 99

ghcn_ts1_last100 <- window(ghcn_ts1, start = time(ghcn_ts1)[start_index])
f1_latent_rw_bivariate_last100 <- window(f1_latent_rw_bivariate, start = time(f1_latent_rw_bivariate)[start_index])
q1_latent_rw_bivariate_last100 <- window(q1_latent_rw_bivariate, start = time(q1_latent_rw_bivariate)[start_index])

ghcn_ts2_last100 <- window(ghcn_ts2, start = time(ghcn_ts2)[start_index])
f2_latent_rw_bivariate_last100 <- window(f2_latent_rw_bivariate, start = time(f2_latent_rw_bivariate)[start_index])
q2_latent_rw_bivariate_last100 <- window(q2_latent_rw_bivariate, start = time(q2_latent_rw_bivariate)[start_index])

par(mfrow = c(2, 1))

ts.plot(ghcn_ts1_last100,
        f1_latent_rw_bivariate_last100,
        f1_latent_rw_bivariate_last100 + 1.96 * q1_latent_rw_bivariate_last100,
        f1_latent_rw_bivariate_last100 - 1.96 * q1_latent_rw_bivariate_last100,
        lty = c(1, 1, 3, 3),
        col = c("lightgrey", "blue", "red", "red"),
        ylab = "Temperature",
        main = "Forecast for T min (Last 100 Observations)")

legend("topleft", legend = c("Observed", "Forecast", "95% CI"),
       lty = c(1, 2, 3), col = c("lightgrey", "blue", "red"), cex = 0.5)

ts.plot(ghcn_ts2_last100,
        f2_latent_rw_bivariate_last100,
        f2_latent_rw_bivariate_last100 + 1.96 * q2_latent_rw_bivariate_last100,
        f2_latent_rw_bivariate_last100 - 1.96 * q2_latent_rw_bivariate_last100,
        lty = c(1, 1, 3, 3),
        col = c("lightgrey", "blue", "red", "red"),
        ylab = "Temperature",
        main = "Forecast for T max (Last 100 Observations)")

legend("topleft", legend = c("Observed", "Forecast", "95% CI"),
       lty = c(1, 2, 3), col = c("lightgrey", "blue", "red"), cex = 0.5)

par(mfrow = c(1, 1))


# Step 7: Checking 

# normality of residuals

res_latent_rw_bivariate <- residuals(filt_latent_rw_bivariate, sd = FALSE) 

#QQ-plot
qqnorm(res_latent_rw_bivariate)
qqline(res_latent_rw_bivariate)

# residuals distribution
hist(res_latent_rw_bivariate) 

# Other diagnostic plots
tsdiag(filt_latent_rw_bivariate)

#Shapiro test with a random subset
set.seed(123)  # for reproducibility
sample_res_latent_rw_bivariate <- sample(res_latent_rw_bivariate, size = 5000)

shapiro.test(sample_res_latent_rw_bivariate)

# lack of serial correlation
# Box-Ljung test
Box.test(res_latent_rw_bivariate[,1], lag = 20, type = "Ljung")


# Step 8: Forecast accuracy measures 

# MAE
mae_latent_rw       <- mean(abs(filt_latent_rw_bivariate$f - ghcn_ts))
mae_latent_rw_min   <- mean(abs(filt_latent_rw_bivariate$f[,1] - ghcn_ts[,1]))
mae_latent_rw_max   <- mean(abs(filt_latent_rw_bivariate$f[,2] - ghcn_ts[,2]))

# MSE
mse_latent_rw       <- mean((filt_latent_rw_bivariate$f - ghcn_ts)^2)
mse_latent_rw_min   <- mean((filt_latent_rw_bivariate$f[,1] - ghcn_ts[,1])^2)
mse_latent_rw_max   <- mean((filt_latent_rw_bivariate$f[,2] - ghcn_ts[,2])^2)

# MAPE
mape_latent_rw      <- mean(abs(filt_latent_rw_bivariate$f - ghcn_ts) / abs(ghcn_ts), na.rm = TRUE)
mape_latent_rw_min  <- mean(abs(filt_latent_rw_bivariate$f[,1] - ghcn_ts[,1]) / abs(ghcn_ts[,1]), na.rm = TRUE)
mape_latent_rw_max  <- mean(abs(filt_latent_rw_bivariate$f[,2] - ghcn_ts[,2]) / abs(ghcn_ts[,2]), na.rm = TRUE)

# Table
metrics_table_latent_rw <- data.frame(
  Metric  = c("MAE", "MSE", "MAPE"),
  Overall = c(mae_latent_rw, mse_latent_rw, mape_latent_rw),
  Min     = c(mae_latent_rw_min, mse_latent_rw_min, mape_latent_rw_min),
  Max     = c(mae_latent_rw_max, mse_latent_rw_max, mape_latent_rw_max)
)

print(metrics_table_latent_rw)
