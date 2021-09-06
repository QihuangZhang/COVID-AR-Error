### This version is forked from Version 1 with following update:

# 1. Change the study endpoint into fatality rate 
# 2. Use the information from He et. al for the asymptomatic cases


### Revision to applied statistics

# 1. Preperation -------------------------------------------------------------

# 1.1 Packages ------------------------------------------------------------

library(tidyr)
require(dplyr)
library(aTSA)
library(forecast)
library(data.table)

library(ggplot2)
require(ggthemes)

library(xtable)

# 2. Data --------------------------------------------------------------------
source(file="code/COVID_functions.R")

# import the cleaned data

source(file="code/COVID_data_Prep.R")

FataRate_Case3_ON_comp <- ON_death$cases / ON_CumCase$Case * 100

FataRate_Case3_BC_comp <- BC_death$cases / BC_CumCase$Case * 100

FataRate_Case3_QC_comp <- QC_death$cases / QC_CumCase$Case * 100

FataRate_Case3_AB_comp <- AB_death$cases / AB_CumCase$Case * 100

FataRate_Case3_ON <- FataRate_Case3_ON_comp[1:(length(FataRate_Case3_ON_comp)-5)]

FataRate_Case3_BC <- FataRate_Case3_BC_comp[1:(length(FataRate_Case3_BC_comp)-5)]

FataRate_Case3_QC <- FataRate_Case3_QC_comp[1:(length(FataRate_Case3_QC_comp)-5)]

FataRate_Case3_AB <- FataRate_Case3_AB_comp[1:(length(FataRate_Case3_AB_comp)-5)]


# 3 Canada Data - Fatality Rate Based Study --------------------------------------------------------------------

# 3.1 British Columbia Data --------------------------------------------------------------------

Date_BC <- BC_Confirm$count

tseries::adf.test(FataRate_Case3_BC)
tseries::adf.test(diff(FataRate_Case3_BC))

produce_pacf_plot(FataRate_Case3_BC, "British Columbia - no differencing", defnumber = 3)
produce_pacf_plot(diff(FataRate_Case3_BC), "British Columbia", defnumber = 3)
# 
# 
find_optimal_lag(diff(FataRate_Case3_BC))
# 
# 
# ### 3.1.2 BC - AR(5) - diff #### 
# 
# Fitted_Models <- diff_fit (tseries = FataRate_Case3_BC, lag = 5, sigma_e = 0.01, sigma_u = 0.002) 
# 
# # xtable(Fitted_Models$Table1)
# 
# Prediction_prov <- prediction_diff(FataRate_Case3_BC, Models = Fitted_Models, Province = "British Columbia", fullseries = FataRate_Case3_BC_comp, origindate = BC_CumCase$Date[1])
# 
# # plot_results(Prediction_prov, defnumber = 3, yrange = c(0,2.5), diff = T)
# 
# Pe_h_est_BC_AR5 <- Evaluations_data(Prediction_prov)
# 
# 
# ### 3.1.2 BC - AR(6) - diff #### 
# 
# Fitted_Models <- diff_fit (tseries = FataRate_Case3_BC, lag = 6, sigma_e = 0.01, sigma_u = 0.002) 
# 
# xtable(Fitted_Models$Table1, digits = 3)
# 
# Prediction_prov <- prediction_diff(FataRate_Case3_BC, Models = Fitted_Models, Province = "British Columbia", fullseries = FataRate_Case3_BC_comp, origindate = BC_CumCase$Date[1])
# 
# plot_results(Prediction_prov, defnumber = 3, yrange = c(0,2.5), diff = T)
# 
# Pe_h_est_BC_AR6 <- Evaluations_data(Prediction_prov)
# 
# 
# 
# ####  (--------------------------------------------------------------)  ####
# # 3.2 Ontario Data --------------------------------------------------------
# 
Date_ON <- ON_Confirm$count

tseries::adf.test(FataRate_Case3_ON)

tseries::adf.test(diff(FataRate_Case3_ON))

find_optimal_lag(FataRate_Case3_ON)

produce_pacf_plot(FataRate_Case3_ON, "Ontario - no differencing", defnumber = 3)
produce_pacf_plot(diff(FataRate_Case3_ON), "Ontario", defnumber = 3)
# 
# ### 3.2.1 ON-AR(2)  ####
# 
# Fitted_Models <- nodiff_fit (tseries = FataRate_Case3_ON, lag = 2, sigma_e = 0.02, sigma_u = 0.001) 
# 
# # xtable(Fitted_Models$Table1)
# 
# Prediction_prov <- prediction_nodiff(FataRate_Case3_ON, Models = Fitted_Models, Province = "Ontario", fullseries = FataRate_Case3_ON_comp, origindate = ON_CumCase$Date[1])
# 
# # plot_results(Prediction_prov, defnumber = 3, yrange = c(0,2.5))
# 
# Pe_h_est_ON_AR2 <- Evaluations_data(Prediction_prov)
# 
# 
# ### 3.2.2 ON-AR(3)  (chosen) ####
# 
# Fitted_Models <- nodiff_fit (tseries = FataRate_Case3_ON, lag = 3, sigma_e = 0.02, sigma_u = 0.001) 
# 
# xtable(Fitted_Models$Table1, digits = 3)
# 
# Prediction_prov <- prediction_nodiff(FataRate_Case3_ON, Models = Fitted_Models, Province = "Ontario", fullseries = FataRate_Case3_ON_comp, origindate = ON_CumCase$Date[1])
# 
# plot_results(Prediction_prov, defnumber = 3, yrange = c(0,2.5))
# 
# Pe_h_est_ON_AR3 <- Evaluations_data(Prediction_prov)
# 
# 
# 
# 
# # !-----------------------------------------             ------------------
# # 3.3 Quebec Data --------------------------------------------------------
# 
Date_QC <- QC_Confirm$count

tseries::adf.test(FataRate_Case3_QC)
tseries::adf.test(diff(FataRate_Case3_QC))

produce_pacf_plot(FataRate_Case3_QC, "Quebec - no differencing", defnumber = 3)
produce_pacf_plot(diff(FataRate_Case3_QC), "Quebec", defnumber = 3)
# 
# find_optimal_lag(diff(FataRate_Case3_QC))
# 
# ### 3.3.1 QC-AR(3) - diff ####
# 
# Fitted_Models <- diff_fit (tseries = FataRate_Case3_QC, lag = 3, sigma_e = 0.002, sigma_u = 0.001) 
# 
# # xtable(Fitted_Models$Table1)
# 
# Prediction_prov <- prediction_diff(FataRate_Case3_QC, Models = Fitted_Models, Province = "Quebec", fullseries = FataRate_Case3_QC_comp, origindate = QC_CumCase$Date[1])
# 
# # plot_results(Prediction_prov, defnumber = 3, yrange = c(0.5,7), diff = T)
# 
# Pe_h_est_QC_AR3 <- Evaluations_data(Prediction_prov)
# 
# ### 3.3.2 QC-AR(4) (chosen)  - diff ####
# 
# Fitted_Models <- diff_fit (tseries = FataRate_Case3_QC, lag = 4, sigma_e = 0.002, sigma_u = 0.001) 
# 
# xtable(Fitted_Models$Table1, digits = 3)
# 
# Prediction_prov <- prediction_diff(FataRate_Case3_QC, Models = Fitted_Models, Province = "Quebec", fullseries = FataRate_Case3_QC_comp, origindate = QC_CumCase$Date[1])
# 
# plot_results(Prediction_prov, defnumber = 3, yrange = c(0.5,7), diff = T)
# 
# Pe_h_est_QC_AR4 <- Evaluations_data(Prediction_prov)
# 
# 
# # !-----------------------------------------             ------------------
# # 3.3 Alberta Data --------------------------------------------------------
Date_AB <- AB_Confirm$count

tseries::adf.test(FataRate_Case3_AB)
tseries::adf.test(diff(FataRate_Case3_AB))

produce_pacf_plot(FataRate_Case3_AB, "Alberta - no differencing", defnumber = 3)
produce_pacf_plot(diff(FataRate_Case3_AB), "Alberta", defnumber = 3)

find_optimal_lag(diff(FataRate_Case3_AB))


### 3.2.1 AB-AR(3)  ####

Fitted_Models <- nodiff_fit (tseries = FataRate_Case3_AB, lag = 3, sigma_e = 0.002, sigma_u = 0.001)

# xtable(Fitted_Models$Table1)

Prediction_prov <- prediction_nodiff(FataRate_Case3_AB, Models = Fitted_Models, Province = "Alberta", fullseries = FataRate_Case3_AB_comp, origindate = AB_CumCase$Date[1])

# plot_results(Prediction_prov, defnumber = 3, yrange = c(0,2))

Pe_h_est_AB_AR3 <- Evaluations_data(Prediction_prov)

### 3.2.2 AB-AR(4)  (chosen)  ####

Fitted_Models <- nodiff_fit (tseries = FataRate_Case3_AB, lag = 4, sigma_e = 0.002, sigma_u = 0.001)

xtable(Fitted_Models$Table1, digits = 3)

Prediction_prov <- prediction_nodiff(FataRate_Case3_AB, Models = Fitted_Models, Province = "Alberta", fullseries = FataRate_Case3_AB_comp, origindate = AB_CumCase$Date[1])

plot_results(Prediction_prov, defnumber = 3, yrange = c(0.3,1.6), diff = T)

Pe_h_est_AB_AR4 <- Evaluations_data(Prediction_prov)


# 6 Summarize Results --------------------------------------------------------------------
## Table of British Columbia

# Table_BC <- get_Peh_table_prov(Pe_h_est_list = list(Pe_h_est_BC_AR5, Pe_h_est_BC_AR6))
# Table_ON <- get_Peh_table_prov(Pe_h_est_list = list(Pe_h_est_ON_AR2, Pe_h_est_ON_AR3))
# Table_QC <- get_Peh_table_prov(Pe_h_est_list = list(Pe_h_est_QC_AR3, Pe_h_est_QC_AR4))
Table_AB <- get_Peh_table_prov(Pe_h_est_list = list(Pe_h_est_AB_AR3, Pe_h_est_AB_AR4))

Tables_Peh <- rbind(Table_AB)


require(xtable)
xtable(Tables_Peh, digits = 3)
