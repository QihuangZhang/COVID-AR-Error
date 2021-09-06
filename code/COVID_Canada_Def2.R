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

FataRate_Case2_ON_comp <- ON_death$cases / c(ON_CumCase$Case[1:10], ON_CumCase$Case[-((length(ON_CumCase$Case)-9):length(ON_CumCase$Case))] ) * 100

FataRate_Case2_BC_comp <- BC_death$cases / c(BC_CumCase$Case[1:10], BC_CumCase$Case[-((length(BC_CumCase$Case)-9):length(BC_CumCase$Case))] ) * 100

FataRate_Case2_QC_comp <- QC_death$cases / c(QC_CumCase$Case[1:10], QC_CumCase$Case[-((length(QC_CumCase$Case)-9):length(QC_CumCase$Case))] ) * 100

FataRate_Case2_AB_comp <- AB_death$cases / c(AB_CumCase$Case[1:10], AB_CumCase$Case[-((length(AB_CumCase$Case)-9):length(AB_CumCase$Case))] ) * 100

FataRate_Case2_ON <- FataRate_Case2_ON_comp[1:(length(FataRate_Case2_ON_comp)-5)]

FataRate_Case2_BC <- FataRate_Case2_BC_comp[1:(length(FataRate_Case2_BC_comp)-5)]

FataRate_Case2_QC <- FataRate_Case2_QC_comp[1:(length(FataRate_Case2_QC_comp)-5)]

FataRate_Case2_AB <- FataRate_Case2_AB_comp[1:(length(FataRate_Case2_AB_comp)-5)]


# 3 Canada Data - Fatality Rate Based Study --------------------------------------------------------------------

# 3.1 British Columbia Data --------------------------------------------------------------------

Date_BC <- BC_Confirm$count

tseries::adf.test(FataRate_Case2_BC)

tseries::adf.test(diff(FataRate_Case2_BC))

produce_pacf_plot(FataRate_Case2_BC, "British Columbia - no differencing", defnumber = 2)

produce_pacf_plot(diff(FataRate_Case2_BC), "British Columbia", defnumber = 2)


find_optimal_lag(diff(FataRate_Case2_BC))


### 3.1.2 BC - AR(1)  #### 

Fitted_Models <- diff_fit (tseries = FataRate_Case2_BC, lag = 1, sigma_e = 0.01, sigma_u = 0.005) 

# xtable(Fitted_Models$Table1)

Prediction_prov <- prediction_diff(FataRate_Case2_BC, Models = Fitted_Models, Province = "British Columbia", fullseries = FataRate_Case2_BC_comp, origindate = BC_CumCase$Date[1])

# plot_results(Prediction_prov, defnumber = 2, diff = T, yrange = c(0,2.5), diff = T)

Pe_h_est_BC_AR1 <- Evaluations_data(Prediction_prov)


### 3.1.2 BC - AR(2) (chosen) #### 

Fitted_Models <- diff_fit (tseries = FataRate_Case2_BC, lag = 2, sigma_e = 0.01, sigma_u = 0.005) 

xtable(Fitted_Models$Table1, digits = 3)

Prediction_prov <- prediction_diff(FataRate_Case2_BC, Models = Fitted_Models, Province = "British Columbia", fullseries = FataRate_Case2_BC_comp, origindate = BC_CumCase$Date[1])

plot_results(Prediction_prov, defnumber = 2, yrange = c(0.5,2.1), diff = T)

Pe_h_est_BC_AR2 <- Evaluations_data(Prediction_prov)





####  (--------------------------------------------------------------)  ####
# 3.2 Ontario Data --------------------------------------------------------

Date_ON <- ON_Confirm$count

tseries::adf.test(FataRate_Case2_ON)
tseries::adf.test(diff(FataRate_Case2_ON))

find_optimal_lag(diff(FataRate_Case2_ON))

produce_pacf_plot(FataRate_Case2_ON, "Ontario - no differencing", defnumber = 2)

produce_pacf_plot(diff(FataRate_Case2_ON), "Ontario", defnumber = 2)

### 3.2.1 ON-AR(1) (chosen) ####

Fitted_Models <- diff_fit (tseries = FataRate_Case2_ON, lag = 1, sigma_e = 0.02, sigma_u = 0.002)

xtable(Fitted_Models$Table1, digits = 3)

Prediction_prov <- prediction_diff(FataRate_Case2_ON, Models = Fitted_Models, Province = "Ontario", fullseries = FataRate_Case2_ON_comp, origindate = ON_CumCase$Date[1])

plot_results(Prediction_prov, defnumber = 2, diff = T, yrange = c(1,4))

Pe_h_est_ON_AR1 <- Evaluations_data(Prediction_prov)


### 3.2.2 ON-AR(2)   ####

Fitted_Models <- diff_fit (tseries = FataRate_Case2_ON, lag = 2, sigma_e = 0.02, sigma_u = 0.002) 

# xtable(Fitted_Models$Table1, digits = 3)

Prediction_prov <- prediction_diff(FataRate_Case2_ON, Models = Fitted_Models, Province = "Ontario", fullseries = FataRate_Case2_ON_comp, origindate = ON_CumCase$Date[1])

# plot_results(Prediction_prov, defnumber = 2, diff = T, yrange = c(0,2.5))

Pe_h_est_ON_AR2 <- Evaluations_data(Prediction_prov)



# !-----------------------------------------             ------------------
# 3.3 Quebec Data --------------------------------------------------------

Date_QC <- QC_Confirm$count

tseries::adf.test(FataRate_Case2_QC)
tseries::adf.test(diff(FataRate_Case2_QC))

produce_pacf_plot(FataRate_Case2_QC, "Quebec - no differencing", defnumber = 2)

produce_pacf_plot(diff(FataRate_Case2_QC), "Quebec", defnumber = 2)

find_optimal_lag(diff(FataRate_Case2_QC))

### 3.3.1 QC-AR(7) ####

Fitted_Models <- diff_fit (tseries = FataRate_Case2_QC, lag = 7, sigma_e = 0.004, sigma_u = 0.01) 

# xtable(Fitted_Models$Table1, digits = 3)

Prediction_prov <- prediction_diff(FataRate_Case2_QC, Models = Fitted_Models, Province = "Quebec", fullseries = FataRate_Case2_QC_comp, origindate = QC_CumCase$Date[1])

# plot_results(Prediction_prov, defnumber = 2, diff = T, yrange = c(0.5,7))

Pe_h_est_QC_AR7 <- Evaluations_data(Prediction_prov)

### 3.3.2 QC-AR(8) (chosen)####

Fitted_Models <- diff_fit (tseries = FataRate_Case2_QC, lag = 8, sigma_e = 0.001, sigma_u = 0.01) 

xtable(Fitted_Models$Table1, digits = 3)

Prediction_prov <- prediction_diff(FataRate_Case2_QC, Models = Fitted_Models, Province = "Quebec", fullseries = FataRate_Case2_QC_comp, origindate = QC_CumCase$Date[1])

plot_results(Prediction_prov, defnumber = 2, diff = T, yrange = c(1.5,5))

Pe_h_est_QC_AR8 <- Evaluations_data(Prediction_prov)




# !-----------------------------------------             ------------------
# 3.3 Alberta Data --------------------------------------------------------
Date_AB <- AB_Confirm$count

tseries::adf.test(FataRate_Case2_AB)
tseries::adf.test(diff(FataRate_Case2_AB))

produce_pacf_plot(FataRate_Case2_AB, "Alberta - no differencing", defnumber = 2)

produce_pacf_plot(diff(FataRate_Case2_AB), "Alberta", defnumber = 2)

find_optimal_lag(diff(FataRate_Case2_AB))


### 3.2.1 AB-AR(1)  ####

Fitted_Models <- diff_fit (tseries = FataRate_Case2_AB, lag = 1, sigma_e = 0.005, sigma_u = 0.005) 

# xtable(Fitted_Models$Table1)

Prediction_prov <- prediction_diff(FataRate_Case2_AB, Models = Fitted_Models, Province = "Alberta", fullseries = FataRate_Case2_AB_comp, origindate = AB_CumCase$Date[1])

# plot_results(Prediction_prov, defnumber = 2, diff = T, yrange = c(0,2))

Pe_h_est_AB_AR1 <- Evaluations_data(Prediction_prov)

### 3.2.2 AB-AR(2)  (chosen)  ####

Fitted_Models <- diff_fit (tseries = FataRate_Case2_AB, lag = 2, sigma_e = 0.005, sigma_u = 0.005) 

xtable(Fitted_Models$Table1, digits = 3)

Prediction_prov <- prediction_diff(FataRate_Case2_AB, Models = Fitted_Models, Province = "Alberta", fullseries = FataRate_Case2_AB_comp, origindate = AB_CumCase$Date[1])

plot_results(Prediction_prov, defnumber = 2, diff = T, yrange = c(0.3,1.6))

Pe_h_est_AB_AR2 <- Evaluations_data(Prediction_prov)


# 6 Summarize Results --------------------------------------------------------------------
## Table of British Columbia

Table_BC <- get_Peh_table_prov(Pe_h_est_list = list(Pe_h_est_BC_AR1, Pe_h_est_BC_AR2))
Table_ON <- get_Peh_table_prov(Pe_h_est_list = list(Pe_h_est_ON_AR1, Pe_h_est_ON_AR2))
Table_QC <- get_Peh_table_prov(Pe_h_est_list = list(Pe_h_est_QC_AR7, Pe_h_est_QC_AR8))
Table_AB <- get_Peh_table_prov(Pe_h_est_list = list(Pe_h_est_AB_AR1, Pe_h_est_AB_AR2))

Tables_Peh <- rbind(Table_BC, Table_ON, Table_QC, Table_AB)


require(xtable)
xtable(Tables_Peh, digits = 3)
