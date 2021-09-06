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

FataRate_Case1_ON_comp <- ON_death$cases / c(ON_CumCase$Case[1:14], ON_CumCase$Case[-((length(ON_CumCase$Case)-13):length(ON_CumCase$Case))] ) * 100

FataRate_Case1_BC_comp <- BC_death$cases / c(BC_CumCase$Case[1:14], BC_CumCase$Case[-((length(BC_CumCase$Case)-13):length(BC_CumCase$Case))] ) * 100

FataRate_Case1_QC_comp <- QC_death$cases / c(QC_CumCase$Case[1:14], QC_CumCase$Case[-((length(QC_CumCase$Case)-13):length(QC_CumCase$Case))] ) * 100

FataRate_Case1_AB_comp <- AB_death$cases / c(AB_CumCase$Case[1:14], AB_CumCase$Case[-((length(AB_CumCase$Case)-13):length(AB_CumCase$Case))] ) * 100

FataRate_Case1_ON <- FataRate_Case1_ON_comp[1:(length(FataRate_Case1_ON_comp)-5)]

FataRate_Case1_BC <- FataRate_Case1_BC_comp[1:(length(FataRate_Case1_BC_comp)-5)]

FataRate_Case1_QC <- FataRate_Case1_QC_comp[1:(length(FataRate_Case1_QC_comp)-5)]

FataRate_Case1_AB <- FataRate_Case1_AB_comp[1:(length(FataRate_Case1_AB_comp)-5)]


# 3 Canada Data - Fatality Rate Based Study --------------------------------------------------------------------

# 3.1 British Columbia Data --------------------------------------------------------------------

Date_BC <- BC_Confirm$count

tseries::adf.test(FataRate_Case1_BC)
tseries::adf.test(diff(FataRate_Case1_BC))

produce_pacf_plot(FataRate_Case1_BC, "British Columbia - no differencing", defnumber = 1)

produce_pacf_plot(diff(FataRate_Case1_BC), "British Columbia", defnumber = 1)

find_optimal_lag(FataRate_Case1_BC)

find_optimal_lag(diff(FataRate_Case1_BC))


### 3.1.1 BC - AR(3)  #### 

Fitted_Models <- diff_fit (tseries = FataRate_Case1_BC, lag = 3, sigma_e = 0.01, sigma_u = 0.002) 

# xtable(Fitted_Models$Table1)

Prediction_prov <- prediction_diff(FataRate_Case1_BC, Models = Fitted_Models, Province = "British Columbia", fullseries = FataRate_Case1_BC_comp, origindate = BC_CumCase$Date[1])

# plot_results(Prediction_prov, defnumber = 1, diff = T, yrange = c(0,2.5))

Pe_h_est_BC_AR3 <- Evaluations_data(Prediction_prov)


### 3.1.2 BC - AR(4)  #### 

Fitted_Models <- diff_fit (tseries = FataRate_Case1_BC, lag = 4, sigma_e = 0.01, sigma_u = 0.002) 

xtable(Fitted_Models$Table1, digits = 3)

Prediction_prov <- prediction_diff(FataRate_Case1_BC, Models = Fitted_Models, Province = "British Columbia", fullseries = FataRate_Case1_BC_comp, origindate = BC_CumCase$Date[1])

plot_results(Prediction_prov, defnumber = 1, diff = T, yrange = c(0.5,2.1))

Pe_h_est_BC_AR4 <- Evaluations_data(Prediction_prov)


### 3.1.3 BC - AR(4)  #### 

Fitted_Models <- nodiff_fit (tseries = FataRate_Case1_BC, lag = 4, sigma_e = 0.01, sigma_u = 0.002) 

xtable(Fitted_Models$Table1, digits = 3)

Prediction_prov <- prediction_nodiff(FataRate_Case1_BC, Models = Fitted_Models, Province = "British Columbia", fullseries = FataRate_Case1_BC_comp, origindate = BC_CumCase$Date[1])

plot_results(Prediction_prov, defnumber = 1, yrange = c(0.5,2.1))

# Pe_h_est_BC_AR4 <- Evaluations_data(Prediction_prov)




####  (--------------------------------------------------------------)  ####
# 3.2 Ontario Data --------------------------------------------------------

Date_ON <- ON_Confirm$count

tseries::adf.test(FataRate_Case1_ON)
tseries::adf.test(diff(FataRate_Case1_ON))

find_optimal_lag(diff(FataRate_Case1_ON))

produce_pacf_plot(FataRate_Case1_ON, "Ontario - no differencing", defnumber = 1)

produce_pacf_plot(diff(FataRate_Case1_ON), "Ontario", defnumber = 1)

### 3.2.1 ON-AR(1)  ####

Fitted_Models <- diff_fit (tseries = FataRate_Case1_ON, lag = 1, sigma_e = 0.005, sigma_u = 0.002) 

# xtable(Fitted_Models$Table1)

Prediction_prov <- prediction_diff(FataRate_Case1_ON, Models = Fitted_Models, Province = "Ontario", fullseries = FataRate_Case1_ON_comp, origindate = ON_CumCase$Date[1])

# plot_results(Prediction_prov, defnumber = 1, diff = T, yrange = c(0,2.5))

Pe_h_est_ON_AR1 <- Evaluations_data(Prediction_prov)


### 3.2.2 ON-AR(2)  (chosen) ####

Fitted_Models <- diff_fit (tseries = FataRate_Case1_ON, lag = 2, sigma_e = 0.005, sigma_u = 0.002) 

xtable(Fitted_Models$Table1, digits = 3)

Prediction_prov <- prediction_diff(FataRate_Case1_ON, Models = Fitted_Models, Province = "Ontario", fullseries = FataRate_Case1_ON_comp, origindate = ON_CumCase$Date[1])

plot_results(Prediction_prov, defnumber = 1, diff = T, yrange = c(1,4))

Pe_h_est_ON_AR2 <- Evaluations_data(Prediction_prov)



# !-----------------------------------------             ------------------
# 3.3 Quebec Data --------------------------------------------------------

Date_QC <- QC_Confirm$count

tseries::adf.test(FataRate_Case1_QC)
tseries::adf.test(diff(FataRate_Case1_QC))

produce_pacf_plot(FataRate_Case1_QC, "Quebec - no differencing", defnumber = 1)

produce_pacf_plot(diff(FataRate_Case1_QC), "Quebec", defnumber = 1)

find_optimal_lag(diff(FataRate_Case1_QC))

### 3.3.1 QC-AR(3) ####

Fitted_Models <- diff_fit (tseries = FataRate_Case1_QC, lag = 3, sigma_e = 0.002, sigma_u = 0.002) 

# xtable(Fitted_Models$Table1)

Prediction_prov <- prediction_diff(FataRate_Case1_QC, Models = Fitted_Models, Province = "Quebec", fullseries = FataRate_Case1_QC_comp, origindate = QC_CumCase$Date[1])

# plot_results(Prediction_prov, defnumber = 1, diff = T, yrange = c(0.5,7))

Pe_h_est_QC_AR3 <- Evaluations_data(Prediction_prov)

### 3.3.2 QC-AR(4) (chosen)####

Fitted_Models <- diff_fit (tseries = FataRate_Case1_QC, lag = 4, sigma_e = 0.002, sigma_u = 0.002) 

xtable(Fitted_Models$Table1, digits = 3)

Prediction_prov <- prediction_diff(FataRate_Case1_QC, Models = Fitted_Models, Province = "Quebec", fullseries = FataRate_Case1_QC_comp, origindate = QC_CumCase$Date[1])

plot_results(Prediction_prov, defnumber = 1, diff = T, yrange = c(1.5,5))

Pe_h_est_QC_AR4 <- Evaluations_data(Prediction_prov)



# !-----------------------------------------             ------------------
# 3.3 Alberta Data --------------------------------------------------------
Date_AB <- AB_Confirm$count

tseries::adf.test(FataRate_Case1_AB)
tseries::adf.test(diff(FataRate_Case1_AB))

produce_pacf_plot(FataRate_Case1_AB, "Alberta - no differencing", defnumber = 1)

produce_pacf_plot(diff(FataRate_Case1_AB), "Alberta", defnumber = 1)

find_optimal_lag(diff(FataRate_Case1_AB))


### 3.2.1 AB-AR(1)   ####

Fitted_Models <- diff_fit (tseries = FataRate_Case1_AB, lag = 1, sigma_e = 0.005, sigma_u = 0.001) 

# xtable(Fitted_Models$Table1, digits = 3)

Prediction_prov <- prediction_diff(FataRate_Case1_AB, Models = Fitted_Models, Province = "Alberta", fullseries = FataRate_Case1_AB_comp, origindate = AB_CumCase$Date[1])

# plot_results(Prediction_prov, defnumber = 1, diff = T, yrange = c(0,2))

Pe_h_est_AB_AR1 <- Evaluations_data(Prediction_prov)

### 3.2.1 AB-AR(2)  (chosen) ####

Fitted_Models <- diff_fit (tseries = FataRate_Case1_AB, lag = 2, sigma_e = 0.005, sigma_u = 0.001) 

xtable(Fitted_Models$Table1, digits = 3)

Prediction_prov <- prediction_diff(FataRate_Case1_AB, Models = Fitted_Models, Province = "Alberta", fullseries = FataRate_Case1_AB_comp, origindate = AB_CumCase$Date[1])

plot_results(Prediction_prov, defnumber = 1, diff = T, yrange = c(0.3,1.6))

Pe_h_est_AB_AR2 <- Evaluations_data(Prediction_prov)

# 6 Summarize Results --------------------------------------------------------------------
## Table of British Columbia

Table_BC <- get_Peh_table_prov(Pe_h_est_list = list(Pe_h_est_BC_AR3, Pe_h_est_BC_AR4))
Table_ON <- get_Peh_table_prov(Pe_h_est_list = list(Pe_h_est_ON_AR1, Pe_h_est_ON_AR2))
Table_QC <- get_Peh_table_prov(Pe_h_est_list = list(Pe_h_est_QC_AR3, Pe_h_est_QC_AR4))
Table_AB <- get_Peh_table_prov(Pe_h_est_list = list(Pe_h_est_AB_AR1, Pe_h_est_AB_AR2))

Tables_Peh <- rbind(Table_BC, Table_ON, Table_QC, Table_AB)


xtable(Tables_Peh, digits = 3)
