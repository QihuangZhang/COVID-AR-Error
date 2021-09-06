# 1.1 Packages ------------------------------------------------------------

# withr::with_libpaths(new = "/u/q259zhan/R/x86_64-pc-linux-gnu-library/3.6/", install_github('GuangchuangYu/nCov2019'))

library(tidyr)
require(dplyr)
library(aTSA)
library(forecast)
library(data.table)


# 2. Data --------------------------------------------------------------------
## Global Parameter
source(file="code/COVID_functions.R")

startdate_overall <- "2020-01-23"
startdate_study <- "2020-04-10"

# 2.1 Case Date
DiagnoseData <- fread("data/covidJHUcanada.csv")
# Transform from cumulative data to dayly data
DiagnoseData_dayly <- apply(DiagnoseData[,3:dim(DiagnoseData)[2]], MARGIN = 1, FUN = function(x){
  x1 <- diff(x)
  x1[x1<0] <- 0
  x2 <- imputezero(x1)
  return (x2)
})

Daly_df <- cbind(DiagnoseData[,1:2], t(DiagnoseData_dayly))
names(Daly_df)[3:dim(Daly_df)[2]] <- as.Date(1:(dim(Daly_df)[2]-2), origin="2020-01-23")
Daly_df_full <- gather(Daly_df, Date, count, "18285":"18853") %>%
  filter ( `Province/State` %in% c("Alberta", "British Columbia", "Ontario", "Quebec")  ) %>%
  mutate(Date = as.Date(as.numeric(Date)-18285, origin = "2020-01-23")) %>%
  mutate ( type = "confirmed" ) %>%
  mutate ( Province = `Province/State`) %>%
  select (c("Province", "Date", "count", "type")) 

I0_ON <- sum( Daly_df_full %>%
                filter(Province == "Ontario") %>%
                filter(type == "confirmed") %>% 
                filter(Date <= as.Date (startdate_study)) %>% select("count"))

I0_BC <- sum( Daly_df_full %>%
                filter(Province == "British Columbia") %>%
                filter(type == "confirmed") %>% 
                filter(Date <= as.Date (startdate_study)) %>% select("count"))

I0_QC <- sum( Daly_df_full %>%
                filter(Province == "Quebec") %>%
                filter(type == "confirmed") %>% 
                filter(Date <= as.Date (startdate_study)) %>% select("count"))

I0_AB <- sum( Daly_df_full %>%
                filter(Province == "Alberta") %>%
                filter(type == "confirmed") %>% 
                filter(Date <= as.Date (startdate_study)) %>% select("count"))


Daly_df_long <- Daly_df_full %>%
  filter ( Date > as.Date (startdate_study)) %>%
  filter ( Date < as.Date ("2021-01-15"))

ON_CumCase <- Daly_df_long   %>%
  filter (Province == "Ontario") %>%
  mutate(Case=cumsum(count) + I0_ON) %>%
  select(c("Date","Case")) 

BC_CumCase <- Daly_df_long   %>%
  filter (Province == "British Columbia") %>%
  mutate(Case=cumsum(count) + I0_BC) %>%
  select(c("Date","Case"))  

QC_CumCase <- Daly_df_long    %>%
  filter (Province == "Quebec") %>%
  mutate(Case=cumsum(count) + I0_QC) %>%
  select(c("Date","Case")) 

AB_CumCase <- Daly_df_long    %>%
  filter (Province == "Alberta") %>%
  mutate(Case=cumsum(count) + I0_AB) %>%
  select(c("Date","Case")) 


ON_Confirm <- Daly_df_long   %>% 
  filter (Province == "Ontario")
  
BC_Confirm <- Daly_df_long %>% 
  filter (Province == "British Columbia")

QC_Confirm <- Daly_df_long  %>%
  filter (Province == "Quebec")

AB_Confirm <- Daly_df_long   %>% 
  filter (Province == "Alberta")


# 2.1 Death Date
DeathData <- fread("data/deathJHUcanada.csv")
# Transform from cumulative data to dayly data

DeathData_dayly <- apply(DeathData[,3:dim(DeathData)[2]], MARGIN = 1, FUN = function(x){
  x1 <- diff(x)
  x1[x1<0] <- 0
  x2 <- imputezero(x1)
  return (x2)
})

Dea_df <- cbind(DeathData[,1:2], t(DeathData_dayly))
names(Dea_df)[3:dim(Dea_df)[2]] <- as.Date(1:(dim(Dea_df)[2]-2), origin="2020-01-23")
Dea_df_long <- gather(Dea_df, Date, count, "18285":"18853") %>%
  filter ( `Province/State` %in% c("Alberta", "British Columbia", "Ontario", "Quebec")  ) %>%
  mutate(Date = as.Date(as.numeric(Date)-18285, origin = "2020-01-23")) %>%
  mutate ( type = "death" ) %>%
  mutate ( Province = `Province/State`) %>%
  select (c("Province", "Date", "count", "type")) %>%
  filter ( Date > as.Date (startdate_study)) %>%
  filter ( Date < as.Date ("2021-01-15"))


ON_death <- Dea_df_long %>%
  filter(Province == "Ontario") %>%
  filter(type == "death") %>%
  mutate( cases = cumsum(count))


BC_death <- Dea_df_long %>%
  filter(Province == "British Columbia") %>%
  filter(type == "death") %>%
  mutate( cases = cumsum(count))

QC_death <- Dea_df_long %>%
  filter(Province == "Quebec") %>%
  filter(type == "death") %>%
  mutate( cases = cumsum(count)) 


AB_death <- Dea_df_long %>%
  filter(Province == "Alberta") %>%
  filter(type == "death") %>%
  mutate( cases = cumsum(count))




# 3. Before selected date ------------------------------------------------------
startdate_study2 <- as.Date (startdate_study) - 14

I0_ON2 <- sum( Daly_df_full %>%
                filter(Province == "Ontario") %>%
                filter(type == "confirmed") %>% 
                filter(Date <= startdate_study2) %>% select("count"))

I0_BC2 <- sum( Daly_df_full %>%
                filter(Province == "British Columbia") %>%
                filter(type == "confirmed") %>% 
                filter(Date <= startdate_study2) %>% select("count"))

I0_QC2 <- sum( Daly_df_full %>%
                filter(Province == "Quebec") %>%
                filter(type == "confirmed") %>% 
                filter(Date <= startdate_study2) %>% select("count"))

I0_AB2 <- sum( Daly_df_full %>%
                filter(Province == "Alberta") %>%
                filter(type == "confirmed") %>% 
                filter(Date <= startdate_study2) %>% select("count"))


Daly_df_long2 <- Daly_df_full %>%
  filter ( Date > startdate_study2) %>%
  filter ( Date < as.Date ("2021-01-15"))

ON_CumCase2 <- Daly_df_long2   %>%
  filter (Province == "Ontario") %>%
  mutate(Case=cumsum(count) + I0_ON) %>%
  select(c("Date","Case")) 

BC_CumCase2 <- Daly_df_long2  %>%
  filter (Province == "British Columbia") %>%
  mutate(Case=cumsum(count) + I0_BC) %>%
  select(c("Date","Case"))  

QC_CumCase2 <- Daly_df_long2    %>%
  filter (Province == "Quebec") %>%
  mutate(Case=cumsum(count) + I0_QC) %>%
  select(c("Date","Case")) 

AB_CumCase2 <- Daly_df_long2    %>%
  filter (Province == "Alberta") %>%
  mutate(Case=cumsum(count) + I0_AB) %>%
  select(c("Date","Case")) 



# 4. Trajactory Plot ------------------------------------------------------


FataRate_Case1_ON_comp <- ON_death$cases / c(ON_CumCase2$Case[1:14], ON_CumCase$Case[-((length(ON_CumCase$Case)-13):length(ON_CumCase$Case))] ) * 100

FataRate_Case1_BC_comp <- BC_death$cases / c(BC_CumCase2$Case[1:14], BC_CumCase$Case[-((length(BC_CumCase$Case)-13):length(BC_CumCase$Case))] ) * 100

FataRate_Case1_QC_comp <- QC_death$cases / c(QC_CumCase2$Case[1:14], QC_CumCase$Case[-((length(QC_CumCase$Case)-13):length(QC_CumCase$Case))] ) * 100

FataRate_Case1_AB_comp <- AB_death$cases / c(AB_CumCase2$Case[1:14], AB_CumCase$Case[-((length(AB_CumCase$Case)-13):length(AB_CumCase$Case))] ) * 100

FataRate_Case1_ON <- FataRate_Case1_ON_comp[1:(length(FataRate_Case1_ON_comp)-5)]

FataRate_Case1_BC <- FataRate_Case1_BC_comp[1:(length(FataRate_Case1_BC_comp)-5)]

FataRate_Case1_QC <- FataRate_Case1_QC_comp[1:(length(FataRate_Case1_QC_comp)-5)]

FataRate_Case1_AB <- FataRate_Case1_AB_comp[1:(length(FataRate_Case1_AB_comp)-5)]



FataRate_Case2_ON_comp <- ON_death$cases / c(ON_CumCase2$Case[5:14], ON_CumCase$Case[-((length(ON_CumCase$Case)-9):length(ON_CumCase$Case))] ) * 100

FataRate_Case2_BC_comp <- BC_death$cases / c(BC_CumCase2$Case[5:14], BC_CumCase$Case[-((length(BC_CumCase$Case)-9):length(BC_CumCase$Case))] ) * 100

FataRate_Case2_QC_comp <- QC_death$cases / c(QC_CumCase2$Case[5:14], QC_CumCase$Case[-((length(QC_CumCase$Case)-9):length(QC_CumCase$Case))] ) * 100

FataRate_Case2_AB_comp <- AB_death$cases / c(AB_CumCase2$Case[5:14], AB_CumCase$Case[-((length(AB_CumCase$Case)-9):length(AB_CumCase$Case))] ) * 100

FataRate_Case2_ON <- FataRate_Case2_ON_comp[1:(length(FataRate_Case2_ON_comp)-5)]

FataRate_Case2_BC <- FataRate_Case2_BC_comp[1:(length(FataRate_Case2_BC_comp)-5)]

FataRate_Case2_QC <- FataRate_Case2_QC_comp[1:(length(FataRate_Case2_QC_comp)-5)]

FataRate_Case2_AB <- FataRate_Case2_AB_comp[1:(length(FataRate_Case2_AB_comp)-5)]



FataRate_Case3_ON_comp <- ON_death$cases / ON_CumCase$Case * 100

FataRate_Case3_BC_comp <- BC_death$cases / BC_CumCase$Case * 100

FataRate_Case3_QC_comp <- QC_death$cases / QC_CumCase$Case * 100

FataRate_Case3_AB_comp <- AB_death$cases / AB_CumCase$Case * 100

FataRate_Case3_ON <- FataRate_Case3_ON_comp[1:(length(FataRate_Case3_ON_comp)-5)]

FataRate_Case3_BC <- FataRate_Case3_BC_comp[1:(length(FataRate_Case3_BC_comp)-5)]

FataRate_Case3_QC <- FataRate_Case3_QC_comp[1:(length(FataRate_Case3_QC_comp)-5)]

FataRate_Case3_AB <- FataRate_Case3_AB_comp[1:(length(FataRate_Case3_AB_comp)-5)]



AllReport_ON <- data.frame(
  Date = ON_death$Date[1:(length(ON_death$Date)-5)],
  ReportDR = c(FataRate_Case1_ON, FataRate_Case2_ON, FataRate_Case3_ON),
  Definition = rep(c("Definition 1", "Definition 2", "Definition 3"), each = length(FataRate_Case1_ON)),
  Province = "Ontario"
)


AllReport_BC <- data.frame(
  Date = BC_death$Date[1:(length(BC_death$Date)-5)],
  ReportDR = c(FataRate_Case1_BC, FataRate_Case2_BC, FataRate_Case3_BC),
  Definition = rep(c("Definition 1", "Definition 2", "Definition 3"), each = length(FataRate_Case1_BC)),
  Province = "British Columbia"
)

AllReport_QC <- data.frame(
  Date = QC_death$Date[1:(length(QC_death$Date)-5)],
  ReportDR = c(FataRate_Case1_QC, FataRate_Case2_QC, FataRate_Case3_QC),
  Definition = rep(c("Definition 1", "Definition 2", "Definition 3"), each = length(FataRate_Case1_QC)),
  Province = "Quebec"
)

AllReport_AB <- data.frame(
  Date = AB_death$Date[1:(length(AB_death$Date)-5)],
  ReportDR = c(FataRate_Case1_AB, FataRate_Case2_AB, FataRate_Case3_AB),
  Definition = rep(c("Definition 1", "Definition 2", "Definition 3"), each = length(FataRate_Case1_AB)),
  Province = "Alberta"
)

AllReport <- rbind(AllReport_ON, AllReport_BC, AllReport_QC, AllReport_AB)

AllReport$Province <- factor(AllReport$Province, levels=c('British Columbia','Ontario','Quebec','Alberta'))

library(ggplot2)
require(ggthemes)
pdf(file = "output/revision/AllReportDeath.pdf",height = 5, width = 16)
strip_color <- "#0A1D37"
AllReportD <- ggplot(AllReport)  +
  geom_line(aes(x=Date, y=ReportDR, color = Definition), size=1.2, alpha= 0.8)+
  # geom_point(aes(x=Date, y=ReportDR, color = Definition, shape= Definition),size=2) + 
  scale_color_manual(values = c("#C4878C", "#469BEC","#6D882B"))  +
  # scale_shape_manual(values = c(15, 17, 19))  +
  labs(x = "Date", y = "Mortality Rate (%)") +
  theme(text=element_text(size=13, family="mono"))+
  labs(color = "Definition of\nDeath Rate", shape = "Definition of\nDeath Rate") +
  facet_wrap(~Province, nrow = 1)  +
  theme_bw()+
  theme(strip.background =element_rect(fill=strip_color,color=strip_color))+ # #535b44
  theme(strip.text = element_text(colour = 'white')) +
  theme(panel.border = element_rect(colour = strip_color))+
  theme(text=element_text(size=16, family="URWHelvetica"))

print(AllReportD)


dev.off()

