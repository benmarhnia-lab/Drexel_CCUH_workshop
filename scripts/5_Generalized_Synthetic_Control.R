#-------------------------------------------------------------------------------#
#---------Quasi-experimental methods for climate epidemiology-------------------# 
#-------------------------Drexel CCUH Workshop----------------------------------#
#-----------------------------Date:03/21/25-------------------------------------#
#--Tarik Benmarhnia (tbenmarhnia@health.ucsd.edu) & Yiqun Ma (yim022@ucsd.edu)--#
#-------------------------------------------------------------------------------#
#--------------------code written by Lara Schwarz-------------------------------#
#paper: Applying a two-stage generalized synthetic control approach to quantify-# 
#----------the heterogeneous health effects of extreme weather events-----------#
#-----------A 2018 large wildfire in California event as a case study-----------#


#-------------------------------------Load packages----------------------
libraries <- c("panelView", "dplyr", "gsynth", "ggplot2", "cowplot", "rstudioapi", 
               "knitr","coefplot", "sjPlot", "svglite","sf", "ggmap", "broom", 
               "ggpubr", "haven", "jtools", "stringr", "metafor", "broom.mixed", 
               "tidyverse", "devtools", "lubridate", "aweek", "padr", "tsibble", 
               "readr", "usethis", "gitcreds","MetBrewer", "choroplethrMaps", 
               "choroplethrZip", "meta")
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
lapply(libraries, install_if_missing)


#----------------------------Loading the fictive data------------------------------
CA_hosp_County <- read_csv(here("data", "CA_hosp_County_2018.csv"))

#---------------------------------Data prep------------------------------
#define study period: from Sept 13 to Dec 5, 2018
CA_hosp_County_week37_49 <- CA_hosp_County[which(CA_hosp_County$week>36 & CA_hosp_County$week<49),]

#define exposure period: from Nov 8 to Dec 5, 2018
CA_hosp_County_week37_49$Exp<-ifelse(CA_hosp_County_week37_49$week<36, 0,
                                     ifelse(CA_hosp_County_week37_49$week>44, 1, NA))

CA_hosp_County_week45_49<-CA_hosp_County[which(CA_hosp_County$week>44 & CA_hosp_County$week<49),]

#define non-exposure period (before wildfires started): from Sept 13 to Nov 7, 2018
CA_hosp_County_week37_44<-CA_hosp_County[which(CA_hosp_County$week>36 & CA_hosp_County$week<45),]

#average smoke exposure during exposure period
Counties_exp <- CA_hosp_County_week45_49 %>%
  group_by(COUNTY_1) %>%
  summarize(smoke_week45_49=mean(smoke))

#average smoke exposure before wildfires started was very low (less <2 on average)
Counties_exp_beforeSmoke <- CA_hosp_County_week37_44 %>%
  group_by(COUNTY_1) %>%
  summarize(smoke_week37_44=mean(smoke))

CA_hosp_County_week37_49 <- left_join(CA_hosp_County_week37_49, Counties_exp,
                                       by = c("COUNTY_1"))

##define smoke exposure level
CA_hosp_County_week37_49$Smokebin<-ifelse(CA_hosp_County_week37_49$smoke_week45_49>=10, 1,
                                          ifelse(CA_hosp_County_week37_49$smoke_week45_49<10, 0, NA))

CA_hosp_County_week37_49_non_exp<-CA_hosp_County_week37_49[which(CA_hosp_County_week37_49$Smokebin==0),]
CA_hosp_County_week37_49_exp<-CA_hosp_County_week37_49[which(CA_hosp_County_week37_49$Smokebin==1),]

# Counties that are exposed
CA_counties_exposed<- CA_hosp_County_week37_49_exp %>%
  group_by(countyname, COUNTY_1) %>%
  summarize(Exposure=mean(Smokebin))

CA_hosp_County_week37_49$week<- as.numeric(CA_hosp_County_week37_49$week)

# create dataset with sum resp by County
CA_hosp_County_week37_49_byweek <- CA_hosp_County_week37_49 %>%
  group_by(week, COUNTY_1, Smokebin, countyname) %>%
  summarize(resp= sum(resp), tmax=mean(tmax), prec=mean(prec), hum=mean(hum), shrtwv_rad=mean(shrtwv_rad), wind=mean(wind))

# creating final dataset for gsynth analysis
df<-CA_hosp_County_week37_49_byweek


#---------------------------------Generalized synthetic control------------------------------
exposed_counties <- unique(CA_counties_exposed$COUNTY_1)

df_final <-data.frame() 

for (i in 1:length(exposed_counties)) {
  ### create new exposure indicator
  ### 0 for all weeks prior to fire week and 1 only if exposed and during fire week
  df$exp_weekly <- ifelse(df$Smokebin==1& df$week>=45,1,0)
  
  ## keep only one exposed county for test
  df2 <- df[which(df$COUNTY_1==exposed_counties[i]|  df$Smokebin==0),]
  
  ### check number of unique zip codes
  length(unique(df2$COUNTY_1))
  
  ### factor zip id
  df2$COUNTY_1 <- as.factor(df2$COUNTY_1)
  
  ### implement generalized synthetic control
  gsynth.out <- gsynth(resp ~ exp_weekly + tmax+prec + hum + shrtwv_rad + wind , data = df2, index = c("countyname","week"), 
                       nboots = 500, inference = "parametric", se = TRUE, parallel = TRUE)
  ### view ATT by period
  ATT<-gsynth.out$att.avg
  SE<-gsynth.out$est.avg[,2]
  
  county<-exposed_counties[i]
  
  result<-cbind(ATT, SE, county)
  
  df_final = rbind(df_final, result)
  
  ### gaps plot - difference between treated and counterfactual
  gap_plot<-plot(gsynth.out, type = "gap" , xlab = "Week", ylab="ATT")
  
  ### treated average and estimated counterfactual average outcomesraw = "none", main="",
  plot(gsynth.out, type = "counterfactual",  xlab = "Week", ylab="Average Hospitalizations")

}

