#------------------------------------------------------------------------------#
#---------Quasi-experimental methods for climate epidemiology------------------# 
#-------------------------Drexel CCUH Workshop---------------------------------#
#-----------------------------Date:03/14/25------------------------------------#
#--Tarik Benmarhnia (tbenmarhnia@health.ucsd.edu) & Yiqun Ma (yim022@ucsd.edu)-#
#------------------------------------------------------------------------------#

#------------------------------------Checking the directory---------------------
getwd()


#-------------------------------------Installing packages-----------------------


if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
} # a nice package to load several packages simultaneously


p_load("tidyverse","magrittr","broom",        #manipulate data
       "cleaR",                               #clear workspace 
       "here",                                #directory managment
       "Synth", "gsynth",                     #synthetic control
       "panelView", "lme4", "estimatr",       #multi-level modeling
       "ggdag",                               #draw a Directed Acylic Diagram (DAG)
       "gtsummary")                           #for tables

#------------------Mise en place (clear everything)----------------------------
clear() #same as remove(list=ls())





#---------------------------------Loading the data------------------------------

mydata <- read_csv(here("data", "sim_data.csv"))

year_policy <- 2000

mydata <- mydata %>% 
  mutate(year_rec = year - year_policy,
         post     = ifelse(year>=year_policy,1,0),
         treated  = ifelse(state %in% c("Alabama",  "Alaska", 
                                        "Arizona", "Arkansas", "California"), 1,0),
         treatedpost = treated*post)


#-----------------------------------------Analysis------------------------------


#--------------------------------------------------------#
#----------------Preliminary notions----------------------
#-----(pre-post designs, controlled pre-post desings and--
#------Interrupted time series designs)-------------------
#This part is to illustrate other related designs
#---------------------------------------------------------#

#-----------------------#
#----Pre-post designs----
#-----------------------#
#Pre-post analysis (before-after) no control group: one state, two time points

#Subset the data to California and the years 1995 and 2005
#Preview the data
#Plot the data
#Fit a linear model to estimate the effect of the policy on the outcome y
#What are potential problems

##Create the data
dt <- mydata %>%
  filter(state=="California",
         year %in% c(1995, 2005)) 


##Preview the data
head(dt)


##Plot the data
dt %>% 
  ggplot(aes(x=year, y=y, group=state, color = state)) + 
  labs(title = paste("Outcome by year"),
       x = "Year", 
       y = "Outcome",
       colour = "Treatment") +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = year_policy, lty=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 

##fit a linear model using the post as the exposure of interest
fit <- lm(y ~ post + xi + xt + xit, data=dt)
tidy(fit)

# Potential problems
#-------------------#
# i.	The years before and after the policy have been implemented are picked arbitrarily
# ii.	It is impossible to get standard errors
# iii.	It is impossible to adjust for covariates
# iv.	This model would be wrong if the parallel trend assumption is violated
# v.	This model would be wrong if the common schock assumption is violated
# Conclusion 
# As expected the estimate is biased and standard errors inestimable


#----------------------------------#
#----Controlled Pre-post designs----
#----------------------------------#
#2 Pre-post analysis (before-after) with a control group: Two states, two time points


#Subset the data to California and Georgia and the years 1995 and 2005
#Preview the data
#Plot the data
#Fit a linear model to estimate the effect of the policy on the outcome y
#What are potential problems

##Create the data
dt1 <- mydata %>% 
  filter(state=="California" | state=="Georgia",
         year %in% c(1995, 2005))


##Preview the data
head(dt1)


##Plot the data
dt1 %>% 
  ggplot(aes(x=year, y=y, group=state, color = state)) + 
  labs(title = paste("Outcome by year"),
       x = "Year", 
       y = "Outcome",
       colour = "Treatment") +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = year_policy, lty=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 

##fit a linear model using the post as the exposure of interest
fit <- lm(y ~ treated*post, data=dt1) 
summary(fit)
tidy(fit)

# Potential problems 
#--------------------#
# i.	The years before and after the policy have been implemented are picked 
# arbitrarily 
# ii.	The control units are picked arbitrarily
# iii.	It is impossible to get standard errors
# iv.	It is impossible to adjust for covariates
# v.	This model would be wrong if the parallel trend assumption is violated
# vi.	This model would be wrong if the common shock assumption is violated

# Conclusion 
# As expected the estimate is biased and standard errors inestimable



#---------------------------------------#
#----Interrupted Time Series Design----
#---------------------------------------#
###Interupted time series: One state, Multiple time points


#Subset the data to California
#Preview the data
#Plot the data
#Fit a linear model to estimate the effect of the policy on the outcome y
#What are potential problems
mydata
##Create the data
dt2 <- mydata %>% 
  filter(state=="California")


##Preview the data
head(dt2)


##Plot the data
dt2_before %>% 
  ggplot(aes(x=year, y=y, group=state, color = state)) + 
  labs(title = paste("Outcome by year"),
       x = "Year", 
       y = "Outcome",
       colour = "Treatment") +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = year_policy, lty=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 

##Method 1: Fit the model
fit <- glm(y ~ year_rec + post + post*year_rec + xi + xt + xit, data=dt2)
tidy(fit)

##Method 2:  Using a two-stage approach

dt2_before <- mydata %>% 
  filter(state=="California") %>% 
  filter(year_rec < 0)

fit_pre <- glm(y ~ year_rec + xt + xit, data=dt2_before)
tidy(fit_pre)


mydata$y_pre <- predict.glm(fit_pre, newdata = mydata, type = "response")


mydata %>% 
  filter(state=="California") %>% 
  ggplot(aes(x=year, group=state, color = state)) + 
  labs(title = paste("Outcome by year"),
       x = "Year", 
       y = "Outcome",
       colour = "Treatment") +
  geom_line(aes(y=y), color="red") +
  geom_point(aes(y=y), color="red") +
  geom_line(aes(y=y_pre), color="blue") +
  geom_point(aes(y=y_pre), color="blue") +
  geom_vline(xintercept = year_policy, lty=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 

res <- mydata %>% 
  filter(state=="California",
         post==1) %>% 
  summarise(y_trt = mean(y),
            y_crtl = mean(y_pre))

res$y_trt - res$y_crtl

# Potential problem (s)
#----------------------#
# This model would be wrong if the common shock assumption is violated 
# (that is if an event occurs in the treated unit at or after the time of 
# the policy that is not attributable to the policy)

# Conclusion 
# As expected the estimate is slightly biased 


#-------------------------------------------------#
#-----Controlled Interrupted Time Series (CITS)----
#----or Difference-in-Difference Designs-----------
#-------------------------------------------------#

# Manual DID
#1)Obtain the average outcome y among the treated in the pre-policy period
#2)Obtain the average outcome y among the treated in the post-policy period
#3)Obtain the average outcome y among the controls in the pre-policy period
#4)Obtain the average outcome y among the controls in the post-policy period

#Estimate the DID by before/after difference in treated - before/after difference in control
#Estimate the DID by control/treated difference in the post period - control/treated difference in the pre period

df00 <- mydata %>% 
  select(treated, post, y) %>% 
  group_by(treated, post) %>% 
  group_modify(~ {.x %>% map_dfr(mean)}) %>% print(n=Inf)

(250 - 75) - (118-55)

(250 - 118) - (75-55)


#---Step 1: Checking the parallel trends assumption-----

#Check the parallel trends assumptions
#To do so
#restrict the data to the period before the policy
#Run a linear regression of the outcome on time-varying
#covariates and on an interaction term between treated indicator and year
#can also try visualizing the trends

p_load("estimatr") #for the lm_robust() function
pretrend_data <- mydata %>% 
  filter(post == 0)

res_pretrend <-  lm(y ~ treated*year + xit + xt + xi, 
                    data = pretrend_data)
tidy(res_pretrend)


#the lm_robust procedure is best because SE are correctly estimated
res_pretrend <- lm_robust(y ~ treated*year + xit, data = pretrend_data,
                          fixed_effects=state,
                          clusters = state, se_type = "stata")
summary(res_pretrend)
#Pay particular attention to the value and p-value value of
# treated:year. If the coefficient of the interaction is close to 0 or 
#the p-value large, then
# can say that the unit-specific confounders do not have time-varying effects
#in other words, the parallel trend assumption might be okay
#recall that this is only true if we assume that the trend is linear

# If okay, now we can implement the DID analysis method
# Fit a linear model with the treated, post and their interaction
# Don't forget to add any other time-varying covariates

year_policy = 2000
#Visual trends
mydata %>% 
  group_by(year, treated) %>% 
  summarise(y=mean(y),.groups="keep") %>% 
  ggplot(aes(x=year, y=y, group=treated, color = factor(treated))) + 
  annotate("rect", fill = "gray", alpha = 0.5,
           xmin = 2000, xmax = Inf,
           ymin = -Inf, ymax = Inf) +
  labs(title = paste("Outcome by year"),
       x = "Year", 
       y = "Outcome",
       colour = "Treatment") +
  geom_line() +
  scale_color_discrete(labels=c("Controls", "Treated")) +
  geom_vline(xintercept = year_policy, lty=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 

#---Step 2: Implementing the DID analysis method-----
# Fit a linear model with the treated, post and their interaction
# Don't forget to add any other time-varying covariates
# What is the regression coefficient of the interaction term


#Method 1: Lm_robust()
p_load("estimatr")
dta <- lm_robust(y ~ treatedpost + factor(year) + xit, 
                 data = mydata,
                 fixed_effects=state,
                 clusters = state, 
                 se_type = "stata")

dta

did <- round(data.frame(ATT     = dta$coefficients["treatedpost"], 
                        se      = dta$std.error["treatedpost"],
                        low_ci  = dta$conf.low["treatedpost"],
                        high_ci = dta$conf.hig["treatedpost"]),2)
did


#Method 2: lmer(): Multilevel
p_load("lmerTest", "lme4")
fit1 <- lmerTest::lmer(y ~  treated + post + treated:post + xt + xi + xit + (1| state) + (1 | year), 
                       data = mydata, 
                       REML = T)

p_load("broom.mixed", "flextable")
broom.mixed::tidy(fit1) %>% mutate_if(is.numeric, round, 2) %>% flextable()





