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
       "did"                                  #for did w/ staggered adoption
)                           

#------------------Mise en place (clear everything)------------------------------
clear() #same as remove(list=ls())




#------------------------------------Load the data------------------------------

sim_data_hte_staggered <- read_csv(here("data", "sim_data_hte_staggered.csv"))



#-------------------------Some info about the data------------------------------

# 50 states
# 15 are treated, 35 are controls
# staggered adoption

# state.name[1:5] ("Alabama", "Alaska", "Arizona", "Arkansas", "California") 
#enacted policy in 2000

# state.name[6:10] ("Colorado", "Connecticut", "Delaware", "Florida", "Georgia") 
#enacted policy in 2003

# state.name[6:10] ("Hawaii", "Idaho", "Illinois", "Indiana", "Iowa") 
#enacted policy in 2006

#-----------------------------------Visualize the data--------------------------
mydata = sim_data_hte_staggered

p_load("panelView")
panelview(y ~ treatedpost, data = mydata,
          index = c("state","year"), 
          pre.post = TRUE) 


panelview(y ~ treatedpost, data = mydata, 
          index = c("state","year"), 
          type = "outcome",  
          by.group = TRUE)


#For all the data
mydata %>% 
  ggplot(aes(x=year, y=y, group=state)) + 
  annotate("rect", fill = "gray", alpha = 0.5,
           xmin = 2000, xmax = Inf,
           ymin = -Inf, ymax = Inf) +
  labs(title = paste("Outcome by year"),
       x = "Year", 
       y = "Outcome",
       color = "Treatment") +
  geom_line(aes(color=factor(treated)), size=0.5) +
  scale_color_discrete(labels=c("Control", "Treated")) +
  geom_vline(xintercept = 2000, lty=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 



#---------------Analyze the data: Callaway and Sant’Anna------------------------
#Tutorial here: https://cran.r-project.org/web/packages/did/vignettes/did-basics.html
#Paper here: https://www.sciencedirect.com/science/article/pii/S0304407620303948?via%3Dihub

#load package
p_load(did)

#Create new variables
mydata1 <- mydata %>% 
  mutate(first_treated = case_when( (state %in% state.name[1:5])   ~ 2000,
                                    (state %in% state.name[6:10])  ~ 2003,
                                    (state %in% state.name[11:15]) ~ 2006,
                                    TRUE~0))



group_time_effects <- att_gt( yname  = "y",
                              tname  = "year",
                              idname = "state_num",
                              gname  = "first_treated",
                              xformla = ~ xit,
                              data = mydata1)


#Group-time effects
summary(group_time_effects)
ggdid(group_time_effects)

#Simple Aggregation
agg.simple <- aggte(group_time_effects, type = "simple")
summary(agg.simple)

#Dynamic Effects and (Event Studies): Effect by length of exposure
agg.es <- aggte(group_time_effects, type = "dynamic")
summary(agg.es)
ggdid(agg.es)

#Group-Specific Effects : : Effect by group
agg.gs <- aggte(group_time_effects, type = "group")
summary(agg.gs)
ggdid(agg.gs)

#Calendar Time effects
agg.ct <- aggte(group_time_effects, type = "calendar")
summary(agg.ct)
ggdid(agg.ct)

# Now that we have seen different effects that could be obtained using the method 
# above, let us see what other models give us.

#-------------------------Analyze the data: Simple DID--------------------------


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


p_load("lmerTest", "lme4")
fit1 <- lmerTest::lmer(y ~ treated:post + 
                         xt + xi + xit + 
                         (1| state) + 
                         (1 | year), 
                       data = mydata, 
                       REML = T)
broom.mixed::tidy(fit1) %>% mutate_if(is.numeric, round, 2) %>% 
  flextable::flextable()
# This quantity above is not giving us the group-specific overall ATT and could also
# be biased.  


#----------------------Analyze the data: Generalized SCM------------------------


y <- gsynth(y ~ treatedpost + xit, 
            data = mydata,  
            EM = F, 
            index = c("state","year"), 
            inference = "parametric", 
            se = TRUE,
            nboots = 100,  #so that it can run faster, default is 200
            r = c(0, 5), 
            CV = TRUE, 
            seed = 123,
            force = "two-way", 
            parallel = FALSE)

y1 <- round(data.frame(y$est.avg),2)


#Period-specific ATT
y$est.att

#average ATT
y$est.avg

#Weights from the generalized synthetic controls
y[["alpha.co"]]

plot(y, type = "counterfactual", raw = "none", main="")

plot(y, type = "counterfactual", raw = "band", main="")
plot(y, type = "counterfactual", raw = "all")
plot(y, type = "ct", raw = "none", main = "", shade.post = FALSE)



#-----Estimate individual state (or time) effects and then Pooled the results----
##Method 1: Estimate single effects and pooled effects using the DID method----


## Group 1: "Early adopters", policy enacted in 2000
data_filter_g1 <- mydata %>% 
  filter(state %in% state.name[1:5] | treated==0)

dta_g1 <- lm_robust(y ~ treatedpost + factor(year) + xit, 
                    data = data_filter_g1,
                    fixed_effects=state,
                    clusters = state, 
                    se_type = "stata")

did_g1 <- data.frame(group ="1-Early_adopters",
                     att     = dta_g1$coefficients["treatedpost"], 
                     se      = dta_g1$std.error["treatedpost"],
                     lowerci = dta_g1$conf.low["treatedpost"],
                     upperci = dta_g1$conf.hig["treatedpost"],
                     row.names = NULL)
did_g1


##Group 2: "Medium adopters", policy enacted in 2003

data_filter_g2 <- mydata %>% 
  filter(state %in% state.name[6:10] | treated==0)

dta_g2 <- lm_robust(y ~ treatedpost + factor(year) + xit, 
                    data = data_filter_g2,
                    fixed_effects=state,
                    clusters = state, 
                    se_type = "stata")

did_g2 <- data.frame(group = "2-Mid-adopters",
                     att     = dta_g2$coefficients["treatedpost"], 
                     se      = dta_g2$std.error["treatedpost"],
                     lowerci = dta_g2$conf.low["treatedpost"],
                     upperci = dta_g2$conf.hig["treatedpost"],
                     row.names = NULL)
did_g2


##Group 3: "ate adopters", policy enacted in 2006


data_filter_g3 <- mydata %>% 
  filter(state %in% state.name[11:15] | treated==0)

dta_g3 <- lm_robust(y ~ treatedpost + factor(year) + xit, 
                    data = data_filter_g3,
                    fixed_effects=state,
                    clusters = state, 
                    se_type = "stata")

did_g3 <- data.frame(group = "3-Late-adopters",
                     att    = dta_g3$coefficients["treatedpost"], 
                     se      = dta_g3$std.error["treatedpost"],
                     lowerci = dta_g3$conf.low["treatedpost"],
                     upperci = dta_g3$conf.hig["treatedpost"],
                     row.names = NULL)
did_g3

combined <- bind_rows(did_g1, did_g2, did_g3)
combined

p_load("metafor")
metaresult <- rma(yi = att, 
                  sei = se, 
                  data = combined, 
                  slab = group, 
                  method = "ML")


combined2 <- combined %>% 
  add_row(group="4-Overall",
          att=as.vector(metaresult$beta),
          se = metaresult$se,
          lowerci = metaresult$ci.lb,
          upperci = metaresult$ci.ub)

combined2



## Method 2: Estimate single effects and pooled effects using the generalized SCM----


#Create a function to estimate the effect of each state
gsynth_meta <- function(states, data, nboots = 200){
  
  
  dat <- data %>%
    filter(state=={{states}} | treated == 0)
  
  y <- gsynth(y ~ treatedpost + xit, 
              data = dat,  
              EM = F, 
              index = c("state","year"), 
              inference = "parametric", 
              se = TRUE,
              nboots = nboots, 
              r = c(0, 5), 
              CV = TRUE, 
              seed = 123,
              force = "two-way", 
              parallel = FALSE)
  
  y1 <- data.frame(y$est.avg)
  
  res <- tibble(
    state = states, 
    ATT=y1$Estimate,
    SE=y1$S.E.,
    lowerci=y1$CI.lower,
    upperci=y1$CI.upper)
  
  res
  
  return(res)
}

#test the function
gsynth_meta(states="California",
            nboots = 100,
            data=mydata)


#create the list and their intersection
states <- mydata %>% 
  filter(treated == 1) %>% 
  select(state) %>% 
  distinct() %>% 
  as.matrix() %>% 
  as.vector()
states


#create list of combinations
list_states <- expand_grid(states) %>% print(n=Inf)
list_states


#Loop through all the states
gsynth_res <- pmap_dfr(list_states, gsynth_meta, nboots=2, data =mydata)
#used nboots = 2 for simplicity. You should aim for at least 200 bootstrap samples
gsynth_res


#estimate the overall effect
metaresult2 <- rma(yi = ATT, 
                   sei = SE, 
                   data = gsynth_res, 
                   slab = state, 
                   method = "ML")


#combine the datasets
combined3 <- gsynth_res %>% 
  add_row(state="Overall",
          ATT=as.vector(metaresult2$beta),
          SE = metaresult2$se,
          lowerci = metaresult2$ci.lb,
          upperci = metaresult2$ci.ub)

combined3

## Method 3: Estimate single effects and pooled effects using the Augmented Synthetic Control Methods----
p_load_gh("ebenmichael/augsynth")
set.seed(123)
augsynth.scm <-multisynth(y ~ treatedpost | xit,
                          unit = state, 
                          time = year, 
                          data = mydata, 
                          
                          fixedeff = T,
                          scm = T,
                          time_cohort = F, #can change this to T if interested in time cohorts instead of unit effects
                          progfunc="Ridge")
#progfunc = None for Traditional SCM,
#progfunc = Ridge for Ridge regression or augmented SCM
#progfunc = GSYN for the Generalized SCM

augsynth.scm 

res <- summary(augsynth.scm)
res

res$att
#is a dataframe that contains all of the point estimates, standard errors, and lower/upper confidence limits. Time = NA denotes the effect averaged across the post treatment periods.

plot(augsynth.scm)
plot(augsynth.scm, levels = "Average")



