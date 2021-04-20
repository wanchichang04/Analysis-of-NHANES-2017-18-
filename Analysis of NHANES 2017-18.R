#### Part 2: Analysis of NHANES 2017-18
#### Wan Chi Chang
#### wanchichang04@tufts.edu
#### Start Date: 2020-11-30
#### Last edit: 2020-12-6, (1)check interpretation of dDev
#                          (2)Conduct stratified analysis, and add plot stratified by sex

##Set up workspace--------------------------------------------
#clear workspace
rm(list=ls())
cat("\014")
dev.off()

#set working directory
setwd("~/Desktop/Grad school/2020 Fall Tufts/Int bio/Final project/Output")

#install packages
library(tidyverse)
library(survey)
library(broom) #Save model regression information
library(car) #CRP plot, vif
library(DescTools) #post-hoc test
library(writexl)
library(lmtest)
library(jtools) #summ function
library(ggthemes)
library(LogisticDx) #Pregibon DBeta
library(forcats)

##Data set import and clean--------------------------------------------
nhanes_pro <- read.csv("nhanes_processed.csv")
summary(nhanes_pro)

nhanes <- nhanes_pro %>%
  mutate(age = RIDAGEYR,
         sex = factor(RIAGENDR, levels = c(1,2), labels = c("male","female")),
         race = factor(RIDRETH3, levels = c(3,4,1,2,6,7), 
                       labels = c("white","black","hispanic","hispanic","asian","other")),
         education = factor(DMDEDUC2, levels = c(1,2,3,4,5), 
                            labels = c("under 9th","9-11th","high school","college","college graduate or above"), NA),
         income = factor(ifelse(INDFMIN2==1|INDFMIN2==2|INDFMIN2==3|INDFMIN2==4|INDFMIN2==13, "under 20000",
                         ifelse(INDFMIN2==5|INDFMIN2==6|INDFMIN2==7|INDFMIN2==12, "20000-44999",
                         ifelse(INDFMIN2==8|INDFMIN2==9|INDFMIN2==10, "45000-74999",
                         ifelse(INDFMIN2==14|INDFMIN2==15, "75000 and above", NA))))),
         psu = SDMVPSU,
         stra = SDMVSTRA,
         weighted = WTDRD1, 
         caff = DR1TCAFF,
         caff_c = factor(ifelse(DR1TCAFF<50,"low",
                         ifelse(DR1TCAFF>=50 & DR1TCAFF<200, "normal",
                         ifelse(DR1TCAFF>=200 & DR1TCAFF<400, "high", 
                         ifelse(DR1TCAFF>=400, "very high", NA)))),
                         levels = c("normal","low","high","very high")),
         fiber = DR1TFIBE,
         fat = DR1TTFAT,
         fat_c = factor(ifelse(DR1TTFAT>80, "high",
                      ifelse(DR1TTFAT<=80 & DR1TTFAT>35, "normal",
                      ifelse(DR1TTFAT<=35, "low", NA))),
                      levels = c("normal","low","high")),
         chol = DR1TCHOL,
         vitc = DR1TVC,
         water = `DR1.320Z`,
         alcohol = factor(ifelse(DR1TALCO<=14, "less than 1 drink",
                          ifelse(DR1TALCO>14&DR1TALCO<=28, "1-2 drinks",
                          ifelse(DR1TALCO>28, "more than 2 drinks", NA))),
                          levels = c("less than 1 drink", "1-2 drinks","more than 2 drinks")),
         wt = BMXWT,
         bmi = BMXBMI,
         bmi_c = factor(ifelse(BMXBMI < 18.5, "Underweight",
                        ifelse(BMXBMI>=18.5 & BMXBMI < 25, "Normal",
                        ifelse(BMXBMI>=25 & BMXBMI<30, "Overweight",
                        ifelse(BMXBMI>=30, "Obese", NA)))),
                        levels = c("Normal","Underweight","Overweight","Obese")),
         waist = BMXWAIST,
         MCQ550 = factor(ifelse(is.na(MCQ550),0,MCQ550)),
         MCQ560 = factor(ifelse(is.na(MCQ560),0,MCQ560)),
         gallstone = factor(MCQ550, levels = c(0,1,2), labels = c(0,"Yes","No"), NA),
         gallsur = factor(MCQ560, levels = c(0,1,2), labels = c(0,"Yes","No"), NA),
         galldis = factor(ifelse(MCQ550==1|MCQ560==1, "Yes",
                          ifelse((MCQ550==0|MCQ550==2) & MCQ560==2, "No",
                          ifelse((MCQ560==0|MCQ560==2) & MCQ550==2, "No",NA))), 
                          levels = c("No","Yes")),
         diabetes = factor(DIQ010,levels = c(1,2), labels = c("Yes","No"),NA),
         smoking = factor(ifelse(SMQ020==1&(SMQ040==1|SMQ040==2), "Yes", 
                          ifelse(SMQ020==2|(SMQ040==3), "No", NA))),
         pregnant = factor(ifelse(RHD143==1&RIAGENDR==2, "Yes", 
                           ifelse(RIAGENDR==1|RHD143==2,"No", NA))),
         lab_chol = LBXTC,
         PAQ605 = ifelse(is.na(PAQ605),0,PAQ605), #recode them as zero bc ifelse will ignore NA
         PAQ620 = ifelse(is.na(PAQ620),0,PAQ620),
         PAD615 = ifelse(is.na(PAD615),0,PAD615),
         PAD630 = ifelse(is.na(PAD630),0,PAD630),
         PAQ650 = ifelse(is.na(PAQ650),0,PAQ650),
         PAQ665 = ifelse(is.na(PAQ665),0,PAQ665),
         PAD660 = ifelse(is.na(PAD660),0,PAD660),
         PAD675 = ifelse(is.na(PAD675),0,PAD675),
         pa_work_c = factor(ifelse(PAQ605==1|PAQ620==1,"Yes",
                            ifelse(PAQ605==2&PAQ620==2, "No", NA))),
         pa_work_t = factor(ifelse((PAD615<1000 & PAD615>=75)|(PAD630 < 1000 & PAD630>=150), "AHA",
                            ifelse(PAD615<75&PAD630<150, "non-AHA", NA))),
         pa_work = factor(ifelse(pa_work_c=="Yes" & pa_work_t=="AHA", "enough",
                          ifelse(pa_work_c=="Yes" & pa_work_t=="non-AHA", "not enough",
                          ifelse(pa_work_c=="No" & pa_work_t=="non-AHA", "don't exercise", NA)))),
         pa_act_c = factor(ifelse(PAQ650==1|PAQ665==1, "Yes",
                          ifelse(PAQ650==2&PAQ665==2, "No", NA))),
         pa_act_t = factor(ifelse((PAD660<1000 & PAD660>=75)|(PAD675<1000 & PAD675>=150), "AHA",
                           ifelse(PAD660<75&PAD675<150, "non-AHA",NA))),
         pa_act = factor(ifelse(pa_act_c=="Yes" & pa_act_t=="AHA", "enough",
                         ifelse(pa_act_c=="Yes" & pa_act_t=="non-AHA", "not enough",
                         ifelse(pa_act_c=="No" & pa_act_t=="non-AHA", "don't exercise", NA)))),
         pa = factor(ifelse(pa_work=="enough"|pa_act=="enough", "enough",
                     ifelse((pa_work=="not enough"|pa_work=="don't exercise")& pa_act=="not enough", "not enough",
                     ifelse((pa_act=="not enough"|pa_act=="don't exercise")& pa_work=="not enough", "not enough",
                     ifelse(pa_work=="don't exercise"& pa_act=="don't exercise", "not enough", NA)))),
                     levels = c("not enough","enough"))) %>%
  select(SEQN, age:education, psu:lab_chol, pa)#decide not to include income based on lit review
summary(nhanes)

nhanes_c <- nhanes %>%
  filter(age>=20 & !(diabetes %in% "Yes" & !(pregnant %in% "Yes"))) %>%
  select(SEQN, age:education, diabetes, psu:bmi_c, galldis, smoking, pa) %>% 
  na.omit(.)%>%
  mutate(keep = 1) 

summary(nhanes_c)

nhanes_subset <- nhanes_c %>% select(SEQN, keep)

nhanes_final <- left_join(nhanes, nhanes_subset, by = "SEQN") %>%
  mutate(keep = factor(ifelse(is.na(keep), 0, keep)))

summary(nhanes_final)

##Descriptive table----------------------------------------------------------
#complex survey design
nhanes_wt <- svydesign(id = ~psu, weights = ~weighted, strata = ~stra, data = nhanes_final, nest = T)
options(survey.lonely.psu = "adjust")
summary(nhanes_wt)

#continuous variables, mean and SE
Table1.con <-svyby(~age+wt+bmi+caff+fiber+fat+chol+vitc+water, ~galldis+keep, nhanes_wt, svymean, ci = T)
Table1.con.tot <- as.data.frame(
  svymean(~age+wt+bmi+caff+fiber+fat+chol+vitc+water, na.rm = T, subset(nhanes_wt, keep==1)))

#categorical variables, unweighted count and percentage
#uwtd count
svyby(~galldis, by = ~galldis+sex, nhanes_wt, unwtd.count)#the same as the original data, so we can use dyplr
Table1.cat.c <- nhanes_c %>%
  select(galldis, sex, race, education, caff_c, alcohol, bmi_c, smoking, pa) %>%
  pivot_longer(names_to = "vars", values_to = "values", cols = sex:pa) %>%
  group_by(galldis, vars, values) %>%
  summarise(n = n())%>%
  arrange(vars)

Table1.cat.c.tot <- nhanes_c %>%
  select(galldis, sex, race, education, caff_c, alcohol, bmi_c, smoking, pa) %>%
  pivot_longer(names_to = "vars", values_to = "values", cols = galldis:pa) %>%
  group_by(vars, values) %>%
  summarise(n = n())%>%
  arrange(vars)

#percentage
Table1.cat <- svyby(~sex+race+education+caff_c+alcohol+bmi_c+smoking+pa, ~galldis+keep, nhanes_wt, svymean)
Table1.cat.tot <- as.data.frame(
  svymean(~sex+race+education+caff_c+alcohol+bmi_c+smoking+pa, subset(nhanes_wt, keep==1))*100)
svymean(~galldis, subset(nhanes_wt, keep==1))*100

#save the tables
write_xlsx(Table1.con, "Table1_con.xlsx")
write_xlsx(Table1.con.tot, "Table1_con_tot.xlsx")
write_xlsx(Table1.cat, "Table1_cat.xlsx")
write_xlsx(Table1.cat.tot, "Table1_cat_tot.xlsx")
write_xlsx(Table1.cat.c, "Table1_cat_c.xlsx")
write_xlsx(Table1.cat.c.tot, "Table1_cat_c_tot.xlsx")

##Figure 1. distribution of gallstones across four caffeine intake level----------
#create the table with needed information
dis_df <- svyby(~galldis, ~caff_c+keep, nhanes_wt, svymean)
dis_df <- dis_df %>%
  slice(5:8) %>%
  select(caff_c, galldisNo, galldisYes) %>%
  pivot_longer(names_to = "Levels", values_to = "values", cols = galldisNo:galldisYes) 

f <- factor(c("low","normal","high","very high"))

percent <- c("87%", "13%","89%","11%","86%","14%","88%","12%")
dis_df <- cbind(dis_df, percent)

I.Plot.Dis <- ggplot(dis_df, aes(x = caff_c, y = values, fill = Levels, label = percent)) + 
  geom_bar(aes(x = caff_c, y = values, fill=Levels), stat = "identity")+
  geom_text(position = position_stack(vjust = 0.9),size = 3, fontface = "bold") + 
  scale_x_discrete(limits = fct_inorder(f),
                   labels =c("Low","Normal","High","Very high")) + 
  scale_fill_manual(values=c("bisque1","bisque2"), labels = c("No", "Yes"))+
  labs(x ="Caffeine Intake Categories", y = "Proportion of subjects", fill = "Gallstones")+
  theme_apa()+
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, face = "bold"))

I.Plot.Dis

##Weighted t-test and chi-squared test--------------------------------------
#t-test looping function
sttloop <- function(vars, design){
  out <- svyttest(as.formula(paste0(vars)), subset(design, keep==1))
  out$data.name <- paste(vars) #paste name on the result
  print(out)
}

#check if the function works
sttloop("caff~galldis", nhanes_wt)
svyttest(caff~galldis, subset(nhanes_wt, keep==1)) #it works!

#create a character vector with the combinations
ttest.col <- c("age~galldis", "wt~galldis", "bmi~galldis","caff~galldis","fiber~galldis",
               "fat~galldis","chol~galldis","vitc~galldis","water~galldis")
  
#t-test
ttest_result <- list()
ttest_result<- lapply(ttest.col, sttloop, nhanes_wt)
cat(capture.output(print(ttest_result), file="ttest_result.txt"))


#chisquare looping function
scsloop <- function(vars, design){
  out <- svychisq(as.formula(paste0( "~" , vars)) , subset(design, keep==1))
  out$data.name <- paste(vars)
  print(out)
           }

#check if the function works
scsloop("sex+galldis", nhanes_wt)
svychisq(~sex+galldis, subset(nhanes_wt, keep==1)) #it works!

#create a character vector with the combinations
chisq.col <- c("sex+galldis", "race+galldis", "education+galldis", "caff_c+galldis", 
               "alcohol+galldis", "bmi_c+galldis", "smoking+galldis", "pa+galldis", "fat_c+galldis")

#chisq
chisq_result <- list()
chisq_result<- lapply(chisq.col, scsloop, nhanes_wt)
cat(capture.output(print(chisq_result), file="chisq_result.txt"))

##Full model--------------------------------------------------------------------
#caffiene as category
glm_full <- svyglm(galldis~caff_c+sex+race+age+bmi+fiber+fat+pa, subset(nhanes_wt, keep==1), 
                   family = binomial()) 
summary(glm_full)
summ(glm_full)#function that gives out AIC and pseudo R2

coef_glm_full <- tidy(glm_full, conf.int = TRUE) %>%
  mutate(OR = exp(estimate), 
         ORL = exp(conf.low),
         ORH = exp(conf.high),
         model = "full")

#anova F-test
glm_pa <- svyglm(galldis~caff_c+sex+race+age+bmi+fiber+fat, subset(nhanes_wt, keep==1), 
                 family = binomial())
summary(glm_pa)
anova(glm_full, glm_pa, method = "LRT") #p=0.47

#write a function for looping
regtest <- function(vars, model){
  out <- regTermTest(glm_full, as.formula(paste0(vars)), method = "LRT")
  out$data.name <- paste(vars)
  print(out)
}

regtest("~pa") #p=0.47, it works!

#create a character vector with the combinations, categorical variable sonly
regtest.col <- c("~caff_c","~sex", "~race","~pa")

#anova F-test, looping
regtest_result <- list()
regtest_result<- lapply(regtest.col, regtest) #caff:0.931; sex:0.027; race: 0.172; pa: 0.479

##Adjusted model--------------------------------------------------------------------
#select caff_c, sex, race, age, bmi, fiber, fat
glm_adj <- svyglm(galldis~caff_c+sex+race+age+bmi+fiber+fat,subset(nhanes_wt, keep==1),
                  family = binomial())
summary(glm_adj)
summ(glm_adj)

coef_glm_ad1 <- tidy(glm_adj, conf.int = TRUE) %>%
  mutate(OR = exp(estimate), 
         ORL = exp(conf.low),
         ORH = exp(conf.high),
         model = "adjusted") 
write_xlsx(coef_glm_ad1, "adj_model.xlsx")

regtest <- function(vars){
  out <- regTermTest(glm_adj, as.formula(paste0(vars)), method = "LRT")
  out$data.name <- paste(vars)
  print(out)
}

regtest.col <- c("~caff_c","~sex", "~race")
regtest_result <- list()
regtest_result<- lapply(regtest.col, regtest)#caff_c:0.943; sex:0.007; race:0.111

##Assumption checking --------------------------------------------------------
#Collinearity, use unweighted sample
uwt_adj <- glm(galldis~caff_c+sex+race+age+bmi+fiber+fat, data = nhanes_c, family = "binomial")
summary(uwt_adj)
pred_adj <- augment(uwt_adj)

vif(uwt_adj) #looks okay

#Linearity check with test+logit plot, use unweighted sample
nhanes_c_lineary <- nhanes_c %>%
  mutate(agesq = age^2,
         bmisq = bmi^2,
         fibersq = fiber^2,
         fatsq = fat^2)

#Age
uwt_age <- glm(galldis~caff_c+sex+race+age+bmi+fiber+fat+agesq, data = nhanes_c_lineary, family = "binomial")
summary(uwt_age) #age^2: p = 0.072

pred_adj %>% 
  ggplot(., aes(y=.fitted, x=age)) +
  geom_point(colour = "gray40") +
  geom_smooth(method = "lm", formula = y~x, se = FALSE, size = 2, colour = "black") +
  geom_smooth(method = "loess", formula = y~x, se = FALSE, size = 2) + 
  theme_bw() #logit plot looks very nice!

#BMI
uwt_bmi <- glm(galldis~caff_c+sex+race+age+bmi+fiber+fat+bmisq, data = nhanes_c_lineary, family = "binomial")
summary(uwt_bmi) #bmi^2: p = 0.21

pred_adj %>% 
  ggplot(., aes(y=.fitted, x=bmi)) +
  geom_point(colour = "gray40") +
  geom_smooth(method = "lm", formula = y~x, se = FALSE, size = 2, colour = "black") +
  geom_smooth(method = "loess", formula = y~x, se = FALSE, size = 2) + 
  theme_bw() #looks great!

#Fiber
uwt_fiber <- glm(galldis~caff_c+sex+race+age+bmi+fiber+fat+fibersq, data = nhanes_c_lineary, family = "binomial")
summary(uwt_fiber)#fiber^2: p = 0.16

pred_adj %>% 
  ggplot(., aes(y=.fitted, x=fiber)) +
  geom_point(colour = "gray40") +
  geom_smooth(method = "lm", formula = y~x, se = FALSE, size = 2, colour = "black") +
  geom_smooth(method = "loess", formula = y~x, se = FALSE, size = 2) + 
  theme_bw() #one point skewed the plot a bit, but acceptable!

#Fat
uwt_fat <- glm(galldis~caff_c+sex+race+age+bmi+fiber+fat+fatsq, data = nhanes_c_lineary, family = "binomial")
summary(uwt_fat) #bmi^2: p = 0.82

pred_adj %>% 
  ggplot(., aes(y=.fitted, x=fat)) +
  geom_point(colour = "gray40") +
  geom_smooth(method = "lm", formula = y~x, se = FALSE, size = 2, colour = "black") +
  geom_smooth(method = "loess", formula = y~x, se = FALSE, size = 2) + 
  theme_bw() #looks great!

#No linearity violation!! Thank goodness....

#Influential point
#check with weighted model
#Leverage
resid.adj <- augment(glm_adj)
resid.adj$h <- hatvalues(glm_adj)

#Leverage plot
Lev_adj <- ggplot(data = resid.adj, aes(x = .std.resid, y = h))+
  geom_point()+
  geom_vline(xintercept = c(-4, 4), linetype = "dashed", col = "#F8766D")+
  #geom_hline(yintercept = 0.016, linetype = "dashed", col = "#F8766D")+
  theme_bw()+
  labs(title = "Leverage vs. standardized residual",
       x = "Standardized residual",
       y = "Leverage") +
  theme(plot.title = element_text(size = 14, face = "bold"))

Lev_adj

#Pregibon DBeta
dbeta_df<- dx(glm_adj, byCov = F)
dbeta_df <- dbeta_df %>%
  select(dDev) %>%
  cbind(nhanes_c) %>%
  arrange(SEQN)

#plot dDev against subject ID
dbeta_adj <- ggplot(data = dbeta_df, aes(x = SEQN, y = dDev))+
  geom_point()+
  #geom_vline(xintercept = c(-4, 4), linetype = "dashed", col = "#F8766D")+
  #geom_hline(yintercept = 0.016, linetype = "dashed", col = "#F8766D")+
  theme_bw()+
  labs(title = "dDev vs. ID",
       x = "Subject IDl",
       y = "dDev") +
  theme(plot.title = element_text(size = 14, face = "bold"))

dbeta_adj #7 data points need to be excluded

#create new var: outlier (Yes/No)
dbeta_df <- dbeta_df %>%
  mutate(outlier = factor(ifelse(dDev>20, "Yes", "No")))

which(grepl("Yes", dbeta_df$outlier))#identify the row no.

nhanes_c <- nhanes_c %>%
  slice(-c(3769,3776,3780,3783,3784,3785,3793)) #then rerun the model

##Unadjusted model----------------------------------------------------------------
glm_uni_c <- svyglm(galldis~caff_c, subset(nhanes_wt, keep==1), family = quasibinomial())
summary(glm_uni_c)
summ(glm_uni_c)

coef_glm_uni_c <- tidy(glm_uni_c, conf.int = TRUE) %>%
  mutate(OR = exp(estimate), 
         ORL = exp(conf.low),
         ORH = exp(conf.high),
         model = "unadjusted") 
write_xlsx(coef_glm_uni_c, "uni_model.xlsx")

##Figure 2. show OR distribution (forest plot)------------------------------------------------
age5_adj <- c("age5", 1.0368343^5, 1.0437951^5, 1.0299200^5, "adjusted")
fiber_adj <- c("fiber5", 0.981034963^5, 0.964411001^5, 0.997945479^5, "adjusted")
fat_adj <- c("fat5", 0.9982671^5, 1.0013947^5, 0.9951492^5, "adjusted")

plotdf <- coef_glm_uni_c %>% 
  rbind(coef_glm_ad1) %>%
  select(term, OR, ORH, ORL, model) %>%   
  rbind(age5_adj, fiber_adj, fat_adj) %>%
  slice(-c(1,5,14, 16, 17)) %>%
  mutate(OR = as.numeric(OR),
         ORH = as.numeric(ORH),
         ORL = as.numeric(ORL)) 

plotdf$term[plotdf$term=="fat5"] <- "Fat (5g)"
plotdf$term[plotdf$term=="fiber5"] <- "Fiber (5g)"
plotdf$term[plotdf$term=="bmi"] <- "BMI"
plotdf$term[plotdf$term=="age5"] <- "Age (5y)"
plotdf$term[plotdf$term=="raceother"] <- "Other races"
plotdf$term[plotdf$term=="raceasian"] <- "Asian"
plotdf$term[plotdf$term=="racehispanic"] <- "Hispanic"
plotdf$term[plotdf$term=="raceblack"] <- "Black"
plotdf$term[plotdf$term=="sexfemale"] <- "Female"
plotdf$term[plotdf$term=="caff_cvery high"] <- "Very high caffeine"
plotdf$term[plotdf$term=="caff_chigh"] <- "High caffeine"
plotdf$term[plotdf$term=="caff_clow"] <- "Low caffeine"


  
f <- factor(c("Fat (5g)","Fiber (5g)","BMI", "Age (5y)","Other races","Asian",
              "Hispanic","Black", "Female","Very high caffeine", "High caffeine",
              "Low caffeine"))

adj = .3
colors <- c("Unadjusted model  p = 0.629\nPseudo R^2 = 0.001"="gray50", 
            "Adjusted model  p < 0.001\nPseudo R^2 = 0.101"="gray50")

II.Plot.OR <- plotdf %>%
  ggplot(., aes(x = OR, y = term)) +
  geom_vline(aes(xintercept = 1), size = .3, linetype = "dashed") +
  geom_errorbarh(data = filter(plotdf, model == "unadjusted"), aes(xmax = ORH, xmin = ORL), 
                 size = .4, height = .2, color = "gray50", position = position_nudge(y = adj)) +
  geom_point(data = filter(plotdf, model == "unadjusted"), aes(color = "Unadjusted model  p = 0.629\nPseudo R^2 = 0.001"), 
             shape = 21, size = 2, fill = "red", position = position_nudge(y = adj)) +
  geom_errorbarh(data = filter(plotdf, model == "adjusted"), aes(xmax = ORH, xmin = ORL), 
                 size = .4, height = .2, color = "gray50") +
  geom_point(data = filter(plotdf, model == "adjusted"), aes(color = "Adjusted model  p < 0.001\nPseudo R^2 = 0.101"),
             shape = 21, size = 2, fill = "pink")+
  theme_apa(remove.y.gridlines = F)+
  scale_y_discrete(limits = fct_inorder(f))+
  scale_x_continuous(breaks = seq(0,6,1) ) +
  coord_trans(x = "log10")+
  labs(x ="Odds ratio (log-scale)", y = "Predictors", color = "Legend") +
  scale_colour_manual(values=colors) + 
  guides(colour = guide_legend(override.aes = list(fill = c("pink","red"))),
         fill = guide_legend(ncol=2, nrow=2, byrow = T))+
  theme(#panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill="white",
                                          size=0.2, 
                                          linetype = "solid",
                                          colour = "black"),
        legend.position = "bottom",
        legend.box="vertical",
        #legend.justification = c("bottom"),
        #legend.box.just = "right",
        legend.margin = margin(0, 3, 0, -10)
  )

II.Plot.OR


##Model with interaction term--------------------------------------------------------
#Does this relationship differ between diff race/ethinicity?
glm_int_race <- svyglm(galldis~caff_c+sex+race+age+bmi+fiber+fat+caff_c*race,subset(nhanes_wt, keep==1),
                  family = quasibinomial())
summary(glm_int_race, df = degf(nhanes_wt)) 
regTermTest(glm_int_race, ~caff_c*race, method = "LRT", df = degf(nhanes_wt)) #p = 0.087
#oevrall caff_c*race isn't significant

#Does this relationship differ between sex, due to hormone?
glm_int_sex <- svyglm(galldis~caff_c+sex+race+age+bmi+fiber+fat+caff_c*sex,subset(nhanes_wt, keep==1),
                      family = quasibinomial())
summary(glm_int_sex, df = degf(nhanes_wt))

coef_glm_int <- tidy(glm_int_sex, conf.int = T) %>%
  mutate(OR = exp(estimate),
         ORL = exp(conf.low),
         ORH = exp(conf.high))

write_xlsx(coef_glm_int, "coef_int.xlsx")

regTermTest(glm_int_sex, ~caff_c*sex, method = "LRT", df = degf(nhanes_wt))#p = 0.014
#overall caff_c*sex is significant
regTermTest(glm_int_sex, ~caff_c, df = degf(nhanes_wt))


##Compare model fit--------------------------------------------------------
regTermTest(glm_uni_c, ~caff_c, method = "LRT") #0.650
regTermTest(glm_adj, ~caff_c+sex+race+age+bmi+fiber+fat, method = "LRT", df = degf(nhanes_wt)) #p<0.001
regTermTest(glm_adj, ~sex+race+age+bmi+fiber+fat, method = "LRT", df = degf(nhanes_wt)) #p<0.001
regTermTest(glm_int_sex, ~caff_c+sex+race+age+bmi+fiber+fat+caff_c*sex, method = "LRT", df = degf(nhanes_wt)) #p<0.001
regTermTest(glm_int_sex, ~caff_c*sex, method = "LRT", df = degf(nhanes_wt)) #p = 0.015

AIC(glm_uni_c, glm_adj)
BIC(glm_uni_c, glm_adj, glm_int_sex, maximal = glm_int_sex)

psrsq(glm_uni_c, type="Nagelkerke")
psrsq(glm_adj, type="Nagelkerke")
psrsq(glm_int_sex, type="Nagelkerke") #model with the interaction term is the best one!

##Stratified analysis by sex--------------------------------------------------------
#create a new data var "sf", those with [keep==1, male]==1, [keep==1,female]==2
nhanes_final_sf <- nhanes_final %>%
  mutate(sf = factor(ifelse(keep==1 & sex=="male", 1,
                     ifelse(keep==1 & sex=="female",2,
                     ifelse(keep==0, 0, NA)))))
summary(nhanes_final_sf)

#complex survey design
nhanes_sf_wt <- svydesign(id = ~psu, weights = ~weighted, strata = ~stra, data = nhanes_final_sf, nest = T)
options(survey.lonely.psu = "adjust")

#adjusted model for men
glm_sf_male <- svyglm(galldis~caff_c+race+age+bmi+fiber+fat, subset(nhanes_sf_wt, sf==1), 
                      family = quasibinomial)
summary(glm_sf_male)

coef_sf_male <- tidy(glm_sf_male, conf.int = T) %>%
  mutate(OR = exp(estimate),
         ORL = exp(conf.low),
         ORH = exp(conf.high),
         model = "Men")

regTermTest(glm_sf_male, ~caff_c, method = "LRT") #0.758

#adjusted model for women
glm_sf_female <- svyglm(galldis~caff_c+race+age+bmi+fiber+fat, subset(nhanes_sf_wt, sf==2), 
                      family = quasibinomial)

summary(glm_sf_female)

coef_sf_female <- tidy(glm_sf_female, conf.int = T) %>%
  mutate(OR = exp(estimate),
         ORL = exp(conf.low),
         ORH = exp(conf.high),
         model = "Women")

regTermTest(glm_sf_female, ~caff_c, method = "LRT") #0.972
##Figure 3. interaction --------------------------------------------------------
#X = 4 different caffiene intake groups, normal group as reference
vec_m <- c("caff_normal", 1,1,1, "Men")
vec_f <- c("caff_normal", 1,1,1, "Women")

sf_plot <- coef_sf_female %>%
  select(term, OR, ORL, ORH, model) %>%
  rbind(coef_sf_male %>% select(term, OR, ORL, ORH, model)) %>%
  slice(2:4, 14:16) %>%
  rbind(vec_m, vec_f) %>%
  mutate(OR = as.numeric(OR),
         ORL = as.numeric(ORL),
         ORH = as.numeric(ORH))
summary(sf_plot)

sf_plot$term[sf_plot$term=="caff_clow"] <- "Low"
sf_plot$term[sf_plot$term=="caff_normal"] <- "Normal"
sf_plot$term[sf_plot$term=="caff_chigh"] <- "High"
sf_plot$term[sf_plot$term=="caff_cvery high"] <- "Very High"

#plotting
caf <- factor(c("Low", "Normal","High","Very High"))
adj = .1
colors_sf <- c("Men" = "gray50",
               "Women" = "gray50")

III.Plot.OR <- sf_plot %>%
  ggplot(., aes(x = term, y = OR)) +
  geom_hline(aes(yintercept = 1), size = .2, linetype = "dashed") +
  geom_errorbar(data = filter(sf_plot, model == "Men"), aes(ymax = ORH, ymin = ORL), 
                 size = .4, width = 0.2, color = "gray50", position = position_nudge(x = adj)) +
  geom_point(data = filter(sf_plot, model == "Men"), aes(color = "Men"), 
            shape = 21, size = 2, fill = "#00A5FF", position = position_nudge(x = adj)) +
  geom_errorbar(data = filter(sf_plot, model == "Women"), aes(ymax = ORH, ymin = ORL), 
                 size = .4, width = 0.2, color = "gray50") +
  geom_point(data = filter(sf_plot, model == "Women"), aes(color = "Women"),
             shape = 21, size = 2, fill = "#FC717F")+
  theme_apa(remove.x.gridlines = F) +
  scale_x_discrete(limits = fct_inorder(caf))+
  scale_colour_manual(values=colors_sf) + 
  guides(colour = guide_legend(override.aes = list(fill = c("#00A5FF","#FC717F"))),
         fill = guide_legend(ncol=2, nrow=2, byrow = T))+
  labs(x = "Caffeine Intake Categories", y = "Odds ratio") +
  theme(#panel.grid.minor = element_blank(),
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold", vjust=-0.5),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 10),
    
  ) 

III.Plot.OR







