#### Part 1: Preparation of the data set
#### Wan Chi Chang
#### wanchichang04@gmail.com
#### Start Date: 2020-10-19
#### Last edit: 2020-11-30, finish creating the final dataset and save as csv

##Set up workspace--------------------------------------------

#clear workspace
rm(list=ls())
cat("\014")
dev.off()

#set working directory
setwd("~/Desktop/2020 Fall Tufts/Int bio/Final project/Output")

#install packages
library(tidyverse)
library(nhanesA) #To download nhanes data set in r
library(broom) #Save model regression information
library(car) #CRP plot
library(DescTools) #post-hoc test

##Data import --------------------------------------------------
nhanesTables("DEMO", 2017)
nhanesTables("DIET", 2017)
nhanesTables("EXAM", 2017)
nhanesTables("LAB", 2017)
nhanesTables("Q",2017)

demo <- nhanes("DEMO_J")
dietary <- nhanes("DR1TOT_J")
sup <- nhanes("DSQTOT_J")
exam <- nhanes("BMX_J")
lab <- nhanes("TCHOL_J")
medical <- nhanes("MCQ_J")
PA <- nhanes("PAQ_J")
DM <- nhanes("DIQ_J")
preg <- nhanes("RHQ_J")
cig <- nhanes("SMQ_J")
names(exam)

##Variable selection and create new data set-----------------------
demo_v <- demo %>%
  select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH3, DMDEDUC2, INDFMIN2, SDMVPSU, SDMVSTRA, WTMEC2YR) %>%
  mutate(SEQN = factor(SEQN))
summary(demo_v)

dietary_v <- dietary %>%
  select(SEQN, WTDRD1, DR1TCAFF, DR1TFIBE, DR1TTFAT, DR1TCHOL, DR1TVC, `DR1.320Z`, DR1TALCO)%>%
  mutate(SEQN = factor(SEQN)) 
summary(dietary_v)

#Not going to include dietary supplementation
sup_v <- sup %>%
  select(SEQN, DSQTCAFF, DSQTFIBE, DSQTVC)%>%
  mutate(SEQN = factor(SEQN)) #too many missing data
summary(sup_v)

exam_v <- exam %>%
  select(SEQN, BMXWT, BMXBMI, BMXWAIST)%>%
  mutate(SEQN = factor(SEQN))
summary(exam_v)
table(exam_v$BMI_c)

#Don't do lab?
lab_v <- lab %>%
  select(SEQN, LBXTC)%>%
  mutate(SEQN = factor(SEQN))
summary(lab_v)

medical_v <- medical %>%
  select(SEQN, MCQ550, MCQ560)%>%
  mutate(SEQN = factor(SEQN))
summary(medical_v)

PA_v <- PA %>%
  select(SEQN, PAQ605, PAD615, PAQ620, PAD630, PAQ650, PAD660, PAQ665, PAD675) %>%
  mutate(SEQN = factor(SEQN)) #too many missing data
summary(PA_v)

DM_v <- DM %>%
  select(SEQN, DIQ010)%>%
  mutate(SEQN = factor(SEQN))
summary(DM_v)

preg_v <- preg %>%
  select(SEQN, RHD143)%>%
  mutate(SEQN = factor(SEQN))

cig_v <- cig %>%
  select(SEQN, SMQ020, SMQ040)%>%
  mutate(SEQN = factor(SEQN))
summary(cig_v)

##Data set create--------------------------------------------
nhanes_re <- list(dietary_v, demo_v, medical_v, DM_v, PA_v, cig_v, exam_v, lab_v, preg_v)%>% #demo, caffiene intake and gallstone
  reduce(left_join, by= "SEQN") 

summary(nhanes_re)
write_csv(nhanes_re, "nhanes_processed.csv")
