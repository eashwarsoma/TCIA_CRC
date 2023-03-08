library(tidyverse)
library(survival)
library(survminer)
library(reshape)
library(stringr)
library(dplyr)
library(ggplot2)
library(readxl)

#Function to linearly scale something to 0 to 1
normalize <- function(col) {
  return ((col - min(col)) / (max(col) - min(col)))
}

#Basic Cox Model Function
univ.Cox <- function(stat, int, vars, dat) {
  #List object to collect data
  list.UR <- list()
  #For loop to run through all the vars
  for (i in 1:length(vars)) {
    #Creating the relevant formula
    form <- as.formula(paste("Surv(", int, ",", stat, ")", "~", vars[i], sep = ""))
    
    #Computing the Cox Model
    res.cox <- coxph(form, data = dat)
    
    # CI
    CI <- exp(confint(res.cox))
    
    #Original HR
    HR <- exp(coef(res.cox))
    
    #Pvalue
    pval <- summary(res.cox)$coefficients[,5]
    
    #Collecting the data and sticking it into the list
    list.UR[[i]] <- cbind(HR, CI, pval)
    
  }
  
  #With the list object, turn it into a data frame
  df.UR <- do.call(rbind, list.UR)
  
  return(df.UR)
  
}

#Multivariable Cox Model Function
multv.Cox <- function(stat, int, vars, dat) {
  #Creating the relevant formula
  form <- as.formula(paste(paste("Surv(", int, ",", stat, ")"), 
                           paste(vars, collapse=" + "), sep=" ~ "))
  
  #Computing the Cox Model
  res.cox <- coxph(form, data = dat)
  
  #lower CI
  CI <- exp(confint(res.cox))
  
  
  
  #Original HR
  HR <- exp(coef(res.cox))
  
  #Pvalue
  pval <- summary(res.cox)$coefficients[,5]
  
  df.mr <- cbind(HR, CI, pval)
  
  return(df.mr)
}

#Reading in excel file
#Sheet 1 for data
clinic <- read_excel("clinical_data.xlsx")

#Good colnames 
names(clinic) <- c("PID", "scout_name", "age", "sex", "is_comorbidity", 
                   "is_prior_chemo", "risk_score", 
                   "risk_cat", "is_extrahep_disease", 
                   "is_steatosis", 
                   "is_sinusoidal_dilatation", 
                   "NASH_score", "is_NASH", 
                   "is_NASH_4_plus", "percent_response", 
                   "percent_necrosis", "percent_fibrosis", 
                   "is_fibrosis_40_percent_plus", 
                   "percent_mucin", "OS_months", "vital_status", 
                   "is_recurrence", 
                   "DFS_months", "vital_status_DFS", 
                   "is_recurrence_in_liver", 
                   "liver_DFS_months", "vital_status_liver_DFS", 
                   "relevant_notes")

####Transforming all the variables####
#Sex
clinic$sex <- case_when(
  clinic$sex == 1 ~ "male",
  clinic$sex == 2 ~ "female"
)

#Comorbidity Presence?
clinic$is_comorbidity <- case_when(
  clinic$is_comorbidity == 1 ~ "yes",
  clinic$is_comorbidity == 0 ~ "no"
)

#Prior Chemo?
clinic$is_prior_chemo <- case_when(
  clinic$is_prior_chemo == 1 ~ "yes",
  clinic$is_prior_chemo == 0 ~ "no"
)

#Clinic Risk Score #-999 values aren't applicable
clinic$risk_score <- case_when(
  clinic$risk_score == -999 ~ NA_real_,
  TRUE ~ clinic$risk_score
)

#Stratified clinic risk score
clinic$risk_cat <- case_when(
  clinic$risk_cat == 1 ~ "high",
  clinic$risk_cat == 0 ~ "low",
  clinic$risk_cat == -999 ~ NA_character_
)


#Extrahepatic Disease
clinic$is_extrahep_disease <- case_when(
  clinic$is_extrahep_disease == 1 ~ "yes",
  clinic$is_extrahep_disease == 0 ~ "no"
)

#Presence of steatosis
clinic$is_steatosis <- case_when(
  clinic$is_steatosis == 1 ~ "yes",
  clinic$is_steatosis == 0 ~ "no"
)

#Presence of dilated sinousoids
clinic$is_sinusoidal_dilatation <- case_when(
  clinic$is_sinusoidal_dilatation == 1 ~ "yes",
  clinic$is_sinusoidal_dilatation == 0 ~ "no"
)

#NASH Score
#No need to change
clinic$NASH_score

#Presence of NASH
clinic$is_NASH <- case_when(
  clinic$is_NASH == 1 ~ "yes",
  clinic$is_NASH == 0 ~ "no"
)

#Is Nash score at least 4
clinic$is_NASH_4_plus <- case_when(
  clinic$is_NASH_4_plus == 1 ~ "yes",
  clinic$is_NASH_4_plus == 0 ~ "no"
)

#Percent of total response
clinic$percent_response <- case_when(
  clinic$percent_response == -999 ~ NA_real_,
  TRUE ~ clinic$percent_response
)

#Percent of tumor necrosis
clinic$percent_necrosis <- case_when(
  clinic$percent_necrosis == -999 ~ NA_real_,
  TRUE ~ clinic$percent_necrosis
)

#Percent of Fibrosis
clinic$percent_fibrosis <- case_when(
  clinic$percent_fibrosis == -999 ~ NA_real_,
  TRUE ~ clinic$percent_fibrosis
)

#Is fibrosis greater than 40 percent
clinic$is_fibrosis_40_percent_plus <- case_when(
  clinic$is_fibrosis_40_percent_plus == 1 ~ "high",
  clinic$is_fibrosis_40_percent_plus == 0 ~ "low",
  clinic$is_fibrosis_40_percent_plus == -999 ~ NA_character_
)

#Is mucin present
clinic$percent_mucin <- case_when(
  clinic$percent_mucin == -999 ~ NA_real_,
  TRUE ~ clinic$percent_mucin
)

#Keeping all the survival and regression metric
#1 indicates event happened (death or progression)
#0 is no event (no recurrence or patient alive)
#Overall survival
clinic$OS_months
#Vital Status 
clinic$vital_status
#Did progression occur
clinic$is_recurrence
#Disease free survival in months
clinic$DFS_months
#Vital status of DFS 
#A 1 means patient either died or had recurrence
clinic$vital_status_DFS
#Progression in liver
clinic$is_recurrence_in_liver
#Months to liver progression
clinic$liver_DFS_months
#Patient status for liver DFS 
#A 1 means patient either died or had liver recurrence
clinic$vital_status_liver_DFS


#####Logistic Model on Regression####













