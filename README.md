# Dissertation-Rene
##-------------------------DISSERTATION 2210245---------------------------------
## TOPIC: Predicting cardiovascular risk in patients with chronic obstructive 
##        pulmonary disease using routine health data in Wales

# Population:Patients with Chronic Obstructive pulmonary disease in Wales

# Intervention: The use of routine health data for cardiovascular risk prediction

# Outcome: Ability to predict cardiovascular risk in patients with COPD using
#          routine health data available in Wales.
--------------------------------------------------------------------------------
# Install and run the packages below.
# install.packages("dplyr")
# install.packages("tidyverse")
# install.packages("caret")
install.packages('e1071')
# install.packages('yardstick')
# install.packages('neuralnet')
# install.packages('xgboost')
# install.packages('gbm')
# install.packages('ggplot2')
install.packages('MLeval')


library(dplyr)
library(tidyverse)
library(lubridate)
library(tcltk)
library(RODBC)
library(ROCR)
library(caret)
library(pROC)
library(caTools)
library(e1071)
library(yardstick)
library(broom)
library(rpart)
library(neuralnet)
library(xgboost)
library(gbm)
library(ggplot2)
--------------------------------------------------------------------------------
## Set working directory and to assess SAIL Databank  use the query below

# setwd("P:/abangr/main-repository-for-1281/renesail/Dissertation")  

source("login.R")

#------------------------------------------------------------------------------
# A. Run Query to import COPD patients and demographics
#    from SAIL1281V project Databank using PEDW_SINGLE_DIAG_20220228 and 
#    PEDW_SINGLE_DIAG_20220228 .
#-----------------------------------------------------------------------------
# Import COPD patients using Read CODES from SAILW1281V

copd <- sqlQuery(channel, "SELECT DISTINCT ALF_PE, min(EVENT_DT) AS DIAG_DT, WOB,

                GNDR_CD AS gender FROM SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 
                
                WHERE EVENT_CD IN (SELECT READ_CODE FROM 
                 
                SAILW1281V.COMORBID_READ_CODES_1281 WHERE CD_TYPE_CATEGORY 
                 
                LIKE 'COPD%') AND EVENT_DT >'2000-12-31'AND ALF_PE IS NOT NULL GROUP BY ALF_PE, 
                 
                 EVENT_DT, WOB, GNDR_CD")

# Calculate age in years and attach copd variable and select age greater than 5

copd <- copd %>% mutate(age = as.numeric(difftime(DIAG_DT, WOB, units = "days"))/
                
                  365.25) %>% mutate(age = round(age)) %>% mutate(copd = 1) %>% 
                  
                 select(-WOB) %>% filter(age >= 30 & age <= 90) %>% arrange(ALF_PE)



# Import CVD patients using Read CODES from SAILW1281V 

cvd <- sqlQuery(channel, "SELECT DISTINCT ALF_PE, min(EVENT_DT) AS D_DT FROM

              SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD IN (
              
              SELECT READ_CODE FROM SAILW1281V.COMORBID_READ_CODES_1281 WHERE 
              
              CD_TYPE_CATEGORY LIKE 'Vascular%' OR CD_TYPE_CATEGORY LIKE 
              
              'Heart failure%') GROUP BY ALF_PE")

# Create CVD variable and attach 1 and remove DIAG_DT

cvd <- cvd %>% mutate(cvd = 1) %>% select(-D_DT) %>% arrange(ALF_PE)



#-------------------------------------------------------------------------------
# B. Import Risk factors from SAIL1281V-GP_EVENT_CLEANSED_20220201

diabetes <- sqlQuery(channel,"SELECT DISTINCT ALF_PE FROM 

                 SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD IN (

                 SELECT READ_CODE FROM SAILW1281V.COMORBID_READ_CODES_1281 WHERE
                      
                      CD_TYPE_CATEGORY LIKE 'Diabetes%') GROUP BY ALF_PE")

diabetes <- diabetes %>% mutate(diabetes = 1) %>% arrange(ALF_PE)


hypertension <- sqlQuery(channel,"SELECT DISTINCT ALF_PE FROM 

                 SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD IN (

                 SELECT READ_CODE FROM SAILW1281V.COMORBID_READ_CODES_1281 WHERE
                      
                      CD_TYPE_CATEGORY LIKE 'Hypertension%') GROUP BY ALF_PE")

hypertension <- hypertension %>% mutate(hypertension = 1) %>% arrange(ALF_PE)

mi <- sqlQuery(channel, "SELECT DISTINCT ALF_PE FROM 

                  SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD IN (
              
                  SELECT READ_CODE FROM SAILW1281V.COMORBID_READ_CODES_1281 
                    
                   WHERE CD_TYPE_CATEGORY LIKE 'MI') GROUP BY ALF_PE")

mi <- mi %>% mutate(mi = 1) %>% arrange(ALF_PE)

smoking <- sqlQuery(channel, "SELECT DISTINCT ALF_PE FROM 

                  SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD IN (
              
                  SELECT READ_CODE FROM SAILW1281V.COMORBID_READ_CODES_1281 
                    
                   WHERE CD_TYPE_CATEGORY IN ('Smoker%', 'Ex-smoker')) 
                    
                  GROUP BY ALF_PE")

smoking <- smoking %>% mutate(smoking = 1) %>% arrange(ALF_PE)


depression <- sqlQuery(channel, "SELECT DISTINCT ALF_PE FROM 
                       
                       SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD 
                       
                       IN (SELECT DISTINCT EVENT_CD FROM SAILW1281V.NR_DEPRESSION nd) 
                       
                       GROUP BY ALF_PE")

depression <- depression %>% mutate(depression = 1) %>% arrange(ALF_PE)


anxiety <- sqlQuery(channel, "SELECT DISTINCT ALF_PE FROM 
                       
                       SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD 
                       
                       LIKE 'E200.' GROUP BY ALF_PE")


anxiety <- anxiety %>% mutate(anxiety = 1) %>% arrange(ALF_PE)


alcohol <- sqlQuery(channel, "SELECT DISTINCT ALF_PE FROM 
                       
                       SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD 
                       
                       IN ('E23..', '1282.', '1462.', '136T.', '136S.', '136V.') 
                    
                       GROUP BY ALF_PE")

alcohol <- alcohol %>% mutate(alcohol = 1) %>% arrange(ALF_PE)


fam_hist_hd <- sqlQuery(channel, "SELECT DISTINCT ALF_PE FROM 
                       
                       SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD 
                       
                       LIKE 'ZV173' GROUP BY ALF_PE")

fam_hist_hd <- fam_hist_hd %>% mutate(fam_hist_hd = 1) %>% arrange(ALF_PE)


emphysema <- sqlQuery(channel, "SELECT DISTINCT ALF_PE FROM 
                       
                       SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD 
                       
                       IN ('H32..', 'SP2y0', 'H581.', 'SK07.', 'HS82.', 'H3121',
                      
                      'Q312..', 'J6502', 'H4640') GROUP BY ALF_PE")

emphysema <- emphysema %>% mutate(emphysema = 1) %>% arrange(ALF_PE)


education <- sqlQuery(channel, "SELECT DISTINCT ALF_PE FROM 
                       
                       SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD 
                       
                       IN('13Z4.', '03...', '03A..') GROUP BY ALF_PE")

education <- education %>% mutate(education = 1) %>% arrange(ALF_PE)


sleep_apnoea <- sqlQuery(channel, "SELECT DISTINCT ALF_PE FROM 

                  SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD IN (
              
                  SELECT READ_CD FROM SAILREFRV.READ_CD WHERE 
                         
                  READ_DESC LIKE '%Sleep apnoea%' OR READ_DESC LIKE 
                         
                  '%sleep apnoea%') GROUP BY ALF_PE")

sleep_apnoea <- sleep_apnoea %>% mutate(sleep_apnoea = 1) %>% arrange(ALF_PE)


bmi <- sqlQuery(channel, "SELECT DISTINCT ALF_PE, min(EVENT_DT) AS EVENT_DT, 
                
                  EVENT_VAL AS bmi FROM SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 
                  
                  WHERE EVENT_VAL > '0' AND (EVENT_CD LIKE '22K%' AND EVENT_CD != '22K9.') 
 
                   GROUP BY ALF_PE, EVENT_DT, EVENT_VAL")

if (!is.numeric(bmi$BMI)) {
 
  bmi$BMI <- as.numeric(bmi$BMI)
}

bmi <- bmi %>% filter(BMI >= 10 & BMI <= 100) %>%  select(ALF_PE, BMI) %>% 
  
     dplyr::distinct(ALF_PE, .keep_all = TRUE)



weight <- sqlQuery(channel, "SELECT DISTINCT ALF_PE, min(EVENT_DT) AS EVENT_DT, 
                   
                   EVENT_VAL AS weight FROM 
                   
                   SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec 
                   
                   WHERE EVENT_VAL > '0' AND (EVENT_CD LIKE '22A%' AND EVENT_CD 
                   
                   NOT IN ('22A7.','22A8.','22A9.','22AA.')) GROUP BY ALF_PE, 
                   
                   EVENT_DT, EVENT_VAL")

if (!is.numeric(weight$WEIGHT)) {
  
  weight$WEIGHT <- as.numeric(weight$WEIGHT)
}

weight <- weight %>% select(ALF_PE, WEIGHT) %>% filter(WEIGHT >= 10 & WEIGHT <= 635) %>% 
  
           dplyr::distinct(ALF_PE, .keep_all = TRUE)



ldl <- sqlQuery(channel, "SELECT DISTINCT ALF_PE, min(EVENT_DT) AS EVENT_DT, EVENT_VAL AS ldl 
            
                FROM SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec WHERE 
                
                EVENT_CD = '44P6.' AND EVENT_VAL IS NOT NULL GROUP BY ALF_PE, 
                
                EVENT_DT, EVENT_VAL")

if (!is.numeric(ldl$LDL)) {
  
  ldl$LDL <- as.numeric(ldl$LDL)
}

ldl <- ldl %>% select(ALF_PE, LDL) %>% filter(LDL > 0) %>%  dplyr::distinct(ALF_PE, .keep_all = TRUE)




hdl <- sqlQuery(channel, "SELECT DISTINCT ALF_PE, min(EVENT_DT) AS EVENT_DT, EVENT_VAL AS hdl

               FROM SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec
                
                WHERE EVENT_CD IN ('44PB.', '44d3.') AND EVENT_VAL IS NOT 
                
                NULL GROUP BY ALF_PE, EVENT_DT, EVENT_VAL")

if (!is.numeric(hdl$HDL)) {
  
  hdl$HDL <- as.numeric(hdl$HDL)
}

hdl <- hdl %>% select(ALF_PE, HDL) %>% dplyr::distinct(ALF_PE, .keep_all = TRUE)


diastolic_bp <- sqlQuery(channel, "SELECT DISTINCT ALF_PE, min(EVENT_DT) AS EVENT_DT,

                        EVENT_VAL AS diastolic_bp FROM 
        
                        SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec 
                        
                        WHERE EVENT_CD IN ('246T.', '246L.', '246R.',
                         
                        '246m.', '246P.', '246f.', '246o1') AND EVENT_VAL IS NOT 
                
                         NULL GROUP BY ALF_PE, EVENT_DT, EVENT_VAL")

if (!is.numeric(diastolic_bp$DIASTOLIC_BP)) {
  
  diastolic_bp$DIASTOLIC_BP <- as.numeric(diastolic_bp$DIASTOLIC_BP)
}

diastolic_bp <- diastolic_bp %>% filter(DIASTOLIC_BP <= 370) %>% select(ALF_PE, DIASTOLIC_BP) %>%
  
                dplyr::distinct(ALF_PE, .keep_all = TRUE)


systolic_bp <- sqlQuery(channel, "SELECT DISTINCT ALF_PE, min(EVENT_DT) AS EVENT_DT,

                      EVENT_VAL AS systolic_bp FROM 
                      
                      SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec WHERE 
                      
                      EVENT_CD IN ('246S.', '246K.', '246Q.', '246N.', 
                        
                      '246I.', '246e.', '246o0') AND EVENT_VAL IS NOT NULL
                
                       GROUP BY ALF_PE, EVENT_DT, EVENT_VAL")

systolic_bp <- systolic_bp %>% mutate(SYSTOLIC_BP = as.numeric(SYSTOLIC_BP))

systolic_bp <- systolic_bp %>% filter(SYSTOLIC_BP <= 370) %>% 
  
    select(ALF_PE, SYSTOLIC_BP) %>% dplyr::distinct(ALF_PE, .keep_all = TRUE)


hba1c <- sqlQuery(channel, "SELECT DISTINCT ALF_PE, min(EVENT_DT) AS EVENT_DT, EVENT_VAL

                  AS hba1c FROM SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec
                  
                  WHERE EVENT_CD IN ('44TB.', '44TC.', '42c..', '42W..', 
                  
                  '42W4.', '42W5.') AND EVENT_VAL IS NOT NULL GROUP BY ALF_PE, 
                  
                   EVENT_DT, EVENT_VAL")

hba1c <- hba1c %>% mutate(HBA1C = as.numeric(HBA1C))

hba1c <- hba1c %>% select(ALF_PE, HBA1C) %>% filter(HBA1C > 0) %>% 
     
             dplyr::distinct(ALF_PE, .keep_all = TRUE)


waist_circ  <- sqlQuery(channel, "SELECT DISTINCT ALF_PE, min(EVENT_DT) AS EVENT_DT, 

                       EVENT_VAL AS waist_circ FROM 
                        
                       SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec WHERE 
                        
                       EVENT_CD LIKE '22N0%' AND EVENT_VAL IS NOT NULL GROUP 
                        
                        BY ALF_PE, EVENT_DT, EVENT_VAL")

waist_circ <- waist_circ %>% mutate(WAIST_CIRC = as.numeric(WAIST_CIRC))

waist_circ <- waist_circ %>% select(ALF_PE, WAIST_CIRC) %>% dplyr::distinct(ALF_PE, .keep_all = TRUE)


creatinine <- sqlQuery(channel, "SELECT DISTINCT ALF_PE, min(EVENT_DT) AS EVENT_DT, 
                   
                   EVENT_VAL AS creatinine FROM 
                   
                   SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec 
                   
                   WHERE EVENT_CD IN ('44J3.', '44J30', '44J31', '44J32', '44J33',
                   
                   '44J3z', '46m7') AND EVENT_VAL IS NOT NULL GROUP BY 
                       
                  ALF_PE, EVENT_DT, EVENT_VAL")

creatinine <- creatinine %>% mutate(CREATININE = as.numeric(CREATININE))


creatinine <- creatinine %>% select(ALF_PE, CREATININE) %>% filter(CREATININE > 0) %>% 
  
                           dplyr::distinct(ALF_PE, .keep_all = TRUE)


fev1 <- sqlQuery(channel,"SELECT DISTINCT ALF_PE, min(EVENT_DT) AS EVENT_DT,

                 EVENT_VAL AS fev1 FROM 

                 SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD IN (

                 select READ_CD from SAILREFRV.READ_CD WHERE READ_DESC LIKE
                 
                 '%FEV1%' OR READ_DESC LIKE '%fev1%') GROUP BY ALF_PE, EVENT_DT,
                 
                 EVENT_VAL")

fev1 <- fev1 %>%  mutate(FEV1 = as.numeric(FEV1))

fev1 <- fev1 %>% select(ALF_PE, FEV1) %>% filter(FEV1 > 0) %>% 
  
         dplyr::distinct(ALF_PE, .keep_all = TRUE)


alanine_amt <- sqlQuery(channel, "SELECT DISTINCT ALF_PE, min(EVENT_DT) AS 
                      
                        EVENT_D, EVENT_VAL AS alanine_amt FROM 
                        
                        SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec WHERE 
                        
                        EVENT_CD IN ('44GB.', '44GA.') AND EVENT_VAL >'0'
                        
                        GROUP BY ALF_PE, EVENT_DT, EVENT_VAL")

alanine_amt <- alanine_amt %>% mutate(ALANINE_AMT = as.numeric(ALANINE_AMT))

alanine_amt <- alanine_amt %>% select(ALF_PE, ALANINE_AMT) %>% 
               
               dplyr::distinct(ALF_PE, .keep_all = TRUE)


fibrinogen <- sqlQuery(channel, "SELECT DISTINCT ALF_PE, min(EVENT_DT) AS 
                      
                        EVENT_D, EVENT_VAL AS fibrinogen FROM 
                        
                        SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec WHERE 
                        
                        EVENT_CD LIKE '42Qn.' AND EVENT_VAL >'0'
                        
                        GROUP BY ALF_PE, EVENT_DT, EVENT_VAL")

fibrinogen <- fibrinogen %>% mutate(FIBRINOGEN = as.numeric(FIBRINOGEN))

fibrinogen <- fibrinogen %>% filter(FIBRINOGEN <= 20) %>% select(ALF_PE, FIBRINOGEN) %>% 
  
              dplyr::distinct(ALF_PE, .keep_all = TRUE)


cr_protein <- sqlQuery(channel, "SELECT DISTINCT ALF_PE, min(EVENT_DT) AS 
                      
                        EVENT_D, EVENT_VAL AS cr_protein FROM 
                        
                        SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec WHERE 
                        
                        EVENT_CD IN ('44CC.', '44CS.') AND EVENT_VAL >'0' 
                       
                        GROUP BY ALF_PE, EVENT_DT, EVENT_VAL")

cr_protein <- cr_protein %>% mutate(CR_PROTEIN = as.numeric(CR_PROTEIN))

cr_protein <- cr_protein %>% select(ALF_PE, CR_PROTEIN) %>% dplyr::distinct(ALF_PE, .keep_all = TRUE)

platelet_c <- sqlQuery(channel, "SELECT DISTINCT ALF_PE, min(EVENT_DT) AS 
                      
                        EVENT_D, EVENT_VAL AS platelet_c FROM 
                        
                        SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec WHERE 
                        
                        EVENT_CD LIKE '42P..' AND EVENT_VAL >'0' 
                       
                        GROUP BY ALF_PE, EVENT_DT, EVENT_VAL")

platelet_c <- platelet_c %>% mutate(PLATELET_C = as.numeric(PLATELET_C))

platelet_c <- platelet_c %>% select(ALF_PE, PLATELET_C) %>% dplyr::distinct(ALF_PE, .keep_all = TRUE)


serum_sodium <- sqlQuery(channel, "SELECT DISTINCT ALF_PE, min(EVENT_DT) AS 
                      
                        EVENT_D, EVENT_VAL AS serum_sodium FROM 
                        
                        SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 wgec WHERE 
                        
                        EVENT_CD LIKE '44I5.' AND EVENT_VAL >'0' 
                       
                        GROUP BY ALF_PE, EVENT_DT, EVENT_VAL")

serum_sodium <- serum_sodium %>% mutate(SERUM_SODIUM = as.numeric(SERUM_SODIUM))

serum_sodium <- serum_sodium %>% select(ALF_PE, SERUM_SODIUM) %>% 
  
            dplyr::distinct(ALF_PE, .keep_all = TRUE)


#==============================================================================
# B. Merge the risk factors and join with COPD

copd_table <- copd %>% left_join(alanine_amt, by = "ALF_PE") %>% left_join(alcohol,
                                                                           
        by = "ALF_PE") %>% left_join(anxiety, by = "ALF_PE") %>% left_join(
            
        cr_protein, by = "ALF_PE") %>% left_join(creatinine, by = "ALF_PE") %>% 
  
        left_join( depression, by = "ALF_PE") %>% left_join(  diabetes, by = "ALF_PE") %>% 
        
        left_join(diastolic_bp, by = "ALF_PE") %>% left_join(education, by = "ALF_PE") %>%
        
        left_join(emphysema, by = "ALF_PE") %>% left_join(fam_hist_hd, by = "ALF_PE") %>% 
  
        left_join(fev1, by = "ALF_PE") %>% left_join(fibrinogen, by = "ALF_PE") %>% 
  
        left_join(hba1c, by = "ALF_PE") %>% left_join(hdl, by = "ALF_PE") %>% left_join(
          
        bmi, by = "ALF_PE") %>% left_join(hypertension, by = "ALF_PE") %>% left_join(
            
        ldl, by = "ALF_PE") %>% left_join(mi, by = "ALF_PE") %>% left_join(platelet_c, 
                                                                           
        by = "ALF_PE") %>% left_join(serum_sodium, by = "ALF_PE") %>% left_join(
          
        sleep_apnoea, by = "ALF_PE") %>% left_join(smoking, by = "ALF_PE") %>% 
  
        left_join(systolic_bp, by = "ALF_PE") %>% left_join(waist_circ, by = "ALF_PE") %>%
  
       left_join(weight, by = "ALF_PE") %>% arrange(ALF_PE)

#===============================================================================
# C. Replace NAs in Binary variables to Zero

copd_table <- copd_table %>% mutate(anxiety = coalesce(anxiety, 0)) %>% mutate(
  
      depression = coalesce(depression, 0)) %>% mutate(diabetes = coalesce(diabetes, 0)) %>% 
      
      mutate(education = coalesce(education, 0)) %>% mutate(emphysema = coalesce(emphysema, 0)) %>% 
  
     mutate(fam_hist_hd = coalesce(fam_hist_hd, 0)) %>% mutate(sleep_apnoea = coalesce(sleep_apnoea, 0)) %>% 
  
     mutate(smoking = coalesce(smoking, 0)) %>% mutate(alcohol = coalesce(alcohol, 0)) %>% mutate(hypertension =
                                                                                                    
     coalesce(hypertension, 0)) %>% mutate(mi = coalesce(mi, 0))

# Join CVD to the table

copd_table <- copd_table %>% left_join(cvd, by = 'ALF_PE') 

copd_table<- copd_table %>% mutate(cvd = coalesce(cvd, 0))



str(copd_table)

# Remove variables with mor the 25% NAs

copd_table1 <- copd_table %>% select(-c(SYSTOLIC_BP, WAIST_CIRC, FIBRINOGEN,
                                        
                           HDL, DIASTOLIC_BP, education, ALF_PE, DIAG_DT, copd))              
                   
sum(!is.na(copd_table1))


#Data visualisation before data input
copd_table1$ALANINE_AMT %>% hist()

copd_table1$CR_PROTEIN %>% hist()

copd_table1[,-23] %>% gather(attributes, value, 1:ncol(copd_table1[,-23])) %>% ggplot(
  
  aes(x = value)) + geom_histogram(fill = 'lightblue', color ='black', bins = 25) +
  
  facet_wrap(~attributes, scales = "free")
#==============================================================================
##Input missing values using mean value for variables with missing data

copd_cvd2 <- copd_table1

copd_cvd2$ALANINE_AMT[which(is.na(copd_cvd2$ALANINE_AMT))] = mean(copd_cvd2$ALANINE_AMT, na.rm = TRUE)

copd_cvd2$CR_PROTEIN[which(is.na(copd_cvd2$CR_PROTEIN))] = mean(copd_cvd2$CR_PROTEIN, na.rm = TRUE)

copd_cvd2$CREATININE[which(is.na(copd_cvd2$CREATININE))] = mean(copd_cvd2$CREATININE, na.rm = TRUE)

copd_cvd2$FEV1[which(is.na(copd_cvd2$FEV1))] = mean(copd_cvd2$FEV1, na.rm = TRUE)

copd_cvd2$HBA1C[which(is.na(copd_cvd2$HBA1C))] = mean(copd_cvd2$HBA1C, na.rm = TRUE)

copd_cvd2$BMI[which(is.na(copd_cvd2$BMI))] = mean(copd_cvd2$BMI, na.rm = TRUE)

copd_cvd2$PLATELET_C[which(is.na(copd_cvd2$PLATELET_C))] = mean(copd_cvd2$PLATELET_C, na.rm = TRUE)

copd_cvd2$SERUM_SODIUM[which(is.na(copd_cvd2$SERUM_SODIUM))] = mean(copd_cvd2$SERUM_SODIUM, na.rm = TRUE)

copd_cvd2$WEIGHT[which(is.na(copd_cvd2$WEIGHT))] = mean(copd_cvd2$WEIGHT, na.rm = TRUE)

copd_cvd2$LDL[which(is.na(copd_cvd2$LDL))] = mean(copd_cvd2$LDL, na.rm = TRUE)

summary(copd_cvd2)
copd_cvd2$ALANINE_AMT %>% hist()
copd_cvd2$CR_PROTEIN %>% hist()

copd_cvd <- copd_cvd2 %>% mutate(cvd = ifelse(cvd == 1, 'yes', 'no')) %>% 
  
                          mutate(cvd = as.factor(cvd)) 

copd_cvd$cvd <- factor(copd_cvd$cvd, levels = c('yes', 'no'))
str(copd_cvd)

b_copd <- copd_cvd %>% select(GENDER, alcohol, anxiety, depression, diabetes, emphysema,
                              
                              fam_hist_hd, hypertension, mi, sleep_apnoea, smoking, cvd)

c_copd <- copd_cvd %>% select(age, ALANINE_AMT, CR_PROTEIN, CREATININE, FEV1, HBA1C,
                              
                              BMI, LDL, PLATELET_C, SERUM_SODIUM, WEIGHT)

# Scale the continuous variable.

c_copd <- scale(c_copd, center = TRUE, scale = TRUE)

#cbind the scaled variable

scaled_copd <- cbind(c_copd, b_copd)

#==============================================================================
## Data Visualization of completed dataset

scaled_copd[,-23] %>% gather(attributes, value, 1:ncol(scaled_copd[,-23])) %>% ggplot(
  
  aes(x = value)) + geom_histogram(fill = 'lightblue', color ='black', bins = 25) +
  
  facet_wrap(~attributes, scales = "free")

#===============================================================================
# Split Data into 70/30 ratio to train the model

set.seed(250)

tr <- createDataPartition(scaled_copd[,'cvd'], p=0.7, list = FALSE)

train_cc <- scaled_copd[tr,]
test_cc <- scaled_copd[-tr,]
#===============================================================================
# Logistics regression with LASSO

# install.packages('glmnet')
library(glmnet)



control <- trainControl(savePredictions = TRUE, classProbs = TRUE, method = 'cv',
                        
                        number = 20, verboseIter = TRUE)

model <- train(cvd~., data = train_cc, method = 'glmnet', 
               
               tuneGrid = expand.grid(alpha = 1, lambda = seq(0.01, 100, length = 30)),
               
               trControl = control, metric = 'Accuracy')

plot(model)

pred <- predict(model, test_cc)
result <- confusionMatrix(pred, test_cc$cvd)
model_auc <- model$pred %>% roc_auc(obs, yes) %>% select(.estimate) %>%  round(2)

#==============================================================================
# Logistics regression for elastic net
control1 <- trainControl(savePredictions = TRUE, classProbs = TRUE, method = 'cv',
                        
                        number = 2, verboseIter = TRUE)

model1 <- train(cvd~., data = train_cc, method = 'glmnet', 
               
               tuneGrid = expand.grid(alpha =seq(0.01, 0.1) , lambda = seq(0.01, 100, length = 30)),
               
               trControl = control1, metric = 'Accuracy')

plot(model1)
pred1 <- predict(model1, test_cc)
result1 <- confusionMatrix(pred1, test_cc$cvd)
model_auc1 <- model1$pred %>% roc_auc(obs, yes) %>% select(.estimate) %>%  round(2)

#===============================================================================
# Logistics regression for ridge
control2 <- trainControl(savePredictions = TRUE, classProbs = TRUE, method = 'cv',
                         
                         number = 5, verboseIter = TRUE)

model2 <- train(cvd~., data = train_cc, method = 'glmnet', 
                
                tuneGrid = expand.grid(alpha = 0, lambda = seq(0.01, 100, length = 30)),
                
                trControl = control2, metric = 'Accuracy')

plot(model2)
pred2 <- predict(model2, test_cc)
result2 <- confusionMatrix(pred2, test_cc$cvd)
model_auc2 <- model2$pred %>% roc_auc(obs, yes) %>% select(.estimate) %>%  round(2)
plot(varImp(model2))
#==============================================================================
# AUC comparison and selection

cat('The AUC of LASSO regression is', model_auc$.estimate,
    '\nThe AUC of Elastic net regression is',  model_auc1$.estimate,
    '\nThe AUC of ridge regression is',  model_auc2$.estimate,
      sep = ' ')

print('Therefore ridge regression was selected')
#==============================================================================
# Select important variables in ridge regression
plot(varImp(model2))

#store the estimates in a vector

model_d <- coef(model$finalModel, model$bestTune$lambda)

model_d <- model_d[, 1] %>% data.frame()

model_d$variables <- rownames(model_d)

model_d <- model_d %>% filter(. !=0) %>% select(variables) %>% unlist(., use.names = FALSE)

# Remove intercept from the list
model_d <- model_d[-1]

# Add outcome variable
model_d[length(model_d) + 1] <- colnames(train_cc)[ncol(train_cc)]
#==============================================================================

#==============================================================================
## Build Model in glm for analysis

copd_cvd1 <- scaled_copd %>% select(GENDER, age, diabetes, hypertension, CREATININE,  
                                    
                                     sleep_apnoea, fam_hist_hd, LDL, emphysema, alcohol,
                                    
                                    depression, BMI, smoking, WEIGHT, HBA1C, PLATELET_C,
                                    
                                    ALANINE_AMT, FEV1, SERUM_SODIUM, anxiety, CR_PROTEIN, cvd)

str(copd_cvd)
write.csv(copd_cvd1, 'copd_cvd1.csv')

model_glm <- glm(cvd~GENDER + age + diabetes + hypertension + CREATININE + anxiety +
                 
                 sleep_apnoea + fam_hist_hd + LDL + emphysema + alcohol +
                 
                 depression + BMI + smoking + WEIGHT + HBA1C + PLATELET_C +
                 
                 ALANINE_AMT + FEV1 + SERUM_SODIUM + CR_PROTEIN, data = copd_cvd1,  
                 
                 family = "binomial")

summary(model_glm)

coef1 <- model_glm %>% tidy(conf.int = TRUE, exp = TRUE)

view(coef1)
model_glm$


# write.csv(coef1, 'statisticaly significant.csv')
# 
# ## Plot coefficient interval

library(readr)
coef2 <- read.csv("coef.csv")


ggplot(coef2, aes(x = Odds_ratio, y = variables, ymin = LowerCI, ymax = UpperCI)) +

  geom_vline(aes(xintercept = 1), size = 0.25, linetype = 'dashed') +

  geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI), linewidth = 0.5, height = 0.2, color = 'grey50') +

  geom_point(size = 2.5, color = 'sky blue') +

  scale_x_continuous(breaks = seq(0.1, 4.0, 0.1)) + 
  
  theme_bw() + 
  
  ylab('') + theme(panel.grid.minor = element_blank()) + coord_cartesian() +
  
  labs(title = 'Odds ratio of statistically significant risk factors for CVD')


install.packages('forcats')
library(forcats)
view(coef2)
#===============================================================================
# Build model to be used 
# Data Partition
set.seed(20)
df <- createDataPartition(copd_cvd1[,'cvd'], p=0.7, list = FALSE)

train <- copd_cvd1[df,]
test <- copd_cvd1[-df,]
set.seed(50)

# Create grid for hyper parameters

control3 <- trainControl(savePredictions = TRUE, classProbs = TRUE, method = 'cv',
                         
                         number = 5, verboseIter = TRUE)

model3 <- train(cvd~., data = train, method = 'glm', family = 'binomial', 
                
                tuneLength = 20,
                
                trControl = control3)

pr <- predict(model3, test)
glmresult <- confusionMatrix(pr, test$cvd)
glm_auc <- model3$pred %>% roc_auc(obs, yes) %>% select(.estimate) %>%  round(2)


#===============================================================================
# Train order important models and compare with AUC

#Random forest

rf <- train(cvd~., data = train, method = 'rf', ntree = 100, trControl = control3
            )

pr1 <- predict(rf, test)
rfresult <- confusionMatrix(pr1, test$cvd)
rf_auc <- rf$pred %>% roc_auc(obs, yes) %>% select(.estimate) %>% round(2)

#===============================================================================
# Support Vector Machine

# controls <- trainControl(method = 'cv', number = 5, savePredictions = TRUE, classProbs = TRUE)
# 
# svm1 <- train(cvd~., data = train, method = "svmRadial", tuneLength = 5,
#             
#             trControl = controls)
# 
# pr2 <- predict(svm1, test)
# svmresult <-confusionMatrix(pr2, test$cvd)
# svm_auc <- lm$pred %>% roc_auc(obs, yes) %>% select(.estimate) %>% round(2)

#===============================================================================
# Neural Network
 nn <- train(cvd~., data = train, method = 'nnet',
             
             trControl = control3, metrics = 'Accuracy')

pr3 <- predict(nn, test)
nnresult<- confusionMatrix(pr3, test$cvd)
nn_auc <- nn$pred %>% roc_auc(obs, yes) %>% select(.estimate) %>% round(2)

#===============================================================================
## Xgboost

xg <- train(cvd~., data = train, method = 'xgbTree', trControl = trainControl(method = 'cv', 
                                                                    
                                  savePredictions = TRUE, classProbs = TRUE))

pr4 <- predict(xg, test)
xgresult <- confusionMatrix(pr4, test$cvd)
xg_auc <- xg$pred %>% roc_auc(obs, yes) %>% select(.estimate) %>% round(2)

#===============================================================================
# Gradient Boosting

gbm_m <- train(cvd~., data = train, method = 'gbm', trControl = trainControl(method = 'cv', summaryFunction = twoClassSummary,
                                                                             
                                      savePredictions = TRUE, classProbs = TRUE))

pr5 <- predict(gbm_m, test)
gbm_m <- confusionMatrix(pr5, test$cvd)
gbm_auc <- gbm_m$pred %>% roc_auc(obs, yes) %>% select(.estimate) %>% round(2)

#==============================================================================
# Plot AUC of Machine models

aa <- model3$pred %>% roc_curve(obs, yes) %>% mutate(.threshold = round(.threshold,2)) %>% 
  
  group_by(.threshold) %>% filter(row_number() == 1)



bb <- rf$pred %>% roc_curve(obs, yes) %>% mutate(.threshold = round(.threshold,2)) %>% 
  
  group_by(.threshold) %>% filter(row_number() == 1)


cc <- nn$pred %>% roc_curve(obs, yes) %>% mutate(.threshold = round(.threshold,2)) %>% 
  
  group_by(.threshold) %>% filter(row_number() == 1)

dd <- xg$pred %>% roc_curve(obs, yes) %>% mutate(.threshold = round(.threshold,2)) %>% 
  
  group_by(.threshold) %>% filter(row_number() == 1)

ee <- gbm_m$pred %>% roc_curve(obs, yes) %>% mutate(.threshold = round(.threshold,2)) %>% 
  
  group_by(.threshold) %>% filter(row_number() == 1)

names(aa) <- c('.threshold', 'glm_specificity', 'glm_sensitivity')

names(bb) <- c('.threshold', 'rf_specificity', 'rf_sensitivity')

names(cc) <- c('.threshold', 'nn_specificity', 'nn_sensitivity')

names(dd) <- c('.threshold', 'xg_specificity', 'xg_sensitivity')

names(ee) <- c('.threshold', 'gbm_specificity', 'gbm_sensitivity')

new <- inner_join(aa, bb, by  = '.threshold')

new <- new %>% inner_join(cc, by = '.threshold')

new <- new %>% inner_join(dd, by = '.threshold')

new <- new %>% inner_join(ee, by = '.threshold')

write.csv(new, "AUC curve.csv")

#Plot the AUC Curve

mina <- new %>% ggplot(.) + geom_step(aes(x = 1 - glm_specificity,
                                      
                                      y = glm_sensitivity), color = 'grey') +
 
  geom_step(aes(x = 1 - rf_specificity,  y = rf_sensitivity), color = 'red') +
  
  geom_step(aes(x = 1 - nn_specificity,  y = nn_sensitivity), color = 'blue') +
  
  geom_step(aes(x = 1 - xg_specificity,  y = xg_sensitivity), color = 'orange') + 
  
  geom_step(aes(x = 1 - gbm_specificity,  y = gbm_sensitivity), color = 'green') +
  
  theme(text = element_text(size = 8)) + aes(label = .threshold) +
  
  geom_text(aes(label = paste('MLR AUC is', glm_auc, '.'), x = 0.7, y = 0.5), color = 'grey') +
  
  geom_text(aes(label = paste('Random Forest AUC is', rf_auc, '.'), x = 0.7, y = 0.46), color = 'red') +
  
  geom_text(aes(label = paste('Neural Network AUC is', nn_auc, '.'), x = 0.7, y = 0.42), color = 'blue') + 
  
  geom_text(aes(label = paste('XGboost AUC is', glm_auc, '.'), x = 0.7, y = 0.38), color = 'orange') + 
  
  geom_text(aes(label = paste('Gradient Boosting AUC is', gbm_auc, '.'), x = 0.7, y = 0.34), color = 'green') +
  
  ggtitle('Comparison of AUC for CVD Risk Models') + labs(y = 'Sensistivity', x = '1 - Sepecificity') +
  
  theme_classic() + geom_abline(intercept = 0, slope = 1, linetype = 'dashed')
  
 pl <- mina  %>% ggplotly()
 install.packages('plotly')
 library(plotly)
#===============================================================================
#Characteristics of Study Population
 library(gtsummary)

xy <- copd_cvd2 %>% mutate(GENDER = ifelse(GENDER == 1, "Female", "Male"),
                           
                          cvd = ifelse(cvd == 1, "CVD", "No CVD"),
                          
                          alcohol = ifelse(alcohol == 1, "Yes", "No"),
                          
                          anxiety = ifelse(anxiety == 1, "Yes", "No"),
                          
                          depression = ifelse(depression == 1, "Yes", "No"),
                          
                          diabetes = ifelse(diabetes == 1, "Yes", "No"),
                          
                          emphysema = ifelse(emphysema == 1, "Yes", "No"),
                          
                          fam_hist_hd = ifelse(fam_hist_hd == 1, "Yes", "No"),
                          
                          hypertension = ifelse(hypertension == 1, "Yes", "No"),
                          
                          sleep_apnoea = ifelse(sleep_apnoea == 1, "Yes", "No"),
                          
                          smoking = ifelse(smoking == 1, "Yes", "No"))

xy <- xy %>% select(-mi)

xy %>% tbl_summary(by = cvd,
                   
                   statistic = list(
                     
                   age ~ "{mean} ({sd})",
                   
                   ALANINE_AMT ~ "{mean} ({sd})",
                   
                   CR_PROTEIN ~ "{mean} ({sd})",
                   
                   CREATININE ~ "{mean} ({sd})",
                   
                   FEV1 ~ "{mean} ({sd})",
                   
                   HBA1C ~ "{mean} ({sd})",
                   
                   BMI ~ "{mean} ({sd})",
                   
                   LDL ~ "{mean} ({sd})",
                   
                   PLATELET_C ~ "{mean} ({sd})",
                   
                   SERUM_SODIUM ~ "{mean} ({sd})",
                   
                   WEIGHT ~ "{mean} ({sd})")) %>% add_overall() %>% 
                   
                   add_stat_label(label = all_continuous() ~ "Mean (IQR)")
  
#==============================================================================
