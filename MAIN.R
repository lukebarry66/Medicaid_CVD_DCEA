#------------------------------------------------------------------------------#
#---------------------Purpose: MSM DCEA Medicaid-CVD---------------------------#
#---------------------------------7-MAIN---------------------------------------#
#------------------------------------------------------------------------------#

clear()

##### Load/install packages ----
library(pacman)
p_load(dplyr
       , magrittr
       , tidyr
       , cleaR
       , tidyverse
       , here
       , heemod
       , labelled
       , officer
       , flextable
       , abind
       , ineq
       , openxlsx
       , gridExtra
       , ggplot2
       , DataExplorer
       , nhanesA
       , colourpicker
       , viridis
       , survey
       , huxtable
       , miceRanger
       , writexl
       , srvyr
       , gtsummary
       , hesim
       , stringr
       , Amelia
       , VIM
       , naniar
       , patchwork
       , haven
       , formattable
       , glue
       , ggpubr
       , cowplot
       , scales)

##### Set environment ----
seed <- 12345
set.seed(seed)

##### PARAMETERS ----
{
  #### Scenario parameters ----
  
  # Set .Deterministic == 1 for .Deterministic and anything else for probabilistic
  .Deterministic    <- 1   
  # Set .one_way == 1 to conduct one-way sensitivity analysis 
  .one_way          <- 1 
  # Set .raw_data == 1 in order to store individual level data
  .raw_data         <- 0 
  
  # Set .ATT == 1 to examine just those receiving treatment (Medicaid expansion) and treatment effects 
  # otherwise whole NHANES sample, not just those who receive medicaid - NB option for distributional analysis
  .ATT              <- 0 
  # Change the strata by which the ASCVD weights are calculated
  .ascvd_cat        <- c("age", "female")                     
  # List of variables for which to calculate summaries
  .grouping_vars    <- c("age_cat", "race", "education", "fam_income", "medicaid_exp")
  # List of equity-relevant variables
  .equity_vars      <- c("Race/Ethnicity", "Education", "Family Income")
  # Set .stop_age (default == 64) according to the maximum starting age of individuals in the model
  # https://www.medicaid.gov/medicaid/eligibility/index.html
  .stop_age         <- 64
  # Set .start_age (default == 19) according to the minimum starting age of individuals in the model
  .start_age        <- 19
  # Set .trunc_age (default == 85) according to the age at which the model should end per person
  .trunc_age        <- 85
  # Set .fpl == 1.38 to match threshold for Medicaid expansion qualification for family income to federal poverty ratio
  .fpl              <- 1.38  
  # Set .fpl == 1.5 to reallocate non-oop costs for those receiving medicaid in the medicaid scenario to those with income >= 1.5 of the FPL
  .fpl_reallocate   <- 1.5                                                       
  
  # Set .half == 1 to examine the effect of the half-cycle correction
  .half             <- 1
  # Set .rate == 1 to calculate rate of MI, stroke etc instead of counts
  .rate             <- 1
  # Set .ir_py == 100000 to set the denominator for the incidence rate 100,000 person-years
  .ir_py            <- 100000
  # Set .excl_death == 1 to exclude deaths from calculation of person_years
  .excl_death       <- 1

  # Set .medicaid_cost == 1 to examine the effect of Medicaid on healthcare costs in the treatment scenario
  .medicaid_cost    <- 1
  # Set .gc == 1 to examine the effect government cost of administering Medicaid in the treatment scenario
  .gc               <- 1                                                          
  # Set .c_absent == 1 to include the cost of workplace absenteeism from an MI or Stroke
  .c_absent         <- 1
  # Set .c_disab == 1 to include the cost of short-term disability from an MI or Stroke
  .c_disab          <- 1
  # Set .prem_mort == 1 to include the cost of premature death in terms of lost annual earnings
  .c_prem_mort      <- 1
  
  # Set .hba1c == 1 for those who newly receive medicaid as part of the medicaid expansion to have lower hba1c and thus diabetes
  .hba1c            <- 1                                                          
  # Set .bp == 1 for those who newly receive medicaid as part of the medicaid expansion to have lower blood pressure
  .bp               <- 1    

  # Re-run Multiple Imputation (.mi == 1); otherwise use previously saved MI datasets (MI datasets are used in both cases)
  .mi               <- 1
  # may cause issues if updated "m" [number of imputed datasets] does not match "m" from previously saved MI
  
  #### Sensitivity analysis parameters  ----
  
  # Define the % swing in baseline parameter values for one-way sensitivity analysis
  delta_one_way <- 0.20  # 20% change
  # Define range of epsilon (inequality aversion) values over which to perform the distributional cost-effectiveness analysis
  epsilon_values  <- seq(0, 10, by = 0.3)
  # Define range of WTP values over which to perform the distributional cost-effectiveness analysis
  wtp_values <- seq(50000, 250000, by = 10000)
  
  #### Progress Display parameters ----
  
  # Set .track_prob == 1 to track whether cycle probabilities sum to 1
  .track_prob       <- 0
  # Set .track_cycle == 1 to track simulation progress as % cycles completed
  .track_cycle      <- 0
  # Set .track_boot == 1 to track bootstrapping progress as % replications completed
  .track_boot       <- 1
  # Set .track_dsa == 1 to track deterministic sensitivity analysis progress as % variables completed
  .track_dsa        <- 1
  # Set timer == 1 to track time of code execution
  .timer            <- 1
  
  #### Model parameters ----
  
  # number of simulated individuals
  n.i          <- 100
  # number of cycles (years) [default to 67; this allows everyone in the model to turn 85, i.e. 19 + 67 = 86]
  n.t          <- 5
  # number of bootstrapped model iterations 
  bootsize     <- 2
  # Number of imputed datasets
  m            <- 2
  # Number of iterations per imputed dataset
  m_iter       <- 2
  # set willingness-to-pay threshold [default should be $100,000 as recommended by the US Institute for Clinical and Economic Research]
  wtp          <- 100000                                                          
  # model state names
  v.n          <- c("No_CVD", 
                    "MI", 
                    "Stroke", 
                    "History_CVD", 
                    "CVD_death", 
                    "Non-CVD_death")  
  # number of states (6)
  n.s          <- length(v.n)  
  # discount rate for QALYs (a.k.a. annual utility); "e" refers more broadly to "effect"
  d.e          <- 0.03           
  # discount rate for annual healthcare costs
  d.c          <- 0.03       
  # strategy names
  v.Trt        <- c("no_medicaid_expansion", 
                    "medicaid_expansion") 
  # Define the duration of each cycle (assuming it's constant)
  cycle_length <- 1                 
  # Calculate the half-cycle correction factor
  half_cycle_correction <- 0.5 * cycle_length 
  # converting 95% CI's to SE's
  ci2se   <- 3.92
  # if no variability for estimate provided, use 20% of mean
  est_se  <- 0.2
  # Inequality aversion parameter
  epsilon <- 0.5                                                               
  
  #### Intervention parameters ----
  
  # absolute reduction in systolic blood pressure after receiving medicaid [−3.03 mmHg; 95% CI, −5.33 mmHg to − 0.73 mmHg] from DOI: 10.1007/s11606-020-06417-6
  delta_sys_bp    <- list(mean = 3.03, lci = 0.73, uci = 5.33)  
  # absolute reduction in hba1c after receiving medicaid [−0.14 percentage points [pp]; 95% CI, −0.24 pp to −0.03 pp;] from DOI: 10.1007/s11606-020-06417-6
  delta_hba1c     <- list(mean = 0.14, lci = 0.03, uci = 0.24)  
  
  #### Transition parameters ----
  
  # "history of CVD" multiplier (see Basu, 2017)
  delta_CVD_history <- list(shape = 3.842, scale = 0.521, mean = 2) 
  
  #### Cost parameters ----
  
  # Annual healthcare costs, as a function of CVD status among other variables, 
  # are estimated per person using the "Costs" function below using
  # the Morey et al 2021 equation (DOI: 10.1161/CIRCOUTCOMES.120.006769)
  
  # Inflation adjustment for Medicaid Admin costs (priced in 2019 USD) using Health PCE [https://apps.bea.gov/]
  pce_hc_2019      <- 1.041
  # Inflation adjustment for Morey et al (2021) healthcare costs (priced in 2017 USD) using Health PCE [https://apps.bea.gov/]
  pce_hc_2017      <- 1.074
  # Inflation adjustment for Song et al (2015) lost productivity costs (priced in 2013 USD) using Overall PCE [https://apps.bea.gov/]
  pce_all_2013     <- 1.141 
  
  # scaling factor ($29,545,128,289 / $645,708,814,950) for the government cost of administering the Medicaid (https://www.medicaid.gov/state-overviews/scorecard/annual-medicaid-chip-expenditures/index.html)
  c_gov_admin      <- 0.046  
  # Medicaid spending per newly eligible adult enrollee (https://www.kff.org/medicaid/state-indicator/medicaid-spending-per-enrollee/)
  c_enrollee       <- 5225
  # Medicaid administration spending per newly eligible adult enrollee in 2021 (https://www.kff.org/medicaid/state-indicator/medicaid-spending-per-enrollee/)
  c_admin_enrollee_mean  <- c_enrollee * c_gov_admin * pce_hc_2019 # probabilistic component, with se = mean*0.2, added in cost function
  c_admin_enrollee_se    <- c_admin_enrollee_mean * est_se
  c_admin_enrollee       <- list(mean = (c_admin_enrollee_mean), shape = ((c_admin_enrollee_mean^2)/(c_admin_enrollee_se^2)), rate = ((c_admin_enrollee_mean)/(c_admin_enrollee_se^2)))
  
  # scaling factor for mi and stroke costs in the year of the event vs the years after the event
  delta_c_y1    <- list(mean = 0.33, shape1 = 16.5, shape2 = 34.3) 
  
  # change in costs from medicaid (oop + other = total) (assumed 2019 prices); https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2809604
  delta_c_total    <- list(mean = 0.105, shape1 = 2.10, shape2 = 17.81) 
  # proportion of total costs that are related to CVD 
  cvd_prop         <- 0.55 # Medicaid spending on heart disease was approx. 3% in 2013 (https://doi.org/10.1097/mlr.0000000000000928)
  # proportion of total costs that are covered by other (e.g. insurer), i.e. not OOP costs 
  non_oop_prop     <- 1.13 
  # Average out-of-pocket costs (USD 2017) per person in pre-expansion states (aged 19-64 and 138% below the FPL); http://dx.doi.org/10.1136/bmj.m40 
  c_oop            <- 429 * pce_hc_2017 # no probabilistic component

  # 2021 US average annual earnings (from OECD 2022) for premature mortality before retirement at age 65
  annual_earnings  <- 74000 
  # Annual per person absenteeism costs
  cost_absenteeism <- list(shape = 3.739, rate = 0.015, mean = 244)
  # Annual per person short-term disability costs
  cost_disability  <- list(shape = 328.882, rate = 0.367, mean = 897)  
  
  # Coefficients for odds of non-zero healthcare costs for Morey et al 2021 equation
  odds_hc <- list(
    int      = list(mean = 1.653, lci = 1.168, uci = 2.339) # Intercept for estimating odds of non-zero healthcare costs 
    , mi       = list(mean = 2.937, lci = 1.798, uci = 4.798) # non-zero healthcare costs for those with an MI from Morey et al (2021)
    , stroke   = list(mean = 0.948, lci = 0.354, uci = 2.540) # non-zero healthcare costs for those with a stroke from Morey et al (2021) - see "Costs" function below   
    , hf       = list(mean = 1.830, lci = 0.402, uci = 8.325) # heart failure
    , cd       = list(mean = 2.389, lci = 1.269, uci = 4.499) # cardiac_dysrhythmia
    , ang      = list(mean = 1.472, lci = 0.566, uci = 3.828) # angina
    , pad      = list(mean = 1.334, lci = 0.617, uci = 2.888) # peripheral_artery_disease
    , diab     = list(mean = 4.382, lci = 3.401, uci = 5.630) # diabetes
    , a_25     = list(mean = 1.171, lci = 1.036, uci = 1.323) # age_cat == "25-44" 
    , a_45     = list(mean = 2.121, lci = 1.834, uci = 2.453) # age_cat == "45-64"
    , a_65     = list(mean = 3.728, lci = 2.076, uci = 6.695) # age_cat == "65+"
    , fem      = list(mean = 2.369, lci = 2.153, uci = 2.607) # female == 1
    , race_w   = list(mean = 1.627, lci = 1.452, uci = 1.824) # race == "white"
    , race_b   = list(mean = 1.013, lci = 0.885, uci = 1.161) # race == "black"
    , race_a   = list(mean = 1.025, lci = 0.848, uci = 1.239) # race == "asian"
    , race_o   = list(mean = 1.461, lci = 1.134, uci = 1.881) # race == "other_race"
    , medicare = list(mean = 1.714, lci = 0.918, uci = 3.2  ) # medicare
    , medicaid = list(mean = 0.959, lci = 0.834, uci = 1.104) # medicaid
    , uninsur  = list(mean = 0.355, lci = 0.313, uci = 0.403) # uninsured
    , inc_np   = list(mean = 0.905, lci = 0.726, uci = 1.128) # fam_income == "near_poor"
    , inc_l    = list(mean = 0.919, lci = 0.789, uci = 1.071) # fam_income == "low"
    , inc_m    = list(mean = 0.960, lci = 0.843, uci = 1.092) # fam_income == "medium"
    , inc_h    = list(mean = 1.344, lci = 1.127, uci = 1.601) # fam_income == "high"
    , ed_gh    = list(mean = 1.178, lci = 1.029, uci = 1.348) # education == "ged_hs"      
    , ed_ab    = list(mean = 1.726, lci = 1.462, uci = 2.038) # education == "associate_bachelor"
    , ed_md    = list(mean = 2.124, lci = 1.684, uci = 2.678) # education == "master_doctorate"
    , bmi_norm = list(mean = 0.985, lci = 0.709, uci = 1.369) # bmi_cat == "normal_weight"
    , bmi_ovw  = list(mean = 1.015, lci = 0.726, uci = 1.417) # bmi_cat == "overweight"  
    , bmi_obe  = list(mean = 1.233, lci = 0.873, uci = 1.74 ) # bmi_cat == "obese"     
    , cci1     = list(mean = 2.037, lci = 1.701, uci = 2.441) # cci == 1
    , cci2     = list(mean = 3.920, lci = 2.08 , uci = 7.387) # cci == 2
    , cci3     = list(mean = 8.287, lci = 3.163, uci = 21.715) # cci > 2
  )
  
  # Coefficients for mean of non-zero healthcare costs for Morey et al 2021 equation
  hc_cost <- list(
    int      = list(mean = 2600.846, lci = 1902.747, uci = 3555.072) # Intercept for estimating odds of non-zero healthcare costs 
    , mi       = list(mean = 1.177, lci = 1.032, uci = 1.343) # non-zero healthcare costs for those with an MI from Morey et al (2021)
    , stroke   = list(mean = 1.016, lci = 0.871, uci = 1.186) # non-zero healthcare costs for those with a stroke from Morey et al (2021) - see "Costs" function below   
    , hf       = list(mean = 0.948, lci = 0.810, uci = 1.110) # heart failure 
    , cd       = list(mean = 1.449, lci = 1.286, uci = 1.634) # cardiac_dysrhythmia
    , ang      = list(mean = 1.442, lci = 1.172, uci = 1.774) # angina
    , pad      = list(mean = 1.425, lci = 1.23 , uci = 1.651) # peripheral_artery_disease
    , diab     = list(mean = 1.676, lci = 1.471, uci = 1.909) # diabetes
    , a_25     = list(mean = 1.403, lci = 1.213, uci = 1.622) # age_cat == "25-44" 
    , a_45     = list(mean = 2.104, lci = 1.815, uci = 2.439) # age_cat == "45-64"
    , a_65     = list(mean = 2.240, lci = 1.711, uci = 2.932) # age_cat == "65+"
    , fem      = list(mean = 1.145, lci = 1.06 , uci = 1.236) # female == 1
    , race_w   = list(mean = 1.180, lci = 1.059, uci = 1.315) # race == "white"
    , race_b   = list(mean = 1.079, lci = 0.953, uci = 1.221) # race == "black"
    , race_a   = list(mean = 1.034, lci = 0.866, uci = 1.235) # race == "asian"
    , race_o   = list(mean = 1.224, lci = 1.02 , uci = 1.469) # race == "other_race"
    , medicare = list(mean = 1.122, lci = 0.898, uci = 1.403) # medicare
    , medicaid = list(mean = 1.308, lci = 1.17 , uci = 1.462) # medicaid
    , uninsur  = list(mean = 0.619, lci = 0.511, uci = 0.749) # uninsured
    , inc_np   = list(mean = 0.958, lci = 0.847, uci = 1.084) # fam_income == "near_poor"
    , inc_l    = list(mean = 1.020, lci = 0.882, uci = 1.179) # fam_income == "low"
    , inc_m    = list(mean = 0.919, lci = 0.824, uci = 1.026) # fam_income == "medium"
    , inc_h    = list(mean = 0.960, lci = 0.853, uci = 1.08 ) # fam_income == "high"
    , ed_gh    = list(mean = 1.076, lci = 0.943, uci = 1.228) # education == "ged_hs"      
    , ed_ab    = list(mean = 1.114, lci = 0.963, uci = 1.289) # education == "associate_bachelor"
    , ed_md    = list(mean = 1.185, lci = 1.008, uci = 1.392) # education == "master_doctorate"
    , bmi_norm = list(mean = 0.822, lci = 0.653, uci = 1.035) # bmi_cat == "normal_weight"
    , bmi_ovw  = list(mean = 0.794, lci = 0.628, uci = 1.004) # bmi_cat == "overweight"  
    , bmi_obe  = list(mean = 0.922, lci = 0.73 , uci = 1.165) # bmi_cat == "obese"     
    , cci1     = list(mean = 1.558, lci = 1.358, uci = 1.788) # cci == 1
    , cci2     = list(mean = 2.582, lci = 1.94 , uci = 3.436) # cci == 2
    , cci3     = list(mean = 2.991, lci = 2.447, uci = 3.655) # cci > 2
  )
  
  #### Utility parameters ----
  
  # Annual utilities, as a function of CVD status among other variables,
  # are estimated per person using the "Effs" function below using
  # the Morey et al 2021 equation (DOI: 10.1161/CIRCOUTCOMES.120.006769)
  
  # Coefficients for mean of utility (SF-6D) for Morey et al 2021 equation
  u <- list(
    int      = list(mean =  0.817, lci =  0.798, uci =  0.835) # Intercept for estimating odds of non-zero healthcare costs 
    , mi       = list(mean = -0.012, lci = -0.021, uci = -0.003) # non-zero healthcare costs for those with an MI from Morey et al (2021)
    , stroke   = list(mean = -0.026, lci = -0.043, uci = -0.009) # non-zero healthcare costs for those with a stroke from Morey et al (2021) - see "Costs" function below   
    , hf       = list(mean = -0.045, lci = -0.066, uci = -0.024) # heart failure
    , cd       = list(mean = -0.024, lci = -0.035, uci = -0.012) # cardiac_dysrhythmia
    , ang      = list(mean = -0.051, lci = -0.071, uci = -0.031) # angina
    , pad      = list(mean = -0.035, lci = -0.055, uci = -0.014) # peripheral_artery_disease
    , diab     = list(mean = -0.033, lci = -0.038, uci = -0.028) # diabetes
    , a_25     = list(mean = -0.026, lci = -0.035, uci = -0.012) # age_cat == "25-44" 
    , a_45     = list(mean = -0.047, lci = -0.071, uci = -0.031) # age_cat == "45-64"
    , a_65     = list(mean = -0.052, lci = -0.055, uci = -0.014) # age_cat == "65+"
    , fem      = list(mean = -0.024, lci = -0.038, uci = -0.028) # female == 1
    , race_w   = list(mean = -0.021, lci = -0.026, uci = -0.016) # race == "white"
    , race_b   = list(mean = -0.003, lci = -0.009, uci =  0.004) # race == "black"
    , race_a   = list(mean = -0.008, lci = -0.016, uci = -0.001) # race == "asian"
    , race_o   = list(mean = -0.030, lci = -0.041, uci = -0.019) # race == "other_race"
    , medicare = list(mean = -0.007, lci = -0.022, uci =  0.007) # medicare
    , medicaid = list(mean = -0.066, lci = -0.073, uci = -0.058) # medicaid
    , uninsur  = list(mean = -0.016, lci = -0.023, uci = -0.009) # medicaid
    , inc_np   = list(mean =  0.009, lci = -0.002, uci =  0.020) # fam_income == "near_poor"
    , inc_l    = list(mean =  0.023, lci =  0.015, uci =  0.032) # fam_income == "low"
    , inc_m    = list(mean =  0.040, lci =  0.033, uci =  0.047) # fam_income == "medium"
    , inc_h    = list(mean =  0.063, lci =  0.055, uci =  0.072) # fam_income == "high"
    , ed_gh    = list(mean =  0.013, lci =  0.007, uci =  0.019) # education == "ged_hs"      
    , ed_ab    = list(mean =  0.016, lci =  0.009, uci =  0.022) # education == "associate_bachelor"
    , ed_md    = list(mean =  0.022, lci =  0.014, uci =  0.030) # education == "master_doctorate"
    , bmi_norm = list(mean =  0.007, lci = -0.009, uci =  0.023) # bmi_cat == "normal_weight"
    , bmi_ovw  = list(mean =  0.006, lci = -0.010, uci =  0.023) # bmi_cat == "overweight"  
    , bmi_obe  = list(mean = -0.015, lci = -0.031, uci =  0.001) # bmi_cat == "obese"     
    , cci1     = list(mean = -0.039, lci = -0.045, uci = -0.034) # cci == 1
    , cci2     = list(mean = -0.052, lci = -0.063, uci = -0.041) # cci == 2
    , cci3     = list(mean = -0.076, lci = -0.093, uci = -0.060) # cci > 2
  )
  
} # Close "PARAMETERS" section

##### Filenames set-up ----

# Set the date scalar
date_scalar <- format(Sys.Date(), "%y%m%d")

# Combine all table filenames to be saved into a list
table_filenames <- c("cea_table_short", "cea_table_long", "ir_table")

# name the version of the model being run
version <- paste0("_trunc", .trunc_age, "_ATT", .ATT,
                  "_", n.i, "_", n.t, "_", bootsize, "_", m,
                  "_Det", .Deterministic, "_seed", seed)

# Create filename for combining tables and figures in one excel file
results_tables          <- paste0("results_tables", "_", date_scalar,"_", version, ".xlsx")
results_graphs          <- paste0("results_graphs", "_", date_scalar,"_", version, ".pptx")
results_combined_output <- paste0("combined_output", "_", date_scalar,"_", version, ".rds")
results_dsa             <- paste0("dsa_output", "_", date_scalar,"_", version, ".rds")

#### DATA ----

# Use this to create non-CVD mortality rates from the CDC and NHANES using nhanesA package
#source(here("scripts", "1-load_data.R"))
nCVD_mort   <- readRDS(here("data", "output_data", "nCVD_mort.rds")) 
nhanes_orig <- readRDS(here("data", "output_data", "appended_nhanes_data_orig.rds")) 

# Create Table 1 and prepare dataset for analysis
#source(here("scripts", "2-clean_data.R"))
nhanes_full <- readRDS(here("data", "output_data", "nhanes.rds")) 

# Identify weighted population counts for 2013 (closest to year of expansion)
source(here("scripts", "3-weighted_pop.R"))

# Prep for imputation
nhanes <- nhanes_full %>%
  select(-WTSAF2YR
         , -SDMVPSU
         , -SDMVSTRA) %>%
  filter(!is.na(WTSAF), #Drop individuals without weights
         age <= .stop_age,
         age >= .start_age) 

# Create imputed datasets
source(here("scripts", "4-imputation.R"))
dataList <- readRDS(here("data", "output_data", "dataList.rds"))

#### FUNCTIONS ----

source(here("scripts", "5-functions.R"))

#### ANALYSIS ----

### Bootstrap the simulation ----

  # Set timer == 1 in "Scenario Analysis" to track model run-time
  if (.timer == 1) {
    # Start timing
    start_time <- Sys.time()
  }
  
  for (dataset_index in 1:m) {
    # Get the dataset for the current iteration
    current_data <- dataList[[paste("Dataset_", dataset_index, sep = "")]]
    
    # Initialize an empty dataframe for the combined output
    output <- data.frame()
    
    # Initialize an empty list for raw data
    raw_data <- list()
    
    for (replication in 1:bootsize) {
      # Using tryCatch to handle errors during the bootstrap function
      result <- tryCatch({
        # Run the Bootstrap_msm function for the current replication on the current dataset
        Bootstrap_msm(original_data = current_data, n.i = n.i)
      }, error = function(e) {
        # Optionally log the error or perform other error handling
        message(paste("Error on dataset", dataset_index, "replication", replication, ":", e$message))
        NULL  # Continue to the next iteration by returning NULL
      })
      
      # Check if result is not NULL before proceeding
      if (!is.null(result)) {
        # Add replication and dataset index variables to the result data frame
        result[[2]]$replication   <- replication
        result[[2]]$dataset_index <- dataset_index
        
        # Append the results to the combined data frame
        output <- bind_rows(output, result[[2]])
        
        if (.raw_data == 1) {
          # Add replication and dataset index variables to the result data frame
          result[[1]]$replication   <- replication
          result[[1]]$dataset_index <- dataset_index
          
          # Append the results to the raw data list
          raw_data <- c(raw_data, list(result[[1]]))
        }
      }
      
      # Display the progress of the simulation
      if (.track_boot == 1) {
        cat('\r', paste(round(replication/bootsize * 100), "% Bootstrap done for Dataset_", dataset_index, sep = " "))
      }
    }
    
    if (.raw_data == 1 && length(raw_data) > 0) {
      # Save the results
      saveRDS(raw_data, here("data", "output_data", glue("raw_data_", "Dataset_", dataset_index, "_{today()}.rds")))
      rm(raw_data)
    }
    
    # Save combined output after processing each dataset
    saveRDS(output, here("data", "output_data", glue("output_", "Dataset_", dataset_index, "_{today()}.rds")))
  

  
  if (.timer == 1) {
    # End timing
    end_time <- Sys.time()
    # Calculate elapsed time
    elapsed_time <- end_time - start_time
    print(elapsed_time)
  }
}  

# Initialize an empty dataframe to store the final combined results
combined_output <- data.frame()

  # Loop through each dataset index
  for (dataset_index in 1:m) {
    # Construct the filename
    filename <- here("data", "output_data", glue("output_", "Dataset_", dataset_index, "_{today()}.rds"))

    # Load the data from each file
    current_output <- readRDS(filename)
    
    # Append the loaded data to the final combined data frame
    combined_output <- bind_rows(combined_output, current_output)
  }
  
  combined_output <- combined_output %>%
    mutate(  diff_mi        = mi_trt - mi_ntrt
             , diff_stroke    = stroke_trt - stroke_ntrt
             , diff_CVDdeath  = CVD_death_trt - CVD_death_ntrt
             , diff_nCVDdeath = nCVD_death_trt - nCVD_death_ntrt
             , inc_cost       = tc_trt - tc_ntrt
             , inc_effect     = te_trt - te_ntrt
             , nhb_ntrt       = (te_ntrt) - ((tc_ntrt) / wtp)
             , nhb_trt        = (te_trt) - ((tc_trt) / wtp)
             , icer           = (tc_trt - tc_ntrt) / (te_trt - te_ntrt)
             , inhb           = (te_trt - te_ntrt) - ((tc_trt - tc_ntrt) / wtp)
             , inmb           = ((te_trt - te_ntrt) * wtp) - (tc_trt - tc_ntrt)
    )
  
  
  ## Save the final combined data frame to a new RDS file
  saveRDS(combined_output, here("data", "output_data", glue("combined_output", "_{today()}.rds")))
  

### Summarize the results  ----

#source the combined_output
combined_output <- readRDS(here("data", "output_data", glue("combined_output", "_{today()}.rds")))

# run results synthesis
source(here("scripts", "6-results.R"))

#### DSA ----

  # Set .Deterministic == 1 for .Deterministic and anything else for probabilistic
  .Deterministic    <- 1   
  # number of simulated individuals
  n.i          <- 100
  # number of cycles (years) [default to 67; this allows everyone in the model to turn 85, i.e. 19 + 67 = 86]
  n.t          <- 10
  # number of bootstrapped model iterations 
  bootsize     <- 1
  # Number of imputed datasets
  m            <- 1
  wtp <-100000

if (.one_way == 1) {
  
  # List of variable names and corresponding base values
  variable_names <- c(
    "delta_CVD_history$mean",
    "delta_hba1c$mean",
    "delta_sys_bp$mean",
    "odds_hc$mi$mean",
    "odds_hc$stroke$mean",
    "hc_cost$mi$mean",
    "hc_cost$stroke$mean",
    "cost_absenteeism$mean",
    "cost_disability$mean",
    "delta_c_total$mean",
    "delta_c_y1$mean",
    "c_admin_enrollee$mean",
    "annual_earnings",
    "u$mi$mean",
    "u$stroke$mean",
    "d.e",
    "d.c"
  )  # Add your variable names
  
  # List of variable names and corresponding base values
  plot_names <- c(
    "CVD history risk mutliplier",
    "Delta HbA1c",
    "Delta systolic BP",
    "Odds non-zero HC-MI",
    "Odds non-zero HC-Stroke",
    "Non-zero HC-MI",
    "Non-zero HC-Stroke",
    "Absenteeism Cost",
    "Disability Cost",
    "Delta HC-Medicaid",
    "Delta y1 HC",
    "Medcaid admin cost/enrollee",
    "Annual earnings",
    "Utility-MI",
    "Utility-Stroke",
    "Discount rate - Utility",
    "Discount rate - Costs"
  )  # Add your variable names
  
  # Define the base_values vector at the start of your script
  base_values <- c(
    delta_CVD_history$mean
    , delta_hba1c$mean
    , delta_sys_bp$mean
    , odds_hc$mi$mean
    , odds_hc$stroke$mean
    , hc_cost$mi$mean
    , hc_cost$stroke$mean
    , cost_absenteeism$mean
    , cost_disability$mean
    , delta_c_total$mean
    , delta_c_y1$mean
    , c_admin_enrollee$mean
    , annual_earnings
    , u$mi$mean
    , u$stroke$mean
    , d.e
    , d.c
  )  # Define your base values
  
  # Initialize a list to store results
  results <- list()
  
  # Loop through each variable
  for (i in seq_along(variable_names)) {
    perturbation_results <- list()  # Initialize a list to store perturbation results
    for (perturbation in c(1 - delta_one_way, 1 + delta_one_way)) {
      perturbed_values <- base_values  # Reset perturbed_values to base_values for each perturbation
      perturbed_values[i] <- base_values[i] * perturbation
      
      result <- dsa(perturbed_values)
      
      # Store the inmb_values for this perturbation
      perturbation_results[[paste(variable_names[i], "Perturbation", perturbation, sep = "_")]] <- result
    }
    
    # display the progress of the DSA
    if (.track_dsa == 1) {
      cat('\r', paste(round(i/length(variable_names) * 100), "% DSA done", sep = " "))                
    }
    
    # Store the perturbation_results for this variable
    results[[variable_names[i]]] <- perturbation_results
  }

  # Calculate INHB values for each variable
  INHB_values <- sapply(variable_names, function(var_name) {
    perturbation_results <- results[[var_name]]
    perturbed_values <- sapply(perturbation_results, function(result) result$inhb_values)
    max_inhb <- max(perturbed_values)
    min_inhb <- min(perturbed_values)
    return(max_inhb - min_inhb)
  })
  
  # Order variables based on the largest differences
  ordered_variables <- variable_names[order(INHB_values, decreasing = TRUE)]
  
  # Create a data frame for plotting without changing factor levels initially
  tornado_data <- data.frame(
    Variable = factor(ordered_variables, levels = ordered_variables),
    Lower_INHB = sapply(ordered_variables, function(var_name) {
      perturbation_results <- results[[var_name]]
      perturbed_values <- sapply(perturbation_results, function(result) result$inhb_values)
      min_inhb <- min(perturbed_values)
      return(min_inhb)
    }),
    Upper_INHB = sapply(ordered_variables, function(var_name) {
      perturbation_results <- results[[var_name]]
      perturbed_values <- sapply(perturbation_results, function(result) result$inhb_values)
      max_inhb <- max(perturbed_values)
      return(max_inhb)
    })
  )
  
  # Calculate the average INMB
  average_INHB <- summarised_output$inhb[summarised_output$category == "Total"]

  # Create a factor with plot_names labels maintaining the original order
  ordered_plot_names <- plot_names[match(ordered_variables, variable_names)]
  tornado_data$Variable <- factor(tornado_data$Variable, levels = unique(ordered_variables))
  levels(tornado_data$Variable) <- ordered_plot_names
  
  # Create the tornado diagram plot with the correct order and labels
  tornado_plot <- ggplot(tornado_data, aes(x = fct_rev(Variable), ymin = Lower_INHB, ymax = Upper_INHB)) +
    geom_linerange(position = position_dodge(width = 0.9), color = "blue") +
    geom_hline(yintercept = average_INHB, linetype = "dashed", color = "red") +
    geom_text(aes(x = 1, y = average_INHB, label = paste("Baseline INHB:", round(average_INHB))), 
              vjust = 3, hjust = 1.05, color = "black", size = 4) +  # Add INHB label
    labs(title = "", y = "Incremental Net Health Benefit", x = "Variable") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_flip()
  
  # Filter the top 5 most impactful variables and reverse the order
  tornado_plot_data <- tornado_data %>%
    #filter(Variable %in% head(ordered_variables, 20)) %>%
    arrange(desc(Variable))
  
  print(tornado_plot + geom_linerange(data = tornado_plot_data, color = "blue"))
  
}
