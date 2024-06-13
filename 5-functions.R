#------------------------------------------------------------------------------#
#---------------------Purpose: MSM DCEA Medicaid-CVD---------------------------#
#-------------------------------5-FUNCTIONS------------------------------------#
#------------------------------------------------------------------------------#


##### FUNCTIONS ----
{
  
  #### Missingness functions ----
  #to plot missing values
  gtmiss <- function(...){
    DataExplorer::profile_missing(...) %>% 
      arrange(pct_missing) %>% 
      print(n=Inf)
  }
  #e.g. gtmiss(dataset)
  
  #to list missing values
  ggmiss <- function(...){
    DataExplorer::plot_missing(...) 
  }
  #e.g. ggmiss(dataset)
  
  
  #### ASCVD weighting function ----
  # The ASCVD function predicts an individuals 10 year ASCVD risk and estimates their average risk within age, sex and race strata
  
  average_ascvd_risk <- function(mydata) {
    
    mydata <- mydata %>%
      
      # Estimate 10 year ASCVD risk
      mutate(ascvd_risk_10yr = case_when(
        # see (https://www.framinghamheartstudy.org/fhs-risk-functions/cardiovascular-disease-10-year-risk/) for details
        (race != "black" & female == 0) ~ (1 - 0.9144) * exp(12.344 * log(age) + 
                                                               0      * log(age) * log(age) + 
                                                               11.853 * log(total_cholesterol) + 
                                                               -2.664 * log(age) * log(total_cholesterol) +
                                                               -7.990 * log(hdl_cholesterol) + 
                                                               1.769  * log(age) * log(hdl_cholesterol) + 
                                                               (1.797 * log(mean_systolic_bp_orig) + 
                                                                  0      * log(age) * log(mean_systolic_bp_orig)) * treated + 
                                                               (1.764 * log(mean_systolic_bp_orig) + 
                                                                  0      * log(age) * log(mean_systolic_bp_orig)) * (1 - treated) + 
                                                               7.837  * smoking +
                                                               -1.795 * log(age) * smoking + 
                                                               0.658  * diabetes_orig + 
                                                               -(61.18)),
        (race != "black" & female == 1) ~ (1 - 0.9665) * exp(-29.799   * log(age) + 
                                                               4.884   * log(age)*log(age) + 
                                                               13.54   * log(total_cholesterol) + 
                                                               -3.114  * log(age) * log(total_cholesterol) +
                                                               -13.578 * log(hdl_cholesterol) + 
                                                               3.149   * log(age) * log(hdl_cholesterol) + 
                                                               (2.019  * log(mean_systolic_bp_orig) + 
                                                                  0       * log(age) * log(mean_systolic_bp_orig)) * treated + 
                                                               (1.957  * log(mean_systolic_bp_orig) + 
                                                                  0       * log(age) * log(mean_systolic_bp_orig)) * (1 - treated) + 
                                                               7.574   * smoking +
                                                               -1.665  * log(age) * smoking + 
                                                               0.661   * diabetes_orig + 
                                                               (-(-29.18))),
        (race == "black" & female == 0) ~ (1 - 0.8954) * exp(2.469    * log(age) + 
                                                               0      * log(age) * log(age) + 
                                                               0.302  * log(total_cholesterol) + 
                                                               0      * log(age) * log(total_cholesterol) +
                                                               -0.307 * log(hdl_cholesterol) + 
                                                               0      * log(age) * log(hdl_cholesterol) + 
                                                               (1.916 * log(mean_systolic_bp_orig) + 
                                                                  0      * log(age) * log(mean_systolic_bp_orig)) * treated + 
                                                               (1.809 * log(mean_systolic_bp_orig) + 
                                                                  0      * log(age) * log(mean_systolic_bp_orig)) * (1 - treated) + 
                                                               0.549  * smoking +
                                                               0      * log(age) * smoking + 
                                                               0.645  * diabetes_orig + 
                                                               -(19.54)),
        (race == "black" & female == 1) ~ (1 - 0.9533) * exp(17.114    * log(age) + 
                                                               0       * log(age) * log(age) + 
                                                               0.940   * log(total_cholesterol) + 
                                                               0       * log(age) * log(total_cholesterol) +
                                                               -18.920 * log(hdl_cholesterol) + 
                                                               4.475   * log(age)*log(hdl_cholesterol) + 
                                                               (29.291 * log(mean_systolic_bp_orig) +
                                                                  -6.432  * log(age) * log(mean_systolic_bp_orig)) * treated + 
                                                               (27.820 * log(mean_systolic_bp_orig) + 
                                                                  -6.087  * log(age) * log(mean_systolic_bp_orig)) * (1 - treated) + 
                                                               0.691   * smoking +
                                                               0       * log(age) * smoking + 
                                                               0.874   * diabetes_orig + 
                                                               -(86.61)),
        TRUE~NA_real_),
        
        ascvd_risk_10yr   = ifelse(ascvd_risk_10yr < 1, 
                                   # Equation can result in probabilities > 1 for very high ages; this puts a limit on it
                                   ascvd_risk_10yr, 0.999999999),               
        # Estimate individual weights to be applied to FHS MI or stroke risk within strata 
        ascvd_rate   = ((-log(1 - ascvd_risk_10yr))/(10)), 
        ascvd_risk   = 1 - exp(-(ascvd_rate) * 1)) %>%
      # Estimate average 10 year ASCVD risk within strata according to age and sex 
      group_by(across(all_of( .ascvd_cat ))) %>%
      mutate(avg_ascvd_risk = mean(ascvd_risk)) %>%
      ungroup() %>%
      select(-ascvd_risk,
             -ascvd_rate,
             -ascvd_risk_10yr)
    
    return(mydata) 
    
  }
  
  #### Probability function ----
  # The Probs function that updates the transition probabilities of every cycle is shown below.
  
  # Predict annual mi risk for those with no history of CVD from FHS monthly risk re-weighted by ASCVD event risk
  Probs <- function(mydata
                    , Trt = FALSE
                    , parameters = parameters
                    ) {
    
    mydata  <- mydata %>%
      # Estimate 10 year ASCVD risk
      mutate(ascvd_risk_10yr = case_when(
        # see (https://www.framinghamheartstudy.org/fhs-risk-functions/cardiovascular-disease-10-year-risk/) for details
        (race != "black" & female == 0) ~ (1 - 0.9144) * exp(12.344 * log(age) + 
                                                               0      * log(age) * log(age) + 
                                                               11.853 * log(total_cholesterol) + 
                                                               -2.664 * log(age) * log(total_cholesterol) +
                                                               -7.990 * log(hdl_cholesterol) + 
                                                               1.769  * log(age) * log(hdl_cholesterol) + 
                                                               (1.797 * log(mean_systolic_bp) + 
                                                                  0      * log(age) * log(mean_systolic_bp)) * treated + 
                                                               (1.764 * log(mean_systolic_bp) + 
                                                                  0      * log(age) * log(mean_systolic_bp)) * (1 - treated) + 
                                                               7.837  * smoking +
                                                               -1.795 * log(age) * smoking + 
                                                               0.658  * diabetes + 
                                                               -(61.18)),
        (race != "black" & female == 1) ~ (1 - 0.9665) * exp(-29.799   * log(age) + 
                                                               4.884   * log(age)*log(age) + 
                                                               13.54   * log(total_cholesterol) + 
                                                               -3.114  * log(age) * log(total_cholesterol) +
                                                               -13.578 * log(hdl_cholesterol) + 
                                                               3.149   * log(age) * log(hdl_cholesterol) + 
                                                               (2.019  * log(mean_systolic_bp) + 
                                                                  0       * log(age) * log(mean_systolic_bp)) * treated + 
                                                               (1.957  * log(mean_systolic_bp) + 
                                                                  0       * log(age) * log(mean_systolic_bp)) * (1 - treated) + 
                                                               7.574   * smoking +
                                                               -1.665  * log(age) * smoking + 
                                                               0.661   * diabetes + 
                                                               (-(-29.18))),
        (race == "black" & female == 0) ~ (1 - 0.8954) * exp(2.469    * log(age) + 
                                                               0      * log(age) * log(age) + 
                                                               0.302  * log(total_cholesterol) + 
                                                               0      * log(age) * log(total_cholesterol) +
                                                               -0.307 * log(hdl_cholesterol) + 
                                                               0      * log(age) * log(hdl_cholesterol) + 
                                                               (1.916 * log(mean_systolic_bp) + 
                                                                  0      * log(age) * log(mean_systolic_bp)) * treated + 
                                                               (1.809 * log(mean_systolic_bp) + 
                                                                  0      * log(age) * log(mean_systolic_bp)) * (1 - treated) + 
                                                               0.549  * smoking +
                                                               0      * log(age) * smoking + 
                                                               0.645  * diabetes + 
                                                               -(19.54)),
        (race == "black" & female == 1) ~ (1 - 0.9533) * exp(17.114    * log(age) + 
                                                               0       * log(age) * log(age) + 
                                                               0.940   * log(total_cholesterol) + 
                                                               0       * log(age) * log(total_cholesterol) +
                                                               -18.920 * log(hdl_cholesterol) + 
                                                               4.475   * log(age)*log(hdl_cholesterol) + 
                                                               (29.291 * log(mean_systolic_bp) +
                                                                  -6.432  * log(age) * log(mean_systolic_bp)) * treated + 
                                                               (27.820 * log(mean_systolic_bp) + 
                                                                  -6.087  * log(age) * log(mean_systolic_bp)) * (1 - treated) + 
                                                               0.691   * smoking +
                                                               0       * log(age) * smoking + 
                                                               0.874   * diabetes + 
                                                               -(86.61)),
        TRUE~NA_real_),
        ascvd_risk_10yr   = ifelse(ascvd_risk_10yr < 1, 
                                   # Equation can result in probabilities > 1 for very high ages; this puts a limit on it
                                   ascvd_risk_10yr, 0.9999999),               
        # Estimate individual weights to be applied to FHS MI or stroke risk within strata 
        ascvd_rate   = ((-log(1 - ascvd_risk_10yr))/(10)), 
        ascvd_risk   = 1 - exp(-(ascvd_rate) * 1),
        
        ## Estimate monthly FHS MI risk, convert to annual rate, apply 10 year ASCVD weights, convert to re-weighted annual MI risk for those without a history of CVD
        # Estimate monthly MI risk as a function of age within categories of gender
        monthly_mi_risk       = ifelse(female == 0,
                                       0.0001   * exp(0.0312 * age),      
                                       0.000008 * exp(0.0599 * age)),
        
        # Equation can result in probabilities > 1 for very high ages; this puts a limit on it
        monthly_mi_risk       = ifelse(monthly_mi_risk < 1, 
                                       monthly_mi_risk, 1),  
        
        # Convert monthly risk to annual rate (assumes constant rate over time)
        annual_mi_rate        = ((-log(1 - monthly_mi_risk))/(1 / 12)),  
        
        # Convert annual rate to annual risk (assumes constant rate over time)  
        annual_mi_risk        = 1 - exp(-(annual_mi_rate) * 1),   
        
        # Using Bayes formula to estimate the probability of an MI given an ASCVD event for age and sex categories 
        annual_ascvd_mi_risk  = annual_mi_risk / avg_ascvd_risk, 
        
        # Use Bayes formula to then estimate the probability of an ASCVD and MI event 
        mi_risk_noCVD         = annual_ascvd_mi_risk * ascvd_risk,          
        
        ## Estimate monthly FHS stroke risk, convert to annual rate, apply 10 year ASCVD weights, convert to re-weighted annual stroke risk for those without a history of CVD
        # Estimate monthly stroke risk as a function of age within categories of gender
        monthly_stroke_risk      = ifelse(female == 0, 
                                          0.000009 * exp(0.0622 * age),             
                                          0.000003 * exp(0.0741 * age)),
        
        # Equation can result in probabilities > 1 for very high ages; this puts a limit on it
        monthly_stroke_risk      = ifelse(monthly_stroke_risk < 1, 
                                          monthly_stroke_risk, 1),  
        
        # Convert monthly risk to annual rate (assumes constant rate over time)
        annual_stroke_rate       = ((-log(1 - monthly_stroke_risk))/(1 / 12)), 
        
        # Convert annual rate to annual risk (assumes constant rate over time) 
        annual_stroke_risk       = 1 - exp(-(annual_stroke_rate) * 1),  
        
        # Using Bayes formula to estimate the probability of an MI given an ASCVD event for age and sex categories 
        annual_ascvd_stroke_risk = annual_stroke_risk / avg_ascvd_risk,
        
        # Use Bayes formula to then estimate the probability of an ASCVD and MI event 
        stroke_risk_noCVD        = annual_ascvd_stroke_risk * ascvd_risk,          
        
        # Estimate MI mortality risk as a function of age within categories of gender
        mi_mortality_risk     = ifelse(female == 0, 
                                       0.0289 * exp(0.0269 * age),        
                                       0.0004 * exp(0.0706 * age)), 
        
        # Equation can result in probabilities > 1 for very high ages; this puts a limit on it
        mi_mortality_risk     = ifelse(mi_mortality_risk < 1, 
                                       mi_mortality_risk, 1),                        
        
        # Estimate stroke mortality risk as a function of age within categories of gender
        stroke_mortality_risk     = ifelse(female == 0, 
                                           0.0003 * exp(0.0782 * age),                
                                           0.0034 * exp(0.0428 * age)),
        
        # Equation can result in probabilities > 1 for very high ages; this puts a limit on it
        stroke_mortality_risk     = ifelse(stroke_mortality_risk < 1, 
                                           stroke_mortality_risk, 1),                
        
        # Convert annual MI risk to rate, apply "history of CVD multiplier" (see Basu, 2017)  and convert back to prob
        mi_risk_CVD            = 1 - exp((-(((-log(1 - mi_risk_noCVD))/(1)) * parameters[["CVD_history_multiplier"]])) * 1), # Conversion assumes constant rate over time
        
        # Convert annual stroke risk to rate, apply "history of CVD multiplier" (see Basu, 2017)  and convert back to prob
        stroke_risk_CVD        = 1 - exp((-(((-log(1 - stroke_risk_noCVD))/(1)) * parameters[["CVD_history_multiplier"]])) * 1) # Conversion assumes constant rate over time
      ) 
    
    # add non-CVD mortality data from CDC WONDER
    mydata <- left_join(mydata, 
                        nCVD_mort, 
                        by = c("age", "female"), 
                        relationship = "many-to-many")

    mydata  <- mydata %>%
      # Estimate non-CVD mortality risk
      mutate(
        # calculate annual non-CVD mortality rate (per person) as difference between all-cause mortality rate from CVD mortality rate in CDC WONDER database
        noncvd_mortality_rate  = case_when(.Deterministic == 1 ~ rate,
                                           .Deterministic != 1 ~ rnorm(n(), rate, se)),

        # Convert non-CVD rate to annual risk
        noncvd_mortality_risk  = 1 - exp((-noncvd_mortality_rate) * 1), # Conversion assumes constant rate over time
        
        # Equation can result in probabilities > 1 for very high ages; this puts a limit on it
        noncvd_mortality_risk  = ifelse(noncvd_mortality_risk < 1, 
                                        noncvd_mortality_risk, 1)                               
      )  %>%  
      # Keep only variables for transitions matrices
      select(stroke_risk_noCVD,
             stroke_risk_CVD,
             stroke_mortality_risk,
             mi_risk_noCVD,
             mi_risk_CVD,
             mi_mortality_risk,
             noncvd_mortality_risk)
    
    return(mydata) 
    
  }
  
  
  #### Cost function ----
  # The Costs function estimates the costs at every cycle
  
  Costs <- function(mydata,
                    cci,
                    cardiac_dysrhythmia,      
                    peripheral_artery_disease,
                    Trt = FALSE,
                    parameters = parameters
                    , t=t) {

    # Below formulas come from Morey et al (2021; DOI: 10.1161/CIRCOUTCOMES.120.006769) using MEPS 2016 data
    mydata  <- mydata %>%
      mutate( # calculate the odds of having non-zero health expenditures 
        # to avoid non-sensical values draw from log values and exponentiate after
        odds_nonzero_hc = ((parameters[["odds_nonzero_hc_int"]]) * # baseline utility                                                 
                             ifelse(mi == 1,                             parameters[["odds_nonzero_hc_mi"]], 1) *
                             ifelse(stroke == 1,                         parameters[["odds_nonzero_hc_stroke"]], 1) *
                             ifelse(heart_failure == 1,                  parameters[["odds_nonzero_hc_hf"]], 1) *
                             # do not have data in nhanes on dc and pad so set == 0 for eveyrone when using function exp
                             ifelse(cardiac_dysrhythmia == 1,            parameters[["odds_nonzero_hc_cd"]]  , 1) *      
                             ifelse(angina == 1,                         parameters[["odds_nonzero_hc_ang"]] , 1) *
                             ifelse(peripheral_artery_disease == 1,      parameters[["odds_nonzero_hc_pad"]] , 1) *
                             ifelse(diabetes == 1,                       parameters[["odds_nonzero_hc_diab"]], 1) *
                             # baseline = "18-24"
                             ifelse((age_cat == "25-44"),                parameters[["odds_nonzero_hc_a_25"]], 1) *       
                             ifelse((age_cat == "45-64"),                parameters[["odds_nonzero_hc_a_45"]], 1) *
                             ifelse((age_cat == "65+"),                  parameters[["odds_nonzero_hc_a_65"]], 1) *
                             ifelse(female == 1,                         parameters[["odds_nonzero_hc_fem"]], 1) *
                             # baseline = "hispanic"
                             ifelse((race == "white"),                   parameters[["odds_nonzero_hc_race_w"]], 1) *       
                             ifelse((race == "black"),                   parameters[["odds_nonzero_hc_race_b"]], 1) *
                             ifelse((race == "asian"),                   parameters[["odds_nonzero_hc_race_a"]], 1) *
                             ifelse((race == "other_race"),              parameters[["odds_nonzero_hc_race_o"]], 1) *
                             # baseline = "Private insurance"
                             ifelse((insurance == "medicare"),      parameters[["odds_nonzero_hc_medicare"]], 1) *      
                             # referred to as "other public" but assumed to be "medicaid" 
                             ifelse((insurance == "medicaid" | 
                                       insurance == "other plan"),  parameters[["odds_nonzero_hc_medicaid"]], 1) *       
                             ifelse((insurance == "uninsured"),     parameters[["odds_nonzero_hc_uninsur"]], 1) *
                             # baseline = "Poor"
                             ifelse((fam_income == "near_poor"),         parameters[["odds_nonzero_hc_inc_np"]], 1) *       
                             ifelse((fam_income == "low"),               parameters[["odds_nonzero_hc_inc_l"]], 1) *
                             ifelse((fam_income == "medium"),            parameters[["odds_nonzero_hc_inc_m"]], 1) *
                             ifelse((fam_income == "high"),              parameters[["odds_nonzero_hc_inc_h"]], 1) *
                             # baseline = "no degree"
                             ifelse((education == "ged_hs"),             parameters[["odds_nonzero_hc_ed_gh"]], 1) *       
                             ifelse((education == "associate_bachelor"), parameters[["odds_nonzero_hc_ed_ab"]], 1) * 
                             ifelse((education == "master_doctorate"),   parameters[["odds_nonzero_hc_ed_md"]], 1) *
                             # baseline = "Underweight" == "1" in NHANES
                             ifelse((bmi_cat == "normal_weight"),        parameters[["odds_nonzero_hc_bmi_norm"]], 1) *        
                             ifelse((bmi_cat == "overweight"),           parameters[["odds_nonzero_hc_bmi_ovw"]], 1) *
                             ifelse((bmi_cat == "obese"),                parameters[["odds_nonzero_hc_bmi_obe"]] , 1) *
                             # baseline = 0
                             ifelse((cci == 1),                          parameters[["odds_nonzero_hc_cci1"]], 1) *       
                             ifelse((cci == 2),                          parameters[["odds_nonzero_hc_cci2"]], 1) *
                             ifelse((cci > 2),                           parameters[["odds_nonzero_hc_cci3"]], 1)),
        
        
        # calculate the probability of having non-zero health expenditures
        prob_nonzero_hc = (odds_nonzero_hc/(1+odds_nonzero_hc)),
        
        # calculate expected healthcare expenditure 
        # baseline utility 
        nonzero_hc_cost = ((parameters[["nonzero_hc_cost_int"]]) *                                        
                             ifelse(mi == 1,                             parameters[["nonzero_hc_cost_mi"]], 1) *
                             ifelse(stroke == 1,                         parameters[["nonzero_hc_cost_stroke"]], 1) *
                             ifelse(heart_failure == 1,                  parameters[["nonzero_hc_cost_hf"]], 1) *
                             # do not have data in nhanes on dc and pad so set == 0 for eveyrone when using function
                             ifelse(cardiac_dysrhythmia == 1,            parameters[["nonzero_hc_cost_cd"]]  , 1) *            
                             ifelse(angina == 1,                         parameters[["nonzero_hc_cost_ang"]] , 1) *
                             ifelse(peripheral_artery_disease == 1,      parameters[["nonzero_hc_cost_pad"]] , 1) *
                             ifelse(diabetes == 1,                       parameters[["nonzero_hc_cost_diab"]], 1) *
                             # baseline = "18-24"
                             ifelse((age_cat == "25-44"),                parameters[["nonzero_hc_cost_a_25"]], 1) *               
                             ifelse((age_cat == "45-64"),                parameters[["nonzero_hc_cost_a_45"]], 1) *        
                             ifelse((age_cat == "65+"),                  parameters[["nonzero_hc_cost_a_65"]], 1) *        
                             ifelse(female == 1,                         parameters[["nonzero_hc_cost_fem"]], 1) * 
                             # baseline = "hispanic"
                             ifelse((race == "white"),                   parameters[["nonzero_hc_cost_race_w"]], 1) *               
                             ifelse((race == "black"),                   parameters[["nonzero_hc_cost_race_b"]], 1) *        
                             ifelse((race == "asian"),                   parameters[["nonzero_hc_cost_race_a"]], 1) *        
                             ifelse((race == "other_race"),              parameters[["nonzero_hc_cost_race_o"]], 1) *
                             # baseline = "Private insurance"
                             ifelse((insurance == "medicare"),      parameters[["nonzero_hc_cost_medicare"]], 1) *               
                             # referred to as "other public"; majority would be "medicaid" 
                             ifelse((insurance == "medicaid" | 
                                       insurance == "other_plan"),  parameters[["nonzero_hc_cost_medicaid"]], 1) *                
                             ifelse((insurance == "uninsured"),     parameters[["nonzero_hc_cost_uninsur"]], 1) *
                             # baseline = "Poor"
                             ifelse((fam_income == "near_poor"),         parameters[["nonzero_hc_cost_inc_np"]], 1) *              
                             ifelse((fam_income == "low"),               parameters[["nonzero_hc_cost_inc_l"]] , 1) *        
                             ifelse((fam_income == "medium"),            parameters[["nonzero_hc_cost_inc_m"]] , 1) *        
                             ifelse((fam_income == "high"),              parameters[["nonzero_hc_cost_inc_h"]] , 1) *
                             # baseline = "no degree"
                             ifelse((education == "ged_hs"),             parameters[["nonzero_hc_cost_ed_gh"]], 1) *              
                             ifelse((education == "associate_bachelor"), parameters[["nonzero_hc_cost_ed_ab"]], 1) *         
                             ifelse((education == "master_doctorate"),   parameters[["nonzero_hc_cost_ed_md"]], 1) *
                             # baseline = "Underweight"
                             ifelse((bmi_cat == "normal_weight"),        parameters[["nonzero_hc_cost_bmi_norm"]], 1) *               
                             ifelse((bmi_cat == "overweight"),           parameters[["nonzero_hc_cost_bmi_ovw"]] , 1) *        
                             ifelse((bmi_cat == "obese"),                parameters[["nonzero_hc_cost_bmi_obe"]] , 1) *
                             # baseline = 0
                             ifelse((cci == 1),                          parameters[["nonzero_hc_cost_cci1"]], 1) *               
                             ifelse((cci == 2),                          parameters[["nonzero_hc_cost_cci2"]], 1) *        
                             ifelse((cci > 2),                           parameters[["nonzero_hc_cost_cci3"]], 1)),
        
        # calculate annual healthcare expenditure times the probability of non-zero healthcare expenditures 
        hc_cost = (prob_nonzero_hc * nonzero_hc_cost),
        
        # inflate to 2021 prices using the Health PCE 
        hc_cost = (hc_cost * pce_hc_2017),
        
        # government administration costs as part of treatment scenario for those receiving medicaid who are under 65 years (at age 65 years everyone switches to Medicare)
        gov_admin_cost = case_when(Trt == TRUE & medicaid_exp == 1 & age < 65 & .gc == 1 ~ parameters[["cost_admin_enrollee"]] * cvd_prop,
                                   Trt == FALSE ~ 0,
                                   TRUE ~ 0),

        # change in total healthcare costs following receiving medicaid in medicaid expansion
        delta_total_hc_cost = case_when((Trt == TRUE & .medicaid_cost == 1) ~ (hc_cost * parameters[["delta_total_cost"]] * medicaid_exp * cvd_prop)
                                        , (Trt == FALSE) ~ 0
                                        , TRUE ~ 0),
        
        delta_other_hc_cost = (delta_total_hc_cost * non_oop_prop),
        delta_oop_hc_cost   = (delta_total_hc_cost - delta_other_hc_cost),
      ) %>% 
      
      # Keep only variables for transitions matrices
      select(fpl
             , fam_income
             , medicaid_exp
             , hc_cost
             , gov_admin_cost
             , delta_total_hc_cost
             , delta_other_hc_cost
             , delta_oop_hc_cost)

    # Calculate number of individuals above or equal to 150% of the FPL and apply non-oop costs to them; http://dx.doi.org/10.1136/bmj.m40 (Supplement, page 3)
    length_non_medicaid <- length(mydata$medicaid_exp[mydata$fpl >= .fpl_reallocate])
    # amount fo gov costs to apportion to all those with income >= 150% FPL
    mean_gov_cost_society <- sum(mydata$gov_admin_cost, na.rm = T) / length_non_medicaid
    # amount fo non-oop costs to apportion to all those with income >= 150% FPL
    mean_noop_cost_society <- sum(mydata$delta_other_hc_cost, na.rm = T) / length_non_medicaid

    mydata <- mydata %>%
      # combining perspective and ATT scenarios with and without treatment
      mutate(redist_gov_cost = case_when(.ATT == 1 & Trt == TRUE ~ ((gov_admin_cost) * medicaid_exp)
                                         , .ATT != 1 & Trt == TRUE ~ (mean_gov_cost_society * (1-medicaid_exp) * (ifelse(fpl >= .fpl_reallocate,1,0)))
                                         , TRUE ~ 0),
             redist_noop_cost = case_when(.ATT == 1 & Trt == TRUE ~ ((delta_other_hc_cost) * medicaid_exp)
                                          , .ATT != 1 & Trt == TRUE ~ (mean_noop_cost_society * (1-medicaid_exp) * (ifelse(fpl >= .fpl_reallocate,1,0)))
                                          , TRUE ~ 0)
             ) %>%
      # Keep only variables for transitions matrices
      select(hc_cost
             , gov_admin_cost
             , delta_total_hc_cost
             , delta_other_hc_cost
             , delta_oop_hc_cost
             , redist_gov_cost
             , redist_noop_cost
             ) %>%
      mutate(delta_nonoop_hc_cost = delta_other_hc_cost) %>%
      select(-delta_other_hc_cost)
    
    return(mydata)
    
  }
  
  #### Utilities function ---- 
  # The Effs function to update the utilities at every cycle
  
  Effs <- function(mydata,
                   cardiac_dysrhythmia,      
                   peripheral_artery_disease,
                   cci, 
                   Trt = FALSE,
                   parameters = parameters) {
    
    # Below formulas come from Morey et al (2021; DOI: 10.1161/CIRCOUTCOMES.120.006769) using MEPS 2016 data
    mydata  <- mydata %>%
      # calculate annual utility (and by default QALYs; cycle length = 1 year) using DOI: 10.1161/CIRCOUTCOMES.120.006769
      mutate(utility = ((parameters[["util_int"]]) +     # baseline utility                                                          
                          mi                                                     * (parameters[["util_mi"]]) +
                          stroke                                                 * (parameters[["util_stroke"]]) +
                          heart_failure                                          * (parameters[["util_hf"]]) +
                          # do not have data in nhanes on dc and pad so set == 0 for eveyrone when using function
                          cardiac_dysrhythmia                                    * (parameters[["util_cd"]]  ) +     
                          angina                                                 * (parameters[["util_ang"]] ) +
                          peripheral_artery_disease                              * (parameters[["util_pad"]] ) +
                          diabetes                                               * (parameters[["util_diab"]]) +
                          (age_cat == "<25")                                     * 0 +                # baseline = "18-24"
                          (age_cat == "25-44")                                   * (parameters[["util_a_25"]]) +
                          (age_cat == "45-64")                                   * (parameters[["util_a_45"]]) +
                          (age_cat == "65+")                                     * (parameters[["util_a_65"]]) +
                          female                                                 * (parameters[["util_fem"]]) +  
                          (race == "hispanic")                                   * 0 +                # baseline = "hispanic"
                          (race == "white")                                      * (parameters[["util_race_w"]]) + 
                          (race == "black")                                      * (parameters[["util_race_b"]]) +
                          (race == "asian")                                      * (parameters[["util_race_a"]]) +
                          (race == "other_race")                                 * (parameters[["util_race_o"]]) +
                          (insurance == "private")                          * 0 +                # baseline = "Private insurance"
                          (insurance == "medicare")                         * (parameters[["util_medicare"]]) +
                          (insurance == "medicaid" | 
                             insurance == "other_plan")                     * (parameters[["util_medicaid"]]) +  # referred to as "other public" but assumed to be "medicaid"
                          (insurance == "uninsured")                        * (parameters[["util_uninsur"]]) +  
                          (fam_income == "poor")                                 * 0 +                # baseline = "Poor"
                          (fam_income == "near_poor")                            * (parameters[["util_inc_np"]]) + 
                          (fam_income == "low")                                  * (parameters[["util_inc_l"]] ) +
                          (fam_income == "medium")                               * (parameters[["util_inc_m"]] ) +
                          (fam_income == "high")                                 * (parameters[["util_inc_h"]] ) +
                          (education == "no_degree")                             * 0 +                # baseline = "no degree"
                          (education == "ged_hs")                                * (parameters[["util_ed_gh"]]) +
                          (education == "associate_bachelor")                    * (parameters[["util_ed_ab"]]) + 
                          (education == "master_doctorate")                      * (parameters[["util_ed_md"]]) +
                          (bmi_cat == "underweight")                             * 0 +                # baseline = "Underweight"
                          (bmi_cat == "normal_weight")                           * (parameters[["util_bmi_norm"]]) + 
                          (bmi_cat == "overweight")                              * (parameters[["util_bmi_ovw"]] ) +
                          (bmi_cat == "obese")                                   * (parameters[["util_bmi_obe"]] )  +
                          (cci == 0)                                             * 0 +                # baseline = 0
                          (cci == 1)                                             * (parameters[["util_cci1"]]) + 
                          (cci == 2)                                             * (parameters[["util_cci2"]]) +
                          (cci > 2)                                              * (parameters[["util_cci3"]]))
      ) %>% 
      
      # Keep only variables for transitions matrices
      select(utility)
    
    return(mydata) 
    
  }
  
  #### Weighted Standard Error function ----
  
  # https://www.alexstephenson.me/post/2022-04-02-weighted-variance-in-r/
  weighted.se.mean <- function(x, w, na.rm = T){
    ## Remove NAs 
    if (na.rm) {
      i <- !is.na(x)
      w <- w[i]
      x <- x[i]
    }
    
    ## Calculate effective N and correction factor
    n_eff <- (sum(w))^2/(sum(w^2))
    correction = n_eff/(n_eff-1)
    
    ## Get weighted variance 
    numerator = sum(w*(x-weighted.mean(x,w))^2)
    denominator = sum(w)
    
    ## get weighted standard error of the mean 
    se_x = sqrt((correction * (numerator/denominator))/n_eff)
    return(se_x)
  }
  
  #### Rubin's rules Standard Error function ----
  
  # Define a function to calculate the standard error using Rubin's rules
  rubin_se <- function(x, se) {
    se <- sqrt((sum(se^2) / length(x)) + ((1 + (1 / length(x))) * var(x)))
    return(se)
  }
  
  #### Transition probability rescaling function ----
  # Function to rescale non-missing values in each row to equal 1
  # This is necessary to rescale transition probabilities when their sum is greater than or equal to 1
  rescale_to_one <- function(row) {
    non_blank_values <- row[!is.na(row)]
    if ((length(non_blank_values) > 0) & (sum(non_blank_values) >= 1)) {
      row[!is.na(row)] <- non_blank_values / sum(non_blank_values)
    }
    return(row)
  }
  
  #### Transition probability residual function ----
  # Function to replace blank values with 1 minus the sum of non-blank values in each row
  replace_blanks <- function(row) {
    non_blank_values <- row[!is.na(row)]
    if (length(non_blank_values) > 0) {
      row[is.na(row)] <- 1 - sum(non_blank_values)
    }
    return(row)
  }
  
  #### Deterministic v Probabilistic values ----
  
  generate_parameters <- function(Deterministic, ci2se = 3.92) {
    parameter_list <- list()
    
    # Set-up intervention parameters
    variables <- c("int", 
                   "mi", "stroke", "hf", "cd", "ang", "pad", "diab", 
                   "a_25", "a_45", "a_65",
                   "fem", 
                   "race_w", "race_b", "race_a", "race_o", 
                   "medicare", "medicaid", "uninsur", 
                   "inc_np", "inc_l", "inc_m", "inc_h", 
                   "ed_gh", "ed_ab", "ed_md", 
                   "bmi_norm", "bmi_ovw", "bmi_obe",
                   "cci1", "cci2", "cci3")
    
    # predicting odds of non-zero costs
    for (variable in variables) {
      odds_nonzero_var <- if (Deterministic == 1) {
        odds_hc[[variable]]$mean
      } else {
        exp(rnorm(1, log(odds_hc[[variable]]$mean), ((log(odds_hc[[variable]]$uci) - log(odds_hc[[variable]]$lci)) / ci2se)))
      }
      
      parameter_list[[paste0("odds_nonzero_hc_", variable)]] <- odds_nonzero_var
    }
    
    # predicting non-zero costs
    for (variable in variables) {
      nonzero_hc_cost_var <- if (Deterministic == 1) {
        hc_cost[[variable]]$mean
      } else {
        exp(rnorm(1, log(hc_cost[[variable]]$mean), ((log(hc_cost[[variable]]$uci) - log(hc_cost[[variable]]$lci)) / ci2se)))
      }
      
      parameter_list[[paste0("nonzero_hc_cost_", variable)]] <- nonzero_hc_cost_var
    }
    
    # predicting utilities
    for (variable in variables) {
      util_var_value <- if (Deterministic == 1) {
        u[[variable]]$mean
      } else {
        (rnorm(1, (u[[variable]]$mean), (((u[[variable]]$uci) - (u[[variable]]$lci)) / ci2se)))
      }
      
      parameter_list[[paste0("util_", variable)]] <- util_var_value
    }
    
    # Additional parameters
    # government administration costs as part of treatment scenario for those receiving medicaid who are under 65 years (at age 65 years everyone switches to Medicare)
    parameter_list[["cost_admin_enrollee"]]    <- ifelse(Deterministic == 1, c_admin_enrollee$mean, rgamma(1, shape = c_admin_enrollee$shape, rate = c_admin_enrollee$rate))
    
    # Annual per person absenteeism costs 
    parameter_list[["c_absenteeism"]]          <- ifelse(Deterministic == 1, cost_absenteeism$mean, rgamma(1, shape = cost_absenteeism$shape, rate = cost_absenteeism$rate)) * pce_all_2013
    
    # Option to switch on and off cost of absenteeism
    parameter_list[["c_absenteeism"]]          <- ifelse(.c_absent == 1, parameter_list[["c_absenteeism"]], 0) 
    
    # Annual per person short-term disability costs 
    parameter_list[["c_disability"]]           <- ifelse(Deterministic == 1, cost_disability$mean, rgamma(1, shape = cost_disability$shape, rate = cost_disability$rate)) * pce_all_2013
    
    # Option to switch on and off cost of short-term disability
    parameter_list[["c_disability"]]           <- ifelse(.c_disab == 1, parameter_list[["c_disability"]] , 0)
    
    # change in total healthcare costs following receiving medicaid in medicaid expansion
    parameter_list[["delta_total_cost"]]       <- ifelse(Deterministic == 1, delta_c_total$mean, rbeta(1, delta_c_total$shape1, delta_c_total$shape2))

    # "history of CVD" multiplier (see Basu, 2017)
    parameter_list[["CVD_history_multiplier"]] <- ifelse(Deterministic == 1, delta_CVD_history$mean, rgamma(1, shape = delta_CVD_history$shape, scale = delta_CVD_history$scale))
    
    # reduction in systolic blood pressure after receiving medicaid [−3.03 mmHg; 95% CI, −5.33 mmHg to − 0.73 mmHg] from DOI: 10.1007/s11606-020-06417-6
    parameter_list[["sysbp_reduction"]]        <- ifelse(Deterministic == 1, delta_sys_bp$mean, rnorm(1, delta_sys_bp$mean, ((delta_sys_bp$uci - delta_sys_bp$lci) / ci2se)))
    
    # reduction in hba1c after receiving medicaid [−0.14 percentage points [pp]; 95% CI, −0.24 pp to −0.03 pp;] from DOI: 10.1007/s11606-020-06417-6
    parameter_list[["hba1c_reduction"]]        <- ifelse(Deterministic == 1, delta_hba1c$mean, rnorm(1, delta_hba1c$mean, ((delta_hba1c$uci - delta_hba1c$lci) / ci2se)))
    
    # change in total healthcare costs in year one of CVD event
    parameter_list[["delta_cost_y1"]]          <- ifelse(Deterministic == 1, delta_c_y1$mean, rbeta(1, delta_c_y1$shape1, delta_c_y1$shape2))
    
    # change in total healthcare costs in year one of CVD event
    parameter_list[["annual_earnings"]]        <- ifelse(.c_prem_mort == 1, annual_earnings, 0)
    
    return(parameter_list)
  }
  
  #### Microsimulation function ----
  # The MicroSim function for the simple microsimulation of the 'Sick-Sicker' model keeps track of what happens to each individual during each cycle. 
  # Adapted from https://doi.org/10.1177%2F0272989X18754513
  
  MicroSim <- function(original_data
                       , n.i
                       , n.t
                       , v.n
                       , d.c
                       , d.e
                       , TR.out = TRUE
                       , TS.out = TRUE
                       , Trt = FALSE
                       , parameters
                       ) {

    # update bootstrapped NHANES data in order to execute probs, cost and effects functions each cycle
    bootdata <- original_data %>%
      # create vector of starting states per individual according to their history of a diagnosis of any CVD conditions in bootstrapped NHANES data
      mutate(cvd              = ifelse((mi == 1 | stroke == 1 | chd == 1 | heart_failure == 1 | angina == 1),
                                       "History_CVD", "No_CVD"),  
             # duplicate age variable so that is can be updated with each cycle
             age_baseline     = age,
             # categorize age into specified categories each cycle that age is updated - Same categories used for cost and utility prediction as for weighting MI/stroke risk by ASCVD risk
             age_cat          = cut(age,                                 
                                    breaks = c(0, 24, 44, 64, Inf), 
                                    labels = c('<25', '25-44', '45-64', '65+')),
             mean_systolic_bp_orig = mean_systolic_bp,
             mean_systolic_bp = case_when(Trt == FALSE ~ mean_systolic_bp, 
                                          Trt == TRUE & medicaid_exp == 0 ~ mean_systolic_bp,
                                          Trt == TRUE & medicaid_exp == 1 & .bp == 0 ~ mean_systolic_bp, 
                                          # option to turn off effect of medicaid on bp in Trt scenario
                                          Trt == TRUE & medicaid_exp == 1 & .bp == 1 ~ (mean_systolic_bp - parameters[["sysbp_reduction"]]),
                                          TRUE ~ mean_systolic_bp), 
             # diagnosed or undiagnosed diabetes https://diabetes.org/diabetes/a1c/diagnosis
             diabetes_orig    = ifelse((diabetes_ever == 1 | (hba1c >= 6.5 & diabetes_ever == 0)), 1, 0),             
             hba1c = case_when(Trt == FALSE ~ hba1c, 
                               Trt == TRUE & medicaid_exp == 0 ~ hba1c,
                               # option to turn off effect of medicaid on hba1c in Trt scenario
                               Trt == TRUE & medicaid_exp == 1 & .hba1c == 0 ~ hba1c, 
                               # apply effect of medicaid on hba1c
                               Trt == TRUE & medicaid_exp == 1 & .hba1c == 1 ~ (hba1c - parameters[["hba1c_reduction"]]),
                               TRUE ~ hba1c),
             # diagnosed or undiagnosed diabetes https://diabetes.org/diabetes/a1c/diagnosis
             diabetes = ifelse((diabetes_ever == 1 | (hba1c >= 6.5 & diabetes_ever == 0)), 1, 0)
      )                                          

    # calculate the cost discount weight based on the discount rate d.c 
    v.dwc <- 1 / (1 + d.c) ^ (0:n.t)          
    # calculate the QALY discount weight based on the discount rate d.e 
    v.dwe <- 1 / (1 + d.e) ^ (0:n.t)                                              
    
    # Create the matrices capturing the state name/costs/health outcomes/ages for all individuals at each time point
    m.M <- m.C <- m.E <- m.age <- matrix(data = 0, nrow = n.i, ncol = n.t + 1, 
                                         dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                                         paste("cycle", 0:n.t, sep = " ")))
    
    # Collect costs across categories
    m.C_hc_cost <- m.C_delta_oop <- m.C_delta_nonoop <- m.C
    m.C_gov_admin <- m.C_disab <- m.C_absent <- m.C_lp <- m.C
    m.C_dist_gov <- m.C_dist_noop <- m.C

    # estimate costs per individual of the initial health state
    costs <- Costs(mydata  = bootdata, 
                   cardiac_dysrhythmia = 0,      
                   peripheral_artery_disease = 0,
                   cci     = 1,
                   Trt     = Trt,
                   parameters = parameters
                   , t = 0)
    
    # estimate effects/QALYs per individual of the initial health state
    effects <- Effs(mydata = bootdata, 
                    cardiac_dysrhythmia = 0,      
                    peripheral_artery_disease = 0,
                    cci    = 1,
                    Trt    = Trt,
                    parameters = parameters) 
    
    # indicate the initial health states as a function of the distribution of MI or stroke (or CHD or heart failure) history
    m.M[, 1] <- bootdata$cvd  
    # assign age at time t to the m.age matrix
    m.age[, 1] <- bootdata$age
    # assign effects/QALYs to the m.E matrix
    m.E[, 1] <- effects$utility    

    # assign cost-breakdown to other cost matrices
    m.C_hc_cost[, 1]      <- costs$hc_cost       
    m.C_delta_oop[, 1]    <- costs$delta_oop_hc_cost
    m.C_delta_nonoop[, 1] <- costs$delta_nonoop_hc_cost
    m.C_gov_admin[, 1]    <- costs$gov_admin_cost
    m.C_disab[, 1]        <- 0
    m.C_absent[, 1]       <- 0
    m.C_lp[, 1]           <- 0
    m.C_dist_gov[, 1]     <- costs$redist_gov_cost
    m.C_dist_noop[, 1]    <- costs$redist_noop_cost

    # assign initial costs to the m.C matrix
    m.C[, 1] <- (m.C_hc_cost[,1] + m.C_delta_oop[,1] +
                   m.C_dist_gov[, 1] + m.C_dist_noop[, 1] +
                   m.C_disab[,1] + m.C_absent[,1] + m.C_lp[,1])     
    
    # Beginning of time point loop
    for (t in 1:n.t) {
      
      bootdata <- bootdata %>%
        select(-age_cat) %>%
        # Update age to age + t with each cycle
        mutate(age            = age_baseline + t,  
               # categorize age into specified categories each cycle that age is updated - Same categories used for cost and utility prediction as for weighting MI/stroke risk by ASCVD risk
               age_cat        = cut(age,                                          
                                    breaks = c(0, 24, 44, 64, Inf), 
                                    labels = c('<25', '25-44', '45-64', '65+')),
               # Update both insurance variables used in the model such that all individual regardless of scenario receive Medicare at age 65
               insurance      = ifelse(age < 65, insurance,      "medicare"), 
               # no more effect of medicaid after age 65 (in particular on costs)
               medicaid_exp = ifelse(age < 65, medicaid_exp, 0),
               # diagnosed or undiagnosed diabetes https://diabetes.org/diabetes/a1c/diagnosis
               diabetes = ifelse((diabetes_ever == 1 | (hba1c >= 6.5 & diabetes_ever == 0)), 1, 0),
               # update bmi categories to correspond
               bmi_cat = case_when(bmi < 18.5               ~ "underweight", 
                                   (bmi >= 18.5 & bmi < 25) ~ "normal_weight",
                                   (bmi >= 25 & bmi < 30)   ~ "overweight",  
                                   bmi >= 30                ~ "obese",
                                   TRUE              ~ bmi_cat)
        )
      
      # Generate ASCVD weights each cycle using values from control group so that weights do not change between treatment scenarios [Otherwise artificially distorts results when BP is adjusted in treatment scenario]
      bootdata <- average_ascvd_risk(mydata = bootdata) 
      
      # Calculate transition probabilities for each cycle
      probs <- Probs(mydata = bootdata
                     , Trt    = Trt
                     , parameters = parameters
                     )
      
      # create vector of state transition probabilities as a function of an individual's characteristics
      # probability of dying from a MI
      p_MI_CVDdeath       <- matrix(probs$mi_mortality_risk,     1, n.i) 
      # probability of dying from a stroke
      p_stroke_CVDdeath   <- matrix(probs$stroke_mortality_risk, 1, n.i) 
      
      # probability of remaining in the dead from other causes state
      p_death_death       <- matrix(1, 1, n.i) 
      # probability of remaining in the dead from CVD state
      p_CVDdeath_CVDdeath <- matrix(1, 1, n.i)                                
      
      # probability of dying from other causes for those with a history of CVD
      p_death             <- matrix(probs$noncvd_mortality_risk, 1, n.i)
      # probability of having MI with a history of CVD
      p_CVD_MI            <- matrix(probs$mi_risk_CVD,           1, n.i)
      # probability of having stroke with a history of CVD
      p_CVD_stroke        <- matrix(probs$stroke_risk_CVD,       1, n.i)
      # probability of having a stroke with no history of CVD
      p_nCVD_MI           <- matrix(probs$mi_risk_noCVD,         1, n.i)
      # probability of having MI with a history of CVD
      p_nCVD_stroke       <- matrix(probs$stroke_risk_noCVD,     1, n.i)      
      
      # create a transition matrix as a function of an individual's characteristics for those in the "No_CVD" state
      m.p.it_nCVD <- matrix(NA, n.s, n.i)   
      # assign names to the vector
      rownames(m.p.it_nCVD) <- v.n                                                
      m.p.it_nCVD[1,] <-  NA # Residual; see function below
      m.p.it_nCVD[2,] <-  p_nCVD_MI * (1-p_death)
      m.p.it_nCVD[3,] <-  p_nCVD_stroke * (1-p_death)
      m.p.it_nCVD[4,] <-  0
      m.p.it_nCVD[5,] <-  0
      m.p.it_nCVD[6,] <-  p_death
      m.p.it_nCVD     <-  t(m.p.it_nCVD)
      
      # create a transition matrix as a function of an individual's characteristics for those in the "MI" state
      m.p.it_MI <- matrix(NA, n.s, n.i) 
      rownames(m.p.it_MI) <- v.n            
      m.p.it_MI[1,] <-  0
      m.p.it_MI[2,] <-  0
      m.p.it_MI[3,] <-  0
      m.p.it_MI[4,] <-  NA # Residual; see function below
      m.p.it_MI[5,] <-  p_MI_CVDdeath * (1-p_death)
      m.p.it_MI[6,] <-  p_death
      m.p.it_MI     <-  t(m.p.it_MI)
      
      # create a transition matrix as a function of an individual's characteristics for those in the "Stroke" state
      m.p.it_stroke <- matrix(NA, n.s, n.i)     
      rownames(m.p.it_stroke) <- v.n          
      m.p.it_stroke[1,] <-  0
      m.p.it_stroke[2,] <-  0
      m.p.it_stroke[3,] <-  0
      m.p.it_stroke[4,] <-  NA # Residual; see function below
      m.p.it_stroke[5,] <-  p_stroke_CVDdeath * (1-p_death)
      m.p.it_stroke[6,] <-  p_death
      m.p.it_stroke     <-  t(m.p.it_stroke)
      
      # create a transition matrix as a function of an individual's characteristics for those in the "History_CVD" state
      m.p.it_CVD <- matrix(NA, n.s, n.i)     
      rownames(m.p.it_CVD) <- v.n            
      m.p.it_CVD[1,] <-  0
      m.p.it_CVD[2,] <-  p_CVD_MI * (1-p_death)
      m.p.it_CVD[3,] <-  p_CVD_stroke * (1-p_death) 
      m.p.it_CVD[4,] <-  NA # Residual; see function below
      m.p.it_CVD[5,] <-  0
      m.p.it_CVD[6,] <-  p_death
      m.p.it_CVD     <-  t(m.p.it_CVD)
      
      # create a transition matrix as a function of an individual's characteristics for those in the "CVD_death" state
      m.p.it_CVD_death <- matrix(NA, n.s, n.i)     
      rownames(m.p.it_CVD_death) <- v.n            
      m.p.it_CVD_death[1,] <-  0
      m.p.it_CVD_death[2,] <-  0
      m.p.it_CVD_death[3,] <-  0
      m.p.it_CVD_death[4,] <-  0
      m.p.it_CVD_death[5,] <-  p_CVDdeath_CVDdeath
      m.p.it_CVD_death[6,] <-  0
      m.p.it_CVD_death     <-  t(m.p.it_CVD_death)
      
      # create a transition matrix as a function of an individual's characteristics for those in the "Non-CVD_death" state
      m.p.it_death <- matrix(NA, n.s, n.i)     
      rownames(m.p.it_death) <- v.n            
      m.p.it_death[1,] <-  0
      m.p.it_death[2,] <-  0
      m.p.it_death[3,] <-  0
      m.p.it_death[4,] <-  0
      m.p.it_death[5,] <-  0
      m.p.it_death[6,] <-  p_death_death
      m.p.it_death     <-  t(m.p.it_death)
      
      # Apply functions to each row to generate the residual transition probability after rescaling transition probabilities whose sum is greater than 1
      for (i in 1:n.i) {
        # Rescale probabilities
        m.p.it_nCVD[i, ]   <- rescale_to_one(m.p.it_nCVD[i, ])
        m.p.it_MI[i, ]     <- rescale_to_one(m.p.it_MI[i, ])
        m.p.it_stroke[i, ] <- rescale_to_one(m.p.it_stroke[i, ])
        m.p.it_CVD[i, ]    <- rescale_to_one(m.p.it_CVD[i, ])
      }
      
      for (i in 1:n.i) {
        # Estimate residual
        m.p.it_nCVD[i, ]   <- replace_blanks(m.p.it_nCVD[i, ]) 
        m.p.it_MI[i, ]     <- replace_blanks(m.p.it_MI[i, ])     
        m.p.it_stroke[i, ] <- replace_blanks(m.p.it_stroke[i, ]) 
        m.p.it_CVD[i, ]    <- replace_blanks(m.p.it_CVD[i, ])
      }
      
      # remove unused variables
      rm(p_MI_CVDdeath,
         p_stroke_CVDdeath,
         p_death_death,
         p_CVDdeath_CVDdeath,
         p_CVD_MI,           
         p_CVD_stroke,       
         p_death,       
         p_nCVD_MI,          
         p_nCVD_stroke)
      
      # Create empty matrix to fill with state transition probabilities
      m.p.it <- matrix(NA, nrow = n.i, ncol = n.s) 
      
      m.p.it <- case_when(
        m.M[, t] == "No_CVD" ~ m.p.it_nCVD,
        m.M[, t] == "MI" ~ m.p.it_MI,
        m.M[, t] == "Stroke" ~ m.p.it_stroke,
        m.M[, t] == "History_CVD" ~ m.p.it_CVD,
        m.M[, t] == "CVD_death" ~ m.p.it_CVD_death,
        m.M[, t] == "Non-CVD_death" ~ m.p.it_death,
        # If none of the conditions match, keep the original m.p.it
        TRUE ~ m.p.it  
      )
      
      rm(m.p.it_nCVD,
         m.p.it_MI,
         m.p.it_stroke,
         m.p.it_CVD,
         m.p.it_CVD_death,
         m.p.it_death)      
      
      if (.track_prob == 1) {
        # print the transition probabilities or produce an error if they don't sum to 1
        ifelse((rowSums(m.p.it) <= 1.0000009) & (rowSums(m.p.it) >= 0.9999990), print("Probabilities sum to 1"), print("ERROR: Probabilities DO NOT sum to 1")) 
      }
      
      # assign names to the vector
      colnames(m.p.it) <- v.n                                                     

      # sample the next health state and store that state in matrix m.M
      m.M[, t + 1] <- rcat(n.i, m.p.it)                                          
      # Mapping values to characters
      m.M[, t + 1] <- v.n[match(m.M[, t + 1], 1:6)]      
      
      # Update bootdata to reflect CVD history
      bootdata <- bootdata %>%
        # Update the MI variable which calculate costs & utilities to reflect whether an individual has had a diagnosis of an MI
        mutate(mi      = ifelse(m.M[, t + 1] == "MI",      
                                1, mi),    
               # Update the Stroke variable which calculate costs & utilities to reflect whether an individual has had a diagnosis of an Stroke
               stroke  = ifelse(m.M[, t + 1] == "Stroke",  
                                1, stroke))

      # estimate costs per individual during cycle t + 1
      costs <- Costs(mydata                    = bootdata,
                     cardiac_dysrhythmia       = 0,      
                     peripheral_artery_disease = 0,
                     cci                       = 1,
                     Trt                       = Trt,
                     parameters                = parameters
                     , t = t)
      
      # estimate effects/QALYs per individual during cycle t + 1
      effects <- Effs(mydata                    = bootdata,
                      cardiac_dysrhythmia       = 0,      
                      peripheral_artery_disease = 0,
                      cci                       = 1,
                      Trt                       = Trt,
                      parameters                = parameters) 

      m.age[, t + 1] <- bootdata$age
      # assign costs during cycle t + 1 to the m.C_oop matrix
      m.C_hc_cost[, t + 1] <- case_when((m.M[, t + 1] == "No_CVD") ~ costs$hc_cost, 
                                (m.M[, t + 1] == "Stroke")         ~ (costs$hc_cost * (1 + parameters[["delta_cost_y1"]])), 
                                (m.M[, t + 1] == "MI")             ~ (costs$hc_cost * (1 + parameters[["delta_cost_y1"]])),
                                (m.M[, t + 1] == "History_CVD")    ~ costs$hc_cost,
                                (m.M[, t + 1] == "CVD_death")      ~ 0,
                                (m.M[, t + 1] == "Non-CVD_death")  ~ 0)
      # assign reallocated societal costs during cycle t + 1 to the m.C_dist_gov matrix
      m.C_dist_gov[, t + 1] <- case_when((m.M[, t + 1] == "No_CVD")    ~ costs$redist_gov_cost, 
                                    (m.M[, t + 1] == "Stroke")         ~ costs$redist_gov_cost, 
                                    (m.M[, t + 1] == "MI")             ~ costs$redist_gov_cost,
                                    (m.M[, t + 1] == "History_CVD")    ~ costs$redist_gov_cost,
                                    (m.M[, t + 1] == "CVD_death")      ~ 0,
                                    (m.M[, t + 1] == "Non-CVD_death")  ~ 0)
      # assign reallocated societal costs during cycle t + 1 to the m.C_dist_noop matrix
      m.C_dist_noop[, t + 1] <- case_when((m.M[, t + 1] == "No_CVD")   ~ costs$redist_noop_cost, 
                                    (m.M[, t + 1] == "Stroke")         ~ costs$redist_noop_cost, 
                                    (m.M[, t + 1] == "MI")             ~ costs$redist_noop_cost,
                                    (m.M[, t + 1] == "History_CVD")    ~ costs$redist_noop_cost,
                                    (m.M[, t + 1] == "CVD_death")      ~ 0,
                                    (m.M[, t + 1] == "Non-CVD_death")  ~ 0)

      # assign change in OOP costs during cycle t + 1 to the m.C_oop matrix
      m.C_delta_oop[, t + 1] <- case_when((m.M[, t + 1] == "No_CVD")  ~ costs$delta_oop_hc_cost, 
                                    (m.M[, t + 1] == "Stroke")        ~ costs$delta_oop_hc_cost, 
                                    (m.M[, t + 1] == "MI")            ~ costs$delta_oop_hc_cost,
                                    (m.M[, t + 1] == "History_CVD")   ~ costs$delta_oop_hc_cost,
                                    (m.M[, t + 1] == "CVD_death")     ~ 0,
                                    (m.M[, t + 1] == "Non-CVD_death") ~ 0)

      # assign change in non-OOP costs during cycle t + 1 to the m.C_delta_nonoop matrix
      m.C_delta_nonoop[, t + 1] <- case_when((m.M[, t + 1] == "No_CVD")     ~ costs$delta_nonoop_hc_cost, 
                                          (m.M[, t + 1] == "Stroke")        ~ costs$delta_nonoop_hc_cost, 
                                          (m.M[, t + 1] == "MI")            ~ costs$delta_nonoop_hc_cost,
                                          (m.M[, t + 1] == "History_CVD")   ~ costs$delta_nonoop_hc_cost,
                                          (m.M[, t + 1] == "CVD_death")     ~ 0,
                                          (m.M[, t + 1] == "Non-CVD_death") ~ 0)

      # assign governemnt administration costs during cycle t + 1 to the m.C_gov_admin matrix
      m.C_gov_admin[, t + 1] <- case_when((m.M[, t + 1] == "No_CVD") ~ costs$gov_admin_cost, 
                                (m.M[, t + 1] == "Stroke")           ~ costs$gov_admin_cost, 
                                (m.M[, t + 1] == "MI")               ~ costs$gov_admin_cost, 
                                (m.M[, t + 1] == "History_CVD")      ~ costs$gov_admin_cost, 
                                (m.M[, t + 1] == "CVD_death")        ~ 0,
                                (m.M[, t + 1] == "Non-CVD_death")    ~ 0)
      # assign short-term disability costs during cycle t + 1 to the m.C_disab matrix
      m.C_disab[, t + 1] <- case_when((m.M[, t + 1] == "No_CVD")  ~ 0,
                                (m.M[, t + 1] == "Stroke")         ~ parameters[["c_disability"]]*(ifelse(bootdata$age <=  65, 1, 0)),
                                (m.M[, t + 1] == "MI")             ~ parameters[["c_disability"]]*(ifelse(bootdata$age <=  65, 1, 0)), 
                                (m.M[, t + 1] == "History_CVD")   ~ 0,
                                (m.M[, t + 1] == "CVD_death")     ~ 0,
                                (m.M[, t + 1] == "Non-CVD_death") ~ 0 )
      # assign short-term disability costs during cycle t + 1 to the m.C_absent matrix
      m.C_absent[, t + 1] <- case_when((m.M[, t + 1] == "No_CVD")       ~ 0,
                                       (m.M[, t + 1] == "Stroke")        ~ parameters[["c_absenteeism"]]*(ifelse(bootdata$age <=  65, 1, 0)),
                                       (m.M[, t + 1] == "MI")            ~ parameters[["c_absenteeism"]]*(ifelse(bootdata$age <=  65, 1, 0)), 
                                      (m.M[, t + 1] == "History_CVD")   ~ 0,
                                      (m.M[, t + 1] == "CVD_death")     ~ 0,
                                      (m.M[, t + 1] == "Non-CVD_death") ~ 0 )
      # assign lost earnings from a premature death (i.e. prior to retirement) during cycle t + 1 to the m.C_lp matrix
      m.C_lp[, t + 1] <- case_when((m.M[, t + 1] == "No_CVD")   ~ 0,
                                (m.M[, t + 1] == "Stroke")      ~ 0,
                                (m.M[, t + 1] == "MI")          ~ 0,
                                (m.M[, t + 1] == "History_CVD") ~ 0,
                                (m.M[, t + 1] == "CVD_death")    ~ parameters[["annual_earnings"]]*(ifelse(bootdata$age <=  65, 1, 0)),
                                (m.M[, t + 1] == "Non-CVD_death")~ parameters[["annual_earnings"]]*(ifelse(bootdata$age <=  65, 1, 0)) 
                                )
      
      # assign total costs during cycle t + 1 to the m.C matrix
      m.C[, t + 1] <- (m.C_hc_cost[,t+1] + m.C_delta_oop[,t+1] +
                         m.C_dist_gov[, t + 1] + m.C_dist_noop[, t + 1] +
                         m.C_disab[,t+1] + m.C_absent[,t+1] + m.C_lp[,t+1])
      
      # assign effects during cycle t + 1 to the m.E matrix
      m.E[, t + 1] <- case_when((m.M[, t + 1] == "No_CVD")        ~ effects$utility,
                                (m.M[, t + 1] == "Stroke")        ~ effects$utility, 
                                (m.M[, t + 1] == "MI")            ~ effects$utility, 
                                (m.M[, t + 1] == "History_CVD")   ~ effects$utility,
                                (m.M[, t + 1] == "CVD_death")     ~ 0,
                                (m.M[, t + 1] == "Non-CVD_death") ~ 0)
      
      # Option to switch on or off the half-cycle correction
      if ( .half == 1) {
        # apply half-cycle corrections across cycles
        m.C[, t] <- (m.C[, t + 1] + m.C[, t]) * half_cycle_correction   
        
        m.C_hc_cost     [, t] <- (m.C_hc_cost     [, t + 1] + m.C_hc_cost     [, t]) * half_cycle_correction   
        m.C_delta_oop   [, t] <- (m.C_delta_oop   [, t + 1] + m.C_delta_oop   [, t]) * half_cycle_correction   
        m.C_delta_nonoop[, t] <- (m.C_delta_nonoop[, t + 1] + m.C_delta_nonoop[, t]) * half_cycle_correction   
        m.C_gov_admin   [, t] <- (m.C_gov_admin   [, t + 1] + m.C_gov_admin   [, t]) * half_cycle_correction   
        m.C_disab       [, t] <- (m.C_disab       [, t + 1] + m.C_disab       [, t]) * half_cycle_correction   
        m.C_absent      [, t] <- (m.C_absent      [, t + 1] + m.C_absent      [, t]) * half_cycle_correction   
        m.C_lp          [, t] <- (m.C_lp          [, t + 1] + m.C_lp          [, t]) * half_cycle_correction   
        m.C_dist_gov    [, t] <- (m.C_dist_gov    [, t + 1] + m.C_dist_gov    [, t]) * half_cycle_correction   
        m.C_dist_noop   [, t] <- (m.C_dist_noop   [, t + 1] + m.C_dist_noop   [, t]) * half_cycle_correction   
        
        m.E[, t] <- (m.E[, t + 1] + m.E[, t]) * half_cycle_correction         
      }
      

      # display the progress of the simulation
      if (.track_cycle == 1) {
        cat('\r', paste(round(t/n.t * 100), "% done", sep = " "))                   
      }
      
    } # close the loop for the time points 
    
    # Remove values above trunc_age (default = 85) is missing
    m.M[m.age >= .trunc_age] <- paste("Age_", .trunc_age, sep = "")
    m.C[m.age >= .trunc_age] <- 0
    m.E[m.age >= .trunc_age] <- 0

    m.C_hc_cost     [m.age >= .trunc_age]<- 0
    m.C_delta_oop   [m.age >= .trunc_age]<- 0
    m.C_delta_nonoop[m.age >= .trunc_age]<- 0
    m.C_gov_admin   [m.age >= .trunc_age]<- 0
    m.C_disab       [m.age >= .trunc_age]<- 0
    m.C_absent      [m.age >= .trunc_age]<- 0
    m.C_lp          [m.age >= .trunc_age]<- 0
    m.C_dist_gov    [m.age >= .trunc_age]<- 0
    m.C_dist_noop   [m.age >= .trunc_age]<- 0

    # Calculate the number of person-years that accumulated for non-elderly adults 
    # Create a subset of v.n containing only the first four values
    if (.excl_death == 1) {
      v.n_subset <- v.n[1:4]
    } else {
      v.n_subset <- v.n[1:6]
    }
    
    # Count occurrences for each row
    person_years <- apply(m.M, MARGIN = 1, FUN = function(row) sum(row %in% v.n_subset))
    
    # total discounted cost per individual
    tc            <- m.C              %*% v.dwc
    hc            <- m.C_hc_cost      %*% v.dwc
    delta_oopc    <- m.C_delta_oop    %*% v.dwc
    delta_nonoopc <- m.C_delta_nonoop %*% v.dwc
    govc          <- m.C_gov_admin    %*% v.dwc
    disabc        <- m.C_disab        %*% v.dwc
    absentc       <- m.C_absent       %*% v.dwc
    lpc           <- m.C_lp           %*% v.dwc 
    r_gov         <- m.C_dist_gov     %*% v.dwc 
    r_noop        <- m.C_dist_noop    %*% v.dwc 
    
    # total discounted QALYs per individual
    te <- m.E %*% v.dwe      
      
    # total number of life-years (before truncated age)
    ly <- apply((matrix(person_years)), 1, sum)
    # Count number of outcomes per person
    mi_count <- apply(m.M == "MI", 1, sum)
    stroke_count <- apply(m.M == "Stroke", 1, sum)
    CVD_death_count <- apply(m.M == "CVD_death", 1
                             , function(row) if (sum(row, na.rm = TRUE) >= 1) 1 else 0)
    NonCVD_death_count <- apply(m.M == "Non-CVD_death", 1
                                , function(row) if (sum(row, na.rm = TRUE) >= 1) 1 else 0)

    # create a  matrix of transitions across states
    if (TS.out == TRUE) {  
      # transitions from one state to the other
      TS <- paste(m.M, cbind(m.M[, -1], NA), sep = "->")                         
      TS <- matrix(TS, nrow = n.i)     
      # name the rows 
      rownames(TS) <- paste("Ind",   1:n.i, sep = " ")  
      # name the columns
      colnames(TS) <- paste("Cycle", 0:n.t, sep = " ")                            
    } else {
      TS <- NULL
    }
    
    # create a trace from the individual trajectories
    if (TR.out == TRUE) { 
      TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE))))
      # create a distribution trace
      TR <- TR / n.i  
      # name the rows 
      rownames(TR) <- paste("Cycle", 0:n.t, sep = " ")  
      # name the columns 
      colnames(TR) <- v.n                                                         
    } else {
      TR <- NULL
    }

    results <- list(
       tc = tc
      , te = te
      , ly = ly
      , hc            = hc           
      , delta_oopc    = delta_oopc   
      , delta_nonoopc = delta_nonoopc
      , govc          = govc         
      , disabc        = disabc       
      , absentc       = absentc      
      , lpc           = lpc 
      , r_gov  = r_gov 
      , r_noop = r_noop
      , person_years = person_years
      , mi_count = mi_count
      , stroke_count = stroke_count
      , CVD_death_count = CVD_death_count
      , NonCVD_death_count = NonCVD_death_count
      , TS = TS
      , TR = TR
      , bootdata = bootdata
    )

    return(results)     
    
  }  # end of the MicroSim function  
  
  #### Microsimulation bootstrapping function ----
  
  # The Bootstrapping function which runs all scenarios across same randomly sampled NHANES dataset
  Bootstrap_msm <- function(original_data 
                            , n.i
                            ) { 
    
    # Do this before Micro-simulation to have the same sample for treatment and control scenarios
    bootdata_original <- original_data %>%
      # categorize age into specified categories each cycle that age is updated - Same categories used for cost and utility prediction as for weighting MI/stroke risk by ASCVD risk
      mutate(age_cat = cut(age,                                            
                           breaks = c(0, 24, 44, 64, Inf), 
                           labels = c('<25', '25-44', '45-64', '65+')),
             # return variables to numeric for operation (following multiple imputation)
             stroke        = as.numeric(stroke),
             mi            = as.numeric(mi),
             heart_failure = as.numeric(heart_failure),
             chd           = as.numeric(chd),
             angina        = as.numeric(angina),
             diabetes_ever = as.numeric(diabetes_ever),
             treated       = as.numeric(treated),
             smoking       = as.numeric(smoking),
             # Create variable which identifies individuals who will receive Medicaid as part of the treatment scenario (under age 65 and family income less than 138% of the federal poverty line)
             medicaid_exp   = ifelse(insurance  == "uninsured" & fpl < .fpl & age < 65, 1, 0),
             medicaid_exp_baseline = ifelse(insurance  == "uninsured" & fpl < .fpl & age < 65, 1, 0),
             # create a family income variable to correspond to MEPS cost and utility calculations and according to the poverty line [see, https://meps.ahrq.gov/survey_comp/hc_technical_notes.shtml]
             # created within function because it follows imputation
             fam_income = case_when(fpl <= 1     ~ "poor",                          
                                    fpl <= 1.25  ~ "near_poor",
                                    fpl <= 2     ~ "low",
                                    fpl <= 4     ~ "medium",
                                    fpl  > 4     ~ "high",
                                    TRUE         ~ NA),
             # same for bmi categories as family income
             bmi_cat = case_when(bmi < 18.5               ~ "underweight", 
                                 (bmi >= 18.5 & bmi < 25) ~ "normal_weight",
                                 (bmi >= 25 & bmi < 30)   ~ "overweight",  
                                 bmi >= 30                ~ "obese",
                                 TRUE              ~ NA) %>%
               set_variable_labels(fam_income = "Family Income")
      )
    
    # Filter the sample to run the model for certain populations
    bootdata_original <- bootdata_original %>%
      filter(# filter to just those receiving the treatment if .ATT == 1, otherwise the whole sample 
        medicaid_exp == 1 | .ATT != 1)

    # sample NHANES data with replacement using survey weights created across waves for fasting subsample
    bootdata_original <- slice_sample(bootdata_original, n = n.i, replace = TRUE, weight_by = WTSAF)
    
    # Create a list of (Deterministic or Probabilistic) parameters to be used in the simulation
    parameters <- generate_parameters(Deterministic = .Deterministic, ci2se = ci2se)
    
    #### Run the simulation ----
    #running the function several times to obtain CI
    sim_no_trt  <- MicroSim(bootdata_original
                            , n.i
                            , n.t
                            , v.n
                            , d.c
                            , d.e
                            , TR.out = TRUE
                            , TS.out = TRUE
                            , Trt = FALSE
                            , parameters
                            ) # run for no treatment
    
    sim_trt     <- MicroSim(bootdata_original
                            , n.i
                            , n.t
                            , v.n
                            , d.c
                            , d.e
                            , TR.out = TRUE
                            , TS.out = TRUE
                            , Trt = TRUE
                            , parameters
                            )  # run for treatment
    
    #### Collate Data ----
    
    # Define the list of variables to be divided by person_years_ntrt
    variables_to_divide_ntrt <- c("mi_ntrt", "stroke_ntrt", "CVD_death_ntrt", "nCVD_death_ntrt")
    variables_to_divide_trt <- c("mi_trt", "stroke_trt", "CVD_death_trt", "nCVD_death_trt")

    combined_data <- data.frame(
      medicaid_exp        = c(sim_no_trt$bootdata$medicaid_exp_baseline)
      , age_baseline      = c(sim_no_trt$bootdata$age_baseline)
      , female            = c(sim_no_trt$bootdata$female)
      , race              = c(sim_no_trt$bootdata$race)
      , insurance         = c(sim_no_trt$bootdata$insurance)
      , education         = c(sim_no_trt$bootdata$education)
      , fam_income        = c(sim_no_trt$bootdata$fam_income)
      , mi_ntrt           = c(sim_no_trt$mi_count)
      , stroke_ntrt       = c(sim_no_trt$stroke_count)
      , CVD_death_ntrt    = c(sim_no_trt$CVD_death_count)
      , nCVD_death_ntrt   = c(sim_no_trt$NonCVD_death_count)
      , person_years_ntrt = c(sim_no_trt$person_years)
      , mi_trt            = c(sim_trt$mi_count)
      , stroke_trt        = c(sim_trt$stroke_count)
      , CVD_death_trt     = c(sim_trt$CVD_death_count)
      , nCVD_death_trt    = c(sim_trt$NonCVD_death_count)
      , person_years_trt  = c(sim_trt$person_years)
      , tc_ntrt           = c(sim_no_trt$tc)
      , tc_trt            = c(sim_trt$tc)
      , te_ntrt           = c(sim_no_trt$te)
      , te_trt            = c(sim_trt$te)
      , ly_ntrt           = c(sim_no_trt$ly)
      , ly_trt            = c(sim_trt$ly)
      , hc_trt            = c(sim_trt$hc)
      , hc_ntrt           = c(sim_no_trt$hc)
      , d_oopc_trt        = c(sim_trt$delta_oopc)
      , d_oopc_ntrt       = c(sim_no_trt$delta_oopc)
      , d_noopc_trt       = c(sim_trt$delta_nonoopc)
      , d_noopc_ntrt      = c(sim_no_trt$delta_nonoopc)
      , govc_trt          = c(sim_trt$govc)
      , govc_ntrt         = c(sim_no_trt$govc)
      , disabc_trt        = c(sim_trt$disabc)
      , disabc_ntrt       = c(sim_no_trt$disabc)
      , absentc_trt       = c(sim_trt$absentc)
      , absentc_ntrt      = c(sim_no_trt$absentc)
      , lpc_trt           = c(sim_trt$lpc)
      , lpc_ntrt          = c(sim_no_trt$lpc)
      , r_gov_trt         = c(sim_trt$r_gov)
      , r_gov_ntrt        = c(sim_no_trt$r_gov)
      , r_noop_trt        = c(sim_trt$r_noop)
      , r_noop_ntrt       = c(sim_no_trt$r_noop)
    ) %>%
      mutate(
        age_cat = cut(age_baseline, breaks = c(0, 24, 44, 64, Inf), labels = c('<25', '25-44', '45-64', '65+'))
      )

    if (.raw_data == 1) {
      raw_data <- c(sim_no_trt, sim_trt)
    } else {
      raw_data <- list()
    }
    
    if(.rate == 1) {
      combined_data <- combined_data %>%
        mutate(
          across(all_of(variables_to_divide_ntrt), ~ . / person_years_ntrt),
          across(all_of(variables_to_divide_trt), ~ . / person_years_trt),
        )
    }
    
    # Initialize an empty list to store results for each variable
    summary_list <- list()
    
    # Define a list of variables for which you want to calculate means
    variables_to_mean <- c("tc_ntrt", "tc_trt", "hc_trt", "hc_ntrt"
                           , "d_oopc_trt", "d_oopc_ntrt"
                           , "d_noopc_trt", "d_noopc_ntrt"
                           , "govc_trt", "govc_ntrt"
                           , "disabc_trt", "disabc_ntrt", "absentc_trt", "absentc_ntrt"
                           , "lpc_trt", "lpc_ntrt"
                           , "r_gov_trt", "r_gov_ntrt", "r_noop_trt", "r_noop_ntrt"
                           , "te_ntrt", "te_trt"
                           , "mi_ntrt", "stroke_ntrt", "CVD_death_ntrt", "nCVD_death_ntrt"
                           , "mi_trt", "stroke_trt", "CVD_death_trt", "nCVD_death_trt"
                           , "ly_ntrt", "ly_trt", "person_years_ntrt", "person_years_trt")
    
    # Loop through each variable
    for (var in .grouping_vars) {
      # Get unique categories within the current variable
      categories <- unique(combined_data %>% pull({{ var }}))
      
      # Initialize an empty list to store results for each category
      category_summary_list <- list()
      
      for (cat in categories) {
        # Calculate the group summary for the current category
        group_summary <- combined_data %>%
          filter(.data[[var]] == cat) %>%
          summarise(
            variable          = var
            , category        = cat
            , num_individuals = n()
            , across(all_of(variables_to_mean), ~ mean(., na.rm = TRUE))
          ) 
        
        # Append the summary for the current category to the list
        category_summary_list[[as.character(cat)]] <- group_summary
      }
      
      # Combine category summaries for the current variable
      var_summary <- do.call(rbind, category_summary_list)
      
      # Append the summary for the current variable to the main list
      summary_list[[var]] <- var_summary
    }
    
    # Combine all variable summaries into a single dataframe
    result <- do.call(rbind, summary_list)
    
    # Calculate the total summary for all categories within each variable
    total_summary <- combined_data %>%
      summarize(
        variable            = "Total"
        , category          = "Total"
        , num_individuals   = n()
        , across(all_of(variables_to_mean), ~ mean(., na.rm = TRUE))
      )
    
    # Append the total summary to the result
    combined_output <- bind_rows(result, total_summary)
    
    # relabel category for female/male
    combined_output <- combined_output %>%
      mutate(category = case_when(category == 1 ~ "Medicaid eligible",
                                  category == 0 ~ "Not Medicaid eligible",
                                  TRUE ~ as.character(category)))   
    
    msm_results <- list(
      raw_data
      , combined_output
    )
    
    return(msm_results)
    
  }
  
  #### Deterministic Sensitivity Analysis ----
  
  # Function to perform Bootstrap_msm with perturbed values
  dsa <- function(perturbed_values) {
    # Copy the perturbed values to avoid overwriting base_values
    modified_values    <- base_values
    modified_values[i] <- perturbed_values[i]
    
    delta_CVD_history$mean <<- perturbed_values[1]
    delta_hba1c$mean       <<- perturbed_values[2]
    delta_sys_bp$mean      <<- perturbed_values[3]
    odds_hc$mi$mean        <<- perturbed_values[4]
    odds_hc$stroke$mean    <<- perturbed_values[5]
    hc_cost$mi$mean        <<- perturbed_values[6]
    hc_cost$stroke$mean    <<- perturbed_values[7]
    cost_absenteeism$mean  <<- perturbed_values[8]
    cost_disability$mean   <<- perturbed_values[9]
    delta_c_total$mean     <<- perturbed_values[10]
    delta_c_y1$mean        <<- perturbed_values[11]
    c_admin_enrollee$mean  <<- perturbed_values[12]
    annual_earnings        <<- perturbed_values[13]
    u$mi$mean              <<- perturbed_values[14]
    u$stroke$mean          <<- perturbed_values[15]
    d.e                    <<- perturbed_values[16]
    d.c                    <<- perturbed_values[17]
    
    #### Bootstrap the simulation ----
    
    # Initialize an empty dataframe for the combined output
    dsa_combined_output <- data.frame()
    
    for (dataset_index in 1:m) {
      # Get the dataset for the current iteration
      current_data <- dataList[[paste("Dataset_", dataset_index, sep = "")]]
      
      for (replication in 1:bootsize) {
        # Run the Bootstrap_msm function for the current replication on the current dataset
        set.seed(seed)
        dsa_result <- Bootstrap_msm(original_data = current_data,
                                    n.i = n.i)
        
        # Add replication and dataset index variables to the result data frame
        dsa_result[[2]]$replication   <- replication
        dsa_result[[2]]$dataset_index <- dataset_index
        
        # Append the results to the combined data frame
        dsa_combined_output <- bind_rows(dsa_combined_output, dsa_result[[2]])
        
      }
      
      # Display the progress of the simulation
      if (.track_boot == 1) {
        cat('\r', paste(round(replication/bootsize * 100), "% Bootstrap done for Dataset_", dataset_index, sep = " "))
      }
    }
    
    #### Summarize the results  ----
    
    # Define a list of variables for which you want to calculate means
    variables_to_operate <- c("tc_ntrt", "tc_trt", 
                              "te_ntrt", "te_trt")
    
    # Initialize an empty list to store summarised outputs
    dsa_impute_outputs <- list()
    
    for (imputation_index in 1:m) {
      # Create the summarised_output for the current dataset
      dsa_summarised_output <- dsa_combined_output %>%
        filter(imputation_index == .data$dataset_index,
               variable == "Total") %>%  # Filter data for the current dataset
        summarize(
          across(all_of(variables_to_operate), ~ (sum(. * num_individuals) / sum(num_individuals)))
          , n = sum(num_individuals)
        ) %>%
        ungroup()
      
      # Store the summarised_output in the list
      dsa_impute_outputs[[paste("summarised_output_", imputation_index, sep = "")]] <- dsa_summarised_output
    }
    
    # Combine all the summarised outputs into one dataset
    dsa_combined_summarised_output <- bind_rows(dsa_impute_outputs) 
    
    # Calculate the mean for selected variables
    dsa_summarised_output <- dsa_combined_summarised_output %>%
      summarize(# Use Rubins rules to combine se's across imputation results
        tc_ntrt = sum(tc_ntrt * n) / sum(n)
        , tc_trt = sum(tc_trt * n) / sum(n)
        , te_ntrt = sum(te_ntrt * n) / sum(n)
        , te_trt = sum(te_trt * n) / sum(n)
        , num_individuals = sum(n)
      ) %>%
      ungroup()
    
    # If only one imputed dataset is used then need "IF" statement to calculate SE's
    if (length(dsa_combined_summarised_output$category[dsa_combined_summarised_output$category == "Total"]) == 1) 
    {
      dsa_summarised_output <- dsa_combined_summarised_output %>%
        mutate(num_individuals = n)
    }
    
    # Calculate mean costs and standard errors for the current category
    v.C_ntrt <- dsa_summarised_output$tc_ntrt
    v.C_trt <- dsa_summarised_output$tc_trt
    v.C <- c(v.C_ntrt, v.C_trt)
    
    # Calculate mean QALYs and standard errors for the current category
    v.E_ntrt <- dsa_summarised_output$te_ntrt
    v.E_trt <- dsa_summarised_output$te_trt
    v.E <- c(v.E_ntrt, v.E_trt)
    
    # calculate incremental costs
    delta.C <- v.C[2] - v.C[1]      
    # calculate incremental QALYs
    delta.E <- v.E[2] - v.E[1]    
    
    # calculate the INMB
    inhb_values <- delta.E - (delta.C / wtp)
    
    # Store the INMB for this perturbation
    result_list <- list(inhb_values = inhb_values, 
                        perturbed_values = modified_values)
    
    return(result_list)
  }
  
  
  #### Cost-effectiveness plane function ----
  
  create_ce_plane <- function(data, x_var, y_var, color_var, title, xlim) {
    ce_plane <- ggplot(data, aes(x = {{x_var}}, y = {{y_var}}, color = factor({{color_var}}))) +
      geom_point(size = 0.5) +
      scale_fill_viridis(discrete = TRUE) +
      geom_abline(intercept = 0, color = "black", slope = wtp, linetype = "dashed") +
      geom_hline(yintercept = 0, color = "grey", linetype = "solid") +
      geom_vline(xintercept = 0, color = "grey", linetype = "solid") +
      labs(title = title,
           x = "Incremental Effect",
           y = "Incremental Cost") +
      theme_minimal() +
      theme(legend.position = c(0.9, 0.9)) +
      coord_cartesian(xlim = xlim) 
    
    return(ce_plane)  # Return the plot object
  }  
  
  #### INHB Bar-chart function ----     
  
  create_grouped_bar_charts <- function(data, variable_values, var, var_title) {
    charts <- list()
    var_sym <- sym(var)
    
    for (variable_value in variable_values) {
      filtered_data <- data %>% filter(variable == variable_value)
      
      # Reorder the 'category' variable based on 'inhb' in descending order
      filtered_data <- filtered_data %>%
        mutate(category = fct_reorder(category, -!!var_sym, .fun = min))
      
      chart <- ggplot(filtered_data, aes(x = category, y = !!var_sym, fill = category)) +
        geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
        scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +  # Turn off scientific notation for y-axis
        labs(title = "",
             x = " ",
             y = var_title,
             fill = "Category") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_viridis(discrete = TRUE) 

      charts[[variable_value]] <- chart
    }
    return(charts)
  }
  
  #### Stacked Bar-chart function ----     
  
  create_grouped_stacked_bar_charts <- function(data, variable_values) {
    charts <- list()
    
    for (variable_value in variable_values) {
      filtered_data <- data %>%
        filter(variable == variable_value) %>%
        pivot_longer(cols = c(hc, oop, redist_gov, redist_noop, disab, abs, lp),
                     names_to = "cost_type",
                     values_to = "value")
      
      # Get the minimum tc value for each category
      min_tc_values <- data %>%
        filter(variable == variable_value) %>%
        group_by(category) %>%
        summarize(min_tc = min(tc, na.rm = TRUE), .groups = 'drop')
      
      # Merge min_tc_values back to filtered_data
      filtered_data <- filtered_data %>%
        left_join(min_tc_values, by = "category") %>%
        mutate(category = fct_reorder(category, min_tc))
      
      chart <- ggplot(filtered_data, aes(x = category, y = value, fill = cost_type)) +
        geom_bar(stat = "identity") +
        labs(
          title = paste("Variable:", variable_value),
          x = "Category",
          fill = "Cost Type"
        ) +
        scale_fill_viridis_d(labels = c("hc" = "MI/Stroke healthcare",
                                        "oop" = "Medicaid (OOP)",
                                        "redist_gov" = "Redist. Medicaid (Gov. Admin)",
                                        "redist_noop" = "Redist. Medicaid (non-OOP)",
                                        "disab" = "Short-term disability",
                                        "abs" = "Absenteeism",
                                        "lp" = "Premature Mortality")) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.text = element_text(size = 9),  # Change legend text size here
              legend.key.size = unit(1, "lines"),  # Change legend key size here
              legend.title = element_text(size = 10))  # Change legend title size here
      
      charts[[variable_value]] <- chart
    }
    
    return(charts)
  }
                      
  #### Costs, QALYs, and Net Health Benefit Bar-chart function ----
  
  create_bar_chart <- function(data, var_prefix, title) {
    
    # Reshape the data from wide to long format
    data_long <- pivot_longer(
      data,
      cols = starts_with(var_prefix),
      names_to = "type",
      values_to = "value"
    ) %>%
      mutate(type = ifelse(str_detect(type, "ntrt$"), "No Medicaid Expansion", "Medicaid Expansion"))
    
    # Extract and order categories based on 'No Medicaid Expansion' values
    no_med_exp_values <-data_long %>%
      filter(type == "No Medicaid Expansion") %>%
      arrange(desc(value))
    
    # Update category levels in bar_data_long based on sorted no_med_exp_values
    data_long$category <- factor(data_long$category, levels = no_med_exp_values$category)
    
    # Create the bar chart
    p <- ggplot(data_long, aes(x = category, y = value, fill = type)) +
      geom_bar(stat = "identity", 
               position = position_dodge(width = 0.8), 
               width = 0.6) +
      labs(x = " ", y = "QALYs", fill = "Treatment Type") +
      scale_fill_viridis_d() +
      theme_minimal() +
      theme(legend.title = element_blank(),
            legend.position = "bottom",
            legend.justification = "center",
            legend.box.just = "center",
            legend.key.size = unit(0.75, "lines"),
            axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
      facet_wrap(~variable, scales = "free_x")
    
    # Set the plot title based on var_prefix
    p <- p + ggtitle(title)
    
    return(p)
  }
  
  #### Cost-Effectiveness Acceptability Curve Function ----
  
  # Function to create CEAC plot for given cost and effect variables
  create_ceac_plot <- function(data, cost_var, effect_var, title_suffix) {
    # Calculate pr_CE for each value of wtp using sapply
    calculate_pr_ce <- function(wtp) {
      pr_CE <- ifelse(((data[[effect_var]] * wtp) - data[[cost_var]]) > 0, 1, 0)
      mean(pr_CE)
    }
    
    pr_CE <- sapply(wtp_range, calculate_pr_ce)
    pr_nCE <- 1 - pr_CE
    
    # Create a data frame for plotting
    ceac_plot_data <- data.frame(wtp_range = wtp_range, pr_CE = pr_CE, pr_nCE = pr_nCE)
    
    # Create the plot using ggplot2
    ggplot(ceac_plot_data, aes(x = wtp_range)) +
      geom_line(aes(y = pr_CE, color = "Medicaid Expansion"), linetype = "solid", linewidth = 1) +
      geom_line(aes(y = pr_nCE, color = "No Medicaid Expansion"), linetype = "solid", linewidth = 1) +
      geom_vline(xintercept = 100000, linetype = "dashed", color = "black") +
      annotate("text", x = 100000, y = 0.05, label = "$100,000", vjust = 0.5, hjust = -0.05, color = "black") +
      labs(title = paste(" ", title_suffix),
           x = "Cost-effectiveness threshold",
           y = "Probability of cost-effectiveness",
           color = " ") +
      scale_color_manual(values = c("Medicaid Expansion" = "blue", "No Medicaid Expansion" = "red")) +
      theme_minimal() +
      theme(legend.position = c(0.5, 0.9)) +
      coord_cartesian(ylim = c(0, 1))
  }
  
  #### Equally Distributed Equivalent Health function ----
  
  calculate_ede <- function(health_scores, epsilon, kp = FALSE, population_size) {

    if (kp == FALSE) {
      if (epsilon == 1) {
        ede <- log(health_scores)
      } else {
        ede_hscore <- ifelse(health_scores >= 0, 
                             health_scores ^ (1 - epsilon),
                             -(-health_scores) ^ (1 - epsilon))
        ede_sum <- sum(ede_hscore * (population_size/sum(population_size)), na.rm = TRUE)
        ede <- ede_sum ^ (1 / (1 - epsilon))
      }
    } else {
      ede_sum <- log(sum(exp(-(health_scores * epsilon)) * (population_size/sum(population_size)), na.rm = TRUE))
      ede <- -(ede_sum / epsilon)
    }
    
    return(ede)
  }
  
  
} # Close "FUNCTIONS" section
