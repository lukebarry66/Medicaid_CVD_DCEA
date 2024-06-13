#------------------------------------------------------------------------------#
#---------------------Purpose: MSM DCEA Medicaid-CVD---------------------------#
#----------------------------------2-CLEAN-------------------------------------#
#------------------------------------------------------------------------------#


#### New sample weights combined across waves ----
# https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx

# Calculate the number of years based on unique values in "wave_year"
num_years <- length(unique(nhanes_orig$wave_year)) * 2

# Loop through the variables you want to divide
variables_to_divide <- c("WTINT2YR", "WTMEC2YR", "WTSAF2YR")

for (variable in variables_to_divide) {
  # Create the new variable name with the updated number
  new_variable_name <- paste0(substr(variable, 1, 5))
  
  # Perform the transformation and assign the result
  nhanes_orig[[new_variable_name]] <- nhanes_orig[[variable]] / num_years
}

#---------------------------------Setting up data------------------------
variable_names <- c("smoking", "stroke_dx", "heart_attack_dx", "chd_dx",
                    "cong_heart_fail_dx", "angina_dx", "meds_hyperten", 
                    "diabetes_dx", "diabetes", "history_cvd")

# remove those older and younger
nhanes_gtsummary <- nhanes_orig %>%
  filter(# no missing weights
    !is.na(WTSAF)
    # non-elderly adults
    , age_years <= .stop_age 
    , age_years >= .start_age
    ) %>%
  mutate(insurance = case_when(health_ins_private == "1"  ~ 1,          
                               medicare           == "1"  ~ 2,
                               medicaid           == "1"  ~ 3,
                               health_ins_cover   == "1"  ~ 4,
                               health_ins_cover   == "0"  ~ 0,
                               TRUE                       ~ NA),
         medicaid_exp   = ifelse(insurance  == 0 & ratio_poverty_inc < .fpl & age_years < 65, 1, 0),
         # create single insurance variable to correspond to MEPS cost and utility calculations and for cleaner distinction of treatment scenario
         education = case_when((edu_lev_5cat == "1" | edu_lev_5cat == "2")  ~ 0
                               , edu_lev_5cat == "3"  ~ 1
                               , edu_lev_5cat == "4"  ~ 2
                               , edu_lev_5cat == "5"  ~ 3
                               , TRUE                 ~ NA),
         # create single race variable to correspond to MEPS cost and utility calculations and for cleaner distinction of treatment scenario
         race = case_when(race_5cat == "1"    ~ 0                    
                          , race_5cat == "2"  ~ 1
                          , race_5cat == "3"  ~ 2
                          , race_5cat == "4"  ~ 3
                          , race_5cat == "5"  ~ 4
                          , TRUE              ~ NA),   
         # create "current smoker" (smoking) variable to correspond to MEPS cost and utility calculations and ASCVD calculation
         smoking = case_when(smoke_3cat         == "1"  ~ 1,                  
                             smoke_3cat         == "2"  ~ 0,
                             smoke_3cat         == "0"  ~ 0,
                             TRUE                       ~ NA),
         # create a family income variable to correspond to MEPS cost and utility calculations and according to the poverty line [see, https://meps.ahrq.gov/survey_comp/hc_technical_notes.shtml]
         fam_income = case_when(ratio_poverty_inc   <= 1     ~ 0                          
                                , ratio_poverty_inc <= 1.25  ~ 1
                                , ratio_poverty_inc <= 2     ~ 2
                                , ratio_poverty_inc <= 4     ~ 3
                                , ratio_poverty_inc  > 4     ~ 4
                                , TRUE         ~ NA),
         # same for bmi categories as family income
         bmi_cat = case_when(BMXBMI < 18.5                 ~ 0
                             , (BMXBMI >= 18.5 & BMXBMI < 25) ~ 1
                             , (BMXBMI >= 25 & BMXBMI < 30)   ~ 2
                             , BMXBMI >= 30                ~ 3
                             , TRUE              ~ NA),
         diabetes = ifelse((diabetes_dx == 1 | (hba1c >= 6.5 & diabetes_dx == 0)), 1, 0),
         wave_year = case_when(  wave_year == "G_11" ~ 0
                                 , wave_year == "H_13" ~ 1
                                 , wave_year == "I_15" ~ 2
                                 , wave_year == "J_17" ~ 3
                                 , TRUE ~ NA),
         history_cvd = ifelse((heart_attack_dx == 1 | stroke_dx == 1 | chd_dx == 1 | cong_heart_fail_dx == 1 | angina_dx == 1),
                              1, 0)
  ) %>%
  set_variable_labels(age_years = "Age in years"
                      , ratio_poverty_inc = "Ratio of family income to FPL"
                      , hdl_chol_mgdl = "HDL-Cholesterol (mg/dL)"
                      , total_chol_mgdl = "Total Cholesterol (mg/dL)"
                      , BMXBMI = "Body Mass Index (kg/m^2)"
                      , bp_systolic = "Systolic blood pressure (mm Hg)"
                      , hba1c = "Glycohemoglobin (%); HbA1c"
                      , heart_attack_dx = "Myocardial infarction"
                      , stroke_dx = "Stroke"
                      , cong_heart_fail_dx = "Congestive heart failure"
                      , chd_dx = "Coronary heart disease"
                      , angina_dx = "Angina pectoris"
                      , diabetes_dx = "Diabetes (doctor diagnosed)"
                      , diabetes = "Diabetes (doctor diagnosed OR HbA1c >= 6.5%)"
                      , meds_hyperten = "Currently taking hypertension medication"
                      , male = "Gender"
                      , insurance = "Insurance Status"
                      , education = "Highest education"
                      , race = "Race/Ethnicity"
                      , smoking = "Current smoker"
                      , medicaid_exp = "Medicaid eligibility under expansion"
                      , fam_income = "Family income"
                      , bmi_cat = "Body Mass Index (kg/m^2) category"
                      , history_cvd = "History of CVD*") %>%
  set_value_labels(
    insurance = c("Uninsured" = 0, "Private" = 1, "Medicare" = 2, "Medicaid" = 3, "Other plan" = 4),
    fam_income = c("Poor" = 0, "Near poor" = 1, "Low" = 2, "Medium" = 3, "High" = 4),
    bmi_cat = c("Underweight" = 0, "Normal weight" = 1, "Overweight" = 2, "Obese" = 3),
    education = c("No degree" = 0, "GED/HS" = 1, "Associate/Bachelor" = 2, "Master/Doctorate" = 3),
    race = c("Hispanic" = 0, "White" = 1, "Black" = 2, "Asian" = 3, "Other race" = 4),
    male = c("Male" = 1, "Female" = 0),
    medicaid_exp = c("Eligible" = 1, "Ineligible" = 0),
    wave_year = c("2011/12" = 0, "2013/14" = 1, "2015/16" = 2, "2017/18" = 3)
  ) %>%
  mutate_at(vars(all_of(variable_names)), ~ set_value_labels(., c("No" = 0, "Yes" = 1))) %>%
  select(-WTSAF
         , -WTINT2YR
         , -WTINT
         , -WTMEC
         , -WTMEC2YR
         , -WTSAF2YR
         , -race_5cat
         , -edu_lev_5cat
         , -ratio_poverty_inc
         , -health_ins_cover
         , -health_ins_private
         , -medicare
         , -medicaid
         , -smoke_3cat
         , -BMXBMI
         , -hba1c
         , -diabetes_dx) 


#### Table 1 ----
{
  
  #---------------------------------Creating functions------------------------
  
  # Define a function to check if a variable is categorical
  is_categorical_var <- function(var){
    # Check if the variable is a factor or has fewer than 10 unique non-missing values
    return(is.factor(var) || (length(unique(var[!is.na(var)])) <= 10))
  }
  
  # Define a function to check if a variable is binary
  is_binary_var <- function(var){
    # Check if the variable is a factor or has exactly 2 unique non-missing values
    return(is.factor(var) || (length(unique(var[!is.na(var)])) == 2))
  }
  
  # Define a function to check if a variable is continuous.
  is_continous_var <- function(var){
    # Check if the variable is labelled (assuming the `is.labelled` function is used) and convert it to a factor if so.
    if (is.labelled(var)) {
      var = as.factor(var)
    }
    # Check if the number of distinct non-missing values is not less than or equal to 2 and if the variable is not a factor.
    return(!(n_distinct(var[!is.na(var)]) <= 2 | is.factor(var)))
  }
  
  # Define a function to replace missing values with a specified value and label them.
  to_label_missing <- function(var, value = 999, label_missing = "N missing") {
    val_label(var, value) <- label_missing
    return(var)
  }
  
  # Define a function to replace missing values with a specified value.
  missing_to_value <- function(variable, value = 0) {
    newx <- replace(variable, which(variable %in% NA), value)
  }
  
  # Define a function to add missing value categories based on variable type.
  add_missing_categories <- function(.data) {
    .data %>% 
      # Replace missing values with 999 and label them for categorical variables.
      mutate_if(is_categorical_var, missing_to_value, value = 999) %>% 
      mutate_if(is_categorical_var, to_label_missing) %>% 
      # Replace missing values with 999 and label them for binary variables.
      mutate_if(is_binary_var, missing_to_value, value = 999) %>% 
      mutate_if(is_binary_var, to_label_missing)
  }
  
  #---------------------------------Creating table 1------------------------------
  
  #With missing in gtsummary---
  # to get table across waves working, need to separate wave_year from data to reuse
  nhanes_wave_year <- nhanes_gtsummary %>%
    select(SEQN,
           wave_year)
  
  # Update the table to include a column for each "wave_year"
  wave_summary_miss <- nhanes_gtsummary %>% 
    select(-wave_year) %>%
    add_missing_categories() %>% 
    modify_if(is.labelled, to_factor) %>% 
    left_join(nhanes_wave_year, by = join_by(SEQN)) %>%
    select(-SEQN) %>%
    tbl_summary(
      missing = "no", 
      by = wave_year,  # Split the table by "wave_year" variable
      type = list(
        all_continuous() ~ "continuous2"),
      statistic = list(
        all_continuous() ~ c("{mean} ({sd})",
                             "{N_miss} ({p_miss}%)"),
        all_categorical() ~ "{n} ({p}%)")
    ) %>%
    modify_header(label = "**Variable**") %>%
    add_stat_label(location = "row") %>%
    modify_spanning_header(starts_with("stat_") ~ "Table 1") %>% 
    bold_labels()
  
  print(wave_summary_miss)
  
  
  #---------------------------------Saving the table in word/pptx------------------------------
  
  wave_summary_miss %>% as_flex_table() %>% 
    flextable::save_as_docx(path=here("tables", "table1_miss.docx"))
  
  wave_summary_miss %>% as_flex_table() %>% 
    flextable::save_as_pptx(path=here("tables", "table1_miss.pptx"))
  
}

rm(nhanes_gtsummary
   , nhanes_wave_year)


#----------------------------NHANES ANALYSIS DATA PREP-------------------------#

#### Load the data ----
## Load the .rds file
nhanes_orig <- nhanes_orig %>%
  as_tibble() %>%
  select(age_years
         , stroke_dx
         , heart_attack_dx
         , chd_dx
         , cong_heart_fail_dx
         , angina_dx
         , race_5cat
         , ratio_poverty_inc
         , edu_lev_5cat
         , hdl_chol_mgdl
         , total_chol_mgdl
         , health_ins_cover
         , health_ins_private
         , medicare
         , medicaid
         , smoke_3cat
         , BMXBMI
         , meds_hyperten
         , bp_systolic
         , diabetes_dx
         , male
         , hba1c
         , wave_year
         , WTINT2YR
         , WTMEC2YR
         , WTSAF2YR
         , SDMVPSU
         , SDMVSTRA
         , SEQN) %>% 
  filter(wave_year == "J_17" | wave_year == "I_15" | wave_year == "H_13" | 
           wave_year == "G_11") %>%
  print(n = 10, width = Inf)

#### New sample weights combined across waves ----
# https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx

# Calculate the number of years based on unique values in "wave_year"
num_years <- length(unique(nhanes_orig$wave_year)) * 2

# Loop through the variables you want to divide
variables_to_divide <- c("WTINT2YR", "WTMEC2YR", "WTSAF2YR")

for (variable in variables_to_divide) {
  # Create the new variable name with the updated number
  new_variable_name <- paste0(substr(variable, 1, 5))
  
  # Perform the transformation and assign the result
  nhanes_orig[[new_variable_name]] <- nhanes_orig[[variable]] / num_years
}

#### Format the data ----

nhanes <- nhanes_orig %>% 
  # rename variables to correspond to those used in the functions
  rename(hdl_cholesterol   = hdl_chol_mgdl,                                     
         total_cholesterol = total_chol_mgdl,
         mean_systolic_bp  = bp_systolic,
         fpl               = ratio_poverty_inc,
         bmi               = BMXBMI
  ) %>%
  mutate(age           = as.numeric(age_years),                                
         stroke        = as.numeric(stroke_dx),                                
         mi            = as.numeric(heart_attack_dx),                          
         heart_failure = as.numeric(cong_heart_fail_dx),                       
         chd           = as.numeric(chd_dx), 
         angina        = as.numeric(angina_dx),
         diabetes_ever = as.numeric(diabetes_dx),                              
         treated       = as.numeric(meds_hyperten),                            
         hba1c         = as.numeric(hba1c),                                    
         female = if_else(male == 1, 0, 1), 
         #create single insurance variable to correspond to MEPS cost and utility calculations and for cleaner distinction of treatment scenario
         insurance = case_when(health_ins_private == "1"  ~ "private",          
                               medicare           == "1"  ~ "medicare",
                               medicaid           == "1"  ~ "medicaid",
                               health_ins_cover   == "1"  ~ "other_plan",
                               health_ins_cover   == "0"  ~ "uninsured",
                               TRUE                       ~ NA),
         # create single insurance variable to correspond to MEPS cost and utility calculations and for cleaner distinction of treatment scenario
         education = case_when(edu_lev_5cat == "1"  ~ "no_degree",              
                               edu_lev_5cat == "2"  ~ "no_degree",
                               edu_lev_5cat == "3"  ~ "ged_hs",
                               edu_lev_5cat == "4"  ~ "associate_bachelor",
                               edu_lev_5cat == "5"  ~ "master_doctorate",
                               TRUE                 ~ NA),
         # create single race variable to correspond to MEPS cost and utility calculations and for cleaner distinction of treatment scenario
         race = case_when(race_5cat == "1"  ~ "hispanic",                       
                          race_5cat == "2"  ~ "white",
                          race_5cat == "3"  ~ "black",
                          race_5cat == "4"  ~ "asian",
                          race_5cat == "5"  ~ "other_race",
                          TRUE              ~ NA),   
         # create "current smoker" (smoking) variable to correspond to MEPS cost and utility calculations and ASCVD calculation
         smoking = case_when(smoke_3cat         == "1"  ~ "1",                  
                             smoke_3cat         == "2"  ~ "0",
                             smoke_3cat         == "0"  ~ "0",
                             TRUE                       ~ NA),
         smoking = as.numeric(smoking)
         ) %>%
  # drop unused variables after formatting
  select(-male
         , -age_years
         , -smoke_3cat
         , -stroke_dx
         , -angina_dx
         , -heart_attack_dx
         , -cong_heart_fail_dx
         , -chd_dx
         , -diabetes_dx
         , -meds_hyperten
         , -race_5cat
         , -edu_lev_5cat
         , -health_ins_cover
         , -health_ins_private
         , -medicare
         , -medicaid
         , -WTINT2YR
         , -WTINT
         , -WTMEC2YR
         , -WTMEC)

### Save R file
# output data before applying exclusion criteria
write_rds(nhanes, here("data", "output_data", "nhanes.rds"))

#### Examine Missingness ----

{
  ### FLOW DIAGRAM DATA ----
  {
    ### SECTION 1 ----
    # number of participants per wave
    length(nhanes$SEQN[nhanes$wave_year == "G_11"])
    length(nhanes$SEQN[nhanes$wave_year == "H_13"])
    length(nhanes$SEQN[nhanes$wave_year == "I_15"])
    length(nhanes$SEQN[nhanes$wave_year == "J_17"])
    
    # number of participants across waves = nhanes
    length(nhanes$SEQN)
    
    # number of participants without fasting blood sample weights
    sum(is.na(nhanes$WTSAF))
    
    # proportion of participants without fasting blood sample weights
    sum(is.na(nhanes$WTSAF))/length(nhanes$SEQN)
    
    ### SECTION 2 ----
    # remove those without fasting blood sample weights
    nhanes_filter <- nhanes %>%
      filter(!is.na(WTSAF))
    
    # number of participants remaining = nhanes_filter
    length(nhanes_filter$SEQN)
    
    # number of participants age 18 years or younger
    length(nhanes_filter$SEQN[nhanes_filter$age < .start_age])
    
    # corresponding proportion
    length(nhanes_filter$SEQN[nhanes_filter$age < .start_age]) / length(nhanes_filter$SEQN)
    
    # number of participants age 65 years or older
    length(nhanes_filter$SEQN[nhanes_filter$age > .stop_age])
    
    # corresponding proportion
    length(nhanes_filter$SEQN[nhanes_filter$age > .stop_age]) / length(nhanes_filter$SEQN)
    
    ### SECTION 3 ----
    # remove those older and younger
    nhanes_filter <- nhanes_filter %>%
      filter(age <= .stop_age, # non-elderly adults
             age >= .start_age)
    # number of participants remaining = nhanes_filter
    length(nhanes_filter$SEQN)
    
    ### SECTION 4 ----
    # Present complete-case data
    # Function to calculate missing values and percentages
    missing_summary <- function(x) {
      missing_count <- sum(is.na(x))
      missing_percentage <- (missing_count / length(x)) * 100
      return(c(missing_count, missing_percentage))
    }
    # Calculate and print missing values and percentages
    missing_summary_df <- nhanes_filter %>%
      summarise(across(everything(), ~ missing_summary(.))) %>%
      t() %>%
      as.data.frame()
    colnames(missing_summary_df) <- c("Missing Count", "Missing Percentage")
    
    # Calculate the total number of observations with missing values
    total_missing_observations <- sum(apply(nhanes_filter, 1, anyNA))
    total_missing_percentage <- (total_missing_observations / nrow(nhanes_filter)) * 100
    
    # Add total missing observations and proportion to the data frame
    total_missing <- c(total_missing_observations, total_missing_percentage)
    missing_summary_df <- rbind(missing_summary_df, total_missing)
    rownames(missing_summary_df)[nrow(missing_summary_df)] <- "Total"
    
    # Assuming missing_summary_df is your data frame
    # To sort in descending order, you can use the 'decreasing' argument
    sorted_missing_summary_df <- missing_summary_df[order(missing_summary_df$`Missing Count`, decreasing = TRUE), ]
    print(sorted_missing_summary_df)
    
    ### SECTION 5 ----
    # remove missing values
    nhanes_miss <- nhanes_filter %>%
      na.omit() %>%
      mutate(m_exp = case_when((insurance == "uninsured" & fpl < 1.38) ~ 1,
                               (insurance != "uninsured" & !is.na(insurance)) | (fpl >= 1.38 & !is.na(fpl)) ~ 0,
                               TRUE ~ NA))
    
    # count after missing values removed
    length(nhanes_miss$SEQN)
    
    # count eligble for medicaid
    table(nhanes_miss$m_exp)
    # proportion eligble for medicaid
    prop_mexp <- sum(nhanes_miss$m_exp[nhanes_miss$m_exp==1]) / length(nhanes_miss$m_exp[!is.na(nhanes_miss$m_exp)])
    print(prop_mexp)
    print(1-prop_mexp)
    
    #remove unnecessary data
    rm(nhanes_miss)
  }
  }
