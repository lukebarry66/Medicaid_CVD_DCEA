#------------------------------------------------------------------------------#
#---------------------Purpose: MSM DCEA Medicaid-CVD---------------------------#
#----------------------------------3-IMPUTE------------------------------------#
#------------------------------------------------------------------------------#

#### Multiple Imputation ----

nhanes_MI <- nhanes %>%
  mutate(stroke        = as.character(stroke),
         mi            = as.character(mi),
         heart_failure = as.character(heart_failure),
         chd           = as.character(chd),
         angina        = as.character(angina),
         diabetes_ever = as.character(diabetes_ever),
         treated       = as.character(treated),
         smoking       = as.character(smoking)) 

numeric_vars <- c("hdl_cholesterol", "total_cholesterol", "mean_systolic_bp", "hba1c", "bmi", "fpl", "age")
cat_vars <- c("stroke", "mi", "heart_failure", "chd", "angina", "diabetes_ever",
              "treated", "education", "smoking", "race", "insurance", "female")
all_vars <- c(numeric_vars, cat_vars)

# Perform mice, return m datasets 
seqTime <- system.time(
  miceObj <- miceRanger(
    nhanes_MI
    , m = m 
    , maxiter = m_iter
    , vars = all_vars
    , returnModels = TRUE
    , verbose = FALSE
    , parallel = FALSE
  )
)
miceObj

# Save m datasets
dataList <- completeData(miceObj)
head(dataList[[1]],m)

#remove unnecessary data
rm(nhanes_MI)

# Save the list to a binary file
saveRDS(dataList, here("data", "output_data", "dataList.rds"))
