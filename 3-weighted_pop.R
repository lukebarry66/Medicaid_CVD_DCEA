#------------------------------------------------------------------------------#
#---------------------Purpose: MSM DCEA Medicaid-CVD---------------------------#
#-----------------------------Weighted counts----------------------------------#
#------------------------------------------------------------------------------#

#### Weighted counts ----

.weight <- "WTSAF2YR"
.psu    <- "SDMVPSU"
.strata <- "SDMVSTRA"
options(survey.lonely.psu = "adjust") 

#--------------------Creating functions--------------------------
#--Create function to display point estimate and confidence interval
ci_fct <- function(pe, low, high, multiplier=1, digits=1, scale="%", count=FALSE){
  if(count==F){
    res <- paste(multiplier*formattable::digits(pe,digits),scale, " " ,"(", 
                 multiplier*formattable::digits(low,digits),",", 
                 multiplier*formattable::digits(high,digits), ")", sep="")
  } else {
    res <- paste(formattable::comma(pe, digits=digits), " " ,"(", 
                 formattable::comma(low, digits=digits),", ", 
                 formattable::comma(high, digits=digits), ")", sep="")
  }
  return(res)
}

#--Test the function
ci_fct(pe=2, low=1, high=3, multiplier=100, scale="%")
ci_fct(pe=200000, low=200000, high=200000, multiplier=NULL, digits=0, scale=NULL, count=T)

#--Create function to calculate the count

count_function <- function(.data, variable){
  .data %>% 
    count(!!sym(variable)) %>% 
    rename(value = variable,
           unweighted_count=n) %>% 
    mutate(variables=variable) %>% 
    relocate(variables, value, unweighted_count)
}

#to create a binary variable
d <- function(expr){
  as.numeric(expr)
}

#--------------------Obtaining weighted and unweighted counts-------------------

##data prep
nhanes_wt <- nhanes_full %>%
  filter(!is.na(.data[[.weight]])
         , age <= .stop_age
         , age >= .start_age
         , wave_year == "H_13") %>% # Restrict to wave relelvant to Medicaid expansion
  # create a family income variable to correspond to MEPS cost and utility calculations and according to the poverty line 
  #see, https://meps.ahrq.gov/survey_comp/hc_technical_notes.shtml
  mutate(medicaid_exp = ifelse(insurance  == "uninsured" & fpl < .fpl & age < 65, 1, 0)) %>%
  #na.omit() %>%
  filter(medicaid_exp == 1 | .ATT != 1) %>%
  select(  .data[[.weight]]
         , .data[[.psu]]
         , .data[[.strata]]) %>%
  mutate(  total = 1) 

##Estimating weighted weighted_count----
weighted_count_data_ci <- nhanes_wt %>% 
  as_survey_design(weights = .data[[.weight]]
                   , strata = .data[[.strata]]
                   , ids = .data[[.psu]]
                   , nest = TRUE
                   ) %>% 
  summarize_all(survey_total, vartype = "ci", na.rm=T) %>% 
  t() %>%
  data.frame(weighted_count=.) %>%
  rownames_to_column(var = "variables")

weighted_count_data_low <- weighted_count_data_ci %>% 
  filter(str_detect(variables, '_low')) %>% 
  mutate(variables=str_replace_all(variables,"_low", "")) %>%
  rename(weighted_count_low=weighted_count)

weighted_count_data_high <- weighted_count_data_ci %>% 
  filter(str_detect(variables, '_upp')) %>% 
  mutate(variables=str_replace_all(variables,"_upp", "")) %>%
  rename(weighted_count_high=weighted_count)

weighted_count_data_pe <- weighted_count_data_ci %>% 
  filter(variables %in% c(
    "total")) %>% 
  rename(weighted_count_pe=weighted_count)

weighted_count_data <- weighted_count_data_pe %>% 
  full_join(weighted_count_data_low, by=c("variables")) %>% 
  full_join(weighted_count_data_high, by=c("variables")) 

table1_us_weighted_count <- weighted_count_data %>% 
  transmute(
    `Pop. characteristic` = case_when(variables== "total" ~ "Total"
                                      , TRUE	~	NA_character_),
    `NHANES weighted_count`=ci_fct(pe=weighted_count_pe, low=weighted_count_low, 
                                  high=weighted_count_high, multiplier=NULL, digits=0, scale=NULL, count=T),
    `NHANES weighted_count pe`= paste(formattable::comma(weighted_count_pe, digits = 0), sep=""))

##Outputing table 1 - weighted ----
table1_us_weighted_count %>% 
  as_hux()

write_xlsx(table1_us_weighted_count, 
           here("tables", "table1_us_overall_weighted_count.xlsx"))

#### Population counts ----

# create weighted population counts to merge
pop_counts <- weighted_count_data_pe %>%
  as_tibble() %>%
  rename(category = variables,
         pop_counts = weighted_count_pe) %>%
  mutate(category = case_when(  category == "total" ~ "Total"),
         category = as_factor(category)) 
