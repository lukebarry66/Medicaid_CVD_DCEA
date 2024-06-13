#------------------------------------------------------------------------------#
#---------------------Purpose: MSM DCEA Medicaid-CVD---------------------------#
#--------------------------------6-RESULTS-------------------------------------#
#------------------------------------------------------------------------------#

#### Summarize the results  ----

# Define a list of variables for which you want to calculate means
variables_to_mean <- c("tc_ntrt", "tc_trt"
                       , "te_ntrt", "te_trt"
                       , "hc_trt", "hc_ntrt"
                       , "d_oopc_trt", "d_oopc_ntrt"
                       , "d_noopc_trt", "d_noopc_ntrt"
                       , "govc_trt", "govc_ntrt"
                       , "disabc_trt", "disabc_ntrt", "absentc_trt", "absentc_ntrt"
                       , "lpc_trt", "lpc_ntrt"
                       , "r_gov_trt", "r_gov_ntrt", "r_noop_trt", "r_noop_ntrt"
                       , "mi_ntrt", "stroke_ntrt", "CVD_death_ntrt", "nCVD_death_ntrt" 
                       , "mi_trt", "stroke_trt", "CVD_death_trt", "nCVD_death_trt"
                       , "diff_mi", "diff_stroke", "diff_CVDdeath", "diff_nCVDdeath"
                       , "inmb", "inhb", "icer"
                       , "inc_cost", "inc_effect", "nhb_ntrt", "nhb_trt"
                       )

# Initialize an empty list to store summarised outputs
impute_outputs <- list()

# Across a large number of iterations, it is possible to get "Inf" values so
# Columns to replace 'Inf' values
cols_to_check <- c(
  "mi_ntrt", "stroke_ntrt", "CVD_death_ntrt", "nCVD_death_ntrt"
  , "mi_trt", "stroke_trt", "CVD_death_trt", "nCVD_death_trt"
)

# Loop through specified columns and replace 'Inf' values with NA
for (col in cols_to_check) {
  combined_output[[col]][is.infinite(combined_output[[col]])] <- NA
}

# Create the summarised_output for the current dataset
combined_output <- combined_output %>%
  group_by(dataset_index, replication, variable, category) %>%
  summarize(
    across(all_of(variables_to_mean), ~ weighted.mean(., num_individuals, na.rm = TRUE))
    , num_individuals = sum(num_individuals)
  ) %>%
  ungroup() 

for (imputation_index in 1:m) {

    # Create the summarised_output for the current dataset
    summarised_output <- combined_output %>%
      filter(imputation_index == dataset_index) %>%  # Filter data for the current dataset
      group_by(variable, category) %>%
      summarize(
        across(all_of(variables_to_mean), ~ weighted.se.mean(., num_individuals, na.rm = TRUE), .names = "{.col}_se")
        , across(all_of(variables_to_mean), ~ weighted.mean(., num_individuals, na.rm = TRUE))
        , n = sum(num_individuals)
      ) %>%
      ungroup()
  
  
  # Store the summarised_output in the list
  impute_outputs[[paste("summarised_output_", imputation_index, sep = "")]] <- summarised_output
}

# Combine all the summarised outputs into one dataset
combined_summarised_output <- bind_rows(impute_outputs) 

# Calculate the mean for selected variables
summarised_output <- combined_summarised_output %>%
  group_by(variable, category) %>%
  summarize(
    # Use Rubin's rules to combine SEs across imputation results
    across(all_of(variables_to_mean),
           ~ rubin_se(.x, get(paste0(cur_column(), "_se"))), .names = "{.col}_se"),
    across(all_of(variables_to_mean), ~ weighted.mean(., n, na.rm = TRUE)),
    num_individuals = sum(n)
  ) %>%
  ungroup()

# If only one imputed dataset is used then need "IF" statement to calculate SE's
if (length(combined_summarised_output$category[combined_summarised_output$category == "Total"]) == 1) 
{
  summarised_output <- combined_summarised_output %>%
    mutate(num_individuals = n)
}


#### Cost-effectiveness analysis ----

# Identify sample & pop sizes
sample <- summarised_output$num_individuals[summarised_output$category == "Total"]
total_pop <- pop_counts$pop_counts[pop_counts$category == "Total"]

### CEA table - Long

  # Merge pop. counts
  cea_output <- summarised_output %>%
    # create single insurance variable to correspond to MEPS cost and utility calculations and for cleaner distinction of treatment scenario
    mutate(
      category = fct_recode(
        category
        , "No degree"          = "no_degree"
        , "GED/HS"             = "ged_hs"
        , "Associate/Bachelor" = "associate_bachelor"
        , "Master/Doctorate"   = "master_doctorate"
        , "Poor"      = "poor"
        , "Near Poor" = "near_poor"
        , "Low"       = "low"
        , "Medium"    = "medium"
        , "High"      = "high"
        , "Asian"      = "asian"
        , "Hispanic"   = "hispanic"
        , "Other race" = "other_race"
        , "Black"      = "black"
        , "White"      = "white"),
      variable = fct_recode(
        variable
        , "Education" = "education"
        , "Family Income" = "fam_income"
        , "Race/Ethnicity" = "race")
    ) 
  
  # Initialize an empty list to store individual table_micro data frames
  table_list <- list()
  
  # Unique category values
  unique_categories <- unique(cea_output$category)
  
  # Loop through each category
  for (category in unique_categories) {
    
    # Subset data for the current category
    category_data <- cea_output[cea_output$category == category, ]
    
    # Calculate mean costs and standard errors for the current category
    v.C_ntrt <- category_data$tc_ntrt
    v.C_trt <- category_data$tc_trt
    v.C <- c(v.C_ntrt, v.C_trt)
    
    se.C_ntrt <- category_data$tc_ntrt_se
    se.C_trt <- category_data$tc_trt_se
    se.C <- c(se.C_ntrt, se.C_trt)
    
    # Calculate mean QALYs and standard errors for the current category
    v.E_ntrt <- category_data$te_ntrt
    v.E_trt <- category_data$te_trt
    v.E <- c(v.E_ntrt, v.E_trt)
    
    se.E_ntrt <-  category_data$te_ntrt_se
    se.E_trt <- category_data$te_trt_se
    se.E <- c(se.E_ntrt, se.E_trt)
    
    # Calculate incremental costs
    delta.C <- v.C_trt - v.C_ntrt
    delta.E <- v.E_trt - v.E_ntrt
    
    # Calculate net health benefits
    nhb_ntrt <- v.E_ntrt - (v.C_ntrt / wtp)
    nhb_trt  <- v.E_trt - (v.C_trt / wtp)
    
    # calculate the ICER
    ICER    <- delta.C / delta.E  
    # calculate the INMB
    INMB    <- (delta.E * wtp) - delta.C
    inmb_alt   <- category_data$inmb
    inmb_alt_se   <- category_data$inmb_se
    
    # calculate the INHB
    INHB    <- delta.E - (delta.C / wtp)
    inhb_alt   <- category_data$inhb
    inhb_alt_se   <- category_data$inhb_se
    
    n <- category_data$num_individuals

    # Create full incremental cost-effectiveness analysis table
    cea_table_micro <- data.frame(
      c(category, "", ""),
      c("nME", "ME", ""),
      c(round(v.C, 0),  ""),                # costs per arm
      c(round(se.C, 0), ""),                # MCSE for costs
      c(round(v.E, 5),  ""),                # health outcomes per arm
      c(round(se.E, 5), ""),                # MCSE for health outcomes
      c("", round(delta.C, 0),   ""),       # incremental costs
      c("", round(delta.E, 5),   ""),       # incremental QALYs 
      c("", "", ""),
      c(round(nhb_ntrt, 5), round(nhb_trt, 5),   ""),      # net health benefits
      c("", "", ""),
      c("", round(ICER, 0),       ""),       # ICER
      c("", round(INMB, 0),       ""),       # INMB
      c("", round(inmb_alt, 0),       ""),       # INMB
      c("", round(inmb_alt_se, 0),       ""),   # INMB
      c("", round(INHB, 5),       ""),       # INHB
      c("", round(inhb_alt, 5),       ""),       # INMB
      c("", round(inhb_alt_se, 5),       ""),   # INMB
      c("", "", ""),
      c("", round(n, 0), "")  # Number of individuals
    )
    
    # Inside the loop, add the table_micro to the list
    table_list[[category]] <- cea_table_micro
    
  }
  
  # After the loop, combine all the individual table_micro data frames into one
  cea_table_long <- do.call(rbind, table_list)
  
  # name the columns
  colnames(cea_table_long) <- c("Category", "Arm",
                                "Costs", "se",  
                                "QALYs", "se", 
                                "Incremental Costs", 
                                "QALYs Gained", 
                                " ",
                                "Net Health Benefit",
                                " ",
                                "ICER", 
                                "INMB", "INMB_alt", "INMB_alt_se",
                                "INHB", "INHB_alt", "INHB_alt_se",
                                " ",
                                "n") 
  
  

  ### CEA table - Short
  
  #Calculate Total INHB and EDEH
  summary_edeh <- cea_output %>%
    group_by(variable) %>%
    summarise(  ede_ntrt    = (calculate_ede(nhb_ntrt, epsilon, kp = FALSE, population_size = num_individuals)) 
                , ede_trt     = (calculate_ede(nhb_trt, epsilon, kp = FALSE, population_size = num_individuals))
                , ede_ntrt_e05    = (calculate_ede(nhb_ntrt, epsilon = 0.5, kp = FALSE, population_size = num_individuals)) 
                , ede_trt_e05     = (calculate_ede(nhb_trt, epsilon = 0.5, kp = FALSE, population_size = num_individuals))
                , ede_ntrt_e15    = (calculate_ede(nhb_ntrt, epsilon = 1.5, kp = FALSE, population_size = num_individuals)) 
                , ede_trt_e15     = (calculate_ede(nhb_trt, epsilon = 1.5, kp = FALSE, population_size = num_individuals))
                , ede_ntrt_e25    = (calculate_ede(nhb_ntrt, epsilon = 2.5, kp = FALSE, population_size = num_individuals)) 
                , ede_trt_e25     = (calculate_ede(nhb_trt, epsilon = 2.5, kp = FALSE, population_size = num_individuals))
                , nhb_ntrt    = weighted.mean(nhb_ntrt, num_individuals, na.rm = TRUE) 
                , nhb_trt     = weighted.mean(nhb_trt, num_individuals, na.rm = TRUE)  
    ) %>%
    ungroup() %>%
    mutate(edeh = (ede_trt - ede_ntrt)
           , edeh_e05 = (ede_trt_e05 - ede_ntrt_e05)
           , edeh_e15 = (ede_trt_e15 - ede_ntrt_e15)
           , edeh_e25 = (ede_trt_e25 - ede_ntrt_e25)
           , Ae_ntrt     = 1 - (ede_ntrt / nhb_ntrt)
           , Ae_trt      = 1 - (ede_trt / nhb_trt)
           , Ae          = -(Ae_trt - Ae_ntrt)
           , inhb        = (nhb_trt - nhb_ntrt)
           , pop_inhb    = inhb * total_pop
           , ineq_burden_trt = nhb_trt * Ae_trt
           , ineq_burden_ntrt = nhb_ntrt * Ae_ntrt
           , ineq_burden = ((ineq_burden_ntrt - ineq_burden_trt)/ineq_burden_ntrt)
           , equity = edeh/inhb
           , relative_equity = Ae_trt / Ae_ntrt
           , pop_ineq_burden_trt = nhb_trt * Ae_trt * total_pop
           , pop_ineq_burden_ntrt = nhb_ntrt * Ae_ntrt * total_pop
           , pop_ineq_burden = ((pop_ineq_burden_ntrt - pop_ineq_burden_trt))
    )
  
  # Initialize an empty list to store individual table_micro data frames
  table_list <- list()
  
  # Extract unique combinations of variable and category
  unique_combinations <- unique(cea_output[c("variable", "category")])
  
  # Loop through each unique combination
  for (i in seq_along(unique_combinations$variable)) {
    variable_name <- unique_combinations$variable[i]
    category <- unique_combinations$category[i]
    
    # Subset data for the current category within the variable
    category_data <- cea_output[cea_output$variable == variable_name & cea_output$category == category, ]
    
    # Order categories by INHB
    ordered_categories <- category_data[order(category_data$inhb, decreasing = TRUE), ]$category
    
    # Loop through ordered categories
    for (category in ordered_categories) {
      category_data <- category_data[category_data$category == category, ]
      
      # Extract values and SEs directly from the data
      delta.C <- category_data$inc_cost
      delta.C_se <- category_data$inc_cost_se
      delta.E <- category_data$inc_effect
      delta.E_se <- category_data$inc_effect_se
      ICER <- delta.C / delta.E
      INHB <- category_data$inhb
      INHB_se <- category_data$inhb_se
      n <- category_data$num_individuals
      pop_frac <- n/sample
      
      # Calculate 95% confidence intervals
      delta.C_CI <- paste0(round(delta.C, 0), " (", round(delta.C - 1.96 * delta.C_se, 0), " to ", round(delta.C + 1.96 * delta.C_se, 0), ")")
      delta.E_CI <- paste0(round(delta.E, 4), " (", round(delta.E - 1.96 * delta.E_se, 3), " to ", round(delta.E + 1.96 * delta.E_se, 3), ")")
      INHB_CI <- paste0(round(INHB, 4), " (", round(INHB - 1.96 * INHB_se, 3), " to ", round(INHB + 1.96 * INHB_se, 3), ")")

      # Total INHB an EDEH
      edeh_e05 <- summary_edeh$edeh_e05[summary_edeh$variable == variable_name]
      edeh_e15 <- summary_edeh$edeh_e15[summary_edeh$variable == variable_name]
      edeh_e25 <- summary_edeh$edeh_e25[summary_edeh$variable == variable_name]
      
      Ae <- summary_edeh$Ae[summary_edeh$variable == variable_name]
      ineq_burden <- summary_edeh$ineq_burden[summary_edeh$variable == variable_name]
      pop_ineq_burden <- summary_edeh$pop_ineq_burden[summary_edeh$variable == variable_name]
      
      # Create incremental cost-effectiveness analysis table row for the category
      cea_table_row <- data.frame(
        Variable = variable_name
        , Category = category
        , Incremental_Costs = delta.C_CI #  round(delta.C, 0)
        , Incremental_QALYs = delta.E_CI #  round(delta.E, 3)
        , ICER = round(ICER, 0)
        , INHB = INHB_CI
        , Pop_frac = round(pop_frac,2)
        , EDEH_e05 = edeh_e05
        , EDEH_e15 = edeh_e15
        , EDEH_e25 = edeh_e25
        , Ae = Ae
        , ineq_burden = round(ineq_burden*100, 1)
        , pop_ineq_burden = round(pop_ineq_burden, 0)
      )
      
      # Add the table row to the list
      table_list[[paste(variable_name, category)]] <- cea_table_row
    }
  }
  
  # Combine all the individual table rows into one table
  cea_table_short <- do.call(rbind, table_list)

  # Revised column naming to match the data frame construction
  colnames(cea_table_short) <- c("Category"
                           , "Group"
                           , "Incremental Cost (USD 2021)"
                           , "Incremental Effect (QALYs)"
                           , "ICER"
                           , "INHB (QALYs)"
                           , "Population Fraction"
                           , "EDEH (QALYs); epsilon = 0.5"
                           , "EDEH (QALYs); epsilon = 1.5"
                           , "EDEH (QALYs); epsilon = 2.5"
                           , "Atkinson Index of Inequality"
                           , "Inequality burden (%)"
                           , "Inequality burden (pop. QALYs")
  
  # Print the note regarding the NHANES population size
  cat("\nThe weighted NHANES population size in 2013/14 was estimated to be", total_pop, "\n")

#### Incidence rates  ----

# Initialize an empty list to store individual table_micro data frames
table_list <- list()

# Unique category values
unique_categories <- unique(summarised_output$category)

# Loop through each category
for (category in unique_categories) {
  
  # Subset data for the current category
  category_data <- summarised_output[summarised_output$category == category, ]
  
  # Incidence of model outcomes
  mi_ir     <- c(category_data$mi_ntrt, 
                 category_data$mi_trt, 
                 category_data$diff_mi)
  mi_ir_se  <- c(category_data$mi_ntrt_se, 
                 category_data$mi_trt_se, 
                 category_data$diff_mi_se)
  
  stroke_ir     <- c(category_data$stroke_ntrt, 
                     category_data$stroke_trt, 
                     category_data$diff_stroke)
  stroke_ir_se  <- c(category_data$stroke_ntrt_se, 
                     category_data$stroke_trt_se, 
                     category_data$diff_stroke_se)
  
  cvd_death_ir     <- c(category_data$CVD_death_ntrt, 
                        category_data$CVD_death_trt, 
                        category_data$diff_CVDdeath)
  cvd_death_ir_se  <- c(category_data$CVD_death_ntrt_se, 
                        category_data$CVD_death_trt_se, 
                        category_data$diff_CVDdeath_se)
  
  ncvd_death_ir     <- c(category_data$nCVD_death_ntrt, 
                         category_data$nCVD_death_trt,
                         category_data$diff_nCVDdeath)
  ncvd_death_ir_se  <- c(category_data$nCVD_death_ntrt_se, 
                         category_data$nCVD_death_trt_se,
                         category_data$diff_nCVDdeath_se)
  
  n <- category_data$num_individuals
  
  # Create full incremental cost-effectiveness analysis table
  ir_table_micro <- data.frame(
    c(category, "", "Difference", ""),
    c("nME", "ME", "", ""),
    c(round(mi_ir      * .ir_py, 5), ""),       
    c(round(mi_ir_se   * .ir_py, 5), ""),  
    c("", "", "", ""),
    c(round(stroke_ir      * .ir_py, 5), ""),   
    c(round(stroke_ir_se   * .ir_py, 5), ""),
    c("", "", "", ""),
    c(round(cvd_death_ir     * .ir_py, 5), ""),   
    c(round(cvd_death_ir_se  * .ir_py, 5), ""),
    c("", "", "", ""),
    c(round(ncvd_death_ir     * .ir_py, 5), ""),   
    c(round(ncvd_death_ir_se  * .ir_py, 5), ""),
    c("", "", "", ""),
    c(round(n, 0), "", "", "")  # Number of individuals
  )
  
  # Inside the loop, add the table_micro to the list
  table_list[[category]] <- ir_table_micro
  
  # After the loop, combine all the individual table_micro data frames into one
  ir_table <- do.call(rbind, table_list)
  
  # name the columns
  colnames(ir_table) <- c("Category", "Arm",
                          "IR MI",
                          "se",
                          " ",
                          "IR Stroke",
                          "se",
                          " ",
                          "IR CVD-death",
                          "se",
                          " ",
                          "IR Non-CVD_death",
                          "se",
                          " ",
                          "N"
  ) 
  
}

#### Cost-Effectiveness Plane ----

# Create dataset for CE planes
ce_data <- combined_output %>% 
  select(dataset_index 
         , replication
         , variable
         , category
         , tc_ntrt
         , tc_trt
         , te_ntrt
         , te_trt
         , inc_effect
         , inc_cost) 

# List of unique variable values
unique_variables <- unique(ce_data$variable)

# Initialize an empty list to store CE planes
ce_planes_list <- list()

# Loop through each unique value of the "variable" variable
for (var in unique_variables) {
  # Filter the data for the current variable
  var_data <- ce_data %>% filter(variable == var)
  
  # Create the CE plane using the function
  ce_plane <- create_ce_plane(var_data, 
                              inc_effect, 
                              inc_cost,
                              category, 
                              paste("Cost-Effectiveness Plane for", var), 
                              c(min(var_data$inc_effect), max(var_data$inc_effect)))
  
  # Add the CE plane to the list
  ce_planes_list[[var]] <- ce_plane
}

#### Cost-Effectiveness Acceptability Curve ----

# Selecting data for the CE planes
ceac <- combined_output %>% 
  select(dataset_index, replication, variable, category,
         inc_cost, inc_effect) %>% 
  filter(category == "Total")

# Define the range of wtp values and the increment
wtp_range <- seq(0, 1000000, by = 100)

# Create and plot each CEAC curve
ceac_plots <- list(
  base = create_ceac_plot(ceac, "inc_cost", "inc_effect", "Base Case")
)

pr_CE_base <- (ceac_plots$base$data$pr_CE[ceac_plots$base$data$wtp_range == 100000])*100

cat("Medicaid expansion over the lifetime has a", pr_CE_base, "% probability of being cost-effectiveat $100,000/WTP", "\n")

#### Threshold Analysis, % Medicaid costs attributable to CVD ----

calculate_inhb <- function(data, new_cvd_prop, cvd_prop, wtp = 100000) {
  adj_prop <- new_cvd_prop / cvd_prop
  result <- data %>%
    filter(variable == "Total") %>%
    mutate(
      tc1 = tc_trt - tc_ntrt,
      hc = hc_trt - hc_ntrt,
      oop = (d_oopc_trt - d_oopc_ntrt) * adj_prop,
      noop = d_noopc_trt - d_noopc_ntrt,
      gov = govc_trt - govc_ntrt,
      redist_gov = (r_gov_trt - r_gov_ntrt) * adj_prop,
      redist_noop = (r_noop_trt - r_noop_ntrt) * adj_prop,
      disab = disabc_trt - disabc_ntrt,
      abs = absentc_trt - absentc_ntrt,
      lp = lpc_trt - lpc_ntrt,
      tc = hc + oop + redist_gov + redist_noop + disab + abs + lp,
      te = te_trt - te_ntrt,
      inhb = (te) - (tc / wtp)
    ) %>%
    pull(inhb)
  return(mean(result, na.rm = TRUE))  # Return mean in case of multiple entries
}

new_cvd_prop_values <- seq(0.03, 0.55, by = 0.01)
closest_zero <- NULL
min_distance <- Inf

for (prop in new_cvd_prop_values) {
  inhb_value <- calculate_inhb(summarised_output, prop, cvd_prop)
  if (abs(inhb_value) < min_distance) {
    min_distance <- abs(inhb_value)
    closest_zero <- prop
  }
}

cat("The proportion of Medicaid costs attributable to the reduction in CVD risk factors where INHB is closest to zero is", closest_zero*100, "%", "\n")


#### Distributional Analysis ----
{
  ### INHB across equity variables ----
  
  # Create dataset for barcharts
  dist_data <- summarised_output %>% 
    filter(variable != "Total"
           , variable != "age_cat"
           , variable != "female") %>%
    # create single insurance variable to correspond to MEPS cost and utility calculations and for cleaner distinction of treatment scenario
    mutate(
      category = fct_recode(
        category
        , "No degree"          = "no_degree"
        , "GED/HS"             = "ged_hs"
        , "Associate/Bachelor" = "associate_bachelor"
        , "Master/Doctorate"   = "master_doctorate"
        , "Poor"      = "poor"
        , "Near Poor" = "near_poor"
        , "Low"       = "low"
        , "Medium"    = "medium"
        , "High"      = "high"
        , "Asian"      = "asian"
        , "Hispanic"   = "hispanic"
        , "Other race" = "other_race"
        , "Black"      = "black"
        , "White"      = "white"),
      variable = fct_recode(
        variable
        , "Education" = "education"
        , "Family Income" = "fam_income"
        , "Race/Ethnicity" = "race")
    )
  
  bar_data <- dist_data %>%
    select(tc_ntrt, tc_trt
           , te_ntrt, te_trt
           , nhb_trt, nhb_ntrt
           , inhb, variable, category
           , inc_effect, inc_cost
           ) %>%
    filter(variable != "medicaid_exp")
  
  # Usage for tc_ variable prefix
  cost_charts <-create_bar_chart(bar_data, "tc_", title = "Costs per person (USD 2021)")
  
  # Usage for te_ variable prefix
  qaly_charts <- create_bar_chart(bar_data, "te_", title = "QALYs per person")
  
  # Usage for nhb_ variable prefix
  nhb_charts <- create_bar_chart(bar_data, "nhb_", title = "Net Health Benefits per person (QALY)")
  
  # Define unique variable values
  unique_variable_values <- .equity_vars
  #unique(bar_data$variable)

  # Create a list of bar charts with different titles
  INHB_bar_charts_list <- create_grouped_bar_charts(bar_data, unique_variable_values, var = "inhb", var_title = "INHB (QALYs)")
  
  # Modify titles for each individual plot
  titles <- c("Race/Ethnicity", "Education", "Family Income")  # Replace with your desired titles for each plot
  
  for (i in seq_along(INHB_bar_charts_list)) {
    # Modify the plot title
    INHB_bar_charts_list[[i]] <- INHB_bar_charts_list[[i]] +
      labs(title = titles[i]) +  # Update the title for each plot
      theme(plot.title = element_text(size = 10)) +  # Change the title size
      if (i > 1) {  # Condition to remove y-axis title from plots 2 and 3
        theme(axis.title.y = element_blank())  # Remove y-axis title for plots 2 and 3
      }
  }
  
  # Combine plots horizontally using patchwork
  combined_INHB_plots <- wrap_plots(INHB_bar_charts_list, ncol = length(INHB_bar_charts_list))
  
  # Display the combined plots
  print(combined_INHB_plots)
  
  ### Total Effect across equity variables ----
  
  # Create a list of bar charts with different titles
  te_bar_charts_list <- create_grouped_bar_charts(bar_data, unique_variable_values, var = "inc_effect", var_title = "QALY")
  
  # Modify titles for each individual plot
  titles <- c("Race/Ethnicity", "Education", "Family Income")  # Replace with your desired titles for each plot
  
  for (i in seq_along(te_bar_charts_list)) {
    # Modify the plot title
    te_bar_charts_list[[i]] <- te_bar_charts_list[[i]] +
      labs(title = titles[i]) +  # Update the title for each plot
      theme(plot.title = element_text(size = 10)) +  # Change the title size
      if (i > 1) {  # Condition to remove y-axis title from plots 2 and 3
        theme(axis.title.y = element_blank())  # Remove y-axis title for plots 2 and 3
      }
  }
  
  # Combine plots horizontally using patchwork
  combined_effect_plots <- wrap_plots(te_bar_charts_list, ncol = length(te_bar_charts_list))
  
  # Display the combined plots
  print(combined_effect_plots)
  
  ### Total Costs across equity variables ----
  
  # Create a list of bar charts with different titles
  tc_bar_charts_list <- create_grouped_bar_charts(bar_data, unique_variable_values, var = "inc_cost", var_title = "USD (2021)")
  
  # Modify titles for each individual plot
  titles <- c("Race/Ethnicity", "Education", "Family Income")  # Replace with your desired titles for each plot
  
  for (i in seq_along(tc_bar_charts_list)) {
    # Modify the plot title
    tc_bar_charts_list[[i]] <- tc_bar_charts_list[[i]] +
      labs(title = titles[i]) +  # Update the title for each plot
      theme(plot.title = element_text(size = 10)) +  # Change the title size
      if (i > 1) {  # Condition to remove y-axis title from plots 2 and 3
        theme(axis.title.y = element_blank())  # Remove y-axis title for plots 2 and 3
      }
  }
  
  # Combine plots horizontally using patchwork
  combined_cost_plots <- wrap_plots(tc_bar_charts_list, ncol = length(tc_bar_charts_list))
  
  # Display the combined plots
  print(combined_cost_plots)
  
  
  ### Cost breakdown across equity variables ----
  
  # Prepare data
  cost_data <- dist_data %>% 
    mutate(tc = tc_trt - tc_ntrt
           , hc = hc_trt - hc_ntrt
           , oop = d_oopc_trt - d_oopc_ntrt
           , noop = d_noopc_trt - d_noopc_ntrt
           , gov = govc_trt - govc_ntrt
           , redist_gov = r_gov_trt - r_gov_ntrt
           , redist_noop = r_noop_trt - r_noop_ntrt
           , disab = disabc_trt - disabc_ntrt
           , abs = absentc_trt - absentc_ntrt
           , lp = lpc_trt - lpc_ntrt
           , check = hc + oop + redist_gov + redist_noop + disab + abs + lp # should equal tc
    ) %>%
    select(variable, category
           , tc, check, hc
           , oop, disab, abs, lp
           , redist_gov, redist_noop)
  
  # Create grouped stacked bar charts
  stacked_list <- create_grouped_stacked_bar_charts(cost_data, unique_variable_values)
  
  # Modify titles for each individual plot
  titles <- c("Race/Ethnicity", "Education", "Family Income")  # Replace with your desired titles for each plot
  
  # Extract the legend from the first plot
  legend_plot <- stacked_list[[1]]
  legend <- get_legend(legend_plot)
  
  # Update plots to remove individual legends, set titles, and modify y-axis labels
  for (i in seq_along(stacked_list)) {
    stacked_list[[i]] <- stacked_list[[i]] +
      labs(title = titles[i]) +
      theme(plot.title = element_text(size = 10),
            legend.position = "none",
            axis.title.x = element_blank())
    
    # Add "USD 2021" label to the y-axis of the first plot only
    if (i == 1) {
      stacked_list[[i]] <- stacked_list[[i]] +
        labs(y = "USD (2021)")
    } else {
      stacked_list[[i]] <- stacked_list[[i]] +
        theme(axis.title.y = element_blank())
    }
  }
  
  # Combine plots horizontally using cowplot
  combined_stack_plots <- plot_grid(plotlist = stacked_list, ncol = 3, align = "h")
  
  # Add the legend to the right of the combined plots
  stacked_cost_plot <- plot_grid(combined_stack_plots, legend, ncol = 2, rel_widths = c(3, 1)) + 
    guides(fill = guide_legend(override.aes = list(size = 10)))  # Adjust the size as needed
  
  # Display the final plot
  print(stacked_cost_plot)
  
  
  #Proportion absenteeism and disablity costs
  subset_cost_data <- summarised_output[summarised_output$category == "Total", ] %>%
    mutate(tc = abs(tc_trt - tc_ntrt)
           , disab = abs(disabc_trt - disabc_ntrt)
           , abs = abs(absentc_trt - absentc_ntrt)) %>%
    select(tc, disab, abs)

  # relative tot total costs
  prop_abs <- (subset_cost_data$abs / subset_cost_data$tc) * 100
  prop_disab <- (subset_cost_data$disab / subset_cost_data$tc) * 100

  cat("The difference in short-term disability and absenteeism costs accounted for"
      , round(prop_disab, 2), "and", round(prop_abs, 2)
      , "% of the difference in total costs respectively", "\n")
  
  
  ### EQUALLY DISTRIBUTED EQUIVALENT HEALTH (EDEH) ----

  ### Using population weights ###
  
  # Calculate the EDEH for treatment and control across bootstrapped replications
  ede_data_pop <- summarised_output %>%
    filter( variable != "female"
            , variable != "age_cat"
            , variable != "Total"
            , variable != "medicaid_exp"
            
    ) %>%
    mutate(  nhb_ntrt = (te_ntrt) - ((tc_ntrt) / wtp)
             , nhb_trt  = (te_trt) - ((tc_trt) / wtp)
             , inhb     = nhb_trt - nhb_ntrt) %>%
    group_by(variable) %>%
    summarise(  ede_ntrt    = (calculate_ede(nhb_ntrt, epsilon, kp = FALSE, population_size = num_individuals)) 
                , ede_trt     = (calculate_ede(nhb_trt, epsilon, kp = FALSE, population_size = num_individuals))
                , nhb_ntrt    = weighted.mean(nhb_ntrt, num_individuals, na.rm = TRUE) 
                , nhb_trt     = weighted.mean(nhb_trt, num_individuals, na.rm = TRUE)  
    ) %>%
    ungroup() %>%
    mutate(edeh = (ede_trt - ede_ntrt)
           , Ae_ntrt     = 1 - (ede_ntrt / nhb_ntrt)
           , Ae_trt      = 1 - (ede_trt / nhb_trt)
           , Ae          = -(Ae_trt - Ae_ntrt)
           , inhb        = (nhb_trt - nhb_ntrt)
           , pop_inhb    = inhb * total_pop
           , ineq_burden_trt = nhb_trt * Ae_trt
           , ineq_burden_ntrt = nhb_ntrt * Ae_ntrt
           , ineq_burden = ((ineq_burden_ntrt - ineq_burden_trt)/ineq_burden_ntrt)
           , equity = edeh/inhb
           , relative_equity = Ae_trt / Ae_ntrt
           , pop_ineq_burden_trt = nhb_trt * Ae_trt * total_pop
           , pop_ineq_burden_ntrt = nhb_ntrt * Ae_ntrt * total_pop
           , pop_ineq_burden = ((pop_ineq_burden_ntrt - pop_ineq_burden_trt))
    ) %>%
    select(variable,
           edeh,
           Ae,
           inhb, pop_ineq_burden, pop_inhb) %>%
    mutate(sample = "population")
  

  # Creating a new column with custom labels
  ede_data_pop$custom_labels <- ifelse(ede_data_pop$variable == "education", "Education", 
                                       ifelse(ede_data_pop$variable == "fam_income", "Family Income", 
                                              ifelse(ede_data_pop$variable == "race", "Race/Ethnicity", NA)))
  
  # Plotting using the custom labels for text
  ede_plot_pop <- ggplot(ede_data_pop, aes(x = pop_ineq_burden, y = pop_inhb, color = variable, shape = variable, label = custom_labels)) +
    geom_point(size = 6) +
    labs(
      x = "Net Equity Impact (QALYs)",
      y = "Incremental Net Health Benefit (QALYs)",
      color = "Equity-relevant variable",
      shape = "Equity-relevant variable"
    ) +
    scale_color_viridis(
      discrete = TRUE,
      breaks = c("education", "fam_income", "race"),
      labels = c("Education", "Family Income", "Race/Ethnicity")
    ) +
    scale_shape_manual(
      values = c("education" = 16, "fam_income" = 17, "race" = 18),
      guide = guide_legend(override.aes = list(size = 4))
    ) +
    geom_vline(xintercept = 0, color = "black") +
    geom_hline(yintercept = 0, color = "black") +
    theme_minimal() +
    theme(
      legend.position = "none"  # Remove the legend
    ) +
    geom_text(
      vjust = -2, hjust = 0.5,  # Adjust vertical and horizontal justification
      size = 3, color = "black"  # Set size and color of the labels
    ) + 
   coord_cartesian(xlim = c((min(ede_data_pop$pop_ineq_burden) + min(ede_data_pop$pop_ineq_burden*0.3)), 
                            (max(ede_data_pop$pop_ineq_burden) + max(ede_data_pop$pop_ineq_burden*0.3))))  # Set x-axis limits
  
  # Adding non-scientific number formatting
  ede_plot_pop <- ede_plot_pop +
    scale_x_continuous(labels = label_number()) +
    scale_y_continuous(labels = label_number())
  
  print(ede_plot_pop)
  
  ### DCEA BY WTP AND INEQUALITY AVERSION ----
  
  # Define a range of values for wtp and epsilon; see sensitivity analysis parameters
  
  # Define the list of variable values
  variable_values <- c("fam_income", "education", "race")
  
  # Create an empty list to store the plots for each variable
  plots_list <- list()
  
  # Loop over the variable values
  for (v in variable_values) {
    # Create an empty matrix to store the edeh values
    results_matrix <- matrix(NA, nrow = length(wtp_values), ncol = length(epsilon_values))
    
    if (v == "fam_income") {
      name <- "Family Income"
    } else if (v == "education"){
      name <- "Highest Education"
    } else if (v == "race") {
      name <- "Race/Ethnicity"
    }
  
    # Loop over wtp values
    for (i in seq_along(wtp_values)) {
      wtp <- wtp_values[i]
      
      # Loop over epsilon values
      for (j in seq_along(epsilon_values)) {
        epsilon <- epsilon_values[j]
        
        # cat("Evaluating combination: wtp =", wtp, "epsilon =", epsilon, "\n")
        
        # Calculate the EDEH for treatment and control
        edeh_data <- summarised_output %>%
          filter(variable == v) %>%
          select(-nhb_ntrt, -nhb_trt, -inhb) %>%
          mutate(nhb_ntrt = (te_ntrt) - ((tc_ntrt) / wtp)
                 , nhb_trt  = (te_trt) - ((tc_trt) / wtp)
                 ) %>%
          group_by(variable) %>%
          summarise(  ede_ntrt    = (calculate_ede(nhb_ntrt, epsilon, kp = FALSE, population_size = num_individuals)) 
                      , ede_trt     = (calculate_ede(nhb_trt, epsilon, kp = FALSE, population_size = num_individuals))
                      , nhb_ntrt    = weighted.mean(nhb_ntrt, num_individuals, na.rm = TRUE) 
                      , nhb_trt     = weighted.mean(nhb_trt, num_individuals, na.rm = TRUE)  
          ) %>%
          mutate(edeh = (ede_trt - ede_ntrt) * total_pop
                 , Ae_ntrt     = 1 - (ede_ntrt / nhb_ntrt)
                 , Ae_trt      = 1 - (ede_trt / nhb_trt)
                 , Ae          = -(Ae_trt - Ae_ntrt)
                 , inhb        = (nhb_trt - nhb_ntrt) * total_pop
          ) %>%
          select(variable,
                 edeh,
                 Ae)
        
        # Check if any non-numeric values exist in the results
        if (any(!is.numeric(edeh_data$edeh) | is.na(edeh_data$edeh) | !is.finite(edeh_data$edeh))) {
          cat("Problematic combination: wtp =", wtp, "epsilon =", epsilon, "\n")
          # Add additional debugging or print statements here to inspect data for that combination
        }
        
        # Store the edeh value in the results matrix
        edeh_data$CE_threshold <- wtp
        edeh_data$Epsilon <- epsilon
        results_matrix[i, j] <- mean(edeh_data$edeh)
      }
    }
    
    # Create a data frame with correct reshaping
    dcea_data <- data.frame(
      Epsilon = rep(epsilon_values, each = length(wtp_values)),
      CE_threshold = rep(wtp_values, times = length(epsilon_values)),
      EDEH = as.vector(results_matrix)
    )
    
    # Subset the data based on the CE_threshold condition
    subset_data <- dcea_data[dcea_data$CE_threshold == 100000, ]
    
    # Calculate the absolute differences between subset_data$EDEH and 0
    abs_diff <- abs(subset_data$EDEH)
    
    # Find the index of the minimum absolute difference
    min_index <- which.min(abs_diff)
    
    # Extract the corresponding Epsilon value
    closest_epsilon <- subset_data$Epsilon[min_index]
    cat("For", name
        , "Medicaid expansion would be cost-effective (at a value of $100,000/WTP) if epsilon was"
        , closest_epsilon, "or higher", "\n")

    ede_epsilon_plot <- ggplot(dcea_data, aes(x = Epsilon, y = EDEH, color = factor(CE_threshold))) +
      geom_line() +  # Use a line plot
      geom_hline(yintercept = 0, color = "black") +
      labs(title = paste0("EDEH vs. Epsilon", " - ", name), x = "Inequality Aversion (Epsilon)", y = "EDEH") +
      theme_minimal()
    
    # Store the plots in the list
    plots_list[[v]] <- list(ede_epsilon_plot = ede_epsilon_plot)
  }
  
  
} # Close "Distributional Analysis" section
# Close "MODEL" section

# Replication stability
ggplot(
  combined_output %>%
    mutate(replication = paste(dataset_index, replication)), 
  aes(x = replication, y = inhb)
) +
  geom_line() + 
  scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 1000)]) +  # Show every 10th label
  labs(x = "Replication", y = "INHB", title = "INHB across Replications and Imputed datasets") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels to manage space


