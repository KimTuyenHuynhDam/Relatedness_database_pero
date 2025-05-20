# Load necessary libraries
library(openxlsx)
library(tidyverse)
library(broom)
library(lubridate)
library(glue)
library(ggpubr)
library(kinship2)
library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)
library(lubridate)
library(tidyr)
library(cowplot)
library(MASS)
library(pROC)

#  Load the data from the Excel file

data <- read_excel ('./data/Number of mice at birth and weaning -BW.xlsx') %>% distinct() %>% 
  filter(!is.na(`STOCK`) & !is.na(`Mating Number`) & 
           !is.na(`Mice born`) & !is.na(`Mice weaned`) ) %>%
  mutate(`Dead mice` = `Mice born` - `Mice weaned`) %>% 
  mutate(`Mice weaned` = ifelse(`Mice weaned` > `Mice born`, `Mice born`, `Mice weaned`)) 


# Capture initial list of breeding pairs
initial_breeding_pairs <- data %>%
  filter(!is.na(`Mating Number`)) %>%
  distinct(`Mating Number`) %>%
  pull(`Mating Number`)

initial_count <- length(initial_breeding_pairs)

cat("Initial number of breeding pairs: ", initial_count, "\n")

# Calculate total mice born, at weaning, and died for each mating cage
mice_summary <- data %>% 
  group_by(`Mating Number`) %>%
  summarise(
    total_mice_born = sum(`Mice born`, na.rm = TRUE),                 # Sum of mice born
    total_mice_weaned = sum(`Mice weaned`, na.rm = TRUE),             # Sum of mice at weaning
    total_mice_died = sum(`Dead mice`, na.rm = TRUE),                 # Sum of dead mice
    min_birthday = min(Birthday, na.rm = TRUE),                       # Earliest birthday in each cage
    max_birthday = max(Birthday, na.rm = TRUE),                       # Latest birthday in each cage
    mating_period_days = as.numeric(max_birthday - min_birthday),     # Mating period in days
    mating_period_years = as.numeric(max_birthday - min_birthday) / 365.25,  # Mating period in years
    total_litters = n_distinct(Birthday)                              # Total number of unique litters
  ) %>%
  mutate(
    avg_litters_per_month = total_litters / (mating_period_days / 30.44),  # Average litters per month (using average days per month)
    avg_litters_per_year = total_litters / mating_period_years             # Average litters per year
  ) %>%
  ungroup()

# Calculate total_records and total_losses
total_records <- data %>%  nrow()

total_losses <- data %>% 
  filter(`Dead mice` > 0) %>%
  nrow()

total_mice_weaned <- sum(data$`Mice weaned`, na.rm = TRUE)

# Calculate the overall probability of losing pups
prob_loss <- total_losses / total_records

#  Perform the analysis for each "Mating Number" (cage) by counting distinct "Birthday" entries (litters)
results_chi_square <- data %>%
  group_by(`Mating Number`) %>%
  summarise(
    num_litters = n_distinct(Birthday),               # Count of distinct litters by Birthday
    num_losses = sum(`Dead mice` > 0, na.rm = TRUE),  # Count of litters where Dead mice > 0
    expected_losses = num_litters * prob_loss         # Expected losses based on overall loss probability
  ) %>%
  filter(num_litters > 0 & num_losses >= 0 & (num_litters - num_losses) >= 0) %>%  # Ensure valid data for chi-sq test
  rowwise() %>%
  mutate(
    test_p_value = chisq.test(c(num_losses, num_litters - num_losses),
                              p = c(prob_loss, 1 - prob_loss))$p.value,
    
    significance_direction = case_when(
      num_losses > expected_losses ~ "Higher than Expected",
      num_losses < expected_losses ~ "Lower than Expected",
      TRUE ~ "Equal to Expected"
    )
  ) %>%
  ungroup()



#################### Calculate relatedness score of Dam and Sire in each Mating Cage

# Clean column names function
clean_column_names <- function(df) {
  names(df) <- gsub("[^[:alnum:]_]", "", names(df))
  return(df)
}

# Clean merged columns
clean_merged_columns <- function(df) {
  base_names <- stringr::str_replace(names(df), "\\.(x|y)$", "")
  unique_base_names <- unique(base_names)
  
  cleaned_df <- data.frame(row.names = row.names(df))
  
  for (base_name in unique_base_names) {
    matching_cols <- grep(paste0("^", base_name, "\\.(x|y)?$"), names(df), value = TRUE)
    if (base_name %in% names(df)) {
      col_to_keep <- base_name
    } else if (paste0(base_name, ".x") %in% names(df)) {
      col_to_keep <- paste0(base_name, ".x")
    } else {
      col_to_keep <- paste0(base_name, ".y")
    }
    cleaned_df[[base_name]] <- df[[col_to_keep]]
  }
  return(cleaned_df)
}

# Function to clean IDs starting with "0"
clean_ids <- function(df, columns) {
  df <- df %>%
    mutate(across(all_of(columns), ~ ifelse(grepl("^0", .), paste0(substr(., 2, nchar(.)), "00000"), .)))
  return(df)
}



# Read and clean data
pero <- read.xlsx( "./data/Peromyscus.xlsx") %>%
  filter(! STOCK == "EPL") %>% 
  dplyr:: select(1:5) %>% distinct() %>% 
  clean_column_names() %>% 
  mutate(Birthday = as.Date(Birthday, origin = "1899-12-30"), 
         BirthMonth = month(Birthday),
         BirthYear = year(Birthday)) %>%
  filter(BirthYear < 2025)


matingcage <- read.xlsx("./data/Mating Records.xlsx") %>%
  filter(! STOCK == "EPL") %>% 
  dplyr:: select(1:5) %>% distinct() %>% 
  clean_column_names() %>%
  mutate(DateofMating = as.Date(DateofMating, origin = "1899-12-30"))

mating_numbers <- unique(data$`Mating Number`)



max_generations = 5


all_stock <- c("BW")

# Function to get all ancestors for given IDs
get_ancestors <- function(ids, ped_data, max_generations = 10) {
  ancestors <- ids
  for (i in 1:max_generations) {
    new_ancestors <- ped_data %>%
      filter(ID %in% ancestors) %>%
      dplyr:: select(Dam, Sire) %>%
      unlist() %>%
      na.omit() %>%
      unique()
    
    if (length(new_ancestors) == 0) break
    
    ancestors <- unique(c(ancestors, new_ancestors))
  }
  return(ancestors)
}

# calculate relatedness score

for (species in all_stock) {
  species_lower <- tolower(species)
  
  # Process IND data
  IND_ori <- pero %>%
    filter(str_detect(tolower(STOCK), species_lower)) %>% 
    distinct() %>% 
    mutate(across(c(MatingNumber, ID), ~str_remove_all(., regex(species_lower, ignore_case = TRUE)) %>% str_replace_all("[^[:alnum:]]", ""))) %>%
    mutate(Sex = str_replace_all(Sex, "[^[:alnum:]]", "") %>% str_to_upper()) %>%
    mutate(Sex = case_when(Sex %in% c("F", "FEMALE", "FEM") ~ "F", Sex %in% c("M", "MALE") ~ "M", TRUE ~ Sex))
  
  # Process DAMSIRE data
  DAMSIRE_ori <- matingcage %>%
    filter(str_detect(tolower(STOCK), species_lower)) %>%
    mutate(across(c(Dam, Sire, MatingNumber), ~str_remove_all(., regex(species_lower, ignore_case = TRUE)) %>% str_replace_all("[^[:alnum:]]", ""))) %>%
    filter(MatingNumber %in% mating_numbers)  # Filter only the cages from the dataset
  
  IND <- clean_ids(IND_ori, c("ID"))
  DAMSIRE <- clean_ids(DAMSIRE_ori, c("Dam", "Sire"))
  
  DAMSIRE <- DAMSIRE %>%
    mutate(across(c(Dam, Sire), as.numeric)) %>%
    filter(!is.na(Dam), !is.na(Sire)) %>% distinct()
  

  
  
  IND <- IND %>%
    mutate(ID = as.numeric(ID)) %>%
    filter(!is.na(ID)) %>% distinct()
  
  # Create ped_data with filtering for rows where either all three columns are NA or all have values
  ped_data <- IND %>% 
    left_join(DAMSIRE, by = "MatingNumber") %>% distinct()
  
  
  # Initialize results list
  results <- list()
  problematic_mating_cages <- list()
  i <- 0
  
  for (mating_number in mating_numbers) {
    i <- i + 1
    tryCatch({
      mating_pair <- DAMSIRE %>% filter(MatingNumber == mating_number)
      
      if (nrow(mating_pair) > 0) {
        StartingAnimalsOfInterest <- as_vector(mating_pair %>% dplyr:: select(Dam, Sire))
        
        dam_id <- mating_pair$Dam
        sire_id <- mating_pair$Sire
        
        # Debugging: Print the current mating pair
        print(glue("Processing MatingNumber: {mating_number}, Dam: {dam_id}, Sire: {sire_id}"))
        
        # Get all relevant IDs
        relevant_ids <- unique(c(dam_id, sire_id))
        all_ancestors <- get_ancestors(relevant_ids, ped_data, max_generations = 6)
        
        # Debugging: Check for NA in all_ancestors
        if (any(is.na(all_ancestors))) {
          print(glue("NA values found in all_ancestors for MatingNumber: {mating_number}"))
        }
        
        # Include all ancestors in the relevant pedigree
        relevant_pedigree <- IND %>%
          filter(ID %in% all_ancestors) %>%
          distinct(ID, .keep_all = TRUE)  # Ensure there are no duplicate IDs
        
        # Debugging: Check for NA in relevant_pedigree
        if (any(is.na(relevant_pedigree$ID))) {
          print(glue("NA values found in relevant_pedigree IDs for MatingNumber: {mating_number}"))
        }
        
        # Ensure parental information is included
        relevant_pedigree <- relevant_pedigree %>%
          left_join(DAMSIRE, by = "MatingNumber") %>%
          distinct(ID, .keep_all = TRUE) %>%
          mutate(Sex = case_when(
            Sex == "M" ~ 1,
            Sex == "F" ~ 2,
            TRUE ~ 0  # Any other value will be converted to 0
          ))
        
        # Debugging: Check for NA after join
        if (any(is.na(relevant_pedigree$ID))) {
          print(glue("NA values found in relevant_pedigree IDs after join for MatingNumber: {mating_number}"))
        }
        
        # Ensure all parents are in the ID list
        parent_ids <- unique(c(relevant_pedigree$Sire, relevant_pedigree$Dam))
        
        # Add missing parents to the relevant pedigree
        missing_ID <- setdiff(parent_ids, relevant_pedigree$ID)
        
        if (length(missing_ID) > 0) {
          missing_parents <- IND %>%
            filter(ID %in% missing_ID) %>%
            distinct(ID, .keep_all = TRUE) %>%
            left_join(DAMSIRE, by = "MatingNumber") %>%
            distinct(ID, .keep_all = TRUE) %>%
            mutate(Sex = case_when(
              Sex == "M" ~ 1,
              Sex == "F" ~ 2,
              TRUE ~ 0  # Any other value will be converted to 0
            ))
          
          relevant_pedigree <- bind_rows(relevant_pedigree, missing_parents)
        }
        
        # Use fixParents to adjust for any missing parent information
        ped2 <- with(relevant_pedigree, fixParents(id = ID, dadid = Sire, momid = Dam, sex = Sex))
        
        # Create pedigree object
        pedALl <- pedigree(id = ped2$id, dadid = ped2$dadid, momid = ped2$momid, sex = ped2$sex)
        
        # Calculate kinship matrix
        kinAll <- kinship(pedALl)
        
        kinUCI <- kinAll[row.names(kinAll) %in% StartingAnimalsOfInterest, colnames(kinAll) %in% StartingAnimalsOfInterest]
        
        print(round(200 * kinUCI, 1))
        
        Relatedness <- round(200 * kinUCI, 1)[1, 2]
        
        print(Relatedness)
        
        
        # Store the result
        results <- append(results, list(data.frame(
          Order = i,
          MatingNumber = mating_number,
          Dam = dam_id,
          Sire = sire_id,
          Relatedness = Relatedness
        )))
      }
    }, error = function(e) {
      cat("Error in mating cage:", i, "- ", e$message, "\n")
      problematic_mating_cages[[length(problematic_mating_cages) + 1]] <- list(
        cage_number = i,
        mating_cage = StartingAnimalsOfInterest,
        error_message = e$message
      )
    })
    
    cat(i, "..")
  }
  
  # Convert results to a data frame
  results_df <- bind_rows(results)
  
  
  combined_results_chi <- results_chi_square %>%
    left_join(results_df %>% dplyr:: select(MatingNumber, Relatedness), by = c(`Mating Number` = "MatingNumber")) %>%
    left_join(mice_summary, by = "Mating Number") %>%  # Join mice summary
    mutate(
      loss_ratio_litter = num_losses / num_litters,            # Loss ratio based on litters
      loss_ratio_mice = total_mice_died / total_mice_born      # Loss ratio based on total mice
    )
  
  
  
  
  ########################### 
  
  
  # Filter relevant data
  relevant_results <- combined_results_chi %>% 
    filter(Relatedness > 0)
  
  # Capture final list of breeding pairs
  final_breeding_pairs <- relevant_results %>%
    filter(!is.na(`Mating Number`)) %>%
    distinct(`Mating Number`) %>%
    pull(`Mating Number`)
  
  # Get the final count
  final_count <- length(final_breeding_pairs)
  
  cat("Final number of breeding pairs: ", final_count, "\n")
  
  # Calculate the number of removed breeding pairs
  removed_count <- initial_count - final_count
  
  # Identify which pairs were removed
  removed_pairs <- setdiff(initial_breeding_pairs, final_breeding_pairs)
  
  cat("Number of breeding pairs removed: ", removed_count, "\n")
  
  # Prepare data for export
  breeding_summary <- data.frame(
    "Description" = c("Initial Breeding Pairs", "Final Breeding Pairs", "Removed Breeding Pairs"),
    "Count" = c(initial_count, final_count, removed_count),
    "Removed Pairs (if any)" = c(NA, NA, paste(removed_pairs, collapse = ", "))
  )
  
  # Write to Excel
  write.xlsx(breeding_summary, "./database/Breeding_Pairs_Summary.xlsx")
  
 
  ####################################

  
  fit_binomial_glm_loss_ratio <- function(data, x, y, total_col) {
    # Fit the binomial GLM using the loss ratio with appropriate weights
    model <- glm(
      formula = as.formula(paste(y, "~", x)),
      data = data,
      family = binomial(link = "logit"),
      weights = data[[total_col]]
    )
    
    # Generate prediction data for the curve
    new_data <- data.frame(x = seq(min(data[[x]], na.rm = TRUE), max(data[[x]], na.rm = TRUE), length.out = 100))
    names(new_data) <- x
    new_data$y_pred <- predict(model, newdata = new_data, type = "response")
    
    # Fit the null model for McFadden's R²
    null_model <- glm(
      formula = as.formula(paste(y, "~ 1")),
      data = data,
      family = binomial(link = "logit"),
      weights = data[[total_col]]
    )
    
    # Calculate McFadden's R²
    mcfadden_r2 <- 1 - (logLik(model) / logLik(null_model))
    
    # Pearson correlation for additional perspective
    cor_test <- cor.test(data[[x]], data[[y]], use = "complete.obs")
    
    # Return a list of results
    return(list(
      model = model,
      predictions = new_data,
      mcfadden_r2 = round(as.numeric(mcfadden_r2), 4),
      glm_p_value = summary(model)$coefficients[2, 4],
      pearson_r = round(cor_test$estimate, 4),
      pearson_p_value = cor_test$p.value
    ))
  }
  
  add_glm_correlation_loss_ratio <- function(plot, x, y, data, total_col, graph_type) {
    # Fit the model and get predictions
    model_info <- fit_binomial_glm_loss_ratio(data, x, y, total_col)
    model <- model_info$model
    new_data <- model_info$predictions
    
    # Add the GLM prediction curve
    plot <- plot + geom_line(data = new_data, aes_string(x = x, y = "y_pred"), 
                             color = "blue", size = 1.2, linetype = "dashed")
    
    # Determine the appropriate N for annotation
    if (graph_type == "litter") {
      n <- sum(data$num_litters, na.rm = TRUE)
      n_label <- "litters"
      
    } else if (graph_type == "mice") {
      n_born <- sum(data$total_mice_born, na.rm = TRUE)
      n_weaned <- sum(data$total_mice_weaned, na.rm = TRUE)
      n <- paste0(n_born, " (born), ", n_weaned, " (weaned)")
      n_label <- "mice"
      
    } else if (graph_type == "mice_lost") {
      n_born <- sum(data$total_mice_born, na.rm = TRUE)
      n_died <- sum(data$total_mice_died, na.rm = TRUE)
      n <- paste0(n_born, " (born), ", n_died, " (died)")
      n_label <- "mice"
      
    } else {
      # Fallback
      n <- sum(complete.cases(data[[x]], data[[y]]))
      n_label <- "observations"
    }
    
    # Calculate safe annotation positions
    max_x <- max(data[[x]], na.rm = TRUE) * 0.95
    max_y <- max(data[[y]], na.rm = TRUE) * 0.95
    
    # Annotate the plot with R² and correlation
    plot + annotate("text", x = max_x, y = max_y,
                    label = paste0(
                      "GLM McFadden's R² = ", model_info$mcfadden_r2,
                      "\nGLM P-value = ", format.pval(model_info$glm_p_value, digits = 4),
                      "\nPearson R = ", model_info$pearson_r,
                      "\nPearson P-value = ", format.pval(model_info$pearson_p_value, digits = 4),
                      "\nN = ", n, " (", n_label, ")"
                    ),
                    hjust = 1, vjust = 1, size = 4.5, color = "black")
  }
  
  # Plot for Loss Ratio (Litter) vs Relatedness
  p_loss_ratio_litter <- ggplot(relevant_results, aes(x = Relatedness, y = loss_ratio_litter)) +
    geom_point(color = "darkgreen", size = 3) +
    labs(title = "Relatedness vs. Loss Ratio (Litter)",
         x = "Relatedness",
         y = "Loss Ratio (Litter)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Add GLM prediction curve and correlation annotation
  p_loss_ratio_litter <- add_glm_correlation_loss_ratio(
    p_loss_ratio_litter, 
    x = "Relatedness", 
    y = "loss_ratio_litter", 
    data = relevant_results, 
    total_col = "num_litters",
    graph_type = "litter"
  )
  
  # Display the plot
  print(p_loss_ratio_litter)
  
  ggsave(paste0("./database/Relatedness_vs_Loss_Ratio_per_Litter_GLM_", species, ".jpeg"), plot = p_loss_ratio_litter, width = 8, height = 6, device = "jpeg")
  
  
  # Plot for Loss Ratio (Mice) vs Relatedness
  p_loss_ratio_mice <- ggplot(relevant_results, aes(x = Relatedness, y = loss_ratio_mice)) +
    geom_point(color = "purple", size = 3) +
    labs(title = "Relatedness vs. Loss Ratio (Mice)",
         x = "Relatedness",
         y = "Loss Ratio (Mice)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Add GLM prediction curve and correlation annotation
  p_loss_ratio_mice <- add_glm_correlation_loss_ratio(
    p_loss_ratio_mice, 
    x = "Relatedness", 
    y = "loss_ratio_mice", 
    data = relevant_results, 
    total_col = "total_mice_born",
    graph_type = "mice"
  )
  
  # Display the plot
  print(p_loss_ratio_mice)
  
  ggsave(paste0("./database/Relatedness_vs_Loss_Ratio_per_Total_Mice_GLM_", species, ".jpeg"), plot = p_loss_ratio_mice, width = 8, height = 6, device = "jpeg")
  
  
  
  


  ##################################################################################
  ########################## interval from Dateofmating to first litter
  
  # Step 1: Find the first litter date for each Mating Number in the litter data
  first_litter_dates <- data %>%
    group_by(`Mating Number`) %>%
    summarise(first_litter_date = min(Birthday, na.rm = TRUE)) %>%
    ungroup()
  
  # Step 2: Merge DAMSIRE data with first litter dates by `Mating Number`
  merged_data <- DAMSIRE %>%
    left_join(first_litter_dates, by = c("MatingNumber" = "Mating Number"))
  
  # Step 3: Calculate the interval between DateofMating and first_litter_date
  merged_data <- merged_data %>%
    mutate(interval_days = as.numeric(difftime(first_litter_date, DateofMating, units = "days")),
           interval_weeks = interval_days / 7,
           interval_years = interval_days / 365)
  
  
  
  # Step 2: Merge to get Relatedness from results_df, and filter for valid Relatedness values
  merged_data_with_relatedness <- merged_data %>%
    left_join(results_df %>% dplyr:: select(MatingNumber, Relatedness), by = c("MatingNumber")) %>%
    filter(Relatedness > 0 & interval_days > 20 & interval_days < 365)
  
  # Step 3: Define a function to add the Pearson correlation annotation on the plot
  add_correlation_top_right <- function(plot, x_var, y_var, data) {
    # Perform correlation test
    cor_test <- cor.test(data[[x_var]], data[[y_var]], use = "complete.obs")
    
    # Create annotation label
    cor_label <- paste("Pearson R =", round(cor_test$estimate, 3),
                       "\nP-value =", format.pval(cor_test$p.value, digits = 3))
    
    # Place the annotation at the top right of the plot area
    plot + annotate("text", x = Inf, y = Inf, label = cor_label,
                    hjust = 1.1, vjust = 1.5, size = 5, color = "black", parse = FALSE)
  }
  
  
  
  # Step 4: Plot interval_days vs Relatedness
  p_interval_vs_relatedness <- ggplot(merged_data_with_relatedness, aes(x = Relatedness, y = interval_days)) +
    geom_point(color = "blue", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "darkgreen", size = 1.2) +
    labs(title = "Interval (Days) Between Mating and First Litter vs. Relatedness",
         x = "Relatedness",
         y = "Interval (Days)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  # Step 5: Add Pearson correlation annotation
  p_interval_vs_relatedness <- add_correlation_top_right(p_interval_vs_relatedness, "Relatedness", "interval_days", merged_data_with_relatedness)
  
  # Display the plot
  print(p_interval_vs_relatedness)
  
  # Save the plot
  ggsave(glue("./database/{species}-Interval_mating_to_1st_litter_vs_Relatedness.jpeg"), plot = p_interval_vs_relatedness, width = 8, height = 6, device = "jpeg")
  
  
  ##################################################################################
  ########################## Delivery interval
  
  
  # Step 1: Sort data by `Mating Number` and `Birthday` and calculate intervals
  delivery_intervals <- data %>%  left_join(DAMSIRE, by = c("Mating Number" = "MatingNumber")) %>%
    arrange(`Mating Number`, Birthday) %>%
    group_by(`Mating Number`) %>%
    mutate(
      # Calculate interval from DateofMating to the first Birthday
      delivery_interval_days = if_else(row_number() == 1,
                                       as.numeric(difftime(Birthday, first(DateofMating), units = "days")),
                                       as.numeric(difftime(Birthday, lag(Birthday), units = "days")))
    ) %>% 
    ungroup()  %>% rename( "MatingNumber" = "Mating Number" )  %>%
    left_join(results_df %>% dplyr:: select(MatingNumber, Relatedness), by = "MatingNumber") %>%
    filter (Relatedness >0 & delivery_interval_days > 20 & delivery_interval_days < 365 )
  
  
  
  
  
  
  # Step 1: Calculate average delivery interval for each Mating Number
  average_delivery_intervals <- delivery_intervals %>%
    # filter(delivery_interval_days > 20) %>%  # Filter out any intervals <= 0
    group_by(MatingNumber) %>%
    summarise(avg_delivery_interval = mean(delivery_interval_days, na.rm = TRUE)) %>%
    ungroup() %>%  
    left_join(results_df %>% dplyr:: select(MatingNumber, Relatedness), by = "MatingNumber")
  
  
  
  
  full_data <- combined_results_chi %>% rename(MatingNumber = `Mating Number`) %>%
    left_join(average_delivery_intervals, by = 'MatingNumber') %>%
    left_join(DAMSIRE_ori,  by = 'MatingNumber' ) %>% clean_merged_columns()
  
  write.xlsx( full_data, glue("./database/mice at born and weaning - {species}- full data.xlsx"))
  
  
  
  # Step 3: Define a function to add Pearson correlation annotation on the plot
  add_correlation_top_right <- function(plot, x_var, y_var, data) {
    # Perform correlation test
    cor_test <- cor.test(data[[x_var]], data[[y_var]], use = "complete.obs")
    
    # Create annotation label
    cor_label <- paste("Pearson R =", round(cor_test$estimate, 3),
                       "\nP-value =", format.pval(cor_test$p.value, digits = 3))
    
    # Place the annotation at the top right of the plot area
    plot + annotate("text", x = Inf, y = Inf, label = cor_label,
                    hjust = 1.1, vjust = 1.5, size = 5, color = "black", parse = FALSE)
  }
  
  
  
  # Step 4: Plot avg_delivery_interval vs Relatedness
  p_avg_delivery_interval_vs_relatedness <- ggplot(average_delivery_intervals, aes(x = Relatedness, y = avg_delivery_interval)) +
    geom_point(color = "blue", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "darkgreen", size = 1.2) +
    labs(title = "Average Delivery Interval vs. Relatedness",
         x = "Relatedness",
         y = "Average Delivery Interval (Days)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  # Step 5: Add Pearson correlation annotation
  p_avg_delivery_interval_vs_relatedness <- add_correlation_top_right(p_avg_delivery_interval_vs_relatedness, "Relatedness", "avg_delivery_interval", average_delivery_intervals)
  
  # Display the plot
  print(p_avg_delivery_interval_vs_relatedness)
  
  
  # Save the plot
  ggsave(glue("./database/{species} Avg_Delivery_Interval_vs_Relatedness.jpeg"), plot = p_avg_delivery_interval_vs_relatedness, width = 8, height = 6, device = "jpeg")
  
  
  #####################
  
  
  
  ##################################################################################
  ########################## Relatedness groups
  
  relevant_results <- relevant_results %>% mutate(
    Relatedness_Group = cut(Relatedness, 
                            breaks = c(-Inf, 25, 30, 35, 40, Inf), 
                            labels = c("<24.9", "25-29.9", "30-34.9", "35-39.9", "≥40"),
                            include.lowest = TRUE)
  )
  
  # Ensure 'Relatedness_Group' is a factor to maintain column order
  relevant_results$Relatedness_Group <- factor(relevant_results$Relatedness_Group, levels = c("<24.9", "25-29.9", "30-34.9", "35-39.9", "≥40"))
  
  
  ####################
  
  # Step 1: Transform the data into long format
  loss_ratio_long <- relevant_results %>% filter(num_litters > 2) %>%
    dplyr:: select(Relatedness_Group, loss_ratio_litter, loss_ratio_mice) %>%
    pivot_longer(
      cols = c(loss_ratio_litter, loss_ratio_mice),
      names_to = "Loss_Type",
      values_to = "Loss_Ratio"
    )
  
  # Step 2: Define a function to create tables with padded columns
  create_loss_ratio_table <- function(data, loss_type) {
    # Filter data for the specified loss type
    filtered_data <- data %>%
      filter(Loss_Type == loss_type) %>%
      dplyr:: select(Relatedness_Group, Loss_Ratio)
    
    # Split data by Relatedness_Group
    split_data <- filtered_data %>%
      group_split(Relatedness_Group)
    
    # Calculate the maximum length among all groups
    max_length <- max(sapply(split_data, nrow))
    
    # Pad each group's Loss_Ratio column to the max_length and convert to a data frame
    padded_data <- map(split_data, ~ {
      tibble(Loss_Ratio = c(.x$Loss_Ratio, rep(NA, max_length - nrow(.x))))
    })
    
    # Bind columns and set names as Relatedness_Group levels
    result <- bind_cols(padded_data) %>%
      setNames(levels(data$Relatedness_Group))
    
    return(result)
  }
  
  # Table for loss_ratio_litter
  loss_ratio_litter_table <- create_loss_ratio_table(loss_ratio_long, "loss_ratio_litter")
  
  
  write.xlsx(loss_ratio_litter_table, "./database/loss_ratio_litter_table(litters gt 2).xlsx")
  
  # Table for loss_ratio_mice
  loss_ratio_mice_table <- create_loss_ratio_table(loss_ratio_long, "loss_ratio_mice")
  
  
  write.xlsx(loss_ratio_mice_table, "./database/loss_ratio_mice_table(litters gt 2).xlsx")
  ###############
  
  # Ensure `Relatedness_Group` is ordered from smallest to largest
  loss_ratio_litter_long <- loss_ratio_litter_table %>%
    pivot_longer(
      cols = everything(),
      names_to = "Relatedness_Group",
      values_to = "Loss_Ratio"
    ) %>%
    drop_na() %>%  # Retain non-missing values
    mutate(Relatedness_Group = factor(Relatedness_Group, levels = c("<24.9", "25-29.9", "30-34.9", "35-39.9", "≥40")))
  
  # Calculate mean and standard error for each Relatedness_Group for the litter data
  plot_data_litter <- loss_ratio_litter_long %>%
    group_by(Relatedness_Group) %>%
    summarise(
      mean_loss_ratio = mean(Loss_Ratio),
      se_loss_ratio = sd(Loss_Ratio) / sqrt(n())
    )
  
  # Perform ANOVA on the data and extract the p-value
  anova_result <- aov(Loss_Ratio ~ Relatedness_Group, data = loss_ratio_litter_long)
  anova_p_value <- summary(anova_result)[[1]][["Pr(>F)"]][1]  # Extract the p-value from ANOVA
  
  # Plot for Loss Ratio Litter with only the ANOVA p-value
  p_litter <- ggplot(plot_data_litter, aes(x = Relatedness_Group, y = mean_loss_ratio, fill = Relatedness_Group)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = mean_loss_ratio - se_loss_ratio, ymax = mean_loss_ratio + se_loss_ratio), width = 0.2) +
    scale_fill_manual(values = c("blue", "red", "green", "purple", "orange", "black")) +
    labs(title = "Loss Ratio per Litter (Litters >2)",
         x = "Relatedness (%)",
         y = "% Matings with Pups Loss") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    # Add the ANOVA p-value as an annotation on the plot
    annotate("text", x = 3, y = max(plot_data_litter$mean_loss_ratio) * 1.1,
             label = paste("ANOVA p-value =", format.pval(anova_p_value, digits = 3)),
             size = 5, color = "black")
  
  # Print and save the litter plot
  print(p_litter)
  ggsave("./database/Loss_Ratio_Litter_Relatedness_Comparison_Litters_over_2.jpeg", plot = p_litter, width = 5, height = 8, device = "jpeg")
  
  ##########################
  
  
  # Ensure `Relatedness_Group` is ordered from smallest to largest
  loss_ratio_mice_long <- loss_ratio_mice_table %>%
    pivot_longer(
      cols = everything(),
      names_to = "Relatedness_Group",
      values_to = "Loss_Ratio"
    ) %>%
    drop_na() %>%  # Retain non-missing values
    mutate(Relatedness_Group = factor(Relatedness_Group, levels = c("<24.9", "25-29.9", "30-34.9", "35-39.9", "≥40")))
  
  # Calculate mean and standard error for each Relatedness_Group for the mice data
  plot_data_mice <- loss_ratio_mice_long %>%
    group_by(Relatedness_Group) %>%
    summarise(
      mean_loss_ratio = mean(Loss_Ratio),
      se_loss_ratio = sd(Loss_Ratio) / sqrt(n())
    )
  
  # Perform ANOVA on the data and extract the p-value
  anova_result_mice <- aov(Loss_Ratio ~ Relatedness_Group, data = loss_ratio_mice_long)
  anova_p_value_mice <- summary(anova_result_mice)[[1]][["Pr(>F)"]][1]  # Extract the p-value from ANOVA
  
  # Plot for Loss Ratio Mice with only the ANOVA p-value
  p_mice <- ggplot(plot_data_mice, aes(x = Relatedness_Group, y = mean_loss_ratio, fill = Relatedness_Group)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = mean_loss_ratio - se_loss_ratio, ymax = mean_loss_ratio + se_loss_ratio), width = 0.2) +
    scale_fill_manual(values = c("blue", "red", "green", "purple", "orange", "black")) +
    labs(title = "Loss Ratio of Mice (Litters >2)",
         x = "Relatedness (%)",
         y = "Loss Ratio (Mice Lost / Mice Born)") +  # Updated y-axis label
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    # Add the ANOVA p-value as an annotation on the plot
    annotate("text", x = 3, y = max(plot_data_mice$mean_loss_ratio) * 1.1,
             label = paste("ANOVA p-value =", format.pval(anova_p_value_mice, digits = 3)),
             size = 5, color = "black")
  
  # Print and save the mice plot
  print(p_mice)
  ggsave("./database/Loss_Ratio_Mice_Relatedness_Comparison_Litters_over_2.jpeg", plot = p_mice, width = 5, height = 8, device = "jpeg")
  
  ##############
  # Add Year column and define the range of years
  full_data_period <- full_data %>% 
    merge(merged_data, by = "MatingNumber") %>%  
    clean_merged_columns() %>% distinct() %>% 
    mutate(Year = year(DateofMating)) %>%
    filter (Year >=1963 )
  
  
  
  
  min_year <- min(full_data_period$Year, na.rm = TRUE)
  max_year <- 2024
  
  
  
  
  #########################################
  ###################
  
  # Summarize data by gender and calculate necessary metrics
  gender_loss_ratios <- data %>%
    group_by(`Mating Number`) %>%
    summarise(
      total_mice_born = sum(`Mice born`, na.rm = TRUE),
      total_F_weaned = sum(`F at weaning`, na.rm = TRUE),
      total_M_weaned = sum(`M at weaning`, na.rm = TRUE),
      num_litters = n_distinct(Birthday),
      num_losses = sum(`Dead mice` > 0, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      total_F_born = total_mice_born / 2,
      total_M_born = total_mice_born / 2,
      F_loss_ratio = pmax(0, (total_F_born - total_F_weaned) / total_F_born),
      M_loss_ratio = pmax(0, (total_M_born - total_M_weaned) / total_M_born)
    )
  
  # Merge the datasets
  merged_data <- gender_loss_ratios %>%
    rename(`MatingNumber` = `Mating Number`) %>%
    inner_join(full_data_period, by = "MatingNumber") %>%
    clean_merged_columns() %>%
    distinct()
  
  # Compute correlations and McFadden's R²
  compute_continuous_correlations_by_gender <- function(data, genders, gender_labels, min_year, max_year) {
    combined_results <- data.frame()
    
    for (i in seq_along(genders)) {
      gender_var <- genders[i]
      gender_label <- gender_labels[i]
      
      gender_results <- data.frame(
        Year = integer(),
        Pearson_R = numeric(),
        Pearson_P = numeric(),
        McFadden_R2 = numeric(),
        McFadden_P = numeric(),
        Gender = character()
      )
      
      for (year in seq(min_year, max_year)) {
        subset_data <- data %>%
          filter(Year >= min_year & Year <= year, !is.na(Relatedness), !is.na(!!sym(gender_var)))
        
        if (nrow(subset_data) >= 5) {
          # Pearson Correlation
          cor_test <- cor.test(subset_data$Relatedness, subset_data[[gender_var]], use = "complete.obs")
          
          # GLM for McFadden's R²
          glm_model <- glm(
            formula = as.formula(paste(gender_var, "~ Relatedness")),
            family = binomial(link = "logit"),
            data = subset_data
          )
          
          # Null model for comparison
          null_model <- glm(
            formula = as.formula(paste(gender_var, "~ 1")),
            family = binomial(link = "logit"),
            data = subset_data
          )
          
          # Compute McFadden's R²
          mcfadden_r2 <- 1 - (logLik(glm_model) / logLik(null_model))
          
          # Likelihood ratio test for McFadden's p-value
          lrt <- anova(null_model, glm_model, test = "Chisq")
          mcfadden_p_value <- ifelse(nrow(lrt) >= 2, lrt$`Pr(>Chi)`[2], NA)
          
          # Append results
          gender_results <- rbind(
            gender_results,
            data.frame(
              Year = year,
              Pearson_R = cor_test$estimate,
              Pearson_P = cor_test$p.value,
              McFadden_R2 = mcfadden_r2,
              McFadden_P = mcfadden_p_value,
              Gender = gender_label
            )
          )
        }
      }
      
      combined_results <- rbind(combined_results, gender_results)
    }
    
    return(combined_results)
  }
  
  # Define genders and labels
  genders <- c("F_loss_ratio", "M_loss_ratio")
  gender_labels <- c("Female", "Male")
  
  # Compute continuous correlations
  continuous_correlations_gender <- compute_continuous_correlations_by_gender(
    data = merged_data,
    genders = genders,
    gender_labels = gender_labels,
    min_year = min(merged_data$Year, na.rm = TRUE),
    max_year = max(merged_data$Year, na.rm = TRUE)
  )
  
  # Plot the four plots
  pearson_plot <- ggplot(continuous_correlations_gender, aes(x = Year, y = Pearson_R, color = Gender)) +
    geom_line(size = 1) +
    labs(
      y = "Pearson R",
      title = "Pearson R by Gender Over Time"
    ) +
    theme_minimal()
  
  pearson_p_plot <- ggplot(continuous_correlations_gender, aes(x = Year, y = -log10(Pearson_P), color = Gender)) +
    geom_line(size = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +
    labs(
      y = "-log10(Pearson P-value)",
      title = "Pearson P-value by Gender Over Time"
    ) +
    theme_minimal()
  
  
  mcfadden_plot <- ggplot(continuous_correlations_gender, aes(x = Year, y = McFadden_R2, color = Gender)) +
    geom_line(size = 1) +
    labs(
      y = "McFadden's R²",
      title = "McFadden's R² by Gender Over Time"
    ) +
    theme_minimal()
  
  mcfadden_p_plot <- ggplot(continuous_correlations_gender, aes(x = Year, y = -log10(McFadden_P), color = Gender)) +
    geom_line(size = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +
    labs(
      y = "-log10(McFadden P-value)",
      title = "McFadden P-value by Gender Over Time"
    ) +
    theme_minimal()
  
  # Combine into one large plot, stacked vertically
  final_plot <- cowplot::plot_grid(
    pearson_plot, 
    pearson_p_plot,
    mcfadden_plot,
    mcfadden_p_plot,
    ncol = 1,
    labels = c("a", "b", "c", "d"),
    label_size = 14
  )
  
  # Add a title to the whole plot
  title <- ggdraw() + draw_label("Gender-Based Loss Ratio Analysis Over Time", fontface = 'bold', size = 15)
  combined_final_plot <- plot_grid(title, final_plot, ncol = 1, rel_heights = c(0.1, 1))
  
  # Save the final plot
  ggsave("./database/Gender_Based_Loss_Ratio_Analysis_Over_Time_Stacked.jpeg", plot = combined_final_plot, width = 10, height = 14)
  
  # Display the final plot
  print(combined_final_plot)
  
  
  
 
 
  
  ##########################
  
  # Function to compute continuous correlations
  compute_continuous_correlations <- function(data, y_var, min_year, max_year) {
    results <- data.frame()
    
    for (year in seq(min_year, max_year)) {
      subset_data <- data %>% filter(Year >= min_year & Year <= year)
      if (nrow(subset_data) >= 5) {  # Ensure at least 5 rows for a valid correlation
        cor_test <- cor.test(subset_data$Relatedness, subset_data[[y_var]], use = "complete.obs")
        results <- rbind(
          results,
          data.frame(
            Year = year,
            R = round(cor_test$estimate, 3),
            P = format.pval(cor_test$p.value, digits = 3),
            y_var = y_var
          )
        )
      }
    }
    return(results)
  }
  
  compute_continuous_correlations <- function(data, y_var, min_year, max_year) {
    results <- data.frame(Year = integer(), Pearson_R = numeric(), P_value = numeric())
    
    for (year in min_year:max_year) {
      subset_data <- data %>% filter(Year <= year)
      
      if (nrow(subset_data) > 2) {  # Ensure enough data points for correlation
        cor_test <- cor.test(subset_data$Relatedness, subset_data[[y_var]], use = "complete.obs")
        results <- rbind(
          results,
          data.frame(
            Year = year,
            Pearson_R = cor_test$estimate,
            P_value = cor_test$p.value
          )
        )
      } else {
        results <- rbind(
          results,
          data.frame(
            Year = year,
            Pearson_R = NA,
            P_value = NA
          )
        )
      }
    }
    return(results)
  }
  
  
  
  plot_continuous_correlation_stacked <- function(data, y_label, title_suffix, species, save_path) {
    # Ensure all years are represented
    all_years <- seq(min(data$Year, na.rm = TRUE), max(data$Year, na.rm = TRUE))
    
    # Fill in missing years with NA for consistency
    data <- data %>%
      complete(Year = all_years, fill = list(Pearson_R = NA, P_value = NA))
    
    # Compute -log10(p-value)
    data <- data %>%
      mutate(logP_transformed = ifelse(P_value > 0, -log10(P_value), NA))
    
    # Pearson R plot
    pearson_plot <- ggplot(data, aes(x = Year)) +
      geom_line(aes(y = Pearson_R, color = "Pearson R"), size = 1) +
      geom_point(aes(y = Pearson_R, color = "Pearson R"), size = 2) +
      labs(
        # title = paste("Pearson R vs Year (", title_suffix, ", Species:", species, ")"),
        x = NULL,  # Remove x-axis title here to avoid redundancy
        y = "Pearson R",
        color = NULL
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),  # Hide x-axis text for the top plot
        axis.title.x = element_blank(),  # Hide x-axis title for the top plot
        legend.position = "top"
      )
    
    # P-value plot
    pvalue_plot <- ggplot(data, aes(x = Year)) +
      geom_line(aes(y = logP_transformed, color = "-log10(P-value)"), linetype = "dashed", size = 1) +
      geom_point(aes(y = logP_transformed, color = "-log10(P-value)"), size = 2, shape = 1) +
      geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +  # Significance threshold
      labs(
        x = "Year",
        y = "-log10(P-value)",
        color = NULL
      ) +
      theme_minimal() +
      theme(
        plot.title = element_blank(),  # No title for the bottom plot
        axis.title.x = element_text(size = 12),
        legend.position = "top"
      )
    
    # Combine both plots vertically
    combined_plot <- cowplot::plot_grid(
      pearson_plot, 
      pvalue_plot, 
      ncol = 1, 
      align = "v", 
      rel_heights = c(1, 1)  # Equal height for both plots
    )
    
    # Add an overarching title for the combined plot
    combined_title <- ggdraw() + 
      draw_label(
        paste(
          "Correlation of Relatedness vs", y_label, 
          "(", title_suffix, ", Species:", species, ")\n",
          "(Period of", min_year, "to Each Year)"
        ),
        fontface = "bold",
        size = 14,
        hjust = 0.5
      )
    
    
    # Combine title and plots
    final_plot <- plot_grid(combined_title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))
    
    # Save the final plot
    ggsave(
      filename = save_path,
      plot = final_plot,
      width = 10,
      height = 8,  # Slightly taller to accommodate stacked plots
      device = "jpeg"
    )
    
    return(final_plot)
  }
  
  # Define variables and their labels
  variables <- list(
    "total_mice_born" = "Total Mice Born",
    "total_mice_weaned" = "Total Mice Weaned",
    "loss_ratio_litter" = "Loss Ratio (By Litters)",
    "loss_ratio_mice" = "Loss Ratio (Total Mice)",
    "avg_litters_per_month" = "Average Litters Per Month",
    "avg_delivery_interval" = "Average Delivery Interval (Days)",
    "interval_days" = "Interval Days (Date of Mating to First Litter)"
  )
  
  
  # Create output folder
  output_folder <- "./database/continuous_year_period"
  dir.create(output_folder, showWarnings = FALSE)
  
  
  
  
  # Generate stacked correlation graphs
  for (litters in c(0, 1, 2, 3, 4)) {
    filtered_data <- full_data_period %>% filter(num_litters > litters) %>%
      filter(Relatedness >0 ) 
    for (var in names(variables)) {
      y_label <- variables[[var]]
      
      continuous_correlations <- compute_continuous_correlations(
        data = filtered_data,
        y_var = var,
        min_year = min_year,
        max_year = max_year
      )
      
      plot <- plot_continuous_correlation_stacked(
        data = continuous_correlations,
        y_label = y_label,
        title_suffix = paste("Litters >", litters),
        species = species,
        save_path = paste0("./database/continuous_year_period/Stacked_Correlation_", gsub(" ", "_", y_label), "_Litters_Gt_", litters, "_", species, ".jpeg")
      )
      print(plot)
    }
  }

  
  ########################################################
  # Function to compute continuous correlations and GLM R-squared
  compute_continuous_glm_correlations <- function(data, y_var, total_col, min_year, max_year) {
    results <- data.frame(Year = integer(), 
                          Pearson_R = numeric(), 
                          Pearson_P = numeric(),
                          McFadden_R2 = numeric(), 
                          GLM_P_value = numeric(),
                          McFadden_P_value = numeric(),
                          AUC = numeric())
    
    for (year in seq(min_year, max_year)) {
      subset_data <- data %>% filter(Year >= min_year & Year <= year)
      
      if (nrow(subset_data) >= 5) {
        
        if (length(unique(subset_data$Relatedness)) > 1 & length(unique(subset_data[[y_var]])) > 1) {
          
          # Drop NA values to ensure same data size for both models
          subset_data <- na.omit(subset_data[, c("Relatedness", y_var, total_col)])
          
          # Pearson Correlation
          cor_test <- tryCatch({
            cor.test(subset_data$Relatedness, subset_data[[y_var]], use = "complete.obs")
          }, error = function(e) {
            list(estimate = NA, p.value = NA)
          })
          
          # Fit the Binomial GLM
          glm_model <- glm(
            formula = as.formula(paste(y_var, "~ Relatedness")),
            family = binomial(link = "logit"),
            weights = subset_data[[total_col]],
            data = subset_data
          )
          
          # Calculate McFadden's R²
          null_model <- glm(
            formula = as.formula(paste(y_var, "~ 1")),
            family = binomial(link = "logit"),
            weights = subset_data[[total_col]],
            data = subset_data
          )
          
          # Apply the LRT safely
          lrt <- tryCatch({
            anova(null_model, glm_model, test = "Chisq")
          }, error = function(e) {
            data.frame(`Pr(>Chi)` = c(NA, NA))
          })
          
          # Extract McFadden's p-value safely
          mcfadden_p_value <- ifelse(nrow(lrt) >= 2, lrt$`Pr(>Chi)`[2], NA)
          
          # AUC Calculation
          binary_response <- ifelse(subset_data[[y_var]] > 0, 1, 0)
          auc_value <- tryCatch({
            auc(roc(binary_response, predict(glm_model, type = "response")))
          }, error = function(e) NA)
          
          # GLM P-value
          if (nrow(summary(glm_model)$coefficients) > 1) {
            r2_mcfadden <- 1 - (logLik(glm_model) / logLik(null_model))
            glm_p_value <- summary(glm_model)$coefficients[2, 4]
          } else {
            r2_mcfadden <- NA
            glm_p_value <- NA
            mcfadden_p_value <- NA
          }
          
          # Store results
          results <- rbind(
            results,
            data.frame(
              Year = year,
              Pearson_R = as.numeric(cor_test$estimate),
              Pearson_P = as.numeric(cor_test$p.value),
              McFadden_R2 = as.numeric(r2_mcfadden),
              GLM_P_value = as.numeric(glm_p_value),
              McFadden_P_value = as.numeric(mcfadden_p_value),
              AUC = auc_value
            )
          )
        }
      }
    }
    return(results)
  }
  
  
  plot_glm_correlation_over_time <- function(data, y_label, title_suffix, species, save_path) {
    all_years <- seq(min(data$Year, na.rm = TRUE), max(data$Year, na.rm = TRUE))
    
    # Fill in missing years with NA for consistency
    data <- data %>%
      complete(Year = all_years, fill = list(Pearson_R = NA, McFadden_R2 = NA, Pearson_P = NA, McFadden_P_value = NA))
    
    # Plot Pearson R and McFadden's R²
    correlation_plot <- ggplot(data, aes(x = Year)) +
      geom_line(aes(y = Pearson_R, color = "Pearson R"), size = 1.2) +
      geom_point(aes(y = Pearson_R, color = "Pearson R"), size = 2) +
      geom_line(aes(y = McFadden_R2, color = "McFadden's R²"), linetype = "dashed", size = 1.2) +
      geom_point(aes(y = McFadden_R2, color = "McFadden's R²"), size = 2) +
      labs(y = "Correlation Coefficients",
           color = NULL) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "top"
      )
    
    # Plot P-values
    pvalue_plot <- ggplot(data, aes(x = Year)) +
      geom_line(aes(y = -log10(Pearson_P), color = "-log10(Pearson P-value)"), linetype = "solid", size = 1.2) +
      geom_point(aes(y = -log10(Pearson_P), color = "-log10(Pearson P-value)"), size = 2) +
      geom_line(aes(y = -log10(McFadden_P_value), color = "-log10(McFadden P-value)"), linetype = "dotted", size = 1.2) +
      geom_point(aes(y = -log10(McFadden_P_value), color = "-log10(McFadden P-value)"), size = 2) +
      labs(y = "-log10(P-value)",
           color = NULL) +
      geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "top"
      )
    
    # Combine both plots vertically
    combined_plot <- cowplot::plot_grid(
      correlation_plot, 
      pvalue_plot, 
      ncol = 1, 
      align = "v", 
      rel_heights = c(1, 1)
    )
    
    # Add a title
    combined_title <- ggdraw() + 
      draw_label(
        paste(
          "Correlation and P-values of Relatedness vs", y_label, 
          "(", title_suffix, ", Species:", species, ")\n",
          "(Period of 1963 to Each Year)"
        ),
        fontface = "bold",
        size = 14,
        hjust = 0.5
      )
    
    # Merge title and plots
    final_plot <- plot_grid(combined_title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))
    
    print(final_plot)
    # Save the plot
    ggsave(filename = save_path, plot = final_plot, width = 10, height = 12, device = "jpeg")
    
    return(final_plot)
  }
  
  
  # Loss Ratio by Litters
  litters_correlation_data <- compute_continuous_glm_correlations(full_data_period, "loss_ratio_litter", "num_litters", 1963, 2024)
  plot_glm_correlation_over_time(litters_correlation_data, "Loss Ratio (By Litters)", "Litters > 0", "BW", "./database/Loss_Ratio_Litters_GLM_Over_Time.jpeg")
  
  # Loss Ratio by Mice
  mice_correlation_data <- compute_continuous_glm_correlations(full_data_period, "loss_ratio_mice", "total_mice_born", 1963, 2024)
  
  plot_glm_correlation_over_time(mice_correlation_data, "Loss Ratio (By Mice)", "Mice > 0", "BW", "./database/Loss_Ratio_Mice_GLM_Over_Time.jpeg")
  
  
  
  
  
}
  