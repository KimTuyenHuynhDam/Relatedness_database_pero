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
library(dplyr)
library(lubridate)
library(tidyr)
library(cowplot)


# Step 1: Load the data from the Excel file

data <- read_excel ('./data/Number of mice at birth and weaning -BW.xlsx') %>% distinct() %>% 
  filter(!is.na(`STOCK`) & !is.na(`Mating Number`) & 
           !is.na(`Mice born`) & !is.na(`Mice weaned`) ) %>%
  mutate(`Dead mice` = `Mice born` - `Mice weaned`) %>% 
  mutate(`Mice weaned` = ifelse(`Mice weaned` > `Mice born`, `Mice born`, `Mice weaned`)) 
  




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

# Step 2: Calculate total_records and total_losses
total_records <- data %>%  nrow()

total_losses <- data %>% 
  filter(`Dead mice` > 0) %>%
  nrow()

total_mice_wean <- sum(data$`Mice weaned`, na.rm = TRUE)

# Calculate the overall probability of losing pups
prob_loss <- total_losses / total_records

# Step 3: Perform the analysis for each "Mating Number" (cage) by counting distinct "Birthday" entries (litters)
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
  
  
   
 
  ########################### producing graphs
 
  
  # Filter relevant data
  relevant_results <- combined_results_chi %>% 
    filter(Relatedness > 0 )
  
  # Helper function to annotate correlation on the plot
  add_correlation <- function(plot, x, y, data) {
    cor_results <- cor.test(data[[x]], data[[y]], use = "complete.obs")
    plot + annotate("text", x = min(data[[x]], na.rm = TRUE), 
                    y = max(data[[y]], na.rm = TRUE),
                    label = paste0("Pearson R = ", round(cor_results$estimate, 4), 
                                   "\nP-value = ", format.pval(cor_results$p.value, digits = 4)),
                    hjust = 0, vjust = 1, size = 5, color = "black")
  }
  
  # 1. Relatedness vs Total Mice Born
  p_mice_born <- ggplot(relevant_results, aes(x = Relatedness, y = total_mice_born)) +
    geom_point(color = "darkblue", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "red", size = 1.2) +
    labs(title = paste("Relatedness vs. Total Mice Born (Species:", species, ")"),
         x = "Relatedness",
         y = "Total Mice Born") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  # Add correlation annotation
  p_mice_born <- add_correlation(p_mice_born, "Relatedness", "total_mice_born", relevant_results)
  print(p_mice_born)  # Print the graph
  
  ggsave(paste0("./database/Relatedness_vs_Total_Mice_Born_", species, ".jpeg"), plot = p_mice_born, width = 8, height = 6, device = "jpeg")
  
  
  
  
  
  # 2. Relatedness vs Loss Ratio (Litter)
  p_loss_ratio_litter <- ggplot(relevant_results, aes(x = Relatedness, y = loss_ratio_litter)) +
    geom_point(color = "darkgreen", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", size = 1.2) +
    labs(title = paste("Relatedness vs. Loss Ratio per Litter (Species:", species, ")"),
         x = "Relatedness",
         y = "Loss Ratio (Litter)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  # Add correlation annotation
  p_loss_ratio_litter <- add_correlation(p_loss_ratio_litter, "Relatedness", "loss_ratio_litter", relevant_results)
  print(p_loss_ratio_litter)  # Print the graph
  ggsave(paste0("./database/Relatedness_vs_Loss_Ratio_per_Litter_", species, ".jpeg"), plot = p_loss_ratio_litter, width = 8, height = 6, device = "jpeg")
  
  # 3. Relatedness vs Loss Ratio (Total Mice)
  p_loss_ratio_mice <- ggplot(relevant_results, aes(x = Relatedness, y = loss_ratio_mice)) +
    geom_point(color = "purple", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "orange", size = 1.2) +
    labs(title = paste("Relatedness vs. Loss Ratio per Total Mice (Species:", species, ")"),
         x = "Relatedness",
         y = "Loss Ratio (Total Mice)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  # Add correlation annotation
  p_loss_ratio_mice <- add_correlation(p_loss_ratio_mice, "Relatedness", "loss_ratio_mice", relevant_results)
  print(p_loss_ratio_mice)  # Print the graph
  ggsave(paste0("./database/Relatedness_vs_Loss_Ratio_per_Total_Mice_", species, ".jpeg"), plot = p_loss_ratio_mice, width = 8, height = 6, device = "jpeg")
  
  # 4. Total Mice Born vs Total Mice Lost
  p_mice_born_vs_lost <- ggplot(relevant_results, aes(x = total_mice_born, y = total_mice_died)) +
    geom_point(color = "darkred", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "green", size = 1.2) +
    labs(title = paste("Total Mice Born vs. Total Mice Lost (Species:", species, ")"),
         x = "Total Mice Born",
         y = "Total Mice Lost") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  # Add correlation annotation
  p_mice_born_vs_lost <- add_correlation(p_mice_born_vs_lost, "total_mice_born", "total_mice_died", relevant_results)
  print(p_mice_born_vs_lost)  # Print the graph
  ggsave(paste0("./database/Total_Mice_Born_vs_Total_Mice_Lost_", species, ".jpeg"), plot = p_mice_born_vs_lost, width = 8, height = 6, device = "jpeg")
  
  
  # 5. Plot: Relatedness vs Total Mice Weaned
  p_total_mice_weaned <- ggplot(relevant_results, aes(x = Relatedness, y = total_mice_weaned)) +
    geom_point(color = "blue", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "darkgreen", size = 1.2) +
    labs(title = paste("Relatedness vs. Total Mice Weaned (Species:", species, ")"),
         x = "Relatedness",
         y = "Total Mice Weaned") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title

  # Add correlation annotation
  p_total_mice_weaned <- add_correlation(p_total_mice_weaned, "Relatedness", "total_mice_weaned", relevant_results)

  # Print the plot and save as JPEG
  print(p_total_mice_weaned)
  ggsave(paste0("./database/Relatedness_vs_Total_Mice_Weaned_", species, ".jpeg"), plot = p_total_mice_weaned, width = 8, height = 6, device = "jpeg")

  ########################
  # Calculate the weaned ratio (total mice weaned / total mice born)
  ratio_data <- relevant_results %>%  filter(num_litters > 2)  %>%
    mutate(weaned_ratio = total_mice_weaned / total_mice_born)
  
  # Define function to add correlation annotation to the plot
  add_correlation2 <- function(plot, x_var, y_var, data) {
    cor_test <- cor.test(data[[x_var]], data[[y_var]], use = "complete.obs")
    cor_label <- paste("Pearson R =", round(cor_test$estimate, 3), 
                       "\nP-value =", format.pval(cor_test$p.value, digits = 3))
    plot + annotate("text", x = -Inf, y = -Inf, label = cor_label, 
                    hjust = -0.1, vjust = -1, size = 5, color = "black", parse = FALSE)
  }
  # Plot: Relatedness vs Weaned Ratio
  p_weaned_ratio <- ggplot(ratio_data, aes(x = Relatedness, y = weaned_ratio)) +
    geom_point(color = "blue", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "darkgreen", size = 1.2) +
    labs(title = paste("Relatedness vs. Weaned Ratio (Species:", species, ")"),
         x = "Relatedness",
         y = "Weaned Ratio (Mice Weaned / Mice Born)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  # Add correlation annotation
  p_weaned_ratio <- add_correlation2(p_weaned_ratio, "Relatedness", "weaned_ratio", ratio_data)
  
  # Print the plot and save as JPEG
  print(p_weaned_ratio)
  ggsave(paste0("./database/Relatedness_vs_Weaned_Ratio_", species, "(litters more than 2).jpeg"), plot = p_weaned_ratio, width = 8, height = 6, device = "jpeg")
  
  #####
  

  # Calculate average litter size at birth and at weaning for each mating cage
  ratio_data2 <- relevant_results %>% 
   mutate(
      avg_litter_size_birth = total_mice_born / num_litters,
      avg_litter_size_wean = total_mice_weaned / num_litters
   )
  
  # Define function to add dynamic correlation annotation based on Pearson R
  add_dynamic_correlation <- function(plot, x_var, y_var, data) {
    cor_test <- cor.test(data[[x_var]], data[[y_var]], use = "complete.obs")
    cor_label <- paste("Pearson R =", round(cor_test$estimate, 3),
                       "\nP-value =", format.pval(cor_test$p.value, digits = 3))

    # Determine annotation position based on Pearson R
    if (cor_test$estimate < 0) {
      plot + annotate("text", x = -Inf, y = -Inf, label = cor_label,
                      hjust = -0.1, vjust = -1, size = 5, color = "black", parse = FALSE)
    } else {
      plot + annotate("text", x = -Inf, y = Inf, label = cor_label,
                      hjust = -0.1, vjust = 1.5, size = 5, color = "black", parse = FALSE)
    }
  }
  
  # Plot: Relatedness vs Average Litter Size at Birth
  p_avg_litter_birth <- ggplot(ratio_data2, aes(x = Relatedness, y = avg_litter_size_birth)) +
    geom_point(color = "blue", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "darkgreen", size = 1.2) +
    labs(title = paste("Relatedness vs. Average Litter Size at Birth (Species:", species, ")"),
         x = "Relatedness",
         y = "Average Litter Size at Birth") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  # Add dynamic correlation annotation for birth plot
  p_avg_litter_birth <- add_correlation(p_avg_litter_birth, "Relatedness", "avg_litter_size_birth", ratio_data2)
  
  
  # Print and save the birth plot
  print(p_avg_litter_birth)
  
  
  ggsave(paste0("./database/Relatedness_vs_Avg_Litter_Size_Birth_", species, ".jpeg"), plot = p_avg_litter_birth, width = 8, height = 6, device = "jpeg")
  
  
  # Plot: Relatedness vs Average Litter Size at Weaning
  p_avg_litter_wean <- ggplot(ratio_data2, aes(x = Relatedness, y = avg_litter_size_wean)) +
    geom_point(color = "purple", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "orange", size = 1.2) +
    labs(title = paste("Relatedness vs. Average Litter Size at Weaning (Species:", species, ")"),
         x = "Relatedness",
         y = "Average Litter Size at Weaning") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  # Add dynamic correlation annotation for weaning plot
  p_avg_litter_wean <- add_correlation(p_avg_litter_wean, "Relatedness", "avg_litter_size_wean", ratio_data2)
  
  # Print and save the weaning plot
  print(p_avg_litter_wean)
  ggsave(paste0("./database/Relatedness_vs_Avg_Litter_Size_Weaning_", species, ".jpeg"), plot = p_avg_litter_wean, width = 8, height = 6, device = "jpeg")
  
  
  ####
  ratio_data3 <- relevant_results %>% filter(num_litters >2)
  
  # Plot: Relatedness vs Average Litters per Month
  p_litters_per_month <- ggplot(ratio_data3, aes(x = Relatedness, y = avg_litters_per_month)) +
    geom_point(color = "blue", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "darkgreen", size = 1.2) +
    labs(title = "Relatedness vs. Average Litters per Month",
         x = "Relatedness",
         y = "Average Litters per Month") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
 
  # Add dynamic correlation annotation for weaning plot
  p_litters_per_month <- add_dynamic_correlation(p_litters_per_month, "Relatedness", "avg_litters_per_month", ratio_data3)
  
  # Print and save the plot for average litters per month
  print(p_litters_per_month)
  
  
  ggsave(glue("./database/{species}- Relatedness_vs_Avg_Litters_per_Month(litters more than 2).jpeg"), plot = p_litters_per_month, width = 8, height = 6, device = "jpeg")
  
 
  
  #####################First and second litters graph
  
 
  # Step 1: Extract the first and second litters for each mating cage
  # and join with results_df to get Relatedness information
  first_two_litters <- data %>%
    group_by(`Mating Number`) %>%
    arrange(Birthday) %>%
    slice(1:2) %>%  # Keep only the first two litters
    mutate(litter_order = row_number()) %>%  # Add a column indicating 1st or 2nd litter
    ungroup() %>%
    left_join(results_df %>% dplyr:: select(MatingNumber, Relatedness), by = c(`Mating Number` = "MatingNumber")) %>%
    filter(Relatedness >0)
  
  # Separate the data for litter size at birth and at weaning for each litter order
  
  # First litter only at birth
  first_litter_birth <- first_two_litters %>%
    filter(litter_order == 1) %>%
    dplyr:: select(`Mating Number`, Relatedness, litter_size_birth = `Mice born`)
  
  # Second litter only at birth
  second_litter_birth <- first_two_litters %>%
    filter(litter_order == 2) %>%
    dplyr:: select(`Mating Number`, Relatedness, litter_size_birth = `Mice born`)
  
  # First litter only at weaning
  first_litter_wean <- first_two_litters %>%
    filter(litter_order == 1) %>%
    dplyr:: select(`Mating Number`, Relatedness, litter_size_wean = `Mice weaned`)
  
  # Second litter only at weaning
  second_litter_wean <- first_two_litters %>%
    filter(litter_order == 2) %>%
    dplyr:: select(`Mating Number`, Relatedness, litter_size_wean = `Mice weaned`)
  
  # Average of the first and second litters at birth
  average_litter_birth <- first_two_litters %>%
    group_by(`Mating Number`, Relatedness) %>%
    summarise(avg_litter_size_birth = mean(`Mice born`, na.rm = TRUE)) %>%
    ungroup()
  
  # Average of the first and second litters at weaning
  average_litter_wean <- first_two_litters %>%
    group_by(`Mating Number`, Relatedness) %>%
    summarise(avg_litter_size_wean = mean(`Mice weaned`, na.rm = TRUE)) %>%
    ungroup()
  
  # Helper function to add Pearson correlation annotation
  
  
  
  # Function to add correlation annotation at the top right of the plot
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
  
  
  # Plot 1: Relatedness vs Litter Size at First Litter (Birth)
  p_first_litter_birth <- ggplot(first_litter_birth, aes(x = Relatedness, y = litter_size_birth)) +
    geom_point(color = "blue", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "darkgreen", size = 1.2) +
    labs(title = paste("Relatedness vs Litter Size at First Litter (Birth, Species:", species, ")"),
         x = "Relatedness",
         y = "Litter Size (First Litter at Birth)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  p_first_litter_birth <- add_correlation_top_right(p_first_litter_birth, "Relatedness", "litter_size_birth", first_litter_birth)
  
  print(p_first_litter_birth)
  
  ggsave(paste0("./database/Relatedness_vs_Litter_Size_First_Litter_Birth_", species, ".jpeg"), plot = p_first_litter_birth, width = 8, height = 6, device = "jpeg")
  
  # Plot 2: Relatedness vs Litter Size at Second Litter (Birth)
  p_second_litter_birth <- ggplot(second_litter_birth, aes(x = Relatedness, y = litter_size_birth)) +
    geom_point(color = "purple", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "orange", size = 1.2) +
    labs(title = paste("Relatedness vs Litter Size at Second Litter (Birth, Species:", species, ")"),
         x = "Relatedness",
         y = "Litter Size (Second Litter at Birth)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  p_second_litter_birth <- add_correlation_top_right(p_second_litter_birth, "Relatedness", "litter_size_birth", second_litter_birth)
  print(p_second_litter_birth)
  
  ggsave(paste0("./database/Relatedness_vs_Litter_Size_Second_Litter_Birth_", species, ".jpeg"), plot = p_second_litter_birth, width = 8, height = 6, device = "jpeg")
  
  # Plot 3: Relatedness vs Average Litter Size (Birth)
  p_avg_litter_birth <- ggplot(average_litter_birth, aes(x = Relatedness, y = avg_litter_size_birth)) +
    geom_point(color = "red", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "darkblue", size = 1.2) +
    labs(title = paste("Relatedness vs Average Litter Size (Birth, Species:", species, ")"),
         x = "Relatedness",
         y = "Average Litter Size (First & Second Litters at Birth)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  p_avg_litter_birth <- add_correlation_top_right(p_avg_litter_birth, "Relatedness", "avg_litter_size_birth", average_litter_birth)
  
  print(p_avg_litter_birth)
  
  ggsave(paste0("./database/Relatedness_vs_Avg_Litter_Size(First & Second Litters at Birth)", species, ".jpeg"), plot = p_avg_litter_birth, width = 8, height = 6, device = "jpeg")
  
  # Plot 4: Relatedness vs Litter Size at First Litter (Weaning)
  p_first_litter_wean <- ggplot(first_litter_wean, aes(x = Relatedness, y = litter_size_wean)) +
    geom_point(color = "blue", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "darkgreen", size = 1.2) +
    labs(title = paste("Relatedness vs Litter Size at First Litter (Weaning, Species:", species, ")"),
         x = "Relatedness",
         y = "Litter Size (First Litter at Weaning)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  p_first_litter_wean <- add_correlation_top_right(p_first_litter_wean, "Relatedness", "litter_size_wean", first_litter_wean)
  print(p_first_litter_wean)
  
  ggsave(paste0("./database/Relatedness_vs_Litter_Size_First_Litter_Wean_", species, ".jpeg"), plot = p_first_litter_wean, width = 8, height = 6, device = "jpeg")
  
  # Plot 5: Relatedness vs Litter Size at Second Litter (Weaning)
  p_second_litter_wean <- ggplot(second_litter_wean, aes(x = Relatedness, y = litter_size_wean)) +
    geom_point(color = "purple", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "orange", size = 1.2) +
    labs(title = paste("Relatedness vs Litter Size at Second Litter (Weaning, Species:", species, ")"),
         x = "Relatedness",
         y = "Litter Size (Second Litter at Weaning)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  p_second_litter_wean <- add_correlation_top_right(p_second_litter_wean, "Relatedness", "litter_size_wean", second_litter_wean)
  print(p_second_litter_wean)
  ggsave(paste0("./database/Relatedness_vs_Litter_Size_Second_Litter_Wean_", species, ".jpeg"), plot = p_second_litter_wean, width = 8, height = 6, device = "jpeg")
  
  # Plot 6: Relatedness vs Average Litter Size (First & Second Litters at Weaning)
  p_avg_litter_wean <- ggplot(average_litter_wean, aes(x = Relatedness, y = avg_litter_size_wean)) +
    geom_point(color = "red", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "darkblue", size = 1.2) +
    labs(title = paste("Relatedness vs Average Litter Size (Weaning, Species:", species, ")"),
         x = "Relatedness",
         y = "Average Litter Size (First & Second Litters at Weaning)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  
  p_avg_litter_wean <- add_correlation_top_right(p_avg_litter_wean, "Relatedness", "avg_litter_size_wean", average_litter_wean)
  print(p_avg_litter_wean)
  
  ggsave(paste0("./database/Relatedness_vs_Avg_Litter_Size(First & Second Litters)_Wean_", species, ".jpeg"), plot = p_avg_litter_wean, width = 8, height = 6, device = "jpeg")
  
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

  
  # Step 1: Calculate summaries, including num_litters and num_losses
  # Summarize data by gender and calculate necessary metrics
  gender_loss_ratios <- data %>%
    group_by(`Mating Number`) %>%
    summarise(
      total_mice_born = sum(`Mice born`, na.rm = TRUE),    # Total mice born
      total_F_weaned = sum(`F at weaning`, na.rm = TRUE), # Total females at weaning
      total_M_weaned = sum(`M at weaning`, na.rm = TRUE), # Total males at weaning
      num_litters = n_distinct(Birthday),                # Number of distinct litters
      num_losses = sum(`Dead mice` > 0, na.rm = TRUE)    # Number of litters with losses
    ) %>%
    ungroup() %>%
    mutate(
      total_F_born = total_mice_born / 2,
      total_M_born = total_mice_born / 2,
      # Ensure loss ratios are non-negative
      F_loss_ratio = pmax(0, (total_F_born - total_F_weaned) / total_F_born),
      M_loss_ratio = pmax(0, (total_M_born - total_M_weaned) / total_M_born)
    )
  
  # Step 3: Merge with Relatedness
  gender_loss_ratios_with_relatedness <- gender_loss_ratios %>%
    left_join(results_df %>% dplyr:: select(MatingNumber, Relatedness), by = c("Mating Number" = "MatingNumber")) %>%
    filter(Relatedness > 0)  # Keep valid relatedness values
  
  # Step 4: Reshape data for plotting
  gender_loss_ratios_long <- gender_loss_ratios_with_relatedness %>%
    dplyr:: select(Relatedness, F_loss_ratio, M_loss_ratio, num_litters) %>%
    pivot_longer(
      cols = c(F_loss_ratio, M_loss_ratio),
      names_to = "Gender",
      values_to = "Loss_Ratio"
    )
  
  
  
  add_correlation_top_right <- function(plot, data, group_var, x_var, y_var) {
    # Calculate Pearson correlation and p-value for each group
    correlations <- data %>%
      group_by(!!sym(group_var)) %>%
      summarise(
        cor_test = list(cor.test(!!sym(x_var), !!sym(y_var), use = "complete.obs"))
      ) %>%
      mutate(
        label = paste0(
          ifelse(!!sym(group_var) == "F_loss_ratio", "Female: ", "Male: "),
          "R = ", round(map_dbl(cor_test, ~ .$estimate), 3),
          ", P = ", map_chr(cor_test, ~ format.pval(.$p.value, digits = 3))
        )
      )
    
    # Determine the annotation positions
    max_x <- max(data[[x_var]], na.rm = TRUE) * 0.75  # Slightly inside the max x-axis
    max_y <- max(data[[y_var]], na.rm = TRUE) * 1.1  # Above the max y-axis
    
    # Add annotations
    for (i in seq_len(nrow(correlations))) {
      y_offset <- if (correlations[[group_var]][i] == "F_loss_ratio") {
        max_y  # Place Female annotation at the very top
      } else {
        max_y * 0.9  # Place Male annotation slightly below
      }
      
      plot <- plot + annotate(
        "text",
        x = max_x,
        y = y_offset,
        label = correlations$label[i],
        hjust = 0,  # Left-align the text
        size = 3.5,
        color = ifelse(correlations[[group_var]][i] == "F_loss_ratio", "red", "blue")
      )
    }
    return(plot)
  }
  
  
  
  
  # Step 5: Plotting function
  plot_gender_loss_ratio <- function(data, filter_condition, title_suffix, species, save_path) {
    # Dynamically evaluate the filter condition
    filtered_data <- data %>%
      filter(!!rlang::parse_expr(filter_condition))
    
    # Create the plot
    plot <- ggplot(filtered_data, aes(x = Relatedness, y = Loss_Ratio, color = Gender)) +
      geom_point(size = 3) +
      geom_smooth(method = "lm", se = FALSE, size = 1.2) +
      scale_color_manual(
        values = c("F_loss_ratio" = "red", "M_loss_ratio" = "blue"),
        labels = c("F_loss_ratio" = "Female Loss Ratio", "M_loss_ratio" = "Male Loss Ratio")
      ) +
      labs(
        title = paste("Relatedness vs Loss Ratio by Gender (", title_suffix, ", Species: ", species, ")", sep = ""),
        x = "Relatedness",
        y = "Loss Ratio",
        color = "Gender"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "right"
      )
    
    # Add Pearson correlation annotations at the top-right
    plot <- add_correlation_top_right(plot, filtered_data, "Gender", "Relatedness", "Loss_Ratio")
    
    # Save the plot
    ggsave(save_path, plot = plot, width = 10, height = 6, device = "jpeg")
    
    # Return the plot for display
    return(plot)
  }
  
  
  
  
  # Step 6: Generate plots for no filter and num_litters > 2
  # Case 1: No filter on num_litters
  p_no_filter <- plot_gender_loss_ratio(
    data = gender_loss_ratios_long,
    filter_condition = "TRUE",  # No filter
    title_suffix = "No Filter on Litters",
    species = species,
    save_path = paste0("./database/Relatedness_vs_Loss_Ratio_by_gender_No_Filter_", species, ".jpeg")
  )
  print(p_no_filter)
  
  
  # Case 2: Filter num_litters > 1
  p_litters_greater_than_1 <- plot_gender_loss_ratio(
    data = gender_loss_ratios_long,
    filter_condition = "num_litters > 1",  # Filter condition as a string
    title_suffix = "Litters > 1",
    species = species,
    save_path = paste0("./database/Relatedness_vs_Loss_Ratio_by_gender_Litters_greater_than_1_", species, ".jpeg")
  )
  print(p_litters_greater_than_1)
  
  # Case 3: Filter num_litters > 2
  p_litters_greater_than_2 <- plot_gender_loss_ratio(
    data = gender_loss_ratios_long,
    filter_condition = "num_litters > 2",  # Filter condition as a string
    title_suffix = "Litters > 2",
    species = species,
    save_path = paste0("./database/Relatedness_vs_Loss_Ratio_by_gender_Litters_greater_than_2_", species, ".jpeg")
  )
  print(p_litters_greater_than_2)
  
  
  # Case 4: Filter num_litters > 3
  p_litters_greater_than_3 <- plot_gender_loss_ratio(
    data = gender_loss_ratios_long,
    filter_condition = "num_litters > 3",  # Filter condition as a string
    title_suffix = "Litters > 3",
    species = species,
    save_path = paste0("./database/Relatedness_vs_Loss_Ratio_by_gender_Litters_greater_than_3_", species, ".jpeg")
  )
  print(p_litters_greater_than_3)
  
  
  
  # Case 5: Filter num_litters > 4
  p_litters_greater_than_4 <- plot_gender_loss_ratio(
    data = gender_loss_ratios_long,
    filter_condition = "num_litters > 4",  # Filter condition as a string
    title_suffix = "Litters > 4",
    species = species,
    save_path = paste0("./database/Relatedness_vs_Loss_Ratio_by_gender_Litters_greater_than_4_", species, ".jpeg")
  )
  print(p_litters_greater_than_4)
  
  #######################
  
  
 
  # Summarize data for average litter size at weaning by gender
  gender_litter_size_weaning <- data %>%
    group_by(`Mating Number`) %>%
    summarise(
      total_F_weaned = sum(`F at weaning`, na.rm = TRUE),  # Total females at weaning
      total_M_weaned = sum(`M at weaning`, na.rm = TRUE),  # Total males at weaning
      num_litters = n_distinct(Birthday),                 # Number of distinct litters
      avg_F_litter_size = total_F_weaned / num_litters,   # Average female litter size
      avg_M_litter_size = total_M_weaned / num_litters    # Average male litter size
    ) %>%
    left_join(results_df %>% dplyr:: select(MatingNumber, Relatedness), by = c("Mating Number" = "MatingNumber")) %>%
    filter(Relatedness > 0 & num_litters > 0)  # Ensure valid relatedness and litter counts
  
  # Convert data to long format for gender-based plotting (include num_litters)
  gender_litter_size_weaning_long <- gender_litter_size_weaning %>%
    dplyr:: select(Relatedness, num_litters, avg_F_litter_size, avg_M_litter_size) %>%
    pivot_longer(
      cols = c(avg_F_litter_size, avg_M_litter_size),
      names_to = "Gender",
      values_to = "Average_Litter_Size"
    ) %>%
    mutate(Gender = recode(Gender, 
                           "avg_F_litter_size" = "Female Litter Size",
                           "avg_M_litter_size" = "Male Litter Size"))
  
  
  
  
  plot_gender_litter_size <- function(data, filter_condition = NULL, title_suffix, species, save_path) {
    # Filter the data if a condition is provided
    filtered_data <- if (!is.null(filter_condition)) {
      data %>% filter(!!rlang::parse_expr(filter_condition))
    } else {
      data
    }
    
    # Calculate correlations by gender
    correlations <- filtered_data %>%
      group_by(Gender) %>%
      summarise(
        cor_test = list(cor.test(Relatedness, Average_Litter_Size, use = "complete.obs"))
      ) %>%
      mutate(
        label = paste0(
          ifelse(Gender == "Female Litter Size", "Female: ", "Male: "),
          "R = ", round(map_dbl(cor_test, ~ .$estimate), 3),
          ", P = ", map_chr(cor_test, ~ format.pval(.$p.value, digits = 3))
        )
      )
    
    # Create the plot
    plot <- ggplot(filtered_data, aes(x = Relatedness, y = Average_Litter_Size, color = Gender)) +
      geom_point(size = 3) +
      geom_smooth(method = "lm", se = FALSE, size = 1.2) +
      scale_color_manual(values = c("Female Litter Size" = "red", "Male Litter Size" = "blue")) +
      labs(
        title = paste("Relatedness vs Average Litter Size at Weaning (", title_suffix, ", Species: ", species, ")", sep = ""),
        x = "Relatedness",
        y = "Average Litter Size",
        color = "Gender"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "right"
      )
    
    # Add Pearson correlation annotations with adjusted positions
    for (i in seq_len(nrow(correlations))) {
      y_offset <- if (correlations$Gender[i] == "Female Litter Size") {
        max(filtered_data$Average_Litter_Size, na.rm = TRUE) * 1.05  # Slightly above for Female
      } else {
        max(filtered_data$Average_Litter_Size, na.rm = TRUE) * 0.95  # Slightly below for Male
      }
      
      plot <- plot + annotate(
        "text",
        x = max(filtered_data$Relatedness, na.rm = TRUE) * 0.8,  # Adjust to 80% of max Relatedness
        y = y_offset,
        label = correlations$label[i],
        hjust = 0,  # Left-align
        size = 3.5,
        color = ifelse(correlations$Gender[i] == "Female Litter Size", "red", "blue")
      )
    }
    
    # Save the plot
    ggsave(save_path, plot = plot, width = 10, height = 6, device = "jpeg")
    
    # Return the plot for display
    return(plot)
  }
  
  
  
  
  # Case 1: No filter on num_litters
  p_no_filter <- plot_gender_litter_size(
    data = gender_litter_size_weaning_long,
    filter_condition = NULL,  # No filter
    title_suffix = "No Filter on Litters",
    species = species,
    save_path = paste0("./database/Relatedness_vs_Avg_Litter_Size_Weaning_By_Gender__No_Filter_", species, ".jpeg")
  )
  print(p_no_filter)
  
  # Case 2: num_litters > 1
  p_litters_greater_than_1 <- plot_gender_litter_size(
    data = gender_litter_size_weaning_long,
    filter_condition = "num_litters > 1",
    title_suffix = "Litters > 1",
    species = species,
    save_path = paste0("./database/Relatedness_vs_Avg_Litter_Size_Weaning_By_Gender__Litters_greater_than_1_", species, ".jpeg")
  )
  print(p_litters_greater_than_1)
  
  # Case 3: num_litters > 2
  p_litters_greater_than_2 <- plot_gender_litter_size(
    data = gender_litter_size_weaning_long,
    filter_condition = "num_litters > 2",
    title_suffix = "Litters > 2",
    species = species,
    save_path = paste0("./database/Relatedness_vs_Avg_Litter_Size_Weaning_By_Gender__Litters_greater_than_2_", species, ".jpeg")
  )
  print(p_litters_greater_than_2)
  
  # Case 4: num_litters > 3
  p_litters_greater_than_3 <- plot_gender_litter_size(
    data = gender_litter_size_weaning_long,
    filter_condition = "num_litters > 3",
    title_suffix = "Litters > 3",
    species = species,
    save_path = paste0("./database/Relatedness_vs_Avg_Litter_Size_Weaning_By_Gender__Litters_greater_than_3_", species, ".jpeg")
  )
  print(p_litters_greater_than_3)
  
  
  ################################
  ##########

  
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
  
  #################
 
  
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
  
  min_year <- min(full_data_period$Year, na.rm = TRUE)
  max_year <- 2024
  
  
  
 
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
  

######################
  # Add Year column and define the range of years
  full_data_period <- full_data %>% 
    merge(merged_data, by = "MatingNumber") %>%  
    clean_merged_columns() %>%
    mutate(Year = year(DateofMating)) %>%
    filter (Year >=1963 )
  

  # Merge the datasets
  
  merged_data <- gender_loss_ratios %>%
    rename(`MatingNumber` = `Mating Number`) %>%  # Standardize column names
    inner_join(full_data_period, by = "MatingNumber") %>% # Merge by 'MatingNumber'
    clean_merged_columns() %>% distinct()
 
  
  
  compute_continuous_correlations_by_gender <- function(data, genders, gender_labels, min_year, max_year) {
    combined_results <- data.frame()
    
    for (i in seq_along(genders)) {
      gender_var <- genders[i]
      gender_label <- gender_labels[i]
      gender_results <- data.frame(Year = integer(), Pearson_R = numeric(), P_value = numeric(), Gender = character())
      
      for (year in seq(min_year, max_year)) {
        subset_data <- data %>%
          filter(Year >= min_year & Year <= year, !is.na(Relatedness), !is.na(!!sym(gender_var)))
        
        if (nrow(subset_data) >= 5) {  # Ensure at least 5 rows for a valid correlation
          cor_test <- cor.test(subset_data$Relatedness, subset_data[[gender_var]], use = "complete.obs")
          gender_results <- rbind(
            gender_results,
            data.frame(
              Year = year,
              Pearson_R = cor_test$estimate,
              P_value = cor_test$p.value,
              Gender = gender_label
            )
          )
        } else {
          gender_results <- rbind(
            gender_results,
            data.frame(
              Year = year,
              Pearson_R = NA,
              P_value = NA,
              Gender = gender_label
            )
          )
        }
      }
      
      combined_results <- rbind(combined_results, gender_results)
    }
    
    return(combined_results)
  }
  
  
  genders <- c("F_loss_ratio", "M_loss_ratio")
  gender_labels <- c("Female", "Male")
  
  # Compute continuous correlations for each gender
  continuous_correlations_gender <- compute_continuous_correlations_by_gender(
    data = merged_data,
    genders = genders,
    gender_labels = gender_labels,
    min_year = min(merged_data$Year, na.rm = TRUE),
    max_year = max(merged_data$Year, na.rm = TRUE)
  )
  
  # Print the results for verification
  #print(continuous_correlations_gender)
  
  

  
  # Stacked Plot Function
  plot_continuous_correlation_by_gender_stacked <- function(data, y_label, title_suffix, species, save_path) {
    # Ensure all years are represented for both genders
    all_years <- seq(min(data$Year, na.rm = TRUE), max(data$Year, na.rm = TRUE))
    data <- data %>%
      complete(Year = all_years, Gender, fill = list(Pearson_R = NA, P_value = NA))
    
    # Compute -log10(p-value)
    data <- data %>%
      mutate(logP_transformed = ifelse(P_value > 0, -log10(P_value), NA))
    
    # Pearson R plot
    pearson_plot <- ggplot(data, aes(x = Year, group = Gender, color = Gender)) +
      geom_line(aes(y = Pearson_R), size = 1) +
      geom_point(aes(y = Pearson_R), size = 2) +
      labs(
        y = "Pearson R",
        color = "Gender"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),  # Hide x-axis text for the top plot
        axis.title.x = element_blank(),  # Hide x-axis title for the top plot
        legend.position = "top"
      )
    
    # P-value plot
    pvalue_plot <- ggplot(data, aes(x = Year, group = Gender, color = Gender)) +
      geom_line(aes(y = logP_transformed), linetype = "dashed", size = 1) +
      geom_point(aes(y = logP_transformed), size = 2, shape = 1) +
      geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +  # Significance threshold
      labs(
        x = "Year",
        y = "-log10(P-value)",
        color = "Gender"
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
          "(From", min(data$Year), "to Each Year)"
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
  
  
  plot_stacked <- plot_continuous_correlation_by_gender_stacked(
    data = continuous_correlations_gender,
    y_label = "Loss Ratio",
    title_suffix = "Male vs Female",
    species = "BW",
    save_path = "./database/stacked_correlation_by_gender.jpeg"
  )
  
  # Print plot
  print(plot_stacked)
  

  }
  