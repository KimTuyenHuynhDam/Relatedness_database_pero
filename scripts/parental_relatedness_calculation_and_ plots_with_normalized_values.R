# Load required libraries
library(openxlsx)
library(tidyverse)
library(broom)
library(lubridate)
library(glue)
library(ggpubr)
library(kinship2)


# Function to clean special characters in column names
clean_column_names <- function(df) {
  names(df) <- gsub("[^[:alnum:]_]", "", names(df))
  return(df)
}

# Function to clean columns after merging multiple dataframes
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

# Function to process IDs starting with "0"
clean_ids <- function(df, columns) {
  df <- df %>%
    mutate(across(all_of(columns), ~ ifelse(grepl("^0", .), paste0(substr(., 2, nchar(.)), "00000"), .)))
  return(df)
}

# Read and clean data
pero <- read.xlsx("./data/Peromyscus.xlsx", detectDates = TRUE) %>% 
  filter(! STOCK == "EPL") %>% 
  dplyr::select(1:5) %>% distinct() %>% 
  clean_column_names() %>% 
  mutate(Birthday = as.Date(Birthday, origin = "1899-12-30"), # Excel's date origin
         BirthMonth = month(Birthday),
         BirthYear = year(Birthday)) %>% filter (BirthYear <= 2024)

matingcage <- read.xlsx("./data/Mating Records.xlsx") %>%
  filter(! STOCK == "EPL") %>% 
  dplyr::select(1:5) %>% distinct() %>% 
  clean_column_names() %>%
  mutate(DateofMating = as.Date(DateofMating, origin = "1899-12-30"))

# Choose max generations:
max_generations = 5

all_stock <- c("BW")

for (species in all_stock) {
  
  # Ensure species is processed consistently in lowercase
  species_lower <- tolower(species)
  
  # Process IND data
  IND_ori <- pero %>%
    filter(str_detect(tolower(STOCK), species_lower)) %>% 
    distinct() %>% 
    mutate(across(c(MatingNumber, ID), ~str_remove_all(., regex(species_lower, ignore_case = TRUE)) %>% str_replace_all("[^[:alnum:]]", ""))) %>%
    mutate(
      Sex = str_replace_all(Sex, "[^[:alnum:]]", "") %>% str_to_upper()
    ) %>%
    mutate(
      Sex = case_when(
        Sex %in% c("F", "FEMALE", "FEM", "FEMALES") ~ "F",
        Sex %in% c("M", "MALE", "MALES") ~ "M",
        TRUE ~ Sex
      )
    ) 
  
  # Process DAMSIRE data
  DAMSIRE_ori <- matingcage %>%
    filter(str_detect(tolower(STOCK), species_lower)) %>%
    mutate(across(c(Dam, Sire, MatingNumber), ~str_remove_all(., regex(species_lower, ignore_case = TRUE)) %>% 
                    str_replace_all("[^[:alnum:]]", ""))) 
  
  common_cols <- setdiff(intersect(names(DAMSIRE_ori), names(IND_ori)), "MatingNumber")
  DAMSIRE <- DAMSIRE_ori %>% dplyr::select(-all_of(common_cols))
  
  # Remove IDs starting with "0" then change to then end "0000" in both datasets
  IND <- clean_ids(IND_ori, c("ID"))
  DAMSIRE <- clean_ids(DAMSIRE, c("Dam", "Sire"))
  
  DAMSIRE <- DAMSIRE %>%
    mutate(across(c(Dam, Sire), as.numeric)) %>%
    filter(!is.na(Dam), !is.na(Sire)) %>% distinct()
  
  IND <- IND %>%
    mutate(ID = as.numeric(ID)) %>%
    filter(!is.na(ID)) %>% distinct()
  
  # Initialize results list
  results <- list()
  problematic_mating_cages <- list()
  i <- 0
  
  # Function to get all ancestors for given IDs
  get_ancestors <- function(ids, ped_data, max_generations = 10) {
    ancestors <- ids
    for (i in 1:max_generations) {
      new_ancestors <- ped_data %>%
        filter(ID %in% ancestors) %>%
        dplyr::select(Dam, Sire) %>%
        unlist() %>%
        na.omit() %>%
        unique()
      
      if (length(new_ancestors) == 0) break
      
      ancestors <- unique(c(ancestors, new_ancestors))
    }
    return(ancestors)
  }
  
  # Create ped_data with filtering for rows where either all three columns are NA or all have values
  ped_data <- IND %>% 
    left_join(DAMSIRE, by = "MatingNumber") %>% distinct()
  #filter((is.na(MatingNumber) & is.na(Dam) & is.na(Sire)) | 
  #         (!is.na(MatingNumber) & !is.na(Dam) & !is.na(Sire))) %>%
  
  
  for (mating_number in unique(DAMSIRE$MatingNumber)) {
    i <- i + 1
    tryCatch({
      mating_pair <- DAMSIRE %>% filter(MatingNumber == mating_number)
      
      if (nrow(mating_pair) > 0) {
        StartingAnimalsOfInterest <- as_vector(mating_pair %>% dplyr::select(Dam, Sire))
        
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
  
  
  
  # Calculate average litter size
  merged_df <- merge(DAMSIRE, IND, by = 'MatingNumber') %>% 
    distinct() %>% 
    mutate(Birthday = as.Date(Birthday, format = "%Y-%m-%d"), # Adjust format as needed
           BirthMonth = month(Birthday),
           BirthYear = year(Birthday))
  
  dam_info <- merged_df %>%
    dplyr::select(ID, Birthday, BirthMonth, BirthYear) %>%
    rename(Birthday_Dam = Birthday, BirthMonth_Dam = BirthMonth, BirthYear_Dam = BirthYear)
  
  sire_info <- merged_df %>%
    dplyr::select(ID, Birthday, BirthMonth, BirthYear) %>%
    rename(Birthday_Sire = Birthday, BirthMonth_Sire = BirthMonth, BirthYear_Sire = BirthYear)
  
  merged_df2 <- merged_df %>%  
    left_join(dam_info, by = c("Dam" = "ID")) %>%
    left_join(sire_info, by = c("Sire" = "ID")) 
  
  litter_size <- merged_df2 %>%
    group_by(BirthYear, BirthMonth, Birthday, MatingNumber) %>%
    summarise(litter_size = n(), .groups = 'drop') %>%
    ungroup()
  
  litter_stats <- litter_size %>%
    group_by(MatingNumber) %>%
    summarise(
      number_of_litters = n(),
      avg_litter_size = mean(litter_size, na.rm = TRUE),
      total_pups = sum(litter_size, na.rm = TRUE),
      mating_time_days = as.numeric(difftime(max(Birthday), min(Birthday), units = "days")),
      mating_time_years = mating_time_days / 365.25,
      total_pups_per_year = total_pups / mating_time_years,
      number_of_litters_per_year = number_of_litters / mating_time_years,
      .groups = 'drop'
    ) %>%
    ungroup()
  
  # Merge average litter size to relatedness score
  mating_info <- results_df %>% 
    dplyr::select(MatingNumber, Relatedness) %>%
    right_join(DAMSIRE_ori, by = 'MatingNumber') %>%
    left_join(litter_stats, by = 'MatingNumber') %>% 
    distinct() %>% 
    filter(! is.na(number_of_litters)) %>% 
    arrange(DateofMating) # Sort by DateofMating from oldest to newest
  
  write.xlsx(mating_info, file = glue("./relatedness vs litter size/{species} - mating info with relatedness score.xlsx"))
  
  # Filter out non-sensible data
  
  mating_info <- mating_info %>%
    mutate(max_litters_per_year = ifelse(species == 'IS', 365.25 / 28, 365.25 / 23)) %>%
    filter(number_of_litters_per_year <= max_litters_per_year) %>%
    dplyr::select(-max_litters_per_year) %>%  # Remove the temporary max_litters_per_year column
    filter(!is.na(Relatedness)) %>% 
    distinct() %>% 
    filter(Relatedness > 0)
  
  
  # Function to create and save a plot with checks for valid data
  create_plot <- function(data, x_var, y_var, x_label, y_label, file_suffix) {
    # Filter out rows with NA or infinite values in the columns of interest
    valid_data <- data %>%
      filter(is.finite(.[[x_var]]) & is.finite(.[[y_var]]))
    
    # Check if there are enough valid observations for correlation and plotting
    if (nrow(valid_data) < 3) {
      warning(glue("Not enough valid observations for {x_var} vs {y_var} to perform correlation and plot."))
      return(NULL)
    }
    
    # Perform correlation test
    cor_test <- cor.test(valid_data[[x_var]], valid_data[[y_var]], method = "pearson")
    lm_model <- lm(reformulate(x_var, y_var), data = valid_data)
    
    plot <- ggplot(valid_data) +
      geom_point(aes_string(x = x_var, y = y_var)) +
      geom_abline(slope = coef(lm_model)[[x_var]], intercept = coef(lm_model)[["(Intercept)"]], col = 'red', linewidth = 1) +
      xlab(x_label) + ylab(y_label) + ggtitle(species) + theme_minimal() +
      geom_label(label = paste("Correlation:", round(cor_test$estimate, 3), "\n", "p-value:", round(cor_test$p.value, 4)), 
                 x = max(valid_data[[x_var]], na.rm = TRUE) * 5/6,
                 y = max(valid_data[[y_var]], na.rm = TRUE) * 5/6, 
                 label.padding = unit(0.55, "lines"), label.size = 0.35,
                 color = "black",
                 fill = "white") +
      theme(
        plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
      )
    
    file_name <- glue("./relatedness vs litter size/{species} - {file_suffix}.jpeg")
    ggsave(file_name, plot, width = 8, height = 5)
  }
  
  ###### Create plots of Relatedness vs other factors
  
  create_plot(mating_info, "avg_litter_size", "Relatedness", "Average litter size", "Relatedness score", "Avg litter size vs relatedness")
  create_plot(mating_info, "total_pups_per_year", "Relatedness", "Total pups per year", "Relatedness score", "Total pups per year vs relatedness")
  create_plot(mating_info, "number_of_litters_per_year", "Relatedness", "Total litters per year", "Relatedness score", "Litters per year vs relatedness")
  
  
  ##### Create plots for relatedness scores and/or average litter size per year
  
  # Convert columns to appropriate types
  mating_info <- mating_info %>%
    mutate(Year = as.factor(year(DateofMating))) %>%
    mutate(Relatedness = as.numeric(Relatedness))
  
  
  # Create the title string separately
  plot_title <- glue("{species} - Relatedness Scores per Year")
  
  # Create the filepath separately
  filepath_relatedness <- sprintf("./relatedness vs litter size/%s - relatedness scores per year.jpeg", species)
  
  # Plot relatedness scores per year with box plot and SEM
  relatedness_plot <- ggplot(mating_info, aes(x = factor(Year), y = Relatedness)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2, color = "red", fill = "red") +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
    xlab("Year") +
    ylab("Relatedness Score") +
    ggtitle(plot_title) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  # Save the plot
  ggsave(filepath_relatedness, plot = relatedness_plot, width = 8, height = 5)
  
  # Create the title string for the next plot
  plot_title_avg_litter <- glue("{species} - Average Litter Size per Year")
  
  # Create the filepath separately
  filepath_avg_litter <- sprintf("./relatedness vs litter size/%s - avg litter size per year.jpeg", species)
  
  # Plot average litter size per year with box plot and SEM
  avg_litter_plot <- ggplot(mating_info, aes(x = factor(Year), y = avg_litter_size)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "red", fill = "red") +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
    xlab("Year") +
    ylab("Average Litter Size") +
    ggtitle(plot_title_avg_litter) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  # Save the plot
  ggsave(filepath_avg_litter, plot = avg_litter_plot, width = 8, height = 5)
  
  
  # Normalize the values for combined plotting
  mating_info <- mating_info %>%
    mutate(Relatedness_norm = (Relatedness - min(Relatedness, na.rm = TRUE)) / (max(Relatedness, na.rm = TRUE) - min(Relatedness, na.rm = TRUE)),
           avg_litter_size_norm = (avg_litter_size - min(avg_litter_size, na.rm = TRUE)) / (max(avg_litter_size, na.rm = TRUE) - min(avg_litter_size, na.rm = TRUE)))
  
  # Prepare data for combined plotting
  relatedness_data <- mating_info %>%
    dplyr::select(Year, Relatedness_norm) %>%
    rename(Value = Relatedness_norm) %>%
    mutate(Type = "Relatedness")
  
  litter_size_data <- mating_info %>%
    dplyr::select(Year, avg_litter_size_norm) %>%
    rename(Value = avg_litter_size_norm) %>%
    mutate(Type = "Average Litter Size")
  
  plot_data <- bind_rows(relatedness_data, litter_size_data)
  
  # Create the title string for the combined plot
  plot_title_combined <- glue("{species} - Relatedness Scores and Average Litter Size per Year")
  
  # Create the filepath separately
  filepath_combined <- sprintf("./relatedness vs litter size/%s - relatedness and avg litter size per year.jpeg", species)
  
  # Plot combined relatedness scores and average litter size per year with box plot and SEM
  combined_plot <- ggplot(plot_data, aes(x = factor(Year), y = Value)) +
    #geom_boxplot(aes(fill = Type), alpha = 0.5, position = position_dodge(width = 0.9), width = 0.4) +
    geom_point(aes(color = Type), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), alpha = 0.5) +
    stat_summary(fun = mean, geom = "line", aes(group = Type, color = Type), size = 1, position = position_dodge(width = 0.9)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(width = 0.9), aes(color = Type)) +
    xlab("Year") +
    ylab("Normalized Value") +
    scale_y_continuous(sec.axis = sec_axis(~ ., name = "Normalized Value")) +
    ggtitle(plot_title_combined) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_color_manual(values = c("Relatedness" = "red", "Average Litter Size" = "purple")) +
    #scale_fill_manual(values = c("Relatedness" = "blue", "Average Litter Size" = "darkgreen")) +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  
  print(combined_plot)
  # Save the combined plot
  ggsave(filepath_combined, plot = combined_plot, width = 15, height = 8)
  
}  
