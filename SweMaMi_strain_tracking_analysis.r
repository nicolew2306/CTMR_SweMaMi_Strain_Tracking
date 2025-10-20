library(ANCOMBC)
library(phyloseq)
library(tidyverse)
library(dplyr)
library(DT)
library(mia)
library(caret)
library(microbiome)
library(microViz)
library(patchwork)
library(ggplot2)
library(mixOmics)
library(devtools)
library(vegan)
library(scutr)
library(BiocManager)
library(decontam)
library(compositions)
library(MASS)
library(MLeval)
library(ggalluvial)
library(ALDEx2)
library(plotly)

#
strain_sheet <- read.csv("/PATH/strain_sample_selection.csv",sep=";")
#
metaphlan_path_fecal <- list.files(path="/PATH/all_output_files", pattern="metaphlan.txt", all.files=T, full.names=T)
#
metaphlan_path_vaginal <- list.files(path="/PATH/all_output_files", pattern="metaphlan.txt", all.files=T, full.names=T)
#
metaphlan_files_rest <- list.files("/PATH/remaining_samples_host_removed",pattern="metaphlan.txt", all.files=T, full.names=T)
#
metached_strains <- read.csv("/PATH/matched_data.csv")
#
metadata <- read.csv("/PATH/swemami_cleaned_20SEPT23.csv",sep=";")
#
species_list <- read.csv("/PATH/species_list.tsv",sep="\t")
#
midas_species_db <- read.csv("/PATH/metadata.tsv",sep="\t")
#
over_3_months_samples <- read.csv("/PATH/months_over_3_all_kit_info.csv",sep=",")
#
kit_information <- read.csv("/PATH/swemami_kits_clean.csv",sep=",")
#
infant_kit_birth <- read.csv("/PATH/infant_birth_kit_ID.csv",sep=";")
#
metadata_selected <- metadata[,c("Studienummer","kit1.vaginal_sample.barcode","kit2.vaginal_sample.barcode","kit3.vaginal_sample.barcode","kit1.faecal_sample.barcode","kit2.faecal_sample.barcode","kit3.faecal_sample.barcode","kit1.salival_sample.barcode","kit2.salival_sample.barcode","kit3.salival_sample.barcode","kit3.infant_faecal_sample.barcode","fstart","deliv_start_GR_Q","Q3_X17..kejsarsnitt")]




metadata_post_delivery <- metadata[,c(1,341,349)]
colnames(metadata_post_delivery) <- c("Studienummer","c-section","breast_fed")

# for infant selection

metadata_2 <- metadata
names(metadata_2)[names(metadata_2) == "Q3_X25..Har.du.kunnat.amma.ditt.dina.barn.efter.färlossningen"] <- "breastfeeding"
names(metadata_2)[names(metadata_2) == "Q3_X17..kejsarsnitt"] <- "c_section"

 metadata_subset <- metadata_post_delivery
colnames(metadata_subset) <- c("Studienummer","c_section","breastfeeding")


# vaginal selection

metadata_subset_for_vaginal <- metadata[,c("Studienummer","Primipara")]


# Step 0: cleaning up the dates
infant_birth_age <- merge(infant_kit_birth[,c(1,2)],over_3_months_samples[,c(3,16)],by="Studienummer")
mother_kit_data <- merge(metadata[,c(1,516)],kit_information[,c(11,28,36)],by="kit1.vaginal_sample.barcode")
mother_infant_submission_data <- merge(mother_kit_data[,c(2:4)],infant_birth_age,by="Studienummer")


# Step 1: Remove the time portion (if needed) and convert to Date
mother_infant_submission_data$cleaned_faecal_date <- as.Date(
  sub("T.*", "", mother_infant_submission_data$kit3.faecal_sample.datetime_taken)
)
mother_infant_submission_data$cleaned_vaginal_date <- as.Date(
  sub("T.*", "", mother_infant_submission_data$kit3.vaginal_sample.datetime_taken)
)
mother_infant_submission_data$birth_date <- as.Date(mother_infant_submission_data$fodelsedatum_barn_register)

# Step 2: Calculate the differences in weeks
mother_infant_submission_data$difference_in_weeks_faecal <- as.numeric(
  difftime(mother_infant_submission_data$cleaned_faecal_date, 
           mother_infant_submission_data$birth_date, 
           units = "weeks")
)
mother_infant_submission_data$difference_in_weeks_vaginal <- as.numeric(
  difftime(mother_infant_submission_data$cleaned_vaginal_date, 
           mother_infant_submission_data$birth_date, 
           units = "weeks")
)

mother_infant_submission_date <- mother_infant_submission_data[,c(1,5,9,10)]
mother_infant_submission_date


midas_species_db_reduced <- midas_species_db[,c(2,5,6)]
species_list_names <- merge(species_list,midas_species_db_reduced,by="species_id")


species_names <- c("Bifidobacterium_longum","Lacticaseibacillus_rhamnosus","Lactobacillus_crispatus",
                  "Bifidobacterium_breve","Ruminococcus_gnavus","Bifidobacterium_adolescentis","Blautia_massiliensis",
                  "Faecalibacterium_prausnitzii","Lactobacillus_gasseri","Lactobacillus_iners","Lactobacillus_jensenii",
                  "Bifidobacterium_vaginale","Fusobacterium_nucleatum")

species_list_names$species <- species_names



FEAST_outputs <- list.files(path="/PATH/all_FEAST_outputs_combined", pattern=".txt", all.files=T, full.names=T)

for (files in FEAST_outputs) {
    
    file_name <-  str_split(files,"/")
    file_name <- data.frame(file_name)
    file_name_2 <- str_split(file_name[12,],".txt")
    file_name_2 <- data.frame(file_name_2)
    name <- paste(file_name_2[1,],"",sep="")
    print(name)
    csv_file <- read.csv(files,sep="\t")
    names(csv_file) <- sub("_.*", "", names(csv_file))
    names(csv_file) <- sub("^X", "", names(csv_file))
    csv_file[1,] <- sub("_.*", "", csv_file[1,])
    csv_file[1,] <- sub("^X", "", csv_file[1,])
    assign(name,csv_file)
}



#Doing this for all the V2s
# List all dataframes in the environment
df_list <- ls(pattern = "^V2_V1_FEAST_output_\\d+_.*_source_contributions_matrix")

# Function to add new column with extracted date part
add_species_column <- function(df_name) {
  # Extract the date part from the dataframe name
  species_part <- str_extract(df_name, "V2_V1_FEAST_output_(\\d+)")
  species_part <- str_remove(species_part, "V2_V1_FEAST_output_")
  
  # Get the dataframe
  df <- get(df_name)
  
  # Add new column with the date part
  df$species <- species_part
  colnames(df) <- c("V2_ID","V1_strain_content","Unknown","species")
  # Assign the updated dataframe back to the original name
  assign(df_name, df, envir = .GlobalEnv)
}

# Apply the function to each dataframe
updated_dfs <- lapply(df_list, add_species_column)
combined_V2_V1_strains_FEAST <- do.call(rbind, updated_dfs)

### merge with species names

V2_V1_strains_FEAST_species_names <- merge(combined_V2_V1_strains_FEAST,species_list_names[,c(1,4)],by.x="species",by.y="species_id")



## let's do for V3_V2_V1 

# List all dataframes in the environment
df_list <- ls(pattern = "V3_V2_V1_FEAST_output_\\d+_.*_source_contributions_matrix")

add_date_column <- function(df_name) {
  # Extract the date part from the dataframe name
  date_part <- str_extract(df_name, "V3_V2_V1_FEAST_output_(\\d+)")
  date_part <- str_remove(date_part, "V3_V2_V1_FEAST_output_")
  
  # Get the dataframe
  df <- get(df_name)
  
  # Add new column with the date part
  df$species <- date_part
  
  return(df)
}

# Apply the function to each dataframe and store in a list
updated_dfs <- lapply(df_list, add_date_column)

# Separate dataframes with fewer than 5 columns
dfs_with_fewer_columns <- list()
dfs_with_five_columns <- list()

for (i in 1:length(updated_dfs)) {
  df <- updated_dfs[[i]]
  original_df <- get(df_list[i])
  if (ncol(original_df) < 4) {
    dfs_with_fewer_columns[[df_list[i]]] <- df
  } else {
    dfs_with_five_columns[[length(dfs_with_five_columns) + 1]] <- df
  }
}

# Define the standard column names for the dataframes
standard_colnames <- c("V3_ID", "strain_V2", "strain_V1", "Unknown", "species")

# Function to standardize column names
standardize_columns <- function(df) {
  colnames(df) <- standard_colnames[1:ncol(df)]
  return(df)
}

# Standardize column names for the dataframes with 5 columns
dfs_with_five_columns <- lapply(dfs_with_five_columns, standardize_columns)

# Combine all dataframes with 5 columns into one
V3_V2_V1_strains_FEAST <- do.call(rbind, dfs_with_five_columns)

### ferew colnames
fewer_column_names <- c("V3_ID","strain","Unknown","species")

standardize_fewer_columns <- function(df) {
  colnames(df) <- fewer_column_names[1:ncol(df)]
  return(df)
}

# dfs with less than five columns
dfs_with_fewer_columns <- lapply(dfs_with_fewer_columns,standardize_fewer_columns)



new_df <- dfs_with_fewer_columns$V3_V2_V1_FEAST_output_100505_S17642f1_source_contributions_matrix 

new_df$strain_V2 <- 0

# Reorder columns to place the new column after the first column
new_df <- new_df %>%
  relocate(strain_V2, .after = V3_ID)
colnames(new_df) <- c("V3_ID", "strain_V2", "strain_V1", "Unknown", "species")


## for the second one

new_df_2 <- dfs_with_fewer_columns$V3_V2_V1_FEAST_output_100505_Sf82011c_source_contributions_matrix

new_df_2$strain_V1 <- 0

# Reorder columns to place the new column after the first column
new_df_2 <- new_df_2 %>%
  relocate(strain_V1, .after = strain)
colnames(new_df_2) <- c("V3_ID", "strain_V2", "strain_V1", "Unknown", "species")


fewer_than_five_combined <- rbind(new_df,new_df_2)

### now combining everything

V3_V2_V1_strains_FEAST <- rbind(V3_V2_V1_strains_FEAST,fewer_than_five_combined)

V3_V2_V1_strains_FEAST_species_names <- merge(V3_V2_V1_strains_FEAST,species_list_names[,c(1,4)],by.x="species",by.y="species_id")


#c Doing the combining for all the F2 F1s
df_list <- ls(pattern = "^F2_F1_FEAST_output_\\d+_.*_source_contributions_matrix")

# Function to add new column with extracted date part
add_species_column <- function(df_name) {
  # Extract the date part from the dataframe name
  species_part <- str_extract(df_name, "F2_F1_FEAST_output_(\\d+)")
  species_part <- str_remove(species_part, "F2_F1_FEAST_output_")
  
  # Get the dataframe
  df <- get(df_name)
  
  # Add new column with the date part
  df$species <- species_part
  colnames(df) <- c("F2_ID","F1_strain_content","Unknown","species")
  # Assign the updated dataframe back to the original name
  assign(df_name, df, envir = .GlobalEnv)
}

# Apply the function to each dataframe
updated_dfs <- lapply(df_list, add_species_column)
combined_F2_F1_strains_FEAST <- do.call(rbind, updated_dfs)

F2_F1_strains_FEAST_species_names <- merge(combined_F2_F1_strains_FEAST,species_list_names[,c(1,4)],by.x="species",by.y="species_id")


## let's do for F3_F2_F1

# List all dataframes in the environment
df_list <- ls(pattern = "F3_F2_F1_FEAST_output_\\d+_.*_source_contributions_matrix")

add_date_column <- function(df_name) {
  # Extract the date part from the dataframe name
  date_part <- str_extract(df_name, "F3_F2_F1_FEAST_output_(\\d+)")
  date_part <- str_remove(date_part, "F3_F2_F1_FEAST_output_")
  
  # Get the dataframe
  df <- get(df_name)
  
  # Add new column with the date part
  df$species <- date_part
  
  return(df)
}

# Apply the function to each dataframe and store in a list
updated_dfs <- lapply(df_list, add_date_column)

# Separate dataframes with fewer than 5 columns
dfs_with_fewer_columns <- list()
dfs_with_five_columns <- list()

for (i in 1:length(updated_dfs)) {
  df <- updated_dfs[[i]]
  original_df <- get(df_list[i])
  if (ncol(original_df) < 4) {
    dfs_with_fewer_columns[[df_list[i]]] <- df
  } else {
    dfs_with_five_columns[[length(dfs_with_five_columns) + 1]] <- df
  }
}

# Define the standard column names for the dataframes
standard_colnames <- c("F3_ID", "strain_F2", "strain_F1", "Unknown", "species")

# Function to standardize column names
standardize_columns <- function(df) {
  colnames(df) <- standard_colnames[1:ncol(df)]
  return(df)
}

# Standardize column names for the dataframes with 5 columns
dfs_with_five_columns <- lapply(dfs_with_five_columns, standardize_columns)

# Combine all dataframes with 5 columns into one
F3_F2_F1_strains_FEAST <- do.call(rbind, dfs_with_five_columns)

### ferew colnames
fewer_column_names <- c("F3_ID","strain","Unknown","species")

standardize_fewer_columns <- function(df) {
  colnames(df) <- fewer_column_names[1:ncol(df)]
  return(df)
}

# dfs with less than five columns
dfs_with_fewer_columns <- lapply(dfs_with_fewer_columns,standardize_fewer_columns)

F3_F2_F1_strains_FEAST_species_names <- merge(F3_F2_F1_strains_FEAST,species_list_names[,c(1,4)],by.x="species",by.y="species_id")



# List all dataframes in the environment
df_list <- ls(pattern = "baby_FEAST_output_\\d+_.*")

add_date_column <- function(df_name) {
  # Extract the date part from the dataframe name
  date_part <- str_extract(df_name, "baby_FEAST_output_(\\d+)")
  date_part <- str_remove(date_part, "baby_FEAST_output_")
  
  # Get the dataframe
  df <- get(df_name)
  
  # Add new column with the date part
  df$species <- date_part
  
  return(df)
}

# Apply the function to each dataframe and store in a list
updated_dfs <- lapply(df_list, add_date_column)

# Separate dataframes with fewer than 5 columns
dfs_with_fewer_columns <- list()
dfs_with_five_columns <- list()

for (i in 1:length(updated_dfs)) {
  df <- updated_dfs[[i]]
  original_df <- get(df_list[i])
  if (ncol(original_df) < 4) {
    dfs_with_fewer_columns[[df_list[i]]] <- df
  } else {
    dfs_with_five_columns[[length(dfs_with_five_columns) + 1]] <- df
  }
}

# Define the standard column names for the dataframes
standard_colnames <- c("infant_ID", "strain_V3", "strain_F3", "Unknown", "species")

# Function to standardize column names
standardize_columns <- function(df) {
  colnames(df) <- standard_colnames[1:ncol(df)]
  return(df)
}

# Standardize column names for the dataframes with 5 columns
dfs_with_five_columns <- lapply(dfs_with_five_columns, standardize_columns)

# Combine all dataframes with 5 columns into one
infant_strains_FEAST <- do.call(rbind, dfs_with_five_columns)

### ferew colnames
fewer_column_names <- c("infant_ID","strain_F3","Unknown","species")

standardize_fewer_columns <- function(df) {
  colnames(df) <- fewer_column_names[1:ncol(df)]
  return(df)
}


dfs_with_fewer_columns <- lapply(dfs_with_fewer_columns,standardize_fewer_columns)
infant_strains_FEAST_only_fecal <- do.call(rbind, dfs_with_fewer_columns)

rownames(infant_strains_FEAST_only_fecal) <- NULL

## Not many here. Gotta get them out of metaphlan. woohoo. So, now combining the five and the four (fixed) ones

rownames(infant_strains_FEAST_only_fecal) <- NULL
infant_strains_FEAST_only_fecal$strain_V3 <- rep(0,31)

infant_strains_FEAST_only_fecal <- infant_strains_FEAST_only_fecal %>% 
  relocate(strain_V3, .after = infant_ID)

df <- dfs_with_five_columns[[1]]

infant_strains_FEAST<- rbind(df,infant_strains_FEAST_only_fecal)

infant_strains_FEAST_species_names <- merge(infant_strains_FEAST,species_list_names[,c(1,4)],by.x="species",by.y="species_id")



metaphlan_input_files <- c(metaphlan_path_fecal, metaphlan_path_vaginal, metaphlan_files_rest)

## combining everything and then removing all the first instances of duplicates
# and then get rid of those in this file and the combine the resequenced removed with the ones that we
# currently have in here.

## files are read in order of flowcell ID

n <- 0


### So, the goal here is to clean up the files and later separate the T1 and T2 fecal analysis files.
for (species_file in metaphlan_input_files){
  
   species_input <- read.csv(species_file,sep = "\t")
   if (n == 0) {combined_species_files <- species_input} else {combined_species_files <- merge(combined_species_files,species_input,by="clade_name",all=T)}
    n <- n + 1  
    print(species_file)
  
}

### getting rid of NAs which occur when certain samples do not have a species in them. We simply set them to zero.
combined_species_files[is.na(combined_species_files)] <- 0

### for breaking down the long names for the samples to only the sample ID
   
sample_full_names <- colnames(combined_species_files)
sample_full_names <- data.frame(sample_full_names)

### this leads to a dataframe of the names split into 3 rows. The second rows is the sample ID and the one we need.
sample_IDs <- lapply(sample_full_names, function(x) str_split(x,"__"))
sample_IDs <- data.frame(sample_IDs)
                                 
sample_IDs_2 <- t(sample_IDs)
sample_IDs_2 <- data.frame(sample_IDs_2)
                    

### Now get rid of the subscripts of the name
                     
sample_IDs_2$X2 <- sapply(sample_IDs_2$X2,function(x) gsub("_1","",x))
sample_IDs_2$X2 <- sapply(sample_IDs_2$X2,function(x) gsub("_2","",x))
sample_IDs_2$X2 <- sapply(sample_IDs_2$X2,function(x) gsub("_3","",x))

                     
### Now here I give the new cleaned column names to my columns
combined_species_with_sample_IDs <- combined_species_files
colnames(combined_species_with_sample_IDs) <- sample_IDs_2$X2
                          
                          
combined_species_with_sample_IDs_2 <- combined_species_with_sample_IDs                           
                          
rownames(combined_species_with_sample_IDs_2) <- combined_species_with_sample_IDs_2$clade_name
combined_species_with_sample_IDs_2 <- combined_species_with_sample_IDs_2[,-c(1)] 
                          
### since the files were read in flowcell order, the final repeat is the most recently sequenced one
## so going to keep that, so will remove the first occurrences of the repeated files.
                          
      
                          
reversed_df <- combined_species_with_sample_IDs_2[, rev(colnames(combined_species_with_sample_IDs_2))]
colnames(reversed_df) <- sub("\\.\\d+$", "", colnames(reversed_df))
reversed_df_filtered <- reversed_df[, !duplicated(colnames(reversed_df))]
combined_samples_unique <- reversed_df_filtered
 
names(combined_samples_unique) <- sub("^X", "", names(combined_samples_unique))                          
                          

# species only

metaphlan_input <- combined_samples_unique
species_metaphlan_input <- NULL

for (i in 1:nrow(metaphlan_input)) {

    if (grepl("s__", rownames(metaphlan_input)[i])== TRUE && grepl("t__",rownames(metaphlan_input)[i]) == FALSE) {

        species_metaphlan_input <- rbind(species_metaphlan_input,metaphlan_input[i,])
    } 

}

# strains only


metaphlan_input <- combined_samples_unique
species_metaphlan_input_strains <- NULL

for (i in 1:nrow(metaphlan_input)) {

    if (grepl("t__",rownames(metaphlan_input)[i]) == TRUE) {

        species_metaphlan_input_strains <- rbind(species_metaphlan_input_strains,metaphlan_input[i,])
    } 

}

# species and their corresponding strains

metaphlan_input <- combined_samples_unique
species_and_strain_metaphlan_input <- NULL

for (i in 1:nrow(metaphlan_input)) {

    if (grepl("s__", rownames(metaphlan_input)[i])== TRUE || grepl("t__",rownames(metaphlan_input)[i]) == TRUE) {

        species_and_strain_metaphlan_input <- rbind(species_and_strain_metaphlan_input,metaphlan_input[i,])
    } 

}


rownames2 <- rownames(species_and_strain_metaphlan_input)
# Extract species names (before the pipe) for all rows
species_in_strain <- sub("\\|.*", "", rownames2)

# Identify rownames without a pipe (species without strains)
no_pipe <- !grepl("\\|", rownames2)

# Check species that appear only without a pipe
species_no_strains <- rownames2[no_pipe][!rownames2[no_pipe] %in% species_in_strain[!no_pipe]]

## removing species with no entries. but I don't think there are any so let's keep them.
species_metaphlan_input <- species_metaphlan_input[rowSums(species_metaphlan_input) != 0,]


# species only
species_names <- lapply(rownames(species_metaphlan_input), function(x) str_split(x,"s__"))
species_names <- data.frame(species_names)

species_names_t <- t(species_names)
species_names_t <- data.frame(species_names_t)
                    
rownames(species_metaphlan_input) <- species_names_t$X2
 
# strains only
                        
strain_names_metaphlan <- lapply(rownames(species_metaphlan_input_strains), function(x) str_split(x,"s__"))
strain_names <- data.frame(strain_names_metaphlan)

strain_names_t <- t(strain_names)
strain_names_t <- data.frame(strain_names_t)
                    
rownames(species_metaphlan_input_strains) <- strain_names_t$X2                       

# species and strains only
                                 
strain_and_species_names_metaphlan <- lapply(rownames(species_and_strain_metaphlan_input), function(x) str_split(x,"s__"))
strain_and_species_names <- data.frame(strain_and_species_names_metaphlan)

strain_and_species_names_t <- t(strain_and_species_names)
strain_and_species_names_t <- data.frame(strain_and_species_names_t)
                    
rownames(species_and_strain_metaphlan_input) <- strain_and_species_names_t$X2
                                      

species_and_strain_metaphlan_t <- t(species_and_strain_metaphlan_input)
species_and_strain_metaphlan_t <- data.frame(species_and_strain_metaphlan_t)
species_strains_F3 <- merge(metached_strains[,c(2,10,16)],species_and_strain_metaphlan_t,by.y=0,
                             by.x="kit3.faecal_sample.barcode")
species_strains_V3 <- merge(metached_strains[,c(2,7,16)],species_and_strain_metaphlan_t,by.y=0,
                             by.x="kit3.vaginal_sample.barcode")
species_strains_inf <- merge(metached_strains[,c(2,4,16)],species_and_strain_metaphlan_t,by.y=0,
                             by.x="kit3.infant_faecal_sample.barcode")

## changing name of first column
colnames(species_strains_F3)[1] <- "kit_number"
colnames(species_strains_V3)[1] <- "kit_number"
colnames(species_strains_inf)[1] <- "kit_number"

## add specifier column
species_strains_F3$specifier <- paste0(species_strains_F3$Studienummer,"_F3")
species_strains_V3$specifier <- paste0(species_strains_V3$Studienummer,"_V3")
species_strains_inf$specifier <- paste0(species_strains_inf$Studienummer,"_zinf")

## merging the time point 3 variables
species_strains_T3 <- rbind(species_strains_F3,species_strains_V3,species_strains_inf)

# Function to extract species-specific strains
extract_species_strains <- function(df, species_name, base_cols = c(1,2,4366)) {
  species_cols <- df[, grep(species_name, colnames(df))]
  result <- cbind(df[, base_cols], species_cols)
  return(result)
}

# List of species 
species_list <- c(
  "Bifidobacterium_longum",
  "Bifidobacterium_breve",
  "Bifidobacterium_adolescentis",
  "Bifidobacterium_vaginale",
  "Lactobacillus_crispatus",
  "Lactobacillus_gasseri",
  "Lactobacillus_iners",
  "Lactobacillus_jensenii",
  "Ruminococcus_gnavus",
  "Blautia_massiliensis",
  "Faecalibacterium_prausnitzii",
  "Lacticaseibacillus_rhamnosus",
  "Fusobacterium_nucleatum",
  "Gardnerella_vaginalis"
)


## picking only the species we're interested in
species_strains_Bifi <- species_strains_T3[,grep("Bifidobacterium",colnames(species_strains_T3))]
species_strains_Fuso <- species_strains_T3[,grep("Fusobacterium",colnames(species_strains_T3))]
species_strains_Lacto <- species_strains_T3[,grep("Lactobacillus",colnames(species_strains_T3))]
species_strains_Blautia <- species_strains_T3[,grep("Blautia",colnames(species_strains_T3))]
species_strains_Feacal <- species_strains_T3[,grep("Faecalibacterium",colnames(species_strains_T3))]
species_strains_Lacticase <- species_strains_T3[,grep("Lacticaseibacillus",colnames(species_strains_T3))]



##combining them
species_strain_bifido_T3 <- cbind(species_strains_T3[,c(1:3,4366)],species_strains_Bifi)
species_strain_fuso_T3 <- cbind(species_strains_T3[,c(1:3,4366)],species_strains_Fuso)
species_strain_lacto_T3 <- cbind(species_strains_T3[,c(1:3,4366)],species_strains_Lacto)
species_strain_blautia_T3 <- cbind(species_strains_T3[,c(1:3,4366)],species_strains_Blautia)
species_strain_faecalo_T3 <- cbind(species_strains_T3[,c(1:3,4366)],species_strains_Feacal)
species_strain_lacticase_T3 <- cbind(species_strains_T3[,c(1:3,4366)],species_strains_Lacticase)
species_strains_of_interest_T3 <- cbind(species_strains_T3[,c(1:3,4366)],species_strains_Bifi,species_strains_Fuso,
                                       species_strains_Lacto,species_strains_Blautia,species_strains_Feacal,species_strains_Lacticase)


# Creating a named list
species_strains_list <- lapply(species_list, function(sp) {
  extract_species_strains(species_strains_T3, sp)
})
names(species_strains_list) <- paste0("species_strains_", 
                                      sapply(strsplit(species_list, "_"), function(x) paste0(substr(x[1],1,1), "_", x[2])), 
                                      "_T3")

# Turn each element of the list into a separate object
list2env(species_strains_list, envir = .GlobalEnv)

species_and_strain_metaphlan_t <- t(species_and_strain_metaphlan_input)
species_and_strain_metaphlan_t <- data.frame(species_and_strain_metaphlan_t)
species_strains_V2 <- merge(metached_strains[,c(2,6,16)],species_and_strain_metaphlan_t,by.y=0,
                             by.x="kit2.vaginal_sample.barcode")
species_strains_V3 <- merge(metached_strains[,c(2,7,16)],species_and_strain_metaphlan_t,by.y=0,
                             by.x="kit3.vaginal_sample.barcode")
species_strains_V1 <- merge(metached_strains[,c(2,5,16)],species_and_strain_metaphlan_t,by.y=0,
                             by.x="kit1.vaginal_sample.barcode")

## changing name of first column
colnames(species_strains_V2)[1] <- "kit_number"
colnames(species_strains_V3)[1] <- "kit_number"
colnames(species_strains_V1)[1] <- "kit_number"

## add specifier column
species_strains_V2$specifier <- paste0(species_strains_V2$Studienummer,"_V2")
species_strains_V3$specifier <- paste0(species_strains_V3$Studienummer,"_V3")
species_strains_V1$specifier <- paste0(species_strains_V1$Studienummer,"_V1")

## merging the time point 3 variables
species_strains_V3_V2_V1 <- rbind(species_strains_V2,species_strains_V3,species_strains_V1)
species_strains_V2_V1 <- rbind(species_strains_V2,species_strains_V1)


# V3 V2 V1
species_strains_list <- lapply(species_list, function(sp) {
  extract_species_strains(species_strains_V3_V2_V1, sp)
})

names(species_strains_list) <- paste0("species_strains_", species_list, "_V3_V2_V1")
list2env(species_strains_list, envir = .GlobalEnv)

# V2 and V1
species_strains_list <- lapply(species_list, function(sp) {
  extract_species_strains(species_strains_V2_V1, sp)
})

names(species_strains_list) <- paste0("species_strains_", species_list, "_V2_V1")
list2env(species_strains_list, envir = .GlobalEnv)

species_and_strain_metaphlan_t <- t(species_and_strain_metaphlan_input)
species_and_strain_metaphlan_t <- data.frame(species_and_strain_metaphlan_t)
species_strains_F2 <- merge(metached_strains[,c(2,9,16)],species_and_strain_metaphlan_t,by.y=0,
                             by.x="kit2.faecal_sample.barcode")
species_strains_F3 <- merge(metached_strains[,c(2,10,16)],species_and_strain_metaphlan_t,by.y=0,
                             by.x="kit3.faecal_sample.barcode")
species_strains_F1 <- merge(metached_strains[,c(2,8,16)],species_and_strain_metaphlan_t,by.y=0,
                             by.x="kit1.faecal_sample.barcode")

## changing name of first column
colnames(species_strains_F2)[1] <- "kit_number"
colnames(species_strains_F3)[1] <- "kit_number"
colnames(species_strains_F1)[1] <- "kit_number"

## add specifier column
species_strains_F2$specifier <- paste0(species_strains_F2$Studienummer,"_F2")
species_strains_F3$specifier <- paste0(species_strains_F3$Studienummer,"_F3")
species_strains_F1$specifier <- paste0(species_strains_F1$Studienummer,"_F1")

## merging the time point 3 variables
species_strains_F3_F2_F1 <- rbind(species_strains_F2,species_strains_F3,species_strains_F1)
species_strains_F2_F1 <- rbind(species_strains_F2,species_strains_F1)


# F3 F2 F1
species_strains_list <- lapply(species_list, function(sp) {
  extract_species_strains(species_strains_F3_F2_F1, sp)
})

names(species_strains_list) <- paste0("species_strains_", species_list, "_F3_F2_F1")
list2env(species_strains_list, envir = .GlobalEnv)


# F2 and F1

species_strains_list <- lapply(species_list, function(sp) {
  extract_species_strains(species_strains_F2_F1, sp)
})

names(species_strains_list) <- paste0("species_strains_", species_list, "_F2_F1")
list2env(species_strains_list, envir = .GlobalEnv)

# Function defined to make the different dataframes the same
reorder_rename <- function(df_0) {
if ("package:mia" %in% search()) {
  detach("package:mia", unload=TRUE)
}
if ("package:ANCOMBC" %in% search()) {
  detach("package:ANCOMBC", unload=TRUE)
}    

    df <- df_0[,c(1:4)]
    colnames(df)[4] <- "abundance"
    name <- colnames(df_0)[4]
    df$species <-  c(rep(name,nrow(df)))
    df <- df %>%
  mutate(sample_type = case_when(
    grepl("F3", specifier) ~ "fecal",
    grepl("V3", specifier) ~ "vaginal",
    grepl("zinf", specifier) ~ "infant"
  ))
 df <- as.data.frame(df)

# Create a dataframe of all possible combinations of Studienummer and sample_type
    all_combinations <- df %>%
      distinct(Studienummer) %>% # Get all unique Studienummer
      expand_grid(sample_type = c("fecal", "vaginal", "infant")) # Create all combinations with sample types
  all_combinations <- as.data.frame(all_combinations)

# Merge the original dataframe with the complete combinations, filling missing values with 0 for abundance
    df_complete <- all_combinations %>%
      left_join(df, by = c("Studienummer", "sample_type")) %>% # Left join to keep all combinations
       mutate(
        abundance = ifelse(is.na(abundance), 0, abundance), # Fill missing abundance with 0
        kit_number = ifelse(is.na(kit_number), "missing", kit_number), # Fill missing kit_number with placeholder
        species = ifelse(is.na(species), name, species),# Fill missing species with placeholder
       )
    df_complete$kit_number <- NULL
    df_complete$specifier <- NULL

    return(df_complete)
}
 

# Combining the species of interest into one sample

species_strains_metaphlan_combined <- rbind(reorder_rename(species_strains_B_longum_T3_ordered),
                                            reorder_rename(species_strains_B_adolesc_T3_ordered),
                                            reorder_rename(species_strains_B_breve_T3_ordered),
                                            reorder_rename(species_strains_B_massilien_T3_ordered),
                                            reorder_rename(species_strains_F_nucl_T3_ordered), 
                                            reorder_rename(species_strains_F_prau_T3_ordered),
                                            reorder_rename(species_strains_L_crispatus_T3_ordered),
                                            reorder_rename(species_strains_L_gasseri_T3_ordered),
                                            reorder_rename(species_strains_L_iners_T3_ordered), 
                                            reorder_rename(species_strains_L_jensenii_T3_ordered),
                                            reorder_rename(species_strains_G_vaginalis_T3_ordered),
                                            reorder_rename(species_strains_L_rhamnosus_T3_ordered),
                                            reorder_rename(species_strains_R_gnavus_T3_ordered))

species_strains_metaphlan_combined <- species_strains_metaphlan_combined %>%
  mutate(sample_type = factor(sample_type, levels = c("vaginal", "fecal", "infant")))
  
species_strains_metaphlan_meta <- merge(species_strains_metaphlan_combined,metadata_subset,by="Studienummer")

# clearning spike in
F_names_2 <- c("Alicyclobacillus_acidiphilus","Alicyclobacillus_acidiphilus.t__SGB30417")

grep("Alicyclobacillus_acidiphilus",rownames(species_metaphlan_input))
species_strains_T3_2_2 <- species_strains_T3[,-c(960,961)]
species_metaphlan_input_2_2 <- species_metaphlan_input[-c(467),]

group.colors <- c("Lactobacillus_crispatus" = "#148F77", 
    "Lactobacillus_iners" = "#A3E4D7", 
    "Lactobacillus_gasseri" ="#6D9F06", 
    "Lactobacillus_jensenii" = "#4E7705", 
    "Gardnerella_vaginalis"= "#CC79A7", 
    "Fannyhessea_vaginae"="#EFB6D6",
    "Fusobacterium_nucleatum"="#F09163",
    "Bifidobacterium_longum"="#56B4E9",
    "Bifidobacterium_breve"="#4991BA",
    "Bifidobacterium_adolescentis"="#154a69",            
    "Faecalibacterium_prausnitzii"="#854E01",
    "Lacticaseibacillus_rhamnosus"="#0a4e57",
    "Ruminococcus_gnavus"="#DD65E6",
    "Blautia_massiliensis"="#e6d765",
    "OtherSpeciesSum"="gray"
    )



species_strains_T3_2 <- species_strains_T3_2_2[, c(1:2, ncol(species_strains_T3_2_2), 3:(ncol(species_strains_T3_2_2) - 1))]
species_strains_T3_species <- species_strains_T3_2[,rownames(species_metaphlan_input_2_2)]
species_strains_T3_species_2 <- cbind(species_strains_T3_2[,(1:4)],species_strains_T3_species)
colnames(species_strains_T3_species_2)[4] <- "c_section"

species_strains_T3_species_3 <- species_strains_T3_species_2 %>%
  mutate(sample_type = case_when(
    grepl("F3", specifier) ~ "fecal",
    grepl("V3", specifier) ~ "avaginal",
    grepl("zinf", specifier) ~ "zinfant"
  ),
    c_section = case_when(
    c_section == 0 ~ "vaginal_delivery",
    c_section == 1 ~ "c_section",
    TRUE ~ as.character(c_section)  # Keeps original values for any other cases
  ))
species_strains_T3_species_3 <- species_strains_T3_species_3[, c(1:2, ncol(species_strains_T3_species_3), 3:(ncol(species_strains_T3_species_3) - 1))]
species_strains_T3_species_3 <- species_strains_T3_species_3[,!names(species_strains_T3_species_3) %in% "Alicyclobacillus_acidiphilus"]
species_strains_T3_species_3 <- species_strains_T3_species_3[order(species_strains_T3_species_3$Studienummer),]
species_strains_T3_species_4 <- species_strains_T3_species_3[,-c(1,4)]

species_names <- c("Bifidobacterium_longum","Lacticaseibacillus_rhamnosus","Lactobacillus_crispatus",
                  "Bifidobacterium_breve","Ruminococcus_gnavus","Bifidobacterium_adolescentis","Blautia_massiliensis",
                  "Faecalibacterium_prausnitzii","Lactobacillus_gasseri","Lactobacillus_iners","Lactobacillus_jensenii",
                  "Fusobacterium_nucleatum","Gardnerella_vaginalis")


# Create a new column for the sum of all other species
species_strains_T3_species_5 <- species_strains_T3_species_4 %>%
  mutate(
    OtherSpeciesSum = rowSums(dplyr::select(., -c(Studienummer, sample_type, c_section, all_of(species_names))))
  ) %>%
  dplyr::select(Studienummer, sample_type, c_section, all_of(species_names), OtherSpeciesSum)


species_strains_T3_species_6 <- merge(mother_infant_submission_date[,c(1,2)],species_strains_T3_species_5,by="Studienummer")

species_strains_T3_long_2 <- species_strains_T3_species_6 %>%
  pivot_longer(
    cols = 5:ncol(species_strains_T3_species_6), # Columns 4 to the end
    names_to = "species", # New column for species names
    values_to = "abundance" # New column for abundance
  )

#species_strains_T3_long


df_normalized <- species_strains_T3_long_2 %>%
  arrange(age_weeks) %>% 
  mutate(Studienummer = factor(Studienummer, levels = unique(Studienummer[order(age_weeks)]))) %>%
  group_by(Studienummer, sample_type)


# Step 2: Create the bar plot, faceted by sample_type (fecal, vaginal, infant)


sample_delivery_type_plot_1 <- ggplot(df_normalized, aes(x = Studienummer, y = abundance, fill = species)) +
 geom_bar(stat = "identity") +
  facet_grid(sample_type ~ c_section,scales = "free_x") +
  theme_minimal() + 
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
  labs(x="samples (youngest to oldest)",y="relative abundance") +  
  scale_fill_manual(values = group.colors) 

#ggsave(filename = "/PATH/sample_delivery_type_plot_all_species_metaphlan_by_age_2.pdf", plot = sample_delivery_type_plot_1, width = 12, height = 8, dpi = 600)



species_strains_metaphlan_meta <- species_strains_metaphlan_meta %>%
  mutate(c_section = case_when(
    c_section == 0 ~ "vaginal_delivery",
    c_section == 1 ~ "c_section",
    TRUE ~ as.character(c_section)  # Keeps original values for any other cases
  ))

# merge on different delivery type
species_strains_metaphlan_meta_2 <- merge(mother_infant_submission_date[,c(1,2)],species_strains_metaphlan_meta,by="Studienummer")

df_normalized <- species_strains_metaphlan_meta_2 %>%
  arrange(age_weeks) %>% 
  mutate(Studienummer = factor(Studienummer, levels = unique(Studienummer[order(age_weeks)]))) %>%
  group_by(Studienummer, sample_type)

# Step 2: Create the bar plot, faceted by sample_type (fecal, vaginal, infant)

sample_delivery_type_plot <- ggplot(df_normalized, aes(x = Studienummer, y = abundance, fill = species)) +
 geom_bar(stat = "identity",position="fill" ) +
  facet_grid(sample_type ~ c_section,scales = "free_x") +
  theme_minimal() + 
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
  labs(x="samples (youngest to oldest)",y="relative abundance") +
  
  scale_fill_manual(values = group.colors) 

#ggsave(filename = "/PATH/sample_delivery_type_plot_species_metaphlan_by_age.pdf", plot = sample_delivery_type_plot, width = 12, height = 8, dpi = 600)
sample_delivery_type_plot


infant_strains_FEAST_studie <- merge(metadata[,c(1,507,341,349)],infant_strains_FEAST_species_names,by.x="kit3.infant_faecal_sample.barcode"
                                           ,by.y="infant_ID")
infant_strains_FEAST_studie <- infant_strains_FEAST_studie[,-c(1,5)]

colnames(infant_strains_FEAST_studie) <- c("Studienummer","c-section","breast_fed","V3","F3","Unknown","species")

infant_strains_FEAST_studie <- infant_strains_FEAST_studie[,c(1,4,5,6,2,3,7)]

#infant_strains_FEAST_studie <- All_strains_combined_infant_ordered

## changing the column name since R doesn't recognize -, replacing it with _
#All_strains_combined_infant_ordered_2 <- All_strains_combined_infant_ordered
#names(All_strains_combined_infant_ordered_2)[5] <- "delivery_type"


All_strains_combined_infant_ordered_2 <- infant_strains_FEAST_studie
names(All_strains_combined_infant_ordered_2)[2] <- "delivery_type"

## changing entries for delivery type
All_strains_combined_infant <- All_strains_combined_infant_ordered_2 %>%
  mutate(delivery_type = case_when(
    delivery_type == 0 ~ "vaginal_delivery",
    delivery_type == 1 ~ "c_section",
    TRUE ~ as.character(delivery_type)  # Keeps original values for any other cases
  ))

## changing names for breastfed type
All_strains_combined_infant <- All_strains_combined_infant %>%
  mutate(breast_fed = case_when(
    breast_fed == "Ja, men slutade" ~ "No",
    breast_fed == "Ja, ammar fortfarande" ~ "Yes",
    breast_fed == "Nej" ~ "No",
    TRUE ~ as.character(breast_fed)  # Keeps original values for any other cases
  ))

## adding an additional column for sum of V3 and F3
All_strains_combined_infant[,c(2,3,4)] <- sapply(All_strains_combined_infant[,c(2,3,4)],function(x) as.numeric(x))
All_strains_combined_infant$V3_F3 <- All_strains_combined_infant$V3+All_strains_combined_infant$F3
                                                 
All_strains_combined_infant <- All_strains_combined_infant[order(All_strains_combined_infant$species),]
                                                

# Choosing the entries for VF
All_strains_combined_infant_VF <- All_strains_combined_infant[,c(1,5,6,7,8)]

# Pivot the data into a wider format
All_strains_combined_infant_VF_wide <- All_strains_combined_infant_VF %>%
  pivot_wider(
    names_from = species,     # Columns created from "species"
    values_from = V3_F3  # Values come from the 'V3_F3' column
  )

## set all NAs to zero since not all samples have all species in them

All_strains_combined_infant_VF_wide[is.na(All_strains_combined_infant_VF_wide)] <- 0

## order by c_section and breast_fed so I can make the heatmap based on them
All_strains_combined_infant_VF_wide <- as.data.frame(All_strains_combined_infant_VF_wide)

All_strains_combined_infant_VF_wide <- All_strains_combined_infant_VF_wide[order(All_strains_combined_infant_VF_wide$delivery_type, All_strains_combined_infant_VF_wide$breast_fed),]

rownames(All_strains_combined_infant_VF_wide) <- All_strains_combined_infant_VF_wide$Studienummer


## now let's make the heatmap  
all_strain_infant_matrix <- as.matrix(All_strains_combined_infant_VF_wide[, -c(which(names(All_strains_combined_infant_VF_wide) == "delivery_type"), 
                                                                             which(names(All_strains_combined_infant_VF_wide) == "breast_fed"), 
                                                                             which(names(All_strains_combined_infant_VF_wide) == "Studienummer"))])

# Create annotation data frame for `A` and `B`
annotation <- All_strains_combined_infant_VF_wide[, c("breast_fed","delivery_type")]

## define colors
colors <- colorRampPalette(c("purple", "pink","white", "orange", "red"))(1000)

# Define breaks to make the heatmap more sensitive to lower values
breaks <- seq(min(all_strain_infant_matrix), max(all_strain_infant_matrix), length.out = 1001)

# Plot heatmap
mother_infant_strain_percent_pheatmap <- pheatmap(all_strain_infant_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,                                         
         annotation_row = annotation,
        show_rownames = FALSE,
        color = colors,
        breaks = breaks,
        main = "Percentage of infant strains from each species from mother")

#ggsave(filename = "/PATH/strain_percent_pheatmap.pdf", plot = mother_infant_strain_percent_pheatmap, width = 12, height = 8, dpi = 600)



# Choosing the entries for unknown
All_strains_combined_infant_unknown <- All_strains_combined_infant[,c(1,4,5,6,7)]

# Pivot the data into a wider format
All_strains_combined_infant_unknown_wide <- All_strains_combined_infant_unknown %>%
  pivot_wider(
    names_from = species,     # Columns created from "species"
    values_from = Unknown  # Values come from the 'Unknown' column
  )

## set all NAs to zero since not all samples have all species in them

All_strains_combined_infant_unknown_wide[is.na(All_strains_combined_infant_unknown_wide)] <- 0

## order by c_section and breast_fed so I can make the heatmap based on them
All_strains_combined_infant_unknown_wide <- as.data.frame(All_strains_combined_infant_unknown_wide)

All_strains_combined_infant_unknown_wide <- All_strains_combined_infant_unknown_wide[order(All_strains_combined_infant_unknown_wide$delivery_type, All_strains_combined_infant_unknown_wide$breast_fed),]

rownames(All_strains_combined_infant_unknown_wide) <- All_strains_combined_infant_unknown_wide$Studienummer


## now let's make the heatmap  
all_strain_infant_matrix <- as.matrix(All_strains_combined_infant_unknown_wide[, -c(which(names(All_strains_combined_infant_unknown_wide) == "delivery_type"), 
                                                                             which(names(All_strains_combined_infant_unknown_wide) == "breast_fed"), 
                                                                             which(names(All_strains_combined_infant_unknown_wide) == "Studienummer"))])

# Create annotation data frame for `A` and `B`
annotation <- All_strains_combined_infant_unknown_wide[, c("breast_fed","delivery_type")]

## define colors
colors <- colorRampPalette(c("purple", "pink","white", "orange", "red"))(1000)

# Define breaks to make the heatmap more sensitive to lower values
breaks <- seq(min(all_strain_infant_matrix), max(all_strain_infant_matrix), length.out = 1001)
# Plot heatmap
unknown_infant_strain_percent_pheatmap <- pheatmap(all_strain_infant_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,                                          
         annotation_row = annotation,
        show_rownames = FALSE,
        color = colors,
        breaks = breaks,
        main = "Percentage of infant strains from each species from sources other than mother")

#ggsave(filename = "/PATH/unknown_infant_strain_percent_pheatmap.pdf", plot = unknown_infant_strain_percent_pheatmap, width = 12, height = 8, dpi = 600)



# Choosing the entries for vaginal
All_strains_combined_infant_vaginal <- All_strains_combined_infant[,c(1,2,5,6,7)]

# Pivot the data into a wider format
All_strains_combined_infant_vaginal_wide <- All_strains_combined_infant_vaginal %>%
  pivot_wider(
    names_from = species,     # Columns created from "species"
    values_from = V3  # Values come from the 'vaginal' column
  )

## set all NAs to zero since not all samples have all species in them

All_strains_combined_infant_vaginal_wide[is.na(All_strains_combined_infant_vaginal_wide)] <- 0

## order by c_section and breast_fed so I can make the heatmap based on them
All_strains_combined_infant_vaginal_wide <- as.data.frame(All_strains_combined_infant_vaginal_wide)

All_strains_combined_infant_vaginal_wide <- All_strains_combined_infant_vaginal_wide[order(All_strains_combined_infant_vaginal_wide$delivery_type, All_strains_combined_infant_vaginal_wide$breast_fed),]

rownames(All_strains_combined_infant_vaginal_wide) <- All_strains_combined_infant_vaginal_wide$Studienummer


## now let's make the heatmap  
all_strain_infant_matrix <- as.matrix(All_strains_combined_infant_vaginal_wide[, -c(which(names(All_strains_combined_infant_vaginal_wide) == "delivery_type"), 
                                                                             which(names(All_strains_combined_infant_vaginal_wide) == "breast_fed"), 
                                                                                   which(names(All_strains_combined_infant_vaginal_wide) == "Studienummer"))])

# Create annotation data frame for `A` and `B`
annotation <- All_strains_combined_infant_vaginal_wide[, c("breast_fed","delivery_type")]

# making color palette more sensitive

colors <- colorRampPalette(c("purple", "pink","white", "orange", "red"))(1000)

# Define breaks to make the heatmap more sensitive to lower values
breaks <- seq(min(all_strain_infant_matrix), max(all_strain_infant_matrix), length.out = 1001)


# Plot heatmap
vaginal_infant_strain_percent_pheatmap <- pheatmap(all_strain_infant_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = annotation,
         color = colors,    # Apply custom color palette
         breaks = breaks,
        show_rownames = FALSE,
        main = "Percentage of infant strains from each species from mother vaginal")

#ggsave(filename = "/PATH/vaginal_infant_strain_percent_pheatmap.pdf", plot = vaginal_infant_strain_percent_pheatmap, width = 12, height = 8, dpi = 600)


# Choosing the entries for motherfecal
All_strains_combined_infant_motherfecal <- All_strains_combined_infant[,c(1,3,5,6,7)]

# Pivot the data into a wider format
All_strains_combined_infant_motherfecal_wide <- All_strains_combined_infant_motherfecal %>%
  pivot_wider(
    names_from = species,     # Columns created from "species"
    values_from = F3  # Values come from the 'motherfecal' column
  )

## set all NAs to zero since not all samples have all species in them

All_strains_combined_infant_motherfecal_wide[is.na(All_strains_combined_infant_motherfecal_wide)] <- 0

## order by c_section and breast_fed so I can make the heatmap based on them
All_strains_combined_infant_motherfecal_wide <- as.data.frame(All_strains_combined_infant_motherfecal_wide)

All_strains_combined_infant_motherfecal_wide <- All_strains_combined_infant_motherfecal_wide[order(All_strains_combined_infant_motherfecal_wide$delivery_type, All_strains_combined_infant_motherfecal_wide$breast_fed),]

rownames(All_strains_combined_infant_motherfecal_wide) <- All_strains_combined_infant_motherfecal_wide$Studienummer


## now let's make the heatmap  
all_strain_infant_matrix <- as.matrix(All_strains_combined_infant_motherfecal_wide[, -c(which(names(All_strains_combined_infant_motherfecal_wide) == "delivery_type"), 
                                                                             which(names(All_strains_combined_infant_motherfecal_wide) == "breast_fed"), 
                                                                                   which(names(All_strains_combined_infant_motherfecal_wide) == "Studienummer"))])

# Create annotation data frame for `A` and `B`
annotation <- All_strains_combined_infant_motherfecal_wide[, c("breast_fed","delivery_type")]

# making color palette more sensitive

colors <- colorRampPalette(c("purple", "pink","white", "orange", "red"))(1000)

# Define breaks to make the heatmap more sensitive to lower values
breaks <- seq(min(all_strain_infant_matrix), max(all_strain_infant_matrix), length.out = 1001)


# Plot heatmap
motherfecal_infant_strain_percent_pheatmap <- pheatmap(all_strain_infant_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = annotation,
         color = colors,    # Apply custom color palette
         breaks = breaks,
        show_rownames = FALSE,
        main = "Percentage of infant strains from each species from mother fecal samples")

#ggsave(filename = "/ceph/projects/010_SweMaMi/analyses/nicole/strain_tracking/plots/motherfecal_infant_strain_percent_pheatmap.pdf", plot = motherfecal_infant_strain_percent_pheatmap, width = 12, height = 8, dpi = 600)

 

## making numeric columns in this case binary..- can use this vars(col1, col2, col3) instead of is.numeric if I don't want to convert all numeric columns
All_strains_combined_infant_binary <- All_strains_combined_infant %>%
  mutate_if(is.numeric, ~ ifelse(. > 0, 1, 0))


all_strains_modified_binary <- All_strains_combined_infant_binary %>%
  mutate(
    mother_present = ifelse(V3_F3 > 0, "mother_present", "mother_absent"),
    vaginal_present = ifelse(V3 > 0, "mother_vaginal_present", "mother_vaginal_absent"),
    fecal_present = ifelse(F3 > 0, "mother_fecal_present", "mother_fecal_absent"),
    unknown_present = ifelse(Unknown > 0, "unknown_present", "unknown_absent"),
    infant_present = rep("infant_present",times=nrow(All_strains_combined_infant_binary))  
  )


all_strains_modified <- All_strains_combined_infant %>%
  mutate(
    mother_present = ifelse(V3_F3 > 0, "mother_present", "mother_absent"),
    vaginal_present = ifelse(V3 > 0, "mother_vaginal_present", "mother_vaginal_absent"),
    fecal_present = ifelse(F3 > 0, "mother_fecal_present", "mother_fecal_absent"),
    unknown_present = ifelse(Unknown > 0, "unknown_present", "unknown_absent"),
    infant_present = rep("infant_present",times=nrow(All_strains_combined_infant))  
  )

### separating by delivery type
c_section_strains_combined_infant_binary <- filter(all_strains_modified_binary, delivery_type %in% c("c_section"))
c_section_strains_combined_infant <- filter(all_strains_modified, delivery_type %in% c("c_section"))

vaginal_delivery_strains_combined_infant_binary <- filter(all_strains_modified_binary, delivery_type %in% c("vaginal_delivery"))
vaginal_delivery_strains_combined_infant <- filter(all_strains_modified, delivery_type %in% c("vaginal_delivery"))

#### 2 column

alluvial_c_section_2_column <- ggplot(c_section_strains_combined_infant_binary, aes(axis1 = Studienummer, axis2 = mother_present, y = Unknown)) +
  geom_alluvium(aes(fill = species)) +
  geom_stratum(width=1/3) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Sample", "present in mother"), expand = c(0.15, 0.02)) +
  ggtitle("c_section infant strain presence in maternal vaginal and fecal samples")

#ggsave(filename = "/PATH/alluvial_c_section_2_column_presence_absence.pdf", plot = alluvial_c_section_2_column, width = 12, height = 8, dpi = 600)


alluvial_vaginal_delivery_2_column <- ggplot(vaginal_delivery_strains_combined_infant_binary, aes(axis1 = Studienummer, axis2 = mother_present, y = Unknown)) +
  geom_alluvium(aes(fill = species)) +
  geom_stratum(width=1/3) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Sample", "present in mother"), expand = c(0.15, 0.02)) +
  ggtitle("vaginal delivery infant strain presence in maternal samples")

#ggsave(filename = "/PATH/alluvial_vaginal_delivery_2_column_presence_absence.pdf", plot = alluvial_vaginal_delivery_2_column, width = 12, height = 8, dpi = 600)




#### 3 column

alluvial_c_section_3_column <- ggplot(c_section_strains_combined_infant, aes(axis1 = vaginal_present, axis2 = Studienummer, axis3 = fecal_present, y = V3_F3)) +
  geom_alluvium(aes(fill = species)) +
  geom_stratum(width=1/2) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("present in mother vaginal","Sample", "present in mother fecal"), expand = c(0.15, 0.02)) +
  ggtitle("c_section infant strain amount in maternal vaginal and fecal samples")

#ggsave(filename = "/PATH/alluvial_c_section_3_column.pdf", plot = alluvial_c_section_3_column, width = 12, height = 8, dpi = 600)


alluvial_vaginal_delivery_3_column <- ggplot(vaginal_delivery_strains_combined_infant, aes(axis1 = vaginal_present, axis2 = Studienummer, axis3 = fecal_present, y = V3_F3)) +
  geom_alluvium(aes(fill = species)) +
  geom_stratum(width=1/2) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("present in mother vaginal","Sample", "present in mother fecal"), expand = c(0.15, 0.02)) +
  ggtitle("vaginal delivery infant strain amount in maternal vaginal and fecal samples")

#ggsave(filename = "/PATH/alluvial_vaginal_delivery_3_column.pdf", plot = alluvial_vaginal_delivery_3_column, width = 12, height = 8, dpi = 600)



# separating the two types of Bifido from the others
strains_infant_filtered_bifi_longum <- filter(All_strains_combined_infant,grepl("longum", species, ignore.case = TRUE))
strains_infant_filtered_bifi_breve <- filter(All_strains_combined_infant,grepl("breve", species, ignore.case = TRUE))
strains_infant_filtered_bifi_breve_longum <- rbind(strains_infant_filtered_bifi_longum,strains_infant_filtered_bifi_breve)

# separating C_section from vaginal delivery

strains_infant_filtered_bifi_breve_longum_c_section <- filter(strains_infant_filtered_bifi_breve_longum,
                                                              delivery_type %in% c("c_section"))

strains_infant_filtered_bifi_breve_longum_vaginal <- filter(strains_infant_filtered_bifi_breve_longum,
                                                              delivery_type %in% c("vaginal_delivery"))


# Permutation testing 
set.seed(123)
n_bootstrap <- 10000
combined <- c(c_section, vaginal_delivery)
group_labels <- c(rep(1, length(c_section)), rep(2, length(vaginal_delivery)))

observed_diff_perm <- mean(c_section) - mean(vaginal_delivery)
perm_diffs <- numeric(n_bootstrap)

for (i in 1:n_bootstrap) {
  shuffled_labels <- sample(group_labels)
  perm_diffs[i] <- mean(combined[shuffled_labels == 1]) - mean(combined[shuffled_labels == 2])
}

p_value_perm <- mean(abs(perm_diffs) >= abs(observed_diff_perm))

cat("Permutation p-value:", p_value_perm, "\n")


### mother strains
infant_filtered_bifi_mother_plot <- ggplot(strains_infant_filtered_bifi_breve_longum, aes(x = delivery_type, y = V3_F3)) +
  geom_violin(fill="yellow")+
#  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="purple") +
  ggtitle("B. longum-breve c_section vs. vaginal mother, permutated p.val = 0.0219")

#ggsave(filename = "/PATH/strain_infant_filtered_bifi_longum_breve_mother_plot_premutated.pdf", plot = infant_filtered_bifi_mother_plot, width = 12, height = 8, dpi = 600)


### unknown strains
infant_filtered_bifi_unknown_plot <- ggplot(strains_infant_filtered_bifi_breve_longum, aes(x = delivery_type, y = Unknown)) +
  geom_violin(fill="yellow")+
#  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="purple") +
  ggtitle("B. longum-breve c_section vs. vaginal unknown, permutated p.val = 0.0207")

#ggsave(filename = "/PATH/strains_infant_filtered_bifi_longum_breve_unknown_plot_permutated.pdf", plot = infant_filtered_bifi_unknown_plot, width = 12, height = 8, dpi = 600)



strains_infant_filtered_non_bifi <- filter(All_strains_combined_infant,!grepl("bifi", species, ignore.case = TRUE))

strains_infant_filtered_non_bifi_c_section <- filter(strains_infant_filtered_non_bifi,
                                                              delivery_type %in% c("c_section"))

strains_infant_filtered_non_bifi_vaginal <- filter(strains_infant_filtered_non_bifi, 
                                          delivery_type %in% c("vaginal_delivery"))
 


# Permutation testing 
set.seed(123)
n_bootstrap <- 10000
combined <- c(c_section, vaginal_delivery)
group_labels <- c(rep(1, length(c_section)), rep(2, length(vaginal_delivery)))

observed_diff_perm <- mean(c_section) - mean(vaginal_delivery)
perm_diffs <- numeric(n_bootstrap)

for (i in 1:n_bootstrap) {
  shuffled_labels <- sample(group_labels)
  perm_diffs[i] <- mean(combined[shuffled_labels == 1]) - mean(combined[shuffled_labels == 2])
}

p_value_perm <- mean(abs(perm_diffs) >= abs(observed_diff_perm))

cat("Permutation p-value:", p_value_perm, "\n") 


c_section_df <- filter(All_strains_combined_infant, delivery_type %in% c("c_section"))
vaginal_df <- filter(All_strains_combined_infant, delivery_type %in% c("vaginal_delivery"))


### mother strains
infant_filtered_non_bifi_mother_plot <- ggplot(strains_infant_filtered_non_bifi, aes(x = delivery_type, y = V3_F3)) +
  geom_violin(fill="lightblue")+
#  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="black") +
  ggtitle("non_bifido c_section vs. vaginal mother, permutated p.val = 0.0601")

#ggsave(filename = "/ceph/projects/010_SweMaMi/analyses/nicole/strain_tracking/plots/strain_infant_filtered_non_bifi_mother_plot_permutated.pdf", plot = infant_filtered_non_bifi_mother_plot, width = 12, height = 8, dpi = 600)
d
### unknown strains
infant_filtered_non_bifi_unknown_plot <- ggplot(strains_infant_filtered_non_bifi, aes(x = delivery_type, y = Unknown)) +
  geom_violin(fill="lightblue")+
#  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="black") +
  ggtitle("non_bifido c_section vs. vaginal unknown, permutated p.val = 0.0544")

#ggsave(filename = "/ceph/projects/010_SweMaMi/analyses/nicole/strain_tracking/plots/strains_infant_filtered_non_bifi_unknown_plot_permutated.pdf", plot = infant_filtered_non_bifi_unknown_plot, width = 12, height = 8, dpi = 600)



# Permutation testing is better for these purposes cause we want to challenge the null hypothesis of
# whether or not there is indeed a difference between the two groups.
set.seed(123)
n_bootstrap <- 10000
combined <- c(c_section, vaginal_delivery)
group_labels <- c(rep(1, length(c_section)), rep(2, length(vaginal_delivery)))

observed_diff_perm <- mean(c_section) - mean(vaginal_delivery)
perm_diffs <- numeric(n_bootstrap)

for (i in 1:n_bootstrap) {
  shuffled_labels <- sample(group_labels)
  perm_diffs[i] <- mean(combined[shuffled_labels == 1]) - mean(combined[shuffled_labels == 2])
}

p_value_perm <- mean(abs(perm_diffs) >= abs(observed_diff_perm))

cat("Permutation p-value:", p_value_perm, "\n")

# All strains
# all strains
strains_infant_filtered_c_section <- filter(All_strains_combined_infant,
                                                              delivery_type %in% c("c_section"))

strains_infant_filtered_vaginal <- filter(All_strains_combined_infant, 
                                          delivery_type %in% c("vaginal_delivery"))
 
### mother strains
infant_filtered_all_mother_plot <- ggplot(All_strains_combined_infant, aes(x = delivery_type, y = F3)) +
  geom_violin(fill="lightblue")+
 # geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="black") +
  ggtitle("all strains c_section vs. vaginal delivery from mother, p.val = 0.0693")

#ggsave(filename = "/PATH/strain_infant_filtered_all_mother_plot_permutated.pdf", plot = infant_filtered_all_mother_plot, width = 12, height = 8, dpi = 600)
infant_filtered_all_mother_plot

### unknown strains
infant_filtered_all_unknown_plot <- ggplot(All_strains_combined_infant, aes(x = delivery_type, y = Unknown)) +
  geom_violin(fill="lightblue")+
 # geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="black") +
  ggtitle("all strains c_section vs. vaginal unknown, permutated p.val = 0.0668")

#ggsave(filename = "/PATH/strains_infant_filtered_all_unknown_plot_permutated.pdf", plot = infant_filtered_all_unknown_plot, width = 12, height = 8, dpi = 600)
infant_filtered_all_unknown_plot


# only get the strains

species_strains_F3_F2_F1_t_names <- grep("\\.t_",names(species_strains_F3_F2_F1),value=TRUE)
species_strains_F3_F2_F1_t <- species_strains_F3_F2_F1[,species_strains_F3_F2_F1_t_names]
species_strains_F3_F2_F1_t <- cbind(species_strains_F3_F2_F1[,c(1:3,ncol(species_strains_F3_F2_F1))],species_strains_F3_F2_F1_t)
# separating the time points for fecal
species_strains_F3 <- filter(species_strains_F3_F2_F1_t,grepl("F3", specifier))
species_strains_F2 <- filter(species_strains_F3_F2_F1_t,grepl("F2", specifier))
species_strains_F1 <- filter(species_strains_F3_F2_F1_t,grepl("F1", specifier))

# getting specific time points, ordering and removing extra rows F1 F2
species_strains_F1_F2 <- rbind(species_strains_F1,species_strains_F2)
species_strains_F1_F2_ordered <- species_strains_F1_F2[order(species_strains_F1_F2$specifier),]
F1_F2_filtered <- species_strains_F1_F2_ordered %>%
  group_by(Studienummer) %>%
  filter(n_distinct(specifier) > 1) %>%
  ungroup()
species_strains_F1_F2_ordered_2 <- F1_F2_filtered[,-c(1:3)]

# getting specific time points, ordering and removing extra rows F1 F3
species_strains_F1_F3 <- rbind(species_strains_F1,species_strains_F3)
species_strains_F1_F3_ordered <- species_strains_F1_F3[order(species_strains_F1_F3$specifier),]
F1_F3_filtered <- species_strains_F1_F3_ordered %>%
  group_by(Studienummer) %>%
  filter(n_distinct(specifier) > 1) %>%
  ungroup()
species_strains_F1_F3_ordered_2 <- F1_F3_filtered[,-c(1:3)]

# getting specific time points, ordering and removing extra rows F2 F3
species_strains_F2_F3 <- rbind(species_strains_F2,species_strains_F3)
species_strains_F2_F3_ordered <- species_strains_F2_F3[order(species_strains_F2_F3$specifier),]
F2_F3_filtered <- species_strains_F2_F3_ordered %>%
  group_by(Studienummer) %>%
  filter(n_distinct(specifier) > 1) %>%
  ungroup()
species_strains_F2_F3_ordered_2 <- F2_F3_filtered[,-c(1:3)]
   

# combination

species_strains_F3_F2_F1_redone <- species_strains_F3_F2_F1[,-c(1,3)]

species_strains_F3_F2_F1_redone <- species_strains_F3_F2_F1_redone %>%
  mutate(timepoint = str_extract(specifier, "(?<=_).*"))


set.seed(123)
# Remove the sample column
F1_F2_strains_only <- species_strains_F1_F2_ordered_2[,-1]

# Calculate Bray-Curtis dissimilarity for all samples
F1_F2_strains_bray_curtis_dissimilarity <- vegdist(F1_F2_strains_only, method = "bray")
dim(F1_F2_strains_bray_curtis_dissimilarity)
dissim_F1_F2 <- as.vector(as.dist(F1_F2_strains_bray_curtis_dissimilarity))


#
samples <- unique(sub("_F[0-9]+", "", species_strains_F1_F2_ordered_2$specifier)) # Extract sample IDs without time points
F1_F2_strains_bray_curtis_dissimilarity_matrix <- as.matrix(F1_F2_strains_bray_curtis_dissimilarity)

F1_F2_dissimilarity_values <- c()
sample_names <- c()

for (sample in samples) {
  # Get the indices of the two time points for each sample
  sample_indices <- grep(sample, species_strains_F1_F2_ordered_2$specifier)
  
  # Check if there are exactly 2 time points for the sample
  if (length(sample_indices) == 2) {
    # Extract the pair of dissimilarities between time points
    F1_F2_dissimilarity <- F1_F2_strains_bray_curtis_dissimilarity_matrix[sample_indices[1], sample_indices[2]]
    
    
    F1_F2_dissimilarity_values <- c(F1_F2_dissimilarity_values, F1_F2_dissimilarity)
    sample_names <- c(sample_names, sample)
  } else {
    print(paste("Sample", sample, "does not have exactly 2 time points. Skipping..."))
  }
}

F1_F2_dissimilarity_df <- data.frame(Studiennumer = sample_names,Time_point = rep("F1_F2",times=100), bray_curtis_dissimilarity = F1_F2_dissimilarity_values)



set.seed(123)
# Remove the sample column
F1_F3_strains_only <- species_strains_F1_F3_ordered_2[,-1]

# Calculate Bray-Curtis dissimilarity for all samples
F1_F3_strains_bray_curtis_dissimilarity <- vegdist(F1_F3_strains_only, method = "bray")
dim(F1_F3_strains_bray_curtis_dissimilarity)
dissim_F1_F3 <- as.vector(as.dist(F1_F3_strains_bray_curtis_dissimilarity))

##

samples <- unique(sub("_F[0-9]+", "", species_strains_F1_F3_ordered_2$specifier)) # Extract sample IDs without time points
F1_F3_strains_bray_curtis_dissimilarity_matrix <- as.matrix(F1_F3_strains_bray_curtis_dissimilarity)

F1_F3_dissimilarity_values <- c()
sample_names <- c()

for (sample in samples) {
  # Get the indices of the two time points for each sample
  sample_indices <- grep(sample, species_strains_F1_F3_ordered_2$specifier)
  
  # Check if there are exactly 2 time points for the sample
  if (length(sample_indices) == 2) {
    # Extract the pair of dissimilarities between time points
    F1_F3_dissimilarity <- F1_F3_strains_bray_curtis_dissimilarity_matrix[sample_indices[1], sample_indices[2]]
    
    
    F1_F3_dissimilarity_values <- c(F1_F3_dissimilarity_values, F1_F3_dissimilarity)
    sample_names <- c(sample_names, sample)
  } else {
    print(paste("Sample", sample, "does not have exactly 2 time points. Skipping..."))
  }
}

F1_F3_dissimilarity_df <- data.frame(Studiennumer = sample_names,Time_point = rep("F1_F3",times=98), bray_curtis_dissimilarity = F1_F3_dissimilarity_values)



set.seed(123)
# Remove the sample column
F2_F3_strains_only <- species_strains_F2_F3_ordered_2[,-1]

# Calculate Bray-Curtis dissimilarity for all samples
F2_F3_strains_bray_curtis_dissimilarity <- vegdist(F2_F3_strains_only, method = "bray")
dim(F2_F3_strains_bray_curtis_dissimilarity)
dissim_F2_F3 <- as.vector(as.dist(F2_F3_strains_bray_curtis_dissimilarity))

##

samples <- unique(sub("_F[0-9]+", "", species_strains_F2_F3_ordered_2$specifier)) # Extract sample IDs without time points
F2_F3_strains_bray_curtis_dissimilarity_matrix <- as.matrix(F2_F3_strains_bray_curtis_dissimilarity)

F2_F3_dissimilarity_values <- c()
sample_names <- c()

for (sample in samples) {
  # Get the indices of the two time points for each sample
  sample_indices <- grep(sample, species_strains_F2_F3_ordered_2$specifier)
  
  # Check if there are exactly 2 time points for the sample
  if (length(sample_indices) == 2) {
    # Extract the pair of dissimilarities between time points
    F2_F3_dissimilarity <- F2_F3_strains_bray_curtis_dissimilarity_matrix[sample_indices[1], sample_indices[2]]
    
    
    F2_F3_dissimilarity_values <- c(F2_F3_dissimilarity_values, F2_F3_dissimilarity)
    sample_names <- c(sample_names, sample)
  } else {
    print(paste("Sample", sample, "does not have exactly 2 time points. Skipping..."))
  }
}

F2_F3_dissimilarity_df <- data.frame(Studiennumer = sample_names,Time_point = rep("F2_F3",times=100), bray_curtis_dissimilarity = F2_F3_dissimilarity_values)


samples <- unique(sub("_F[0-9]+", "", species_strains_F2_F3_ordered_2$specifier)) # Extract sample IDs without time points
F2_F3_strains_bray_curtis_dissimilarity_matrix <- as.matrix(F2_F3_strains_bray_curtis_dissimilarity)

F2_F3_dissimilarity_values <- c()
sample_names <- c()

for (sample in samples) {
  # Get the indices of the two time points for each sample
  sample_indices <- grep(sample, species_strains_F2_F3_ordered_2$specifier)
  
  # Check if there are exactly 2 time points for the sample
  if (length(sample_indices) == 2) {
    # Extract the pair of dissimilarities between time points
    F2_F3_dissimilarity <- F2_F3_strains_bray_curtis_dissimilarity_matrix[sample_indices[1], sample_indices[2]]
    
    
    F2_F3_dissimilarity_values <- c(F2_F3_dissimilarity_values, F2_F3_dissimilarity)
    sample_names <- c(sample_names, sample)
  } else {
    print(paste("Sample", sample, "does not have exactly 2 time points. Skipping..."))
  }
}

F2_F3_dissimilarity_df <- data.frame(Studiennumer = sample_names,Time_point = rep("F2_F3",times=100), bray_curtis_dissimilarity = F2_F3_dissimilarity_values)



### finding p_values between them
set.seed(123)
F1_F2_with_F1_F3_p_val <- t.test(F1_F2_dissimilarity_df$bray_curtis_dissimilarity,F1_F3_dissimilarity_df$bray_curtis_dissimilarity)

F1_F2_with_F2_F3_p_val <- t.test(F1_F2_dissimilarity_df$bray_curtis_dissimilarity,F2_F3_dissimilarity_df$bray_curtis_dissimilarity)

F1_F3_with_F2_F3_p_val <- t.test(F1_F3_dissimilarity_df$bray_curtis_dissimilarity,F2_F3_dissimilarity_df$bray_curtis_dissimilarity)


dissim_data_F <- rbind(F2_F3_dissimilarity_df,F1_F3_dissimilarity_df,F1_F2_dissimilarity_df)


Fecal_dissimilarity_plot <- ggplot(dissim_data_F, aes(x = Time_point, y = bray_curtis_dissimilarity)) +
  geom_violin(fill="lightpink")+
#  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_boxplot(width=0.1,color="black") +
  ggtitle("Fecal dissimilarity plot \n p_val F1_F2 vs F2_F3 = 0.014, \n F1_F2 vs F1_F3 = 0.02, \n F2_F3 vs F1_F3 = 0.712)")


#ggsave(filename="/PATH/Fecal_dissimilarity_plot.pdf", plot=Fecal_dissimilarity_plot, width = 12, height = 8, dpi = 600)


species_strains_F3_F2_F1_redone <- species_strains_F3_F2_F1[,-c(1,3)]

species_strains_F3_F2_F1_redone <- species_strains_F3_F2_F1_redone %>%
  mutate(timepoint = str_extract(specifier, "(?<=_).*"))

# cleaning

F_names_2 <- c("Alicyclobacillus_acidiphilus","Alicyclobacillus_acidiphilus.t__SGB30417")

species_strains_F3_F2_F1_redone <- species_strains_F3_F2_F1_redone[,! colnames(species_strains_F3_F2_F1_redone) %in% F_names_2]

species_strains_F3_F2_F1_redone_2 <- species_strains_F3_F2_F1_redone
F_colnames <- colnames(species_strains_F3_F2_F1_redone)[grepl("t__", colnames(species_strains_F3_F2_F1_redone))]
species_strains_F3_F2_F1_redone<- species_strains_F3_F2_F1_redone[, colnames(species_strains_F3_F2_F1_redone) %in% F_colnames]

## adonis three time points
set.seed(1234)

adonis_file_species <- species_strains_F3_F2_F1_redone#[,c(2:(ncol(species_strains_F3_F2_F1_redone)-2))] ### the section that only contains the species + the metadata we want to look into
adonis_file_species_2 <- adonis_file_species + 0.0000001
adonis_file_species_clr <- clr(adonis_file_species_2)
adonis_file_species_clr <- as.data.frame(adonis_file_species_clr)
adonis_stra_timepoint <- adonis2(adonis_file_species ~ timepoint , data = species_strains_F3_F2_F1_redone_2, na.action = na.omit, permutations = 999,method="bray")

## two time points
set.seed(12345)

species_strains_F2_F1_redone <- filter(species_strains_F3_F2_F1_redone_3,timepoint %in% c("F1","F2"))
adonis_file_species <- species_strains_F2_F1_redone[,c(2:(ncol(species_strains_F2_F1_redone)-2))] ### the section that only contains the species + the metadata we want to look into
adonis_file_species_2 <- adonis_file_species + 0.0000001
adonis_file_species_clr <- clr(adonis_file_species_2)
adonis_file_species_clr <- as.data.frame(adonis_file_species_clr)

adonis_stra_timepoint <- adonis2(adonis_file_species ~ timepoint , data = species_strains_F2_F1_redone, na.action = na.omit, permutations = 999,method="bray")

# F1 and F3

set.seed(1234)
species_strains_F3_F1_redone <- filter(species_strains_F3_F2_F1_redone,timepoint %in% c("F1","F3"))
adonis_file_species <- species_strains_F2_F1_redone[,c(2:(ncol(species_strains_F3_F1_redone)-2))] ### the section that only contains the species + the metadata we want to look into
adonis_file_species_2 <- adonis_file_species + 0.0000001
adonis_file_species_clr <- clr(adonis_file_species_2)
adonis_file_species_clr <- as.data.frame(adonis_file_species_clr)
adonis_stra_timepoint <- adonis2(adonis_file_species ~ timepoint , data = species_strains_F2_F1_redone, na.action = na.omit, permutations = 9999,method="bray")

adonis_stra_timepoint

dist_matrix <- vegdist(adonis_file_species, method = "bray")

# Perform PCoA
pcoa_result <- cmdscale(dist_matrix)

# Convert to dataframe
pcoa_df <- data.frame(PC1 = pcoa_result[,1], 
                      PC2 = pcoa_result[,2],
                      Group = species_strains_F2_F1_redone$timepoint)

# Plot
plot <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  stat_ellipse(type = "t") +
  theme_minimal() +
  labs(title = "PCoA Plot (Bray-Curtis)", x = "PC1", y = "PC2")

#ggsave(filename="/PATH/fecal_PCA_plot.pdf", plot=plot, width = 12, height = 8, dpi = 600)

plot
dispersion <- betadisper(dist_matrix, species_strains_F2_F1_redone$timepoint)
#pdf("/PATH/fecal_all_dispersion.pdf")

boxplot(dispersion, main = "Beta Dispersion of fecal groups")

#dev.off()



species_strains_F2_F3_redone <- filter(species_strains_F3_F2_F1_redone_2,timepoint %in% c("F2","F3"))
metadata_for_aldex <- species_strains_F2_F3_redone[,c("specifier","timepoint")]
rownames(metadata_for_aldex) <- metadata_for_aldex$specifier
metadata_for_aldex <- metadata_for_aldex[,2,drop=FALSE]

aldex.clr <- aldex.clr(pseudo_counts, conds = metadata_for_aldex$timepoint, mc.samples = 128, denom = "all")

grouping <- metadata_for_aldex$timepoint  # e.g., Case vs. Control

# Run ALDEx2
aldex_result_F2_F3 <- aldex(pseudo_counts, 
                      conditions = grouping, 
                      mc.samples = 128,  # Monte Carlo Dirichlet instances
                      test = "t",  # "t" for Welch's t-test, "wilcox" for Wilcoxon
                      effect = TRUE,  # Calculate effect sizes
                      denom = "iqlr")  # Use IQLR normalization

# plot

effect <- ggplot(aldex_result_F1_F2_2, 
                 aes(x = diff.btw, y = -log10(we.ep))) +
  geom_point(aes(color = (abs(diff.btw) > 1.5 & we.ep < 0.1))) +
  theme_minimal() +
  labs(title = "ALDEx2 F1-F2 Volcano Plot", 
       x = "Log fold change", 
       y = "Welch's p-value")
#ggsave(filename = "/PATH/effect_V2_V3_aldex.pdf", plot = effect, width = 12, height = 8, dpi = 600)


F3_F2_F1_strains_FEAST_studie <- merge(metadata[,c(1,508)],F3_F2_F1_strains_FEAST_species_names,by.x="kit3.faecal_sample.barcode"
                                           ,by.y="F3_ID")
F3_F2_F1_strains_FEAST_studie <- F3_F2_F1_strains_FEAST_studie[,-c(1,3)]

colnames(F3_F2_F1_strains_FEAST_studie) <- c("Studienummer","F2","F1","Unknown","species")
F3_F2_F1_strains_FEAST_studie[,c(2,3,4)] <- sapply(F3_F2_F1_strains_FEAST_studie[,c(2,3,4)],function(x) as.numeric(x))
                                                 

All_strains_combined_F <- F3_F2_F1_strains_FEAST_studie[,-c(3,4)]
All_strains_combined_F <- All_strains_combined_F[order(All_strains_combined_F$species,decreasing = TRUE),]

# Pivot the data into a wider format
All_strains_combined_F_wide <- All_strains_combined_F %>%
  pivot_wider(
    names_from = species,     # Columns created from "species"
    values_from = F2  # Values come from the 'V3_F3' column
  )
All_strains_combined_F_wide[is.na(All_strains_combined_F_wide)] <- 0

## order by c_section and breast_fed so I can make the heatmap based on them
All_strains_combined_F_wide <- as.data.frame(All_strains_combined_F_wide)

#All_strains_combined_F_wide <- All_strains_combined_F_wide[order(All_strains_combined_F_wide$delivery_type, All_strains_combined_F_wide$breast_fed),]

rownames(All_strains_combined_F_wide) <- All_strains_combined_F_wide$Studienummer

## now let's make the heatmap  
all_strain_F_matrix <- as.matrix(All_strains_combined_F_wide[,-1])


## define colors
colors <- colorRampPalette(c("purple", "pink","white", "orange", "red"))(1000)

# Define breaks to make the heatmap more sensitive to lower values
breaks <- seq(min(all_strain_F_matrix), max(all_strain_F_matrix), length.out = 1001)

# Plot heatmap
F3_F2_strain_percent_pheatmap <- pheatmap(all_strain_F_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,                                         
      #   annotation_row = annotation,
        show_rownames = FALSE,
        color = colors,
        breaks = breaks,
        main = "Fecal transfer to TP3 from TP2")

#ggsave(filename = "/PATH/F3_from_F2_strain_percent_pheatmap.pdf", plot = F3_F2_strain_percent_pheatmap, width = 12, height = 8, dpi = 600)


All_strains_combined_F <- F3_F2_F1_all_strains_ordered[,-c(2,4)]
All_strains_combined_F <- All_strains_combined_F[order(All_strains_combined_F$species,decreasing = TRUE),]

# Pivot the data into a wider format
All_strains_combined_F_wide <- All_strains_combined_F %>%
  pivot_wider(
    names_from = species,     # Columns created from "species"
    values_from = F1  # Values come from the 'V3_F3' column
  )
All_strains_combined_F_wide[is.na(All_strains_combined_F_wide)] <- 0

## order by c_section and breast_fed so I can make the heatmap based on them
All_strains_combined_F_wide <- as.data.frame(All_strains_combined_F_wide)


rownames(All_strains_combined_F_wide) <- All_strains_combined_F_wide$Studienummer

## now let's make the heatmap  
all_strain_F_matrix <- as.matrix(All_strains_combined_F_wide[,-1])


## define colors
colors <- colorRampPalette(c("purple", "pink","white", "orange", "red"))(1000)

# Define breaks to make the heatmap more sensitive to lower values
breaks <- seq(min(all_strain_F_matrix), max(all_strain_F_matrix), length.out = 1001)

# Plot heatmap
F3_F1_strain_percent_pheatmap <- pheatmap(all_strain_F_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,                                         
      #   annotation_row = annotation,
        show_rownames = FALSE,
        color = colors,
        breaks = breaks,
        main = "Fecal transfer to TP3 from TP1")

#ggsave(filename = "/PATH/F3_from_F1_strain_percent_pheatmap.pdf", plot = F3_F1_strain_percent_pheatmap, width = 12, height = 8, dpi = 600)



F_strains_wide <- F_species_strains_metaphlan_combined %>%
  pivot_wider(names_from = sample_type, values_from = abundance)


F_strains_wide_presence <- F_strains_wide %>%
    mutate(F_TP1_presence = ifelse(F_TP1>0,1,0),
           F_TP2_presence = ifelse(F_TP2>0,1,0),
           F_TP3_presence = ifelse(F_TP3>0,1,0))

#F_strains_wide_presence

F_strain_data_upset <- F_strains_wide_presence %>%
  mutate(T1_T2 = ifelse(F_TP1_presence == 1 & F_TP2_presence == 1, 1, 0),
         T2_T3 = ifelse(F_TP2_presence == 1 & F_TP3_presence == 1, 1, 0)) 


#library(UpSetR)

# Function 
make_upset_plots <- function(df, species_name, species_label, outdir = "/PATH") {
  # Filter the data
  df_species <- filter(df, species %in% c(species_name))
  
  # --- 3TP plot ---
  pdf(file.path(outdir, paste0("F_strain_data_upset_", gsub(" ", "_", tolower(species_label)), "_3TP.pdf")))
  upset(as.data.frame(df_species),
        sets = c("F_TP1_presence", "F_TP2_presence", "F_TP3_presence"),
        order.by = "freq", keep.order = TRUE, main.bar.color = "steelblue",
        mainbar.y.label = paste(species_label, "presence"),
        text.scale = c(1.3, 1.3, 1, 1, 2, 2))
  dev.off()
  
  # --- 2TP plot ---
  pdf(file.path(outdir, paste0("F_strain_data_upset_", gsub(" ", "_", tolower(species_label)), "_2TP.pdf")))
  upset(as.data.frame(df_species),
        sets = c("T1_T2", "T2_T3"),
        order.by = "freq", keep.order = TRUE, main.bar.color = "steelblue",
        mainbar.y.label = paste(species_label, "presence"),
        text.scale = c(1.3, 1.3, 1, 1, 2, 2))
  dev.off()
}

# -----------------
#calls:
make_upset_plots(F_strain_data_upset, "Bifidobacterium_longum", "B. longum")
make_upset_plots(F_strain_data_upset, "Bifidobacterium_breve", "B. breve")
make_upset_plots(F_strain_data_upset, "Bifidobacterium_adolescentis", "B. adolescentis")
make_upset_plots(F_strain_data_upset, "Blautia_massiliensis", "B. massiliensis")
make_upset_plots(F_strain_data_upset, "Faecalibacterium_prausnitzii", "F. prausnitzii")
make_upset_plots(F_strain_data_upset, "Fusobacterium_nucleatum", "F. nucleatum")
make_upset_plots(F_strain_data_upset, "Gardnerella_vaginalis", "G. vaginalis")
make_upset_plots(F_strain_data_upset, "Lacticaseibacillus_rhamnosus", "L. rhamnosus")
make_upset_plots(F_strain_data_upset, "Lactobacillus_crispatus", "L. crispatus")
make_upset_plots(F_strain_data_upset, "Lactobacillus_gasseri", "L. gasseri")
make_upset_plots(F_strain_data_upset, "Lactobacillus_iners", "L. iners")
make_upset_plots(F_strain_data_upset, "Lactobacillus_jensenii", "L. jensenii")
make_upset_plots(F_strain_data_upset, "Ruminococcus_gnavus", "R. gnavus")

species_strains_F3_F2_F1_2 <- species_strains_F3_F2_F1[, c(1:2, ncol(species_strains_F3_F2_F1), 3:(ncol(species_strains_F3_F2_F1) - 1))]
species_strains_F3_F2_F1_species <- species_strains_F3_F2_F1_2[,rownames(species_metaphlan_input)]
species_strains_F3_F2_F1_species_2 <- cbind(species_strains_F3_F2_F1_2[,(1:4)],species_strains_F3_F2_F1_species)
colnames(species_strains_F3_F2_F1_species_2)[4] <- "c_section"


species_strains_F3_F2_F1_species_3 <- species_strains_F3_F2_F1_species_2 %>%
  mutate(sample_type = case_when(
    grepl("F3", specifier) ~ "F_TP3",
    grepl("F2", specifier) ~ "F_TP2",
    grepl("F1", specifier) ~ "F_TP1"
  ),
    c_section = case_when(
    c_section == 0 ~ "vaginal_delivery",
    c_section == 1 ~ "c_section",
    TRUE ~ as.character(c_section)  # Keeps original values for any other cases
  ))
species_strains_F3_F2_F1_species_3 <- species_strains_F3_F2_F1_species_3[, c(1:2, ncol(species_strains_F3_F2_F1_species_3), 3:(ncol(species_strains_F3_F2_F1_species_3) - 1))]

species_strains_F3_F2_F1_species_3 <- species_strains_F3_F2_F1_species_3[order(species_strains_F3_F2_F1_species_3$Studienummer),]
species_strains_F3_F2_F1_species_4 <- species_strains_F3_F2_F1_species_3[,-c(1,4)]

species_names <- c("Bifidobacterium_longum","Lacticaseibacillus_rhamnosus","Lactobacillus_crispatus",
                  "Bifidobacterium_breve","Ruminococcus_gnavus","Bifidobacterium_adolescentis","Blautia_massiliensis",
                  "Faecalibacterium_prausnitzii","Lactobacillus_gasseri","Lactobacillus_iners","Lactobacillus_jensenii",
                  "Fusobacterium_nucleatum","Gardnerella_vaginalis")

# Create a new column for the sum of all other species
species_strains_F3_F2_F1_species_5 <- species_strains_F3_F2_F1_species_4 %>%
  mutate(
    OtherSpeciesSum = rowSums(dplyr::select(., -c(Studienummer, sample_type, c_section, all_of(species_names))))
  ) %>%
  dplyr::select(Studienummer, sample_type, c_section, all_of(species_names), OtherSpeciesSum)




species_strains_F3_F2_F1_species_6 <- merge(mother_infant_submission_date[,c(1,3)],species_strains_F3_F2_F1_species_5,by="Studienummer")

species_strains_F3_F2_F1_long_2 <- species_strains_F3_F2_F1_species_6 %>%
  pivot_longer(
    cols = 5:ncol(species_strains_F3_F2_F1_species_6), # Columns 4 to the end
    names_to = "species", # New column for species names
    values_to = "abundance" # New column for abundance
  )



df_normalized <- species_strains_F3_F2_F1_long_2 %>%
  arrange(difference_in_weeks_faecal) %>% 
  mutate(Studienummer = factor(Studienummer, levels = unique(Studienummer[order(difference_in_weeks_faecal)]))) %>%
  group_by(Studienummer, sample_type)



sample_delivery_type_plot_1 <- ggplot(df_normalized, aes(x = Studienummer, y = abundance, fill = species)) +

 geom_bar(stat = "identity",position="fill" ) +
  facet_grid(sample_type ~ .) +
  theme_minimal() + 
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
  labs(x="samples (earliest to latest submisssion)",y="relative abundance") +
  
  scale_fill_manual(values = group.colors) 

#ggsave(filename = "/PATH/sample_F3_F2_F1_plot_all_species_metaphlan_age_2.pdf", plot = sample_delivery_type_plot_1, width = 12, height = 8, dpi = 600)

df_normalized_2 <- filter(df_normalized, !species %in% c("OtherSpeciesSum"))

sample_delivery_type_plot_1 <- ggplot(df_normalized_2, aes(x = Studienummer, y = abundance, fill = species)) +
 geom_bar(stat = "identity",position="fill") +
  facet_grid(sample_type ~ .) +
  theme_minimal() + 
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
  labs(x="samples (earliest to latest submisssion)",y="relative abundance") +
  
  scale_fill_manual(values = group.colors) 

#ggsave(filename = "/PATH/sample_F3_F2_F1_plot_all_species_metaphlan_age_earliest_to_latest_no_fill.pdf", plot = sample_delivery_type_plot_1, width = 12, height = 8, dpi = 600)



F2_F1_strains_FEAST_b_longum <- filter(F2_F1_strains_FEAST_studie,species %in% c("Bifidobacterium_longum"))
F2_F1_strains_FEAST_b_adolesc <- filter(F2_F1_strains_FEAST_studie,species %in% c("Bifidobacterium_adolescentis"))
F2_F1_strains_FEAST_r_gnavus <- filter(F2_F1_strains_FEAST_studie,species %in% c("Ruminococcus_gnavus"))
F2_F1_strains_FEAST_l_rhamnosus <- filter(F2_F1_strains_FEAST_studie,species %in% c("Lacticaseibacillus_rhamnosus"))
F2_F1_strains_FEAST_bl_massi <- filter(F2_F1_strains_FEAST_studie,species %in% c("Blautia_massiliensis"))
F2_F1_strains_FEAST_f_praus <- filter(F2_F1_strains_FEAST_studie,species %in% c("Faecalibacterium_prausnitzii"))
F2_F1_strains_FEAST_l_crisp <- filter(F2_F1_strains_FEAST_studie,species %in% c("Lactobacillus_crispatus"))


All_strains_combined_F <- F2_F1_all_strains_ordered[,-c(3)]
All_strains_combined_F <- All_strains_combined_F[order(All_strains_combined_F$species,decreasing = TRUE),]

# Pivot the data into a wider format
All_strains_combined_F_wide <- All_strains_combined_F %>%
  pivot_wider(
    names_from = species,     # Columns created from "species"
    values_from = F1  # Values come from the 'V3_F3' column
  )
All_strains_combined_F_wide[is.na(All_strains_combined_F_wide)] <- 0

## order by c_section and breast_fed so I can make the heatmap based on them
All_strains_combined_F_wide <- as.data.frame(All_strains_combined_F_wide)

#All_strains_combined_F_wide <- All_strains_combined_F_wide[order(All_strains_combined_F_wide$delivery_type, All_strains_combined_F_wide$breast_fed),]

rownames(All_strains_combined_F_wide) <- All_strains_combined_F_wide$Studienummer

## now let's make the heatmap  
all_strain_F_matrix <- as.matrix(All_strains_combined_F_wide[,-1])

# Create annotation data frame for `A` and `B`
#annotation <- All_strains_combined_F_wide[, c("breast_fed","delivery_type")]

## define colors
colors <- colorRampPalette(c("purple", "pink","white", "orange", "red"))(1000)

# Define breaks to make the heatmap more sensitive to lower values
breaks <- seq(min(all_strain_F_matrix), max(all_strain_F_matrix), length.out = 1001)

# Plot heatmap
F3_F1_strain_percent_pheatmap <- pheatmap(all_strain_F_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,                                         
      #   annotation_row = annotation,
        show_rownames = FALSE,
        color = colors,
        breaks = breaks,
        main = "Fecal transfer to TP2 from TP1")

#ggsave(filename = "/PATH/F2_from_F1_strain_percent_pheatmap.pdf", plot = F3_F1_strain_percent_pheatmap, width = 12, height = 8, dpi = 600)


## for general, I'm going to separate them into parity and calculate the Bray-Curtis change between the timepoint.
species_strains_V3_V2_V1_md <- merge(metadata_subset_for_vaginal,species_strains_V3_V2_V1,by="Studienummer")
species_strains_V3_V2_V1_p0 <- filter(species_strains_V3_V2_V1_md,Primipara == 0)
species_strains_V3_V2_V1_p0 <- species_strains_V3_V2_V1_p0[,-c(2)]
species_strains_V3_V2_V1_p1 <- filter(species_strains_V3_V2_V1_md,Primipara == 1)
species_strains_V3_V2_V1_p1 <- species_strains_V3_V2_V1_p1[-c(2)]

# only get the strains

## changed below
species_strains_V3_V2_V1_t_names <- grep("\\.t_",names(species_strains_V3_V2_V1_p1),value=TRUE)
species_strains_V3_V2_V1_t <- species_strains_V3_V2_V1_p1[,species_strains_V3_V2_V1_t_names]
species_strains_V3_V2_V1_t <- cbind(species_strains_V3_V2_V1_p1[,c(1:3,ncol(species_strains_V3_V2_V1_p1))],species_strains_V3_V2_V1_t)

#separating the time points for vaginal
species_strains_V3 <- filter(species_strains_V3_V2_V1_t,grepl("V3", specifier))
species_strains_V2 <- filter(species_strains_V3_V2_V1_t,grepl("V2", specifier))
species_strains_V1 <- filter(species_strains_V3_V2_V1_t,grepl("V1", specifier))

# getting specific time points, ordering and removing extra rows V1 V2
species_strains_V1_V2 <- rbind(species_strains_V1,species_strains_V2)
species_strains_V1_V2_ordered <- species_strains_V1_V2[order(species_strains_V1_V2$specifier),]
V1_V2_filtered <- species_strains_V1_V2_ordered %>%
  group_by(Studienummer) %>%
  filter(n_distinct(specifier) > 1) %>%
  ungroup()
species_strains_V1_V2_ordered_2 <- V1_V2_filtered[,-c(1:3)]

# getting specific time points, ordering and removing extra rows V1 V3
species_strains_V1_V3 <- rbind(species_strains_V1,species_strains_V3)
species_strains_V1_V3_ordered <- species_strains_V1_V3[order(species_strains_V1_V3$specifier),]
V1_V3_filtered <- species_strains_V1_V3_ordered %>%
  group_by(Studienummer) %>%
  filter(n_distinct(specifier) > 1) %>%
  ungroup()
species_strains_V1_V3_ordered_2 <- V1_V3_filtered[,-c(1:3)]

# getting specific time points, ordering and removing extra rows V2 V3
species_strains_V2_V3 <- rbind(species_strains_V2,species_strains_V3)
species_strains_V2_V3_ordered <- species_strains_V2_V3[order(species_strains_V2_V3$specifier),]
V2_V3_filtered <- species_strains_V2_V3_ordered %>%
  group_by(Studienummer) %>%
  filter(n_distinct(specifier) > 1) %>%
  ungroup()
species_strains_V2_V3_ordered_2 <- V2_V3_filtered[,-c(1:3)]
   

set.seed(123)
# Remove the sample column
V1_V2_strains_only <- species_strains_V1_V2_ordered_2[,-1]

# Calculate Bray-Curtis dissimilarity for all samples
V1_V2_strains_bray_curtis_dissimilarity <- vegdist(V1_V2_strains_only, method = "bray")
dim(V1_V2_strains_bray_curtis_dissimilarity)
dissim_V1_V2 <- as.vector(as.dist(V1_V2_strains_bray_curtis_dissimilarity))



samples <- unique(sub("_V[0-9]+", "", species_strains_V1_V2_ordered_2$specifier)) # Extract sample IDs without time points
V1_V2_strains_bray_curtis_dissimilarity_matrix <- as.matrix(V1_V2_strains_bray_curtis_dissimilarity)

V1_V2_dissimilarity_values <- c()
sample_names <- c()

for (sample in samples) {
  # Get the indices of the two time points for each sample
  sample_indices <- grep(sample, species_strains_V1_V2_ordered_2$specifier)
  
  # Check if there are exactly 2 time points for the sample
  if (length(sample_indices) == 2) {
    # Extract the pair of dissimilarities between time points
    V1_V2_dissimilarity <- V1_V2_strains_bray_curtis_dissimilarity_matrix[sample_indices[1], sample_indices[2]]
    
    
    V1_V2_dissimilarity_values <- c(V1_V2_dissimilarity_values, V1_V2_dissimilarity)
    sample_names <- c(sample_names, sample)
  } else {
    print(paste("Sample", sample, "does not have exactly 2 time points. Skipping..."))
  }
}

V1_V2_dissimilarity_df <- data.frame(Studiennumer = sample_names,Time_point = rep("V1_V2",times=40), bray_curtis_dissimilarity = V1_V2_dissimilarity_values)



set.seed(123)
# Remove the sample column
V1_V3_strains_only <- species_strains_V1_V3_ordered_2[,-1]

# Calculate Bray-Curtis dissimilarity for all samples
V1_V3_strains_bray_curtis_dissimilarity <- vegdist(V1_V3_strains_only, method = "bray")
dim(V1_V3_strains_bray_curtis_dissimilarity)
dissim_V1_V3 <- as.vector(as.dist(V1_V3_strains_bray_curtis_dissimilarity))

#
samples <- unique(sub("_V[0-9]+", "", species_strains_V1_V3_ordered_2$specifier)) # Extract sample IDs without time points
V1_V3_strains_bray_curtis_dissimilarity_matrix <- as.matrix(V1_V3_strains_bray_curtis_dissimilarity)

V1_V3_dissimilarity_values <- c()
sample_names <- c()

for (sample in samples) {
  # Get the indices of the two time points for each sample
  sample_indices <- grep(sample, species_strains_V1_V3_ordered_2$specifier)
  
  # Check if there are exactly 2 time points for the sample
  if (length(sample_indices) == 2) {
    # Extract the pair of dissimilarities between time points
    V1_V3_dissimilarity <- V1_V3_strains_bray_curtis_dissimilarity_matrix[sample_indices[1], sample_indices[2]]
    
    
    V1_V3_dissimilarity_values <- c(V1_V3_dissimilarity_values, V1_V3_dissimilarity)
    sample_names <- c(sample_names, sample)
  } else {
    print(paste("Sample", sample, "does not have exactly 2 time points. Skipping..."))
  }
}

V1_V3_dissimilarity_df <- data.frame(Studiennumer = sample_names,Time_point = rep("V1_V3",times=40), bray_curtis_dissimilarity = V1_V3_dissimilarity_values)



set.seed(123)
# Remove the sample column
V2_V3_strains_only <- species_strains_V2_V3_ordered_2[,-1]

# Calculate Bray-Curtis dissimilarity for all samples
V2_V3_strains_bray_curtis_dissimilarity <- vegdist(V2_V3_strains_only, method = "bray")
dim(V2_V3_strains_bray_curtis_dissimilarity)
dissim_V2_V3 <- as.vector(as.dist(V2_V3_strains_bray_curtis_dissimilarity))

##

samples <- unique(sub("_V[0-9]+", "", species_strains_V2_V3_ordered_2$specifier)) # Extract sample IDs without time points
V2_V3_strains_bray_curtis_dissimilarity_matrix <- as.matrix(V2_V3_strains_bray_curtis_dissimilarity)

V2_V3_dissimilarity_values <- c()
sample_names <- c()

for (sample in samples) {
  # Get the indices of the two time points for each sample
  sample_indices <- grep(sample, species_strains_V2_V3_ordered_2$specifier)
  
  # Check if there are exactly 2 time points for the sample
  if (length(sample_indices) == 2) {
    # Extract the pair of dissimilarities between time points
    V2_V3_dissimilarity <- V2_V3_strains_bray_curtis_dissimilarity_matrix[sample_indices[1], sample_indices[2]]
    
    
    V2_V3_dissimilarity_values <- c(V2_V3_dissimilarity_values, V2_V3_dissimilarity)
    sample_names <- c(sample_names, sample)
  } else {
    print(paste("Sample", sample, "does not have exactly 2 time points. Skipping..."))
  }
}

V2_V3_dissimilarity_df <- data.frame(Studiennumer = sample_names,Time_point = rep("V2_V3",times=42), bray_curtis_dissimilarity = V2_V3_dissimilarity_values)


species_metaphlan_input <- NULL

for (i in 1:nrow(metaphlan_input)) {

    if (grepl("s__", rownames(metaphlan_input)[i])== TRUE && grepl("t__",rownames(metaphlan_input)[i]) == FALSE) {

        species_metaphlan_input <- rbind(species_metaphlan_input,metaphlan_input[i,])
    } 

}

# V3 V2 V1
species_strains_V3_V2_V1_redone <- species_strains_V3_V2_V1_md[,-c(1,3)] 

species_strains_V3_V2_V1_redone <- species_strains_V3_V2_V1_redone %>%
  mutate(timepoint = str_extract(specifier, "(?<=_).*"))


V_colnames <- colnames(species_strains_V3_V2_V1_redone)[grepl("t__", colnames(species_strains_V3_V2_V1_redone))]
adonis_file_species <- species_strains_V3_V2_V1_redone[, colnames(species_strains_V3_V2_V1_redone) %in% V_colnames]
adonis_file_species_2 <- adonis_file_species + 0.0000001
adonis_file_species_clr <- clr(adonis_file_species_2)
adonis_file_species_clr <- as.data.frame(adonis_file_species_clr)
adonis_stra_timepoint <- adonis2(adonis_file_species ~ timepoint , data = species_strains_V3_V2_V1_redone, na.action = na.omit, permutations = 9999,method="bray")

# V2 V1
species_strains_V2_V1_redone <- filter(species_strains_V3_V2_V1_redone,timepoint %in% c("V1","V2"))

V_colnames <- colnames(species_strains_V2_V1_redone)[grepl("t__", colnames(species_strains_V2_V1_redone))]
adonis_file_species <- species_strains_V2_V1_redone[, colnames(species_strains_V2_V1_redone) %in% V_colnames]

adonis_file_species_2 <- adonis_file_species + 0.0000001
adonis_file_species_clr <- clr(adonis_file_species_2)
adonis_file_species_clr <- as.data.frame(adonis_file_species_clr)
adonis_stra_timepoint <- adonis2(adonis_file_species ~ timepoint , data = species_strains_V2_V1_redone, na.action = na.omit, permutations = 9999,method="bray")

# V3 V1

species_strains_V3_V1_redone <- filter(species_strains_V3_V2_V1_redone,timepoint %in% c("V1","V3"))

V_colnames <- colnames(species_strains_V3_V1_redone)[grepl("t__", colnames(species_strains_V3_V1_redone))]
adonis_file_species <- species_strains_V3_V1_redone[, colnames(species_strains_V3_V1_redone) %in% V_colnames]

adonis_file_species_2 <- adonis_file_species + 0.0000001
adonis_file_species_clr <- clr(adonis_file_species_2)
adonis_file_species_clr <- as.data.frame(adonis_file_species_clr)
adonis_stra_timepoint <- adonis2(adonis_file_species ~ timepoint , data = species_strains_V2_V1_redone, na.action = na.omit, permutations = 9999,method="bray")

# V3 V2

species_strains_V3_V2_redone <- filter(species_strains_V3_V2_V1_redone,timepoint %in% c("V2","V3"))

V_colnames <- colnames(species_strains_V3_V2_redone)[grepl("t__", colnames(species_strains_V3_V2_redone))]
adonis_file_species <- species_strains_V3_V2_redone[, colnames(species_strains_V3_V2_redone) %in% V_colnames]

adonis_file_species_2 <- adonis_file_species + 0.0000001
adonis_file_species_clr <- clr(adonis_file_species_2)
adonis_file_species_clr <- as.data.frame(adonis_file_species_clr)
adonis_stra_timepoint <- adonis2(adonis_file_species ~ timepoint , data = species_strains_V2_V1_redone, na.action = na.omit, permutations = 9999,method="bray")


dist_matrix <- vegdist(adonis_file_species_clr, method = "euclidean")

# Perform PCoA
pcoa_result <- cmdscale(dist_matrix)

# Convert to dataframe
pcoa_df <- data.frame(PC1 = pcoa_result[,1], 
                      PC2 = pcoa_result[,2],
                      Group = species_strains_V2_V1_redone$timepoint)

# Plot
plot <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
stat_ellipse(type = "t") + 
  theme_minimal() +
  labs(title = "PCoA Plot (Bray-Curtis) P1", x = "PC1", y = "PC2")
#ggsave(filename = "/PATH/PCOA_all_vaginal_p1.pdf", plot = plot, width = 12, height = 8, dpi = 600)

plot

dispersion <- betadisper(dist_matrix, species_strains_V2_V1_redone$timepoint)
#pdf("/PATH/vaginal_all_dispersion_P1.pdf")

boxplot(dispersion, main = "Beta Dispersion of Groups P1")

#dev.off()

species_strains_V2_V3_redone <- filter(species_strains_V3_V2_V1_redone,timepoint %in% c("V2","V3"))
metadata_for_aldex <- species_strains_V2_V3_redone[,c("specifier","timepoint")]
rownames(metadata_for_aldex) <- metadata_for_aldex$specifier
metadata_for_aldex <- metadata_for_aldex[,2,drop=FALSE]

grouping <- metadata_for_aldex$timepoint  # e.g., Case vs. Control

# Run ALDEx2
aldex_result_V3_V2 <- aldex(pseudo_counts, 
                      conditions = grouping, 
                      mc.samples = 128,  # Monte Carlo Dirichlet instances
                      test = "t",  # "t" for Welch's t-test, "wilcox" for Wilcoxon
                      effect = TRUE,  # Calculate effect sizes
                      denom = "iqlr")  # Use IQLR normalization

## doing it for the third time point

V3_V2_V1_all_strains_ordered_metadata <- merge(V3_V2_V1_all_strains_ordered,metadata_subset_for_vaginal,
                                            by="Studienummer")

V3_V2_V1_all_strains_ordered_metadata_2 <- V3_V2_V1_all_strains_ordered_metadata %>% 
                                                mutate(V1_V2 = V1+V2,
                                                      V1_present =ifelse(V1>0,"present_in_V1","absent_in_V1"),
                                                      V2_present = ifelse(V2>0,"present_in_V2","absent_in_V2"),
                                                      V_present = ifelse(V1_V2>0,"present_in_V","absent_in_V")) 


V3_V2_V1_all_strains_ordered_metadata_p1 <- filter(V3_V2_V1_all_strains_ordered_metadata_2,Primipara == 1)
V3_V2_V1_all_strains_ordered_metadata_p0 <- filter(V3_V2_V1_all_strains_ordered_metadata_2,Primipara == 0)


All_strains_combined_V <- V3_V2_V1_all_strains_ordered_metadata_2[,-c(2,4,6:10,12)]
All_strains_combined_V <- All_strains_combined_V[order(All_strains_combined_V$species,decreasing = TRUE),]
# Pivot the data into a wider format
All_strains_combined_V_wide <- All_strains_combined_V %>%
  pivot_wider(
    names_from = species,     # Columns created from "species"
    values_from = V1  # Values come from the 'V3_V3' column
  )
All_strains_combined_V_wide[is.na(All_strains_combined_V_wide)] <- 0

## order by c_section and breast_fed so I can make the heatmap based on them
All_strains_combined_V_wide <- as.data.frame(All_strains_combined_V_wide)

All_strains_combined_V_wide <- All_strains_combined_V_wide[order(All_strains_combined_V_wide$Primipara_chr),]


rownames(All_strains_combined_V_wide) <- All_strains_combined_V_wide$Studienummer

## now let's make the heatmap  
all_strain_V_matrix <- as.matrix(All_strains_combined_V_wide[,-c(1,2)])

# Create annotation data frame for `A` and `B`
annotation <- All_strains_combined_V_wide[, c("Primipara_chr"),drop=FALSE]

## define colors
colors <- colorRampPalette(c("purple", "pink","white", "orange", "red"))(1000)

# Define breaks to make the heatmap more sensitive to lower values
breaks <- seq(min(all_strain_V_matrix), max(all_strain_V_matrix), length.out = 1001)

# Plot heatmap
V3_V1_strain_percent_pheatmap <- pheatmap(all_strain_V_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,                                         
         annotation_row = annotation,
        show_rownames = FALSE,
        color = colors,
        breaks = breaks,
        main = "Vaginal transfer to TP3 from TP1")

#ggsave(filename = "/PATH/V3_from_V1_strain_percent_pheatmap.pdf", plot = V3_V1_strain_percent_pheatmap, width = 12, height = 8, dpi = 600)





V3_V2_V1_all_strains_ordered_metadata_metaphlan_2 <- V_strains_wide %>% 
                                                mutate(V1_V2 = V_TP1 + V_TP2,
                                                     V1_present =ifelse(V_TP1>0,"present_in_V1","absent_in_V1"),
                                                      V2_present = ifelse(V_TP2>0,"present_in_V2","absent_in_V2"),
                                                      V3_present = ifelse(V_TP3>0,"present_in_V3","absent_in_V3"),
                                                      V_present = ifelse(V1_V2>0,"present_earlier","absent_earlier"),
                                                      V1_bool =ifelse(V_TP1>0,1,0),
                                                      V2_bool =ifelse(V_TP2>0,1,0),
                                                      V3_bool =ifelse(V_TP3>0,1,0)) 

V3_V2_V1_all_strains_ordered_metadata_metaphlan_2 <- merge(V3_V2_V1_all_strains_ordered_metadata_metaphlan_2,metadata_subset_for_vaginal,
                                            by="Studienummer")
V3_V2_V1_all_strains_ordered_metadata_metaphlan_p1 <- filter(V3_V2_V1_all_strains_ordered_metadata_metaphlan_2,Primipara.x == 1)
V3_V2_V1_all_strains_ordered_metadata_metaphlan_p0 <- filter(V3_V2_V1_all_strains_ordered_metadata_metaphlan_2,Primipara.x == 0)


V3_V2_V1_all_strains_ordered_metadata_2$Primipara_chr <- as.character(V3_V2_V1_all_strains_ordered_metadata_2$Primipara)
V3_V2_V1_all_strains_ordered_metadata_2$parity_presence <- paste(V3_V2_V1_all_strains_ordered_metadata_2$V_present,
                                                                V3_V2_V1_all_strains_ordered_metadata_2$Primipara_chr,
                                                                sep="_")

V2_V1_all_strains_ordered_metadata_2$Primipara_chr <- as.character(V2_V1_all_strains_ordered_metadata_2$Primipara)

V_species_strains_metaphlan_combined <- rbind(reorder_rename_fecal(species_strains_B_longum_V3_V2_V1_ordered),
                                            reorder_rename_fecal(species_strains_B_adolesc_V3_V2_V1_ordered),
                                            reorder_rename_fecal(species_strains_B_breve_V3_V2_V1_ordered),
                                            reorder_rename_fecal(species_strains_B_massilien_V3_V2_V1_ordered),
                                            reorder_rename_fecal(species_strains_F_nucl_V3_V2_V1_ordered), 
                                            reorder_rename_fecal(species_strains_F_prau_V3_V2_V1_ordered),
                                            reorder_rename_fecal(species_strains_L_crispatus_V3_V2_V1_ordered),
                                            reorder_rename_fecal(species_strains_L_gasseri_V3_V2_V1_ordered),
                                            reorder_rename_fecal(species_strains_L_iners_V3_V2_V1_ordered), 
                                            reorder_rename_fecal(species_strains_L_jensenii_V3_V2_V1_ordered),
                                            reorder_rename_fecal(species_strains_G_vaginalis_V3_V2_V1_ordered),
                                            reorder_rename_fecal(species_strains_L_rhamnosus_V3_V2_V1_ordered),
                                            reorder_rename_fecal(species_strains_R_gnavus_V3_V2_V1_ordered))

V_species_strains_metaphlan_combined <- V_species_strains_metaphlan_combined %>%
  mutate(sample_type = factor(sample_type, levels = c("V_TP1", "V_TP2", "V_TP3")))
  
species_strains_V3_V2_V1_parity <- merge(species_strains_V3_V2_V1,metadata_subset_for_vaginal,
                                            by="Studienummer")



species_strains_V3_V2_V1_parity_2 <- species_strains_V3_V2_V1_parity[, c(1:2, ncol(species_strains_V3_V2_V1_parity), 3:(ncol(species_strains_V3_V2_V1_parity) - 1))]

species_strains_V3_V2_V1_parity_2 <- species_strains_V3_V2_V1_parity_2[,-4]


species_strains_V3_V2_V1_2 <- species_strains_V3_V2_V1_parity_2[, c(1:2, ncol(species_strains_V3_V2_V1_parity_2), 3:(ncol(species_strains_V3_V2_V1_parity_2) - 1))]
species_strains_V3_V2_V1_species <- species_strains_V3_V2_V1_2[,rownames(species_metaphlan_input)]
species_strains_V3_V2_V1_species_2 <- cbind(species_strains_V3_V2_V1_2[,(1:4)],species_strains_V3_V2_V1_species)


species_strains_V3_V2_V1_species_3 <- species_strains_V3_V2_V1_species_2 %>%
  mutate(sample_type = case_when(
    grepl("V3", specifier) ~ "V_TP3",
    grepl("V2", specifier) ~ "V_TP2",
    grepl("V1", specifier) ~ "V_TP1"
  ),
    Primipara = case_when(
    Primipara == 0 ~ "Parity_0",
    Primipara == 1 ~ "Parity_1",
    TRUE ~ as.character(Primipara)  # Keeps original values for any other cases
  ))
species_strains_V3_V2_V1_species_3 <- species_strains_V3_V2_V1_species_3[, c(1:2, ncol(species_strains_V3_V2_V1_species_3), 3:(ncol(species_strains_V3_V2_V1_species_3) - 1))]
species_strains_V3_V2_V1_species_3 <- species_strains_V3_V2_V1_species_3[,!names(species_strains_V3_V2_V1_species_3) %in% "Alicyclobacillus_acidiphilus"]

species_strains_V3_V2_V1_species_3 <- species_strains_V3_V2_V1_species_3[order(species_strains_V3_V2_V1_species_3$Studienummer),]
species_strains_V3_V2_V1_species_4 <- species_strains_V3_V2_V1_species_3[,-c(2,4)]

species_names <- c("Bifidobacterium_longum","Lacticaseibacillus_rhamnosus","Lactobacillus_crispatus",
                  "Bifidobacterium_breve","Ruminococcus_gnavus","Bifidobacterium_adolescentis","Blautia_massiliensis",
                  "Faecalibacterium_prausnitzii","Lactobacillus_gasseri","Lactobacillus_iners","Lactobacillus_jensenii",
                  "Fusobacterium_nucleatum","Gardnerella_vaginalis")

# Create a new column for the sum of all other species
species_strains_V3_V2_V1_species_5 <- species_strains_V3_V2_V1_species_4 %>%
  mutate(
    OtherSpeciesSum = rowSums(dplyr::select(., -c(Studienummer, sample_type, Primipara, all_of(species_names))))
  ) %>%
  dplyr::select(Studienummer, sample_type, Primipara, all_of(species_names), OtherSpeciesSum)


species_strains_V3_V2_V1_long <- species_strains_V3_V2_V1_species_5 %>%
  pivot_longer(
    cols = 4:ncol(species_strains_V3_V2_V1_species_5), # Columns 4 to the end
    names_to = "species", # New column for species names
    values_to = "abundance" # New column for abundance
  )

#species_strains_V3_V2_V1_long


species_strains_V3_V2_V1_species_6 <- merge(mother_infant_submission_date[,c(1,4)],species_strains_V3_V2_V1_species_5,by="Studienummer")

species_strains_V3_V2_V1_long_2 <- species_strains_V3_V2_V1_species_6 %>%
  pivot_longer(
    cols = 5:ncol(species_strains_V3_V2_V1_species_6), # Columns 4 to the end
    names_to = "species", # New column for species names
    values_to = "abundance" # New column for abundance
  )

#species_strains_T3_long


df_normalized <- species_strains_V3_V2_V1_long_2 %>%
  arrange(difference_in_weeks_vaginal) %>% 
  mutate(Studienummer = factor(Studienummer, levels = unique(Studienummer[order(difference_in_weeks_vaginal)]))) %>%
  group_by(Studienummer, sample_type)

# Step 2: Create the bar plot, faceted by sample_type (fecal, vaginal, infant)

sample_delivery_type_plot_1 <- ggplot(df_normalized, aes(x = Studienummer, y = abundance, fill = species)) +
# geom_col() +
 geom_bar(stat = "identity",position="fill" ) +
  facet_grid(sample_type ~  Primipara,scales = "free_x") +
  theme_minimal() + 
 # theme(axis.text.x = element_text(size=2,angle = 90, hjust = 1))+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
  labs(x="samples (earliest to latest submission)",y="relative abundance") +
  
  scale_fill_manual(values = group.colors) 

#ggsave(filename = "/PATH/sample_V3_V2_V1_plot_all_species_metaphlan_parity_age_2.pdf", plot = sample_delivery_type_plot_1, width = 12, height = 8, dpi = 600)



species_strains_V3_V2_V1_2 <- species_strains_V3_V2_V1_parity_2[, c(1:2, ncol(species_strains_V3_V2_V1_parity_2), 3:(ncol(species_strains_V3_V2_V1_parity_2) - 1))]
species_strains_V3_V2_V1_species <- species_strains_V3_V2_V1_2[,species_names]
species_strains_V3_V2_V1_species_2 <- cbind(species_strains_V3_V2_V1_2[,(1:4)],species_strains_V3_V2_V1_species)
#colnames(species_strains_V3_V2_V1_species_2)[4] <- "c_section"


species_strains_V3_V2_V1_species_3 <- species_strains_V3_V2_V1_species_2 %>%
  mutate(sample_type = case_when(
    grepl("V3", specifier) ~ "V_TP3",
    grepl("V2", specifier) ~ "V_TP2",
    grepl("V1", specifier) ~ "V_TP1"
  ),
    Primipara = case_when(
    Primipara == 0 ~ "Parity_0",
    Primipara == 1 ~ "Parity_1",
    TRUE ~ as.character(Primipara)  # Keeps original values for any other cases
  ))
species_strains_V3_V2_V1_species_3 <- species_strains_V3_V2_V1_species_3[, c(1:2, ncol(species_strains_V3_V2_V1_species_3), 3:(ncol(species_strains_V3_V2_V1_species_3) - 1))]

species_strains_V3_V2_V1_species_3 <- species_strains_V3_V2_V1_species_3[order(species_strains_V3_V2_V1_species_3$Studienummer),]
species_strains_V3_V2_V1_species_4 <- species_strains_V3_V2_V1_species_3[,-c(2,4)]

species_names <- c("Bifidobacterium_longum","Lacticaseibacillus_rhamnosus","Lactobacillus_crispatus",
                  "Bifidobacterium_breve","Ruminococcus_gnavus","Bifidobacterium_adolescentis","Blautia_massiliensis",
                  "Faecalibacterium_prausnitzii","Lactobacillus_gasseri","Lactobacillus_iners","Lactobacillus_jensenii",
                  "Fusobacterium_nucleatum","Gardnerella_vaginalis")

# Create a new column for the sum of all other species
species_strains_V3_V2_V1_species_5 <- species_strains_V3_V2_V1_species_4 


species_strains_V3_V2_V1_long <- species_strains_V3_V2_V1_species_5 %>%
  pivot_longer(
    cols = 4:ncol(species_strains_V3_V2_V1_species_5), # Columns 4 to the end
    names_to = "species", # New column for species names
    values_to = "abundance" # New column for abundance
  )



df_normalized <- species_strains_V3_V2_V1_long %>%
  group_by(Studienummer, sample_type)
# Step 2: Create the bar plot, faceted by sample_type (fecal, vaginal, infant)

sample_delivery_type_plot_1 <- ggplot(df_normalized, aes(x = Studienummer, y = abundance, fill = species)) +
# geom_col() +
 geom_bar(stat = "identity",position="fill" ) +
  facet_grid(sample_type ~ Primipara,scales = "free_x") +
  theme_minimal() + 
 # theme(axis.text.x = element_text(size=2,angle = 90, hjust = 1))+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
  labs(x="samples",y="relative abundance") +
  
  scale_fill_manual(values = group.colors) 

#ggsave(filename = "/PATH/sample_V3_V2_V1_plot_species_of_interest_metaphlan_parity.pdf", plot = sample_delivery_type_plot_1, width = 12, height = 8, dpi = 600)



V2_V1_strains_FEAST_l_rhamnosus <- filter(V2_V1_strains_FEAST_studie,species %in% c("Lacticaseibacillus_rhamnosus"))
V2_V1_strains_FEAST_l_crisp <- filter(V2_V1_strains_FEAST_studie,species %in% c("Lactobacillus_crispatus"))
V2_V1_strains_FEAST_l_iners <- filter(V2_V1_strains_FEAST_studie,species %in% c("Lactobacillus_iners"))
V2_V1_strains_FEAST_l_gass <- filter(V2_V1_strains_FEAST_studie,species %in% c("Lactobacillus_gasseri"))
V2_V1_strains_FEAST_l_jens <- filter(V2_V1_strains_FEAST_studie,species %in% c("Lactobacillus_jensenii"))


V2_V1_all_strains <- rbind(V2_df_new_g_vag_new_2,l_gass_V2_strains,l_jens_V2_strains,
                          V2_df_new_f_nucl_new_2,l_iners_V2_strains,l_crisp_V2_strains,
                          V2_df_new_b_breve_new_2,V2_df_new_l_rham_new_2,V2_df_new_b_longum_new_2)

V2_V1_all_strains_ordered <- V2_V1_all_strains[order(V2_V1_all_strains$Studienummer, decreasing = FALSE),]


V2_V1_all_strains_ordered_metadata <- merge(V2_V1_all_strains_ordered,metadata_subset_for_vaginal,
                                            by="Studienummer")
V2_V1_all_strains_ordered_metadata_2 <- V2_V1_all_strains_ordered_metadata
V2_V1_all_strains_ordered_metadata_2$Primipara_chr <- as.character(V2_V1_all_strains_ordered_metadata_2$Primipara)



V_species_strains_metaphlan_combined_2 <- merge(V_species_strains_metaphlan_combined, metadata_subset_for_vaginal,
                                             by="Studienummer")

V_strains_wide <- V_species_strains_metaphlan_combined_2 %>%
  pivot_wider(names_from = sample_type, values_from = abundance)


V_strains_wide_presence <- V_strains_wide %>%
    mutate(V_TP1_presence = ifelse(V_TP1>0,1,0),
           V_TP2_presence = ifelse(V_TP2>0,1,0),
           V_TP3_presence = ifelse(V_TP3>0,1,0))



#V_strains_wide_presence

V_strain_data_upset <- V_strains_wide_presence %>%
  mutate(T1_T2 = ifelse(V_TP1_presence == 1 & V_TP2_presence == 1, 1, 0),
         T2_T3 = ifelse(V_TP2_presence == 1 & V_TP3_presence == 1, 1, 0)) 

V_strain_data_upset_p1 <- filter(V_strain_data_upset,Primipara==1)
V_strain_data_upset_p0 <- filter(V_strain_data_upset,Primipara==0)

## Bifidobacterium longum
V_strain_data_upset_l_crispatus_p1 <- filter(V_strain_data_upset_p1,species %in% c("Lactobacillus_crispatus"))



upset(as.data.frame(V_strain_data_upset_l_crispatus_p1), sets = c("V_TP1_presence", "V_TP2_presence","V_TP3_presence")
     ,order.by = "freq",keep.order = TRUE,main.bar.color = "purple",mainbar.y.label = "P1 L. crispatus presence",
      text.scale = c(1.3, 1.3, 1, 1, 2, 2))



V_strain_data_upset_l_crispatus_p0 <- filter(V_strain_data_upset_p0,species %in% c("Lactobacillus_crispatus"))



upset(as.data.frame(V_strain_data_upset_l_crispatus_p0), sets = c("V_TP1_presence", "V_TP2_presence","V_TP3_presence")
     ,order.by = "freq",keep.order = TRUE,main.bar.color = "purple",mainbar.y.label = "P0 L. crispatus presence",
      text.scale = c(1.3, 1.3, 1, 1, 2, 2))



All_strains_combined_V <- V3_V2_V1_all_strains_ordered_metadata_2[,-c(2,4,6:10,12)]
All_strains_combined_V <- All_strains_combined_V[order(All_strains_combined_V$species,decreasing = TRUE),]
# Pivot the data into a wider format
All_strains_combined_V_wide <- All_strains_combined_V %>%
  pivot_wider(
    names_from = species,     
    values_from = V1  
All_strains_combined_V_wide[is.na(All_strains_combined_V_wide)] <- 0

All_strains_combined_V_wide <- as.data.frame(All_strains_combined_V_wide)

All_strains_combined_V_wide <- All_strains_combined_V_wide[order(All_strains_combined_V_wide$Primipara_chr),]


rownames(All_strains_combined_V_wide) <- All_strains_combined_V_wide$Studienummer

## now let's make the heatmap  
all_strain_V_matrix <- as.matrix(All_strains_combined_V_wide[,-c(1,2)])

# Create annotation data frame for `A` and `B`
annotation <- All_strains_combined_V_wide[, c("Primipara_chr"),drop=FALSE]

## define colors
colors <- colorRampPalette(c("purple", "pink","white", "orange", "red"))(1000)

# Define breaks to make the heatmap more sensitive to lower values
breaks <- seq(min(all_strain_V_matrix), max(all_strain_V_matrix), length.out = 1001)

# Plot heatmap
V3_V1_strain_percent_pheatmap <- pheatmap(all_strain_V_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,                                         
         annotation_row = annotation,
        show_rownames = FALSE,
        color = colors,
        breaks = breaks,
        main = "Vaginal transfer to TP3 from TP1")

#ggsave(filename = "/PATH/V3_from_V1_strain_percent_pheatmap.pdf", plot = V3_V1_strain_percent_pheatmap, width = 12, height = 8, dpi = 600)


ll_strains_combined_V <- V2_V1_all_strains_ordered_metadata_2[,-c(3,5,6)]
All_strains_combined_V <- All_strains_combined_V[order(All_strains_combined_V$species,decreasing = TRUE),]

# Pivot the data into a wider format
All_strains_combined_V_wide <- All_strains_combined_V %>%
  pivot_wider(
    names_from = species,     
    values_from = V1  
  )
All_strains_combined_V_wide[is.na(All_strains_combined_V_wide)] <- 0

All_strains_combined_V_wide <- as.data.frame(All_strains_combined_V_wide)

All_strains_combined_V_wide <- All_strains_combined_V_wide[order(All_strains_combined_V_wide$Primipara_chr),]


rownames(All_strains_combined_V_wide) <- All_strains_combined_V_wide$Studienummer

## now let's make the heatmap  
all_strain_V_matrix <- as.matrix(All_strains_combined_V_wide[,-c(1,2)])

# Create annotation data frame for `A` and `B`
annotation <- All_strains_combined_V_wide[, c("Primipara_chr"),drop=FALSE]

## define colors
colors <- colorRampPalette(c("purple", "pink","white", "orange", "red"))(1000)

# Define breaks to make the heatmap more sensitive to lower values
breaks <- seq(min(all_strain_V_matrix), max(all_strain_V_matrix), length.out = 1001)

# Plot heatmap
V2_V1_strain_percent_pheatmap <- pheatmap(all_strain_V_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,                                         
         annotation_row = annotation,
        show_rownames = FALSE,
        color = colors,
        breaks = breaks,
        main = "Vaginal transfer to TP2 from TP1")

#ggsave(filename = "/PATH/V2_from_V1_strain_percent_pheatmap.pdf", plot = V2_V1_strain_percent_pheatmap, width = 12, height = 8, dpi = 600)

