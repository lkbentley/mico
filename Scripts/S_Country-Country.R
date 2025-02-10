########### Calculate connections between countries #############
#### Libraries
library(tidyverse)
library(reshape2)

##how many refs?
sites <- read_csv("Data/site_details_multiple-species V1.csv", 
                  col_types = cols(Ver. = col_character()))
sites%>%select(ReviewID) %>%unique() %>% count()

species <- read_csv("Data/species.csv")
groups <- read_csv("Data/groups_UN.csv")
regions <- c("Northern_America", "Central_America", "Caribbean", "South_America",
             "Southern_Africa", "Middle_Africa", "Western_Africa", "Northern_Africa", "Eastern_Africa",
             "Southern_Europe", "Western_Europe", "Northern_Europe", "Eastern_Europe",
             "Western_Asia", "Central_Asia", "Southern_Asia", "Eastern_Asia", "SE_Asia",
             "Melanesia", "Micronesia", "Polynesia", "Australia_NZ", "Antarctica", "ABNJ")

linkcol <- c(
  Northern_America = "#f7fcb9",
  Central_America = "#78c679",
  Caribbean = "#d9f0a3",
  South_America = "#238443",
  Southern_Africa = "#67000d",
  Middle_Africa = "#a50f15",
  Western_Africa = "#d7301f",
  Northern_Africa = "#ef6548",
  Eastern_Africa = "#fc8d59",
  Southern_Europe = "#7fa8c1",
  Western_Europe = "#557b9e",
  Northern_Europe = "#2f4b7c",
  Eastern_Europe = "#003f5c",
  Western_Asia = "#FFDAB9",
  Central_Asia = "#FFD700",
  Southern_Asia = "#FFA500",
  Eastern_Asia = "#FF8C00",
  SE_Asia = "#FF4500",
  Melanesia = "#54278f",
  Micronesia = "#756bb1",
  Polynesia =  "#9e9ac8",
  Australia_NZ = "#cbc9e2",
  Antarctica = "#702963",
  ABNJ = "grey50"
)


ALLroutesconnections <- read_csv("Data/mico_metasites_crossing_oceans_all.csv")

ALLroutesconnections <- ALLroutesconnections %>% 
  rename(source_eezID = `eez_id...6`,
         source_eezName = `eez_name...7`,
         target_eezID = `eez_id...13`,
         target_eezName = `eez_name...14`) %>%
  mutate(source_eezName = gsub(" ", "_", source_eezName)) %>%
  mutate(target_eezName = gsub(" ", "_", target_eezName)) %>%
  mutate(source_ABNJ = case_when(is.na(source_eezName) & ocean_basin_source != "Unknown" ~ "T")) %>%
  mutate(target_ABNJ = case_when(is.na(target_eezName) & ocean_basin_target != "Unknown" ~ "T")) %>%
  mutate(source_eezName = if_else((is.na(source_eezName) & source_ABNJ=="T"), "ABNJ", source_eezName)) %>% 
  mutate(target_eezName = if_else((is.na(target_eezName) & target_ABNJ=="T"), "ABNJ", target_eezName)) %>% 
  left_join(species, by="common_name") %>%
  left_join(groups, by=join_by("source_eezName" == "value")) %>%
  rename(source_group = group) %>%
  rename(source_country = country) %>%
  left_join(groups, by=join_by("target_eezName" == "value")) %>%
  rename(target_group = group) %>%
  rename(target_country = country) %>%
  filter(!is.na(source_eezName),
         !is.na(target_eezName)) %>%
  mutate(taxon = case_when(
    taxon == "bird" ~ "Seabirds",
    taxon == "fish" ~ "Fishes",
    taxon == "mammal" ~ "Marine Mammals",
    taxon == "seaturtle" ~ "Sea Turtles"
  )) %>%
  mutate(source_group = factor(source_group, levels = names(linkcol))) %>%
  mutate(target_group = factor(target_group, levels = names(linkcol))) %>%
  mutate(source_region = case_when(
    grepl("Africa", source_group) ~ "Africa",
    grepl("America", source_group) ~ "Americas",
    grepl("Caribbean", source_group) ~ "Americas",
    grepl("Asia", source_group) ~ "Asia",
    grepl("Europe", source_group) ~ "Europe",
    grepl("esia", source_group) ~ "Oceania",
    grepl("alia", source_group) ~ "Oceania",
    grepl("Antarctica", source_group) ~ "Antarctica",
    grepl("ABNJ", source_group) ~ "ABNJ")) %>%
  mutate(target_region = case_when(
    grepl("Africa", target_group) ~ "Africa",
    grepl("America", target_group) ~ "Americas",
    grepl("Caribbean", target_group) ~ "Americas",
    grepl("Asia", target_group) ~ "Asia",
    grepl("Europe", target_group) ~ "Europe",
    grepl("esia", target_group) ~ "Oceania",
    grepl("alia", target_group) ~ "Oceania",
    grepl("Antarctica", target_group) ~ "Antarctica",
    grepl("ABNJ", target_group) ~ "ABNJ"))


CountryPairs <- ALLroutesconnections %>%
  select(source_country, target_country) %>%
  distinct()

identify_unique_pairs <- function(data) {
  sorted_data <- t(apply(data, 1, function(x) sort(x)))
  unique_rows <- unique(sorted_data)
  original_order <- apply(unique_rows, 1, function(x) data.frame(col1 = x[1], col2 = x[2]))
  return(original_order)
}

unique_pairs <- identify_unique_pairs(CountryPairs)

# list into a data frame
unique_pairs_df <- do.call(rbind, unique_pairs)

unique_pairs_df <- unique_pairs_df %>%
  filter(!(col1==col2)) 

# the number of connections a country has will be the number of times it appears in the list -1
all_countries <- c(unique_pairs_df$col1, unique_pairs_df$col2)

# Count the occurrences of each country
country_counts <- table(all_countries)

# Convert the counts to a data frame
country_counts_df <- data.frame(letter = names(country_counts), count = as.numeric(country_counts))

# Display the summary
print(country_counts_df)

# average number of connections?
country_counts_df %>%
  summarise(mean_countries=mean(count), sd=sd(count), median_countries=median(count))
