#### Figures for age, sex, life stage, observations and years of data ####

### load packages ###
library(readr)
library(tidyverse)
library(jsonlite)
library(scales)
library(ggrepel)
library(patchwork)

### import data ###
site_demographicsRAW <- read_csv("Data/mico_metasites.csv")
species <- read_csv("Data/species.csv")
groups <- read_csv("Data/groups_UN.csv")
regions <- c("Northern_America", "Central_America", "Caribbean", "South_America",
             "Southern_Africa", "Middle_Africa", "Western_Africa", "Northern_Africa", "Eastern_Africa",
             "Southern_Europe", "Western_Europe", "Northern_Europe", "Eastern_Europe",
             "Western_Asia", "Central_Asia", "Southern_Asia", "Eastern_Asia", "SE_Asia",
             "Melanesia", "Micronesia", "Polynesia", "Australia_NZ", "Antarctica", "ABNJ")

site_demographics <- site_demographicsRAW %>%
  left_join(species, by="common_name") %>%
  mutate(taxon = case_when(
    taxon == "bird" ~ "Seabirds",
    taxon == "fish" ~ "Fishes",
    taxon == "mammal" ~ "Marine Mammals",
    taxon == "seaturtle" ~ "Sea Turtles"))


taxon_coloursdark <- c("Seabirds" = "#FF7F0E",
                       "Fishes" = "#1F77B4",
                       "Marine Mammals" = "#D62728",
                       "Sea Turtles" = "#2CA02C")

sx_colours <-c("Female" = "#FFD54F",
               "Male" = "#7986CB",
               "Unknown" = "#B0BEC5")
  
ls_colours <- c("Breeding Adult" = "#80CBC4",
                "Juvenile" = "#FFCC80",
                "Non-breeding Adult" = "#A5D6A7",
                "Unknown" = "#B0BEC5")

iucn_colours <- c("DD" = "#B0BEC5",
                "LC" = "#80CBC4",
                "NT" = "#AEEE9C",
                "VU" = "#FFD54F",
                "EN" = "#FFB74D",
                "CR" = "#FF8A65")


## number of metasites per taxon
site_demographics %>%
  select(metasite_id, species_code, taxon) %>%
  distinct() %>% 
  group_by(taxon) %>%
  summarise(n())


### make a set of site demographics to correct (those with >100 to start)
site_demographics %>%
  arrange(-n_ndvdl)# %>%
 # write_csv("sites_by_nIDs.csv")





### Plot for SEXES ###
site_demographics$Sex_Male <- NA
site_demographics$Sex_Unknown <- NA
site_demographics$Sex_Female <- NA

# repair the problem where people wrote "UNK" instead of "U"
replace_unk <- function(input_string) {
  return(gsub("UNK", "U", input_string))
}

site_demographics$sex <- replace_unk(site_demographics$sex)


# Convert "sex" JSON column to list and extract counts for "M", "U", and "F"
for (i in seq_along(site_demographics$sex)) {
  tryCatch({
    json_data <- jsonlite::fromJSON(site_demographics$sex[i])
    site_demographics$Sex_Male[i] <- ifelse("M" %in% names(json_data), json_data[["M"]], NA)
    site_demographics$Sex_Unknown[i] <- ifelse("U" %in% names(json_data), json_data[["U"]], NA)
    site_demographics$Sex_Female[i] <- ifelse("F" %in% names(json_data), json_data[["F"]], NA)
  }, error = function(e) {
    # Handle invalid JSON strings here (e.g., print a warning or log the row numbers)
    warning(paste("Invalid JSON string in row", i, "- Setting M, U, F to NA"))
  })
}

site_demographics2<-site_demographics

site_demographics_sex <- site_demographics2 %>%
  select(metasite_id, species_code, common_name, latitude, longitude, radius, n_ndvdl, sex, Sex_Female, Sex_Unknown, Sex_Male) %>%
  mutate(Sex_Female = ifelse(is.na(Sex_Female), 0, Sex_Female)) %>%
  mutate(Sex_Unknown = ifelse(is.na(Sex_Unknown), 0, Sex_Unknown)) %>%
  mutate(Sex_Male = ifelse(is.na(Sex_Male), 0, Sex_Male)) %>%
  pivot_longer(cols = starts_with("Sex_"),
               names_to = "Sex",
               values_to = "N_Indiv") %>%
  mutate(Sex = case_when(Sex == "Sex_Male" ~ "Male",
                         Sex == "Sex_Female" ~ "Female",
                         Sex == "Sex_Unknown" ~ "Unknown")) %>%
  left_join(species, by="common_name") %>%
  mutate(taxon = case_when(
    taxon == "bird" ~ "Seabirds",
    taxon == "fish" ~ "Fishes",
    taxon == "mammal" ~ "Marine Mammals",
    taxon == "seaturtle" ~ "Sea Turtles")) %>%
  filter(!species_code == "HAGR",
         !species_code == "URMA"#,
       #  !species_code == "THCR",
      #   !species_code == "DIEX",
        # !species_code == "MANI",
       #  !species_code == "THOB"
      )



# Calculate the summaries of each section within the fill categories
site_demographics_sexSUM <- site_demographics_sex %>%
  group_by(taxon, Sex # unhash for numbers 
           ) %>%
  summarise(N_Indiv = sum(N_Indiv))

## plot proportions
ggplot(site_demographics_sex, aes(x=taxon, y = N_Indiv, fill = Sex)) +
  geom_bar(position = "fill", stat = "identity") +
 # geom_text(data = site_demographics_sexSUM, aes(label = scales::comma(N_Indiv))) +
  scale_fill_manual(values = sx_colours) +
  scale_y_continuous(name="Proportion of Individuals") +
  scale_x_discrete("Taxonomic Group") +
  theme_bw()

## plot absolute values ### USE THIS FOR PAPER 
ggplot(site_demographics_sex, aes(x=taxon, y = N_Indiv, fill = Sex)) +
  geom_col(position = "stack") +
 # geom_text_repel(data = site_demographics_sexSUM, aes(label = scales::comma(N_Indiv)),direction = "y", position = position_stack(vjust=0.5)) +
  scale_y_continuous(label=comma, name = "Number of Individuals") +
  scale_x_discrete("Taxonomic Group") +
  scale_fill_manual(values = sx_colours) +
  theme_bw() +
  theme(legend.position = c(0.85,0.85)) 
  

#### plot for life stage ####

codes_life_stage <- read_csv("Data/codes_life_stage.csv")

site_demographics <- site_demographicsRAW %>%
  left_join(species, by="common_name") %>%
  mutate(taxon = case_when(
    taxon == "bird" ~ "Seabirds",
    taxon == "fish" ~ "Fishes",
    taxon == "mammal" ~ "Marine Mammals",
    taxon == "seaturtle" ~ "Sea Turtles")) %>%
  filter(metasite_id != "METASITE_HAGR_AJ") # GREY SEALS

site_demographics$life_stage <- replace_unk(site_demographics$life_stage)

site_demographics$A <- NA
site_demographics$AB <- NA
site_demographics$ABG <- NA
site_demographics$`ABG-BG` <- NA
site_demographics$`ABG-CB` <- NA
site_demographics$`ABG-CR` <- NA
site_demographics$ACR <- NA
site_demographics$AI <- NA
site_demographics$AIN <- NA
site_demographics$APG <- NA
site_demographics$AS <- NA
site_demographics$`AS-PB` <- NA
site_demographics$`AS-S` <- NA
site_demographics$C <- NA
site_demographics$CA <- NA
site_demographics$COY <- NA
site_demographics$I <- NA
site_demographics$J <- NA
site_demographics$NB <- NA
site_demographics$NH <- NA
site_demographics$SUB <- NA
site_demographics$UNK <- NA
site_demographics$YL <- NA

# Convert "life_stage" JSON column to list and extract counts for values
for (i in seq_along(site_demographics$life_stage)) {
  tryCatch({
    json_data <- jsonlite::fromJSON(site_demographics$life_stage[i])
    site_demographics$A[i] <- ifelse("A" %in% names(json_data), json_data[["A"]], NA)
    site_demographics$AB[i] <- ifelse("AB" %in% names(json_data), json_data[["AB"]], NA)
    site_demographics$ABG[i] <- ifelse("ABG" %in% names(json_data), json_data[["ABG"]], NA)
    site_demographics$`ABG-BG`[i] <- ifelse("ABG-BG" %in% names(json_data), json_data[["ABG-BG"]], NA)
    site_demographics$`ABG-CB`[i] <- ifelse("ABG-CB" %in% names(json_data), json_data[["ABG-CB"]], NA)
    site_demographics$`ABG-CR`[i] <- ifelse("ABG-CR" %in% names(json_data), json_data[["ABG-CR"]], NA)
    site_demographics$ACR[i] <- ifelse("ACR" %in% names(json_data), json_data[["ACR"]], NA)
    site_demographics$AI[i] <- ifelse("AI" %in% names(json_data), json_data[["AI"]], NA)
    site_demographics$AIN[i] <- ifelse("AIN" %in% names(json_data), json_data[["AIN"]], NA)
    site_demographics$APG[i] <- ifelse("APG" %in% names(json_data), json_data[["APG"]], NA)
    site_demographics$AS[i] <- ifelse("AS" %in% names(json_data), json_data[["AS"]], NA)
    site_demographics$`AS-PB`[i] <- ifelse("AS-PB" %in% names(json_data), json_data[["AS-PB"]], NA)
    site_demographics$`AS-S`[i] <- ifelse("AS-S" %in% names(json_data), json_data[["AS-S"]], NA)
    site_demographics$C[i] <- ifelse("C" %in% names(json_data), json_data[["C"]], NA)
    site_demographics$CA[i] <- ifelse("CA" %in% names(json_data), json_data[["CA"]], NA)
    site_demographics$COY[i] <- ifelse("COY" %in% names(json_data), json_data[["COY"]], NA)
    site_demographics$I[i] <- ifelse("I" %in% names(json_data), json_data[["I"]], NA)
    site_demographics$J[i] <- ifelse("J" %in% names(json_data), json_data[["J"]], NA)
    site_demographics$NB[i] <- ifelse("NB" %in% names(json_data), json_data[["NB"]], NA)
    site_demographics$NH[i] <- ifelse("NH" %in% names(json_data), json_data[["NH"]], NA)
    site_demographics$SUB[i] <- ifelse("SUB" %in% names(json_data), json_data[["SUB"]], NA)
    site_demographics$UNK[i] <- ifelse("U" %in% names(json_data), json_data[["U"]], NA)
    site_demographics$YL[i] <- ifelse("YL" %in% names(json_data), json_data[["YL"]], NA)
  }, error = function(e) {
    # Handle invalid JSON strings here (e.g., print a warning or log the row numbers)
    warning(paste("Invalid JSON string in row", i, "- Setting col to NA"))
  })
}

site_demographics_lifestage <- site_demographics %>%
  select(metasite_id, species_code, common_name, latitude, longitude, radius, n_ndvdl, taxon, life_stage, c(14:36)) %>%
  mutate_at(vars(10:32), ~replace_na(., 0)) %>%
  pivot_longer(cols = (10:32),
               names_to = "Life_Stage",
               values_to = "N_Indiv") %>%
  filter(!(N_Indiv == 0)) %>%
  left_join(codes_life_stage, by=c("Life_Stage" = "code")) %>%
  filter(!species_code == "HAGR",
         !species_code == "URMA"#,
        # !species_code == "THCR",
        # !species_code == "DIEX",
       #  !species_code == "MANI",
        # !species_code == "THOB"
       )

  
# Calculate the summaries of each section within the fill categories
site_demographics_lifestageSUM <- site_demographics_lifestage %>%
  group_by(taxon, LS_Simple) %>%
  summarise(N_Indiv = sum(N_Indiv))

## plot proportions
ggplot(site_demographics_lifestage, aes(x=taxon, y = N_Indiv, fill = LS_Simple)) +
  geom_bar(position = "fill", stat = "identity") +  
  #geom_text(data = site_demographics_lifestageSUM, aes(label = scales::comma(N_Indiv)), position = position_fill(vjust = 0.5)) +
geom_text(aes(label = scales::comma(after_stat(y)), group = taxon), 
    stat = 'summary', fun = sum, position = position_fill(vjust = -0.02)) +
  scale_fill_manual(values = ls_colours, name = "Life Stage") +
  scale_y_continuous(label=comma, name = "Proportion of Individuals") +
  scale_x_discrete("Taxonomic Group") +
  theme_bw()

## plot absolute values
ggplot(site_demographics_lifestage, aes(x=taxon, y = N_Indiv, fill = LS_Simple)) +
  geom_col(position = "stack") +  
  #geom_text_repel(data = site_demographics_lifestageSUM, aes(label = scales::comma(N_Indiv)), direction="y", position = position_stack(vjust = 0.5)) +
  scale_y_continuous(label=comma, name = "Number of Individuals") +
  scale_x_discrete("Taxonomic Group") +
  scale_fill_manual(values = ls_colours, name = "Life Stage") +
  theme_bw()  +
  theme(legend.position = c(0.85,0.85)) 

## plot proportions (number of individuals)
ggplot(filter(site_demographics_iucn), aes(x=taxon, y = n_ndvdl, fill = IUCN)) +
  geom_bar(position = "fill", stat = "identity") +  
  geom_text(aes(label = scales::comma(after_stat(y)), group = taxon), 
            stat = 'summary', fun = sum, position = position_fill(vjust = -0.02)) +
  scale_fill_manual(values = iucn_colours, name = "IUCN Red List Category") +
  scale_y_continuous(label=comma, name = "Proportion of Individuals") +
  scale_x_discrete("Taxonomic Group") +
  theme_bw()


## plot proportions of SPECIES
ggplot(filter(unique(select(site_demographics_iucn, common_name, taxon, IUCN))), aes(x=taxon, fill = IUCN)) +
  geom_bar(position = "fill") +  
  scale_fill_manual(values = iucn_colours, name = "IUCN Red List Category") +
  scale_y_continuous(label=comma, name = "Proportion of Species") +
  scale_x_discrete("Taxonomic Group") +
  theme_bw()

## plot NUMBERS OF SPECIES ## I think this is the most interesting one?
ggplot(filter(unique(select(site_demographics_iucn, common_name, taxon, IUCN))), aes(x=taxon, fill = IUCN)) +
  geom_bar(position = "stack") +  
  scale_fill_manual(values = iucn_colours, name = "Red List Category") +
  scale_y_continuous(label=comma, name = "Number of Species") +
  scale_x_discrete("Taxonomic Group") +
  theme_bw(base_size = 14) +
  theme(legend.position = c(0.6,0.7)) 



## plots for paper

paper_plotS <- ggplot(site_demographics_sex, aes(x=taxon, y = N_Indiv, fill = Sex)) +
  geom_col(position = "stack") +
  # geom_text_repel(data = site_demographics_sexSUM, aes(label = scales::comma(N_Indiv)),direction = "y", position = position_stack(vjust=0.5)) +
  scale_y_continuous(label=comma, name = "Number of Individuals") +
  scale_x_discrete("Taxonomic Group") +
  scale_fill_manual(values = sx_colours) +
  theme_bw(base_size = 14) +
  theme(legend.position = c(0.75,0.75)) 

paper_plotLS <- ggplot(site_demographics_lifestage, aes(x=taxon, y = N_Indiv, fill = LS_Simple)) +
  geom_col(position = "stack") +  
  #geom_text_repel(data = site_demographics_lifestageSUM, aes(label = scales::comma(N_Indiv)), direction="y", position = position_stack(vjust = 0.5)) +
  scale_y_continuous(label=comma, name = "Number of Individuals") +
  scale_x_discrete("Taxonomic Group") +
  scale_fill_manual(values = ls_colours, name = "Life Stage") +
  theme_bw(base_size = 14) +
  theme(legend.position = c(0.75,0.75)) 


### patchworking

paper_plotS + paper_plotLS + plot_annotation(tag_levels = "A")

### supp table of species in these plots
#site_demographics_lifestage %>%
#  select(taxon, species_code, common_name) %>%
#  arrange(taxon) %>%
#  distinct() %>%
#  write_csv("Table_S1.csv")

## select only telemetry sites from master file, count tags

telem_sites <- site_demographicsRAW %>%
  filter(Method=="T") %>%
  mutate(NumIDs = as.double(NumIDs)) %>%
  select(-`Ver.`) %>%
  left_join(focus_n10_fixed) %>%
  mutate(plot_NumID = ifelse(is.na(NumIDsTelem), NumIDs, NumIDsTelem)) %>%
  filter(!Species == "HAGR",
         !Species == "URMA")

tags_deployed_telem_sites <- telem_sites %>%
  filter(SiteID == "A") %>%
  left_join(species, by= c("CommonName"="common_name")) %>%
  mutate(taxon = case_when(
    taxon == "bird" ~ "Seabirds",
    taxon == "fish" ~ "Fishes",
    taxon == "mammal" ~ "Marine Mammals",
    taxon == "seaturtle" ~ "Sea Turtles")) %>% 
  dplyr::select(Species, taxon, CommonName, ReviewID, Method, SiteID, Location, Ocean, Long, Lat, Years,
                Radius, Sampling, plot_NumID, NumIDs, Numroutes, Numconns, Activities, UniqueIDbySPP) 

# number of tags deployed across species 
tags_deployed_telem_sites %>%
  filter(!is.na(plot_NumID)) %>% 
  filter(!str_detect(Sampling, "radio")) %>% # remove radio tags
  filter(!str_detect(ReviewID, "S3BMTSF9")) %>%# drop the tancell ref as it duplicates tags 
  group_by(Species, Location) %>%
  slice_max(order_by = plot_NumID) %>%
  group_by(taxon) %>% 
  summarise(count=sum(plot_NumID))
