###### Figures for Map Figure in MiCO Paper #####

#### THIS SCRIPT CREATES FIGURE 3 ####

#### Load Libraries ####
library(readr)
library(tidyverse)
library(ggbreak)
library(scales)
library(ggOceanMaps)
library(ggExtra)
library(rnaturalearth)


#### Extra data and colours ####
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

#### Metasite Data ####

### import data ###
site_demographicsRAW <- read_csv("Data/mico_metasitesdemogrFIXED.csv")

site_demographics <- site_demographicsRAW %>%
  left_join(species, by="common_name") %>%
  mutate(taxon = case_when(
    taxon == "bird" ~ "Seabirds",
    taxon == "fish" ~ "Fishes",
    taxon == "mammal" ~ "Marine Mammals",
    taxon == "seaturtle" ~ "Sea Turtles"))


all_nodes <- site_demographicsRAW %>%
  select(metasite_id, species_code, common_name, 
         latitude, longitude) %>%
  left_join(species) %>%
  mutate(Taxon = case_when(
    taxon == "bird" ~ "Seabirds",
    taxon == "fish" ~ "Fishes",
    taxon == "mammal" ~ "Marine Mammals",
    taxon == "seaturtle" ~ "Sea Turtles")) %>%
  unique()


## how many metasites per taxon
all_nodes %>%
  group_by(Taxon) %>%
  summarise(metasites_per_taxon=n())

## how many species per taxon
all_nodes %>%
  select(Taxon, common_name) %>%
  group_by(Taxon) %>%
  unique() %>%
  summarise(species_per_taxon=n())


### map of metasites globally

taxon_colours <- c("Seabirds" = "#FFDDBB",
                   "Fishes" = "#BBDDFF",
                   "Marine Mammals" = "#FFBBBB",
                   "Sea Turtles" = "#BBFFBB")

taxon_coloursdark <- c("Seabirds" = "#FF7F0E",
                       "Fishes" = "#1F77B4",
                       "Marine Mammals" = "#D62728",
                       "Sea Turtles" = "#2CA02C")

land <- ne_download(scale = "medium", 
                    type = "land", 
                    category = "physical",
                    returnclass = "sf")

globalmetas<- ggplot(land,
        grid.col = NA,
        land.col = "grey80") +
  geom_sf(fill = "grey90", color = "black")+
  geom_point(data = all_nodes, 
             aes(x = longitude, y = latitude, color=Taxon),
             size=.8)+
  scale_color_manual(values=taxon_coloursdark) +
  geom_hline(yintercept = 0,linewidth = 0.25) +
  geom_hline(yintercept=23.43, lty=2, linewidth = 0.25) +
  geom_hline(yintercept=-23.43, lty=2, linewidth=0.25) +
  theme_minimal(base_size = 8) +
  theme(legend.position = "bottom",
    plot.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8))

marg <- ggMarginal(globalmetas, groupColour = TRUE, groupFill = TRUE,
           margins = "both", alpha=0.3,size = 5)

ggsave("Figures/Fig_3.pdf", marg, 
       width = (21.0 - 2.54 * 2), 
       units = "cm", dpi = 300)

###########################################################






