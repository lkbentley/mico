#### Summaries of Connectivity ####

#### Libraries
library(tidyverse)
library(ggbreak)
library(scales)
library(ggOceanMaps)

species <- read_csv("species.csv")
groups <- read_csv("groups_UN.csv")
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
  left_join(groups, by=join_by("target_eezName" == "value")) %>%
  rename(target_group = group) %>%
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


#### Summary Plots ####

### Source/Target EEZs Histogram 

ggplot(ALLroutesconnections, aes(source_group, fill=source_group)) +
  geom_bar(stat="count", position="dodge") +
  facet_wrap(~taxon, nrow = 2) +
  coord_flip(ylim=c(0, 600)) +
  geom_text(stat = "count", aes(label = (ifelse(after_stat(count) > 600, after_stat(count), ""))),
            size = 3,
            position = position_dodge(width = -0.8), vjust = 0, y = 600) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
  scale_fill_manual(values=linkcol)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "off") 

ggplot(ALLroutesconnections, aes(target_group, fill=target_group)) +
  geom_bar(stat="count", position="dodge") +
  facet_wrap(~taxon, nrow = 2) +
  coord_flip(ylim=c(0, 600)) +
  geom_text(stat = "count", aes(label = (ifelse(after_stat(count) > 600, after_stat(count), ""))),
            size = 3,
            position = position_dodge(width = -0.8), vjust = 0, y = 600) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
  scale_fill_manual(values=linkcol)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "off") 



### Number of Individuals DEPARTING and ARRIVING to each region (biases in measurement)

taxon_coloursdark <- c("Seabirds" = "#FF7F0E",
                   "Fishes" = "#1F77B4",
                   "Marine Mammals" = "#D62728",
                   "Sea Turtles" = "#2CA02C")

taxon_colours <- c("Seabirds" = "#FFDDBB",
                       "Fishes" = "#BBDDFF",
                       "Marine Mammals" = "#FFBBBB",
                       "Sea Turtles" = "#BBFFBB")


y_breaks <- seq(0, 1000, by = 100)

ggplot(ALLroutesconnections, aes(source_region, fill=taxon)) +
  geom_bar(stat="count", position=position_dodge(width=-0.8)) +
  #facet_wrap(~taxon, nrow = 2) +
  #scale_y_break(c(1500, 5000)) +
  #scale_y_break(c(5100, 11200)) +
  coord_flip(ylim=c(0,1000)) +
  theme_bw()+
  geom_text(stat = "count", aes(label = (ifelse(after_stat(count) > 1000, after_stat(count), ""))),
            size = 3,
            position = position_dodge(width = -0.8), vjust = 0, y = 1000) +
  labs(x = "Source Region", y = "Individuals (N)", fill = "Taxonomic Group") +
  scale_fill_manual(values=taxon_colours, drop=FALSE) +
  scale_y_continuous(breaks = y_breaks) +
  theme(legend.position = "right", axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5))


ggplot(ALLroutesconnections, aes(target_region, fill=taxon)) +
  geom_bar(stat="count", position=position_dodge(width=-0.8)) +
  #facet_wrap(~taxon, nrow = 2) +
  #scale_y_break(c(1500, 5000)) +
  #scale_y_break(c(5100, 11200)) +
  coord_flip(ylim=c(0,1000)) +
  theme_bw()+
  geom_text(stat = "count", aes(label = (ifelse(after_stat(count) > 1000, after_stat(count), ""))),
            size = 3,
            position = position_dodge(width = -0.8), vjust = 0, y = 1000) +
  labs(x = "Target Region", y = "Individuals (N)", fill = "Taxonomic Group") +
  scale_fill_manual(values=taxon_colours, drop=FALSE) +
  scale_y_continuous(breaks = y_breaks) +
  theme(legend.position = "right", axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5))


### Counting nodes

all_nodes <- ALLroutesconnections %>%
  select(species_code, common_name, source, 
         start_x, start_y, source_radius) %>%
  mutate(source_target = "source") %>%
  rename(metasite = source,
         x = start_x,
         y = start_y,
         radius = source_radius) %>%
  bind_rows(ALLroutesconnections %>% select(species_code, common_name, target, end_x, end_y, target_radius) %>%
            transmute(metasite = target,
                      species_code = species_code,
                      common_name = common_name,
                      x = end_x,
                      y = end_y,
                      radius = target_radius)) %>%
  mutate(source_target = ifelse(is.na(source_target), "target", source_target)) %>% 
  left_join(species) %>%
  mutate(Taxon = case_when(
    taxon == "bird" ~ "Seabirds",
    taxon == "fish" ~ "Fishes",
    taxon == "mammal" ~ "Marine Mammals",
    taxon == "seaturtle" ~ "Sea Turtles")) %>%
  unique()


## all nodes from demog datasheet
site_demographicsRAW <- read_csv("Data/mico_metasitesdemogrFIXED.csv")

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
basemap(all_nodes,
        grid.col = NA,
        land.col = "grey80") +
  geom_point(data = all_nodes, aes(x = longitude, y = latitude, color=Taxon))+
  scale_color_manual(values=taxon_coloursdark) +
  geom_hline(yintercept = 0,linewidth = 0.25) +
  geom_hline(yintercept=23.43, lty=2, linewidth = 0.25) +
  geom_hline(yintercept=-23.43, lty=2, linewidth=0.25) 

## zoom in on caribbean
basemap(all_nodes,limits = c(-100, -45, 40, -10),
        grid.col = NA,
        land.col = "grey80") +
  geom_point(data = all_nodes, aes(x = x, y = y, color=Taxon))+
  scale_color_manual(values=taxon_coloursdark) 


### DENSITY PLOT LATITUDE
ggplot(all_nodes, aes(y=latitude, colour=Taxon, fill=Taxon)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = taxon_coloursdark) +
  scale_fill_manual(values = taxon_colours) +
  theme_classic() +
  geom_hline(yintercept=0, lty=1, linewidth = 0.25) +
  geom_hline(yintercept=23.43, lty=2, linewidth = 0.25) +
  geom_hline(yintercept=-23.43, lty=2, linewidth=0.25) +
  ylab("Latitude (Decimal Degrees)") 

### density plot longitude?
ggplot(all_nodes, aes(x=longitude, colour=Taxon, fill=Taxon)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = taxon_coloursdark) +
  scale_fill_manual(values = taxon_colours) +
  theme_classic() +
  xlab("Longitude (Decimal Degrees)") 

### assembling together (this is broken)
library(patchwork)

plotmap <-basemap(all_nodes,
                  grid.col = NA,
                  land.col = "grey80") +
  geom_point(data = all_nodes, aes(x = longitude, y = latitude, color=Taxon))+
  scale_color_manual(values=taxon_coloursdark) +
  geom_hline(yintercept = 0,linewidth = 0.25) +
  geom_hline(yintercept=23.43, lty=2, linewidth = 0.25) +
  geom_hline(yintercept=-23.43, lty=2, linewidth=0.25) +
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        legend.position = "off")

plotlat <- ggplot(all_nodes, aes(y=latitude, colour=Taxon, fill=Taxon)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = taxon_coloursdark) +
  scale_fill_manual(values = taxon_colours) +
  theme_classic() +
  geom_hline(yintercept=0, lty=1, linewidth = 0.25) +
  geom_hline(yintercept=23.43, lty=2, linewidth = 0.25) +
  geom_hline(yintercept=-23.43, lty=2, linewidth=0.25) +
  ylab("Latitude (Decimal Degrees)") +
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length.y = unit(0, "pt"),
        legend.position = "off")
  
plotlong <- ggplot(all_nodes, aes(x=longitude, colour=Taxon, fill=Taxon)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = taxon_coloursdark) +
  scale_fill_manual(values = taxon_colours) +
  theme_classic() +
  xlab("Longitude (Decimal Degrees)") +
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length.x = unit(0, "pt"),
        legend.position = "off")


plotlong + plot_spacer() + plot_layout(widths = c(4,1))
  
(plotmap + plotlat) + plot_layout(widths = c(6,1))


plot_grid(plotlong,plotmap, plotlat, ncol = 1, align = "hv")





