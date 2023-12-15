###### Data Exploration and Preliminary Visualisations ######

### Load Libraries
library(readr)
library(tidyverse)
library(circlize)


### Import datasheet

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
  Southern_Europe = "#003f5c",
  Western_Europe = "#2f4b7c",
  Northern_Europe = "#557b9e",
  Eastern_Europe = "#7fa8c1",
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

linkcolS <- c(
  "Northern America" = "#f7fcb9",
  "Central America" = "#78c679",
  "Caribbean" = "#d9f0a3",
  "South America" = "#238443",
  "Southern Africa" = "#67000d",
  "Middle Africa" = "#a50f15",
  "Western Africa" = "#d7301f",
  "Northern Africa" = "#ef6548",
  "Eastern Africa" = "#fc8d59",
  "Southern Europe" = "#003f5c",
  "Western Europe" = "#2f4b7c",
  "Northern Europe" = "#557b9e",
  "Eastern Europe" = "#7fa8c1",
  "Western Asia" = "#FFDAB9",
  "Central Asia" = "#FFD700",
  "Southern Asia" = "#FFA500",
  "Eastern Asia" = "#FF8C00",
  "SE Asia" = "#FF4500",
  "Melanesia" = "#54278f",
  "Micronesia" = "#756bb1",
  "Polynesia" =  "#9e9ac8",
  "Australia NZ" = "#cbc9e2",
  "Antarctica" = "#702963",
  "ABNJ" = "grey50"
)

routesconnections <- read_csv("Data/mico_metasites_crossing_oceans_all.csv")

routesconnections <- routesconnections %>% 
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
  filter(!is.na(source_eezName),
         !is.na(target_eezName)) %>%
  left_join(groups, by=join_by(source_eezName == value)) %>%
  rename(source_eezGroup = group) %>%
  left_join(groups, by=join_by(target_eezName == value)) %>%
  rename(target_eezGroup = group)


#### identify REGIONAL CONNECTIONS for each metasite. (reduce duplication by country)

routesconnectionsMINI <- routesconnections %>%
  distinct(species_code, common_name, scientific_name, type, 
           source, start_x, start_y, source_radius, source_eezGroup, 
           target, end_x, end_y, target_radius, target_eezGroup,
           taxon) 

# split into taxonomic groups
seaturtleRC <- routesconnectionsMINI %>% filter(taxon=="seaturtle")
mammalRC <- routesconnectionsMINI %>% filter(taxon=="mammal")
fishRC <- routesconnectionsMINI %>% filter(taxon=="fish")
birdRC <- routesconnectionsMINI %>% filter(taxon=="bird")


sourceRCs <- list(seaturtleRC=seaturtleRC, 
                  mammalRC = mammalRC, 
                  fishRC = fishRC, 
                  birdRC = birdRC)

# Set the radius for the labels
radius <- 1.2

# Set the starting angle for the labels
start_angle <- 0

image_resolution <- 300


for (i in 1:length(sourceRCs)) {
  dataset <- sourceRCs[[i]]
  name <- names(sourceRCs)[i]
  
  dataset <- dataset %>%
    mutate(target_eezGroup = gsub("_", " ", target_eezGroup),
           source_eezGroup = gsub("_", " ", source_eezGroup))
  
  # pull nodes out
  nodes <- unique(c(dataset$source_eezGroup, dataset$target_eezGroup))
  
  # empty matrix
  adjacency_matrix <- matrix(0, nrow = length(nodes), ncol = length(nodes), dimnames = list(nodes, nodes))
  
  # Iterate over the pair data and update the adjacency matrix
  for (i in 1:nrow(dataset)) {
    node1 <- dataset$source_eezGroup[i]
    node2 <- dataset$target_eezGroup[i]
    adjacency_matrix[node1, node2] <- adjacency_matrix[node1, node2] + 1
  }
  
  # make list of regions in this subset
  subset_regions <- nodes[nodes %in% gsub("_", " ", regions)]
  subset_regions <- subset_regions[match(gsub("_", " ", regions), subset_regions)]
  subset_regions <- na.omit(subset_regions)
  
  linkcol2 <- linkcolS[names(linkcolS) %in% subset_regions]
  
  # Save the chord diagram as an image with high resolution
  png(filename = paste0("ChordDiagrams/Paper/", name, "_region_CD_mini.png"), width = 5, height = 5, units = "in", res = 300)
  plot.new()
  par(cex = 0.45, mar = c(4, 4, 4,4))
  circos.clear()
  
  chordDiagram(adjacency_matrix, 
               order = subset_regions,
               grid.col = linkcol2,
               transparency = 0.3,
               annotationTrack = "grid", #group = group, 
               #directional = 1,
               big.gap = 3, small.gap = 1,
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(adjacency_matrix))))))
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    xplot = get.cell.meta.data("xplot")
    sector.name = get.cell.meta.data("sector.index")
    if(abs(xplot[2] - xplot[1]) < 12) {
      circos.text(CELL_META$xcenter, ylim[1] + cm_h(2.5), 
                  sector.name, facing = "clockwise",
                  niceFacing = TRUE, adj = c(0, 0.5), cex = 1.9, font = 1)
    } else {
      circos.text(CELL_META$xcenter, ylim[1] + cm_h(2.5), 
                  sector.name, facing = "inside",
                  niceFacing = TRUE, adj = c(0.5, 0), cex = 1.9, font = 1)
    }
    circos.axis(h = "bottom",
                labels.cex = .8,
                sector.index = sector.name
    )
  }, bg.border = NA)
  dev.off()
  
  cat("Chord diagram for", name, "created and saved as", paste0(name, "_region_CD_mini.png"), "\n\n")
}













######################################################################################### unused

## chord diagrams for all countries within 4 groups
for (i in 1:length(sourceRCs)) {
  dataset <- sourceRCs[[i]]
  name <- names(sourceRCs)[i]
  
  # pull nodes out
  nodes <- unique(c(dataset$source_eezName, dataset$target_eezName))
  
  # empty matrix
  adjacency_matrix <- matrix(0, nrow = length(nodes), ncol = length(nodes), dimnames = list(nodes, nodes))
  
  # Iterate over the pair data and update the adjacency matrix
    for (i in 1:nrow(dataset)) {
      node1 <- dataset$source_eezName[i]
      node2 <- dataset$target_eezName[i]
      adjacency_matrix[node1, node2] <- adjacency_matrix[node1, node2] + 1
    }
  
  # Find nodes with no connections
  no_connections <- which(rowSums(adjacency_matrix) == 0)
  
  # Remove nodes with no connections from the adjacency matrix
  adjacency_matrix <- adjacency_matrix[-no_connections, -no_connections]
  
  # create grouping variable
  nm <- unique(unlist(dimnames(adjacency_matrix)))
  groups.sub <- as.data.frame(nm) %>% left_join(groups, by=join_by("nm"== "value"))
  group <- structure(groups.sub$group, names = nm)
  
  # Save the chord diagram as an image with high resolution
  png(filename = paste0("ChordDiagrams/", name, "_chord_diagram.png"), width = 10, height = 10, units = "in", res = image_resolution)
  plot.new()
  par(cex = 0.5, mar = c(0, 0, 0, 0))
  circos.clear()

  chordDiagram(adjacency_matrix, annotationTrack = "grid", group = group, big.gap = 1, small.gap = 0.5,
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(adjacency_matrix))))))
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, 
                ylim[1] + cm_h(2), 
                sector.name, 
                facing = "clockwise",
                niceFacing = TRUE, 
                adj = c(0, 0.5),
                cex = 1,
                font = 2)
    circos.axis(h = "bottom",
                labels.cex = .6,
                sector.index = sector.name
    )
  }, bg.border = NA)
  dev.off()
  
  cat("Chord diagram for", name, "created and saved as", paste0(name, "_chord_diagram.png"), "\n\n")
}


###### chord diagrams for regions ######

for (i in 1:length(sourceRCs)) {
  dataset <- sourceRCs[[i]]
  name <- names(sourceRCs)[i]
  
  dataset <- dataset %>%
    left_join(groups, by=join_by(source_eezName == value)) %>%
    rename(source_eezGroup = group) %>%
    left_join(groups, by=join_by(target_eezName == value)) %>%
    rename(target_eezGroup = group)
  
  # pull nodes out
  nodes <- unique(c(dataset$source_eezGroup, dataset$target_eezGroup))
  
  # empty matrix
  adjacency_matrix <- matrix(0, nrow = length(nodes), ncol = length(nodes), dimnames = list(nodes, nodes))
  
  # Iterate over the pair data and update the adjacency matrix
  for (i in 1:nrow(dataset)) {
    node1 <- dataset$source_eezGroup[i]
    node2 <- dataset$target_eezGroup[i]
    adjacency_matrix[node1, node2] <- adjacency_matrix[node1, node2] + 1
  }
  
  # make list of regions in this subset
  subset_regions <- nodes[nodes %in% regions]
  subset_regions <- subset_regions[match(regions, subset_regions)]
  subset_regions <- na.omit(subset_regions)
  
  linkcol2 <- linkcol[names(linkcol) %in% subset_regions]
  
  # Find nodes with no connections
 # no_connections <- which(rowSums(adjacency_matrix) == 0)
  
  # Remove nodes with no connections from the adjacency matrix
  #adjacency_matrix <- adjacency_matrix[-no_connections, -no_connections]
  
  # create grouping variable
  #nm <- unique(unlist(dimnames(adjacency_matrix)))
 # groups.sub <- as.data.frame(nm) %>% left_join(groups, by=join_by("nm"== "value"))
 # group <- structure(groups.sub$group, names = nm)
  
  # Save the chord diagram as an image with high resolution
  png(filename = paste0("ChordDiagrams/Group_Level/", name, "_region_chord_diagram.png"), width = 10, height = 10, units = "in", res = image_resolution)
  plot.new()
  par(cex = 0.5, mar = c(0, 0, 0, 0))
  circos.clear()
  
  chordDiagram(adjacency_matrix, 
               order = subset_regions,
               grid.col = linkcol2,
               transparency = 0.3,
               annotationTrack = "grid", #group = group, big.gap = 1, small.gap = 0.5,
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(adjacency_matrix))))))
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, 
                ylim[1] + cm_h(2), 
                sector.name, 
                facing = "clockwise",
                niceFacing = TRUE, 
                adj = c(0, 0.5),
                cex = 1.3,
                font = 2)
    circos.axis(h = "bottom",
                labels.cex = .8,
                sector.index = sector.name
    )
  }, bg.border = NA)
  dev.off()
  
  cat("Chord diagram for", name, "created and saved as", paste0(name, "region_chord_diagram.png"), "\n\n")
}

##### BY REGION BY SPECIES FOR TOP 50 MOST COMMON SPECIES #####
top50 <- routesconnections %>%
       group_by(common_name) %>%
       summarise(n=n()) %>%
       arrange(desc(n)) %>%
       slice(1:50) %>%
       select(common_name)

for (i in 1:nrow(top50)) {
  common_nameT <- as.vector(top50[i,1])
  dataset <- routesconnections %>%
    filter(common_name == common_nameT)
  
  dataset <- dataset %>%
    left_join(groups, by=join_by(source_eezName == value)) %>%
    rename(source_eezGroup = group) %>%
    left_join(groups, by=join_by(target_eezName == value)) %>%
    rename(target_eezGroup = group)
  
  # pull nodes out
  nodes <- unique(c(dataset$source_eezGroup, dataset$target_eezGroup))
  
  # empty matrix
  adjacency_matrix <- matrix(0, nrow = length(nodes), ncol = length(nodes), dimnames = list(nodes, nodes))
  
  # Iterate over the pair data and update the adjacency matrix
  for (i in 1:nrow(dataset)) {
    node1 <- dataset$source_eezGroup[i]
    node2 <- dataset$target_eezGroup[i]
    adjacency_matrix[node1, node2] <- adjacency_matrix[node1, node2] + 1
  }
  
  # make list of regions in this subset
  subset_regions <- nodes[nodes %in% regions]
  subset_regions <- subset_regions[match(regions, subset_regions)]
  subset_regions <- na.omit(subset_regions)
  
  linkcol2 <- linkcol[names(linkcol) %in% subset_regions]
  
  # Save the chord diagram as an image with high resolution
  png(filename = paste0("ChordDiagrams/Indiv_Species/", common_nameT, "_chord_diagram.png"), width = 10, height = 10, units = "in", res = image_resolution)
  plot.new()
  par(cex = 0.5, mar = c(0, 0, 0, 0))
  circos.clear()
  
  chordDiagram(adjacency_matrix, 
               order = subset_regions,
               grid.col = linkcol2,
               transparency = 0.3,
               annotationTrack = "grid", #group = group, big.gap = 1, small.gap = 0.5,
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(adjacency_matrix))))))
  title(paste0(common_nameT))
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, 
                ylim[1] + cm_h(2), 
                sector.name, 
                facing = "clockwise",
                niceFacing = TRUE, 
                adj = c(0, 0.5),
                cex = 1.3,
                font = 2)
    circos.axis(h = "bottom",
                labels.cex = .8,
                sector.index = sector.name
    )
  }, bg.border = NA)
  dev.off()
}



##### BY SPECIES BY COUNTRY #####
for (i in 1:nrow(top50)) {
  common_nameT <- as.vector(top50[i,1])
  dataset <- routesconnections %>%
    filter(common_name == common_nameT)
  
  dataset <- dataset %>%
    left_join(groups, by=join_by(source_eezName == value)) %>%
    rename(source_eezGroup = group) %>%
    left_join(groups, by=join_by(target_eezName == value)) %>%
    rename(target_eezGroup = group)
  
  # pull nodes out
  nodes <- unique(c(dataset$source_eezName, dataset$target_eezName))
  
  # empty matrix
  adjacency_matrix <- matrix(0, nrow = length(nodes), ncol = length(nodes), dimnames = list(nodes, nodes))
  
  # Iterate over the pair data and update the adjacency matrix
  for (i in 1:nrow(dataset)) {
    node1 <- dataset$source_eezName[i]
    node2 <- dataset$target_eezName[i]
    adjacency_matrix[node1, node2] <- adjacency_matrix[node1, node2] + 1
  }
  
  # create grouping variable
  nm <- unique(unlist(dimnames(adjacency_matrix)))
  groups.sub <- as.data.frame(nm) %>% left_join(groups, by=join_by("nm"== "value"))
  group <- structure(groups.sub$group, names = nm)
  
  # Save the chord diagram as an image with high resolution
  png(filename = paste0("ChordDiagrams/Indiv_Species_Countries/", common_nameT, "_chord_diagram.png"), width = 10, height = 10, units = "in", res = image_resolution)
  plot.new()
  par(cex = 0.5, mar = c(0, 0, 0, 0))
  circos.clear()
  
  chordDiagram(adjacency_matrix, 
               transparency = 0.5,
               annotationTrack = "grid", group = group, big.gap = 2, small.gap = 0.8,
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(adjacency_matrix))))))
  title(paste0(common_nameT))
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, 
                ylim[1] + cm_h(2), 
                sector.name, 
                facing = "clockwise",
                niceFacing = TRUE, 
                adj = c(0, 0.5),
                cex = 1.3,
                font = 2)
    circos.axis(h = "bottom",
                labels.cex = .8,
                sector.index = sector.name
    )
  }, bg.border = NA)
  dev.off()
}





##### This doesn't work -- just reproduces the originals which are not scaled. 
##### Working in this section to adjust the labelling etc. 
### Comparing two chord diagrams (same scale across both -- relative values) 
linkcolS <- c(
  "Northern America" = "#f7fcb9",
  "Central America" = "#78c679",
  "Caribbean" = "#d9f0a3",
  "South America" = "#238443",
  "Southern Africa" = "#67000d",
  "Middle Africa" = "#a50f15",
  "Western Africa" = "#d7301f",
  "Northern Africa" = "#ef6548",
  "Eastern Africa" = "#fc8d59",
  "Southern Europe" = "#003f5c",
  "Western Europe" = "#2f4b7c",
  "Northern Europe" = "#557b9e",
  "Eastern Europe" = "#7fa8c1",
  "Western Asia" = "#FFDAB9",
  "Central Asia" = "#FFD700",
  "Southern Asia" = "#FFA500",
  "Eastern Asia" = "#FF8C00",
  "SE Asia" = "#FF4500",
  "Melanesia" = "#54278f",
  "Micronesia" = "#756bb1",
  "Polynesia" =  "#9e9ac8",
  "Australia NZ" = "#cbc9e2",
  "Antarctica" = "#702963",
  "ABNJ" = "grey50"
)
## Make the "big" one first and then the smaller one, for species level the biggest is fish so set the others to that
dataset <- fishRC
name <- names(fishRC)

dataset_FISH <- dataset %>%
  left_join(groups, by=join_by(source_eezName == value)) %>%
  rename(source_eezGroup = group) %>%
  left_join(groups, by=join_by(target_eezName == value)) %>%
  rename(target_eezGroup = group) %>%
  mutate(target_eezGroup = gsub("_", " ", target_eezGroup),
         source_eezGroup = gsub("_", " ", source_eezGroup))

# pull nodes out
nodesFISH <- unique(c(dataset_FISH$source_eezGroup, dataset_FISH$target_eezGroup))

# empty matrix
adjacency_matrixFISH <- matrix(0, nrow = length(nodesFISH), ncol = length(nodesFISH), dimnames = list(nodesFISH, nodesFISH))

# Iterate over the pair data and update the adjacency matrix
for (i in 1:nrow(dataset_FISH)) {
  node1 <- dataset_FISH$source_eezGroup[i]
  node2 <- dataset_FISH$target_eezGroup[i]
  adjacency_matrixFISH[node1, node2] <- adjacency_matrixFISH[node1, node2] + 1
}

# make list of regions in this subset
subset_regionsFISH <- nodesFISH[nodesFISH %in% gsub("_", " ", regions)]
subset_regionsFISH <- subset_regionsFISH[match(gsub("_", " ", regions), subset_regionsFISH)]
subset_regionsFISH <- na.omit(subset_regionsFISH)

linkcolFISH <- linkcol[names(linkcol) %in% subset_regionsFISH]

# Save the chord diagram as an image with high resolution
png(filename = paste0("ChordDiagrams/Scaled/", "FISH", "_region_CD_scaled.png"), width = 10, height = 10, units = "in", res = image_resolution)
plot.new()
par(cex = 0.5, mar = c(0, 0, 0, 0))
circos.clear()

chordDiagram(adjacency_matrixFISH, 
             order = subset_regionsFISH,
             grid.col = linkcolFISH,
             transparency = 0.3,
             annotationTrack = "grid", #group = group, 
             big.gap = 3, small.gap = 1,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(adjacency_matrixFISH))))))

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  xplot = get.cell.meta.data("xplot")
  sector.name = get.cell.meta.data("sector.index")
  if(abs(xplot[2] - xplot[1]) < 10) {
    circos.text(CELL_META$xcenter, ylim[1] + cm_h(1), 
                sector.name, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), cex = 1.3, font = 1)
  } else {
    circos.text(CELL_META$xcenter, ylim[1] + cm_h(1), 
                sector.name, facing = "inside",
                niceFacing = TRUE, adj = c(0.5, 0), cex = 1.3, font = 1)
  }
  circos.axis(h = "bottom",
              labels.cex = .8,
              sector.index = sector.name
  )
}, bg.border = NA)
dev.off()

### Now make the smaller ones, accounting for the gap based on the fish one
for (i in 1:length(sourceRCs)) {
  dataset <- sourceRCs[[i]]
  name <- names(sourceRCs)[i]
  
  dataset <- dataset %>%
    left_join(groups, by=join_by(source_eezName == value)) %>%
    rename(source_eezGroup = group) %>%
    left_join(groups, by=join_by(target_eezName == value)) %>%
    rename(target_eezGroup = group) %>%
    mutate(target_eezGroup = gsub("_", " ", target_eezGroup),
           source_eezGroup = gsub("_", " ", source_eezGroup))
  
  # pull nodes out
  nodes <- unique(c(dataset$source_eezGroup, dataset$target_eezGroup))
  
  # empty matrix
  adjacency_matrix <- matrix(0, nrow = length(nodes), ncol = length(nodes), dimnames = list(nodes, nodes))
  
  # Iterate over the pair data and update the adjacency matrix
  for (i in 1:nrow(dataset)) {
    node1 <- dataset$source_eezGroup[i]
    node2 <- dataset$target_eezGroup[i]
    adjacency_matrix[node1, node2] <- adjacency_matrix[node1, node2] + 1
  }
  
  # make list of regions in this subset
  subset_regions <- nodes[nodes %in% gsub("_", " ", regions)]
  subset_regions <- subset_regions[match(gsub("_", " ", regions), subset_regions)]
  subset_regions <- na.omit(subset_regions)
  
  linkcol2 <- linkcolS[names(linkcolS) %in% subset_regions]
  
  # Save the chord diagram as an image with high resolution
  png(filename = paste0("ChordDiagrams/Scaled/", name, "_region_CD_scaled.png"), width = 10, height = 10, units = "in", res = image_resolution)
  plot.new()
  par(cex = 0.5, mar = c(0, 0, 0, 0))
  circos.clear()
  
  chordDiagram(adjacency_matrix, 
               order = subset_regions,
               grid.col = linkcol2,
               transparency = 0.3,
               annotationTrack = "grid", #group = group, 
               #directional = 1,
               big.gap = 3, small.gap = 1,
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(adjacency_matrix))))))
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    xplot = get.cell.meta.data("xplot")
    sector.name = get.cell.meta.data("sector.index")
    if(abs(xplot[2] - xplot[1]) < 10) {
      circos.text(CELL_META$xcenter, ylim[1] + cm_h(2.5), 
                  sector.name, facing = "clockwise",
                  niceFacing = TRUE, adj = c(0, 0.5), cex = 1.9, font = 1)
    } else {
      circos.text(CELL_META$xcenter, ylim[1] + cm_h(2.5), 
                  sector.name, facing = "inside",
                  niceFacing = TRUE, adj = c(0.5, 0), cex = 1.9, font = 1)
    }
    circos.axis(h = "bottom",
                labels.cex = .8,
                sector.index = sector.name
    )
  }, bg.border = NA)
  dev.off()
  
  cat("Chord diagram for", name, "created and saved as", paste0(name, "_region_CD_mini.png"), "\n\n")
}



### Trying out edgebundleR
library(edgebundleR)
library(igraph)

dataset <- routesconnections
name <- names(routesconnections)

# pull nodes out
nodes <- unique(c(dataset$source_eezName, dataset$target_eezName))

# empty matrix
adjacency_matrix <- matrix(0, nrow = length(nodes), ncol = length(nodes), dimnames = list(nodes, nodes))

# Iterate over the pair data and update the adjacency matrix
for (i in 1:nrow(dataset)) {
  node1 <- dataset$source_eezName[i]
  node2 <- dataset$target_eezName[i]
  adjacency_matrix[node1, node2] <- adjacency_matrix[node1, node2] + 1
}

# Find nodes with no connections
no_connections <- which(rowSums(adjacency_matrix) == 0)

# Remove nodes with no connections from the adjacency matrix
adjacency_matrix <- adjacency_matrix[-no_connections, -no_connections]

g <- graph_from_adjacency_matrix(adjacency_matrix)
edgebundle(g,
           tension = 0.9,
           padding = 150,
           fontsize = 10)

