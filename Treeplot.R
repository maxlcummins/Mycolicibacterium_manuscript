# Load necessary libraries
library(tidyverse)     # For data manipulation and visualization
library(ggtree)        # For visualizing phylogenetic trees
library(ggtreeExtra)   # For enhancing ggtree plots with additional layers

# Read in FastANI pairwise comparison results (tab-delimited, no column headers)
fastani <- read_delim("D:/sequence_runs/2025/March/NTM_25012024/tree_manuscript/fastANI_output.txt", col_names = FALSE)

# Rename columns to something meaningful
fastani <- fastani %>% rename('name' = X1, 'query' = X2, 'ANI' = X3)

# Clean up file paths and extensions from 'name' column
fastani <- fastani %>% mutate(name = str_replace(name, 'assemblies/', ''))
fastani <- fastani %>% mutate(name = str_replace(name, '.fna', ''))

# Clean up file paths and extensions from 'query' column
fastani <- fastani %>% mutate(query = str_replace(query, 'assemblies/', ''))
fastani <- fastani %>% mutate(query = str_replace(query, '.fna', ''))

# Read the phylogenetic tree in Newick format
tree <- read.tree("D:/sequence_runs/2025/March/NTM_25012024/tree_manuscript/tree.dnd")

# Read mapping of genome IDs to species names2
species_ids <- read_delim("D:/sequence_runs/2025/March/NTM_25012024/tree_manuscript/species_ids.txt", delim = " ", col_names = FALSE)

# Clean up genome IDs and rename columns
species_ids <- species_ids %>% 
  mutate(X2 = str_replace(X2, '.fna.gz', '')) %>% 
  rename('name' = X2)

# Extract tip labels from the tree and convert to a data frame
tree_names <- tree$tip.label %>% as.data.frame() %>% rename('name' = '.')

# Join species IDs to tree tip labels; label missing species as "Unknown"
species_ids <- full_join(tree_names, species_ids, by = 'name') %>% 
  rename('species' = X1) %>% 
  mutate(species = replace_na(species, 'Unknown'))

# Manually edit species IDs for unknowns
species_ids <- species_ids %>%
  mutate(species = if_else(name == "AP023287", "sp.", species)) %>% 
  mutate(species = if_else(name == "GCF_014218295.1", "sp.", species)) %>% 
  mutate(species = if_else(name == "TPG_MSP_01", "sp.", species)) %>% 
  mutate(species = if_else(name == "NZ_AP022586", "litorale", species)) %>% 
  mutate(species = if_else(name == "GCF_017352375.1", "sp.", species)) %>% 
  mutate(species = if_else(name == "GCF_023015925.1", "sp.", species))

# Bind in the genus abbreviation
species_ids <- species_ids %>% mutate(species = paste0("M. ", species))

# Create a supp table
Supp_table_1 <- left_join(fastani, species_ids, by = 'name') %>% select(-X4, -X5)

# Write our supp table to file
write_delim(Supp_table_1, "D:/sequence_runs/2025/March/NTM_25012024/tree_manuscript/Supplementary_table_1", delim = " ")

# Italicise our species names
species_ids <- species_ids %>% mutate(species = paste0("italic('", species, "')"))

# Reroot the tree
rerooted_tree <- ape::root(tree, outgroup = "GCF_022179545.1", resolve.root = TRUE)

# Create base ggtree plot with species labels
p1 <- ggtree(rerooted_tree) %<+% species_ids +
  geom_tiplab(align = TRUE, size = 4) +  # genome labels
  geom_tiplab(aes(label = species), parse = TRUE, size = 4, offset = 0.1, align = TRUE, linetype = NULL) +  # species labels
  xlim(c(NA, 0.25))  # adjust x-axis limits

# Extract the order of tips from the plotted tree
tip_order <- p1$data %>% 
  as.data.frame() %>% 
  filter(isTip == TRUE) %>% 
  arrange(desc(y)) %>% 
  pull(label)

# Reorder 'query' factor levels in fastani data to match tree tip order
fastani <- fastani %>% mutate(query = factor(query, levels = tip_order))

# Add a heatmap of ANI values as a tile layer beside the tree
p1 +
  geom_fruit(
    data = fastani,
    geom = geom_tile,
    mapping = aes(x = query, y = name, fill = ANI),
    offset = 0.95,
    pwidth = 3,
    axis.params = list(
      axis = "x",
      text.angle = 90,
      text.size = 3,
      line.size = 0,
      vjust = 0.5,
      hjust = 1
    )
  ) +
  # Set custom color gradient for ANI values
  scale_fill_gradientn(
    colours = rev(c("#d73027", "#fc8d59", "#fee08b", "#91bfdb", "#4575b4")),
    values = scales::rescale(c(79, 85, 90, 95, 100)),
    limits = c(79, 100),
    name = "ANI"
  ) +
  # Adjust y- and x-axis limits to accommodate extra plot elements
  ylim(c(-2, NA)) +
  xlim(c(NA, 0.7)) + theme(
    legend.margin = margin(0, 10, 0, -20),     
    legend.box.spacing = unit(0, "pt"),       
    plot.margin = unit(c(1,1,1,1), "mm")    
  )
