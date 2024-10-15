library(ggplot2)
library(scater)
library(scran)
library(here)
library(RColorBrewer)

## save directories
plot_dir = here("plots", "07_annotation_and_characterization")
processed_dir = here("processed-data","07_annotation_and_characterization")

# load sce
sce <- readRDS(here(processed_dir, "sce_excit_final_subclusters_annotated.rds"))


# ======= Macaque samples ========

#subset to different species
sce_macaque <- sce[,which(colData(sce)$species == "macaque")]
colnames(colData(sce_macaque))
sce_macaque$Sample_num

# read in new sample labels
new_labels <- read.csv(here("raw-data","sampleinfo", "relabeling_macaque_metadata.csv"))

# add new Sample_num coldata that's Sample1 - Sample 35 --- to new_labels
new_labels$Sample_num <- paste0("Sample", 1:35)

# make new coldata for "$Sample_num _ ONPRC.ID _ Region _ Dorsal.Ventral" to new_labels
new_labels$Final <- paste0(new_labels$Sample_num, "_", new_labels$ONPRC.ID, "_", new_labels$Region, "_", new_labels$Dorsal.Ventral)


# now merge this to sce_macaque coldata using new_label$Current.Label.in.Metadata to match with sce_macaque$Sample
sce_macaque$Sample <- factor(sce_macaque$Sample)
new_labels$Current.Label.in.Metadata <- factor(new_labels$Current.Label.in.Metadata)
new_coldata <- merge(colData(sce_macaque), new_labels, by.x = "Sample", by.y = "Current.Label.in.Metadata", all=TRUE)

# now replace coldata in sce_macaque with new_coldata
colData(sce_macaque) <- new_coldata


library(dplyr)

# Create the table
celltype_table <- table(sce_macaque$fine_celltype, sce_macaque$Final, sce_macaque$Subregion)

# Convert the table to a data frame
df <- as.data.frame(celltype_table)

# Rename columns for clarity
colnames(df) <- c("CellType", "Sample", "Subregion", "Count")

# Calculate proportions within each Sample and Subregion
df <- df %>%
  group_by(Sample, Subregion) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  filter(!is.na(Proportion))

# Define custom colors for cell types
nb.cols <- length(unique(sce$fine_celltype))
colors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)

celltype_colors <- setNames(colors, unique(df$CellType))

# Create the stacked bar plot with facet wrap for subregion
png(here(plot_dir, "stackedbars_macaque_sample_by_excit_subregion.png"), width = 8, height = 10, units = "in", res = 300)
ggplot(df, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Proportion", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12)) +  # Customize facet labels
  coord_flip() +
  scale_fill_manual(values = celltype_colors) +
  facet_grid(Subregion~., scales="free", space="free")  # Adjust the number of columns in the facet wrap
dev.off()










# ======= Baboon samples ========

#subset to different species
sce_baboon <- sce[,which(colData(sce)$species == "baboon")]

# Create the table
celltype_table <- table(sce_baboon$fine_celltype, sce_baboon$Sample, sce_baboon$Subregion)

# Convert the table to a data frame
df <- as.data.frame(celltype_table)

# Rename columns for clarity
colnames(df) <- c("CellType", "Sample", "Subregion", "Count")

# Calculate proportions within each Sample and Subregion
df <- df %>%
  group_by(Sample, Subregion) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  filter(!is.na(Proportion))

# Define custom colors for cell types
nb.cols <- length(unique(sce$fine_celltype))
colors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)

celltype_colors <- setNames(colors, unique(df$CellType))

# Create the stacked bar plot with facet wrap for subregion
png(here(plot_dir, "stackedbars_baboon_sample_by_cexcit_subregion.png"), width = 8, height = 6, units = "in", res = 300)
ggplot(df, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Proportion", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12)) +  # Customize facet labels
  coord_flip() +
  scale_fill_manual(values = celltype_colors) +
  facet_grid(Subregion~., scales="free", space="free")  # Adjust the number of columns in the facet wrap
dev.off()




# ======= Human samples ========
#subset to different species
sce_human <- sce[,which(colData(sce)$species == "human")]

# Create the table
celltype_table <- table(sce_human$fine_celltype, sce_human$Sample)

# Convert the table to a data frame
df <- as.data.frame(celltype_table)

# Rename columns for clarity
colnames(df) <- c("CellType", "Sample", "Count")

# Calculate proportions within each Sample and Subregion
df <- df %>%
  group_by(Sample) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  filter(!is.na(Proportion))

# Define custom colors for cell types
nb.cols <- length(unique(sce$fine_celltype))
colors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)

celltype_colors <- setNames(colors, unique(df$CellType))

# Create the stacked bar plot with facet wrap for subregion
png(here(plot_dir, "stackedbars_human_sample_by_excit.png"), width = 6, height = 5, units = "in", res = 300)
ggplot(df, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Proportion", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12)) +  # Customize facet labels
  coord_flip() +
  scale_fill_manual(values = celltype_colors) 
dev.off()