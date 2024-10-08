# The chemosynthetic microbiome of deep-sea hydrothermal vents across space and time

## Metatranscriptomic analysis

# Load packages
library(ggplot2)
library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)
library(ggforce)
library(readxl)
library(RColorBrewer)

# LVP files
# Read in excel spreadsheets
lvp4_gff <- read_excel(here("../dataset/raw_data_for_mod_gene_abundance/lvp4_gff_modified_split.xlsx"))
lvp4_tpm <- read_excel(here("../dataset/raw_data_for_mod_gene_abundance/lvp4_tpm.xlsx"))
lvp4_prod_names <- read_excel(here("../dataset/raw_data_for_mod_gene_abundance/lvp4_product_names.xlsx"))
lvp4_phylodist <- read_excel(here("../dataset/raw_data_for_mod_gene_abundance/lvp4_phylodist.xlsx"))

# First merge .gff and .tpm by gene_ID
merged_df_lvp <- merge(lvp4_gff, lvp4_tpm, by = "gene_ID", all = TRUE)

#Then merge that data frame with prod_names by locus_tag
merged_df_lvp4 <- merge(merged_df_lvp, lvp4_prod_names, by = "locus_tag", all = TRUE)

# Change column header cv88_tpm
colnames(merged_df_lvp4)[colnames(merged_df_lvp4) == "lvp4_tpm"] <- "tpm"

# Add Sample_ID column
merged_df_lvp4$Sample_ID <- "lvp4"

# The taxonomy can also be added (at x resolution) to this as well using the phylodist file
merged_df_lvp4p <- merge(merge(merged_df_lvp4, lvp4_prod_names, by = "locus_tag", all = TRUE), lvp4_phylodist, by = "locus_tag", all = TRUE)

# Clean up the CV88 files
# Read in excel spreadsheets
cv88_gff <- read_excel(here("../dataset/raw_data_for_mod_gene_abundance/cv88_gff_modified.xlsx"))
cv88_tpm <- read_excel(here("../dataset/raw_data_for_mod_gene_abundance/cv88_tpm.xlsx"))
cv88_prod_names <- read_excel(here("../dataset/raw_data_for_mod_gene_abundance/cv88_product_names.xlsx"))
cv88_phylodist <- read_excel(here("../dataset/raw_data_for_mod_gene_abundance/cv88_phylodist.xlsx"))

# Clean up the CV88 .gff file
cv88_to_split <- read_excel(here("../dataset/raw_data_for_mod_gene_abundance/cv88_gff_tofix.xlsx"))
cv88_exon <- subset(cv88_to_split, c == "exon")
cv88_rRNA <- subset(cv88_to_split, c == "rRNA")
cv88_tRNA <- subset(cv88_to_split, c == "tRNA")
cv88_CDS <- subset(cv88_to_split, c == "CDS")
cv88_repeat_region <- subset(cv88_to_split, c == "repeat_region")

# Keep: extract locus tag and gene ID
# Extract 'locus_tag' from column 'd' and create a new column
cv88_rRNA <- cv88_rRNA %>%
  mutate(locus_tag = str_extract(d, "locus_tag=(\\S+?);") %>% str_remove_all("locus_tag=|;"))
cv88_rRNA <- cv88_rRNA %>%
  mutate(gene_ID = str_extract(d, "ID=(\\S+?);") %>% str_remove_all("ID=|;"))
cv88_rRNA <- cv88_rRNA %>% select(-a, -b, -c, -d)

cv88_tRNA <- cv88_tRNA %>%
  mutate(locus_tag = str_extract(d, "locus_tag=(\\S+?);") %>% str_remove_all("locus_tag=|;"))
cv88_tRNA <- cv88_tRNA %>%
  mutate(gene_ID = str_extract(d, "ID=(\\S+?);") %>% str_remove_all("ID=|;"))
cv88_tRNA <- cv88_tRNA %>% select(-a, -b, -c, -d)

cv88_CDS <- cv88_CDS %>%
  mutate(locus_tag = str_extract(d, "locus_tag=(\\S+?);") %>% str_remove_all("locus_tag=|;"))
cv88_CDS <- cv88_CDS %>%
  mutate(gene_ID = str_extract(d, "ID=(\\S+?);") %>% str_remove_all("ID=|;"))
cv88_CDS <- cv88_CDS %>% select(-a, -b, -c, -d)

save.image()

# Concatonate and merge CV88 files
# Concatonate the dataframes and clean up
concatenated_df <- rbind(cv88_rRNA, cv88_tRNA, cv88_CDS)

# First merge .gff and .tpm by gene_ID
merged_df <- merge(concatenated_df, cv88_tpm, by = "gene_ID", all = TRUE)

#Then merge that data frame with prod_names by locus_tag
merged_df_cv88 <- merge(merged_df, cv88_prod_names, by = "locus_tag", all = TRUE)

# Change column header cv88_tpm
colnames(merged_df_cv88)[colnames(merged_df_cv88) == "cv88_tpm"] <- "tpm"

# Add Sample_ID column
merged_df_cv88$Sample_ID <- "cv88"

# The taxonomy can also be added (at x resolution) to this as well using the phylodist file
merged_df_cv88p <- merge(merge(merged_df_cv88, cv88_prod_names, by = "locus_tag", all = TRUE), cv88_phylodist, by = "locus_tag", all = TRUE)

merged_df_cv88p_subset <- merged_df_cv88p[, 1:14]

# Concatonate CV88 and LVP4 files
all_df1 <- rbind(merged_df_cv88p_subset, merged_df_lvp4p)

# Subset prokaryotic reads
all_df <- subset(all_df1, Domain %in% c("Archaea", "Bacteria"))

# Extract genes of interest
# Extract rows by name
# ATP citrate lyase (aclAB)
ATP_citrate_lyase_aclAB <- c("ATP-citrate lyase alpha-subunit", "ATP-citrate lyase beta-subunit")

# Type-IV hydrogenase (hyfBEF)
TypeIV_hydrogenase_hyfBEF <- c("hydrogenase-4 component B", "hydrogenase-4 component E", "hydrogenase-4 component F")

# Hydrogenase maturation proteins (hypABCDEF)
Hydrogenase_maturation_proteins_hypABCDEF <- c("hydrogenase nickel incorporation protein HypA/HybF","hydrogenase nickel incorporation protein HypA/HybF/putative two-component system protein, hydrogenase maturation factor HypX/HoxX","Zn finger protein HypA/HybF (possibly regulating hydrogenase expression)", "hydrogenase nickel incorporation protein HypB","hydrogenase expression/formation protein HypC","hydrogenase expression/formation protein HypD","hydrogenase expression/formation protein HypE","Hydrogenase maturation factor HypF (carbamoyltransferase)","hydrogenase maturation protein HypF")

# thiol peroxidase
thiol_peroxidase <- c("thiol peroxidase, atypical 2-Cys peroxiredoxin")

# superoxide reductase
superoxide_reductase <- c("superoxide reductase")

# superoxide dismutase
superoxide_dismutase <- c("superoxide dismutase, Fe-Mn family")

# cytochrome c peroxidase
cytochrome_c_peroxidase <- c("cytochrome c peroxidase")

# catalase-peroxidase
catalase_peroxidase <- c("catalase-peroxidase")

# alkyl hydroperoxide reductase subunit C (ahp)
alkyl_hydroperoxide_reductase_subunit_C_ahp <- c("peroxiredoxin (alkyl hydroperoxide reductase subunit C)")

# Sulfur-oxidizing protein (soxABCXYZ)
Sulfur_oxidizing_protein_soxABCXYZ <- c("sulfur-oxidizing protein SoxA","sulfur-oxidizing protein SoxB","sulfane dehydrogenase subunit SoxC","sulfur-oxidizing protein SoxX","sulfur-oxidizing protein SoxY","sulfur-oxidizing protein SoxZ")

# Sulfide:quinone oxidoreductase (sqr)
Sulfide_quinone_oxidoreductase_sqr <- c("sulfide:quinone oxidoreductase")

# Polysulfide reductase (psrA)
Polysulfide_reductase_psrA <- c("thiosulfate reductase / polysulfide reductase chain A")

# Periplasmic nitrate reductase (napABCDFGH)
Periplasmic_nitrate_reductase_napABCDFGH <- c("nitrate reductase NapA","cytochrome c-type protein NapB","cytochrome c-type protein NapC","nitrate reductase NapD","	
ferredoxin-type protein NapF","ferredoxin-type protein NapG","ferredoxin-type protein NapH")

# Nitric oxide reductase (norBCDEQW)
Nitric_oxide_reductase_norBCDEQW <- c("nitric oxide reductase subunit B","nitric oxide reductase subunit C","nitric oxide reductase NorD protein","nitric oxide reductase NorE protein","nitric oxide reductase NorQ protein","nitric oxide reductase FlRd-NAD(+) reductase")

# nitrite reductase (nirBD)
#group15 <- c("")??? Not sure, set aside
  
# membrane-bound nitrate reductase (narGHJ)
#group16 <-c("")??? Not sure, set aside
  
# Formate-dependent nitrite reductase (nrfAD) # only D
Formate_dependent_nitrite_reductase_nrfAD <- c("Formate-dependent nitrite reductase, membrane component NrfD")

# Assimilatory nitrate reductase catalytic subunity (nas)
Assimilatory_nitrate_reductase_catalytic_subunite_nas <- c("assimilatory nitrate reductase catalytic subunit")

# Hydroxylamine reductase (har)
Hydroxylamine_reductase_har <- c("nitrite reductase (NO-forming) / hydroxylamine reductase", "hydroxylamine reductase")

# Cytochrome d ubiquinol oxidase
Cytochrome_d_ubiquinol_oxidase <- c("cytochrome d ubiquinol oxidase subunit I","cytochrome d ubiquinol oxidase subunit II")

# Cytochrome c oxidase cbb3-type
Cytochrome_c_oxidase_cbb3_type <- c("Cbb3-type cytochrome oxidase, cytochrome c subunit","Cbb3-type cytochrome oxidase, subunit 1","	Cbb3-type cytochrome oxidase, subunit 3","cytochrome c oxidase cbb3-type subunit 1","cytochrome c oxidase cbb3-type subunit 1/cytochrome c oxidase cbb3-type subunit I/II","cytochrome c oxidase cbb3-type subunit 2","cytochrome c oxidase cbb3-type subunit 3","cytochrome c oxidase cbb3-type subunit 4","Cytochrome C oxidase, cbb3-type, subunit III")

# Cytochrome oxidase
Cytochrome_oxidase <- c("cytochrome c oxidase subunit 1","cytochrome c oxidase subunit 2","cytochrome c oxidase subunit 3","Cytochrome C oxidase subunit IV")

# L-lactate dehydrogenase complex protein LldEFG
L_lactate_dehydrogenase_complex_protein_LldEFG <- c("L-lactate dehydrogenase complex protein LldE","L-lactate dehydrogenase complex protein LldF","L-lactate dehydrogenase complex protein LldG")

# Methylmalonyl-CoA mutase (mcm)
Methylmalonyl_CoA_mutase_mcm <- c("methylmalonyl-CoA mutase","methylmalonyl-CoA mutase, C-terminal domain","	
methylmalonyl-CoA mutase, N-terminal domain")

# Acetyl-CoA carboxylase (accABCD) #maybe add more
Acetyl_CoA_carboxylase_accABCD <- c("Acetyl-CoA carboxylase alpha subunit","	Acetyl-CoA carboxylase beta subunit")

# 3-HP-4HB or DC-4HB
#group26 <- c("")?? Not sure, set aside
  
# Methylenetetrahydrofolate reductase (metF)
Methylenetetrahydrofolate_reductase_metF <- c("5,10-methylenetetrahydrofolate reductase")

# Carbon-monoxide dehydrogenase large subunit (codH)
Carbon_monoxide_dehydrogenase_large_subunit_codH <- c("carbon-monoxide dehydrogenase large subunit")

# Formate dehydrogenase (fdh)
Formate_dehydrogenase_fdh <- c("formate dehydrogenase","formate dehydrogenase beta subunit","formate dehydrogenase iron-sulfur subunit","formate dehydrogenase major subunit","	
formate dehydrogenase subunit gamma")

# Pyruvate:ferredoxin oxidoreductase (pfor)
Pyruvate_ferredoxin_oxidoreductase_pfor <- c("Pyruvate:ferredoxin oxidoreductase","Pyruvate:ferredoxin oxidoreductase, beta subunit","Pyruvate:ferredoxin oxidoreductase, gamma subunit","pyruvate ferredoxin oxidoreductase gamma subunit","pyruvate ferredoxin oxidoreductase delta subunit","pyruvate ferredoxin oxidoreductase beta subunit","pyruvate ferredoxin oxidoreductase alpha subunit")

# 2-oxoglutarate ferredoxin oxidoreductase (korABGD)
Two_oxoglutarate_ferredoxin_oxidoreductase_korABGD <- c("2-oxoglutarate ferredoxin oxidoreductase subunit alpha","2-oxoglutarate ferredoxin oxidoreductase subunit beta","2-oxoglutarate ferredoxin oxidoreductase subunit delta","2-oxoglutarate ferredoxin oxidoreductase subunit gamma")

# Define group names
group_names <- c(
  "ATP_citrate_lyase_aclAB",
  "TypeIV_hydrogenase_hyfBEF",
  "Hydrogenase_maturation_proteins_hypABCDEF",
  "thiol_peroxidase",
  "superoxide_reductase",
  "superoxide_dismutase",
  "cytochrome_c_peroxidase",
  "catalase_peroxidase",
  "alkyl_hydroperoxide_reductase_subunit_C_ahp",
  "Sulfur_oxidizing_protein_soxABCXYZ",
  "Sulfide_quinone_oxidoreductase_sqr",
  "Polysulfide_reductase_psrA",
  "Periplasmic_nitrate_reductase_napABCDFGH",
  "Nitric_oxide_reductase_norBCDEQW",
  "Formate_dependent_nitrite_reductase_nrfAD",
  "Assimilatory_nitrate_reductase_catalytic_subunite_nas",
  "Hydroxylamine_reductase_har",
  "Cytochrome_d_ubiquinol_oxidase",
  "Cytochrome_c_oxidase_cbb3_type",
  "Cytochrome_oxidase",
  "L_lactate_dehydrogenase_complex_protein_LldEFG",
  "Methylmalonyl_CoA_mutase_mcm",
  "Acetyl_CoA_carboxylase_accABCD",
  "Methylenetetrahydrofolate_reductase_metF",
  "Carbon_monoxide_dehydrogenase_large_subunit_codH",
  "Formate_dehydrogenase_fdh",
  "Pyruvate_ferredoxin_oxidoreductase_pfor",
  "Two_oxoglutarate_ferredoxin_oxidoreductase_korABGD"
)

# Define group definitions
group_definitions <- lapply(group_names, function(group_name) get(group_name))

# Subset and assign gene categories
# Create an empty data frame to store the results
goi <- data.frame()

# Loop through each group definition
for (i in seq_along(group_definitions)) {
  # Subset the data based on the current group definition
  subset_data <- subset(all_df, product_name.x %in% group_definitions[[i]])
  
  # Add a column for gene categories and assign the current group name
  subset_data$gene_categories <- group_names[i]
  
  # Append the subset data to the goi data frame
  goi <- rbind(goi, subset_data)
}
# Remove NAs
goi <- goi[complete.cases(goi$tpm), ]

# Check number of gene categories
num_unique <- goi %>% 
  distinct(gene_categories) %>% 
  nrow()

print(num_unique)

# Sum gene categories by sample
library(dplyr)

# Step 1: Split the dataframe based on Sample_ID, keeping only tpm and gene_categories
split_data <- split(select(goi, Sample_ID, tpm, gene_categories), goi$Sample_ID)

# Step 2: Calculate the summed value for tpm of each gene_categories
summed_data <- lapply(split_data, function(df) {
  df %>%
    group_by(gene_categories) %>%
    summarise(total_tpm = sum(tpm, na.rm = TRUE))
})

# Step 3: Add Sample_ID back as a column to each split dataframe
for (i in seq_along(summed_data)) {
  summed_data[[i]]$Sample_ID <- names(summed_data)[i]
}

# Step 4: Merge the split dataframes back together, excluding product_name.x
merged_data <- do.call(rbind, summed_data)

# Print the merged dataframe
print(merged_data)

# To plot: Bubbleplot of gene abundances

# Convert gene_categories to a factor with custom levels
merged_data$gene_categories <- factor(merged_data$gene_categories, levels = c(
  "Two_oxoglutarate_ferredoxin_oxidoreductase_korABGD",
  "Pyruvate_ferredoxin_oxidoreductase_pfor",
  "ATP_citrate_lyase_aclAB",
  "Formate_dehydrogenase_fdh",
  "Carbon_monoxide_dehydrogenase_large_subunit_codH",
  "Methylenetetrahydrofolate_reductase_metF",
  "Acetyl_CoA_carboxylase_accABCD",
  "Methylmalonyl_CoA_mutase_mcm",
  "L_lactate_dehydrogenase_complex_protein_LldEFG",
  "Cytochrome_oxidase",
  "Cytochrome_c_oxidase_cbb3_type",
  "Cytochrome_d_ubiquinol_oxidase",
  "Hydroxylamine_reductase_har",
  "Assimilatory_nitrate_reductase_catalytic_subunite_nas",
  "Formate_dependent_nitrite_reductase_nrfAD",
  "Nitric_oxide_reductase_norBCDEQW",
  "Periplasmic_nitrate_reductase_napABCDFGH",
  "Polysulfide_reductase_psrA",
  "Sulfide_quinone_oxidoreductase_sqr",
  "Sulfur_oxidizing_protein_soxABCXYZ",
  "alkyl_hydroperoxide_reductase_subunit_C_ahp",
  "catalase_peroxidase",
  "cytochrome_c_peroxidase",
  "superoxide_dismutase",
  "superoxide_reductase",
  "thiol_peroxidase",
  "Hydrogenase_maturation_proteins_hypABCDEF",
  "TypeIV_hydrogenase_hyfBEF"
))

# Set factor layer
merged_data <- merged_data %>%
  mutate(metabolic_group = case_when(
    gene_categories %in% c("Two_oxoglutarate_ferredoxin_oxidoreductase_korABGD", "Pyruvate_ferredoxin_oxidoreductase_pfor",
                           "ATP_citrate_lyase_aclAB", "Formate_dehydrogenase_fdh", "Carbon_monoxide_dehydrogenase_large_subunit_codH",
                           "Methylenetetrahydrofolate_reductase_metF", "Acetyl_CoA_carboxylase_accABCD", "Methylmalonyl_CoA_mutase_mcm",
                           "L_lactate_dehydrogenase_complex_protein_LldEFG") ~ "carbon",
    gene_categories %in% c("Cytochrome_oxidase", "Cytochrome_c_oxidase_cbb3_type", "Cytochrome_d_ubiquinol_oxidase") ~ "oxygen",
    gene_categories %in% c("Hydroxylamine_reductase_har", "Assimilatory_nitrate_reductase_catalytic_subunite_nas",
                           "Formate_dependent_nitrite_reductase_nrfAD", "Nitric_oxide_reductase_norBCDEQW",
                           "Periplasmic_nitrate_reductase_napABCDFGH") ~ "nitrogen",
    gene_categories %in% c("Polysulfide_reductase_psrA", "Sulfide_quinone_oxidoreductase_sqr", "Sulfur_oxidizing_protein_soxABCXYZ") ~ "sulfur",
    gene_categories %in% c("alkyl_hydroperoxide_reductase_subunit_C_ahp", "catalase_peroxidase", "cytochrome_c_peroxidase",
                           "superoxide_dismutase", "superoxide_reductase", "thiol_peroxidase") ~ "ROS defense",
    gene_categories %in% c("Hydrogenase_maturation_proteins_hypABCDEF", "TypeIV_hydrogenase_hyfBEF") ~ "hydrogen",
    TRUE ~ NA_character_
  ))

summary(merged_data$total_tpm)

# Plot
bubble_plot <- ggplot(merged_data, aes(x = Sample_ID, y = gene_categories, size = total_tpm, fill = metabolic_group)) +
  geom_point(alpha = 0.7,shape=21) +  # Add points with transparency
  # scale_size_continuous(range = c(3, 10)) +  # Adjust the range of bubble sizes
  scale_size_continuous(limits = c(0, 5500), range = c(1,18), breaks = c(0,100,1000,5000)) +
  labs(title = "Gene transcript abundances", x = "Sample ID", y = "Product Name", size = "TPM") +  # Add labels
  scale_fill_manual(labels = c("Carbon", "Hydrogen","Nitrogen","Oxygen","ROS defense","Sulfur"), 
                    values = c("lightblue", "gold","green3","violet","red3","navy")) +
  theme_dpml() +
  theme(legend.position = "right",
    text = element_text(family = "Helvetica", size = 10))

bubble_plot
#ggsave("../plot/bubble_plot.svg",width=14,height=16)

write.csv(bubble_plot$data,"../plot/bubbleplot_data_table.csv")

save.image()

# Plot pie graphs
#Warning, scatterpie will conflict with other packages, use it and detach it.
library(scales)
library(scatterpie)

theme_dpml <- function() {
  theme(
    text = element_text(family = "Helvetica", size = 20),
    axis.title = element_text(size = rel(1.2)),
    axis.text = element_text(size = rel(1)),
    
    legend.title = element_text(size = rel(1.2)),
    legend.text = element_text(size = rel(1)),
    
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    
    legend.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white"),
    strip.background = element_rect(fill = "black"),
    
    # Custom element for a rectangle in the plot area
    panel.border = element_rect(
    fill = NA,
    colour = 'black')
      )
}
options(repr.plot.width=16, repr.plot.height=12)


save.image()

# Prepare dataframe for pie charts

# Group by gene category, sum, and take unique values for each gene category
goi2_MS <- goi %>%
  group_by(Sample_ID, gene_categories,Genus) %>%
  summarise(tpm = sum(tpm)) # Summarize tpm by Sample_ID and gene_categories

goi2_MS

# Merge goi2 with goi to include taxonomy data
goi3_MS <- goi2_MS %>%
  left_join(goi, by = c("Sample_ID", "gene_categories","Genus"))
goi3_MS

# Add a numeric column for x-axis 
goi4_MS <- goi3_MS %>%
  mutate(x = as.numeric(factor(Sample_ID)), 
         x = scales::rescale(x, to = range(Sample_ID)))

# Subsetting columns from goi4
subset_goi4_MS <- goi4_MS[, c("Sample_ID", "gene_categories", "tpm.x", "Genus")]

# Keep only unique rows
unique_subset_goi4_MS <- unique(subset_goi4_MS)

# Ungroup the dataframe
ungrouped_unique_subset_goi4_MS <- unique_subset_goi4_MS %>% 
  ungroup()

ungrouped_unique_subset_goi4_MS

# Pie Charts: Make Sample_ID and product_name.x numeric column equivalent
# Convert Sample_ID to factor with levels sorted
goi_n_MS <- ungrouped_unique_subset_goi4_MS %>%
  mutate(sampleID_num = as.numeric(factor(Sample_ID)),  # Convert to factor and then to numeric
    gene_num = as.numeric(factor(gene_categories)))

# Sanity check
length(unique(unique_subset_goi4_MS$gene_categories)) #28
length(unique(goi_n_MS$gene_num)) #28

length(goi_n_MS$Sample_ID)
length(unique(goi_n_MS$gene_num))
length(unique(goi_n_MS$gene_categories))

goi_n_MS

# Adjust the data frame for plotting
goi_n_long_MS <- goi_n_MS %>%
  mutate(x = as.numeric(factor(Sample_ID, levels = unique(Sample_ID))),
         x = scales::rescale(x, to = range(sampleID_num)),
         value = tpm.x) %>%
  group_by(Genus, gene_num) %>%
  mutate(r = 2 * sqrt(sum(value) / 2 / pi)) %>%
  ungroup()

#update taxonomy from silva138
silva138_ID<-read.csv('~/database/silva/silva138_fasta_ids.csv',header=TRUE,sep=',')
silva138_ID<-subset(silva138_ID, Kingdom!='Eukaryota')

silva138_ID

goi_n_long_MS_taxa<-data.frame(goi_n_long_MS)

goi_n_long_MS_taxa2 <-goi_n_long_MS_taxa %>%
  left_join(silva138_ID, by = c("Genus"),relationship ="many-to-many")
  
goi_n_long_MS_taxa2 <-unique(goi_n_long_MS_taxa2)

goi_n_long_MS_taxa2

print(colnames(goi_n_long_MS_taxa2))

write.csv(goi_n_long_MS_taxa2[,c(1,2,3,10:14,4)],"metatrascriptomic_tpm_taxonomy.csv")

save.image()

write.csv(goi_n_long_MS_taxa2,'../plot/goi_n_long_MS_taxa2.csv')

goi_n_long_MS_taxa3<-read.csv('../plot/goi_n_long_MS_taxa2.csv',header=TRUE,sep=',')

# Sort the data by 'sort_col' in descending order (modify order as needed)
sorted_data <- arrange(goi_n_long_MS_taxa3, Phylum)

#svg("../plot/pie_MS3.svg",width=14,height=16)
    ggplot() +
  geom_scatterpie(aes(x = x, y = gene_num,fill=Phylum, r = .25),
                  data = sorted_data, cols = "Genus", long_format = TRUE) +

  #scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(8,"Set3"))(55))) +

  scale_x_continuous(
    breaks = unique(sorted_data$x),
    labels = unique(sorted_data$Sample_ID)) +
  scale_y_continuous(
    breaks = unique(sorted_data$gene_num),
    labels = unique(sorted_data$gene_categories)) +
  coord_equal() + theme_dpml() +

  theme(legend.position = "right",
    text = element_text(family = "Helvetica", size = 10))
#dev.off()



save.image()




