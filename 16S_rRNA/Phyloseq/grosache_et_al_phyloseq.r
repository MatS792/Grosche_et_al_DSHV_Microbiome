# The chemosynthetic microbiome of deep-sea hydrothermal vents across space and time

## Deepsea Hydrothermal vent  microbiome investigation: 16S rRNA based

#loading libraries
library(phyloseq)
library(ape)
library(microbiome)
library(microbiomeutilities)
library(vegan)
library(repr)
library(tidyverse)
library(ggthemes) # additional themes fro ggplot2
library(ggpubr)
library(ggpmisc) #to use stat_poly_eq
library(RColorBrewer) # nice color options
library(gridExtra) # gridding plots
library(viridis)
library(ggrepel)
library(rioja) # plotting poackages for tabular bubbleplots
library(reshape2) 
library(patchwork)
options(repr.plot.width=16, repr.plot.height=12)
set.seed(100000000)

# DeepseaMicrobiologyLab theme
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

#Uploading seqtab_nochim and taxa .rds file
seqtab_nochim<-readRDS("dada2_data/rds/seqtabnochim_run13dec23.rds")
taxa<-readRDS("dada2_data/rds/taxa_run13dec23.rds")
tree <- read.tree(file = "dada2_data/tree/alignment_13dec23.tree")

#Uploading env dataset
env_data<-read.csv("../dataset/Biofilms_metadata.csv",header=TRUE,sep=",")
row.names(env_data)<-env_data[,1]
head(env_data)

# Representing physico-chemical data

library(polygonPlot)

#uploading datasets
alv_polyplot_df<-read.csv("../dataset/polygonplot/alvinella_polyplot.csv",header=TRUE,sep=",")
muss_polyplot_df<-read.csv("../dataset/polygonplot/mussel_polyplot.csv",header=TRUE,sep=",")
rif_polyplot_df<-read.csv("../dataset/polygonplot/riftia_polyplot.csv",header=TRUE,sep=",")
fil_polyplot_df<-read.csv("../dataset/polygonplot/filaments.csv",header=TRUE,sep=",")
bas_polyplot_df<-read.csv("../dataset/polygonplot/bare_basalt.csv",header=TRUE,sep=",")
seawat_polyplot_df<-read.csv("../dataset/polygonplot/seawater_polyplot.csv",header=TRUE,sep=",")

alv_polygPlot<-polygonplot(alv_polyplot_df, shape = 3, fillcolor = "green", linecolor = "green3",labels_axis = c("Temperature (°C)", "H2S", "pH"), title="Alvinella")
alv_polygPlot

muss_polygPlot<-polygonplot(muss_polyplot_df, shape = 3, fillcolor = "seagreen1", linecolor = "seagreen3",labels_axis = c("Temperature (°C)", "H2S", "pH"), title="Mussel")
muss_polygPlot

rif_polygPlot<-polygonplot(rif_polyplot_df, shape = 3, fillcolor = "orange", linecolor = "orange3",labels_axis = c("Temperature (°C)", "H2S", "pH"), title="Rfitia")
rif_polygPlot

fil_polygPlot<-polygonplot(fil_polyplot_df, shape = 3, fillcolor = "gold", linecolor = "gold3",labels_axis = c("Temperature (°C)", "H2S", "pH"), title="Filaments")
fil_polygPlot

bas_polygPlot<-polygonplot(bas_polyplot_df, shape = 3, fillcolor = "grey", linecolor = "grey2",labels_axis = c("Temperature (°C)", "H2S", "pH"), title="Bare basalt")
bas_polygPlot

seawat_polygPlot<-polygonplot(seawat_polyplot_df, shape = 3, fillcolor = "purple", linecolor = "grey2",labels_axis = c("Temperature (°C)", "H2S", "pH"), title="Seawater")
seawat_polygPlot

ggarrange(alv_polygPlot,muss_polygPlot,rif_polygPlot,fil_polygPlot,bas_polygPlot,seawat_polygPlot,
          ncol=3,nrow=2) #using INKSCAPE to adjust the colors and other details

# Starting phyloseq analysis

# Phyloseq object
prok_data_raw<-phyloseq(sample_data(env_data),
                        otu_table(as.matrix(seqtab_nochim), taxa_are_rows=FALSE),
                        tax_table(as.matrix(taxa)),
                       phy_tree(tree))
prok_data_raw
sum(readcount(prok_data_raw))

# Cleaning up unwanted sequences from Eukarya, mitochrondria and chloroplast

prok_data_euk <- subset_taxa(prok_data_raw,  (Kingdom != "Eukaryota") | is.na(Kingdom))
prok_data_euk
sum(readcount(prok_data_euk))
1 - (sum(readcount(prok_data_euk))/sum(readcount(prok_data_raw))) # sanity check percentage of reads after Eukaryota removal

prok_data_chl <- subset_taxa(prok_data_euk, (Order!="Chloroplast") | is.na(Order))
prok_data_chl
sum(readcount(prok_data_chl))
1- (sum(readcount(prok_data_chl))/sum(readcount(prok_data_raw))) # sanity check percentage of reads after Chloroplast removal

prok_data_mit <- subset_taxa(prok_data_chl, (Family!="Mitochondria") | is.na(Family))
prok_data_mit
sum(readcount(prok_data_mit))
1- (sum(readcount(prok_data_mit))/sum(readcount(prok_data_raw))) # sanity check percentage of reads after Chloroplast removal

## Removing the known DNA Extraction contaminants from Sheik et al., 2018
# Assuming you are starting from a phyloseq object called prok_data_raw
prok_data_contam <- subset_taxa(prok_data_mit,  (Genus != "Afipia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Aquabacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Asticcacaulis") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Aurantimonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Beijerinckia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Bosea") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Bradyrhizobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Brevundimonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Caulobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Craurococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Devosia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Hoefleae") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Mesorhizobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Methylobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Novosphingobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Ochrobactrum") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Paracoccus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Pedomicrobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Phyllobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Rhizobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Sphingobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Sphingomonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Sphingopyxis") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Acidovorax") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Azoarcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Azospira") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Burkholderia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Comamonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Cupriavidus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Curvibacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Delftiae") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Duganella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Herbaspirillum") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Janthinobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Kingella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Leptothrix") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Limnobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Massilia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Methylophilus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Methyloversatilis") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Oxalobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Pelomonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Polaromonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Neisseria") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Ralstonia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Schlegelella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Sulfuritalea") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Undibacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Variovorax") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Acinetobactera") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Enhydrobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Enterobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Escherichia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Nevskia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Pasteurella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Pseudoxanthomonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Psychrobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Stenotrophomonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Xanthomonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Aeromicrobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Actinomyces") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Arthrobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Beutenbergia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Brevibacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Corynebacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Curtobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Dietzia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Janibacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Kocuria") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Microbacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Micrococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Microlunatus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Patulibacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Propionibacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Rhodococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Tsukamurella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Chryseobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Dyadobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Flavobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Hydrotalea") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Niastella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Parabacteroides") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Pedobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Prevotella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Wautersiella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Deinococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Abiotrophia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Bacillus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Brevibacillus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Brochothrix") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Facklamia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Olivibacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Lactobacillus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Paenibacillus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Ruminococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Staphylococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Streptococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Veillonella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Fusobacterium") | is.na(Genus))
prok_data_contam
message("reads saved from prok_data_raw after contaminants removal")
sum(readcount(prok_data_contam))/sum(readcount(prok_data_raw)) # reads saved from prok_data_raw after contaminants removal 
message("reads saved from prok_data_prune_mit after Eukaryota, Chloroplast, and Mitochondria removal")
sum(readcount(prok_data_contam))/sum(readcount(prok_data_mit)) # reads saved from prok_data_prune_mit after Eukaryota, Chloroplast, and Mitochondria removal

# Removing the potential human pathogens and contaminants
# This step needs to be evaluated with attention since many of these genera
# might be relevant in many environmental settings. Usually it is better to compare
# before-after removal to see what and how much you are removing. Feel free to
# experiment with the different groups and evaluate the results.

prok_data_humpath <- subset_taxa(prok_data_contam, (Genus != "Abiotrophia") |  is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Achromobacter") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Acinetobacter") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Actinobacillus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Arcanobacterium") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Babesia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Bifidobacterium") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Bartonella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Bordetella") |  is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Borrelia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Brodetella") |  is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Brucella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Capnocytophaga") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Chlamydia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Citrobacter") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Comamonas") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Corynebacterium_1") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Corynebacterium") | is.na(Genus)) # anche marino
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Coxiella") | is.na(Genus)) # anche marino
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Cronobacter") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Cutibacterium") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Dermatophilus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Ehrlichia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Order != "Enterobacteriales") | is.na(Order))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Enterococcus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Erysipelothrix") |  is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Escherichia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Escherichia/Shigella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Francisella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Gardnerella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Granulicatella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Haemophilus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Hafnia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Helicobacter") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Klebsiella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Kocuria") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Lactococcus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Lactobacillus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Lawsonia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Legionella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Leptospira") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Listeria") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Merkel_cell") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Micrococcus") | is.na(Genus)) # anche marino
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Morganella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Mycoplasma") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Neisseria") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Nocardia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Pasteurella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Plesiomonas") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Propionibacterium") | is.na(Genus))
#prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Proteus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Providencia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Pseudomonas") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Rhodococcus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Rickettsiae") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Roseomonas") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Rothia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Salmonella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Serratia") | is.na(Genus))
#prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Shewanella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Shigella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Sphaerophorus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Staphylococcus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Stenotrophomonas") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Streptococcus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Treponema") | is.na(Genus))
#prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Vibrio") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Yersinia") | is.na(Genus))
message("reads saved from prok_data_raw after contaminants removal")
sum(readcount(prok_data_humpath))/sum(readcount(prok_data_raw)) # reads saved from prok_data_raw after contaminants removal 
message("reads saved from prok_data_contam after DNA KIT EXTRAC. contaminants removal")
sum(readcount(prok_data_humpath))/sum(readcount(prok_data_contam)) # reads saved from prok_data_contam after DNA KIT EXTRAC. contaminants removal
prok_data_humpath
sum(readcount(prok_data_humpath))

(1 - (sum(readcount(prok_data_contam))/sum(readcount(prok_data_mit))))*100

cleaning_table1<-t(data.frame(
c((1 - (sum(readcount(prok_data_euk))/sum(readcount(prok_data_raw))))*100,'%'),
c((1 - (sum(readcount(prok_data_chl))/sum(readcount(prok_data_raw))))*100,'%'),
c((1 - (sum(readcount(prok_data_mit))/sum(readcount(prok_data_raw))))*100,'%'),
c((1 - (sum(readcount(prok_data_contam))/sum(readcount(prok_data_mit))))*100,'%'),
c((1 - (sum(readcount(prok_data_humpath))/sum(readcount(prok_data_contam))))*100,'%')    
    )
 )
rownames(cleaning_table1)<-c("Eukaryota","Chloroplast","Mitochondria","Kit contam.","Human contam.")
cleaning_table1 #tiny organization of the contamination results

#sanity check
prok_data2 = filter_taxa(prok_data_humpath, function(x) sum(x) > 0, TRUE)
prok_data2
sum(readcount(prok_data2))
100*(1-(sum(readcount(prok_data2))/sum(readcount(prok_data_raw))))

library(randomcoloR)
n <- 20
palette <- distinctColorPalette(n)
palette<-as.vector(palette)
scales::show_col(palette, labels=TRUE)

#svg("../results/plot/rarefaction_curves/rarefaction_curve.svg", width=12,height=12)
rarecurve(as.matrix(data.frame(otu_table(prok_data2))), col=palette, step=100 , lwd=3, ylab="ASVs", label=T)
#dev.off() #using INKSCAPE to adjust the colors and other details

# Alpha diversity by Bioregime

alpha_shannon<- 
plot_richness(prok_data2, x="Bioregime",measures=c("Shannon"), sortby = "Shannon") +
geom_boxplot(aes(fill=Bioregime)) +  ggtitle("Alpha diversity with Shannon index")+
scale_fill_brewer("Bioregime",palette = "Dark2",na.value="black",direction = -1) +
geom_jitter(size=3) +
theme_dpml() +
theme(axis.text.x = element_text(angle = -90))
alpha_shannon
#ggsave("../results/plot/alpha_diversity/alpha_diversity_shannon_bioregime.svg",width=16,height=12) #using INKSCAPE to adjust the colors and other details

alpha_observed<- 
plot_richness(prok_data2, x="Bioregime",measures=c("Observed"), sortby = "Observed") +
geom_boxplot(aes(fill=Bioregime)) +  ggtitle("Alpha diversity with Shannon index")+
scale_fill_brewer("Bioregime",palette = "Dark2",na.value="black",direction = -1) +    geom_jitter(size=3) +
ggtitle("Alpha diversity with Observed index")+
theme_dpml() +
theme(axis.text.x = element_text(angle = -90))
alpha_observed
#ggsave("../results/plot/alpha_diversity/alpha_diversity_observed_bioregime.svg",width=16,height=12) #using INKSCAPE to adjust the colors and other details

shannon_age_filt<-
plot_richness(subset_samples(prok_data2,! SampleID %in% c("SW1","SW2","NM2","NM3")),x="Age",measures=c("Shannon"), sortby = "Shannon") +
geom_boxplot(aes(fill=Age)) +  
geom_jitter(size=3) +
theme_dpml() +
theme(axis.text.x = element_text(angle = -90))
shannon_age_filt_df<-shannon_age_filt$data
shannon_age_filt
#ggsave("../results/plot/alpha_diversity/alpha_diversity_age_filt.svg",width=16,height=12) #using INKSCAPE to adjust the colors and other details

adonis_test_shannon_age_filt<-adonis2(shannon_age_filt_df$value~shannon_age_filt_df$Age,data=shannon_age_filt_df)
adonis_test_shannon_age_filt

#Alphaa diversity multi-index comparing
alpha_div_summary<-cbind(
    data.frame(sample_sums(prok_data2)),
   (estimate_richness(prok_data2, split = TRUE, measures = "Shannon")),
    (estimate_richness(prok_data2, split = TRUE, measures = "Chao1")),
    (estimate_richness(prok_data2, split = TRUE, measures = "Observed"))
                        )
alpha_div_summary <- alpha_div_summary[,-4]
alpha_div_summary<-cbind(alpha_div_summary,alpha_div_summary[,4]/alpha_div_summary[,3])
colnames(alpha_div_summary)<-c("Reads","Shannon","Chao1","Observed","Coverage")
alpha_div_summary

#write.csv(alpha_div_summary, "../results/plot/tables/alpha_div_summary.csv")

# Normalize the counts across the different samples by converting the abundance to relative abundance 
# and multiply by the median library size
prok_ndata <- transform_sample_counts(prok_data2, function(x) ((x / sum(x))*median(readcount(prok_data2))))
prok_ndata
sum(readcount(prok_ndata))

# Transform normalized abundance to relative abundances for plotting and some stats
prok_ra = transform_sample_counts(prok_ndata, function(x){x / sum(x)})
prok_ra

## Agglomerate at a specific taxonomic level at the Genus level
prok_ra_genus = tax_glom(prok_ra, "Genus", NArm = FALSE)
prok_ra_family = tax_glom(prok_ra, "Family", NArm = FALSE)
prok_ra_order = tax_glom(prok_ra, "Order", NArm = FALSE)
prok_ra_class = tax_glom(prok_ra, "Class", NArm = FALSE)
prok_ra_phyla = tax_glom(prok_ra, "Phylum", NArm = FALSE)

# Clustering and dendogram
# Compute distance matrix
dist.wunifrac <- UniFrac(prok_ndata, weighted=TRUE, normalized=TRUE)
# Clustering using Ward linkage
clustering <- hclust(dist.wunifrac, method = "ward.D2") #plot(clustering, hang = -1, cex = 0.6) 
# Convert clustering analysis to dendrogram and plot
hcd <- as.dendrogram(clustering)
sample_data <- as.data.frame(sample_data(prok_ndata))
sample_data2 <- data.frame(sample_data)
# Define nodePar
nodePar <- list(lab.cex = 1, pch = c(NA, 19), 
                cex = 3, col = "blue")

#svg("../results/plot/dendogram/dendrogram.svg")
  par(mar = c(2,2,2,10))
  hcd %>% 
      plot(main="16S Metabarcoding Weighted Unifrac: Ward Linkage ", horiz = TRUE, nodePar=nodePar)
 #dev.off() #using INKSCAPE to adjust the colors and other details

library(fantaxtic)
library(MicEco)
library(see) #viasualization and half-violin plots

# Using only taxa with a relative abundance > 1%
pursed_pseq.phylum_perc1 <- ps_prune(prok_ra_phyla, min.abundance = 0.01)
pursed_pseq.class_perc1 <- ps_prune(prok_ra_class, min.abundance = 0.01)
pursed_pseq.gen_perc1 <- ps_prune(prok_ra_genus, min.abundance = 0.01)

#Other and NAs are mixed up - Manual correction
tax_table(pursed_pseq.phylum_perc1)
tax_table(pursed_pseq.phylum_perc1)[3,2]<-"Kingdom Bacteria - Phylum NA"
tax_table(pursed_pseq.phylum_perc1) #sanity check

tax_table(pursed_pseq.class_perc1)
tax_table(pursed_pseq.class_perc1)[3,3]<-"Kingdom Bacteria - Class NA"
tax_table(pursed_pseq.class_perc1) #sanity check

tax_table(pursed_pseq.gen_perc1)
tax_table(pursed_pseq.gen_perc1)[1,6] <- "Family Bacteroidetes BD2-2 - Genus NA"
tax_table(pursed_pseq.gen_perc1)[3,6] <- "Kingdom Bacteria - Genus NA"
tax_table(pursed_pseq.gen_perc1)[5,6] <- "Family Thiotrichaceae - Genus NA"
tax_table(pursed_pseq.gen_perc1)[7,6] <- "Class Gammaproteobacteria - Genus NA"
tax_table(pursed_pseq.gen_perc1)[11,6] <- "Family Arcobacteraceae - Genus NA"
tax_table(pursed_pseq.gen_perc1) #sanity check

#Barplot: Phylum
phylum_1perc<-plot_bar(pursed_pseq.phylum_perc1, x="SampleID", fill ="Phylum") +
  ggtitle("Phylum Abundance of >1% ASVs") +
  theme(axis.text.x = element_text(face="bold",size=10))+
  theme_dpml() + theme(legend.position = "bottom") +
  scale_y_continuous(labels=c(0,25,50,75,100)) + scale_fill_brewer("Phylum", palette = "PuBuGn",na.value="black") +
  coord_flip() 

phylum_1perc$data$SampleID <- factor(phylum_1perc$data$SampleID, levels = c('RIF1','MUS1','ALV4','MUS2','MUS4','MUS3','RIF5','RIF4','NM3',
                                                                        'NM2','RIF3','ALV2','ALV3','ALV1','RIF6','NM1','RIF2','SW2','SW1'))
phylum_1perc + scale_x_discrete(limits = rev(levels(phylum_1perc$data$SampleID)))
#ggsave("../results/plot/taxonomy/last/phylum_1perc.svg", width=16,height=12) #using INKSCAPE to adjust the colors and other details

#Barplot: Class
class_1perc<-plot_bar(pursed_pseq.class_perc1, x="SampleID", fill ="Class") +
  ggtitle("Class Abundance of >1% ASVs") +
  theme(axis.text.x = element_text(face="bold",size=10))+
  theme_dpml() + theme(legend.position = "bottom") +
  scale_y_continuous(labels=c(0,25,50,75,100)) + scale_fill_brewer("Class", palette = "Spectral",na.value="black") +
  coord_flip() 

class_1perc$data$SampleID <- factor(class_1perc$data$SampleID, levels = c('RIF1','MUS1','ALV4','MUS2','MUS4','MUS3','RIF5','RIF4','NM3',
                                                                        'NM2','RIF3','ALV2','ALV3','ALV1','RIF6','NM1','RIF2','SW2','SW1'))
class_1perc + scale_x_discrete(limits = rev(levels(class_1perc$data$SampleID)))
#ggsave("../results/plot/taxonomy/last/class_1perc.svg", width=16,height=12) #using INKSCAPE to adjust the colors and other details

#Barplot: Genus
genus_1perc<-plot_bar(pursed_pseq.gen_perc1, x="SampleID", fill ="Genus") +
  ggtitle("Genus Abundance of >1% ASVs") +
  theme(axis.text.x = element_text(face="bold",size=10))+
  theme_dpml() + theme(legend.position = "bottom") +
  scale_y_continuous(labels=c(0,25,50,75,100)) + 
#scale_fill_brewer("Genus", palette = "RdBu",na.value="black") +
scale_fill_manual(na.value="black","Genus",values = colorRampPalette(brewer.pal(11, "RdBu"))(16)) + 
coord_flip() 

genus_1perc$data$SampleID <- factor(genus_1perc$data$SampleID, levels = c('RIF1','MUS1','ALV4','MUS2','MUS4','MUS3','RIF5','RIF4','NM3',
                                                                        'NM2','RIF3','ALV2','ALV3','ALV1','RIF6','NM1','RIF2','SW2','SW1'))
genus_1perc + scale_x_discrete(limits = rev(levels(genus_1perc$data$SampleID))) 
#ggsave("../results/plot/taxonomy/last/genus_1perc.svg", width=16,height=12) #using INKSCAPE to adjust the colors and other details

#Saving a .csv with combined taxonomy and ASVs abundances for the selected taxonomic level
summary_object <-prok_ra
summary_otu <- as.matrix(t(otu_table(summary_object)))
summary_tax <- as.matrix(tax_table(summary_object))
summary_comb <- cbind(summary_otu, summary_tax)
write.csv(summary_comb, "../results/summary/summary_prok_ra.csv")

#PCoA
.JACCARD W
.VECTOR FITTING
.ORDISURF
.ORDINATION AS TAXA
.ORDINATIONS (SAMPLES-TAXA) OVERLAP with VECTOR FITTING

prok_ndata_noSW<-subset_samples(prok_ndata,Substrate!='Seawater') #Seawater samples will be excluded from the following ordination

prok_dist_wjac_noSW <- phyloseq::distance(prok_ndata_noSW, method = "jaccard")

# https://rdrr.io/bioc/phyloseq/man/distance.html
prok_mds_jw_noSW <- ordinate(prok_ndata_noSW,prok_dist_wjac_noSW, method = "MDS")

#comparing weighted vs unweighted 
plot_ordination(prok_ndata_noSW,prok_mds_jw_noSW) +
geom_text(aes(label= SampleID), size=5, hjust=-0.5,vjust=1) +
geom_point(aes(fill=Bioregime),shape=21,size=6,color="black",stroke=0.3) +
scale_fill_brewer("Bioregime",palette = "Dark2",na.value="black",direction = -1) +
scale_shape_manual(values=c(21:25)) +
theme_dpml()

mds_samples_noSW<-(plot_ordination(prok_ndata_noSW, prok_mds_jw_noSW, type="samples"))$data
mds_taxa_noSW<-(plot_ordination(prok_ndata_noSW, prok_mds_jw_noSW, type="taxa"))$data

#Vector fitting using env variables
env_jaccardW_noSW <-envfit(mds_samples_noSW[,c(1:2)], mds_samples_noSW[,c(8:9)], perm = 9999, na.rm = T)
env_jaccardW_noSW

#Vector fitting
env.scores_noSW <- as.data.frame(scores(env_jaccardW_noSW, display = "vectors")) #extracts relevant scores from envifit
env.scores_noSW <- cbind(env.scores_noSW, env.variables = rownames(env.scores_noSW)) #and then gives them their names

env.scores_noSW <- cbind(env.scores_noSW, pval = env_jaccardW_noSW$vectors$pvals) # add pvalues to dataframe
sig.env.scrs_noSW <- subset(env.scores_noSW, pval<=0.1) #subset data to show variables significant at 0.05

head(sig.env.scrs_noSW)

#ordisurf

env_ordsurf_temp_noSW<-ordisurf(mds_samples_noSW[,1:2] ~ Temperature, mds_samples_noSW[,c(1,2,8)], family=gaussian, plot = T,knots = 2)
env_ordsurf_h2s_noSW<-ordisurf(mds_samples_noSW[,1:2] ~ Hydrogen_sulfide, mds_samples_noSW[,c(1,2,9)], family=gaussian, plot = T,knots = 2)

summary(env_ordsurf_temp_noSW)

env_ordsurf_temp_noSW
env_ordsurf_h2s_noSW

library(BiodiversityR) # also loads vegan
library(readxl)
library(ggsci)
library(ggrepel)
library(ggforce)

ordisurf_temp.grid_noSW <- ordisurfgrid.long(env_ordsurf_temp_noSW)
ordisurf_h2s.grid_noSW <- ordisurfgrid.long(env_ordsurf_h2s_noSW)

print(colnames(mds_samples_noSW))

ggplot(data = mds_samples_noSW, aes(x = Axis.1, y = Axis.2)) +

geom_contour(data=ordisurf_temp.grid_noSW, 
                        aes(x=x, y=y, z=z,colour=factor(after_stat(level))),linewidth=1) +
scale_colour_viridis_d(option = 'plasma') +

geom_text(aes(label= SampleID), size=8, hjust=1,vjust=-1) +
geom_point(size=7,color="black",stroke=0.3,aes(fill=Bioregime,shape=Substrate)) +
scale_fill_brewer("Bioregime",palette = "Dark2",na.value="black",direction = -1) +
scale_shape_manual(values = c(22:25)) +
theme_dpml() 
#ggsave("../../13Dec2023/results/plot/ordinations/mds_JW_noSW_ordisurf_temp.svg", width=14, height=14) #using INKSCAPE TO ADJUST COLORS AND DETAILS


ggplot(data = mds_samples_noSW, aes(x = Axis.1, y = Axis.2)) +

geom_contour(data=ordisurf_h2s.grid_noSW, 
                        aes(x=x, y=y, z=z,colour=factor(after_stat(level))),linewidth=1) +
scale_colour_viridis_d(option = 'viridis',) +

geom_text(aes(label= SampleID), size=8, hjust=1,vjust=-1) +
geom_point(size=7,color="black",stroke=0.3,aes(fill=Bioregime,shape=Substrate)) +
scale_fill_brewer("Bioregime",palette = "Dark2",na.value="black",direction = -1) +

scale_shape_manual(values = c(22:25)) +
theme_dpml() 
#ggsave("../../13Dec2023/results/plot/ordinations/mds_JW_noSW_ordisurf_h2s.svg", width=14, height=14)  #using INKSCAPE TO ADJUST COLORS AND DETAILS

library(pairwiseAdonis)

# Adonis
adonis2(prok_dist_wjac_noSW~Bioregime+Substrate,permutations=999,
    data = mds_samples_noSW)

adonis2(prok_dist_wjac_noSW~Bioregime+Substrate+Site,permutations=999,
    data = mds_samples_noSW)

pairwise.adonis(prok_dist_wjac_noSW, mds_samples_noSW$Bioregime, reduce='Mussels')

###TAXA###
mds_taxa_noSW_TAXA<-ggplot(data = mds_taxa_noSW, aes(x = Axis.1, y = Axis.2)) +

geom_point(shape=21,size=8,color="black",stroke=0.25,fill="grey",alpha=0,
           data=subset(mds_taxa_noSW,! Class %in% c("Campylobacteria","Gammaproteobacteria"))) +

geom_point(shape=21,size=8,color="black",stroke=0.25,fill="orange",alpha=0.5,
           data=subset(mds_taxa_noSW,Class %in% c("Campylobacteria"))) +

geom_point(shape=21,size=8,color="black",stroke=0.25,fill="navy",alpha=0.5,
           data=subset(mds_taxa_noSW,Class %in% c("Gammaproteobacteria"))) + 
theme_dpml() + theme(legend.position = "right") + xlab('Axes.1 (13.2%)') + ylab('Axes.2 (10.7%)')

mds_taxa_noSW_TAXA + geom_segment(data = sig.env.scrs_noSW, aes(x = 0, y = 0, xend = Axis.1/2, yend = Axis.2/2), 
        linewidth =.8, alpha = 0.8, colour = "grey") + 
geom_text(data = sig.env.scrs_noSW, aes(x = Axis.1/2, y = Axis.2/2), colour = "grey30", 
       fontface = "bold", size=8, label = row.names(sig.env.scrs_noSW))
 

#ggsave("../results/plot/ordinations/mds_taxa_vectors.svg",width=16,height=12) #using INKSCAPE to adjust the colors and other details

#Saving a .csv with combined taxonomy and ASVs abundances for the selected taxonomic level
summary_object <- prok_ra_class
summary_otu <- as.matrix(t(otu_table(summary_object)))
summary_tax <- as.matrix(tax_table(summary_object))
summary_comb <- cbind(summary_otu, summary_tax)
gamma_epsilon_df<-subset(summary_comb, summary_comb[,22] %in% c('Gammaproteobacteria'))
gamma_epsilon_df<-data.frame(gamma_epsilon_df[,-c(20:21,23:25)])
colnames(gamma_epsilon_df)<-gamma_epsilon_df[20,]
gamma_epsilon_df<-rownames_to_column(data.frame(gamma_epsilon_df),var='SampleID')
gamma_epsilon_df<-gamma_epsilon_df[-20,]
gamma_epsilon_df

#Saving a .csv with combined taxonomy and ASVs abundances for the selected taxonomic level
summary_object <- prok_ra_class
summary_otu <- as.matrix(t(otu_table(summary_object)))
summary_tax <- as.matrix(tax_table(summary_object))
summary_comb <- cbind(summary_otu, summary_tax)
campyl_epsilon_df<-subset(summary_comb, summary_comb[,22] %in% c('Campylobacteria'))
campyl_epsilon_df<-data.frame(campyl_epsilon_df[,-c(20:21,23:25)])
colnames(campyl_epsilon_df)<-campyl_epsilon_df[20,]
campyl_epsilon_df<-rownames_to_column(data.frame(campyl_epsilon_df),var='SampleID')
campyl_epsilon_df<-campyl_epsilon_df[-20,]
campyl_epsilon_df

epsilon_df<-cbind(campyl_epsilon_df,gamma_epsilon_df)
epsilon_df

epsilon_df2<-list(mds_samples_noSW,epsilon_df[,c(1,2,4)]) %>%
    reduce(left_join, by="SampleID")
epsilon_df2

epsilon_vs_gamma<-read.csv("csv/prok_ra_campy_gamma.csv",header=TRUE,sep=',')
epsilon_vs_gamma$ra<-epsilon_vs_gamma$ra*100
epsilon_vs_gamma$ra

p_cor<- ggplot(data = epsilon_vs_gamma,aes(x=age, y=ra)) +
geom_boxplot(aes(fill=group),position = position_dodge(width = 0.8),alpha=0.5) +
scale_fill_manual(values=c("orange","navy")) + labs(x="",y="Relative Abundance (%)",fill='Class') +
theme_dpml()

df_p_val_new <- rstatix::wilcox_test(subset(epsilon_vs_gamma,age=='Newly-formed'),ra~group) %>% 
	rstatix::add_xy_position()
df_p_val_est <- rstatix::wilcox_test(subset(epsilon_vs_gamma,age=='Established'),ra~group) %>% 
	rstatix::add_xy_position()

p_cor <-p_cor + stat_pvalue_manual(df_p_val_new, label = "p", tip.length = 1, hide.ns = FALSE) + 
                stat_pvalue_manual(df_p_val_est, label = "p", tip.length = 1, hide.ns = FALSE)

p_cor
#ggsave("../results/plot/taxonomy/boxplot_campylobacteria_gammaprotebacteria.svg",width=16,height=14) #using INKSCAPE to adjust the colors and other details

#Venn diagrams

library(VennDiagram)
library(gplots)
library(venn)
library(venneuler)
library(UpSetR)
library(MicrobiotaProcess)
library(ggvenn)

prok_ra_venn_list<-get_vennlist(subset_samples(prok_ra, ! Bioregime %in% c('Bare Basalt','Fluids','Filaments')),
                                      factorNames="Bioregime")

ggvenn(
  prok_ra_venn_list, 
  fill_color = c("#990f0fff", "#260f99ff", "#6b990fff"),
  stroke_size = 0.5, set_name_size = 10, text_size=8
  ) + 
theme_dpml()
#ggsave("../results/plot/venn_diagram/venn_bioregimes.svg", width=14,height=16) #using INKSCAPE to adjust the colors and other details

#sanity check
message('% unique ASV Riftia')
round((1703/(218+217+298+1703+405+2+281))*100)
message('% unique ASV Mussels')
round((405/(218+217+298+1703+405+2+281))*100)
message('% unique ASV Alvinella')
round((281/(218+217+298+1703+405+2+281))*100)

#Bioregime specifc

#Riftia

prok_ra_rftia_list<-get_vennlist(subset_samples(prok_ra, Bioregime %in% 'Riftia'),
                                      factorNames="Age")

ggvenn(
  prok_ra_rftia_list, 
  fill_color = c("#990f0ffd", "#260f99ff"),
  stroke_size = 0.5, set_name_size = 10, text_size=8
  ) + theme_dpml()
#ggsave("../results/plot/venn_diagram/venn_riftia_age.svg", width=14,height=16) #using INKSCAPE to adjust the colors and other details

#sanity check
message('established-Newly-formed shared %')
round((565/(565+1583+288))*100)
message('established unique %')
round((1583/(565+1583+288))*100)
message('newly-formed unique %')
round((288/(565+1583+288))*100)

#Mussel

prok_ra_mussel_list<-get_vennlist(subset_samples(prok_ra, Bioregime %in% 'Mussels'),
                                      factorNames="Age")

ggvenn(
  prok_ra_mussel_list, 
  fill_color = c("#990f0ffd", "#260f99ff"),
  stroke_size = 0.5, set_name_size = 10, text_size=8
  ) + theme_dpml()
#ggsave("../results/plot/venn_diagram/venn_mussels_age.svg", width=14,height=16) #using INKSCAPE to adjust the colors and other details

message('established-Newly-formed shared %')
round(100*(171/(534+171+218)))
message('established unique %')
round(100*(534/(534+171+218)))
message('newly-formed unique %')
round(100*(218/(534+171+218)))

#Alvinella

prok_ra_alvinella_list<-get_vennlist(subset_samples(prok_ra, Bioregime %in% 'Alvinella'),
                                      factorNames="Age")

ggvenn(
  prok_ra_alvinella_list, 
  fill_color = c("#990f0ffd", "#260f99ff"),
  stroke_size = 0.5, set_name_size = 10, text_size=8
  ) + theme_dpml()
#ggsave("../results/plot/venn_diagram/venn_alvinella_age.svg", width=14,height=16) #using INKSCAPE to adjust the colors and other details

message('established-Newly-formed shared %')
100*(221/(287+221+210))
message('established unique %')
100*(210/(287+221+210))
message('newly-formed unique %')
100*(287/(287+221+210))

#Donuts plots

library(wesanderson)

riftia_donuts<-read.csv('csv/donuts/1percent/riftia_est_unique_genus_1perc_labels.csv',header=TRUE, sep=',')

mycolors <-wes_palette("Rushmore1", 11, type = "continuous")

ggplot(riftia_donuts, aes(x = hsize, y = Established, fill = Label)) +
geom_col() + 
scale_fill_brewer("Riftia Established Biofilms Unique Taxa", palette = "PRGn",na.value="black") +
coord_polar(theta = "y") +
xlim(c(0.2, hsize + 0.5)) + 
xlab('') + ylab('') +
theme_dpml() #using INKSCAPE to adjust the colors and other details

#Mussel

mussel_donuts<-read.csv('csv/donuts/1percent/mussel_est_unique_genus_1perc_labels.csv',header=TRUE, sep=',')

mycolors <-wes_palette("Royal1", 11, type = "continuous")

ggplot(mussel_donuts, aes(x = hsize, y = Established, fill = Label)) +
geom_col() + 
geom_col() + scale_fill_brewer("Mussels Established Biofilms Unique Taxa", palette = "RdBu",na.value="black") +
coord_polar(theta = "y") +
xlim(c(0.2, hsize + 0.5)) + 
xlab('') + ylab('') +
theme_dpml() #using INKSCAPE to adjust the colors and other details

#Alvinella

alvinella_donuts<-read.csv('csv/donuts/1percent/alvinella_est_unique_genus_1perc_labels.csv',header=TRUE, sep=',')

colourCount =length(unique(alvinella_donuts$Established))
colourCount

ggplot(alvinella_donuts, aes(x = hsize, y = Established, fill = Label)) +
geom_col() + 
scale_fill_manual("Alvinella Established Biofilms Unique Taxa",values = colorRampPalette(brewer.pal(11, "BrBG"))(colourCount)) + 
coord_polar(theta = "y") +
xlim(c(0.2, hsize + 0.5)) + 
xlab('') + ylab('') +
theme_dpml() #using INKSCAPE to adjust the colors and other details

save.image()
