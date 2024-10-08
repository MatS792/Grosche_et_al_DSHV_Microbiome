# Using Metaxa2 to investigate the taxonomic content
## from "Metagenomic Annotation Workshop" 
### https://metagenomics-workshop.readthedocs.io/en/2014-5/annotation/metaxa2.html

library(tidyverse)
library(viridis)
library(RColorBrewer)
library(phyloseq)
library(ggpubr)
options(repr.plot.width=16, repr.plot.height=12)

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

count_table<-read.csv("../16S_tables/count_table.csv", sep=',',header=TRUE,row.names = 1)
tax_table<-read.csv("../16S_tables/taxa_table.csv", sep=',',header=TRUE,row.names = 1)

head(count_table)
head(tax_table)

count_table$cv88_tpm<-count_table$CV88_Count/2266097 *1000000 #16S rRNA counts for the different genus per million reads
count_table$cv88_ra<-count_table$CV88_Count/sum(count_table$CV88_Count)
head(count_table)

count_table$lvp4_tpm<-count_table$LVP4_Count/17299 *1000000 #16S rRNA counts for the different genus per million reads
count_table$lvp4_ra<-count_table$LVP4_Count/sum(count_table$LVP4_Count)
head(count_table)

count_table_ra<-as.matrix(count_table[,c(4,6)])
tax_table<-as.matrix(tax_table)

sam_data<-read.csv("../16S_tables/sample_data.csv",header=TRUE,sep=",",row.names = 1)
sam_data

physeq.obj<-phyloseq(sample_data(sam_data),
                     otu_table(count_table_ra,taxa_are_rows = TRUE),
                     tax_table(tax_table))
physeq.obj

physeq.obj_class = tax_glom(physeq.obj, "Class", NArm = FALSE)

physeq.obj_class

#Saving a .csv with combined taxonomy and ASVs abundances for the selected taxonomic level
summary_object <-physeq.obj
summary_otu <- as.matrix((otu_table(summary_object)))
summary_tax <- as.matrix(tax_table(summary_object))
summary_comb <- cbind(summary_otu, summary_tax)
write.csv(summary_comb, "../summary_prok_ra_genus_metatrascript.csv")

library(fantaxtic)
library(MicEco)
library(see) #viasualization and half-violin plots

pseq.class_perc05 <- ps_prune(physeq.obj_class, min.abundance = 0.005)
pseq.class_perc05

metatrascr_class_bp<-plot_bar(pseq.class_perc05, fill ="Class") +
  ggtitle("") +
#  scale_fill_brewer("Genus", palette = "PiYG",na.value="black") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "BrBG"))(7)) +
  scale_x_discrete(labels=c('RIF3', 'LVP4')) +
  labs(y="Relative Abundance (%)", x='') +
  theme_dpml()  + theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0, size=20)) +  
  scale_y_continuous(labels=c(0,25,50,75,100)) + theme(aspect.ratio = 2.5)
metatrascr_class_bp

#ggsave('../metatrascr_class_barplot.svg',width=10,height=14)

pseq.gen_perc05 <- ps_prune(physeq.obj, min.abundance = 0.005)
pseq.gen_perc05

metatrascr_genus_bp<-plot_bar(pseq.gen_perc05, fill ="Genus") +
  ggtitle("") +
#  scale_fill_brewer("Genus", palette = "PiYG",na.value="black") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "PiYG"))(16)) +
  scale_x_discrete(labels=c('RIF3', 'LVP4')) +
  labs(y="Relative Abundance (%)", x='') +
  theme_dpml()  + theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0, size=20)) +
  scale_y_continuous(labels=c(0,25,50,75,100)) + theme(aspect.ratio = 2.5)
metatrascr_genus_bp

#ggsave('../metatrascr_genus_barplot.svg',width=10,height=14)



save.image()
