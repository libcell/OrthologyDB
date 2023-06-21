################################################################################
### Exploring consistency between gene tree, protein tree and species tree   ###
### Author: Bo Li, et al., June 15, 2023                                     ###
################################################################################

### ****************************************************************************
### code chunk number 01: Selecting the model organisms.
### ****************************************************************************

### ------------------------------------------------------------------------ ###
### Step-01. Selecting sixteen most popular model organisms from GEO database.

### Url for information on 16 model organisms. 
### https://www.ncbi.nlm.nih.gov/geo/summary/?type=tax

### Checking the information of 16 model organisms in NCBI genome. 

### End of Step-01.
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-02. Obtaining orthologs for these 16 model organisms via KEGGREST.

### 1) Reading the information on 16 model organism. 

species_list <- openxlsx::read.xlsx("./data/raw/Selected sixteen organims.xlsx", 
                                    sheet = 1)

species_abbr <- species_list$Organism

### 2) Counting the information on CDSs of 16 model organisms individuals. 

library(KEGGREST)

species_gene <- NULL

Gene_list <- list()

i <- 0

for (s in species_abbr) {
  
  i <- i + 1
  
  ## returns the entire list of human genes
  genes <- keggList(s) 
  
  genes <- genes[genes == "CDS"]
  
  species_gene <- c(species_gene, length(genes))
  
  Gene_list[[i]] <- genes
  
  print(s)
  
  print(species_gene)
  
  Sys.sleep(10)
}

names(species_gene) <- species_abbr

names(Gene_list) <- species_abbr

saveRDS(species_gene, "./data/species_gene.rds")
saveRDS(Gene_list, "./data/Gene_list.rds")






