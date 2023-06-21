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

### 3) Converting gene id to KO id for 16 model organisms. 

species_gene <- readRDS("species_gene.rds")
Gene_list <- readRDS("Gene_list.rds")

Gene_list_ko <- as.list(rep(NA, 16))

len <- length(Gene_list_ko) <- length(Gene_list)

for (s in 1:len) {
  
  print("#####################################################################")
  print(s)
  
  gl <- names(Gene_list[[s]])
  
  for (g in gl) {
    
    repeat {
      query <- try(names(keggGet(g)[[1]]$ORTHOLOGY), silent = TRUE)
      if (class(query) != "try-error") {
        break
      }
    }
    
    Gene_list_ko[[s]] <- c(Gene_list_ko[[s]], query)
    
    print(query)
    
    # print(Gene_list_ko[[s]])
    
    Sys.sleep(0.01)
    
  }
  
  Gene_list_ko[[s]] <- unique(Gene_list_ko[[s]])
  
  Gene_list_ko[[s]] <- as.vector(na.omit(Gene_list_ko[[s]]))
  
  print("#####################################################################")
  
}

names(Gene_list_ko) <- names(Gene_list)

for (i in 1:16) {
  
  Gene_list_ko[[i]] <- unique(Gene_list_ko[[i]])
  
}

for (i in 1:16) {
  
  Gene_list_ko[[i]] <- Gene_list_ko[[i]][-1]
  
}

genes <- Reduce(intersect, Gene_list_ko)

orth_16_species <- list(Geneset = Gene_list_ko, 
                        Comgene = genes)

saveRDS(orth_16_species, file = "./data/orth_16_species.rds")

gene_stats <- keggList("hsa")

table(gene_stats)

### 4) Statistics for all genes in individual species.

spe <- names(Gene_list_ko)
ko_genes <- matrix(NA, nrow = 176, ncol = 16)
for (rw in 1:length(genes)) {
  txt_1 <- paste("row number is", rw, sep = " ")
  print("#####################################################################")
  print(txt_1)
  spe_gene <- keggGet(genes[rw])[[1]]$GENES
  for (cl in 1:16) {
    txt_2 <- paste("column number is", cl, sep = " ")
    print(txt_2)
    text <- paste0("^", toupper(spe[cl]), ":")
    loc <- grep(text, spe_gene)
    gs <- spe_gene[loc]
    gs <- strsplit(gs, " ")[[1]]
    gs <- paste0(tolower(gs[1]), gs[2])vv
    gs <- strsplit(gs, "\\(")[[1]][1]
    ko_genes[rw, cl] <- gs
  }
}
ko_genes <- as.data.frame(ko_genes)

rownames(ko_genes) <- genes

colnames(ko_genes) <- names(Gene_list_ko)

DT::datatable(ko_genes)

library("Biostrings")

all.ntseq <- as.list(rep(NA, 176))
all.aaseq <- as.list(rep(NA, 176))
names(all.ntseq) <- names(all.aaseq) <- genes

for (r in 1:176) {
  print(r)
  ntseq <- DNAStringSet(NULL)
  aaseq <- AAStringSet(NULL)
  for (c in 1:16) {
    xl_nt <- keggGet(as.character(ko_genes[r, c]), "ntseq")
    xl_aa <- keggGet(as.character(ko_genes[r, c]), "aaseq")
    ntseq <- c(ntseq, xl_nt)
    aaseq <- c(aaseq, xl_aa)
  }
  all.ntseq[[r]] <- ntseq
  all.aaseq[[r]] <- aaseq
  Sys.sleep(5)
}

saveRDS(all.ntseq, file = "./data/all.ntseq.rds")
saveRDS(all.aaseq, file = "./data/all.aaseq.rds")

### 5) Statistics for all genes in individual species.







