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

### 5) Building UPGMA tree using DNA & protein sequences for orthologs.

################################################################################
### Start here!!!

all.ntseq <- readRDS("./data/all.ntseq.rds")
all.aaseq <- readRDS("./data/all.aaseq.rds")

organism <- species_abbr


ko <- names(all.ntseq)

coef.mat <- NULL

for (k in ko) {
  
  # --------- Construct evolution tree using DNA sequence ------------------ ###
  
  dna <- all.ntseq[[k]]
  
  # sequence alignment
  names(dna) <- species_abbr
  
  library(msa)
  align.seq <- msa(dna,
                   method = "ClustalW")
  library(ggmsa)
  consv.motif <- msaConsensusSequence(align.seq)
  
  library(bios2mds)
  DNA <- msaConvert(align.seq,
                    type = "bios2mds::align")
  export.fasta(DNA,
               outfile = "./output/sequences/ntK00858.fas",
               ncol = 60,
               open = "w")
  #
  library(adegenet) 
  dna <- fasta2DNAbin(file = "./output/sequences/ntK00858.fas")
  
  library(ape)
  D <- dist.dna(dna, model = "TN93") #替代模型用TN93
  #UPGMA tree
  D[is.na(D)] <- 0
  D[is.nan(D)] <- 0
  sum(is.infinite(D))  # THIS SHOULD BE 0
  nth_cluster <- hclust(D, method = "average") 
  
  g.tree.name <- paste("Gene tree constructed by DNA sequences of", 
                       k, sep = " ")
  g.tree.file <- paste0("./output/plots/gene_tree/Gene_tree_of_", k, ".pdf")
  
  pdf(g.tree.file)
  plot(nth_cluster, main = g.tree.name)
  dev.off()
  
  
  # -------- Construct evolution protein using protein sequence ------------ ###
  
  aa <- all.aaseq[[k]]
  # sequence alignment
  names(aa) <- species_abbr
  library(msa)
  align.seq <- msa(aa, method = "Muscle")
  library(ggmsa)
  consv.motif <- msaConsensusSequence(align.seq, type = "upperlower")
  
  AA <- msaConvert(align.seq,
                   type = "bios2mds::align")
  library(bios2mds)
  AA <- export.fasta(AA,
                     outfile = "./output/sequences/aaK00858.fas",
                     ncol = 60,
                     open = "w")
  library(seqinr)
  aa <- read.alignment("./output/sequences/aaK00858.fas", format = "fasta")
  # Computing the distance between DNA sequences.
  Da <- dist.alignment(aa)
  #UPGMA tree
  Da[is.na(Da)] <- 0
  Da[is.nan(Da)] <- 0
  sum(is.infinite(Da))  # THIS SHOULD BE 0
  # method = average is used for UPGMA
  aah_cluster <- hclust(Da, method = "average") 
  
  p.tree.name <- paste("Protein tree constructed by AA sequences of", 
                       k, sep = " ")
  p.tree.file <- paste0("./output/plots/protein_tree/Protein_tree_of_", 
                        k, ".pdf")
  
  pdf(p.tree.file)
  plot(aah_cluster, main = p.tree.name)
  dev.off()
  
  
  ### 6) Obtaining the taxa tree using taxize package.
  
  ### 7) Comparison between gene tree, protein tree and species tree.
  
  ##------------ Dendrograms comparison ------------
  
  library(dendextend)
  library(treeio)
  library(ape)
  # load tree
  tree.species <- ape::read.nexus("./output/objects/tree16abbr.nexus")
  tree.aaseq <- as.phylo(aah_cluster)
  tree.ntseq <- as.phylo(nth_cluster)
  
  # Create two dendrograms
  dend.nt <- as.dendrogram (tree.ntseq)
  dend.aa <- as.dendrogram (tree.aaseq)
  dend.sp <- as.dendrogram (tree.species)
  
  # Create a list to hold dendrograms
  
  dend_list <- dendlist(dend.nt, dend.sp)
  r.g.s <- round(entanglement(dend_list), digits = 3)
  tanglegram(dend.nt, dend.sp,  
             main = paste("entanglement", 
                          r.g.s, 
                          sep = "="), 
             sub = "L: Gene; R: Species") 
  
  nt.tangle.sp <- paste("Tanglegram of", 
                       k, sep = " ")
  nt.tangle.sp.file <- paste0("./output/plots/gene_species/Tanglegram of_", 
                        k, ".pdf")
  
  pdf(nt.tangle.sp.file, width = 9, height = 6)
  tanglegram(dend.nt, dend.sp,  
             main = paste("entanglement", 
                          r.g.s, 
                          sep = "="), 
             sub = "L: Gene; R: Species") 
  dev.off()
  
  dend_list <- dendlist(dend.aa, dend.sp)
  r.p.s <- round(entanglement(dend_list), digits = 3)
  tanglegram(dend.aa, dend.sp, 
             main = paste("entanglement", 
                          r.p.s, 
                          sep = "="), 
             sub = "L: Protein; R: Species") 
  
  aa.tangle.sp <- paste("Tanglegram of", 
                        k, sep = " ")
  aa.tangle.sp.file <- paste0("./output/plots/protein_species/Tanglegram of_", 
                              k, ".pdf")
  
  pdf(aa.tangle.sp.file, width = 9, height = 6)
  tanglegram(dend.aa, dend.sp, 
             main = paste("entanglement", 
                          r.p.s, 
                          sep = "="), 
             sub = "L: Protein; R: Species") 
  dev.off()
  
  
  dend_list <- dendlist(dend.nt, dend.aa)
  r.g.p <- round(entanglement(dend_list), digits = 3)
  tanglegram(dend.nt, dend.aa,  
             main = paste("entanglement", 
                          r.g.p, 
                          sep = "="), 
             sub = "L: Gene; R: Protein") 
  
  nt.tangle.aa <- paste("Tanglegram of", 
                        k, sep = " ")
  nt.tangle.aa.file <- paste0("./output/plots/gene_protein/Tanglegram of_", 
                              k, ".pdf")
  
  pdf(nt.tangle.aa.file, width = 9, height = 6)
  tanglegram(dend.nt, dend.aa,  
             main = paste("entanglement", 
                          r.g.p, 
                          sep = "="), 
             sub = "L: Gene; R: Protein") 
  dev.off()
  
  ## ------------- Correlation matrix between a list dendrograms ------------- ###
  
  ## A). Computing the entanglement coefficient between trees. 
  
  r.seq <- c(r.g.s, r.p.s, r.g.p)
  
  coef.mat <- rbind(coef.mat, r.seq)
  
  ## B). Computing the correlation coefficient between trees. 
  
  #. tree_list <- dendlist(dend.nt, dend.aa, dend.sp)
  #. cor.dendlist(tree_list, method = "cophenetic")
  #. cor.dendlist(tree_list, method = "baker")
  #. cor.dendlist(tree_list, method = "common_nodes")
  
}

rownames(coef.mat) <- ko
colnames(coef.mat) <- c("Gene-Species", "Protein-Species", "Gene-Protein")
coef.mat <- as.data.frame(coef.mat)
coef.mat <- 1-coef.mat

coef.mat <- as.matrix(coef.mat)

saveRDS(coef.mat, file = "./output/objects/coef.mat.rds")

### 8) Comparison between gene tree, protein tree and species tree.

d <- dist(coef.mat)
h <- hclust(d)
plot(h)

library(pheatmap)
pheatmap(coef.mat)
pheatmap(coef.mat, kmeans_k = 3)

# Extracting the gene symbols. 
gene.info <- readRDS("./data/gene.info.rds ")

head(gene.info)

kid.sym <- openxlsx::read.xlsx("./data/Common_orthologs.xlsx", 
                               sheet = "All") [, 1]

### End of Step-02.
### ------------------------------------------------------------------------ ###
