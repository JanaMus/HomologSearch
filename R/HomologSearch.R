# Tool for identic genes from model and nonmodel organism searching. 
# It uses rBLAST package to create a local databse from model organism sequence and subsequently to blast (find similar genes/proteins/sequences) non-model organism seqeunce to local database (model organism sequence). Finally the tool selects the most similar genes/proteins/sequences.
# Needed tools: BLAST+ (files: blastdbcmd.exe, blastp.exe/blastn.exe, makeblastdb.exe); R packages: rBLAST, Biostrings

## inputs:
# model_organism = fasta file with proteins/nucleotides sequences
# nonmodel_organism = fasta file with proteins/nucleotides sequences
# fasta_type = protein/nucleotide

# Example: HomologSearch("model_organism.fasta", "nonmodel_organism.fasta", "protein")

HomologSearch <- function(model_organism, nonmodel_organism, fasta_type){
library(rBLAST)
library(Biostrings)

# Data import, local database creation
if (fasta_type == "protein"){
  model_seq <- readAAStringSet(model_organism)
  nonmodel_seq <- readAAStringSet(nonmodel_organism)
  makeblastdb(model_organism, dbtype = 'prot') # creates a local database - working directory: NCBI/blast-2.x.x/bin
  local_database <- blast(model_organism, type = "blastp")
  print("protein files imported, database created")
  
} else if (fasta_type == "nucleotide"){
  model_seq <- readDNAStringSet(model_organism)
  nonmodel_seq <- readDNAStringSet(nonmodel_organism)
  makeblastdb(model_organism, dbtype = "nucl") # creates a local database - working directory: NCBI/blast-2.x.x/bin
  local_database <- blast(model_organism, type = "blastn")
  print("nucleotide files imported, database created")
  
} else {
  stop("Incorrect fasta_type. Write 'protein' or 'nucleotide'.")
}


# BLAST query and selection of pairs with Perc.Identity > 50, E-value < 0.001 / 10^-10, and Bits > 50
blast_similarity <- predict(local_database, nonmodel_seq)

if (fasta_type == "protein"){
  blast_similarity_selected <- blast_similarity[(blast_similarity$Perc.Ident > 50 & blast_similarity$E < 0.001 & blast_similarity$Bits > 50),]
} else if (fasta_type == "nucleotide"){
  blast_similarity_selected <- blast_similarity[(blast_similarity$Perc.Ident > 50 & blast_similarity$E < 10^-10 & blast_similarity$Bits > 50),]
}

# Filtering duplicate sequences (deleting the sequence with lower Bits value)
blast_similarity_selected_filtered <- blast_similarity_selected[!duplicated(blast_similarity_selected$QueryID), ]


# Extraction of input sequences names
nonmodel_seq_names <- sub(" \\[.*", "", names(nonmodel_seq))
model_seq_names <- sub(" \\[.*", "", names(model_seq))


# Selection of identic locus_tags and gene names, export to table locus_tags_genes
locus_tags_genes <- matrix(data=NA, nrow = dim(blast_similarity_selected_filtered)[1], ncol = 4)
colnames(locus_tags_genes) = c("nonmodel_locus_tag", "model_locus_tag", "nonmodel_gene", "model_gene")


for (i in 1:dim(blast_similarity_selected_filtered)[1]){
  QueryID <- blast_similarity_selected_filtered[i,1] # name of the nonmodel homolog sequence
  SubjectID <- blast_similarity_selected_filtered[i,2] # name of the model homolog sequence
  
  # export nonmodel_locus_tag name
  sub_nonmodel <- sub(".*locus_tag=", "", names(nonmodel_seq[which(QueryID == nonmodel_seq_names)]))
  locus_tags_genes[i,1] <- sub("].*", "", sub_nonmodel) 
  
  # export model_locus_tag name
  sub_model <- sub(".*locus_tag=", "", names(model_seq[which(SubjectID == model_seq_names)]))
  locus_tags_genes[i,2] <- sub("].*", "", sub_model) 
  
  # IF exists, export nonmodel_gene name
  if (gregexpr(pattern="gene", names(nonmodel_seq[which(QueryID == nonmodel_seq_names)]))[[1]][1] > 0){
    sub_gene_nonmodel <- sub(".*gene=", "", names(nonmodel_seq[which(QueryID == nonmodel_seq_names)]))
    locus_tags_genes[i,3] <- sub("].*", "", sub_gene_nonmodel)
  }
  
  # IF exists, export model_gene name
  if (gregexpr(pattern="gene", names(model_seq[which(SubjectID == model_seq_names)]))[[1]][1] > 0){
    sub_gene_model <- sub(".*gene=", "", names(model_seq[which(SubjectID == model_seq_names)]))
    locus_tags_genes[i,4] <- sub("].*", "", sub_gene_model)
  }
}

write.table(locus_tags_genes, "homologs.csv", sep = ";", row.names = F)
print(paste0("Output table homologs.csv is created with ", dim(blast_similarity_selected_filtered)[1], " homologs."))

}