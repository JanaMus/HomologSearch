\name{HomologSearch}
\alias{HomologSearch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Package for identic genes from model and nonmodel organism searching. }
\description{
The package computes the similarity between two sequences in the fasta files. The similarity is based on the BLAST method. It uses rBLAST package to create a local databse from model organism sequence and subsequently to blast (find similar genes/proteins/sequences) non-model organism seqeunce to local database (model organism sequence). Finally the tool selects the most similar genes/proteins/sequences and save locus tags and genes (if available) into the output .csv file.
}
\usage{
HomologSearch(model_organism, nonmodel_organism, fasta_type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model_organism}{
%%     ~~Describe \code{model_organism} here~~
}
  \item{nonmodel_organism}{
%%     ~~Describe \code{nonmodel_organism} here~~
}
  \item{fasta_type}{
%%     ~~Describe \code{fasta_type} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
Jana Musilova
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
HomologSearch("model_organism.fasta", "nonmodel_organism.fasta", "protein")

## The function is currently defined as
function (model_organism, nonmodel_organism, fasta_type)
{
    library(rBLAST)
    library(Biostrings)
    if (fasta_type == "protein") {
        model_seq <- readAAStringSet(model_organism)
        nonmodel_seq <- readAAStringSet(nonmodel_organism)
        makeblastdb(model_organism, dbtype = "prot")
        local_database <- blast(model_organism, type = "blastp")
        print("protein files imported, database created")
    }
    else if (fasta_type == "nucleotide") {
        model_seq <- readDNAStringSet(model_organism)
        nonmodel_seq <- readDNAStringSet(nonmodel_organism)
        makeblastdb(model_organism, dbtype = "nucl")
        local_database <- blast(model_organism, type = "blastn")
        print("nucleotide files imported, database created")
    }
    else if (fasta_type == "example") {
        makeblastdb(model_organism, dbtype = "prot")
        local_database <- blast(model_organism, type = "blastp")
        print("example files imported, database created")
    }
    else {
        stop("Incorrect fasta_type. Write 'protein' or 'nucleotide'.")
    }
    blast_similarity <- predict(local_database, nonmodel_seq)
    if (fasta_type == "protein") {
        blast_similarity_selected <- blast_similarity[(blast_similarity$Perc.Ident >
            50 & blast_similarity$E < 0.001 & blast_similarity$Bits >
            50), ]
    }
    else if (fasta_type == "nucleotide") {
        blast_similarity_selected <- blast_similarity[(blast_similarity$Perc.Ident >
            50 & blast_similarity$E < 10^-10 & blast_similarity$Bits >
            50), ]
    }
    else if (fasta_type == "example") {
        blast_similarity_selected <- blast_similarity[(blast_similarity$Perc.Ident >
            50 & blast_similarity$E < 0.001 & blast_similarity$Bits >
            50), ]
    }
    blast_similarity_selected_filtered <- blast_similarity_selected[!duplicated(blast_similarity_selected$QueryID),
        ]
    nonmodel_seq_names <- sub(" \\[.*", "", names(nonmodel_seq))
    model_seq_names <- sub(" \\[.*", "", names(model_seq))
    locus_tags_genes <- matrix(data = NA, nrow = dim(blast_similarity_selected_filtered)[1],
        ncol = 4)
    colnames(locus_tags_genes) = c("nonmodel_locus_tag", "model_locus_tag",
        "nonmodel_gene", "model_gene")
    for (i in 1:dim(blast_similarity_selected_filtered)[1]) {
        QueryID <- blast_similarity_selected_filtered[i, 1]
        SubjectID <- blast_similarity_selected_filtered[i, 2]
        sub_nonmodel <- sub(".*locus_tag=", "", names(nonmodel_seq[which(QueryID ==
            nonmodel_seq_names)]))
        locus_tags_genes[i, 1] <- sub("].*", "", sub_nonmodel)
        sub_model <- sub(".*locus_tag=", "", names(model_seq[which(SubjectID ==
            model_seq_names)]))
        locus_tags_genes[i, 2] <- sub("].*", "", sub_model)
        if (gregexpr(pattern = "gene", names(nonmodel_seq[which(QueryID ==
            nonmodel_seq_names)]))[[1]][1] > 0) {
            sub_gene_nonmodel <- sub(".*gene=", "", names(nonmodel_seq[which(QueryID ==
                nonmodel_seq_names)]))
            locus_tags_genes[i, 3] <- sub("].*", "", sub_gene_nonmodel)
        }
        if (gregexpr(pattern = "gene", names(model_seq[which(SubjectID ==
            model_seq_names)]))[[1]][1] > 0) {
            sub_gene_model <- sub(".*gene=", "", names(model_seq[which(SubjectID ==
                model_seq_names)]))
            locus_tags_genes[i, 4] <- sub("].*", "", sub_gene_model)
        }
    }
    write.table(locus_tags_genes, "homologs.csv", sep = ";",
        row.names = F)
    print(paste0("Output table homologs.csv is created with ",
        dim(blast_similarity_selected_filtered)[1], " homologs."))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
