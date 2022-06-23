library(stringr)

# Script to generate Gencode.annotation.RData file (features_gene, features_transcript and tx2gene dataframes) from Gencode gtf files.
# Used to generate Gencode.v40.annotation.RData and Gencode.v19.annotation.RData in https://github.com/BHKLAB-Pachyderm/Annotations.git
# gencode.v40.annotation.gtf (GRCh38) obtained from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz
# gencode.v19.annotation.gtf (GRCh37) obtained from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

work_dir <- "~/Documents/PredictIO/curation/gene_annotation"
gtf_filename <- "gencode.v19.annotation.gtf"
rdata_filename <- "Gencode.v19.annotation.RData"

get_features_df <- function(gtf_df, feature_name, features_colnames){
  features_df <- gtf_df[gtf_df$V3 == feature_name, ]
  colnames(features_df) <- c("seqnames", "source", "type", "start", "end", "strand", "other")
  features_df$score <- NA
  features_df$phase <- NA
  other_values <- lapply(features_df$other, function(line){
    split_trimmed <- str_trim(str_split(line, ";")[[1]])
    split_trimmed <- split_trimmed[split_trimmed != ""]
    key_vals = lapply(split_trimmed, function(token){
      key_val <- str_split(token, "\\s")[[1]]
      l <- list()
      l[[key_val[1]]] <- key_val[2]
      return(l)
    })
    return(unlist(key_vals))
  })
  
  col_vals <- lapply(features_colnames, function(colname){
    return(
      unlist(lapply(other_values, function(row){
        if(colname %in% names(row)){
          return(row[[colname]])
        }else{
          return(NA)
        }
      }))
    )
  })
  names(col_vals) <- features_colnames
  
  for(colname in names(col_vals)){
    features_df[[colname]] <- col_vals[[colname]]
  }
  return(features_df)
}

features_gene_cols <- c(
  "gene_id", "gene_type", "gene_name", "level", "havana_gene", "tag"
)
features_transcript_cols <- c(
  "gene_id", "gene_type", "gene_name", "level", "hgnc_id", "havana_gene", "transcript_id", "transcript_type", "transcript_name", "transcript_support_level", 
  "tag", "havana_transcript", "ont", "protein_id", "ccdsid"
)

gtf_data <- read.table(file.path(work_dir, gtf_filename), header=FALSE, sep="\t")
gtf_data <- gtf_data[, !(names(gtf_data) %in% c("V6", "V8"))]

features_gene <- get_features_df(gtf_data, "gene", features_gene_cols)
features_gene <- features_gene[, !names(features_gene) %in% c("other")]
rownames(features_gene) <- features_gene$gene_id

features_transcript <- get_features_df(gtf_data, "transcript", features_transcript_cols)
features_transcript <- features_transcript[, !names(features_transcript) %in% c("other")]
rownames(features_transcript) <- features_transcript$transcript_id

tx2gene <- data.frame(transcripts=rownames(features_transcript), genes=features_transcript$gene_id)

save(features_gene, features_transcript, tx2gene, file=file.path(work_dir, rdata_filename))