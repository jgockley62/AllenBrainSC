BiocManager::install("snpStats")
install.packages("GenomicTools.fileHandler")
BiocManager::install("Biostrings")
BiocManager::install("BSgenome")

source('AllenBrainSC/R/functions.R')

gtf_stats <- function( i, data, genome){
  gene_model <- data[data$gene_symbol == i,]
  
  gr <- GenomicRanges::reduce(
    GenomicRanges::GRanges(
      paste0(
        names(genome),':',
        gene_model$V4, '-',gene_model$V5, ':',
        gene_model$V7
      )))
  
  # calculate the length
  len <- sum(GenomicRanges::width(gr))
  
  ## Get Exon Sequence
  seq <- BSgenome::getSeq(genome, gr)
  
  ## Get Base Frequency
  alphafreq <- BSgenome::alphabetFrequency(seq)
  totalfreq <- colSums(alphafreq)
  
  ## get G/C content as a fraction
  gc <- (sum(totalfreq[c('G','C')]) /
           sum(totalfreq[c('A','G','C','T')])) * 100
  return(c(i, gc, len))
}

synapser::synLogin()

gtfID <-'syn26434482'
syns_used <- c(gtfID)
gtf <- synapser::synGet(gtfID)


feature_gtf <- GenomicTools.fileHandler::importGTF(
  file=gtf$path,
  level='gene',
  features=c("gene_id",
             "gene_symbol",
             'transcript_id')
)
biom <- feature_gtf[ , c("V1", "gene_id", "gene_symbol", "transcript_id")]

mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
test = biomaRt::getBM(attributes = c('ensembl_gene_id', 
                                   "hgnc_symbol", "entrezgene_id"), 
                    filters = "hgnc_symbol", 
                    values = biom$gene_symbol,
                    mart = mart)
#########
# Fasta:

fastaID <-'syn26485910'
syns_used <- c(syns_used, fastaID)
genome <- synapser::synGet(fastaID)

dna <- Biostrings::readDNAStringSet(genome$path)
foo <- do.call(rbind,strsplit(names(dna),  ' '))
names(dna) <- foo[,1]

chroms <- list()
for(contig in names(table(biom$V1))){
  chroms[[contig]] <- dna[grep(contig,names(dna),value=T)][1]
}

# Translation for fasta chrom name alterations
region_names <- rep(NA, length(names(chroms)))
names(region_names) <- names(chroms)
for(nam in names(region_names)){
  region_names[nam] <- names(chroms[[nam]])
}

# Gene named vector with of chromosome designations
#gene_chr <- biom$V1
#names(gene_chr) <- biom$gene_id
gene_chr <- biom$gene_symbol


cl <- snow::makeCluster(15, outfile = "log.log")
stats <- matrix(NA,0,3)
for (element in names(table(biom$V1))) {
  calcs <- do.call(rbind, parallel::parLapply(
    cl = cl,
    as.list(biom[biom$V1 == element, ]$gene_symbol),
    gtf_stats,
    data = as.data.frame(feature_gtf)[feature_gtf$V1 == element,],
    genome = chroms[[element]]
  ))
  stats <- rbind(stats,calcs)
}
snow::stopCluster(cl)
rm(cl)

stats <- as.data.frame(stats)
colnames(stats) <- c('gene_symbol', "percentage_gene_gc_content",
                     "gene_length")
#colnames(biom) <- c("chromosome_name", filters, "hgnc_symbol",
#                    "gene_biotype")
biomart_results <- feature_gtf %>%
  dplyr::full_join(stats, by = "gene_symbol")

biomart_results <- biomart_results[ ,c('V1', 'gene_symbol', 'percentage_gene_gc_content', 'gene_length')]
colnames(biomart_results)[1] <- 'chromosome_name'

gene_biotype
entrezgene_id

mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
test = biomaRt::getBM(attributes = c('gene_biotype', 
                                     "hgnc_symbol"), 
                      filters = "hgnc_symbol", 
                      values = biomart_results$gene_symbol,
                      mart = mart)
test <- test[!duplicated(test$hgnc_symbol),]
colnames(test) <- c('gene_biotype', 'gene_symbol')

biomart_results <- biomart_results %>%
  dplyr::full_join(test, by = "gene_symbol")
biomart_results$hgnc_symbol <- biomart_results$gene_symbol
biomart_results <- biomart_results[, c('gene_symbol', 'hgnc_symbol', 'percentage_gene_gc_content',
                 'gene_biotype', 'chromosome_name', 'gene_length')]

chrs <- c( "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
           "14", "15", "16", "17", "18", "19", "20", "21",  "22", "X", "Y", "MT") 

names(chrs) <- c("NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000004.12",
                 "NC_000005.10", "NC_000006.12", "NC_000007.14", "NC_000008.11",
                 "NC_000009.12", "NC_000010.11", "NC_000011.10", "NC_000012.12",
                 "NC_000013.11", "NC_000014.9",  "NC_000015.10", "NC_000016.10",
                 "NC_000017.11", "NC_000018.10", "NC_000019.10", "NC_000020.11",
                 "NC_000021.9",  "NC_000022.11", "NC_000023.11", "NC_000024.10", 
                 "NC_012920.1")

biomart_results$chromosome_name <- as.character(chrs[biomart_results$chromosome_name])

write.table(biomart_results, file = 'test_biomart_object.tsv', row.names = F, 
            col.names = T, sep = '\t', quote = F)

##############
syns_used <- c(syns_used,'syn25883034')
meta <- as.data.frame(
  data.table::fread(
    synapser::synGet('syn25883034')$path
  )
)

# Gene expression Matrix
syns_used <- c(syns_used,'syn25883506')
gene_exp <- as.data.frame(
  data.table::fread(
    synapser::synGet('syn25883506')$path
  )
)
# Assign row names
row.names(gene_exp) <- gene_exp$sample_name
gene_exp <- as.data.frame(gene_exp[,2:dim(gene_exp)[2]])
# transpose - assigns gene features to rows, cell IDs to columns
gene_exp <- as.data.frame(t(gene_exp))
# remove gena features with no counts in any cell
gene_exp <- gene_exp[ rowSums(gene_exp != 0) > 0,] 

#Ad_At <- filter_test(exp = gene_exp, met = meta, is.broad = FALSE, pcnt = .5, 
#                     value = 0, cell_t = 'Exc_L2-3_LINC00507_RPL9P17')

####################
# Clean the Metadata
#Remove leading/trailing white space and replace spaces with '_'
meta <- as.data.frame(apply( meta, 2, rep_space ))

#Identify Broad Cell Types
meta$broad_cell <- do.call(
  rbind,
  stringr::str_split(
    meta$cluster_label,
    '_', 
    n = 2, 
    simplify = FALSE
  )
)[,1]

meta$broad_cell[is.na(meta$broad_cell)] <- 'Unk'

meta[ is.na(meta$cell_type_alias_label),]$cell_type_alias_label <- 
  meta[ is.na(meta$cell_type_alias_label),]$outlier_type

meta <- meta[ meta$outlier_call == 'FALSE',]

meta <- meta[,c('sample_name','class_label', 'donor_sex_label', 'region_label',
                'cortical_layer_label', 'cell_type_alias_label', 
                'external_donor_name_label', 'broad_cell')]

write.csv(meta, file = 'test_sageseqr_meta.csv', row.names = F, quote = F)

gene_exp$gene_symbol <- row.names(gene_exp)
gene_exp <- gene_exp[,c('gene_symbol', meta$sample_name)]
write.table(gene_exp, file = 'test_sageseqr_exp.tsv', row.names = F,
            col.names = F, sep='\t', quote = F)

################################################################################
# Setting Synapse ID's
parentID = 'syn26485899'

this_repo <- githubr::getRepo(
  repository = 'jgockley62/AllenBrainSC', 
  ref="branch", 
  refName='function_dev'
)
prov <- githubr::getPermlink(
  repository = this_repo,
  repositoryPath = 'scripts/sageseqR_format.R'
)

#Biomart Object
bm <- synapser::synStore( synapser::File(
    path = 'test_biomart_object.tsv',
    name = 'Input Biomart Object',
    parentId = parentID),
  used = syns_used,
  executed = prov,
  activityName = 'Input Biomart Object',
  activityDescription = 'Object containing biomart-object information'
)

#Covariates Object
meta <- synapser::synStore( synapser::File(
  path = 'test_sageseqr_meta.csv',
  name = 'Input Metadata Object',
  parentId = parentID),
  used = syns_used,
  executed = prov,
  activityName = 'Input Metadata Object',
  activityDescription = 'Object containing Metadata information'
)

#Expression Object
counts <- synapser::synStore( synapser::File(
  path = 'test_sageseqr_exp.tsv',
  name = 'Input Expression Object',
  parentId = parentID),
  used = syns_used,
  executed = prov,
  activityName = 'Input Expression Object',
  activityDescription = 'Object containing Expression information'
)



