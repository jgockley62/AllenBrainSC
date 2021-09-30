library(dplyr)

synapser::synLogin()
seqs <- read.table(synapser::synGet('syn26178508')$path, header = T, sep='\t')
  
#Meta Vars: csv
meta <-  as.data.frame(data.table::fread(synapser::synGet('syn26161132')$path, sep=','))
r_ind <- apply(meta, 1, function(x) all(is.na(x) | x == ''))
c_ind <- apply(meta, 2, function(x) all(is.na(x) | x == ''))
meta <- meta[!r_ind, !c_ind]

# Remove non-variant columns
c_ind <-vapply(meta, function(x) length(unique(x)) > 1, logical(1L))
meta <- meta[ , c_ind]

# Remove useless columns
remove <- c('Run', 'Bases', 'BioSample', 'Bytes', 'DATASTORE filetype', 'Experiment',
            'GEO_Accession (exp)', 'Sample Name', 'V30', 'V31')
meta <- meta[,colnames(meta)[!(colnames(meta) %in% remove)]]

# Remove heavily cell-type confounded variables
meta <- meta[,colnames(meta)[!(colnames(meta) %in% 
                                 c('AvgSpotLen','cell_subtype','source_name')
                              )]]
# Clean age
meta$AGE <- gsub(' years old', '', meta$AGE)

# clean column names
colnames(meta) <- c('sample', 'age', 'sex', 'cell_type')

#Pull seqs
seqs <- seqs[ , c('sample', 'AlignmentSummaryMetrics__PCT_PF_READS_ALIGNED',
          'RnaSeqMetrics__PCT_INTRONIC_BASES',	'RnaSeqMetrics__PCT_INTERGENIC_BASES',
          'RnaSeqMetrics__PCT_CODING_BASES')]
colnames(seqs) <- gsub('__', '_', colnames(seqs))

meta <- meta %>%
  dplyr::left_join(seqs)

internal_parentid <- 'syn25784097'
folder_loc <- 'syn25792226'

#Set Activity
activity <- synapser::synGet('syn26161132')

##Set Annotations:
all.annotations = list(
  dataType = c('clinical','geneExpression'),
  resourceType = 'metadata',
  metadataType = 'analytical covariates',
  isModelSystem = 'FALSE',
  isMultiSpecimen = 'TRUE',
  fileFormat = 'csv',
  grant = 'U01AG046152',
  species = 'Human',
  organ = 'brain',
  tissue = c('hippocampus', 
             'temporal cortex'
  ),
  study = c('Allen_Brain_Single_Cell_Analysis'), 
  consortium = 'AMP-AD',
  assay = 'rnaSeq'
)

#Github Code Pull
thisRepo <- githubr::getRepo(
  repository = "jgockley62/AllenBrainSC",
  ref="branch",
  refName='main'
)
thisFile <- githubr::getPermlink(
  repository = thisRepo,
  repositoryPath=paste0('scripts/','metadata_preprocessing.R'
  )
)

activityName = 'Metadata Input for Sageseqr'
activityDescription = 'Cleaned Metadata for SageSeqR Normalization'

write.csv(meta,
          file = "../validation_metadata_sageseqr.csv",
          row.names = F,
          quote = F)

ENRICH_OBJ <- synapser::synStore( synapser::File( 
  path='../validation_metadata_sageseqr.csv',
  name = 'Metadata Input for Sageseqr',
  parentId=activity$properties$id ),
  used = synids_used,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
file.remove("../validation_metadata_sageseqr.csv")
