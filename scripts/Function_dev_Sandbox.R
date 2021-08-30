######################################################
#### The Purpose of this sandbox is to test the   #### 
#### functions for handeling the single cell data ####
#### from the allen brain institute. Data is      ####
#### SMART-Seq2 Data                              ####
######################################################
# Source the functions file
library(dplyr)
source('AllenBrainSC/R/functions.R')
synapser::synLogin()

# Load the Allen Brain Data
meta <- as.data.frame(
  data.table::fread(
    synapser::synGet('syn25883034')$path
  )
)

# Gene expression Matrix
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


Ad_At <- filter_test(exp = gene_exp, met = meta, is.broad = FALSE, pcnt = .5, 
                     value = 0, cell_t = 'Exc_L2-3_LINC00507_RPL9P17')

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
 
####################
############################################
# Important Questions:
#1.1 What has more variability, Region or Donor?
#1.2 Based on different strategies, what is our missingness
#1.3 How does Braod or Specific Celltype variability compare accross 1.1?
#1.4 Best Practices for measuring variability. -> PCA?
############################################
####################
  ##  ##  1.1 What has more variability, Region or Donor?  ##  ##  

  ####Look at tables by Cell Type
  ## Donor Based Tables
    # Donor By Tissue Table
tissue_donor <- kableExtra::kable(
  table(meta$external_donor_name_label, meta$region_label)
) %>%
  kableExtra::kable_styling()                 
tissue_donor

    # Donor By celltype Table
donor_celltype <- kableExtra::kable(
  table(meta$cell_type_alias_label,meta$external_donor_name_label)
) %>%
  kableExtra::kable_styling()                 
donor_celltype

    # Donor By broad celltype Table
donor_broadcelltype <- kableExtra::kable(
  table(meta$broad_cell, meta$external_donor_name_label)
) %>%
  kableExtra::kable_styling()                 
donor_broadcelltype

  ## Region Based Tables
    # Region By Celltype Table
tissue_celltype <- kableExtra::kable(
  table(meta$cell_type_alias_label, meta$region_label)
) %>%
  kableExtra::kable_styling()                 
tissue_celltype

    # Region By Broad Celltype Table
tissue_broadcelltype <- kableExtra::kable(
  table(meta$broad_cell, meta$region_label)
) %>%
  kableExtra::kable_styling()                 
tissue_broadcelltype

  #Exc_L2-3_LINC00507_RPL9P17 cell pop
# Region By Broad Celltype Table
vig_celltype <- kableExtra::kable(
  table(
    meta[meta$cell_type_alias_label %in% 'Exc_L2-3_LINC00507_RPL9P17',]$external_donor_name_label, 
    meta[meta$cell_type_alias_label %in% 'Exc_L2-3_LINC00507_RPL9P17',]$region_label)
) %>%
  kableExtra::kable_styling()                 
vig_celltype

######################################

###### Tests on: Exc_L2-3_LINC00507_RPL9P17
Ad_At <- filter_test(exp = gene_exp, met = meta, is.broad = FALSE, pcnt = .5, 
            value = 0, cell_t = 'Exc_L2-3_LINC00507_RPL9P17')
exc_neuron <- filter_test(exp = gene_exp, met = meta, is.broad = TRUE, pcnt = .5, 
                     value = 0, cell_t = 'Exc')
all_data <- filter_test(exp = gene_exp, met = meta, pcnt = .5, value = 0)

#Build Output
subtype_tests_pairwise <- list()
feature_ret <- as.data.frame(matrix(
  0,
  length(table(meta$external_donor_name_label)),
  length(table(meta$region_label))
))
row.names(feature_ret) <- names(table(meta$external_donor_name_label))
colnames(feature_ret) <- names(table(meta$region_label))
data_ret <- feature_ret
donorxtissue_features <- NULL
donorxtissue_featurelist <- list()
#Run Donor x Tissue Type
for(don in names(table(meta$external_donor_name_label))){
  subtype_tests_pairwise[[don]] <- list()
  for(tis in names(table(meta$region_label))){
    subtype_tests_pairwise[[don]][[tis]] <- 
      filter_test(exp = gene_exp, met = meta,  pcnt = .5, value = 0, 
                cell_t = 'Exc_L2-3_LINC00507_RPL9P17', is.broad = FALSE,
                region = tis, donor = don
              )
    feature_ret[ don,tis ] <- subtype_tests_pairwise[[don]][[tis]]$features
    data_ret[ don,tis ] <- subtype_tests_pairwise[[don]][[tis]]$count_coverage
    donorxtissue_features <- c(
      donorxtissue_features, row.names(subtype_tests_pairwise[[don]][[tis]]$exp)
    )
    donorxtissue_featurelist[[paste0(don,'_',tis)]]<-
      row.names(subtype_tests_pairwise[[don]][[tis]]$exp)
  }
}

data_both <- kableExtra::kable(
  data_ret
) %>%
  kableExtra::kable_styling()                 
data_both

features_both <- kableExtra::kable(
  feature_ret
) %>%
  kableExtra::kable_styling()                 
features_both

 #Examine Histogram:
hist(table(donorxtissue_features), las = 1, 
     main='Feature overlap of Exc_L2-3_LINC00507_RPL9P17: Tissue X Donor',
     xlab = 'Frequency retained',
     ylab = 'Number of Genes')

 #UpSet Plot - Largely Useless with this many comparisons
upset_obj <- ComplexHeatmap::make_comb_mat(
  ComplexHeatmap::list_to_matrix(
    donorxtissue_featurelist
  )
)
ComplexHeatmap::UpSet(upset_obj)

####
  #Look at Donors only
donor_tests <- list()
don_feature_ret <- as.data.frame(matrix(
  0,
  length(table(meta$external_donor_name_label)),
  1
))
row.names(don_feature_ret) <- names(table(meta$external_donor_name_label))
colnames(don_feature_ret) <- 'Features'
don_data_ret <- don_feature_ret
donor_features <- NULL
donor_featurelist <- list()
#Run Donor x Tissue Type
for(don in names(table(meta$external_donor_name_label))){
  donor_tests[[don]] <- 
    filter_test(exp = gene_exp, met = meta,  pcnt = .5, value = 0, 
                cell_t = 'Exc_L2-3_LINC00507_RPL9P17', is.broad = FALSE,
                donor = don
    )
  don_feature_ret[ don,1 ] <- donor_tests[[don]]$features
    don_data_ret[ don,1 ] <- donor_tests[[don]]$count_coverage
    donor_features <- c(
      donor_features, row.names(donor_tests[[don]]$exp)
    )
  donor_featurelist[[don]]<-
    row.names(donor_tests[[don]]$exp)
}

data_don <- kableExtra::kable(
  don_data_ret
) %>%
  kableExtra::kable_styling()                 
data_don

features_don <- kableExtra::kable(
  don_feature_ret
) %>%
  kableExtra::kable_styling()                 
features_don

#Examine Histogram:
hist(table(donor_features), las = 1, 
     main='Feature overlap of Exc_L2-3_LINC00507_RPL9P17: Donor',
     xlab = 'Frequency retained',
     ylab = 'Number of Genes')

#UpSet Plot - Largely Useless with this many comparisons
upset_obj <- ComplexHeatmap::make_comb_mat(
  ComplexHeatmap::list_to_matrix(
    donor_featurelist
  )
)
ComplexHeatmap::UpSet(upset_obj,
                      comb_order = order(ComplexHeatmap::comb_size(upset_obj)))


####
#Look at Tissues only
tissue_tests <- list()
tis_feature_ret <- as.data.frame(matrix(
  0,
  1,
  length(table(meta$region_label)),
))
colnames(tis_feature_ret)<- names(table(meta$region_label))
row.names(tis_feature_ret)  <- 'Features'
tis_data_ret <- tis_feature_ret
tissue_features <- NULL
tissue_featurelist <- list()
#Run tissue x Tissue Type
for(tis in names(table(meta$region_label))){
  tissue_tests[[tis]] <- 
    filter_test(exp = gene_exp, met = meta,  pcnt = .5, value = 0, 
                cell_t = 'Exc_L2-3_LINC00507_RPL9P17', is.broad = FALSE,
                region = tis
    )
  tis_feature_ret[ 1,tis ] <- tissue_tests[[tis]]$features
  tis_data_ret[ 1,tis ] <- tissue_tests[[tis]]$count_coverage
  tissue_features <- c(
    tissue_features, row.names(tissue_tests[[tis]]$exp)
  )
  tissue_featurelist[[tis]]<-
    row.names(tissue_tests[[tis]]$exp)
}

data_tis <- kableExtra::kable(
  tis_data_ret
) %>%
  kableExtra::kable_styling()                 
data_tis

features_tis <- kableExtra::kable(
  tis_feature_ret
) %>%
  kableExtra::kable_styling()                 
features_tis

#Examine Histogram:
hist(table(tissue_features), las = 1, 
     main='Feature overlap of Exc_L2-3_LINC00507_RPL9P17: Tissue',
     xlab = 'Frequency retained',
     ylab = 'Number of Genes')

#UpSet Plot - Largely Useless with this many comparisons
upset_obj <- ComplexHeatmap::make_comb_mat(
  ComplexHeatmap::list_to_matrix(
    tissue_featurelist
  )
)
ComplexHeatmap::UpSet(upset_obj,
                      comb_order = order(ComplexHeatmap::comb_size(upset_obj))
                     )
############################
###### Tests on: All Excitatory Neurons
exc_neuron <- filter_test(exp = gene_exp, met = meta, is.broad = TRUE, pcnt = .5, 
                          value = 0, cell_t = 'Exc')

#Build Output
subtype_tests_pairwise <- list()
feature_ret <- as.data.frame(matrix(
  0,
  length(table(meta$external_donor_name_label)),
  length(table(meta$region_label))
))
row.names(feature_ret) <- names(table(meta$external_donor_name_label))
colnames(feature_ret) <- names(table(meta$region_label))
data_ret <- feature_ret
donorxtissue_features <- NULL
donorxtissue_featurelist <- list()
#Run Donor x Tissue Type
for(don in names(table(meta$external_donor_name_label))){
  subtype_tests_pairwise[[don]] <- list()
  for(tis in names(table(meta$region_label))){
    subtype_tests_pairwise[[don]][[tis]] <- 
      filter_test(exp = gene_exp, met = meta,  pcnt = .5, value = 0, 
                  cell_t = 'Exc', is.broad = TRUE,
                  region = tis, donor = don
      )
    feature_ret[ don,tis ] <- subtype_tests_pairwise[[don]][[tis]]$features
    data_ret[ don,tis ] <- subtype_tests_pairwise[[don]][[tis]]$count_coverage
    donorxtissue_features <- c(
      donorxtissue_features, row.names(subtype_tests_pairwise[[don]][[tis]]$exp)
    )
    donorxtissue_featurelist[[paste0(don,'_',tis)]]<-
      row.names(subtype_tests_pairwise[[don]][[tis]]$exp)
  }
}

data_both <- kableExtra::kable(
  data_ret
) %>%
  kableExtra::kable_styling()                 
data_both

features_both <- kableExtra::kable(
  feature_ret
) %>%
  kableExtra::kable_styling()                 
features_both

#Examine Histogram:
hist(table(donorxtissue_features), las = 1, 
     main='Feature overlap of Excitatory: Tissue X Donor',
     xlab = 'Frequency retained',
     ylab = 'Number of Genes')

#UpSet Plot - Largely Useless with this many comparisons
upset_obj <- ComplexHeatmap::make_comb_mat(
  ComplexHeatmap::list_to_matrix(
    donorxtissue_featurelist
  )
)
ComplexHeatmap::UpSet(upset_obj)

####
#Look at Donors only
donor_tests <- list()
don_feature_ret <- as.data.frame(matrix(
  0,
  length(table(meta$external_donor_name_label)),
  1
))
row.names(don_feature_ret) <- names(table(meta$external_donor_name_label))
colnames(don_feature_ret) <- 'Features'
don_data_ret <- don_feature_ret
donor_features <- NULL
donor_featurelist <- list()
#Run Donor x Tissue Type
for(don in names(table(meta$external_donor_name_label))){
  donor_tests[[don]] <- 
    filter_test(exp = gene_exp, met = meta,  pcnt = .5, value = 0, 
                cell_t = 'Exc', is.broad = TRUE,
                donor = don
    )
  don_feature_ret[ don,1 ] <- donor_tests[[don]]$features
  don_data_ret[ don,1 ] <- donor_tests[[don]]$count_coverage
  donor_features <- c(
    donor_features, row.names(donor_tests[[don]]$exp)
  )
  donor_featurelist[[don]]<-
    row.names(donor_tests[[don]]$exp)
}

data_don <- kableExtra::kable(
  don_data_ret
) %>%
  kableExtra::kable_styling()                 
data_don

features_don <- kableExtra::kable(
  don_feature_ret
) %>%
  kableExtra::kable_styling()                 
features_don

#Examine Histogram:
hist(table(donor_features), las = 1, 
     main='Feature overlap of Excitatory: Donor',
     xlab = 'Frequency retained',
     ylab = 'Number of Genes')

#UpSet Plot - Largely Useless with this many comparisons
upset_obj <- ComplexHeatmap::make_comb_mat(
  ComplexHeatmap::list_to_matrix(
    donor_featurelist
  )
)
ComplexHeatmap::UpSet(upset_obj,
                      comb_order = order(ComplexHeatmap::comb_size(upset_obj)))


####
#Look at Tissues only
tissue_tests <- list()
tis_feature_ret <- as.data.frame(matrix(
  0,
  1,
  length(table(meta$region_label)),
))
colnames(tis_feature_ret)<- names(table(meta$region_label))
row.names(tis_feature_ret)  <- 'Features'
tis_data_ret <- tis_feature_ret
tissue_features <- NULL
tissue_featurelist <- list()
#Run tissue x Tissue Type
for(tis in names(table(meta$region_label))){
  tissue_tests[[tis]] <- 
    filter_test(exp = gene_exp, met = meta,  pcnt = .5, value = 0, 
                cell_t = 'Exc', is.broad = TRUE,
                region = tis
    )
  tis_feature_ret[ 1,tis ] <- tissue_tests[[tis]]$features
  tis_data_ret[ 1,tis ] <- tissue_tests[[tis]]$count_coverage
  tissue_features <- c(
    tissue_features, row.names(tissue_tests[[tis]]$exp)
  )
  tissue_featurelist[[tis]]<-
    row.names(tissue_tests[[tis]]$exp)
}

data_tis <- kableExtra::kable(
  tis_data_ret
) %>%
  kableExtra::kable_styling()                 
data_tis

features_tis <- kableExtra::kable(
  tis_feature_ret
) %>%
  kableExtra::kable_styling()                 
features_tis

#Examine Histogram:
hist(table(tissue_features), las = 1, 
     main='Feature overlap of Excitatory: Tissue',
     xlab = 'Frequency retained',
     ylab = 'Number of Genes')

#UpSet Plot - Largely Useless with this many comparisons
upset_obj <- ComplexHeatmap::make_comb_mat(
  ComplexHeatmap::list_to_matrix(
    tissue_featurelist
  )
)
ComplexHeatmap::UpSet(upset_obj,
                      comb_order = order(ComplexHeatmap::comb_size(upset_obj))
)

##############################################################################
####################
##  ##  1.2 Based on different strategies, what is our missingness?  ##  ## 

broad <- c('Astro', 'Endo','Exc', 'Inh', 'Micro', 'Oligo',  'OPC', 'Unk' )
specific <- c('Astro', 'Exc', 'Inh', 'Micro', 'Oligo',  'OPC', 'Unk' )

b_meta_filt <- meta[ meta$broad_cell %in% broad, ]
s_meta_filt <- meta[ meta$broad_cell %in% broad, ]

hist( table(b_meta_filt$broad_cell), breaks = 100, las = 1, 
      xlab = 'Number of Cells', ylab = 'Number of Populations',
      main = 'Broad Cell Types No Filter' )


hist( table(b_meta_filt$cell_type_alias_label), breaks = 100, las = 1, 
      xlab = 'Number of Cells', ylab = 'Number of Populations',
      main = 'Broad Cell Types No Filter' )
# Dropped in all comparisons due to low numbers: VLMC and Peri







############################################
# Important Questions:
  #1.1 What has more variability, Region or Donor?
  #1.2 How does Braod or Specific Celltype variability compare accross 1.1?
  #1.3 Best Practices for measuring variability. -> PCA?