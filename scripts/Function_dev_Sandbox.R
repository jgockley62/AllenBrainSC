######################################################
#### The Purpose of this sandbox is to test the   #### 
#### functions for handeling the single cell data ####
#### from the allen brain institute. Data is      ####
#### SMART-Seq2 Data                              ####
######################################################
# Source the functions file
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
row.names(gene_exp) <- gene_exp$sample_name
gene_exp <- as.data.frame(gene_exp[,2:dim(gene_exp)[2]])

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

####################
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


####################
# Principle Components of the entire data
geneexp_pca <- prcomp(gene_exp, center = TRUE,scale. = TRUE)



############################################
# Important Questions:
  #1.1 What has more variability, Region or Donor?
  #1.2 How does Braod or Specific Celltype variability compare accross 1.1?
  #1.3 Best Practices for measuring variability. -> PCA?