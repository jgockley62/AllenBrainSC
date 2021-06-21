# Allen Brain Data
# https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq

# Libraries ----------------------------------------------------------------

library(tidyverse)
library(synapser)

# Load Data -----------------------------------------------------------------

# Metadata for Allen Brain dataset
meta <- read.csv("rstudio/metadata.csv",header = TRUE, sep = ',')

# Gene expression Matrix
gene_exp <- read.csv("rstudio/matrix.csv", header = TRUE, sep = ',')

# Analysis -----------------------------------------------------------------

# Counts Per Million (CPM) Normalization by sample 
CPM <- function(x){
  den <- sum(x)
  x/den*(10^6)
}
cpm_exp <- as.data.frame(t(apply(gene_exp[,-1], 1, CPM)))
cpm_exp <- cbind(gene_exp[,1], cpm_exp)
colnames(cpm_exp)[1] = "sample_name"

remove(gene_exp) # Removing un normalized expression matrix to save memory

# Counts of Broad Cell Types based on Brain Region
regional_class <- as.data.frame.matrix(xtabs(formula = ~region_label+class_label, meta))
colnames(regional_class)[1] <- "Unlabelled"

# Counts of Narrow Cell Types based on Brain Region
regional_subclass <- as.data.frame.matrix(xtabs(formula = ~region_label+subclass_label, meta))
colnames(regional_subclass)[1] <- "Unlabelled"

# Counts of Narrow Cell Types for Broad Cell Types
class_sub_rel <- as.data.frame.matrix(xtabs(formula = ~subclass_label+class_label, meta))

# Counts of the most specific distinction vs other two 
cell_type <- as.data.frame.matrix(xtabs(formula = ~cell_type_accession_label+class_label, meta))
cell_type_sub <- as.data.frame.matrix(xtabs(formula = ~cell_type_accession_label+subclass_label, meta))


# Broad Cell Type Subsets Metadata
unlabelled <- subset(meta, class_label == "") # Unlabelled samples
inhibitory <- subset(meta, class_label == "GABAergic") # Inhibitory samples
excitatory <- subset(meta, class_label == "Glutamatergic") # Excitatory samples
non_neuronal <- subset(meta, class_label == "Non-neuronal") # Non-neuronal samples

# Double checking that all the data is accounted for
nrow(meta) == nrow(unlabelled) + nrow(inhibitory) + nrow(excitatory) + nrow(non_neuronal)

unlab_data <- semi_join(cpm_exp, unlabelled, by = 'sample_name') # Unlabelled Gene exp
inhib_data <- semi_join(cpm_exp, inhibitory, by = 'sample_name') # Inhibitory Gene exp
excit_data <- semi_join(cpm_exp, excitatory, by = 'sample_name') # Excitatory Gene exp
non_data <- semi_join(cpm_exp, non_neuronal, by = 'sample_name') # Non-neuronal Gene exp

# Pushing data to synapse -----------------------------------------------------------

synLogin(authToken = "")

# Setting Synapse ID's
parentID = 'syn25881694'
folder_loc = 'syn25881691'

# Set Activity
activity <- synGet('syn25881691')
activityName = 'Allen Brain Data Analysis'
activityDescription = 'Single Cell analysis of Allen Brain Data'


# Counts Table for Broad Cell type per Brain Regions
write.csv(regional_class,
          file = 'Counts_Broad_Type_vs_Brain_Region.csv',
          quote = FALSE
)

CountsTable <- synStore( File(
  path = 'Counts_Broad_Type_vs_Brain_Region.csv',
  name = 'Counts of Broad Cell Types based on Brain Region',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
#synapser::synSetAnnotations(CountsTable, annotations = all.annotations)
file.remove('Counts_Broad_Type_vs_Brain_Region.csv')


# Inhibitory Cell Types
write.csv(inhib_data,
          file = 'Inhibitory_Cells.csv',
          quote = FALSE
)

Inhib <- synStore( File(
  path = 'Inhibitory_Cells.csv',
  name = 'Subset of Inhibitory Cells from Allen Brain data',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
#synapser::synSetAnnotations(Inhib, annotations = all.annotations)
file.remove('Inhibitory_Cells.csv')


# Excitatory Cell Types
write.csv(excit_data,
          file = 'Excitatory_Cells.csv',
          quote = FALSE
)

Excit <- synStore( File(
  path = 'Excitatory_Cells.csv',
  name = 'Subset of Excitatory Cells from Allen Brain data',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
#synapser::synSetAnnotations(Excit, annotations = all.annotations)
file.remove('Excitatory_Cells.csv')


# Non Neuronal Cell Types
write.csv(non_data,
          file = 'Non_Neuronal_Cells.csv',
          quote = FALSE
)

Non <- synStore( File(
  path = 'Non_Neuronal_Cells.csv',
  name = 'Subset of Non Neuronal Cells from Allen Brain data',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
#synapser::synSetAnnotations(Non, annotations = all.annotations)
file.remove('Non_Neuronal_Cells.csv')


# Unlabelled Cells
write.csv(unlab_data,
          file = 'Unlabelled_Cells.csv',
          quote = FALSE
)

Unlab <- synStore( File(
  path = 'Unlabelled_Cells.csv',
  name = 'Subset of Unlabelled Cells from Allen Brain data',
  parentId = activity$properties$id),
  activityName = activityName,
  activityDescription = activityDescription
)
#synapser::synSetAnnotations(Unlab, annotations = all.annotations)
file.remove('Unlabelled_Cells.csv')



