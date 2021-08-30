#' Replace Spaces
#' 
#' This function replaces the spaces in a charcter vector with a user
#' specified character. Specifically ment to be used in an apply loop over
#' a data frame's columns
#' 
#' @export
#' @usage 
#' 
#' @param input a character vector to replace spaces in
#' @param repchar a character to replace the spaces with default = '_'
#' 
#' @examples 
#' example <- c('abc', 'a b c', 'a  b   c', '    ab c ', '', '  ')
#' rep_space( input=example, repchar='_' )
rep_space <- function( input, repchar='_' ){
  #Clean leading spaces
  output <- trimws(input,which = c("both"),whitespace = "[ \t\r\n]")
  
  #Replace in character spaces
  output <- gsub('\\s+',repchar,output)
  
  #Replace empty values with NA
  output <- dplyr::na_if(output, '')
  
  return(output)
}
####################
#' Expression Filter Test
#' Tests a filter schema on a user defined cut of the data
#' 
#' @usage 
#' 
#' @param exp An expression data set of raw counts
#' @param met A metadata object with columnl lables of; region_label, 
#' external_donor_name_label, broad_cell, cell_type_alias_label, and 
#' sample_name. These correspond to brain region, donor ID, broad cell type 
#' label, specific cell type label, and sample name. 
#' @pcnt A numeric between 0 and 1 corresponding to a the rank value of the
#' feature's expression to apply the filter of `value`. default = 0.5 (median)
#' @value A numeric value to apply the filter to the `pcnt` value
#' @cell_t character vector of meta data column values to select for broad OR 
#' specific cell type values. If values are broad cell type(s) `is.broad` must
#' be set to TRUE (Optional)
#' @region a character vector of brain region(s)/tissue(s) to subset the data 
#' by (Optional)
#' @donor_id value of donor ID's to subset data by. (Optional)
#' @is.broad if cell_t is assigned then is.broad must be set to TRUE 
#' (cell types are broad cell types) or FALSE (cell types are specific cell 
#' types)
#' 
#' @return a list object of filtered counts matrix (list$exp), filtered meta data
#' (list$met), number of features remaining (list$features), and percent of counts 
#' retained after filtering counts based on percentile cutoff in the filtered 
#' cell population (list$count_coverage)
#' 
#' @export
#' 
#' @examples 
#'set.seed(42)
#'s_meta <- as.data.frame(cbind(
#'  sample_name = sample(random::randomStrings(
#'    n=1000, len=6, digits=F, upperalpha=TRUE, loweralpha=F, unique=T), 1000, 
#'    replace = F),
#'  cell_type_alias_label = sample( random::randomStrings(
#'    n=17, len=6, digits=F, upperalpha=TRUE, loweralpha=F, unique=F), 1000, 
#'    replace = T),
#'  region_label = sample( random::randomStrings(
#'    n=7, len=3, digits=F, upperalpha=TRUE, loweralpha=F, unique=F), 1000, 
#'    replace = T),
#'  external_donor_name_label = sample( random::randomStrings(
#'    n=4, len=6, digits=F, upperalpha=TRUE, loweralpha=F, unique=F), 1000, 
#'    replace = T),
#'  broad_cell = sample( random::randomStrings(
#'    n=7, len=5, digits=F, upperalpha=TRUE, loweralpha=F, unique=F), 1000, 
#'    replace = T)
#'))
#'
#'s_exp <- (rbind(
#'  rnbinom(n=1000, 20, p=.75 ),
#'  rnbinom(n=1000, 20, p=.95 ),
#'  rnbinom(n=1000, 20, p=.25 ),
#'  rnbinom(n=1000, 4, p=.25 ),
#'  rnbinom(n=1000, 4, p=.95 ),
#'  rnbinom(n=1000, 4, p=.75 ),
#'  rnbinom(n=1000, 8, p=.25 ),
#'  rnbinom(n=1000, 8, p=.95 ),
#'  rnbinom(n=1000, 8, p=.75 ),
#'  rnbinom(n=1000, 14, p=.5 )
#'))
#'s_exp <- s_exp[rep(1:nrow(s_exp), times = 100), ]
#'colnames(s_exp) <- s_meta$sample_name
#'example <- filter_test(met = s_meta, exp = s_exp, pcnt = .9, value = 19)

filter_test <- function( exp, met, pcnt = .5, value = 1, 
                         cell_t = NULL, donor = NULL, region = NULL, 
                         is.broad = NULL) {
  if(missing(cell_t)){
    cell_t <- NULL
  }
  if(missing(donor)){
    donor <- NULL
  }
  if(missing(region)){
    region <- NULL
  }
  if(missing(is.broad)){
    is.broad <- NULL
  }
  #Confirm is.broad is set correctly if cell_t is specified
  if(is.null(cell_t)) {
  }else{
    if(!is.null(cell_t) & is.null(is.broad)) {
      stop(
        paste0('Must specify is.broad as TRUE or FALSE if cell_t is specified')
      )
    }else{
      if(!is.null(cell_t) & (isTRUE(is.broad) | isFALSE(is.broad))) {
        
      }else{
        stop(
          paste0('is.broad must be set to TRUE or FALSE')
        )
      }
    }
  }
  # Filter meta data
  if(!is.null(region)){
    met <- met[ met$region_label %in% region, ]
  }
  if(!is.null(donor)){
    met <- met[ met$external_donor_name_label %in% donor, ]
  }
  if(!is.null(cell_t)) {
    if(is.broad == TRUE){
      met <- met[ met$broad_cell %in% cell_t, ]
    }else{
      met <- met[ met$cell_type_alias_label %in% cell_t, ]
    }
  }
  
  # Filter gene expression object by cells
  exp <- exp[ ,met$sample_name ]
  
  # Filter gene expression object by expression parameters
  initial_counts <- sum(exp)
  
  ## Filter out genes of X percent zeros
  present <- function(x, percent = pcnt) mean(x > value) >= percent
  filt_exp <- exp[apply(exp, 1, present), ]
  
  features <- dim(filt_exp)[1]
  count_coverage <- signif( 100*(sum(filt_exp)/initial_counts), 4)
  
  return( list( exp = filt_exp,
                met = met,
                features = features,
                count_coverage = count_coverage
  )
  )
  
}

