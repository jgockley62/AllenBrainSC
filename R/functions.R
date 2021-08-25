#' Replace Spaces
#' 
#' This function replaces the spaces in a charcter vector with a user
#' specified character. Specifically ment to be used in an apply loop over
#' a dataframe's columns
#' 
#' @param input a character vector to replace spaces in
#' @param repchar a character to replace the spaces with default = '_'
#' 
#' @examples 
#' example <- c('abc', 'a b c', 'a  b   c')
#' rep_space( input=example, repchar='_' )
rep_space <- function( input, repchar='_' ){
  output <- gsub('\\s+', repchar,input)
  return(output)
}
#' Replace Spaces
#' 
#' This function replaces the spaces in a charcter vector with a user
#' specified character. Specifically ment to be used in an apply loop over
#' a dataframe's columns
#' 
#' @param input a character vector to replace spaces in
#' @param repchar a character to replace the spaces with default = '_'
#' 
#' @examples 
#' example <- c('abc', 'a b c', 'a  b   c')
#' rep_space( input=example, repchar='_' )
rep_space <- function( input, repchar='_' ){
  output <- gsub('\\s+', repchar,input)
  return(output)
}
#' Replace Spaces
#' 
#' This function replaces the spaces in a charcter vector with a user
#' specified character. Specifically ment to be used in an apply loop over
#' a dataframe's columns
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

