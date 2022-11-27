# calculates parsimony informative sites given a DNAbin object
# calc = absolute means return count of parsimony informative sites
# calc = relative means return the % of parsimony informative sites

#' Count Parsimony informative sites
#'
#' @description
#' Counts the parsimony informative sites in a single alignment. The input must be an alignment in the
#' form of [ape]'s DNAbin object. This is primarily an internal function used by [UCETools::getPisCounts].
#' Note that ambiguity codes are ignored! The only bases considered are the standard DNA bases: A, G, C, T
#'
#' @param DNAbinObj An alignment in the form of [ape]'s DNAbin
#' @param calc A string indicating if the calculation should return the absolute count of parsimony informative sites ('absolute') or the percent of parsimony informative sites relative to the length of the alignment ('relative')
#'
#' @returns If calc is 'absolute', returns an integer count of parsimony informative sites in the alignment. If calc is 'relative', returns a numeric 0 <= n <= 1 representing the percent of parsimony informative sites relative to the length of the alignment
#'
#' @seealso [UCETools::getPisCounts] [UCETools::getAlignmentsByPIS]
#'
#' @import ape
#'
#' @export

countParsimonyInformative <- function(DNAbinObj, calc = "absolute"){
    ## Check input
    if(!inherits(DNAbinObj, "DNAbin")){
        stop("Input must be of class 'DNAbin'")
    }
    if(!calc %in% c("absolute", "relative")){
        stop("'calc' argument must be 'absolute' or 'relative'.")
    }
    ## my own slightly different implementation of par.inf from inside the pis() function:
    ## https://rdrr.io/cran/ips/src/R/pis.R
    ## My version handles '?' properly by only handling the 4 dna bases
    ## Mine is also very slightly faster on average (2% faster was not worth the time I spent thinking about it...)
    pars.inf <- function(charCol){
        validChars <- c('a', 'g', 'c', 't')
        charTable <- table(charCol)
        # Drop characters with exactly 1 occurrence (apomorphic characters)
        charTable <- charTable[charTable > 1]
        return(length(charTable[names(charTable) %in% validChars]) > 1)
    }
    alignLength <- dim(DNAbinObj)[2]
    # Convert DNAbin to character; can't get around this as far as I can tell
    # The binary representation of DNAbin is described here: http://ape-package.ird.fr/misc/BitLevelCodingScheme_20April2007.pdf
    DNAbinObj <- as.character(DNAbinObj)
    countPIS <- sum(apply(DNAbinObj, 2, pars.inf), na.rm = TRUE)
    if(calc == "relative"){
        return((countPIS/alignLength)*100)
    }
    return(countPIS)
}
