# calculates parsimony informative sites given a DNAbin object
# calc = absolute means return count of parsimony informative sites
# calc = relative means return the % of parsimony informative sites

#' Count Parsimony informative sites
#'
#' @description
#' Counts the parsimony informative sites in a single alignment. The input must be an alignment in the
#' form of [ape]'s DNAbin object. This is primarily an internal function used by [UCETools::getPisCounts].
#' Note that ambiguity codes are ignored by default. The only bases considered are the standard DNA bases: A, G, C, T
#' I have added an experimental feature for including ambiguity codes in the counts following the comparison methods
#' described in [the DNAbin description](https://emmanuelparadis.github.io/misc/BitLevelCodingScheme_20April2007.pdf).
#' This allows for using this function with, for example, RY coded alignments.
#'
#' @param DNAbinObj An alignment in the form of [ape]'s DNAbin
#' @param calc A string indicating if the calculation should return the absolute count of parsimony informative sites ('absolute') or the percent of parsimony informative sites relative to the length of the alignment ('relative')
#' @param useAmbiguity A logical indicating whether to account for IUPAC ambiguity codes (e.g. 'R' for A or G) when counting parsimony informative sites. If FALSE, they are ignored. If TRUE, ambiguity codes are accounted for (e.g. a column 'YYRR' is parsimony informative).\bold{This is an experimental feature. Please check that the values returned are what you expect.}
#'
#' @returns If calc is 'absolute', returns an integer count of parsimony informative sites in the alignment. If calc is 'relative', returns a numeric 0 <= n <= 1 representing the percent of parsimony informative sites relative to the length of the alignment
#'
#' @seealso [UCETools::getPisCounts] [UCETools::getAlignmentsByPIS]
#'
#' @import ape
#' @importFrom RcppAlgos comboGeneral
#'
#' @export

countParsimonyInformative <- function(DNAbinObj, calc = "absolute", useAmbiguity = FALSE){
    ## Check input
    if(!inherits(DNAbinObj, "DNAbin")){
        stop("Input must be of class 'DNAbin'")
    }
    if(!calc %in% c("absolute", "relative")){
        stop("'calc' argument must be 'absolute' or 'relative'.")
    }

    alignLength <- dim(DNAbinObj)[2]

    if(useAmbiguity){
        # DNAbin is described in this document: https://emmanuelparadis.github.io/misc/BitLevelCodingScheme_20April2007.pdf
        # handling the raw binary representation makes taking ambiguity codes into account possible
        pars.inf <- function(rawCol){
            # Throw out characters that are 'N', alignment gap ('-'), or unknown ('?')
            # Implicitly, we're considering gaps strictly uninformative
            rawCol <- rawCol[rawCol != as.raw(240) & rawCol != as.raw(4) & rawCol != as.raw(2)]
            ##TODO: Below necessary?
            # if(length(rawCol)<4){
            #     # Can't be parsimony informative? What about 'NNGAA'?
            #     return(FALSE)
            # }

            # Build table(); R won't run this on raw, so have to convert to numeric
            numTable <- table(as.numeric(rawCol))
            lenNumTable <- length(numTable)

            # If there are less than 2 uniques in a column (table < 2) OR if one of the uniques has a count of length(col) - 1, col can't be parsimony informative
            if(lenNumTable < 2 || (length(rawCol) - 1) %in% numTable){
                return(FALSE)
            }

            # If func hasn't returned at this point, col *might* be informative, but we need to compare chars
            # I'm sure there's some better way of doing this...
            pairs <- RcppAlgos::comboGeneral(lenNumTable, m = 2) # All possible pairs of unique chars (referenced by indices of numTable) in the col
            compRes <- vector("logical", nrow(pairs))
            for(i in 1:nrow(pairs)){
                # Ugly, but works; also, seems to work just calling as.raw() on the string from names rather than making integer first
                # This is my R implementation of DifferentBase(a,b) described in the appendix of the DNAbin document here: https://emmanuelparadis.github.io/misc/BitLevelCodingScheme_20April2007.pdf
                # TRUE implies the characters are surely different
                compRes[i] <- (as.raw(names(numTable[pairs[i,1]])) & as.raw(names(numTable[pairs[i,2]]))) < 16
                # For my own sanity if I need edit this later: Pairs contains pairs of indices in numTable.
                # numTable is a vector with counts at those indices and each index has a name. That name is the
                # integer form of the DNAbin encoding... This is awful
            }

            # If all comparisons are false; col can't be parsimony informative
            if(!any(compRes)){
                return(FALSE)
            }

            # If any character with a count >=2 is surely different than any other character with a count >=2, the site is parsimony informative
            # First remove pairs that aren't different from each other
            pairs <- pairs[compRes,, drop = FALSE]
            # Then, if both members of a pair have counts greater than 2, return TRUE
            for(row in 1:nrow(pairs)){
                if((numTable[pairs[row,][1]] > 1) && (numTable[pairs[row,][2]] > 1)){
                    return(TRUE)
                }
            }
            ##TODO: Should, for example, 'RGYYT' be parsimony informative? I don't think so... It isn't in the above

            # Finally, return FALSE if reaching this point (no pair had both members with counts >2)
            return(FALSE)
        }
        DNAbinObj <- t(sapply(X = DNAbinObj, FUN = as.raw))
        countPIS <- sum(apply(DNAbinObj, 2, pars.inf), na.rm = TRUE)
    }else{
        ## my own slightly different implementation of par.inf from inside the pis() function:
        ## https://rdrr.io/cran/ips/src/R/pis.R
        ## My version ignores '?' properly by only handling the 4 dna bases
        ## Mine is also very slightly faster on average (2% faster was not worth the time I spent thinking about it...)
        pars.inf <- function(charCol){
            validChars <- c('a', 'g', 'c', 't')
            charTable <- table(charCol)
            # Drop characters with exactly 1 occurrence
            charTable <- charTable[charTable > 1]
            return(length(charTable[names(charTable) %in% validChars]) > 1)
        }
        # Convert DNAbin to character; can't get around this as far as I can tell
        # The binary representation of DNAbin is described here: http://ape-package.ird.fr/misc/BitLevelCodingScheme_20April2007.pdf
        DNAbinObj <- as.character(DNAbinObj)
        countPIS <- sum(apply(DNAbinObj, 2, pars.inf), na.rm = TRUE)
    }
    if(calc == "relative"){
        return((countPIS/alignLength)*100)
    }
    return(countPIS)
}
