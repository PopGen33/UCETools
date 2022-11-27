# Uses the PISdata matrix to retrieve the top <count> or top <percent>% of alignments and
# puts the alignments in a new directory. Also returns the post-filtering pisData matrix if requested

#' Get top alignments by count of parsimony informative sites
#'
#' @description
#' Uses the PISdata matrix output by [UCETools::getPisCounts] to create a directory containing copies of the top alignments
#' as ranked by their counts of parsimony informative sites.
#'
#' @param pisData A PISdata matrix produced by [UCETools::getPisCounts]
#' @param cutoff A numeric 0 < n < 1 or an integer n >= 1. If cutoff is < 1, it is interpreted as requesting the top
#' loci above the specified percentile as ranked by PIS (e.g. cutoff = 0.8 requests loci with counts of parsimony
#' informative sites greater than or equal to the 80th percentile). If cutoff is an integer >= 1, it is interpreted as
#' requesting the top &lt;cutoff&gt; loci (e.g. cutoff = 20 requests the top 20 loci). NOTE: requesting an integer
#' number of loci returns EXACTLY that many loci even if there are ties (e.g. if you request 20 loci, you get 20 loci
#' even if loci 21-30 are tied with the 20th ranked locus by counts of parsimony informative sites).
#' @param outputDir The directory to which the requested subset of alignments will be copied
#' @param returnPisData Boolean. If TRUE, the function returns a PISdata matrix containing only the subset
#' of loci remaining after filtering. If FALSE, the function returns nothing.
#'
#' @returns If returnPisData is TRUE, returns a matrix of class 'PISdata' for only the subset of loci remaining after filtering.
#' The first column contains paths to alignment files. The second column
#' contains DNAbin representations of the alignments. The third column contains counts of parsimony informative sites.
#' The fourth column contains alignment lengths. If returnPisData is false, the function returns nothing and only copies
#' the subset of alignments remaining after filtering to outputDir.
#'
#' @seealso [UCETools::countParsimonyInformative] [UCETools::getPisCounts]
#'
#' @importFrom stats quantile
#'
#' @export

getAlignmentsByPIS <- function(pisData, cutoff, outputDir, returnPisData = FALSE){
    ## Check inputs
    if(!inherits(pisData, "PISdata")){
        stop("'pisData' argument must be of class 'PISdata'. 'PISdata' objects are produced by the 'getPisCounts()' function.")
    }
    if(cutoff <= 0){
        stop("'cutoff' must be positive.")
    }
    if(!dir.exists(outputDir)){
        cat("'outputDir' does not exist. Attempting to create new output directory...")
        dir.create(outputDir)
        cat(paste("Done.\nCreated new output directory: ", outputDir, "\n"))
    } else if(length(list.files(outputDir, all.files = TRUE, include.dirs = TRUE, no.. = TRUE)) > 0){
        # warn if outputDir isn't empty
        warning("WARNING: 'outputDir' is not an empty directory! Make sure it contains what you're expecting!")
    }
    # Check if sorted by PIS. This is required for the subsetting to work as intended
    # is.unsorted checks ascending order, but this descending, so need to rev() it
    if(is.unsorted(rev(unlist(pisData[,3])))){
        stop("'pisData' is not sorted. Did this input come from getPisCounts()?")
    }
    if(cutoff >= 1){
        # Test if cutoff is non-integer (must be integer if greater than 1)
        if(cutoff %% 1 != 0){
            stop("If 'cutoff' is greater than 1, it must be an integer.")
        }
        # Handle integer cutoff; this is the easier case if the matrix is sorted
        # This is actually the only reason I need to sort that matrix; could maybe avoid it?
        pisData <- pisData[1:cutoff,]
    }else{
        # Handle percentile cutoff
        # A bit hard to read but, essentially "top loci above <cutoff> percentile as ranked by PIS"
        pisData <- pisData[pisData[,3] >= quantile(unlist(pisData[,3]), probs = cutoff),]
    }

    ## Go get the alignment files we need and put them in outputDir
    # Could maybe parallel this but I'm not sure it would be faster b/c
    # all of the parallel copy functions would need access to the disk
    cat(paste0("Copying selected alignments to ", outputDir, "..."))
    for(alignFile in pisData[,1]){
        file.copy(from = alignFile, to = outputDir, overwrite = TRUE)
    }
    cat("Done.\n\n")
    ## Return pisData if requested
    if(returnPisData){
        cat("Returning subset of 'pisData'.\n")
        # I don't know if this is actually necessary?
        # I think RStudio might not be displaying the class correctly in the global environment...
        pisData <- as.matrix(pisData)
        class(pisData) <- c(class(pisData), "PISdata")
        return(pisData)
    }
}
