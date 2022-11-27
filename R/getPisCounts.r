# Creates the matrix containing file paths to alignments, their DNAbin representations,
# and counts of parsimony informative sites (using countParsimonyInformative()).
# This matrix has a special class ('PISdata') and is used as input for getAlignmentsByPIS()
# Uses the parallel library to parallelize things in a way that should work on Windows or Unix-like

#' Get counts of parsimony informative sites
#'
#' @description
#' Counts the parsimony informative sites in multiple alignments (input as a path to a directory
#' containing alignments). Returns a matrix of class 'PISdata' that's used by [UCETools::getAlignmentsByPIS]
#' to subset alignments.
#'
#' @param alignDir Path to the directory containing alignment files
#' @param ext The file extension of the alignment files (e.g. "file.nexus" is matched by ext = "nexus"). If not set, all files in the directory are used.
#' @param calc A string indicating if the calculation should return the absolute count of parsimony informative sites ('absolute') or the percent of parsimony informative sites relative to the length of the alignment ('relative')
#' @param threads Number of simultaneous processes used for counting parsimony infomative sites. By
#' default, this uses the number of cores returned by [parallel::detectCores].
#' @param format Either 'nexus' or 'fasta' indicating what format the alignments are stored in. The default is 'fasta'
#'
#' @returns A matrix of class 'PISdata'. The first column contains paths to alignment files. The second column
#' contains DNAbin representations of the alignments. The third column contains counts of parsimony informative sites.
#' The fourth column contains alignment lengths.
#'
#' @seealso [UCETools::countParsimonyInformative] [UCETools::getAlignmentsByPIS]
#'
#' @importFrom  ape read.FASTA read.nexus.data nexus2DNAbin
#' @importFrom  parallel makePSOCKcluster detectCores parLapply stopCluster
#'
#' @export

getPisCounts <- function(alignDir, ext, calc = "absolute", threads = detectCores(), format = "fasta"){
    ## Check inputs
    if(!dir.exists(alignDir)){
        stop("Alignment dirctory does not exist or is not a directory.")
    }
    if(!format %in% c("fasta", "nexus")){
        stop("'format' must be 'fasta' or 'nexus'.")
    }
    if(!calc %in% c("absolute", "relative")){
        stop("'calc' argument must be 'absolute' or 'relative'.")
    }
    ## Get file paths
    if(!missing(ext)){
        alignFile <- list.files(alignDir, full.names = TRUE, pattern = paste0("*.", ext))
    }else{
        alignFile <- list.files(alignDir, full.names = TRUE)
    }
    ## Parallelism starts here
    cl <- makePSOCKcluster(threads)
    cat("Reading alignments...")
    if(format == "fasta"){
        DNAbin <- parLapply(cl, alignFile, read.FASTA)
    }else{
        # we're just going to assume it's nexus since we check if it was one of these before
        DNAbin <- parLapply(cl, alignFile, read.nexus.data)
        DNAbin <- parLapply(cl, DNAbin, nexus2DNAbin)
    }
    # DNAbin objects can have different underlying types (in this case need to convert list to matrix)
    DNAbin <- parLapply(cl, DNAbin, as.matrix)

    cat("Done.\n")

    cat("Counting Parsimony Informative Sites...")
    pisCount <- parLapply(cl, DNAbin, countParsimonyInformative, calc = calc)
    cat("Done.\n")

    cat("Getting Alignment Lengths...")
    alignmentLength <- parLapply(cl, DNAbin, ncol)
    cat("Done.\n")

    stopCluster(cl)
    ## Parallelism ends here
    # accessing the DNAbin is weird after this; use like 'pisMatrix[row,2]$DNAbin'
    pisMatrix <- as.matrix(cbind(alignFile, DNAbin, pisCount, alignmentLength))

    ## Sort the results
    ## ***MAKE SURE TO CHECK IF IT'S STILL SORTED IN THE OTHER SCRIPT. NO ASSUMING***
    ## highest pisCount go on top
    cat("Sorting...")
    pisMatrix <- pisMatrix[order(unlist(pisMatrix[,3]), decreasing = TRUE),]
    cat("Done.\n\n")

    # Print some info
    cat(paste("Read", length(alignFile), "alignments.\n"))
    cat("Summary of Parsimony Informative Site Counts:\n\n")
    print(summary(unlist(pisCount)))
    # Make it a special class so I can check this in the next function
    class(pisMatrix) <- c(class(pisMatrix), "PISdata")
    return(pisMatrix)
}
