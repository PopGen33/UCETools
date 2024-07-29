# Read a directory of genetrees and return a multiphylo object where the name of each genetree is the file name it was read from

#' Read a directory of genetrees in newick format
#'
#' @description
#' Reads a directory of genetrees in newick format where each file contains one genetree and returns a [ape] multiphylo object where each
#' genetree is named the same as the file it was read from (e.g. when reading a file 'UCEtree01.treefile', the corresponding tree
#' in the multiphylo object will be named 'UCEtree01.treefile). This function mainly exists to make using <futureFunction> easier.
# [UCETools::addTreesToPISData] ### Add this to the above line once implemented
#'
#' @param treesdir A path the the directory containing the genetrees where each file contains exactly one genetree in newick format
#' @param ext The file extension of the alignment files (e.g. "file.treefile" is matched by ext = "treefile"). If not set, all files in the directory are used.
#' @param threads Number of simultaneous processes used for reading trees. By
#' default, this uses the number of cores returned by [parallel::detectCores].
#'
#' @returns An [ape] multiPhylo object containing the trees that were read from the directory where each tree is named as the file it was read from
#'
#' @details The way [ape]'s read.tree() is implemented, text preceding the first '(' is treated as the name of the tree. However, the way
#' I've written this, I assume there is no text preceding the first '('. If there is, this will result in trees being named like '<filename> <precedingText>'.
#'
#' @seealso [ape::read.tree]
# [UCETools::addTreesToPISData] ### Add this to the above line once implemented
#'
#' @importFrom  ape read.tree
#' @importFrom  parallel makePSOCKcluster detectCores parSapply stopCluster
#'
#' @export

readTreesDir <- function(treesdir, ext, threads = detectCores()){
    ## Check inputs
    if(!dir.exists(treesdir)){
        stop("Treesdir dirctory does not exist or is not a directory.")
    }
    ## Get file paths
    if(!missing(ext)){
        treefiles <- list.files(treesdir, pattern = paste0("*.", ext))
    }else{
        treefiles <- list.files(treesdir)
    }
    ## custom function to pass to parSapply
    customReadTrees <- function(treefile, trees_dir){
        return(read.tree(text = paste(treefile, readChar(file.path(trees_dir, treefile), file.info(file.path(trees_dir, treefile))$size)), keep.multi = TRUE))
    }
    ## Parallelism starts here
    cl <- makePSOCKcluster(threads)
    cat("Reading trees...")
    trees_multiphylo <- parSapply(cl, treefiles, FUN = customReadTrees, USE.NAMES = FALSE, trees_dir = treesdir)
    cat("Done.\n")
    cat(paste0("Read ", length(trees_multiphylo), " trees.\n"))
    stopCluster(cl)
    ## Parallelism ends here
    class(trees_multiphylo) <- c("multiPhylo")
    return(trees_multiphylo)
}
