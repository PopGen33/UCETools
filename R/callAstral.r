# call ASTRAL given at least a set of trees and an output file
# checks global option for "astralPath" if not passed in as an argument

#' Run Astral
#'
#' @description
#' Calls Astral on a set of input genetrees. Astral returns the unrooted species that is
#' most compatible with quartets sampled from the input genetrees. See [Astral's github](https://github.com/smirarab/ASTRAL)
#' for more info on the program.
#'
#' @param help Boolean. If TRUE, all other inputs are ignored and the function simply prints Astral's help message
#' @param trees Either a path to a file containing newick formatted genetrees or an ape multiPhylo object containing genetrees. Used as input for Astral
#' @param outfile Path to file where output tree from Astral will be stored
#' @param astralPath Path to the Astral .jar file; if not set, checks the global astralPath option (set by [UCETools::getAstral()]) for its location
#' @param annotation Sets Astral '-t' option. The type of annotation Astral performs on its output tree; See Astral' help message for more details
#' @param loadOutputTree Boolean. If TRUE, function returns the output Astral tree as an [ape] phylo object. If FALSE, the function returns nothing
#'
#' @returns If help is TRUE, the function prints Astral's help message and returns nothing. If loadOutputTree is FALSE, the function returns nothing and simply runs Astral. If loadOutputTree is TRUE, the function returns the tree output from Astral as an [ape] phylo object.
#'
#' @seealso [UCETools::getAstral]
#'
#' @import ape
#'
#' @export

callAstral <- function(trees, outfile, astralPath, annotation = 3, loadOutputTree = FALSE, help = FALSE){
    ## Help overrides all other options. Print astral --help and return.
    if(help){
        system2(command = "java", args = c("--help"), wait = TRUE)
        return()
    }
    ## validate inputs
    if(missing(trees)){
        stop("'trees' argument is required. either not set or is not a path to a file that exists.")
    }
    if(!inherits(trees, "multiPhylo")){
        # if it isn't a mutliPhylo object, we should check if the file exists; can type error
        if(!file.exists(trees)){
            stop("'trees' argument is not class 'multiPhylo' AND is not a path to a file that exists.")
        }
    }
    if(annotation %% 1 != 0 | annotation < 0){
        stop("Annotation should be integer and positive. See Astral manual.")
    }
    if(missing(outfile)){
        stop("'outfile' argument is required.")
    }
    # maybe should add test of outfile being a valid file path? How though?
    if(missing(astralPath)){
        if(!is.null(getOption("astralPath"))){
            astralPath = getOption("astralPath")
        } else {
            stop("'astralPath' argument not set AND astralPath global option not set.")
        }
    }

    ## Test if trees inherits multiPhylo
    # Need a temp file if input is multiPhylo
    # If it's just a path to a file with newick trees in it, no prep is needed
    if(inherits(trees, "multiPhylo")){
        # If multiphylo, need to prep temp file for astral input containing the trees
        treesLoc <- tempfile()
        write.tree(phy = trees, file = treesLoc)
        trees <- treesLoc
    }

    ## Finally, call ASTRAL
    # shQuote preps args for the terminal
    system2(command = "java", args = c("-jar", shQuote(astralPath), "-i", shQuote(trees), "-o", shQuote(outfile), "-t", annotation), wait = TRUE)

    ## Load tree if user requested it using ape's read.tree
    if(loadOutputTree){
        return(read.tree(outfile))
    }
}
