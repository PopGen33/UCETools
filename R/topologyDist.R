# Gets distances between trees using TreeDist::TreeDistance() and returns a vector of those distances

#' TreeDistance Wrapper function
#'
#' @description
#' Internal [TreeDist::TreeDistance] Wrapper function for use with parApply; might be avoidable? Matrix of trees rather than indices?
#' @param combo a vector containing exactly 2 items: the indices of a pair of trees. The distance between this pair is
#' evaluated by [TreeDist::TreeDistance]
#' @param treeSet an ape multiPhylo object containing the trees to which the indices in combo refer
#' @returns description Returns a distance between the pair of trees
#' @importFrom TreeDist TreeDistance
#' @keywords internal
#' @export
TreeDistWrapper <- function(combo, treeSet){
    return(TreeDist::TreeDistance(treeSet[combo[1]], treeSet[combo[2]]))
}

#' Get Names For Distances function
#'
#' @description
#' Internal function getting names for the vector returned by toplogyDist for use with parApply
#' @param combo a vector containing exactly 2 items: the indices of a pair of trees. Used to select the trees from treeSet and extract their names
#' @param treeSet an ape multiPhylo object containing the trees to which the indices in combo refer
#' @returns description Returns a distance between the pair of trees
#' @keywords internal
#' @export
getNamesForDistances <- function(combo, treeSet){
    return(paste0(names(treeSet[combo[1]]), "|", names(treeSet[combo[2]])))
}

#' Get distances between trees using [TreeDist::TreeDistance]
#'
#' @description
#' Returns a vector containing generalized Robinson-Foulds distances between trees within a set of trees or between a single
#' 'target' tree and a set of trees using the [TreeDist::TreeDistance] function.This metric is slightly different
#' from how Robinson-Foulds distances were originally defined, but it works with trees that don't necessarily
#' contain exactly the same taxa. Please see the description for the [TreeDist::TreeDistance] function for a more
#' thorough explanation. This matrix returned by this function can be used to visualize the effects of filtering trees
#' by parsimony informative sites
#'
#' @param treeSet A set of trees as an ape multiPhylo object. See, for example, the [ape::read.tree] documentation.
#' @param threads Number of simultaneous processes used for calculating tree distances. By
#' default, this uses the number of cores returned by [parallel::detectCores].
#' @param targetTree Target tree (as an ape 'phylo' object) to which the treeSet is compared. If set, this
#' changes the behaviour of the function substantially. Instead of returning distances between pairs of trees
#' in treeSet the function will
#' return distances between the trees in treeSet and the targetTree.
#' @param samples Integer. If not set (default), distance between all possible sets of trees are returned. If set,
#' a random sample of those distances (with replacement) of size 'samples' is returned instead.
#' @param withNames Boolean. If true, the vector of distances is returned with names such that "tree1|tree2" is the name
#' for the distance from tree1 to tree2. Note that this is unused if the targetTree argument is used (all distances are
#' to the targetTree)
#'
#' @returns a vector containing generalized Robinson-Foulds distances
#'
#'
#' @importFrom  TreeDist TreeDistance
#' @importFrom  parallel makePSOCKcluster detectCores parApply stopCluster clusterEvalQ
#'
#' @export
topologyDist <- function(treeSet, threads = detectCores(), targetTree, samples, withNames = FALSE){
    if(!missing(targetTree)){
        if(!inherits(targetTree, "phylo")){
            stop("'targetTree' must be an object of class 'phylo'.")
        }
        # Get distances to the targetTree
        distances <- TreeDistance(targetTree, treeSet)
        # If samples is set, return that many distances; otherwise return all
        if(!missing(samples)){
            if(!inherits(samples, "numeric") || samples %% 1 != 0){
                stop("'samples' argument must be integer if set.")
            }
            distances <- sample(distances, size = samples, replace = TRUE)
        }
        return(distances)
    }else{
        # Make cluster (for parallelism that works on Windows)
        cl <- makePSOCKcluster(threads)
        clusterEvalQ(cl = cl, library(TreeDist)) # Because the virtual nodes spawn with only base packages loaded

        # All possible tree combinations (expressed as combinations of indices)
        combos <- comboGeneral(length(treeSet), m = 2, nThreads = threads)

        # If samples isn't zero, sample <samples> rows from the combos matrix with replacement
        if(!missing(samples)){
            if(!inherits(samples, "numeric") || samples %% 1 != 0){
                stop("'samples' argument must be integer if set.")
            }
            combos <- combos[sample(nrow(combos), size = samples, replace = TRUE),]
        }

        # Get the distances
        result <- parApply(cl, FUN = TreeDistWrapper, MARGIN = 1, X = combos, treeSet)

        # Name the distances (if requested)
        if(withNames){
            # name the results vector like "tree1|tree2" = dist(tree1, tree2)
            namesVect <- parApply(cl, FUN = getNamesForDistances, MARGIN = 1, X = combos, treeSet)
            names(result) <- namesVect
        }

        # Stop cluster
        stopCluster(cl)

        return(result)
    }
}
