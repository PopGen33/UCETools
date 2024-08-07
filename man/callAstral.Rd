% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/callAstral.r
\name{callAstral}
\alias{callAstral}
\title{Run Astral}
\usage{
callAstral(
  trees,
  outfile,
  astralPath,
  annotation = 3,
  loadOutputTree = FALSE,
  help = FALSE
)
}
\arguments{
\item{trees}{Either a path to a file containing newick formatted genetrees or an ape multiPhylo object containing genetrees. Used as input for Astral}

\item{outfile}{Path to file where output tree from Astral will be stored}

\item{astralPath}{Path to the Astral .jar file; if not set, checks the global astralPath option (set by \code{\link[=getAstral]{getAstral()}}) for its location}

\item{annotation}{Sets Astral '-t' option. The type of annotation Astral performs on its output tree; See Astral' help message for more details}

\item{loadOutputTree}{Boolean. If TRUE, function returns the output Astral tree as an \link{ape} phylo object. If FALSE, the function returns nothing}

\item{help}{Boolean. If TRUE, all other inputs are ignored and the function simply prints Astral's help message.
Note that the vast majority of Astral's functions aren't available through callAstral(). The ability to pass custom options may be added in the future if requested.}
}
\value{
If help is TRUE, the function prints Astral's help message and returns nothing. If loadOutputTree is FALSE, the function returns nothing and simply runs Astral. If loadOutputTree is TRUE, the function returns the tree output from Astral as an \link{ape} phylo object.
}
\description{
Calls Astral on a set of input genetrees. Astral returns the unrooted species that is
most compatible with quartets sampled from the input genetrees. See \href{https://github.com/smirarab/ASTRAL}{Astral's github}
for more info on the program.
}
\seealso{
\link{getAstral}
}
