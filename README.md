# UCETools
 Simplified UCE filtering and ASTRAL interface from R

This R package facilitates counting parsimony informative sites (including allowing for IUPAC ambiguity codes) in alignments of DNA data and then filtering loci based on those counts. 
This is especially helpful for UCEs where counts of parsimony informative sites for some loci are near zero and lead to incorrectly infered gene trees that can mislead some species 
tree inferences methods (e.g. ASTRAL). This package also provides a (very) simple ASTRAL interface from R.

## Setup and basic tutorial
This is a stopgap until I can add more explicit usage instructions here. In addition to the below, please also call `help("functionName")` in R to see full documentation for the functions. You 
Can also see all available functions in the package by calling `help(package = "UCETools")`.

### Install
To install the package from github, simply load the devtools library (install it first if you need to), and the use the install_github() function to unstall it.
```
library(devtools)
install_github("PopGen33/UCETools")
library(UCETools)
help(package = "UCETools")
```
### Basic Workflow
Alignments are expected to be in a directory with each locus in a seperate file. These files can be in either nexus or fasta format (this is controlled by the format option; 
please see this functions help() for usage). To read these files and count parsimony informative sites, call the `getPisCounts()` function. Afterwards, you can use the 
`getAlignmentsByPIS()` function to filter the alignemnts and export them to a new directory.
```
# Get alignments and count parsimony informative sites; see help for options
help("getPisCounts")
pisCounts <- getPisCounts(alignDir = "path/to/your/alignments", format = "fasta")
# If you'd like to load some example alignments to play around with, you can call:
pisCounts_example <- getPisCounts(alignDir = system.file("extdata", "UCE_aligns", package = "UCETools"))

# Filter alignments alignments with the most parsimony informative sites (in this case, the top 50%)
help("getAlignmentsByPIS")
pisCounts_filtered_top50percent <- getAlignmentsByPIS(pisData = pisCounts, cutoff = 0.5, outputDir = "path/to/your/outputAlignments", returnPisData = TRUE) 
```
At this point, you can find the filtered alignments in the output directory you specified. When called as above (returnPisData = TRUE), it also returns a new 
PISMatrix with just the loci that passed the filter. Then... you need to use that information to retrieve the genetrees for these loci. When I made these scripts, 
I used an external shell script to retrieve the genetrees based on the files in the output directory (just mataching names), but that's more than I'd like to ask of 
people using this package, so I'm working on a way to associate genetrees with the PISMatrix. That isn't implemented...yet. If you have genetrees, you can use the 
simple ASTRAL interface I provide (assuming you have java installed!):
```
# Automatically retrieve the ASTRAL jar file; see help() to control where this downloads to
getAstral()
allUCEtrees <- ape::read.tree("path/to/your/trees/trees.treefile") # multiple trees in one file, or...
allUCEtrees <- readTreesDir("path/to/your/trees") # when you have many files each containing one newick formatted tree in a directory; again, see help("readTreesDir")
allUCEtrees_ASTRALtree <- callAstral(allUCEtrees, outfile = "path/to/your/astralOutput.nwk", annotation = 1, loadOutputTree = TRUE)
plot(allUCEtrees_ASTRALtree, show.node.label = TRUE, main = "ASTRAL Tree: All UCEs")
```

## Citation
For now, I'd recommend using the R package `grateful` and simply calling `citation("UCETools")` with the package installed. That will generate a bibtex citation. 
This citation style will change in the future. This will fail to return a date, but please use "2024". Again, this should be more straightforward in the future.
```
library(grateful)
citation("UCETools")
```