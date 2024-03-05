# Various tests for the UCETools package

# These assume the 'testData' directory is in the R current working directory (and is decompressed, of course)

# libraries
# automatically download uninstalled packages from the registry
# Not sure if I'll need these for testing, but just in case
libNames <- c("devtools","roxygen2")
for(i in 1:length(libNames)){
    if (!libNames[i] %in% rownames(installed.packages())){
        install.packages(libNames[i], dependencies = TRUE)
    }
    library(libNames[i], character.only=TRUE)
}

## Set Working directory (if you need to)
setwd(r"(C:\Users\Conrad\OneDrive-SIU\OneDrive - Southern Illinois University\Zool536\Project)")

## Install the package (it's should be in the current working directory)
#install("UCETools")
# If the github install below doesn't work, use the above to install the local copy

## Install from Github using token to access private repository (expires after 90 days from Nov 27 2022)
install_github("PopGen33/UCETools")
library(UCETools)

## Test help pages
help("getAstral")
help("callAstral")
help("countParsimonyInformative")
help("getPisCounts")
help("getAlignmentsByPIS")

## View tree on which I simulated the testing UCE dataset
simTree <- read.tree("./testData/simulation_tree.nwk")
plot(simTree)
add.scale.bar(x=0,y=10)

## Visualize what one of the alignments looks like (this is just image() called on ape's DNAbin)
exampleAlign <- read.FASTA("./testData/UCE_aligns/simUCE8.fasta")
image(exampleAlign)
#DNAbin must be represented as matrix; this is handled internally by getPisCounts()
countParsimonyInformative(as.matrix(exampleAlign))

## View tree and support for IQTREE concatenated (and partitioned by locus) tree inference
## Suport is ultrafast bootstrap and SH-aLRT (SH-like approximate likelihood ratio test)
## IQTREE call:
## iqtree -p UCE_aligns/ -T auto -m TEST -B 1000 --alrt 1000
iqtree_simTree <- read.tree("testData/UCE_aligns.treefile")
plot(root(iqtree_simTree, outgroup = "A"), show.node.label = TRUE, main = "IQTREE (partitioned by locus)")
add.scale.bar(x=0,y=10)

## Also ran IQTREE to produce individual genetrees for each alignment
## IQTREE call:
## iqtree -S UCE_aligns/ -T 4 -m TEST --prefix UCE_aligns_genetrees/UCE_aligns_genetrees

## Let's read those in and run Astral on them
## TESTS: getAstral() and callAstral()
getAstral()
allUCEtrees <- read.tree("testData/UCE_aligns_genetrees/UCE_aligns_genetrees.treefile")
allUCEtrees_ASTRALtree <- callAstral(allUCEtrees, outfile = "testData/UCE_aligns_genetrees/UCE_aligns_genetrees_ASTRAL.nwk", annotation = 1, loadOutputTree = TRUE)
plot(allUCEtrees_ASTRALtree, show.node.label = TRUE, main = "ASTRAL Tree: All UCEs")

## Let's see if we can improve those quartet scores by filtering
## TESTS: getPisCounts() and getAlignmentsByPIS()
## Also tests countParsimonyInformative(), which is called inside getPisCounts
pisMatrix <- getPisCounts("testData/UCE_aligns")
# Visualize PIS counts
hist(unlist(pisMatrix[,3]), xlab = "PIS Count", main = "Unfiltered Dataset")
# Filter to just loci above or equal to the 60th percentile and visualize that
pisMatrix_filtered <- getAlignmentsByPIS(pisMatrix, 0.6, "testData/UCE_aligns_PIS60thPercentileAndAbove", returnPisData = TRUE)
hist(unlist(pisMatrix_filtered[,3]), xlab = "PIS Count", main = "60th Percentile and Above")

## Now go run IQTREE on those to get new genetrees (yes, you *could* pull a subset of the 
## genetrees already run on the unfiltered dataset; I have a shell script that did that when
## I was doing this with my own data)
## IQTREE call:
## iqtree -S UCE_aligns_PIS60thPercentileAndAbove/ -T 2 -m TEST --prefix UCE_aligns_PIS60thPercentileAndAbove_genetrees/UCE_aligns_PIS60thPercentileAndAbove_genetrees

## Read those genetrees in and run Astral on them
filteredUCEtrees <- read.tree("testData/UCE_aligns_PIS60thPercentileAndAbove_genetrees/UCE_aligns_PIS60thPercentileAndAbove_genetrees.treefile")
filteredUCEtrees_ASTRALtree <- callAstral(filteredUCEtrees, outfile = "testData/UCE_aligns_PIS60thPercentileAndAbove_genetrees/UCE_aligns_PIS60thPercentileAndAbove_genetrees_ASTRAL.nwk", annotation = 1, loadOutputTree = TRUE)
plot(filteredUCEtrees_ASTRALtree, show.node.label = TRUE, main = "ASTRAL Tree: 60th Percentile and Above UCEs")

