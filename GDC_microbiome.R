
## Clean/reset environment
rm(list = ls()) 

## Load packages
install.packages("phyloseq")
install.packages("microbiome")
install.packages("igraph")
install.packages("qgraph")
install.packages("SpiecEasi")
install.packages("MCL")
install.packages("igraphdata")
install.packages("visNetwork")
install.packages("ggnetwork")




library(phyloseq)
library(microbiome)
library(RColorBrewer)
library(igraph)
library(qgraph)
library(SpiecEasi)
library(vegan)
library(MCL)
library(igraphdata)
library(visNetwork)
library(ggnetwork)


## Create new working directory
dir.create("Food"); setwd("Food")

## Download data (R data image file)
food.url <- "https://www.gdc-docs.ethz.ch/MDA/data/Chaillou2015.zip"
utils::download.file(food.url, destfile = "Chaillou2015.zip")
unzip("Chaillou2015.zip")
file.remove("Chaillou2015.zip")
list.files()
# "chaillou.biom"
# "Chaillou2015b.pdf"
# "otu_table.tsv"
# "sample_data.tsv"
# "sequences.fasta"
# "tax_table.tsv"    
# "tree.nwk"  
setwd("../")

## Function twee - a plain text listing of directories
## similar to tree in linux
dir.create("Scripts"); setwd("Scripts")
twee.url <- "https://www.gdc-docs.ethz.ch/MDA/scripts/twee.R"
utils::download.file(twee.url, destfile = "twee.R")
source("twee.R")
setwd("../")

## Show diretories
twee("MDA/")

## Import from biom
biomfile <- "chaillou.biom"
treefile <- "tree.nwk"
food <- import_biom("C:/Users/michi/Documents/PhD/R tutorials/food/chaillou.biom", 
                    parseFunction = parse_taxonomy_greengenes)
phy_tree(food) <- read_tree("C:/Users/michi/Documents/PhD/R tutorials/food/tree.nwk")
food
## Verify Import
food
# otu_table()   OTU Table:         [ 508 taxa and 64 samples ]
# sample_data() Sample Data:       [ 64 samples by 3 sample variables ]
# tax_table()   Taxonomy Table:    [ 508 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 508 tips and 507 internal nodes ]

## Summary
microbiome::summarize_phyloseq(food)

## Variables
sample_data(food)

## French to English
dictionary = c("BoeufHache"      = "Ground_Beef", 
               "VeauHache"       = "Ground_Veal", 
               "MerguezVolaille" = "Poultry_Sausage", 
               "DesLardons"      = "Bacon_Dice", 
               "SaumonFume"      = "Smoked_Salmon", 
               "FiletSaumon"     = "Salmon_Fillet", 
               "FiletCabillaud"  = "Cod_Fillet", 
               "Crevette"        = "Shrimp")
env_type <- sample_data(food)$EnvType
sample_data(food)$EnvType <- factor(dictionary[env_type], levels = dictionary)

## Add Sample ID
sample_data(food)$SID <- sample_names(food)

## Save / Load
save.image("food.Rdata")


### ============================================= ###
###             Network Graphs Examples           ###
### --------------------------------------------- ###
###                 Version 200106                ###
### --------------------------------------------- ###
###               Jean-Claude Walser              ###
### ============================================= ###

#### Help: terms and definitions  ----

# network theory: study of graphs as a representation of either symmetric relations or asymmetric relations between discrete objects
# graph theory: study of graphs, which are mathematical structures used to model pairwise relations between objects
# vertex (plural vertices) == node -> V
# edge == link -> E
# bridge == isthmus == cut-edge == cut arc is an edge of a graph whose deletion increases its number of connected components
# sparse network: fewer links than the possible maximum number of links within that network
# degree: the number of edges connecting each node to the rest of the network

#### Get Ready ----

### Example Food Data ----

# Load food data
load("food.Rdata")  ## Load food data

## Verify Import
food
# otu_table()   OTU Table:         [ 508 taxa and 64 samples ]
# sample_data() Sample Data:       [ 64 samples by 3 sample variables ]
# tax_table()   Taxonomy Table:    [ 508 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 508 tips and 507 internal nodes ]
class(food)


#### Data Filtering ----

# Prepare data by filtering rare OTUs or standardize data.
# These are fairly arbitary steps. Adjust to you "believes".

## Are there any OTUs in the dataset that have no counted reads?
sum(taxa_sums(food) == 0)

## Save original data
food.original <- food

## Remove "unobserved" OTUs
food <- prune_taxa(taxa_sums(food) > 0, food)
any(taxa_sums(food) == 0)

## OTU counts
food.readsums <- data.frame(nreads = sort(taxa_sums(food), TRUE), sorted = 1:ntaxa(food), type = "OTU")
ggplot(food.readsums, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")

## Remove OTUs with low counts
summary(taxa_sums(food))
food.p <-  prune_taxa(taxa_sums(food) > 50, food)

## Standardize abundances to the median sequencing depth
mrs     <- median(sample_sums(food.p))
standf  <- function(x, t = mrs) round(t * (x / sum(x)))
food.s  <- transform_sample_counts(food.p, standf)

#### Network Analysis ----

### |A| (Dis)-Similarity- or Distance-Based Method ----

# A wide range of methods, with varying levels of efficiency and accuracy,
# have been used to construct networks based on microbiome data.
# The simplest methods are (dis)similarity- or distance-based techniques.

### Netwoks with Phyloseq ----

#  Functions for phyloseq objects
#
#  - [plot_net](https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_net)
#  - [make_network](https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/make_network)
#
#  - [plot_network](https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_network)

## (A) phyloseq::make_network & phyloseq::plot_network

# igraph object
distance.cutoff <- 0.7
food.bc.ig <- make_network(food.s, dist.fun = "bray", max.dist = distance.cutoff)
class(food.bc.ig)

# plot igraph object
plot_network(food.bc.ig, food, color = "EnvType", shape = "FoodType", line_weight = 0.5, label = NULL)

# Identify important nodes unsing different methods

sort(degree(food.bc.ig))
sort(betweenness(food.bc.ig))
sort(closeness(food.bc.ig))
ec <- eigen_centrality(food.bc.ig)
sort(ec$vector)

par(mfrow = c(2,2))
plot(degree(food.bc.ig), typ = "h", col = "tomato", main = "Degree")
plot(betweenness(food.bc.ig), typ = "h", col = "tomato", main = "Betweenness Centrality")
plot(closeness(food.bc.ig), typ = "h", col = "tomato", main = "Closeness Centrality of vertices")
plot(ec$vector, typ = "h", col = "tomato", main = "Eigenvector Centrality Scores")
par(mfrow = c(1,1))

## (B) phyloseq::plot_net

# Note: plot_net function does not require a separate make_network or an igraph object. 

phyloseq::plot_net(food.s, distance = "bray", maxdist = distance.cutoff, point_label = "SID", color = "EnvType", shape = "FoodType", laymeth = "auto")

# Distance Methods:
# unlist(distanceMethodList)

# Note: point_label - the variable name in physeq covariate data to map to vertex labels
# Either create a new sample ID: 
# sample_data(food)$SID <- sample_names(food)
# or use
# point_label = factor(sample_names(food)))

## Compare distance methods

p.wu <- plot_net(food.s, maxdist = distance.cutoff, distance = "wunifrac", color = "EnvType", shape = "FoodType", title = "Weighted-Unifrac")
p.un <- plot_net(food.s, maxdist = distance.cutoff, distance = "unifrac",  color = "EnvType", shape = "FoodType", title = "Unifrac")
p.ja <- plot_net(food.s, maxdist = distance.cutoff, distance = "jaccard",  color = "EnvType", shape = "FoodType", title = "Jaccard")
p.br <- plot_net(food.s, maxdist = distance.cutoff, distance = "bray",     color = "EnvType", shape = "FoodType", title = "Bray-Curtis")
ggpubr::ggarrange(p.wu, p.un, p.ja, p.br)

### Netwoks without Phyloseq ----

## Create dissimilarity matrix
food.bc <- phyloseq::distance(food.s, method = "bray")
min(food.bc); max(food.bc); class(food.bc)

## Explore matrix
boxplot(food.bc, col = "blue", notch = TRUE)
stripchart(food.bc, vertical = TRUE, method = "jitter", add = TRUE, pch = 1, col = "gray")

# Define threshold
distance.cutoff <- 0.7
# Data loss (%)
round(sum(food.bc < distance.cutoff) / length(food.bc) * 100, 1)

## Threshold 
food.bc <- as.matrix(food.bc)
class(food.bc); boxplot(food.bc)

## Convert dissimilarity matrix to adjacency matrix (n x n binary matrix)
food.bc.adj <- ifelse(food.bc <= distance.cutoff, 1, 0)

## Construct network from adjacency matrix
food.bc.adj.ig <- graph.adjacency(food.bc.adj, mode = "undirected", diag = FALSE)
class(food.bc.adj.ig)

## Quick and ugly plot
plot.igraph(food.bc.adj.ig, vertex.color = get_variable(food, "EnvType"))

## A nicer plot

## Fruchterman-Reingold layout algorithm
co <- layout_with_fr(food.bc.adj.ig, niter = 2000)
# niter: Integer scalar, the number of iterations to perform

plot.igraph(food.bc.adj.ig, layout = co, asp = 0, 
            main = "Food: Meat/Seafood", sub = "Bray-Curtis with 0.7 cutoff",
            ## nodes =======================================
            vertex.label  = as.character(get_variable(food.s, "EnvType")),
            vertex.color  = as.factor(get_variable(food.s, "EnvType")),
            vertex.size   = 5,
            ## edges =======================================
            edge.color = "lightgray",
            edge.width = 5,               ## default = 1
            edge.lty = "solid",           ## linetype: blank, solid, dashed, dotted,
            ## dotdash, longdash, or twodash
            edge.curved = 0.05            ## 0 to 1 or TRUE (0.5)
)

## Interactive 2D plot (using the tcltk package)
tkplot(food.bc.adj.net)

## Alternaitve plot (using ggplot)

V(food.bc.adj.ig)$name    <- get_variable(food.s, "SID")
V(food.bc.adj.ig)$envc    <- get_variable(food.s, "EnvType")
V(food.bc.adj.ig)$envn    <- as.character(get_variable(food.s, "EnvType"))
V(food.bc.adj.ig)$foodn   <- get_variable(food.s, "FoodType")
V(food.bc.adj.ig)$foodc   <- as.factor(get_variable(food.s, "FoodType"))

food.bc.adj.ig.gg <- ggnetwork(food.bc.adj.ig, layout = "fruchtermanreingold", niter = 2000)

ggplot(food.bc.adj.ig.gg, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black") +
  geom_nodes(aes(x = x, y = y), size = 5, alpha = 0.8, color = "orange") +
  geom_nodelabel(aes(label = envn), size = 4, color = "blue")


### |B| Correlation- Based Network ----

food.otu <- t(data.frame(phyloseq::otu_table(food.s),
                         check.names = FALSE)
)
class(food.otu)

food.otu.pearson <- cor(food.otu, method = "pearson")
boxplot(abs(food.otu.pearson))

# Convert correlation matrix to binary adjacency matrix
cor.cutoff <- 0.6
food.otu.cor.adj <- ifelse(abs(food.otu.pearson) >= cor.cutoff, 1, 0)

# Construct microbiome network from adjacency matrix 
food.otu.cor.adj.ig <- graph.adjacency(food.otu.cor.adj, mode = "undirected", diag = FALSE)

plot.igraph(food.otu.cor.adj.ig, layout = co, asp = 0, 
            main = "Food: Meat/Seafood", sub = "Correlation: Pearson with 0.6 cutoff",
            ## nodes =======================================
            vertex.label  = as.character(get_variable(food.s, "EnvType")),
            vertex.color  = as.factor(get_variable(food.s, "EnvType")),
            vertex.size   = 5,
            ## edges =======================================
            edge.color = "lightgray",
            edge.width = 5,               ## default = 1
            edge.lty = "solid",           ## linetype: blank, solid, dashed, dotted,
            ## dotdash, longdash, or twodash
            edge.curved = 0.05            ## 0 to 1 or TRUE (0.5)
)

### SparCC Networks ----

# SparCC (Sparse Correlations for Compositional data), is a technique that uses linear correlations 
# between the log-transformed components to infer associations in compositional data (Friedman 2012)

# food.otu <- t(data.frame(phyloseq::otu_table(food.s),) check.names = FALSE)

food.sparcc <- sparcc(food.otu, iter = 20, inner_iter = 10, th = 0.1)
boxplot(food.sparcc)

sparcc.cutoff   <- 0.3
sparcc.adj <- ifelse(abs(food.sparcc$Cor) >= sparcc.cutoff, 1, 0)

## Add OTU names to rows and columns
rownames(sparcc.adj) <- colnames(food.otu)
colnames(sparcc.adj) <- colnames(food.otu)

## Build network from adjacency
sparcc.adj.ig <- graph.adjacency(sparcc.adj, mode = "undirected", diag = FALSE)

plot.igraph(sparcc.adj.ig, layout = co, asp = 0, 
            main = "Food: Meat/Seafood", sub = "Correlation: SPARCC with 0.3 cutoff",
            ## nodes =======================================
            vertex.label  = as.character(get_variable(food.s, "EnvType")),
            vertex.color  = as.factor(get_variable(food.s, "EnvType")),
            vertex.size   = 5,
            ## edges =======================================
            edge.color = "lightgray",
            edge.width = 5,               ## default = 1
            edge.lty = "solid",           ## linetype: blank, solid, dashed, dotted,
            ## dotdash, longdash, or twodash
            edge.curved = 0.05            ## 0 to 1 or TRUE (0.5)
)


### SPIEC-EASI Networks ----

# Sparse InversE Covariance estimation for Ecological Association and Statistical Inference
#
# SPIEC-EASI (SParse InversE Covariance Estima- tion for Ecological Association Inference) is a statistical method 
# for the inference of microbial ecological networks that combines data transformations developed for compositional data 
# analysis with a graphical model inference framework with the assumption that the underlying ecological association network is spars.

# Main SPIEC-EASI pipeline: 
# se <- spiec.easi(X, method = 'mb', lambda.min.ratio = 1e-2, nlambda = 15)
# Data transformation, sparse inverse covariance estimation and model selection

food.se <- spiec.easi(
  food.otu,
  method = "glasso",
  lambda.min.ratio = 1e-2,
  nlambda = 20,
  icov.select.params = list(rep.num = 50)
)

# food.se <- spiec.easi(
#     food.otu,
#     method = 'mb',
#     lambda.min.ratio = 1e-2,
#     nlambda = 20,
#     pulsar.params = list(rep.num = 50)
#   )

# Add OTU names to rows and columns
rownames(SpiecEasi.matrix$refit) <- colnames(amgut1.filt) # Build network from adjacency
SpiecEasi.net <- graph.adjacency(SpiecEasi.matrix$refit,
                                 mode = "undirected", diag = FALSE)

food.se.ig <- adj2igraph(getRefit(food.se), vertex.attr = list(name = taxa_names(food.s)))

plot_network(ig.mb.food, food, type = "taxa", color = "Family")


### Correlation Network (Almost)-All-In-One-Function ----

# Careful: This can get messy. You need a few more packages and not all functions are
#          supported anymore.

source("scripts/CorrelationNetwork.R")

amp_network(food.s, threshold.cor = 0.7, show.label = TRUE, add.tax.info = FALSE, scale.cor = TRUE, show.top = 50)
amp_network(food.s, show.label = TRUE, add.tax.info = FALSE, scale.cor = TRUE, tax.display = "Family")

