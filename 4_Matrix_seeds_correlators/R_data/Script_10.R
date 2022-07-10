library(igraph)
library(ape)

setwd("C:/Users/MARIA DEL MAR/Desktop/TFM_WINDOWS/TFM_WINDOWS_1")

# Read matrix of seeds versus correlators
m2 <- as.matrix(read.csv("posneg_correlators_matrix.txt", header = T, sep = '\t', row.names = 1))

################## GeneNetwork of seeds and correlators ########################

# Define seeds and correlators from files
seeds <- readLines("seeds.tsv")
corr_pos <- readLines("correlators_pos.txt")
corr_neg <- readLines("correlators_neg.txt")
corr_posneg <- readLines("correlators_posneg.txt")

# Convert to a bipartitie network
bg = igraph::graph.incidence(m2)
bg

# Label color
V(bg)$label.color <- "white"
for (g in rownames(m2)) {
  edg <- "black"
  try(V(bg)[g]$label.color <- edg, silent = TRUE)
}

# Shape
V(bg)$shape <- "circle"
for (g in rownames(m2)) {
  shap <- "circle"
  try(V(bg)[g]$shape <- shap, silent = TRUE)
}

# Shape color
V(bg)$color <- "black"
for (g in rownames(m2)) {
  col <- "black"
  if (g %in% corr_pos) { col <- "tomato3"}
  if (g %in% corr_neg) { col <- "ivory"}
  if (g %in% corr_posneg) { col <- "moccasin"}
  try(V(bg)[g]$color <- col, silent = TRUE)
}

# Plot
pdf("GeneNetwork.pdf", width=10, height=10)
plot(bg, edge.width = 1, edge.color = "grey", vertex.size = 12, vertex.label.cex = 0.7, 
    asp = 0, cex.main = 5, main = "Seeds and correlators Network", )
legend("bottomleft", title="Type of correlator", c("Positive","Negative","Both"), 
       fill= c("tomato3","ivory","moccasin"), cex=1)
dev.off()

################## GeneNetwork of seeds and correlators ########################

# Transpose matrix
m3 <- t(m2)

# Compute distances and hierarchical clustering. 
# The hclust function in R uses the complete linkage method for hierarchical clustering by default
dd <- dist(m3, method = "manhattan")
hc <- hclust(dd)

# Plot fan
plot(as.phylo(hc), type = "fan")

# Cut the dendrogram into 5 clusters
colors = c("red", "blue", "green", "black", "orange")
clus5 = cutree(hc, 5)

# Plot
pdf("Dendrogram.pdf", width=10, height=10)
par(mar = c(6, 2, 6, 2))
plot(as.phylo(hc), type = "fan", tip.color = colors[clus5],
     label.offset = 0, cex = 1, main = "Dendrogram of seeds", cex.main = 1.5)
dev.off()