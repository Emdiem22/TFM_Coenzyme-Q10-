library(viridis)
library(here)        
library(readr)      
library(dplyr)      
library(tidyr)       
library(tibble)      
library(ggplot2) 

setwd("C:/Users/MARIA DEL MAR/Desktop/TFM_WINDOWS/TFM_WINDOWS_1")

############ Barplot of number of correlators per seed ################

# Read the data frame of seeds versus seeds as matrix
m <- as.matrix(read.csv("posneg_seeds_matrix.txt", header = T, sep = '\t', row.names = 1))

# Define variables of the data frame
gene <- c("COQ2", "COQ4", "COQ5", "COQ6", "COQ7", "COQ8a", "COQ8b", "COQ9", "PDSS1", "PDSS2")
correlators <- c(m[1,1],m[2,2],m[3,3],m[4,4],m[5,5],m[6,6],m[7,7],m[8,8],m[9,9],m[10,10])
df <- data.frame(gene, correlators)
col <- plasma(10)

# Plot
pdf("Barplot.pdf", width=10, height=10)
barplot(df$correlators, main = "Number of correlators per seed", cex.main = 1.6, axis.lty=1, 
        col = col, ylim=c(0,500), xlab = "Seed", names = df$gene, cex.names = 0.7, 
                                ylab = "Number of correlators")
dev.off()

############ Heatmap of common correlators between seeds ################

# Read the data frame of seeds versus seeds
(dat <- readr::read_delim(here::here("posneg_seeds_matrix.txt")))

# Convert the data frame to a matrix
(mat <- dat |> column_to_rownames("Identifier") |> as.matrix())

# Remove half of the matrix
mat[upper.tri(mat, diag = TRUE)] <- NA

# Convert it back to a data frame
(mat_df <- as.data.frame(mat) |> rownames_to_column("Gene1"))

# Pivote it to have one column for each axis of the plot, and remove the NA:
(mat_df <- mat_df |> pivot_longer(-Gene1, names_to = "Gene2") |> drop_na(value))

# Plot
pdf("Heatmap.pdf", width=10, height=10)
ggplot(mat_df, aes(x = Gene1, y = Gene2)) +
  geom_tile(aes(fill = value), colour = "white", size = rel(1.2)) + 
  geom_text(aes(label = value), size = rel(4.5)) +
  scale_fill_gradient(low = "white", high = "firebrick4", 
    na.value = "white", trans = "pseudo_log") +
  scale_y_discrete(position = "right") +
  labs(title = "Number of common correlators between seeds", x = "", y = "") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20), legend.position = "none", panel.grid.major = element_blank())
dev.off()