library(vcfR)
library(adegenet)

vcf <- read.vcfR("olfusa.filtered.vcf")

x <- vcfR2genind(vcf)

# Dealing with NA values
y <- tab(x, NA.method="mean")

# Performing a PCA
y.pca = dudi.pca(y, scannf=FALSE, nf=20)

library(factoextra)
library(ggpubr)
library(ggplot2)

fviz_eig(y.pca)

scree_plot <- fviz_eig(y.pca, addlabels = TRUE,
                       barfill = "#fdae61",
                       barcolor = "#f46d43",
                       linecolor = "black") + 
  theme(text = element_text(family="Times New Roman")) + labs(title = NULL)


pop <- read.table(file="popmap.olfusa.clean.tsv", header=F)
groups <- as.factor(pop$V2)

## with dim1 vs dim2
fviz_pca_ind(y.pca,
             habillage = groups,
             axes = c(1,2),
             geom = c("point"), 
             mean.point = F,
             repel = F,
             addEllipses = F,
             label="none",
             pointsize=5,
             palette=c("#fee090","#a50026", "#fee090", "#f46d43", "#a50026",
                       "#f46d43", "#74add1", "#a50026", "#f46d43", "#fee090", 
                       "#313695", "#f46d43", "#fee090", "#fee090", "#fee090", 
                        "#a50026", "#fee090")) +
scale_shape_manual(values=c(23,22,4,22,19,19,19,24,24,21,19,21,22,19,24,21,25)) + 
labs(title = NULL, x = "PC1 (20.3%)", y = "PC2 (8.2%)") + 
theme(text = element_text(family="Times New Roman"))


## with dim3 vs dim4
fviz_pca_ind(y.pca,
             habillage = groups,
             axes = c(3,4),
             geom = c("point"), 
             mean.point = F,
             repel = F,
             addEllipses = F,
             label="none",
             pointsize=5,
             palette=c("#fee090","#a50026", "#fee090", "#f46d43", "#a50026",
                       "#f46d43", "#74add1", "#a50026", "#f46d43", "#fee090", 
                       "#313695", "#f46d43", "#fee090", "#fee090", "#fee090", 
                       "#a50026", "#fee090")) +
scale_shape_manual(values=c(23,22,4,22,19,19,19,24,24,21,19,21,22,19,24,21,25)) + 
  labs(title = NULL, x = "PC3 (4.0%)", y = "PC4 (3.5%)") + 
  theme(text = element_text(family="Times New Roman"))


## with dim5 vs dim6
fviz_pca_ind(y.pca,
             habillage = groups,
             axes = c(5,6),
             geom = c("point"), 
             mean.point = F,
             repel = F,
             addEllipses = F,
             label="none",
             pointsize=5,
             palette=c("#fee090","#a50026", "#fee090", "#f46d43", "#a50026",
                       "#f46d43", "#74add1", "#a50026", "#f46d43", "#fee090", 
                       "#313695", "#f46d43", "#fee090", "#fee090", "#fee090", 
                       "#a50026", "#fee090")) +
scale_shape_manual(values=c(23,22,4,22,19,19,19,24,24,21,19,21,22,19,24,21,25)) + 
  labs(title = NULL, x = "PC5 (2.7%)", y = "PC6 (2.3%)") + 
  theme(text = element_text(family="Times New Roman"))