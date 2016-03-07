library(limma)
library(edgeR)

manfred <- read.table("e_coli_normalized_counts.out", header = T) # Read Manfred Moose data (4293 genes)
torgny <- read.table("e_coli_counts.txt", header = T) # Read in raw data from Torgny (4336 genes)

# Create DGE lists that are suitable for RNA-seq count data from the count table

DGEmanfred <- DGEList(counts = manfred[,2:10], genes = manfred[,c("gene","reflen")], group = c(rep("A", 3), rep("B", 3), rep("C", 3)))
DGEtorgny <- DGEList(counts = torgny[,2:10], genes = torgny[,c("gene","reflen")], group = c(rep("A", 3), rep("B", 3), rep("C", 3)))

# Evaluate the role of using normfactors in the analysis (might be a bad idea for the already normalised data)

DGEtorgny_norm <- calcNormFactors(DGEtorgny)
DGEmanfred_norm <- calcNormFactors(DGEmanfred)

# Run voom on the objects, needs a designmatrix so first that is created

des.manfred <- model.matrix(~0 + DGEmanfred$sample$group)
des.torgny <- model.matrix(~0 + DGEtorgny$sample$group)
colnames(des.manfred) <- levels(DGEmanfred$sample$group)
colnames(des.torgny) <- levels(DGEtorgny$sample$group)

manfred_voom <- voom(DGEmanfred, des.manfred)
manfred_norm_voom <- voom(DGEmanfred_norm, des.manfred)

torgny_voom <- voom(DGEtorgny, des.torgny)
torgny_norm_voom <- voom(DGEtorgny_norm, des.torgny)

# Test for differential expression
contM <- makeContrasts(B-A, C-A, C-B, levels = des.manfred)

fit.manfred <- lmFit(manfred_voom, des.manfred)
fit.manfred2 <- contrasts.fit(fit.manfred, contM)
fit.manfred2 <- eBayes(fit.manfred2)

fit.manfred_norm <- lmFit(manfred_norm_voom, des.manfred)
fit.manfred_norm2 <- contrasts.fit(fit.manfred_norm, contM)
fit.manfred_norm2 <- eBayes(fit.manfred_norm2)

fit.torgny <- lmFit(torgny_voom, des.torgny)
fit.torgny2 <- contrasts.fit(fit.torgny, contM)
fit.torgny2 <- eBayes(fit.torgny2)

fit.torgny_norm <- lmFit(torgny_norm_voom, des.torgny)
fit.torgny_norm2 <- contrasts.fit(fit.torgny_norm, contM)
fit.torgny_norm2 <- eBayes(fit.torgny_norm2)

# Extract the number of significant genes after fdr correction eg adjusted p-values less than 0.05

summary(res.manfred2 <- decideTests(fit.manfred2))
summary(res.manfred_norm2 <- decideTests(fit.manfred_norm2))
summary(res.torgny2 <- decideTests(fit.torgny2))
summary(res.torgny_norm2 <- decideTests(fit.torgny_norm2))

# Save results of the analysis

write.fit(file = "DE_moose2.txt", results = res.manfred2, digits = 25, sep = "\t", fit = fit.manfred2)
write.fit(file = "DE_moose2_norm.txt", results = res.manfred_norm2, digits = 25, sep = "\t", fit = fit.manfred_norm2)
write.fit(file = "DE_torgny_counts.txt", results = res.torgny2, digits = 25, sep = "\t", fit = fit.torgny2)
write.fit(file = "DE_torgny_counts_norm.txt", results = res.torgny_norm2, digits = 25, sep = "\t", fit = fit.torgny_norm2)
