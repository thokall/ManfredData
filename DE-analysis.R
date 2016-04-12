library(limma)
library(edgeR)

manfred <- read.table("e_coli_normalized_counts.out", header = T) # Read Manfred Moose data (4293 genes)
# torgny <- read.table("e_coli_counts.txt", header = T) # Read in raw data from Torgny (4336 genes)

# Create DGE lists that are suitable for RNA-seq count data from the count table

DGEmanfred <- DGEList(counts = manfred[,2:10], genes = manfred[,c("gene","reflen")], group = c(rep("A", 3), rep("B", 3), rep("C", 3)))
# DGEtorgny <- DGEList(counts = torgny[,2:10], genes = torgny[,c("gene","reflen")], group = c(rep("A", 3), rep("B", 3), rep("C", 3)))

# Evaluate the role of using normfactors in the analysis (might be a bad idea for the already normalised data)

# DGEtorgny_norm <- calcNormFactors(DGEtorgny)
# DGEmanfred_norm <- calcNormFactors(DGEmanfred)

# Make data suitable for linear modelling
# Voom is the preferred method for most data set, estimates a weight for each observation
# that is incorporated in the linear modelling
# The other option is to not use weigths for each observation
# and instead take care of the mean-variance trend at gene level.
# Uses a designmatrix to keep track of samples so first that is created

des.manfred <- model.matrix(~0 + DGEmanfred$sample$group)
# des.torgny <- model.matrix(~0 + DGEtorgny$sample$group)
colnames(des.manfred) <- levels(DGEmanfred$sample$group)
# colnames(des.torgny) <- levels(DGEtorgny$sample$group)
contM <- makeContrasts(B-A, C-A, C-B, levels = des.manfred)

manfred_voom <- voom(DGEmanfred, des.manfred) # typical use for RNA-seq counts
# manfred_voom_noweight <- manfred_voom
# manfred_voom_noweight$weights <- NULL
# manfred_norm_voom <- voom(DGEmanfred_norm, des.manfred)
manfred_no_voom <- cpm(DGEmanfred$counts, log = TRUE)

# torgny_voom <- voom(DGEtorgny, des.torgny)
# torgny_norm_voom <- voom(DGEtorgny_norm, des.torgny)

# Test for differential expression using voom that uses individual weights

fit.manfred <- lmFit(manfred_voom, des.manfred)
fit.manfred2 <- contrasts.fit(fit.manfred, contM)
fit.manfred2 <- eBayes(fit.manfred2)

#fit.manfred.nw <- lmFit(manfred_voom_noweight, des.manfred)
#fit.manfred2.nw <- contrasts.fit(fit.manfred.nw, contM)
#fit.manfred2.nw <- eBayes(fit.manfred2.nw, trend = TRUE)


#fit.manfred_norm <- lmFit(manfred_norm_voom, des.manfred)
#fit.manfred_norm2 <- contrasts.fit(fit.manfred_norm, contM)
#fit.manfred_norm2 <- eBayes(fit.manfred_norm2)

#fit.torgny <- lmFit(torgny_voom, des.torgny)
#fit.torgny2 <- contrasts.fit(fit.torgny, contM)
#fit.torgny2 <- eBayes(fit.torgny2)

#fit.torgny_norm <- lmFit(torgny_norm_voom, des.torgny)
#fit.torgny_norm2 <- contrasts.fit(fit.torgny_norm, contM)
#fit.torgny_norm2 <- eBayes(fit.torgny_norm2)

# Test for differential expression using Limma trend
fit.manfred.LT <- lmFit(manfred_no_voom, des.manfred)
fit.manfred.LT2 <- contrasts.fit(fit.manfred.LT, contM)
fit.manfred.LT2 <- eBayes(fit.manfred.LT2, trend = TRUE)


# Extract the number of significant genes after fdr correction eg adjusted p-values less than 0.05

summary(res.manfred2 <- decideTests(fit.manfred2))
summary(res.manfred.LT2 <- decideTests(fit.manfred.LT2))
# summary(res.manfred.nw <- decideTests(fit.manfred2.nw))
# summary(res.torgny2 <- decideTests(fit.torgny2))
# summary(res.torgny_norm2 <- decideTests(fit.torgny_norm2))



# Save results of the analysis

write.fit(file = "DE_moose2voom.txt", results = res.manfred2, digits = 25, sep = "\t", fit = fit.manfred2)
write.fit(file = "DE_moose2_novoom.txt", results = res.manfred.LT2, digits = 25, sep = "\t", fit = fit.manfred.LT2)

# write.fit(file = "DE_moose2_norm.txt", results = res.manfred_norm2, digits = 25, sep = "\t", fit = fit.manfred_norm2)
# write.fit(file = "DE_torgny_counts.txt", results = res.torgny2, digits = 25, sep = "\t", fit = fit.torgny2)
# write.fit(file = "DE_torgny_counts_norm.txt", results = res.torgny_norm2, digits = 25, sep = "\t", fit = fit.torgny_norm2)

