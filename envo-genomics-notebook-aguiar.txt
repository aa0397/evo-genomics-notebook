# evo-genomics-notebook-aguiar
# 

#### Initialization & Dependencies ----
{
rm(list=ls())

required_packages <- c("tidyverse","vcfR","adegenet","hierfstat","ggsci","ggplot2")

require2 <- function(packs){
  assign('b',c())
  for (i in seq_along(packs)) {
    assign('b',rbind(b,paste0("require('",packs[[i]],"')")))
  }
  eval(parse(text=b))}

require2(required_packages)
}

####

#### Examples with Visualization ----

data(iris)
head(iris)

#plot of sepal width by sepal length for each species in the data set.
p_iris <- ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species, shape = Species)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw(base_size = 11) +
  labs(
    title = "Sepal Dimensions across Iris Species",
    x = "Sepal Length (cm)",
    y = "Sepal Width (cm)"
  ) +
  theme(panel.grid.minor = element_blank())

print(p_iris)

# Save plot to working directory
# ggsave("iris_sepal_scatter.pdf", p_iris, width = 6, height = 4.5)

####

#### Ho, He, Fis, and PCA in R with Scutellaria saxatilis ----

# Load Metadata
meta <- read.csv("./data/master_metadata_scut_2026_07_27.csv", stringsAsFactors = FALSE)

vcf <- read.vcfR("./Tassel_filtered_2023_12_06_noindels.vcf", verbose = FALSE)

#Convert VCF object to genind structure
genind_obj <- vcfR::vcfR2genind(vcf)

sample_names <- adegenet::indNames(genind_obj)

# Align metadata sample order with VCF genind order
meta <- meta[match(sample_names, meta$name_in_vcf), ]

# Assign population strata
adegenet::strata(genind_obj) <- data.frame(
  State = meta$state,
  Pop = meta$pop,
  Subpop = meta$subpop
)
adegenet::pop(genind_obj) <- meta$state


# Generating the diversity metrics

# Convert genind to hierfstat structure
hf_obj <- hierfstat::genind2hierfstat(genind_obj, pop = meta$state)

# Compute basic statistics (Ho, Hs/He, Fis per locus and population)
pop_stats <- hierfstat::basic.stats(hf_obj)

# Summarize per-population metrics
diversity_summary <- data.frame(
  State = names(colMeans(pop_stats$Ho, na.rm = TRUE)),
  Ho    = colMeans(pop_stats$Ho, na.rm = TRUE),
  He    = colMeans(pop_stats$Hs, na.rm = TRUE),
  Fis   = colMeans(pop_stats$Fis, na.rm = TRUE)
)

print(diversity_summary)

####

#### Step 3: Principal Component Analysis (PCA) & Outlier Filtering ----

# Impute missing locus values with mean allele frequency for initial PCA pass
gen_imputed <- adegenet::tab(genind_obj, NA.method = "mean")

# Pass 1: Perform initial PCA to identify extreme outliers
pca_pass1 <- ade4::dudi.pca(gen_imputed, scannf = FALSE, nf = 2)

# Identify outliers based on PC2 coordinates (PC2 > 25 or PC2 < -10)
outliers <- rownames(pca_pass1$li)[pca_pass1$li$Axis2 > 25 | pca_pass1$li$Axis2 < -10]

# Repeat filtering loop to catch secondary outliers driven by updated axes
genind_clean <- genind_obj[!indNames(genind_obj) %in% outliers, ]
gen_imputed_pass2 <- adegenet::tab(genind_clean, NA.method = "mean")
pca_pass2 <- ade4::dudi.pca(gen_imputed_pass2, scannf = FALSE, nf = 2)

# Identify any additional outlier after re-scaling
outliers_pass2 <- rownames(pca_pass2$li)[pca_pass2$li$Axis2 > 25]
all_outliers <- unique(c(outliers, outliers_pass2))

# Filter genind and metadata objects for final clean run
genind_clean <- genind_obj[!indNames(genind_obj) %in% all_outliers, ]
meta_clean   <- meta[!meta$name_in_vcf %in% all_outliers, ]

# Final PCA on cleaned dataset
gen_imputed_final <- tab(genind_clean, NA.method = "mean")
pca_res <- dudi.pca(gen_imputed_final, scannf = FALSE, nf = 4)

# Calculate variance explained
percent_var <- (pca_res$eig / sum(pca_res$eig)) * 100

# Store clean PC coordinates
meta_clean$PC1 <- pca_res$li$Axis1
meta_clean$PC2 <- pca_res$li$Axis2

# Print removed sample names to console for record-keeping
cat("Removed Outlier Samples:", paste(all_outliers, collapse = ", "), "\n")

####

#### Step 4: Plotting the PCA ----

p_pca <- ggplot(meta_clean, aes(x = PC1, y = PC2, fill = subpop, shape = state)) +
  geom_point(size = 3.5, stroke = 0.3, color = "black") +
  scale_shape_manual(values = c("OH" = 21, "VA" = 22, "WV" = 24)) +
  scale_fill_igv() +
  labs(
    title = "Principal Component Analysis (PCA - Outliers Removed)",
    x = paste0("PC1 (", round(percent_var[1], 1), "% Var)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "% Var)"),
    fill = "Subpopulation",
    shape = "State"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank()) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 3.5, color = "black"))
  )

print(p_pca)

# Save figure output
#ggsave("Figure1_PCA_Ordination_Cleaned.pdf", p_pca, width = 8, height = 6)


# Saving an R session
# save.image(file = "my_workspace.RData")
#loading an r session
# load("my_workspace.RData")

####

#### Questions ---- 

#Notes:
# Fis = (He - Ho)/He.

{"
# Q1: How would you characterize these populations at the state/regional level, in terms of Ho, He, and Fis?

The heterozygosity for WV is considerably higher than VA and OH and is much closer to its expected heterozygosity He of 0.064. The inbreeding Fis is also considrably lower in 
WV than in OH and VA. This is very interesting and suggests that two distinct populations exisist in OH and VA and WV is a region of inbreeding between the two populations 
which increasing the heterzogosity. Alternatively, WV could be a source of two distinct populations which are for some reason expanding in different regions with limited inbreeding. I suppose this is a much less likely scenario. 
This Ho & He values seem low overall, but maybe it is because these are plants and we looked at mammals in previous classes. I would be curious to know if this plant reproduces clonally some times or how many chromosomes it has.

I am curious why the expercte dheterozygosity and observed heterozygosity are so far apart. Perhaps this is an artifact of the founder effect and the limited population sizes as well as
the limited evolutionary time. Maybe if these were bacteria we would see closer concordance between the observed and expected heterozygosity values. 

# Q2: Can you come up with some R code, on your own (use any resources at your disposal), to calculate the ratio of Ho/He? Can you then use an R command to add those values as an extra column in the table below?

diversity_summary$Proportion_Heterozygosity <- diversity_summary$Ho/diversity_summary$He
diversity_summary

# Q3: Examine the PCA output. How would you describe the population structure of S. saxatilis, if there is any? Are there any patterns of distinctness among states or localities (subpopulations)?

We see VA and OH subpopulations separating out along PC2 and WV ppulations separating along PC1. Overall a low amount of variation is explained. Overall we see pretty good
clustering of subpopulations, except for SR-P1 and whatever outliers we removed. We can sort of visualize the cross-breeding individuals at the centroid, there,
where OH and VA populations dip into WV populations around 0,0.

"}


