#!/usr/bin/env Rscript
#Script to generate Figure 3
#Loading the required packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(wesanderson))
suppressPackageStartupMessages(library(taxize))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggVennDiagram))
suppressPackageStartupMessages(library(patchwork))

#Importing sponge scaffold contamination results from file
sponge_scaff_res_dt <- fread("scaffold_contamination_info_sponges.tsv")
sponge_scaff_res_dt[, "species" := sub("_scaffold_contamination_info.tsv", "", filename, fixed = T)]
sponge_scaff_res_dt[, "species" := sub("_", " ", species, fixed = T)]

#Calculating the proportion of contaminant scaffolds per assembly
sponge_scaff_contam_N <- sponge_scaff_res_dt[, .("scaff_N" = .N), by = c("contamination", "species")]
sponge_scaff_contam_N[, "scaff_prop" := scaff_N / sum(scaff_N), by = "species"]
species_fact_levels <- sponge_scaff_contam_N[contamination == F][order(-scaff_prop), species]
sponge_scaff_contam_N[, "species" := factor(species, levels = species_fact_levels)]

bar_colors <- c("#4FA3BCFF", "#F2AD00")
bar_colors <- c(alpha(bar_colors[1], 0.9), alpha(bar_colors[2], 0.8))

#Plotting the percentage of contaminant scaffolds per assembly (Figure 3a)
fig_3a <- ggplot(data = sponge_scaff_contam_N, aes(x = species, y = scaff_N, fill = contamination)) +
 geom_bar(stat = "identity", position = "fill", width = 0.8, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = bar_colors, labels = c("Sponge", "Contamination")) +
 theme(legend.position = "right", legend.margin = margin(0, 0, 0, 0)) +
 theme(legend.title = element_blank()) +
 theme(legend.text = element_text(size = 5)) +
 theme(legend.key.size = unit(0.27, 'cm')) +
 theme(panel.grid = element_blank()) +  
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 ylab("% scaffolds") +
 scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "20", "40", "60", "80", "100")) +
 theme(axis.title.x = element_blank()) +
 theme(axis.title.y = element_text(size = 7.5)) +
 theme(axis.text.x = element_text(size = 5.6, angle = 55, hjust = 0.9, face = "italic")) +
 theme(axis.text.y = element_text(size = 4.5)) 

ggsave(fig_3a, file = "fig_3a.svg", height = 2.5, width = 6)

#Calculating the proportion of contaminant genome length per assembly
sponge_scaff_contam_len <- sponge_scaff_res_dt[, .("scaff_len_sum" = sum(scaffold_len)), by = c("contamination", "species")]
sponge_scaff_contam_len[, "species" := factor(species, levels = species_fact_levels)]

#Plotting the percentage of contaminant genome length per assembly (Figure 3b)
fig_3b <- ggplot(data = sponge_scaff_contam_len, aes(x = species, y = scaff_len_sum, fill = contamination)) +
 geom_bar(stat = "identity", position = "fill", width = 0.8, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = bar_colors, labels = c("Sponge", "Contamination")) +
 theme(legend.position = "right", legend.margin = margin(0, 0, 0, 0)) +
 theme(legend.title = element_blank()) +
 theme(legend.text = element_text(size = 5)) +
 theme(legend.key.size = unit(0.27, 'cm')) +
 theme(panel.grid = element_blank()) +  
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 ylab("% scaffold length") +
 scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "20", "40", "60", "80", "100")) +
 theme(axis.title.x = element_blank()) +
 theme(axis.title.y = element_text(size = 7.5)) +
 theme(axis.text.x = element_text(size = 5.6, angle = 55, hjust = 0.9, face = "italic")) +
 theme(axis.text.y = element_text(size = 4.5)) 

ggsave(fig_3b, file = "fig_3b.svg", height = 2.5, width = 6)

#Analysing contaminant scaffold lengths
sponge_scaff_res_dt_cont <- sponge_scaff_res_dt[contamination == T]
sponges_wo_cont <- sponge_scaff_contam_N[contamination == F & scaff_prop == 1, species]
sponge_scaff_res_dt_cont <- rbindlist(list(sponge_scaff_res_dt_cont, sponge_scaff_res_dt[species %in% sponges_wo_cont]))
sponge_scaff_res_dt_cont[species %in% sponges_wo_cont, "scaffold_len" := 0]
sponge_scaff_res_dt_cont[, "species" := factor(species, levels = species_fact_levels)]

#Plotting contaminant scaffold length distribution (Figure 3c)
fig_3c <- ggplot(data = sponge_scaff_res_dt_cont, aes(x = species, y = scaffold_len, fill = species)) +
 geom_boxplot(width = 0.8, size = 0.2, outlier.shape = NA, color = "gray20", fill = "#4FA3BCFF") +
 theme_minimal() +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +  
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 scale_y_continuous(breaks = c(1000, 50000, 100000, 150000, 200000), labels = c(1, 50, 100, 150, 200), limits = c(0, max(sponge_scaff_res_dt_cont[, scaffold_len]))) +
 ylab("Scaffold length (kb)") +
 theme(axis.title.x = element_blank()) +
 theme(axis.title.y = element_text(size = 7.5)) +
 theme(axis.text.x = element_text(size = 5.6, angle = 55, hjust = 0.9, face = "italic")) +
 theme(axis.text.y = element_text(size = 4.5)) 
  
ggsave(fig_3c, file = "fig_3c.svg", height = 2.5, width = 5.25)

#Importing contaminant scaffolds from file
contam_scaff_taxonomy <- fread("contaminant_scaffolds_taxonomy_sponges.tsv")
contam_scaff_taxonomy[, "species" := sub("_contaminant_scaffolds_taxonomy.tsv", "", filename, fixed = T)]
contam_scaff_taxonomy[, "species" := sub("_", " ", species)]

#Assigning higher taxonomic categories to contaminant taxa
taxa_to_assign_higher_tax <- unique(contam_scaff_taxonomy[is.na(assigned_tax) == F, assigned_tax])
taxa_to_assign_higher_tax_classfied <- classification(taxa_to_assign_higher_tax, db = "ncbi", rows = 1)

taxa_to_assign_higher_tax_domain <- sapply(taxa_to_assign_higher_tax_classfied, function(x) {
 if (is.data.frame(x)) {
  val <- x$name[x$rank == "domain"]
  if (length(val) == 0) NA else val
  } else {
  NA
 }
})

taxa_to_assign_higher_tax_phylum <- sapply(taxa_to_assign_higher_tax_classfied, function(x) {
 if (is.data.frame(x)) {
  val <- x$name[x$rank == "phylum"]
  if (length(val) == 0) NA else val
  } else {
  NA
 }
})

taxa_to_assign_higher_tax <- as.data.table(cbind(taxa_to_assign_higher_tax, taxa_to_assign_higher_tax_domain, taxa_to_assign_higher_tax_phylum))
setnames(taxa_to_assign_higher_tax, c("assigned_tax", "domain", "phylum"))
contam_scaff_taxonomy <- merge(contam_scaff_taxonomy, taxa_to_assign_higher_tax, by = "assigned_tax", all.x = T)
contam_scaff_taxonomy <- contam_scaff_taxonomy[is.na(phylum) == F]
contam_scaff_taxonomy <- contam_scaff_taxonomy[phylum %chin% c("Porifera", "Platyhelminthes", "Mollusca", "Arthropoda", "Echinodermata", "Chordata") == F]
contam_scaff_taxonomy[, "phylum" := sub("Candidatus ", "", phylum, fixed = T)]

#Calculating the proportion of contaminant scaffold phyla across sponge assemblies
contam_scaff_taxonomy_N <- contam_scaff_taxonomy[, .("scaff_N" = .N), by = c("phylum", "species")]
non_contam_scaffold_N <- sponge_scaff_contam_N[contamination == F, .(species, scaff_N)]
non_contam_scaffold_N[, "phylum" := "Non-contaminants"]
setcolorder(non_contam_scaffold_N, c(3, 1, 2))
contam_scaff_taxonomy_N <- rbindlist(list(contam_scaff_taxonomy_N, non_contam_scaffold_N))
contam_scaff_taxonomy_N[, "scaff_prop" := scaff_N / sum(scaff_N), by = "species"]
contam_scaff_taxonomy_N <- contam_scaff_taxonomy_N[phylum != "Non-contaminants"]

phyla_to_keep <- contam_scaff_taxonomy_N[, .("species_N" = uniqueN(species)), by = "phylum"][species_N > 1, phylum]
contam_scaff_taxonomy_N <- contam_scaff_taxonomy_N[phylum %chin% phyla_to_keep]

domain_phylum_dt <- unique(contam_scaff_taxonomy[, .(domain, phylum)])
setnames(domain_phylum_dt, c("Domain", "phylum"))
contam_scaff_taxonomy_N <- dcast(contam_scaff_taxonomy_N, species ~ phylum, value.var = "scaff_prop", fun.aggregate = sum, fill = 0)

contam_scaff_taxonomy_N_mat <- as.matrix(contam_scaff_taxonomy_N[, -1])
contam_scaff_taxonomy_N_mat <- log10(contam_scaff_taxonomy_N_mat + 1e-4)
rownames(contam_scaff_taxonomy_N_mat) <- contam_scaff_taxonomy_N[, species]

anno_col <- data.frame(Domain = domain_phylum_dt$Domain[match(colnames(contam_scaff_taxonomy_N_mat), domain_phylum_dt$phylum)], row.names = colnames(contam_scaff_taxonomy_N_mat))
domain_levels <- sort(unique(anno_col$Domain))
domain_colors <- setNames(wes_palette("Darjeeling1", 5, type = "discrete")[c(3, 2, 5)][seq_along(domain_levels)], domain_levels)
domain_colors <- list(Domain = domain_colors)
heatmap_colors <- c("#4FA3BCFF", "#77B7CFFF", "#9EC3A9FF", "#CFCB63FF", "#E3AE12FF", "#F06C0BFF")

#Plotting the proportion of contaminant scaffold phyla across sponge assemblies (Figure 3d)
pdf(file = "fig3d.pdf", width = 7.1, height = 4)
pheatmap(contam_scaff_taxonomy_N_mat, cluster_rows = T, cluster_cols = T, color = colorRampPalette(heatmap_colors)(100), border_color = NA, fontsize = 7, fontsize_row = 5.6, fontsize_col = 5.1, angle_col = 45, legend_breaks = c(-4, -3, -2, -1, -0.7), legend_labels = c("0  ","0.1 ", "1  ","10  ", "% scaffolds  "), annotation_col = anno_col, annotation_colors = domain_colors)
dev.off()

#Plotting Venn diagram of compositional, protein-based, and nucleotide-based scores per assembly (Figure 3e)
scores_venn <- function(sponge_scaff_res_dt, species) {
 sets <- list(C = sponge_scaff_res_dt[compositional_score == 1, scaffold_name], P = sponge_scaff_res_dt[prot_taxonomy_score == 1, scaffold_name], N  = sponge_scaff_res_dt[nucl_taxonomy_score == 1, scaffold_name])

 ggVennDiagram(sets,label = "count", label_alpha = 0, set_size = 1.9, label_size = 2.2, edge_size = 0.2) +
  theme_void() +
  scale_fill_gradient(low = "#F7FBFF", high = "#4FA3BCFF") +
  theme(legend.position = "none") +
  theme(plot.margin = margin(2, 2, 2, 2)) +
  labs(title = species) +
  theme(plot.title = element_text(hjust = 0.5, face = "italic", size = 5.2)) 
}

venn_list <- lapply(species_fact_levels, function(sp) {
 sponge_scaff_res_dt_species <- sponge_scaff_res_dt[species == sp] 
 scores_venn(sponge_scaff_res_dt_species, sp)})

fig_3e <- wrap_plots(plotlist = venn_list, nrow = 5) 
ggsave(fig_3e, file = "fig_3e.svg", height = 5.5, width = 7.9)

#Importing UMAP projection coordinates of k-mer-based PCs from file
umap_res_dt <- fread("mc_umap_dt_sponges.tsv")
umap_res_dt[, "species" := sub("_mc_umap_dt.tsv", "", filename, fixed = T)]
umap_res_dt[, "species" := sub("_", " ", species, fixed = T)]
umap_res_dt <- merge(umap_res_dt, sponge_scaff_res_dt[, .(scaffold_name, species, prot_taxonomy_score, nucl_taxonomy_score)], by = c("scaffold_name", "species"))
umap_res_dt[, "prot_nucl_tax_score" := paste(prot_taxonomy_score, nucl_taxonomy_score)]
umap_res_dt[, "species" := factor(species, levels = species_fact_levels)]

#Plotting UMAP projection coordinates of k-mer-based PCs from file (Figure 3f)
fig_3f <-  ggplot(data = umap_res_dt, aes(UMAP1, UMAP2, color = prot_nucl_tax_score)) +
 geom_point(size = 0.003, shape = 16, alpha = 0.75) +
 theme_minimal() +
 scale_color_manual(values = wes_palette("Darjeeling1", 5, type = "discrete")[c(5, 2, 3, 4)], name = "Taxonomy scores", labels = c("None", expression(italic("P")), expression(italic("N")), expression(italic("P") ~ "+" ~ italic("N")))) +
 theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.05, "cm")) +      
 theme(legend.title = element_text(size = 5.2)) +
 theme(legend.text = element_text(size = 4.2)) +
 guides(color = guide_legend(override.aes = list(size = 1, alpha = 1))) +
 facet_wrap(. ~ species, scales = "free") +
 theme(strip.text = element_text(size = 5, face = "italic")) +
 theme(panel.grid = element_blank())  +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 theme(axis.title = element_text(size = 6)) +
 theme(axis.text = element_text(size = 4.3))

ggsave(fig_3f, file = "fig_3f.svg", height = 5.3, width = 6.1)

#Importing BUSCO results from file
orig_genome_busco_res <- fread("busco_original_genome_complete_merged.tsv")
orig_genome_busco_res[, "category" := "Original"]
decom_genome_busco_res <- fread("busco_decontaminated_genome_complete_merged.tsv")
decom_genome_busco_res[, "category" := "Decontaminated"]

#Merging BUSCO results
busco_res <- rbindlist(list(orig_genome_busco_res, decom_genome_busco_res))
setnames(busco_res, c("N_BUSCO", "species", "genome_category"))
busco_res[, "species" := sub("^(([^_]*_){1}[^_]*)_.*$", "\\1", species)]
busco_res[, "species" := sub("_", " ", species, fixed = T)]
busco_res[, "species" := factor(species, levels = species_fact_levels)]
busco_res[, "genome_category" := factor(genome_category, levels = c("Original", "Decontaminated"))] 

#Plotting the number of complete BUSCOs in original and decontaminated assemblies (Figure 3g)
fig_3g <- ggplot(data = busco_res,aes(x = species, y = N_BUSCO, color = genome_category)) +
 geom_point(position = position_dodge(width = 0), size = 1, shape = 16, alpha = 0.7) +
 theme_minimal() +
 scale_color_manual(values = c("#F06C0BFF", "#4FA3BCFF")) +
 theme(legend.position = "right", legend.justification = "center", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-6, -6, -3, -6)) +
 theme(legend.title = element_blank()) +
 theme(legend.text = element_text(size = 4.3)) +
 theme(legend.key.size = unit(0.3, 'cm')) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 ylim(c(0, 672)) +
 ylab("N complete BUSCOs") +
 theme(axis.title.x = element_blank()) +
 theme(axis.title.y = element_text(size = 7.5)) +
 theme(axis.text.x = element_text(size = 5.6, angle = 55, hjust = 0.9, face = "italic")) +
 theme(axis.text.y = element_text(size = 4.5)) 

ggsave(fig_3g, file = "fig_3g.svg", height = 2.5, width = 6)
