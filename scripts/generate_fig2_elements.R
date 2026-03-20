#!/usr/bin/env Rscript
#Script to generate Figure 2
#Loading the required packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(wesanderson))
suppressPackageStartupMessages(library(ComplexUpset))
suppressPackageStartupMessages(library(pheatmap))

#Importing validation dataset scaffold contamination results from file
validation_dataset_decon_results <- fread("scaffold_contamination_info_validation.tsv")
validation_dataset_decon_results[, "species" := sub(".*\\.(?:[1-9]|1[0-2])_", "", sub("_plasmid.*", "", sub("_longest", "", sub("_mito", "", sub("\\.fa.*", "", scaffold_name)))))]

#Calculating decontamination pipeline performance metrics
validation_dataset_decon_results[, "contamination" := as.numeric(contamination)]
validation_dataset_decon_results[, "contamination_truth" := fifelse(grepl("contaminant_scaffold", scaffold_name, fixed = T), 1, 0)]

calculate_perf_metrics <- function(input_dt, truth_col, pred_col) {
 TN <- sum(input_dt[[pred_col]] == 0 & input_dt[[truth_col]] == 0)
 TP <- sum(input_dt[[pred_col]] == 1 & input_dt[[truth_col]] == 1)
 FN <- sum(input_dt[[pred_col]] == 0 & input_dt[[truth_col]] == 1)
 FP <- sum(input_dt[[pred_col]] == 1 & input_dt[[truth_col]] == 0)

 precision <- TP / (TP + FP)
 recall    <- TP / (TP + FN)
 accuracy  <- (TP + TN) / (TP + TN + FP + FN)

 data.table(predictor = pred_col, precision = precision, recall = recall, accuracy = accuracy)
}

predictor_cols <- c("compositional_score",  "prot_taxonomy_score", "nucl_taxonomy_score", "contamination")
perf_metrics <- rbindlist(lapply(predictor_cols, function(col) calculate_perf_metrics(validation_dataset_decon_results, truth_col = "contamination_truth", pred_col = col)))
perf_metrics <- melt(perf_metrics, id.vars = "predictor", measure.vars = c("precision", "recall", "accuracy"), variable.name = "metric", value.name = "value")
perf_metrics[, "metric" := factor(metric, levels = c("accuracy", "precision", "recall"))]
perf_metrics[, "predictor" := factor(predictor, levels = c("compositional_score", "prot_taxonomy_score", "nucl_taxonomy_score", "contamination"))]

#Plotting accuracy, precision, and recall for each individual scoring component and for the integrated contamination score (Figure 2a)
fig_2a <- ggplot(data = perf_metrics, aes(x = predictor, y = value, fill = metric)) +
 geom_bar(stat = "identity", position = "dodge", width = 0.8, color = "gray40", alpha = 0.7) +
 theme_minimal() +
 scale_fill_manual(values = wes_palette("Darjeeling1", 5, type = "discrete")[c(1, 3, 5)], labels = c("Accuracy", "Precision", "Recall")) +
 theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0)) +
 theme(legend.title = element_blank()) +
 theme(legend.text = element_text(size = 5))  +
 theme(legend.key.size = unit(0.27, 'cm')) +
 theme(panel.grid = element_blank()) +  
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Score") +
 scale_x_discrete(labels = c("C", "P", "N", "S")) +
 scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
 theme(axis.title = element_blank()) +
 theme(axis.text.x = element_text(size = 6, face = "italic")) +
 theme(axis.text.y = element_text(size = 6)) 

ggsave(fig_2a, file = "fig_2a.svg", height = 1.8, width = 3.5)

#Preparing for intersection analysis
validation_dataset_decon_results[contamination_truth == 0 & contamination == 0, "class_result" := "TN"]
validation_dataset_decon_results[contamination_truth == 1 & contamination == 1, "class_result" := "TP"]
validation_dataset_decon_results[contamination_truth == 1 & contamination == 0, "class_result" := "FN"]
validation_dataset_decon_results[contamination_truth == 0 & contamination == 1, "class_result" := "FP"]

upset_dt <- validation_dataset_decon_results[, .("C" = compositional_score == 1, "P" = prot_taxonomy_score == 1, "N" = nucl_taxonomy_score == 1, "Classification" = class_result, "scaffold_length" = scaffold_len)]
upset_dt[, "Classification" := factor(Classification, levels = c("TN", "TP", "FN", "FP"))]

#Plotting UpSet plot of scoring component results (Figure 2b)
fig_1b <- upset(upset_dt, intersect = c("C", "P", "N"), set_sizes = F, height_ratio = 0.4, stripes = 'white', 
 themes = upset_modify_themes(list("intersections_matrix" = theme(axis.text.y = element_text(size = 6.5, face = "italic"), axis.title.x = element_blank()))), 
 base_annotations = list("Intersection size" = intersection_size(mapping = aes(fill = Classification), text = list(size = 0), alpha = 0.7) +
  scale_fill_manual(values = wes_palette("Darjeeling1", 5, type = "discrete")[c(5, 4, 2, 1)]) +
  theme(legend.position = "top", legend.direction = "horizontal", legend.margin = margin(0, 0, 0, 0)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 6)) +
  theme(legend.key.size = unit(0.28, 'cm')) +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_rect(fill = NA, color = "grey70")) +
  theme(axis.title.y = element_text(size = 7)) +
  theme(axis.text.y = element_text(size = 6))), 
 annotations = list("Scaffold length" = (
  ggplot(mapping = aes(y = log10(scaffold_length))) + 
   geom_boxplot(width = 0.8, size = 0.2, outlier.shape = NA, color = "gray40", fill = "#5BBCD6", alpha = 0.7) +
   ylab("log(sequence length)") +
   theme(panel.grid = element_blank() +
   theme(panel.background = element_rect(fill = NA, color = "grey70") +
   theme(axis.title.y = element_text(size = 6)) +
   theme(axis.text.y = element_text(size = 5.5)))))))

ggsave(fig_1b, file = "fig_1b.svg", height = 3.2, width = 3.5)

#Analysing recall per phylum
external_contaminant_taxonomy_dt <- fread("taxonomy_info_validation.csv")
external_contaminant_taxonomy_dt[external_contaminant_taxonomy_dt == ""] <- NA
validation_dataset_decon_results <- merge(validation_dataset_decon_results, external_contaminant_taxonomy_dt, by = "species", all.x = T)

validation_dataset_decon_results_cont <- validation_dataset_decon_results[contamination_truth == 1]
validation_dataset_decon_results_cont <- melt(validation_dataset_decon_results_cont, measure.vars = predictor_cols, variable.name = "score_type", value.name = "prediction")
validation_dataset_decon_results_cont[, "class_result" := fifelse(prediction == 1, "TP", "FN")]
validation_dataset_decon_results_cont <- validation_dataset_decon_results_cont[, .N, by = c("score_type", "phylum", "domain", "class_result")]
validation_dataset_decon_results_cont <- dcast(validation_dataset_decon_results_cont, score_type + phylum + domain ~ class_result, value.var = "N", fill = 0)

validation_dataset_decon_results_cont <- validation_dataset_decon_results_cont[, .("recall" = TP / (TP + FN)), by = c("score_type", "phylum", "domain")]
validation_dataset_decon_results_cont <- dcast(validation_dataset_decon_results_cont, score_type ~ phylum, value.var = "recall")
validation_dataset_decon_results_cont <- as.matrix(validation_dataset_decon_results_cont[, -1])
rownames(validation_dataset_decon_results_cont) <- c("C", "P", "N", "S")

phylum_domain_dt <- unique(validation_dataset_decon_results[, .(phylum, domain)])
ann_col <- as.data.frame(phylum_domain_dt[match(colnames(validation_dataset_decon_results_cont), phylum), .(domain)])
colnames(ann_col) <- "Domain"
rownames(ann_col) <- colnames(validation_dataset_decon_results_cont)

domain_colors <-  wes_palette("Darjeeling1", 5, type = "discrete")[c(1, 2, 3, 5)]
names(domain_colors) <- c("Viruses", "Bacteria", "Archaea", "Eukaryota")
domain_colors <- list(Domain = domain_colors)
heatmap_colors <- c("#4FA3BCFF", "#77B7CFFF", "#9EC3A9FF", "#CFCB63FF", "#E3AE12FF", "#F06C0BFF")

#Plotting recall for each individual scoring component and for the integrated contamination score per phylum (Figure 2c)
pdf(file = "fig_2c.pdf", width = 5.5, height = 4)
pheatmap(validation_dataset_decon_results_cont, cluster_rows = T, cluster_cols = T, color = heatmap_colors, border_color = "gray20", cellheight = 12, treeheight_row = 20, treeheight_col = 25, fontsize = 5.5, legend_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.95, 1), legend_labels = c("0  ", "0.2  ", "0.4  ", "0.6  ", "0.8  ", "1  ", "Recall  "), annotation_col = ann_col, annotation_colors = domain_colors)
dev.off()

#Importing contaminant scaffold taxonomy from file
contam_scaff_taxonomy <- fread("contaminant_scaffolds_taxonomy_validation.tsv")
contam_scaff_taxonomy[, "species_true" := sub(".fa.*", "", scaffold_name)]
contam_scaff_taxonomy[, "species_true" := sub("_mito", "", species_true, fixed = T)]
contam_scaff_taxonomy[, "species_true" := sub("_longest", "", species_true, fixed = T)]
contam_scaff_taxonomy[, "species_true" := sub("_plasmid.*", "", species_true)]
contam_scaff_taxonomy[, "species_true" := sub(".*\\.(?:[1-9]|1[0-2])_", "", species_true)]

colnames(external_contaminant_taxonomy_dt) <- paste(colnames(external_contaminant_taxonomy_dt), "true", sep = "_")
contam_scaff_taxonomy <- merge(contam_scaff_taxonomy, external_contaminant_taxonomy_dt, by = "species_true", all.x = T)
contam_scaff_taxonomy <- contam_scaff_taxonomy[phylum_true != "Porifera"]
contam_scaff_taxonomy[, "species_true" := sub("_", " ", species_true, fixed = T)]

#Analysing contaminant scaffold taxonomy calls
contam_scaff_taxonomy[assigned_tax == genus_true, "assignment" := "Genus"]
contam_scaff_taxonomy[assigned_tax == family_true, "assignment" := "Family"]
contam_scaff_taxonomy[assigned_tax == order_true, "assignment" := "Order"]
contam_scaff_taxonomy[assigned_tax == class_true, "assignment" := "Class"]
contam_scaff_taxonomy[assigned_tax == phylum_true, "assignment" := "Phylum"]
contam_scaff_taxonomy[assigned_tax == domain_true, "assignment" := "Domain"]
contam_scaff_taxonomy[is.na(assigned_tax), "assignment" := "Unknown"]
contam_scaff_taxonomy[is.na(assignment), "assignment" := "Incorrect"]

contam_scaff_taxonomy_N <- contam_scaff_taxonomy[, .("assignment_N" = .N), by = c("assignment", "phylum_true", "domain_true")]
contam_scaff_taxonomy_N[, "assignment_perc" := assignment_N / sum(assignment_N), by = "phylum_true"]

phylum_factor_levels <- contam_scaff_taxonomy_N[assignment == "Genus"][order(-assignment_perc)][, phylum_true]
all_phyla <- unique(contam_scaff_taxonomy_N[, phylum_true])
phylum_factor_levels <- c(phylum_factor_levels, all_phyla[all_phyla %chin% phylum_factor_levels == F])
contam_scaff_taxonomy_N[, "phylum_true" := factor(phylum_true, levels = phylum_factor_levels)]

setorder(contam_scaff_taxonomy_N, phylum_true)
x_axis_colors <- domain_colors$Domain[as.character(unique(contam_scaff_taxonomy_N, by = "phylum_true")[, domain_true])]
contam_scaff_taxonomy_N[, "assignment" := factor(assignment, levels = c("Genus", "Family", "Order", "Class", "Phylum", "Domain", "Incorrect", "Unknown"))]

bar_colors <- c("#4FA3BCFF", colorRampPalette(wes_palette("Darjeeling1", 5, type = "discrete"))(6)[c(5:3, 1, 2)])
bar_colors <- c(alpha(bar_colors[1], 0.9), alpha(bar_colors[2], 0.7), alpha(bar_colors[3], 0.7), alpha(bar_colors[4], 0.7), alpha(bar_colors[5], 0.7), alpha(bar_colors[6], 0.7))

#Plotting contaminant scaffold taxonomy calls per phylum (Figure 2d)
fig_2d <- ggplot(data = contam_scaff_taxonomy_N, aes(x = phylum_true, y = assignment_N, fill = assignment)) +
 geom_bar(stat = "identity", position = "fill", width = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = bar_colors, name = "Assigned taxon level") +
 theme(legend.position = "right", legend.margin = margin(0, 0, 0, 0)) +
 theme(legend.title = element_text(size = 5.5))  +
 theme(legend.text = element_text(size = 5.2))  +
 theme(legend.key.size = unit(0.28, 'cm')) +
 theme(panel.grid = element_blank()) +  
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Phylum") +
 ylab("% of sequences") +
 scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 20, 40, 60, 80, 100)) +
 theme(axis.title = element_text(size = 6.5)) +
 theme(axis.text.x = element_text(size = 5.2, angle = 40, hjust = 0.9, color = x_axis_colors)) +
 theme(axis.text.y = element_text(size = 6)) 

ggsave(fig_2d, file = "fig_2d.svg", height = 1.8, width = 4.53)