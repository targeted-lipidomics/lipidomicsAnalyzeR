library(tidyverse)
library(magrittr)
library(tibble)
library(structToolbox)
library(openxlsx)

sapply(list.files("C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/processing_lipidomics/lipidomicsAnalyzeR/R", full.names = T), source)

# Set filenames

setwd("C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/Projects/GP_tissue_extract.PRO/analysis")

out_name_oct = "octadecanoids_peakArea.xlsx"
out_name_GOM = "GOM_peakArea.xlsx"

fdata_fn_oct <- "Octadecanoids_compoundInfo.csv"
metadata_fn_oct = "Octadecanoids_metadata.csv"
lcms_fns_oct <- c("Octadecanoids_lcms.csv") # In batch order

fdata_fn_GOM <- "GOM_compoundInfo.csv"
metadata_fn_GOM = "GOM_metadata.csv"
lcms_fns_GOM <- c("GOM_lcms.csv") # In batch order

#Load feature data info to include per Compound (Concentration or Ratio), generate a vector with compound names
fdata_oct <- read.csv(fdata_fn_oct, header=TRUE, encoding="UTF-8")
fdata_GOM <- read.csv(fdata_fn_GOM, header=TRUE, encoding="UTF-8")

# Load Metadata
metadata_oct <- read.csv(metadata_fn_oct, header=TRUE) %>% filter(Include == "YES")
metadata_GOM <- read.csv(metadata_fn_GOM, header=TRUE) %>% filter(Include == "YES")


###############################################################
# Octadecanoids - Generate datamatrices of selected features and response types
octa_area_NAs = generateDataMatrix(fdata = fdata_oct,
                                   metadata = metadata_oct,
                                   lcms_fns = lcms_fns_oct,
                                   datatype = "Area", # "Response"  "Conc."  "Area", "S/N"
                                   tl_headers = c("ID","Name","Area", "S.N"),
                                   use_type ="Ratio",
                                   snr = 5,
                                   blank_filter = 5,
                                   replaceNA = F,
                                   batch_correction = F,
                                   qc_label = "Sample",
                                   factor_name = "Sample_type")
# Save xlsx

wb_oct = createWorkbook()
# variable metadata
addWorksheet(wb_oct, "feature_metadata")
writeData(wb_oct, "feature_metadata", octa_area_NAs$variable_meta, colNames=TRUE, rowNames = FALSE, keepNA=FALSE)
# sample metadata
addWorksheet(wb_oct, "sample_metadata")
writeData(wb_oct, "sample_metadata", octa_area_NAs$sample_meta, colNames=TRUE, rowNames = FALSE, keepNA=FALSE)
# peak area matrix
addWorksheet(wb_oct, "peak_area_matrix")
writeData(wb_oct, "peak_area_matrix",
          octa_area_NAs$data %>%
            select(where(function(x) any(!is.na(x)))),
          colNames=TRUE, rowNames = TRUE, keepNA=FALSE)
# save workbook
saveWorkbook(wb_oct,
             file = out_name_oct,
             overwrite = TRUE)


# Save RData
saveRDS(octa_area_NAs, file = "octadecanoids_peakArea.RDS")


# Do some stats
apply_pmp_wf = model_apply(pmp_wf, octa_area_NAs)

scores_p = pca_scores_plot(factor_name="Sample_type2",
                           points_to_label = 'none')

chart_plot(scores_p,apply_pmp_wf[8])


loadings_df = apply_pmp_wf[8]@loadings %>%
  tibble::rownames_to_column("feature")

ggplot(loadings_df, aes(x=PC1, y=PC2)) +
  geom_point(size=2, alpha=1, shape=23) +
  theme_Publication() +
  stat_ellipse(type='norm', level=0.9) +
  ggrepel::geom_text_repel(size=2.5, label.padding = unit(0.2, "lines"),
                           col="red", aes(label=feature), nudge_x = 0, nudge_y = 0.01)



###############################################################
# GOM - Generate datamatrices of selected features and response types
GOM_area_NAs = generateDataMatrix(fdata = fdata_GOM,
                                  metadata = metadata_GOM,
                                  lcms_fns = lcms_fns_GOM,
                                  datatype = "Area", # "Response"  "Conc."  "Area", "S/N"
                                  tl_headers = c("ID","Name","Area", "S.N"),
                                  use_type ="Ratio",
                                  snr = 5,
                                  blank_filter = 5,
                                  replaceNA = F,
                                  batch_correction = F,
                                  qc_label = "Sample",
                                  factor_name = "Sample_type")
# Save xlsx

wb_GOM = createWorkbook()
# variable metadata
addWorksheet(wb_GOM, "feature_metadata")
writeData(wb_GOM, "feature_metadata", GOM_area_NAs$variable_meta, colNames=TRUE, rowNames = FALSE, keepNA=FALSE)
# sample metadata
addWorksheet(wb_GOM, "sample_metadata")
writeData(wb_GOM, "sample_metadata", GOM_area_NAs$sample_meta, colNames=TRUE, rowNames = FALSE, keepNA=FALSE)
# peak area matrix
addWorksheet(wb_GOM, "peak_area_matrix")
writeData(wb_GOM, "peak_area_matrix",
          GOM_area_NAs$data %>%
            select(where(function(x) any(!is.na(x)))),
          colNames=TRUE, rowNames = TRUE, keepNA=FALSE)
# save workbook
saveWorkbook(wb_GOM,
             file = out_name_GOM,
             overwrite = TRUE)

# Save RData
saveRDS(GOM_area_NAs, file = "GOM_peakArea.RDS")


# Do some stats
apply_pmp_wf = model_apply(pmp_wf, GOM_area_NAs)

scores_p = pca_scores_plot(factor_name="Sample_type2",
                    points_to_label = 'none')

chart_plot(scores_p,apply_pmp_wf[8])


loadings_df = apply_pmp_wf[8]@loadings %>%
  tibble::rownames_to_column("feature")

ggplot(loadings_df, aes(x=PC1, y=PC2)) +
  geom_point(size=2, alpha=1, shape=23) +
  theme_Publication() +
  stat_ellipse(type='norm', level=0.9) +
  ggrepel::geom_text_repel(size=2.5, label.padding = unit(0.2, "lines"),
                           col="red", aes(label=feature), nudge_x = 0, nudge_y = 0.01)
