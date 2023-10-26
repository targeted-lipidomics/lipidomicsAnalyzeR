# Stats workflow

pmp_wf =
  filter_smeta(mode = "include",
               levels = c("Lung"), #, "Trachea"),
               factor_name = "Sample_type1") +
  mv_feature_filter(threshold = 49,
                    method = "across",
                    factor_name = "Sample_type1") +
  mv_sample_filter(mv_threshold = 49) +
  structToolbox::pqn_norm(qc_label = "Sample",
                          factor_name = "Sample_type",
                          qc_frac = 0,
                          sample_frac = 0,
                          ref_method = "mean") +
  #structToolbox::vec_norm() +
  structToolbox::knn_impute(neighbours=5,
                            sample_max=50,
                            feature_max=50) +
  #structToolbox::pareto_scale() +
  structToolbox::log_transform() +
  structToolbox::mean_centre()+
  PCA(number_components=5)
