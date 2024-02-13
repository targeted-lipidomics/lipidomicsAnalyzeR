# Function to combine data matrices

#- TargetLynx data were combined and processed (SNR, missing values, batch correction).

#In the case of selecting data type concentration:
#  - The concentration predicted from TargetLynx using the calibration curve (ng/mL or Conc.) is normalised such that:

#  conc = conc * (std_vol / sample_vol) * (sample_IS / std_IS)


generateDataMatrix = function(fdata,
                              metadata,
                              lcms_fns,
                              datatype = "Area",
                              tl_headers = c("ID", "Name","Area", "ng/mL", "Response", "S/N"),
                              use_type ="Ratio",
                              snr = 5,
                              blank_filter = 5,
                              replaceNA = F,
                              batch_correction = F,
                              bc_qc_label = "Sample",
                              bc_factor_name = "Sample_type",
                              bc_header = "extract_batch",
                              keep_samples = c("Blank", "Sample", "QC"),
                              keep_sample_head = "Sample_type"){

  # Select feature names to extract
  fnames <- fdata %>% filter(Report_peak == "YES") %>% select("Compound")

  #Load lcms data
  lcms <- do.call(rbind, lapply(lcms_fns, read.csv, header=T)) %>%
    janitor::row_to_names(row_number = 6, remove_row = F, remove_rows_above = F) %>%
    rename(ID = 1) %>%
    select(tl_headers)

  # Generate table from TargetLynx output and trim
  lcms_table <- extractTable(lcms) %>%
    subset(Name != "")

  out_table =  lcms_table %>%
    distinct() %>%
    select("ID","Name", datatype) %>%         # Select the column you want to pivot longer
    group_by(ID) %>%                        # These two lines to avoid repeated values
    pivot_wider(names_from = "ID", values_from = datatype) %>%
    select(Name, fnames[,1]) %>%
    tibble::column_to_rownames("Name") %>%
    mutate_at(fnames[,1], as.numeric)

  if(snr != FALSE){
    snr_table = lcms_table %>%
      distinct() %>%
      select("ID","Name", "S/N") %>%         # Select the column you want to pivot longer
      group_by(ID) %>%                        # These two lines to avoid repeated values
      pivot_wider(names_from = "ID", values_from = "S/N") %>%
      select(Name, fnames[,1]) %>%
      tibble::column_to_rownames("Name") %>%
      mutate_at(fnames[,1], as.numeric)

    out_table[snr_table < snr] = NA

  }

  # Empty final matrix to add to
  out_matrix = data.frame(matrix(nrow = 0, ncol = length(fnames[,1])))
  colnames(out_matrix) = fnames[,1]


  for(batch in unique(metadata[[bc_header]])){

    # Filter metadata by batch and Include column
    metadata_batch <- metadata %>%
      subset(.[[bc_header]] == batch) %>%
      filter(Include == "YES")

    if(blank_filter != FALSE){

      blank_samples = metadata_batch$Name[which(metadata_batch$Sample_type1 == "Blank")]
      out_table_blank = out_table[rownames(out_table) %in% blank_samples, ]


      ext_blank_samples = metadata_batch$Name[which(metadata_batch$Sample_type1 == "ExtractBlank")]
      out_table_ext_blank = out_table[rownames(out_table) %in% ext_blank_samples, ]


      blank_df = data.frame(feature = names(out_table_blank),
                            blank = colMeans(out_table_blank, na.rm=T),
                            ext_blank = colMeans(out_table_ext_blank, na.rm=T)) %>%
        mutate(max_val = pmax(blank, ext_blank, na.rm=T))

      blank_df$blank_thresh = blank_df$max_val * blank_filter

      for(col_idx in which(!is.na(blank_df$blank_thresh))){

        thresh = blank_df$blank_thresh[col_idx]

        temp_col = out_table[,col_idx]
        temp_col = ifelse(temp_col < thresh, NA, temp_col)

        out_table[,col_idx] = temp_col
      }
    }

    # Keep only samples of interest
    sample_names = metadata_batch$Name[which(metadata_batch[[keep_sample_head]] %in% keep_samples)]

    out_table_batch = subset(out_table, row.names(out_table) %in% sample_names)


    # Calculate the concentration of lipid in each sample - TO BE COMPLETED
        # conc from cal curve (ng/mL) - based on known conc of cal mixes and adjustment factor for each std group
    # To find conc in each sample
      # conc = conc * (std_vol / sample_vol) * (sample_IS / std_IS)

    if(datatype %in% c("Conc.", "ng/mL")){

      sample_IS_vol_uL = unique(metadata_batch$sample_IS_vol_uL)
      cal_IS_vol_uL = unique(metadata_batch$cal_IS_vol_uL)
      sample_volume_uL = unique(metadata_batch$sample_volume_uL)
      cal_vol_uL = unique(metadata_batch$cal_vol_uL)

      if( all(c(length(sample_IS_vol_uL), length(cal_IS_vol_uL), length(sample_volume_uL), length(cal_vol_uL)) == 1)){

        IS_fact = sample_IS_vol_uL / cal_IS_vol_uL
        vol_fact = cal_vol_uL / sample_volume_uL

        out_table_batch = out_table_batch * vol_fact * IS_fact

      } else{
        return("FAIL - variable volumes, implement sample wise (not written yet)")
      }
    }


    ###################################################################
    ###### Replace missing values
    # NA to 0.2 * Lowest detected value)

    # First calculate the minimum values that will replace NA (20% of the minimum value reported, including all samples)

    if(replaceNA == T){
      min_vals = apply(out_table_batch, 2, FUN = min, na.rm=T)

      for(col_idx in 1:length(min_vals)){

        min_val = min_vals[col_idx]

        if(!is.infinite(min_val)){

          temp_col = out_table_batch[,col_idx]
          temp_col = ifelse(is.na(temp_col), (min_val * 0.2), temp_col)

          out_table_batch[,col_idx] = temp_col
        }
      }
    }

    # Concatenate files
    out_matrix <- rbind(out_matrix,
                        out_table_batch)

  }

  fdata_output = fdata %>%
    subset(Compound %in% colnames(out_matrix)) #%>%
    #select(Export, Compound, PUFAS, Pathways, Enzymes, Precursors, Double_bond,
    #       RT, Polarity, SRM1, SRM2)


  # Create dataset experiment
  lcms_experiment = DatasetExperiment(data=out_matrix,
                                      sample_meta=metadata[which(metadata[[keep_sample_head]] %in% keep_samples), ],
                                      variable_meta=fdata_output)

  if(batch_correction == T){

    bc_wf =
      batch_correct(qc_label = bc_qc_label,
                    factor_name = bc_factor_name)

    lcms_experiment = model_apply(bc_wf, lcms_experiment)
  }


  return(lcms_experiment)
}
