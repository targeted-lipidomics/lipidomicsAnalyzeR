# Function to combine data matrices
generateDataMatrix = function(fdata,
                              metadata,
                              lcms_fns,
                              datatype = "Area",
                              tl_headers = c("ID","Name","Area", "Conc.", "Response"),
                              use_type ="Ratio",
                              snr = 5,
                              blank_filter = 5,
                              replaceNA = F,
                              batch_correction = F,
                              qc_label = "Sample",
                              factor_name = "Sample_type"){

  # Select feature names to extract
  fnames <- fdata %>% filter(Include == "YES" & Use_type == use_type) %>% select("Compound")

  # Empty final matrix to add to
  out_matrix = data.frame(matrix(nrow = 0, ncol = length(fnames[,1])))
  colnames(out_matrix) = fnames[,1]

  for(fn_ind in 1:length(lcms_fns)){

    # Select fn from vectors
    lcms_fn = lcms_fns[fn_ind]

    #Load lcms data
    lcms <- read.csv(lcms_fn, header = T) %>%
      select(tl_headers)

    # Filter metadata by batch and Include column
    metadata_batch <- metadata %>%
      subset(Chrom_Batch == fn_ind) %>%
      filter(Include == "YES")


    # Generate table from TargetLynx output and trim
    lcms_table <- extractTable(lcms) %>%
      subset(Name != "")

    out_table =  lcms_table %>%
      select("ID","Name", datatype) %>%         # Select the column you want to pivot longer
      group_by(ID) %>%                        # These two lines to avoid repeated values
      pivot_wider(names_from = "ID", values_from = datatype) %>%
      select(Name, fnames[,1]) %>%
      tibble::column_to_rownames("Name") %>%
      mutate_at(fnames[,1], as.numeric)



    if(snr != FALSE){
      snr_table = lcms_table %>%
        select("ID","Name", "S.N") %>%         # Select the column you want to pivot longer
        group_by(ID) %>%                        # These two lines to avoid repeated values
        pivot_wider(names_from = "ID", values_from = "S.N") %>%
        select(Name, fnames[,1]) %>%
        tibble::column_to_rownames("Name") %>%
        mutate_at(fnames[,1], as.numeric)

      out_table[snr_table < snr] = NA

    }

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

    # Remove everything but real samples
    sample_names = metadata_batch$Name[which(metadata$Sample_type == "Sample")]

    out_table = subset(out_table, row.names(out_table) %in% sample_names)


    # Calculate the concentration of lipid in each sample - TO BE COMPLETED
    if(datatype == "Conc."){

      return("FAIL - not implemented code to calculate concentration in biological samples")
      # conc based on cal curve is ng/mL
      # Step 1 = Multiply by volume of cal curve injected (mL)
          # Gives the amount of lipid per sample (ng) ???

      # Step 2 = Divide by the change in volume (extracted to analysed): Amount of lipid / vol of plasma
          # Gives the conc of lipid in each sample prior to sample prep (ng / mL)
    }


    ###################################################################
    ###### Replace missing values
    # NA to 0.2 * Lowest detected value)

    # First calculate the minimum values that will replace NA (20% of the minimum value reported, including all samples)

    if(replaceNA == T){
      min_vals = apply(out_table, 2, FUN = min, na.rm=T)

      for(col_idx in 1:length(min_vals)){

        min_val = min_vals[col_idx]

        if(!is.infinite(min_val)){

          temp_col = out_table[,col_idx]
          temp_col = ifelse(is.na(temp_col), (min_val * 0.2), temp_col)

          out_table[,col_idx] = temp_col
        }
      }
    }

    # Concatenate files
    out_matrix <- rbind(out_matrix,
                        out_table)

  }

  fdata_output = fdata %>%
    select(Compound, PUFAS, Pathways, Enzymes, Precursors, Double_bond,
           RT, Polarity, SRM1, SRM2)


  # Create dataset experiment
  lcms_experiment = DatasetExperiment(data=out_matrix,
                                      sample_meta=metadata %>% subset(Sample_type == "Sample"),
                                      variable_meta=fdata_output)

  if(batch_correction == T){

    bc_wf =
      batch_correct(qc_label = qc_label,
                    factor_name = factor_name)

    lcms_experiment = model_apply(bc_wf, lcms_experiment)
  }


  return(lcms_experiment)
}
