#' @eval get_description('batch_correct')
#' @export batch_correct
#' @examples
#' D = iris_DatasetExperiment()
#' M = batch_correct(factor_name='Species',qc_label='all')
#' M = model_apply(M,D)
#'
batch_correct = function(
    qc_label='QC',
    factor_name,
    qc_frac=0,
    sample_frac=0,
    ...) {

  out=struct::new_struct('batch_correct',
                         qc_label=qc_label,
                         factor_name=factor_name,
                         qc_frac=qc_frac,
                         sample_frac=sample_frac,
                         ...)
  return(out)
}


.batch_correct<-setClass(
  "batch_correct",
  contains = c('model'),
  slots=c(
    qc_label='entity',
    factor_name='entity',
    qc_frac='entity',
    sample_frac='entity',
    corrected='entity',
    coeffs='entity'
  ),
  prototype=list(name = 'Batch correction',
                 description = paste0(
                   'Batch correction ',
                   'to correct for differences in performance in different analytical batches. '),
                 type = 'Batch correction',
                 predicted = 'corrected',
                 .params=c('qc_label','factor_name','qc_frac','sample_frac'),
                 .outputs=c('corrected','coeffs'),

                 qc_label=entity(name = 'QC label',
                                 description = 'The label used to identify QC samples.',
                                 value = 'QC',
                                 type='character',
                                 max_length = 1
                 ),

                 qc_frac=entity(name = 'QC fraction',
                                description=paste0(
                                  "A value between 0 and 1 to indicate the minimum proportion ",
                                  "of QC samples a feature must be present in for it to be ",
                                  "included when computing the reference. Default qc_frac = 0. "
                                ),
                                type='numeric',
                                value=0,
                                max_length = 1
                 ),

                 sample_frac=entity(name = 'Sample fraction',
                                    description=paste0(
                                      "A value between 0 and 1 to indicate the minimum proportion ",
                                      "of samples a feature must be present in for it to be ",
                                      "considered when computing the normalisation coefficients. "
                                    ),
                                    type='numeric',
                                    value=0,
                                    max_length = 1
                 ),
                 corrected=entity(name = 'Batch corrected DatasetExperiment',
                                   description = 'A DatasetExperiment object containing the batch corrected data.',
                                   type='DatasetExperiment',
                                   value=DatasetExperiment()
                 ),
                 coeffs=entity(name = 'Batch coefficients',
                              description = 'The batch specific coefficients per feature',
                              type='data.frame',
                              value=data.frame()
                 ),

                 factor_name=entity(name = 'Factor name',
                                    description = 'The name of the factor with QC label.',
                                    value = 'Sample_type',
                                    type='character',
                                    max_length = 1
                 )
  )
)

#' @export
#' @template model_train
setMethod(f="model_train",
          signature=c("batch_correct","DatasetExperiment"),
          definition=function(M,D)
          {
            opt=param_list(M)

            smeta=D$sample_meta
            x=D$data
            fdata = D$variable_meta

            coeffs = data.frame(matrix(nrow=0, ncol=(1+nrow(fdata))))
            colnames(coeffs) = c("Batch", fdata$Compound)

            # Reference samples (QCs or all samples)
            ref_samples = smeta$Name[which(smeta[,opt$factor_name] == opt$qc_label)]
            ref_df = x[which(rownames(x) %in% ref_samples),]

            # Grand median intenisty (of median QC or all samples) for each feature
            grand_med = as.numeric(apply(ref_df, 2, FUN = median, na.rm=T))


            for(i in unique(smeta$Chrom_Batch)){

              # Batch specific reference samples (QCs or all samples in specific batch)
              sample_names = smeta$Name[which(smeta$Chrom_Batch == i)]
              sample_names[which(sample_names %in% ref_samples)]
              temp_df = x[which(rownames(x) %in% sample_names),]

              # Batch median intensity for each feature
              batch_med = as.numeric(apply(temp_df, 2, FUN = median, na.rm=T))

              # Coefficient used to normale (grand / batch intenisty)
              coeff = batch_med / grand_med

              # Correct each feature by batch specific coeff
              x[which(rownames(x) %in% sample_names),] = t( apply(temp_df, 1, function(x) x / coeff) )

              # Store coefficients
              coeffs[i,] = c(i, coeff)

            }

            D$data = x

            output_value(M,'corrected') = D
            output_value(M,'coeffs') = coeffs

            #A=attributes(corrected)$processing_history$batch_correct

            return(M)
          }
)


#' @export
#' @template model_predict
setMethod(f="model_predict",
          signature=c("batch_correct","DatasetExperiment"),
          definition=function(M,D)
          {
            opt=param_list(M)

            smeta=D$sample_meta
            x=D$data
            fdata = D$variable_meta

            coeffs = data.frame(matrix(nrow=0, ncol=(1+nrow(fdata))))
            colnames(coeffs) = c("Batch", fdata$Compound)

            # Reference samples (QCs or all samples)
            ref_samples = smeta$Name[which(smeta[,opt$factor_name] == opt$qc_label)]
            ref_df = x[which(rownames(x) %in% ref_samples),]

            # Grand median intenisty (of median QC or all samples) for each feature
            grand_med = as.numeric(apply(ref_df, 2, FUN = median, na.rm=T))


            for(i in unique(smeta$Chrom_Batch)){

              # Batch specific reference samples (QCs or all samples in specific batch)
              sample_names = smeta$Name[which(smeta$Chrom_Batch == i)]
              sample_names[which(sample_names %in% ref_samples)]
              temp_df = x[which(rownames(x) %in% sample_names),]

              # Batch median intensity for each feature
              batch_med = as.numeric(apply(temp_df, 2, FUN = median, na.rm=T))

              # Coefficient used to normale (grand / batch intenisty)
              coeff = batch_med / grand_med

              # Correct each feature by batch specific coeff
              x[which(rownames(x) %in% sample_names),] = t( apply(temp_df, 1, function(x) x / coeff) )

              # Store coefficients
              coeffs[i,] = c(i, coeff)

            }

            D$data = x

            output_value(M,'corrected') = D
            output_value(M,'coeffs') = coeffs

            #A=attributes(corrected)$processing_history$batch_correct

            return(M)
          }
)
