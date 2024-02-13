setClassUnion("character_NULL", c('character', 'NULL'))

parameterClass = setClass('parameterClass',
                          representation("VIRTUAL"))


pathsClass = setClass('pathsClass',
                      slots = c(path_project = 'character_NULL',
                                path_alternative = 'character_NULL',
                                path_data = 'character_NULL',
                                path_mzml = 'character_NULL',
                                path_raw = 'character_NULL',
                                path_workflows = 'character_NULL',
                                path_envs = 'character_NULL',
                                path_results = 'character_NULL',
                                path_metadata = 'character_NULL',
                                path_annotations = 'character_NULL',
                                path_reports = 'character_NULL',
                                path_rdata = 'character_NULL',
                                path_xlsx = 'character_NULL',
                                path_scripts = 'character_NULL',
                                path_results_TPs = 'character_NULL',
                                path_mspurity = 'character_NULL',
                                path_resources = 'character_NULL',
                                path_compounds_df = 'character_NULL',
                                path_metfrag ='character_NULL',
                                nist_database = 'character_NULL'
                      ))

setGeneric("validate_paths", function(object, ...) standardGeneric("validate_paths"))

#' Method to validate pathsClass
#' @export
#' @template validate_paths
setMethod(f="validate_paths", signature="pathsClass",  definition= function(object){
    
    nms <- slotNames(object)
    
    paths = sapply(nms, FUN=function(nm){
        path = slot(object, nm)
        
        exists=T
        
        if(is.null(path)){
            print(sprintf("No path provided for %s, passing with caution.", path))
        } else if(!dir.exists(path) & !file.exists(path)){
            print(sprintf("Path for %s does not exist, CHECK %s path.", path, path))
            exists=F}
        return(exists)
    })
    
    if(any(paths ==F)){
        stop("Check paths!!!!")
    }
    
    return(object)
})

mspurity_parsClass = setClass('mspurity_parsClass',
                                       slots = c(
                                           sample_group = 'character',
                                           sample_level = 'character'
                                       ))

sygma_parsClass = setClass('sygma_parsClass',
                            slots = c(
                                phase1_cycles = 'numeric',
                                phase2_cycles = 'numeric',
                                score_threshold = 'numeric',
                                path_substances = 'character_NULL',
                                ppm = 'numeric',
                                rerun = 'logical',
                                isomers_together = 'logical',
                                filter_small_fragments = 'logical',
                                substance_list = 'vector',
                                col_comp = 'character',
                                beamspy_isotopes_filename = 'character',
                                path_peaklist =  'character',
                                beamspy_adducts_filename =  'character'
                            ))

setGeneric("validate_sygma", function(object, ...) standardGeneric("validate_sygma"))

#' Method to validate sygma_parsClass
#' @export
#' @template validate_sygma
setMethod(f="validate_sygma", signature="sygma_parsClass",  definition= function(object){
    
    if(object@ppm > 10){
        print(sprintf("ppm_tol of %s exceeds usual values. CHeck this is correct", object@ppm))
    }
    if(object@phase1_cycles > 10){
        print(sprintf("Number of phase1 cycles of %s exceeds usual values. CHeck this is correct", object@phase1_cycles))
    }
    
    return(object)
})


filterFragSpectra_parsClass = setClass('filterFragSpectra_parsClass',
                           slots = c(
                               ilim = 'numeric',
                               plim = 'numeric',
                               ra = 'numeric',
                               snr = 'numeric',
                               rmp = 'logical',
                               snmeth = 'character'
                           ))

averageIntraFragSpectra_parsClass = setClass('averageIntraFragSpectra_parsClass',
                                       slots = c(
                                           minfrac = 'numeric',
                                           minnum = 'numeric',
                                           ppm = 'numeric',
                                           snr = 'numeric',
                                           ra = 'numeric',
                                           av = 'character',
                                           sumi = 'logical',
                                           rmp = 'logical'
                                       ))

averageInterFragSpectra_parsClass = setClass('averageInterFragSpectra_parsClass',
                                             slots = c(
                                                 minfrac = 'numeric',
                                                 minnum = 'numeric',
                                                 ppm = 'numeric',
                                                 snr = 'numeric',
                                                 ra = 'numeric',
                                                 av = 'character',
                                                 sumi = 'logical',
                                                 rmp = 'logical'
                                             ))



chroms_parsClass = setClass('chroms_parsClass',
                            slots = c(
                                col_class = 'character',
                                sample_group = 'character',
                                sample_level = 'vector',
                                peaklist = 'character',
                                ppm_tol = 'numeric',
                                rt_tol = 'numeric',
                                chroms_suffix = 'character'
                            ))

setGeneric("validate_chroms", function(object, ...) standardGeneric("validate_chroms"))

#' Method to validate chroms_parsClass
#' @export
#' @template validate_chroms
setMethod(f="validate_chroms", signature="chroms_parsClass",  definition= function(object){
    
    if(object@ppm_tol > 10){
        print(sprintf("ppm_tol of %s exceeds usual values. CHeck this is correct", object@ppm_tol))
    }
    if(object@rt_tol > 45){
        print(sprintf("rt_tol of %s exceeds usual values. CHeck this is correct", object@rt_tol))
    }
    
    return(object)
})



readMSData_parsClass = setClass('readMSData_parsClass',
                                slots = c(
                                    sample_col1 = 'character',
                                    sample_levels1 = 'vector',
                                    sample_col2 = 'character',
                                    sample_levels2 = 'vector',
                                    path_rdata = 'character',
                                    suffix = 'character'
                                ))


findChromPeaks_parsClass = setClass('findChromPeaks_parsClass',
                                    slots = c(
                                        ppm_tol = 'numeric',
                                        mzdiff = 'numeric',
                                        peakwidth = 'vector',
                                        snthresh = 'numeric',
                                        prefilter = 'vector',
                                        fitgauss = 'logical',
                                        noise = 'numeric',
                                        verboseColumns = 'logical',
                                        firstBaselineCheck = 'logical',
                                        polarity = 'character',
                                        path_xlsx = 'character',
                                        path_rdata = 'character'
                                    ))

targeted_findChromPeaks_parsClass = setClass('targeted_findChromPeaks_parsClass',
                                             slots = c(
                                                 snthreshIsoROIs = 'numeric',
                                                 maxCharge = 'numeric',
                                                 maxIso = 'numeric',
                                                 mzIntervalExtension = 'logical',
                                                 ppm_tol_targeted = 'numeric',
                                                 rt_tol_targeted = 'numeric',
                                                 path_lib = 'character_NULL'
                                             ))


adjustRtime_parsClass = setClass('adjustRtime_parsClass',
                                 slots = c(
                                     centerSampleName = 'character_NULL',
                                     method = 'character',
                                     #binSize: 1.0
                                     binsize_obiwarp = 'numeric',
                                     distFun = 'character',
                                     response = 'numeric',
                                     gapInit = 'numeric',
                                     gapExtend = 'numeric',
                                     localAlignment = 'logical',
                                     wrtcor_suffix = 'character',
                                     rtcor_suffix = 'character',
                                     path_xlsx = 'character',
                                     path_rdata = 'character'
                                 ))


groupChromPeaks_parsClass = setClass('groupChromPeaks_parsClass',
                                     slots = c(
                                         wrtcor_suffix = 'character',
                                         rtcor_suffix = 'character',
                                         method = 'character',
                                         bw = 'numeric',
                                         minFraction = 'numeric',
                                         binsize_peakdensity = 'numeric',
                                         maxFeatures = 'numeric',
                                         group_minFrac = 'character_NULL',
                                         path_xlsx = 'character',
                                         path_rdata = 'character'
                                     ))


mzcor_parsClass = setClass('mzcor_parsClass',
                           slots = c(
                               correction_apply = 'vector',
                               min_frac = 'numeric',
                               ppm_toll_search = 'numeric',
                               rt_toll_search = 'numeric',
                               rtcor = 'logical',
                               injection_cor = 'numeric',
                               mzcor_level = 'character_NULL',
                               error_tol = 'character_NULL',
                               apply_mzcor = 'logical',
                               LS_mzcor = 'logical',
                               MTox_mzcor = 'logical',
                               CD_mzcor = 'logical'
                           ))

setGeneric("validate_mzcor", function(object, ...) standardGeneric("validate_mzcor"))

#' Method to validate mzcor_parsClass
#' @export
#' @template validate_mzcor
setMethod(f="validate_mzcor", signature="mzcor_parsClass",  definition= function(object){
    
    if(is.null(object@mzcor_level)){
        print("No value provided for mzcor_level. Add later to correct annotations.")
    } else if(!object@mzcor_level %in% c("sample", "sample_mz", "mzOnly")){
        stop("Unexpected value for mzcor_level!!")
    }
    
    if(is.null(object@error_tol)){
        print("No value provided for error_tol. Add later to correct annotations.")
    } else if(!object@error_tol %in% c("", "_10ppm", "_5ppm")){
        print("Unusual value for error_tol, check the value. Passing with caution.")
    }
    
    return(object)
})

xcms_report_parsClass = setClass('xcms_report_parsClass',
                                 slots = c(
                                     rtcor_suffix = 'character',
                                     rtcor = 'logical',
                                     mzcor = 'logical',
                                     rt_lib = 'character',
                                     rt_prefix = 'character',
                                     rt_toll_search = 'numeric',
                                     rt_toll = 'numeric',
                                     ppm_toll_search = 'numeric',
                                     ppm_factor = 'numeric',
                                     path_reports = 'character_NULL'
                                 ))

#######################################################################################
#######################################################################################

setClassUnion("tool_classes", c('NULL', 'mspurity_parsClass', 'sygma_parsClass', 'filterFragSpectra_parsClass', 'averageIntraFragSpectra_parsClass', 'averageInterFragSpectra_parsClass', 'chroms_parsClass', 'readMSData_parsClass', 'findChromPeaks_parsClass',
                                'adjustRtime_parsClass', 'groupChromPeaks_parsClass', 'targeted_findChromPeaks_parsClass',
                                'xcms_report_parsClass', 'mzcor_parsClass'))

toolClass = setClass('toolClass',
                     slots = c(
                         name = 'character',
                         execute = 'logical',
                         index = 'numeric',
                         parameters = 'tool_classes'
                     ))

setClassUnion("tool_list", c('list', 'toolClass'))

workflowClass = setClass('workflowClass',
                         slots = c(category = 'character',
                                   version_script = 'numeric',
                                   execute = 'logical',
                                   tools = 'tool_list'
                         ))

comb_workflowsClass = setClass("comb_workflowsClass", contains="parameterClass",
                               slots=c(xeno = 'workflowClass',
                                       chroms = 'workflowClass',
                                       xcms = 'workflowClass',
                                       xcms_report = 'workflowClass'))

projectClass = setClass("projectClass", contains="parameterClass",
                        slots=c(project_name = 'character',
                                paths = 'pathsClass',
                                workflows = 'comb_workflowsClass'))


#######################################################################################
#######################################################################################




#' Function to load YAML and concatenate paths using anchors
#' @export
#' @param path_yaml Full path to YAML file
#' @return params list from the YAML
load_yaml = function(path_yaml, assay){
    params = read_yaml(path_yaml,
                       handlers = list(
                           concat_path = function(x){
                               path = paste(x, collapse="/")
                               print(path)
                               return(path)},
                           concat_str = function(x) {
                               path = paste(x, collapse="_")
                               print(path)
                               return(path)}
                       )
    )}


#' Function to extract specific sections from YAML
#' @export
#' @param yaml_params parameters list as read in from YAML file, at the level to subset
#' @param level the level in the YAML to subset i.e. "assays", "workflows", "tools"
#' @param name the name of the level to extract the parameters from i.e. "HILIC_POS", "xcms" "findChromPeaks"
#' @return params list subset to the desired level and name
extract_yaml_chunk = function(yaml_params, level, name){
    
    names = unlist(lapply(yaml_params[[level]], FUN=function(x){
        return(x$name)}))
    
    name_idx = which(names == name)
    
    if(length(name_idx) == 0){
        print(sprintf("%s not specified. Returning NULL", level))
        return(NULL)
    }
    
    return(yaml_params[[level]][[name_idx]])
}


#' Non-recursive version of function modifyList
#' @export
#' @param x A named list, possibly empty.
#' @param val A named list with components to replace corresponding components in x.
#' @return params list from the YAML
modify_list_nonrecursive <- function (x, val, keep.null = FALSE) {
    if(!is.list(x)){
        print("'x' not a list, returning 'val'")
        return(val)
    } else if(!is.list(val)){
        print("'val' not a list, returning 'x'")
        return(x)
    }
    
    xnames <- names(x)
    idx = which(sapply(xnames, FUN= function(xname){
        typeof(x[[xname]])}) == "list")
    if(as.numeric(length(idx)) > 0){
        xnames = xnames[-as.numeric(idx)]
    }
    
    vnames <- names(val)
    vnames <- vnames[nzchar(vnames)]
    idx = which(sapply(vnames, FUN= function(vname){
        typeof(val[[vname]])}) == "list")
    if(as.numeric(length(idx)) > 0){
        vnames = vnames[-as.numeric(idx)]
    }
    
    for (v in vnames) {
        x[[v]] <- val[[v]]
    }
    
    return(x)
}

#' Function to merge parameters from specific tools given by the default and config parameters.
#' Specifically default parameters are overwritten by the config parameters.
#' @export
#' @param default_params list read in from YAML containing default parameters for a specific workflow and tool.
#' @param config_params list read in from YAML containing assay specific parameters for a specific workflow and tool.
#' @param wf string defining the name of the workflow.
#' @param tool string defining the name of the tool.
#' @return list of updated workflow and tool parameters.
combine_wf_tools = function(default_params, config_params, wf, tool){
    
    # Workflow
    default_wf = extract_yaml_chunk(yaml_params = default_params$project,
                                    level = "workflows",
                                    name = wf)
    
    assay_wf = extract_yaml_chunk(yaml_params = config_params,
                                  level = "workflows",
                                  name = wf)
    
    wf = modify_list_nonrecursive(
        default_wf,
        assay_wf)
    
    # Tool
    default_tool = extract_yaml_chunk(yaml_params = default_wf,
                                      level = "tools",
                                      name = tool)
    
    assay_tool = extract_yaml_chunk(yaml_params = assay_wf,
                                    level = "tools",
                                    name = tool)
    
    tool = modify_list_nonrecursive(
        default_tool,
        assay_tool)
    
    
    # Chroms tool pars
    
    tool_params = modify_list_nonrecursive(
        x = default_tool$parameters,
        val = assay_tool$parameters)
    
    return(list(wf = wf,
                tool = tool,
                tool_params = tool_params))
}