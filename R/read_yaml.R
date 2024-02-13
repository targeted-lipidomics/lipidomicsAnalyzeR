.libPaths(c("E:/R-4.0.3/library")) # local

library(yaml)
library(optparse)
library(dplyr)
library(BiocParallel)

#######################################################
#######################################################
########### Read in YAML ##############################
#######################################################
#######################################################


# Get arguments (single argument to *yml)
option_list <- list(
    make_option(c("-y", "--path_yaml"),
                type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

path_yaml = opt$path_yaml

print(sprintf("%s/general_functions/params_class_functions.R", path_scripts))
source(sprintf("%s/general_functions/params_class_functions.R", path_scripts))

default_yaml = "workflows_config_defaults.yml"
project_yaml = "workflows_project_config.yml"


# Load default paramaters
default_params = load_yaml(path_yaml = sprintf("%s/%s", path_yaml, default_yaml),
                           assay = assay)

# Load project parameters
project_params = load_yaml(path_yaml = sprintf("%s/%s", path_yaml, project_yaml))


# Extract project parameters for specific assay
assay_params = extract_yaml_chunk(yaml_params = project_params$project,
                                  level = "assays",
                                  name = assay)

#######################################################
#######################################################
########### Project and path pars #####################
#######################################################
#######################################################

#############
#### Project
# Update non-recursive
project = modify_list_nonrecursive(
    default_params$project, project_params$project)


#############
#### Project
# Update non-recursive
paths = modify_list_nonrecursive(
    default_params$project$paths, project_params$project$paths)

paths = append(paths,
               list(path_mzml = assay_params[["paths"]][["path_mzml"]],
                    path_raw = assay_params[["paths"]][["path_raw"]]))

project_name = project[["name"]]

###################
#### Into S4 class

paths_obj = pathsClass(path_project = paths[["path_project"]],
                       #path_alternative = paths[["path_alternative"]],
                       path_data = paths[["path_data"]],
                       path_mzml = paths[["path_mzml"]],
                       path_raw = paths[["path_raw"]],
                       path_workflows = paths[["path_workflows"]],
                       path_envs = paths[["path_envs"]],
                       path_results = paths[["path_results"]],
                       path_metadata = paths[["path_metadata"]],
                       path_annotations = paths[["path_annotations"]],
                       path_reports = paths[["path_reports"]],
                       path_rdata = paths[["path_rdata"]],
                       path_xlsx = paths[["path_xlsx"]],
                       path_scripts = paths[["path_scripts"]],

                       path_results_TPs = paths[["path_results_TPs"]],
                       path_mspurity = paths[["path_mspurity"]],
                       path_resources = paths[["path_resources"]],
                       path_compounds_df = paths[["path_compounds_df"]],
                       path_metfrag = paths[["path_metfrag"]],
                       nist_database = paths[["nist_database"]]
)

#paths_obj <- validate_paths(paths_obj)

#######################################################
#######################################################
########### Workflow (inc tool) pars ##################
#######################################################
#######################################################

###################################################
###################################################
########### xeno_annotation workflow ##############
###################################################

# mspurity

xeno_annotation_wf_tools = combine_wf_tools(default_params = default_params,
                                            config_params = assay_params,
                                            wf = "xeno_annotation",
                                            tool = "mspurity")

xeno_annotation_wf = xeno_annotation_wf_tools[["wf"]]

mspurity_tool = xeno_annotation_wf_tools[["tool"]]
mspurity_tool_params = xeno_annotation_wf_tools[["tool_params"]]

mspurity_pars = mspurity_parsClass(
    sample_group = mspurity_tool_params[["sample_group"]],
    sample_level = mspurity_tool_params[["sample_level"]]
)

mspurity_tool = toolClass(
    name = mspurity_tool[["name"]],
    execute = mspurity_tool[["execute"]],
    index = mspurity_tool[["index"]],
    parameters = mspurity_pars )

# sygma

xeno_annotation_wf_tools = combine_wf_tools(default_params = default_params,
                                            config_params = assay_params,
                                            wf = "xeno_annotation",
                                            tool = "sygma")

xeno_annotation_wf = xeno_annotation_wf_tools[["wf"]]

sygma_tool = xeno_annotation_wf_tools[["tool"]]
sygma_tool_params = xeno_annotation_wf_tools[["tool_params"]]

sygma_pars = sygma_parsClass(
    phase1_cycles = sygma_tool_params[["phase1_cycles"]],
    phase2_cycles = sygma_tool_params[["phase2_cycles"]],
    score_threshold = sygma_tool_params[["score_threshold"]],
    path_substances = sygma_tool_params[["path_substances"]],
    ppm = sygma_tool_params[["ppm"]],
    rerun = sygma_tool_params[["rerun"]], #if False, avoid re-running annotation of peaklist isotopes/adducts
    isomers_together = sygma_tool_params[["isomers_together"]], #True = run metfrag on isomers at once
    filter_small_fragments = sygma_tool_params[["filter_small_fragments"]], #boolean, when True excludes metabolites with less than 15% of original atoms
    substance_list = sygma_tool_params[["substance_list"]],
    col_comp = sygma_tool_params[["col_comp"]],
    beamspy_isotopes_filename = sygma_tool_params[["beamspy_isotopes_filename"]],
    path_peaklist = sygma_tool_params[["path_peaklist"]],
    beamspy_adducts_filename = sygma_tool_params[["beamspy_adducts_filename"]]
)

sygma_pars = validate_sygma(sygma_pars)

sygma_tool = toolClass(
    name = sygma_tool[["name"]],
    execute = sygma_tool[["execute"]],
    index = sygma_tool[["index"]],
    parameters = sygma_pars )


# filterFragSpectra
xeno_annotation_wf_tools = combine_wf_tools(default_params = default_params,
                                            config_params = assay_params,
                                            wf = "xeno_annotation",
                                            tool = "filterFragSpectra")

xeno_annotation_wf = xeno_annotation_wf_tools[["wf"]]


filterFragSpectra_tool = xeno_annotation_wf_tools[["tool"]]
filterFragSpectra_tool_params = xeno_annotation_wf_tools[["tool_params"]]

filterFragSpectra_pars = filterFragSpectra_parsClass(
    ilim = filterFragSpectra_tool_params[["ilim"]],
    plim = filterFragSpectra_tool_params[["plim"]],
    ra = filterFragSpectra_tool_params[["ra"]],
    snr = filterFragSpectra_tool_params[["snr"]],
    rmp = filterFragSpectra_tool_params[["rmp"]],
    snmeth = filterFragSpectra_tool_params[["snmeth"]]
)

filterFragSpectra_tool = toolClass(
    name = filterFragSpectra_tool[["name"]],
    execute = filterFragSpectra_tool[["execute"]],
    index = filterFragSpectra_tool[["index"]],
    parameters = filterFragSpectra_pars )


# averageIntraFragSpectra
xeno_annotation_wf_tools = combine_wf_tools(default_params = default_params,
                                            config_params = assay_params,
                                            wf = "xeno_annotation",
                                            tool = "averageIntraFragSpectra")

xeno_annotation_wf = xeno_annotation_wf_tools[["wf"]]


averageIntraFragSpectra_tool = xeno_annotation_wf_tools[["tool"]]
averageIntraFragSpectra_tool_params = xeno_annotation_wf_tools[["tool_params"]]

averageIntraFragSpectra_pars = averageIntraFragSpectra_parsClass(
    minfrac = averageIntraFragSpectra_tool_params[["minfrac"]],
    minnum = averageIntraFragSpectra_tool_params[["minnum"]],
    ppm = averageIntraFragSpectra_tool_params[["ppm"]],
    snr = averageIntraFragSpectra_tool_params[["snr"]],
    ra = averageIntraFragSpectra_tool_params[["ra"]],
    av = averageIntraFragSpectra_tool_params[["av"]],
    sumi = averageIntraFragSpectra_tool_params[["sumi"]],
    rmp = averageIntraFragSpectra_tool_params[["rmp"]]
)

averageIntraFragSpectra_tool = toolClass(
    name = averageIntraFragSpectra_tool[["name"]],
    execute = averageIntraFragSpectra_tool[["execute"]],
    index = averageIntraFragSpectra_tool[["index"]],
    parameters = averageIntraFragSpectra_pars )

# averageInterFragSpectra
xeno_annotation_wf_tools = combine_wf_tools(default_params = default_params,
                                            config_params = assay_params,
                                            wf = "xeno_annotation",
                                            tool = "averageInterFragSpectra")

xeno_annotation_wf = xeno_annotation_wf_tools[["wf"]]

averageInterFragSpectra_tool = xeno_annotation_wf_tools[["tool"]]
averageInterFragSpectra_tool_params = xeno_annotation_wf_tools[["tool_params"]]

averageInterFragSpectra_pars = averageInterFragSpectra_parsClass(
    minfrac = averageInterFragSpectra_tool_params[["minfrac"]],
    minnum = averageInterFragSpectra_tool_params[["minnum"]],
    ppm = averageInterFragSpectra_tool_params[["ppm"]],
    snr = averageInterFragSpectra_tool_params[["snr"]],
    ra = averageInterFragSpectra_tool_params[["ra"]],
    av = averageInterFragSpectra_tool_params[["av"]],
    sumi = averageInterFragSpectra_tool_params[["sumi"]],
    rmp = averageInterFragSpectra_tool_params[["rmp"]]
)

averageInterFragSpectra_tool = toolClass(
    name = averageInterFragSpectra_tool[["name"]],
    execute = averageInterFragSpectra_tool[["execute"]],
    index = averageInterFragSpectra_tool[["index"]],
    parameters = averageInterFragSpectra_pars )


# Generate wf
xeno_annotation_wf = workflowClass(
    category = xeno_annotation_wf[["category"]],
    version_script = xeno_annotation_wf[["version_script"]],
    execute = xeno_annotation_wf[["execute"]],
    tools = list(
        mspurity = mspurity_tool,
        sygma = sygma_tool,
        filterFragSpectra = filterFragSpectra_tool,
        averageIntraFragSpectra = averageIntraFragSpectra_tool,
        averageInterFragSpectra = averageInterFragSpectra_tool)
)

#View(chroms_wf)


##########################################
##########################################
########### chroms workflow ##############
##########################################

chroms_wf_tools = combine_wf_tools(default_params = default_params,
                                   config_params = assay_params,
                                   wf = "chroms",
                                   tool = "chroms")

chroms_wf = chroms_wf_tools[["wf"]]
chroms_tool = chroms_wf_tools[["tool"]]
chroms_tool_params = chroms_wf_tools[["tool_params"]]

chroms_pars = chroms_parsClass(
    col_class = chroms_tool_params[["col_class"]],
    sample_group = chroms_tool_params[["sample_group"]],
    sample_level = chroms_tool_params[["sample_level"]],
    peaklist = chroms_tool_params[["peaklist"]],
    ppm_tol = chroms_tool_params[["ppm_tol"]],
    rt_tol = chroms_tool_params[["rt_tol"]],
    chroms_suffix = chroms_tool_params[["chroms_suffix"]] )

chroms_pars = validate_chroms(chroms_pars)

chroms_tool = toolClass(
    name = chroms_tool[["name"]],
    execute = chroms_tool[["execute"]],
    index = chroms_tool[["index"]],
    parameters = chroms_pars )

chroms_wf = workflowClass(
    category = chroms_wf[["category"]],
    version_script = chroms_wf[["version_script"]],
    execute = chroms_wf[["execute"]],
    tools = list(
        chroms = chroms_tool )
)

#View(chroms_wf)



##########################################
##########################################
########### xcms workflow ################
##########################################

# readMSData

xcms_wf_readMSData = combine_wf_tools(default_params = default_params,
                                      config_params = assay_params,
                                      wf = "xcms",
                                      tool = "readMSData")

xcms_wf = xcms_wf_readMSData[["wf"]]
readMSData_tool = xcms_wf_readMSData[["tool"]]
readMSData_params = xcms_wf_readMSData[["tool_params"]]

readMSData_pars = readMSData_parsClass(
    sample_col1 = readMSData_params[["sample_col1"]],
    sample_levels1 = readMSData_params[["sample_levels1"]],
    sample_col2 = readMSData_params[["sample_col2"]],
    sample_levels2 = readMSData_params[["sample_levels2"]],
    path_rdata = paths_obj@path_rdata,
    suffix = readMSData_params[["suffix"]] )

readMSData_tool = toolClass(
    name = readMSData_tool[["name"]],
    execute = readMSData_tool[["execute"]],
    index = readMSData_tool[["index"]],
    parameters = readMSData_pars )


# findChromPeaks

xcms_wf_findChromPeaks = combine_wf_tools(default_params = default_params,
                                          config_params = assay_params,
                                          wf = "xcms",
                                          tool = "findChromPeaks")

findChromPeaks_tool = xcms_wf_findChromPeaks[["tool"]]
findChromPeaks_params = xcms_wf_findChromPeaks[["tool_params"]]

findChromPeaks_pars = findChromPeaks_parsClass(
    ppm_tol = findChromPeaks_params[["ppm_tol"]],
    mzdiff = findChromPeaks_params[["mzdiff"]],
    peakwidth = findChromPeaks_params[["peakwidth"]],
    snthresh = findChromPeaks_params[["snthresh"]],
    prefilter = findChromPeaks_params[["prefilter"]],
    fitgauss = findChromPeaks_params[["fitgauss"]],
    noise = findChromPeaks_params[["noise"]],
    verboseColumns = findChromPeaks_params[["verboseColumns"]],
    firstBaselineCheck = findChromPeaks_params[["firstBaselineCheck"]],
    polarity = findChromPeaks_params[["polarity"]],
    path_xlsx = paths_obj@path_xlsx,
    path_rdata = paths_obj@path_rdata )

findChromPeaks_tool = toolClass(
    name = findChromPeaks_tool[["name"]],
    execute = findChromPeaks_tool[["execute"]],
    index = findChromPeaks_tool[["index"]],
    parameters = findChromPeaks_pars )

# targeted_findChromPeaks

xcms_wf_targeted_findChromPeaks = combine_wf_tools(default_params = default_params,
                                                   config_params = assay_params,
                                                   wf = "xcms",
                                                   tool = "targeted_findChromPeaks")

targeted_findChromPeaks_tool = xcms_wf_targeted_findChromPeaks[["tool"]]
targeted_findChromPeaks_params = xcms_wf_targeted_findChromPeaks[["tool_params"]]

targeted_findChromPeaks_pars = targeted_findChromPeaks_parsClass(
    snthreshIsoROIs = targeted_findChromPeaks_params[["snthreshIsoROIs"]],
    maxCharge = targeted_findChromPeaks_params[["maxCharge"]],
    maxIso = targeted_findChromPeaks_params[["maxIso"]],
    mzIntervalExtension = targeted_findChromPeaks_params[["mzIntervalExtension"]],
    ppm_tol_targeted = targeted_findChromPeaks_params[["ppm_tol_targeted"]],
    rt_tol_targeted = targeted_findChromPeaks_params[["rt_tol_targeted"]],
    path_lib = targeted_findChromPeaks_params[["path_lib"]] )

targeted_findChromPeaks_tool = toolClass(
    name = targeted_findChromPeaks_tool[["name"]],
    execute = targeted_findChromPeaks_tool[["execute"]],
    index = targeted_findChromPeaks_tool[["index"]],
    parameters = targeted_findChromPeaks_pars )

# adjustRtime

xcms_wf_adjustRtime = combine_wf_tools(default_params = default_params,
                                       config_params = assay_params,
                                       wf = "xcms",
                                       tool = "adjustRtime")

adjustRtime_tool = xcms_wf_adjustRtime[["tool"]]
adjustRtime_params = xcms_wf_adjustRtime[["tool_params"]]

adjustRtime_pars = adjustRtime_parsClass(
    centerSampleName = adjustRtime_params[["centerSampleName"]],
    method = adjustRtime_params[["method"]],
    binsize_obiwarp = adjustRtime_params[["binsize_obiwarp"]],
    distFun = adjustRtime_params[["distFun"]],
    response = adjustRtime_params[["response"]],
    gapInit = adjustRtime_params[["gapInit"]],
    gapExtend = adjustRtime_params[["gapExtend"]],
    localAlignment = adjustRtime_params[["localAlignment"]],
    wrtcor_suffix = adjustRtime_params[["wrtcor_suffix"]],
    rtcor_suffix = adjustRtime_params[["rtcor_suffix"]],
    path_xlsx = paths_obj@path_xlsx,
    path_rdata = paths_obj@path_rdata )

adjustRtime_tool = toolClass(
    name = adjustRtime_tool[["name"]],
    execute = adjustRtime_tool[["execute"]],
    index = adjustRtime_tool[["index"]],
    parameters = adjustRtime_pars )


# groupChromPeaks

xcms_wf_groupChromPeaks = combine_wf_tools(default_params = default_params,
                                           config_params = assay_params,
                                           wf = "xcms",
                                           tool = "groupChromPeaks")

groupChromPeaks_tool = xcms_wf_groupChromPeaks[["tool"]]
groupChromPeaks_params = xcms_wf_groupChromPeaks[["tool_params"]]

groupChromPeaks_pars = groupChromPeaks_parsClass(
    wrtcor_suffix = groupChromPeaks_params[["wrtcor_suffix"]],
    rtcor_suffix = groupChromPeaks_params[["rtcor_suffix"]],
    method = groupChromPeaks_params[["method"]],
    bw = groupChromPeaks_params[["bw"]],
    minFraction = groupChromPeaks_params[["minFraction"]],
    binsize_peakdensity = groupChromPeaks_params[["binsize_peakdensity"]],
    maxFeatures = groupChromPeaks_params[["maxFeatures"]],
    group_minFrac = groupChromPeaks_params[["group_minFrac"]],
    path_xlsx = paths_obj@path_xlsx,
    path_rdata = paths_obj@path_rdata )

groupChromPeaks_tool = toolClass(
    name = groupChromPeaks_tool[["name"]],
    execute = groupChromPeaks_tool[["execute"]],
    index = groupChromPeaks_tool[["index"]],
    parameters = groupChromPeaks_pars )



# mz_correct

xcms_wf_mzcor = combine_wf_tools(default_params = default_params,
                                 config_params = assay_params,
                                 wf = "xcms",
                                 tool = "mzcor")

mzcor_tool = xcms_wf_mzcor[["tool"]]
mzcor_params = xcms_wf_mzcor[["tool_params"]]

mzcor_pars = mzcor_parsClass(
    correction_apply = mzcor_params[["correction_apply"]],
    min_frac = mzcor_params[["min_frac"]],
    ppm_toll_search = mzcor_params[["ppm_toll_search"]],
    rt_toll_search = mzcor_params[["rt_toll_search"]],
    rtcor = mzcor_params[["rtcor"]],
    injection_cor = mzcor_params[["injection_cor"]],
    mzcor_level = mzcor_params[["mzcor_level"]],
    error_tol = mzcor_params[["error_tol"]],
    apply_mzcor = mzcor_params[["apply_mzcor"]],
    LS_mzcor = mzcor_params[["LS_mzcor"]],
    MTox_mzcor = mzcor_params[["MTox_mzcor"]],
    CD_mzcor = mzcor_params[["CD_mzcor"]],)

mzcor_pars = validate_mzcor(mzcor_pars)

mzcor_tool = toolClass(
    name = mzcor_tool[["name"]],
    execute = mzcor_tool[["execute"]],
    index = mzcor_tool[["index"]],
    parameters = mzcor_pars )


###############################
###############################


xcms_wf = workflowClass(
    category = xcms_wf[["category"]],
    version_script = xcms_wf[["version_script"]],
    execute = xcms_wf[["execute"]],
    tools = list(readMSData = readMSData_tool,
                 findChromPeaks = findChromPeaks_tool,
                 targeted_findChromPeaks = targeted_findChromPeaks_tool,
                 adjustRtime = adjustRtime_tool,
                 groupChromPeaks = groupChromPeaks_tool,
                 mzcor = mzcor_tool) )

#View(xcms_wf)

##########################################
##########################################
########### xcms report ##################
##########################################

xreport_wf_tools = combine_wf_tools(default_params = default_params,
                                    config_params = assay_params,
                                    wf = "xcms_report",
                                    tool = "xcms_report")

xreport_wf = xreport_wf_tools[["wf"]]
xreport_tool = xreport_wf_tools[["tool"]]
xreport_tool_params = xreport_wf_tools[["tool_params"]]

xreport_pars = xcms_report_parsClass(
    rtcor_suffix = xreport_tool_params[["rtcor_suffix"]],
    rtcor = xreport_tool_params[["rtcor"]],
    mzcor = xreport_tool_params[["mzcor"]],
    rt_lib = xreport_tool_params[["rt_lib"]],
    rt_prefix = xreport_tool_params[["rt_prefix"]],
    rt_toll_search = xreport_tool_params[["rt_toll_search"]],
    rt_toll = xreport_tool_params[["rt_toll"]],
    ppm_toll_search = xreport_tool_params[["ppm_toll_search"]],
    ppm_factor = xreport_tool_params[["ppm_factor"]],
    path_reports = xreport_tool_params[["path_reports"]] )

xreport_tool = toolClass(
    name = xreport_tool[["name"]],
    execute = xreport_tool[["execute"]],
    index = xreport_tool[["index"]],
    parameters = xreport_pars )

xreport_wf = workflowClass(
    category = xreport_wf[["category"]],
    version_script = xreport_wf[["version_script"]],
    execute = xreport_wf[["execute"]],
    tools = list(
        xcms_report = xreport_tool ))

#View(xreport_wf)


#######################################################
#######################################################
##### Combine workflows and paths to 1 S4 object ######
#######################################################
#######################################################

workflows = comb_workflowsClass(
    xeno = xeno_annotation_wf,
    chroms = chroms_wf,
    xcms = xcms_wf,
    xcms_report = xreport_wf
)

project_pars = projectClass(
    project_name = project_name,
    paths = paths_obj,
    workflows = workflows
)


#View(project_pars)

###############################
#### Write to YAML

save(project_pars, file = sprintf("%s/%s_parameters.RData", project_pars@paths@path_rdata, assay))

print(sprintf("Saved parameter object for %s.", assay))

sink(sprintf("%s/parameters/%s_parameters.txt", project_pars@paths@path_results, assay))

print(project_pars)
sink()
