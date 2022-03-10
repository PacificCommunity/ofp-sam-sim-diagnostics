# Diagnostic Part of Simulation Project

As part of the collaborative [Spatial Stock Assessment Simulation
Experiment](https://aaronmberger-nwfsc.github.io/Spatial-Assessment-Modeling-Workshop/)
([repo](https://github.com/aaronmberger-nwfsc/Spatial-Assessment-Modeling-Workshop)),
Nicholas sent an email (29 Jan 2022) to Arni and Thom, inviting them to
participate in running model diagnostics.

The email has three attachments:

* `launch.2021_diagnostic_retrospective_w_LLprof.r`
* `launch.2021_hindcast.r`
* `launch.2021_newExtIdx_dataUpdate_5kg_PeatmanExtComp_selexExplore.r`

The general purpose of this repository is to have a platform where we can
collaborate on developing and running model diagnostics.

It might be worthwhile to 'upgrade' this analysis to make it run on Arni's and
Thom's computers. With that goal in mind, the first steps could be:

1. Identify all dependencies to run the scripts, in terms of software and data
   files.

2. Make the data files available in a location that allows us to run the scripts
   on different computers without modifying the scripts.

Below are some notes related to objectives 1 and 2.

## Required packages beyond CRAN

* [frqit](https://github.com/PacificCommunity/ofp-sam-frqit) -
  `install_github("PacificCommunity/ofp-sam-frqit")`

## Required input files

**launch.2021_diagnostic_retrospective_w_LLprof.r** uses these directories:

Variable          | Path
----------------- | ---------------------------------------------------------------------------------------------
`dir.model_runs`  | `C:/Users/nicholasd/HOME/SPC/SPC_SAM/2021_SC/SWO/Assessment/Model_Runs/`
`dir.condor`      | `C:/Users/nicholasd/HOME/SPC/SPC_SAM/2021_SC/SWO/Assessment/condor_files/`
`dir.stem`        | `C:/Users/nicholasd/HOME/SPC/SPC_SAM/2021_SC/SWO/Assessment/Model_Runs/2021_diagnostic_case/`
`dir.frq`         | `C:/Users/nicholasd/HOME/SPC/SPC_SAM/2021_SC/SWO/Assessment/Data_Prep/frq_file/`
`dir.launch`      | `/home/nicholasd/`
`dir.mfcl.launch` | `/home/nicholasd/mfcl/`
