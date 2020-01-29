# ThompsonLab

This repository contains the various scripts used in the Thompson lab. Currently, user-ready version include:

  1) **SirMixaPlot.R** - Current version 8.0 (Danger! High voltage) developed for the YggData ImageJ extraction macro
  2) **YggData.ijm** - Current version 1.0 developed for nuclear and whole cell iamge extraction for SirMixaPloyv8
  3) **jetData.ijm** - Fiji ImageJ macro for extracted ROI information across any number of channel pictures for SirMixaPlotv7
  4) **joinR** - contains an altered jetData.ijm that will generally extract ROI data without a pre-determined ROI. This folder also contains the R script joinR, which will use a 'nearest neighbor' approach to assign each ROI to a specific cell number (**No longer supported**)
  5) **RunnR** - Current version 2.0 (GraphParty). RunnR is an R package for graphing and analyzing line analyses
  
For version specific changes/program use, consult the annotation within each script
