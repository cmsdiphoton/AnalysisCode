# AnalysisCode

** Analysis **
This directory contains the necessary code to make the analysis histograms from ggNtuples.

* Diempt.C will take the input gg and ff or ee histograms of diempt and calculate the diempt ratios which can be used to reweight the MET distributions. This macro also makes a nice plot of the diempt ratio and calculates the appropriate errors.
* FinalAnalysis.C takes the output of the event loop and calculates the final background estimates and errors and makes several plots of the final quantities.
* ratioMethod.C performs the fit to the gg/ff ratio.



** Limits **
Code to make limit plots from the datacards.

* PlotsSMSFiles is a directory that contains all the files you need to edit to use the PlotsSMS package to make nice looking limit plots.



** NtupleCode **
Code to manipulate the ggNtuples.

* Condor_skimNtuples contains example code to submit condor jobs on the LPC to skim ggNtuples and only select interesting events.



** Utils **
This directory contains random useful pieces of code

* FinalPlot.C uses the output of FinalAnalysis.C to produce the final pretty stack plot.
* GridDiemptContamination.C will calculate the signal contamination in the gg diempt distribution from each of the signal mass bins. 
* PrintBinContents.C is useful for debugging - it will print the contents of a histogram or its errors.
* accXeff.C can be used to make the acceptance x efficiency plots for the T5gg or T6gg grids.
