# AnalysisCode

** Analysis **
This directory contains the necessary code to make the analysis histograms from ggNtuples.

* AllieCode/ana.C defines what files will be used.
* AllieCode/ggEventTree.h defines the variables in the ggNtuples that will be used by the analyzer
* AllieCode/SusyEventAnalyzer.cc is the analysis code that sorts the events into signal and control regions and makes all the histograms.
* Diempt.C will take the input gg and ff or ee histograms of diempt and calculate the diempt ratios which can be used to reweight the MET distributions. This macro also makes a nice plot of the diempt ratio and calculates the appropriate errors.
* FinalAnalysis.C takes the output of the event loop and calculates the final background estimates and errors and makes several plots of the final quantities.
* GGM_Combination/finishHists_*C calculates the uncertainties that we implemented for the GGM combination paper.
* lester_mt2_bisect.h is the code to calculate MT2.
* ratioMethod.C performs the fit to the gg/ff ratio.



** Limits **
Code to make limit plots from the datacards.

* MakeCountingFiles.ipynb is a jupyter notebook for creation of the datacards and scripts needed to run combine
* MakeLimitPlots.ipynb is a jupyter notebook for creation of limit plots from the counting files created by combine
* mlfit.ipynb and .py creates diagnostic plots of the limit calculation
* PlotsSMSFiles is a directory that contains all the files you need to edit to use the PlotsSMS package to make nice looking limit plots.



** NtupleCode **
Code to manipulate the ggNtuples.

* Condor_skimNtuples contains example code to submit condor jobs on the LPC to skim ggNtuples and only select interesting events.
* sortSignalByMass.py will take the events in a ggNtuple and sort it by squark/gluino vs neutralino mass.


** Utils **
This directory contains random useful pieces of code

* FinalPlot.C uses the output of FinalAnalysis.C to produce the final pretty stack plot.
* GridDiemptContamination.C will calculate the signal contamination in the gg diempt distribution from each of the signal mass bins. 
* PrintBinContents.C is useful for debugging - it will print the contents of a histogram or its errors.
* accXeff.C can be used to make the acceptance x efficiency plots for the T5gg or T6gg grids.
