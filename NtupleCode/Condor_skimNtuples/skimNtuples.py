#Allie Reinsvold-Hall
#May 2016
#Script to skim ggNtuples and only keep the events that pass the criteria specified
#Usage: adjust input file names, starting_event, and process_n_event parameters below, then run using python skimNtuples.py

import os, sys, ROOT, math

sw = ROOT.TStopwatch()
sw.Start()


chain_in = ROOT.TChain("ggNtuplizer/EventTree") #name of tree and directory the tree is in
#Add the names of the files you want to skim here. Can use wildcard * character.
chain_in.Add("root://cmsxrootd.fnal.gov//store/group/phys_higgs/itopsisg/lisa/crab/crab3/80X_bkg/GJ20to40_2/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/crab_GJ20to40_2/160518_083428/0000/ggtree_mc_1.root");


starting_event   = 0
#Number of events to process. This is a hack, since this script isn't set up to run using condor
#For large files, I run this script several times, changing starting_event each time, 
#and then at the end merge the resulting output files
process_n_events =  1000
if (starting_event + process_n_events) > chain_in.GetEntries():
    last_event = chain_in.GetEntries()
else:
    last_event = starting_event + process_n_events
print chain_in.nPho
#define output name here
file_out = ROOT.TFile("test.root", "recreate")
dir_out = file_out.mkdir("ggNtuplizer") #make a directory in the output file so the resulting tree has the same structure
dir_out.cd() # tree will be made in this directory
# Number tells it how many events to clone into the new tree. We use 0 here so it makes the same structure 
#but doesn't copy in any of the events. Putting no number there would copy all of the events
tree_out = chain_in.CloneTree(0)

n_events_saved = 0
n_TnP = 0

for j_entry in range(starting_event,last_event):
    i_entry = chain_in.LoadTree(j_entry)

    if i_entry < 0:
        break

    nb = chain_in.GetEntry(j_entry)
    if nb <= 0:
        continue
    
    if j_entry % 10000 ==0:
        print "Processing entry " + str(j_entry) + " of " + str(last_event)

    save_event = False
    #Define criteria for saving the event here. In this case it is set to save tag and probe pairs with 
    #invariant mass greater than 70 GeV 
    ele_tags = []
    pho_probes = []

    for i in range(chain_in.nEle):
        if(chain_in.eleIDbit[i]>>1&1 and chain_in.elePt[i] > 30):
            ele_tags.append(i)

    for i in range(chain_in.nPho):
        if(chain_in.phoIDbit[i]>>1&1):
            pho_probes.append(i)

    for it_tag in ele_tags:
        for it_probe in pho_probes:
            invmass = math.sqrt(2*chain_in.phoEt[it_probe]*chain_in.elePt[it_tag]*(math.cosh(chain_in.phoEta[it_probe]-chain_in.eleEta[it_tag]) - math.cos(chain_in.phoPhi[it_probe]-chain_in.elePhi[it_tag]) ) )
            if(invmass >70):
                save_event = True
                n_TnP += 1

#Save all interesting events
    if save_event:
        tree_out.Fill()
        n_events_saved += 1

file_out.Write() #Write everything we've created (including the output tree and any histograms) to the file
file_out.Close() 
sw.Stop()

print "Real Time = " + str(sw.RealTime() / 60.0 ) + " minutes."
print "CPU Time = " + str(sw.CpuTime() / 60.0 ) + " minutes."
print "Saved " + str(n_events_saved) + " events."
print "Total of " + str(n_TnP) + " tag and probe pairs."
