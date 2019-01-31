import os, sys, ROOT, math
#from array import *
sw = ROOT.TStopwatch()
sw.Start()


chain_in = ROOT.TChain("ggNtuplizer/EventTree") #name of tree and directory the tree is in
#chain_in.Add("root://eoscms.cern.ch//eos/cms/store/group/phys_smp/ggNtuples/13TeV/mc/job_spring15_gjet_pt20_MGG_40to80_25ns.root") #file names. Can use wildcard * !
chain_in.Add("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/itopsisg/lisa/crab/crab3/T5Wg_signal/SMS-T5Wg_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_T5Wg/161024_164050/0000/ggtree_mc_*.root")

#Trimming
#Don't save any of the tau variables
#chain_in.SetBranchStatus("tau*",0)
#chain_in.SetBranchStatus("pfTau*",0)
#chain_in.SetBranchStatus("",0)

starting_event =           0000
process_n_events =     20000000
if (starting_event + process_n_events) > chain_in.GetEntries():
    last_event = chain_in.GetEntries()
else:
    last_event = starting_event + process_n_events
#print chain_in.nPho
#define output file names
base = "./T5gg_sorted_ntuples/T5gg_mGlu_"
mid = "_mNeu_"
end = ".root"
#define mass grid here
min_gluino_mass = 1800
max_gluino_mass = 1851
gluino_mass_step = 50
min_neu_mass = 100
neu_mass_step = 100

#Define output files
file_num = 0
file_list = []
dir_list = []
tree_list = []
for glu_mass in xrange(min_gluino_mass,max_gluino_mass,gluino_mass_step):
    for neu_mass in range(min_neu_mass,glu_mass+1,neu_mass_step):
        out_file_name = base + str(glu_mass) +mid+str(neu_mass) +end
        print "File number " + str(file_num) + ": " + out_file_name
        temp_file = ROOT.TFile(out_file_name, "recreate")
        file_list.append(temp_file)
        temp_dir = file_list[file_num].mkdir("ggNtuplizer")
        dir_list.append(temp_dir)
        dir_list[file_num].cd()
        temp_tree = chain_in.CloneTree(0)
        tree_list.append(temp_tree)
        file_num=file_num+1
total_num_files = file_num-1
#Histogram example
#h1 = TH1F("parameters","title",10,0,100)

for j_entry in range(starting_event,last_event):
    i_entry = chain_in.LoadTree(j_entry)

    if i_entry < 0:
        break

    nb = chain_in.GetEntry(j_entry)
    if nb <= 0:
        continue
    
    if j_entry % 10000 ==0:
        print "Processing entry " + str(j_entry) + " of " + str(last_event)
    true_glu_mass = 0
    true_neu_mass = 0
    tag = chain_in.EventTag
    str_tag = str(tag)

    true_glu_mass = int(str_tag[13:17])
#    print true_glu_mass
    true_neu_mass = int(str_tag[18:])
#    print true_neu_mass

    if true_glu_mass ==0:
        print "Gluino mass = 0. Something is wrong!"
 #   if true_glu_mass < 1380 or true_glu_mass > 1570:
 #       print "Gluino mass = " + str(true_glu_mass) + " which is out of bounds!"
    if true_neu_mass ==0:
        print "Neutralino mass = 0. Something is wrong!"
#    if true_neu_mass > 1520:
#        print "Uh oh. Neutralino mass = " + str(true_neu_mass)
    file_it = -1
    file_num =0
    for glu_mass in xrange(min_gluino_mass,max_gluino_mass,gluino_mass_step):
        if(file_it == -1):  #only need to continue if file_it has not already been set to something
            for neu_mass in range(min_neu_mass,glu_mass+1,neu_mass_step):
                if(true_glu_mass < (glu_mass + 20) and true_glu_mass > (glu_mass - 20) and
                   true_neu_mass < (neu_mass + 49) and true_neu_mass > (neu_mass - 49) ):
                    file_it = file_num
                    break
                else:
                    file_num+=1

    if file_it != -1:
        dir_list[file_it].cd()
        tree_list[file_it].Fill()
    else:
        print "No match! You done messed up. Gluino mass = " + str(true_glu_mass) + " and neutralino mass = " + str(true_neu_mass)

for i in xrange(0,total_num_files+1):
    file_list[i].Write()
    file_list[i].Close()

sw.Stop()

#print automatically appends a new line character
print "Real Time = " + str(sw.RealTime() / 60.0 ) + " minutes."
print "CPU Time = " + str(sw.CpuTime() / 60.0 ) + " minutes."
