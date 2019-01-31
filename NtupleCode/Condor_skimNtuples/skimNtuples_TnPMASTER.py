import os, sys, ROOT, math

sw = ROOT.TStopwatch()
sw.Start()


chain_in = ROOT.TChain("ggNtuplizer/EventTree") #name of tree and directory the tree is in
chain_in.Add("ZZZZ");

starting_event   =  XXXX
process_n_events =  YYYY
if (starting_event + process_n_events) > chain_in.GetEntries():
    last_event = chain_in.GetEntries()
else:
    last_event = starting_event + process_n_events
#define output name here
file_out = ROOT.TFile("OOOO", "recreate")
dir_out = file_out.mkdir("ggNtuplizer") #make a directory in the output file so the resulting tree has the same structure
dir_out.cd() # tree will be made in this directory

tree_out = chain_in.CloneTree(0)

n_events_saved = 0
n_GoodPho = 0

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
    n_GoodPho = 0
    for i in range(chain_in.nPho):
        if(chain_in.phoIDbit[i]>>1&1 and abs(chain_in.phoEta[i]) < 1.4442):
            n_GoodPho += 1

    if n_GoodPho > 1:
        save_event = True


#Save all interesting TnP events
    if save_event:
        tree_out.Fill()
        n_events_saved += 1

file_out.Write() #Write everything we've created (including the output tree and any histograms) to the file
file_out.Close() 
sw.Stop()

#print automatically appends a new line character
print "Real Time = " + str(sw.RealTime() / 60.0 ) + " minutes."
print "CPU Time = " + str(sw.CpuTime() / 60.0 ) + " minutes."
print "Saved " + str(n_events_saved) + " events."
