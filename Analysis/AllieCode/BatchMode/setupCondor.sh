#! /bin/bash                                                                            
rm Analysis.tar.gz

mkdir Analysis
mkdir Analysis/src
mkdir Analysis/macro
cp ../src/SusyNtuplizer_LinkDef.h Analysis/src/SusyNtuplizer_LinkDef.h
cp ../src/ggEventTree.h Analysis/src/ggEventTree.h
cp golden_ReReco.json Analysis/macro
cp SusyEventAnalyzer.h Analysis/macro
cp SusyEventAnalyzer.cc Analysis/macro
cp SusyEventAnalyzer_cc.so Analysis/macro
cp SusyEventAnalyzer_cc_ACLiC_dict_rdict.pcm Analysis/macro
cp Event_dict_rdict.pcm Analysis/macro
cp libSusyEvent.so Analysis/macro
cp reweights.root Analysis/macro
cp NJet_reweights.root Analysis/macro
cp powhegPU.root Analysis/macro
cp lester_mt2_bisect.h Analysis/macro
cp MyDataPileupHistogram.root Analysis/macro

tar -cvzf Analysis.tar.gz Analysis

rm Analysis/macro/*
rm Analysis/src/*
rmdir Analysis/macro
rmdir Analysis/src
rmdir Analysis