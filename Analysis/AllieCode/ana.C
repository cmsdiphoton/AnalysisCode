// Come back to this directory and do
// $ make
// $ root -b -q -l ana.C
// will produce filename

//Need to first compile SusyEventAnalyzer.cc with any changes with
//root -l
//.L SusyEventAnalyzer.cc+


{
  cout << "Here we go." << endl;
  bool GGM = false;
  bool GGM3 = false;
  bool data = false;
  bool powheg=false;
  bool amcatnlo=false;
  bool trig2016 = false;
  bool trig2017 = false;
  bool ZZ = false;
  bool ZGG = false;
  bool DY = false;
  bool DY_mu = false;
  bool data_mu = false;
  bool T5Wg = false;
  bool T5Wg_ext = false;
  bool T6Wg = false;
  bool T6Wg_ext = false;
  bool GJet_LowPt = false;
  bool GJet_HighPt = false;
  bool QCD_LowPt = false;
  bool QCD_HighPt = false;
  bool WGToLNuG = false;
  
  //Will write to output file hist_filename_anaysis.root
  //  TString filename = "powheg_reject7";
  //  TString filename = "data_testFF";
  //  TString filename = "T5Wg_ext_full";
  //  TString filename ="GGM_M1M3_J";//_dEta2Med";
  //  TString filename = "data_reject7_jetpt40";
  //TString filename = "T5Wg_phoJetDR7_Iso1_Med_phoDR6";//1700_500_newcleaning";
  TString filename = "T5gg_2000_100_testMix";
  TStopwatch ts;
  ts.Start();
  
  R__LOAD_LIBRARY(libSusyEvent.so);
  R__LOAD_LIBRARY(SusyEventAnalyzer_cc.so);

  // chain of inputs
   TChain* chain = new TChain("ggNtuplizer/EventTree");
  
  std::cout<<"Grabbing Files..."<<std::endl;

  //  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/T6Wg_ntuplesv2/SMS-T6Wg_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_T6Wgv2/170630_145433/0000/ggtree_mc_46.root");
  //  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/T6Wg_ntuples_ext2/SMS-T6Wg_mSq1850To2150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_T6Wg_extv2/170630_145547/0000/ggtree_mc_19.root");
  //  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/T6Wg_ntuplesv2/SMS-T6Wg_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_T6Wgv2/170630_145433/0000/ggtree_mc_46.root");

  // chain->Add("root://cms-xrd-global.cern.ch//store/group/phys_susy/gpaspala/GGM/GGM_GravitinoLSP_M1-200to1500_M2-200to1500_TuneCUETP8M1_13TeV_pythia8/crab_GGM_M1M2/180219_145514/0000/ggtree_mc_1.root");
  //  chain->Add("/uscms/home/ahall/nobackup/Diphoton/T5Wg_Files/T5gg_2000_1000.root");
  chain->Add("/uscms/home/ahall/nobackup/Diphoton/T5Wg_Files/T5Wg_mGlu_2000_mNeu_100.root"); 
  //  chain->Add("root://xrootd.unl.edu//store/group/phys_susy/gpaspala/Spring_80X_bkg/job_spring16_WGToLNuG/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_job_spring16_WGToLNuG/170524_113730/0000/ggtree_mc_1.root");

  if(ZGG){
    /*    stringstream ss;
    for(int i=1; i<=25; i++){
      string r = ".root";
      ss.str("");
      string s = "root://xrootd.unl.edu//store/group/phys_susy/gpaspala/Spring_80X_bkg/ZGGToNuNuGG/ZGGToNuNuGG_5f_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_job_spring16_ZGGToNuNuGG/171006_094702/0000/ggtree_mc_"; //1 to 25
    ss << s;
    ss << i;
    ss << r;
    const char* infilename = ss.str().c_str();
    TString num = "";
    num.Form("%d",i);
    TString test = TString(s) + Form("%d",i) +TString(r);
    chain->Add(test);
    }*/

    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/ZGGToNuNuGG/Output_Skim_0.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/ZGGToNuNuGG/Output_Skim_1.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/ZGGToNuNuGG/Output_Skim_2.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/ZGGToNuNuGG/Output_Skim_3.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/ZGGToNuNuGG/Output_Skim_4.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/ZGGToNuNuGG_large/Output_Skim_0.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/ZGGToNuNuGG_large/Output_Skim_1.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/ZGGToNuNuGG_large/Output_Skim_2.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/ZGGToNuNuGG_large/Output_Skim_3.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/ZGGToNuNuGG_large/Output_Skim_4.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/ZGGToNuNuGG_large/Output_Skim_5.root");
  }
  if(GGM3) { //goes from 0 to 476
    for(int i=450; i<=476; i++){
      string s = "root://cmseos.fnal.gov//store/user/areinsvo/2016MC/GGM_M1M3_Skims/Output_Skim_";
      string r = ".root";
      TString test = TString(s) + Form("%d",i) +TString(r);
      chain->Add(test);
    }
  }

  if(GGM) {
    for(int i=0; i<=140; i++){
      string s = "root://cmseos.fnal.gov//store/user/areinsvo/2016MC/GGM_M1M2_Skims/Output_Skim_";
      string r = ".root";
      TString test = TString(s) + Form("%d",i) +TString(r);
      chain->Add(test);
    }
  }
  if(T6Wg){
    stringstream ss;
    for(int i=1; i<=46; i++){
      string r = ".root";
      ss.str("");
    string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/T6Wg_ntuplesv2/SMS-T6Wg_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_T6Wgv2/170630_145433/0000/ggtree_mc_";
    ss << s;
    ss << i;
    ss << r;
    const char* infilename = ss.str().c_str();
    TString num = "";
    num.Form("%d",i);
    TString test = TString(s) + Form("%d",i) +TString(r);
    chain->Add(test);
    }
  }
  if(T6Wg_ext){
    stringstream ss;
    for(int i=1; i<=19; i++){
      string r = ".root";
      ss.str("");
    string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/T6Wg_ntuples_ext2/SMS-T6Wg_mSq1850To2150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_T6Wg_extv2/170630_145547/0000/ggtree_mc_";
    ss << s;
    ss << i;
    ss << r;
    const char* infilename = ss.str().c_str();
    TString num = "";
    num.Form("%d",i);
    TString test = TString(s) + Form("%d",i) +TString(r);
    chain->Add(test);
    }
  }
  if(T5Wg){//goes from 1 to 80
    stringstream ss;
    for(int i=1; i<=80; i++){
      string r = ".root";
      ss.str("");
    string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/T5Wg_ntuples/SMS-T5Wg_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_T5Wg/170116_182044/0000/ggtree_mc_";
    ss << s;
    ss << i;
    ss << r;
    const char* infilename = ss.str().c_str();
    TString num = "";
    num.Form("%d",i);
    TString test = TString(s) + Form("%d",i) +TString(r);
    //    cout << test << endl;
    chain->Add(test);
    }
  }
  if(T5Wg_ext){
    stringstream ss;
    for(int i=1; i<=16; i++){
      string r = ".root";
      ss.str("");
      string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/T5Wg_ntuples_ext2/SMS-T5Wg_mGo2150To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_T5Wg_extv2/170630_145320/0000/ggtree_mc_";
      ss << s;
      ss << i;
      ss << r;
      const char* infilename = ss.str().c_str();
      TString num = "";
      num.Form("%d",i);
      TString test = TString(s) + Form("%d",i) +TString(r);
      chain->Add(test);
    }
  }
  if(DY_mu){
    stringstream ss;
    for(int i=0; i<=12; i++){
      string r = ".root";
      ss.str("");
      string s = "root://cms-xrd-global.cern.ch//store/user/areinsvo/2016MC/DYJets_MuSkim/even/Output_Skim_";
      ss << s;
      ss << i;
      ss << r;
      const char* infilename = ss.str().c_str();
      chain->Add(infilename);
      s = "root://cms-xrd-global.cern.ch//store/user/areinsvo/2016MC/DYJets_MuSkim/odd/Output_Skim_";
      ss.str("");
      ss << s;
      ss << i;
      ss << r;
      infilename = ss.str().c_str();
      chain->Add(infilename);
    }
  }
  if(DY){
    stringstream ss;
    for(int i=1; i<=395; i++){
      string r = ".root";
      ss.str("");
      string s = "root://cms-xrd-global.cern.ch//store/group/phys_higgs/itopsisg/lisa/crab/crab3/Spring16_80X/job_spring16_DYJetsToLL_m50_2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_job_spring16_DYJetsToLL_m50/170124_105022/0000/ggtree_mc_";
      ss << s;
      ss << i;
      ss << r;
      const char* infilename = ss.str().c_str();
      chain->Add(infilename);
    }
  }
  if(ZZ){
    stringstream ss;
    for(int i=0; i<=23; i++){
      string r = ".root";
      ss.str("");
      string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/ZZToLLNuNu/Output_Skim_";
      ss << s;
      ss << i;
      ss << r;
      const char* infilename = ss.str().c_str();
      chain->Add(infilename);
    }
  }
  if(GJet_HighPt){
    /* stringstream ss;
    for(int i=0; i<=1998; i++){
      if(i > 324 && i < 506) continue;
      string r = ".root";
      ss.str("");
      string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet40_diEM/Output_Skim_";
      ss << s;
      ss << i;
      ss << r;
      string intermediate = ss.str();
      const char* infilename = intermediate.c_str();    
      cout << infilename << "*************" << intermediate <<  endl;
      chain->Add(infilename);
    }
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet40_diEM/Output_Skim_326.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet40_diEM/Output_Skim_327.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet40_diEM/Output_Skim_340.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet40_diEM/Output_Skim_377.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet40_diEM/Output_Skim_381.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet40_diEM/Output_Skim_428.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet40_diEM/Output_Skim_444.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet40_diEM/Output_Skim_446.root");
    stringstream ss2;
    for(int i=0; i<=172; i++){
      string r = ".root";
      ss2.str("");
      string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet40_diEM/Recovery/Output_Skim_";
      ss2 << s;
      ss2 << i;
      ss2 << r;
      string intermediate = ss2.str();
      const char* infilename = intermediate.c_str();
      cout << infilename << " -------- " << intermediate << endl;
      chain->Add(infilename);
      }*/
    for(int i=100; i<=108; i++){ 
      string r = ".root";
      string s = "root://xrootd.unl.edu//store/group/phys_susy/gpaspala/Summer16_bkg/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/crab_GJet_Pt-40toInf_DoubleEMEnriched/180114_174243/0000/ggtree_mc_";
      TString test = TString(s) + Form("%d",i) +TString(r);
      chain->Add(test);
    }
  }
  if(GJet_LowPt){
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet20to40_diEM/GJetLowPtATakeSix.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet20to40_diEM/GJetLowPtATakeFour.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet20to40_diEM/GJetLowPtB.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet20to40_diEM/GJetLowPtBTakeTwo.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet20to40_diEM/GJetLowPtC.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet20to40_diEM/GJetLowPtCTakeTwo.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/GJet20to40_diEM/GJetLowPtD.root");
  }
  if(QCD_HighPt){
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/QCD40_diEM/QCDHighPtA.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/QCD40_diEM/QCDHighPtB.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/QCD40_diEM/QCDHighPtARecovery.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/QCD40_diEM/QCDHighPtC.root");

  }
  if(QCD_LowPt){
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/QCD30to40_diEM/QCDLowPtA.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/QCD30to40_diEM/QCDLowPtB.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/QCD30to40_diEM/QCDLowPtBRecovery.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/QCD30to40_diEM/QCDLowPtC.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/QCD30to40_diEM/QCDLowPtCRecovery.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/QCD30to40_diEM/QCDLowPtD.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/QCD30to40_diEM/QCDLowPtDRecovery.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/QCD30to40_diEM/QCDLowPtE.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/QCD30to40_diEM/QCDLowPtF.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/QCD30to40_diEM/QCDLowPtFRecovery.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/QCD30to40_diEM/QCDLowPtG.root");
    chain->Add("root://cmseos.fnal.gov//store/user/areinsvo/2016MC/Summer16/QCD30to40_diEM/QCDLowPtGRecovery.root");

  }
  if(WGToLNuG){
    stringstream ss;
    for(int i=0; i<=14; i++){
      string r = ".root";
      ss.str("");
      string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/2016MC/Summer16/WGToLNuG/Original/Output_Skim_";
      ss << s;
      ss << i;
      ss << r;
      string intermediate = ss.str();
      const char* infilename = intermediate.c_str();
      chain->Add(infilename);
    }
    stringstream ss2;
    for(int i=0; i<=61; i++){
      string r = ".root";
      ss2.str("");
      string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/2016MC/Summer16/WGToLNuG/Ext/Output_Skim_";
      ss2 << s;
      ss2 << i;
      ss2 << r;
      string intermediate = ss2.str();
      const char* infilename = intermediate.c_str();
      chain->Add(infilename);
    }
  }
  if(amcatnlo){
    stringstream ss;
    for(int i=0; i<=32; i++){
      string r = ".root";
      ss.str("");
      string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_amcatnlo/Output_Skim_";
      ss << s;
      ss << i;
      ss << r;
      const char* infilename = ss.str().c_str();
      chain->Add(infilename);
    }
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_amcatnlo_leftover/Output_Skim_0.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_amcatnlo_leftover/Output_Skim_1.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_amcatnlo_leftover/Output_Skim_2.root");
  }

  //  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims/Output.root");

  if(powheg){
    stringstream ss;                                                      
    for(int i=0; i<=39; i++){
      string r = ".root"; 
      ss.str("");
      string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims/Output_Skim_";
      ss << s;
      ss << i;
      ss << r; 
      const char* infilename = ss.str().c_str();  
      chain->Add(infilename);
    }

    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_Leftovers/Output_Skim_0.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_Leftovers/Output_Skim_1.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_Leftovers/Output_Skim_2.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_Leftovers/Output_Skim_3.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_Leftovers/Output_Skim_4.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_Leftovers/Output_Skim_5.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_Leftovers/Output_Skim_6.root");
  }
  //  chain->Add("ggtree_mc_3.root");

  //  chain->Add("root://xrootd.unl.edu//store/group/phys_higgs/itopsisg/lisa/crab/crab3/Spring16_80X/job_spring16_TTJets_powheg/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_job_spring16_TTJets_powheg/170815_152458/0000/ggtree_mc_3.root");

  //  chain->Add("/uscms_data/d3/areinsvo/T5Wg/T5gg_2000_1000.root");     
  //  chain->Add("/uscms_data/d3/areinsvo/T5Wg/T5Wg_mGlu_2000_mNeu_100.root");
  //  chain->Add("/uscms_data/d3/areinsvo/T5Wg/T5Wg_mGlu_2000_mNeu_1900.root");

  //    chain->Add("/uscms_data/d1/areinsvo/T5Wg/T5gg_2000_1000.root");
  // chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/T5Wg_mGlu_2000_mNeu_1000");

  if(trig2017){
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/2017TrigSkims/Output_Skim_0.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/2017TrigSkims/Output_Skim_1.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/2017TrigSkims/Output_Skim_2.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/2017TrigSkims/Output_Skim_3.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/2017TrigSkims/Output_Skim_0.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/2017TrigSkims/Output_Skim_1.root");
  }
  if(trig2016){
    stringstream ss; 
    for(int i=1; i<=42; i=i+2){
      if(i == 16) continue;
      string r = ".root";     
      ss.str("");             
      string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_";
      ss << s;
      ss << i;
      ss << r;
      const char* infilename = ss.str().c_str();
      chain->Add(infilename);
    }
  }
  /*  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_0.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_1.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_2.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_7.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_8.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_9.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_10.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_11.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_13.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_15.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_16.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_17.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_18.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_20.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_21.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_24.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_26.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_28.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_29.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_30.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_36.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_38.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_39.root");
  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TriggerSkims/SingleElectron/RunCDEF/Output_Skim_40.root");;*/


  if(data_mu){
    stringstream ss;
    for(int i=0; i<=21; i++){
      string r = ".root";
      ss.str("");
      string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/2016Data/DoubleMuon/Run2016B/0000/Output_Skim_";
      ss << s;
      ss << i;
      ss << r;
      const char* infilename = ss.str().c_str();
      chain->Add(infilename);
    }
    //    chain->Add("./ggtree_data.root");                                                                                                                      
  }
  if(data){   
    /*    stringstream ss;
    for(int i=0; i<=103; i++){ 
      //for(int i=91; i<=91; i++){  
      //  re mini AOD data goes from 0 to 103 (inclusive)
      string r = ".root";
      ss.str("");
      string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/ReMINIAOD_Skims/DoubleEG/RequireTrig/Output_Skim_";
      //      string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/ReMINIAOD_Skims/DoubleEG/Output_Skim_";
      ss << s;
      ss << i; 
      ss << r; 
      const char* infilename = ss.str().c_str();
      chain->Add(infilename);
      }*/
    for(int i=0; i<=103; i++){
      //  re mini AOD data goes from 0 to 103 (inclusive)  
      string r = ".root";
      //      string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/Output_Skim_";
      string s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/ReMINIAOD_Skims/DoubleEG/RequireTrig/Output_Skim_";
      TString test = TString(s) + Form("%d",i) +TString(r);
      chain->Add(test);
    }
    //    chain->Add("./ggtree_data.root");
  }
  //  chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/ReMINIAOD_Skims/DoubleEG/Output_Skim_4.root");

   
  std::cout<<"Files grabbed.  Now to define sea:"<<std::endl;
  SusyEventAnalyzer* sea = new SusyEventAnalyzer(chain);
  std::cout<<"sea defined"<<std::endl;
  // configuration parameters
  // any values given here will replace the default values

    sea->SetDoSignal(true);

  if(WGToLNuG) sea->SetPUFileName("WGPU.root");
  if(powheg)   sea->SetPUFileName("powhegPU.root");
  if(amcatnlo) sea->SetPUFileName("amcatnloPU.root");
  if(ZZ)       sea->SetPUFileName("ZZPU.root");
  if(ZGG)      sea->SetPUFileName("ZggPU.root");
  if(GJet_LowPt || GJet_HighPt || QCD_LowPt || QCD_HighPt) sea->SetPUFileName("GJetPU.root");
  if(GGM || GGM3) sea->SetDoSignal(true); 

  if(T5Wg or T5Wg_ext ){
    sea->SetDoSignal(true);
    sea->SetDoT5Wg(true);
  }
  if(T6Wg or T6Wg_ext){
    sea->SetDoSignal(true);
    sea->SetDoT5Wg(false);
  }

  sea->SetDataset(filename);        // dataset name

  bool tt = (powheg || amcatnlo);
  std::cout<< "doTTBar? " << tt << std::endl;
  sea->SetDoTTBar(tt);
  sea->SetDoToys(false); //Generate toy histograms for diempt systematic
  sea->SetPrintInterval(2.5e5);             // print frequency
  sea->SetPrintLevel(0);                  // print level for event contents
  sea->SetOutputEventNumbers(false);      // print run and event numbers
  sea->SetUseTrigger(true);
  sea->SetUseJSON(true);
  sea->DoRhoCorrection(true);
  sea->DoNvertexCorrection(false);
  sea->SetDR03Rho25Corr(0.081,0.022,0.011);//Ecal,Hcal.Track
  sea->SetPFisoRho25Corr(0.031,0.013,0.078);//chargedHadronIso,neutralHadronIso,photonIso
  sea->SetFilter(false);                  // filter events passing final cuts
  sea->isFastSim(true);    //only matters for HiggsAna()
  sea->isFullSim(false);   //only matters for HiggsAna()
  sea->SetProcessNEvents(-1);//9000000 );//5000000);     //(set to -1 to process all events) number of events to be processed
  std::cout<<"Done seting variables"<<endl;
  // HLT trigger path names
  //sea->AddHltName("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v");
  //sea->AddHltName("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v");
  //sea->AddHltName("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v");
  //sea->AddHltName("HLT_Photon36_R9Id85_Photon22_R9Id85_v");
  

  //fully combined json
  //  sea->IncludeAJson("/uscms_data/d3/areinsvo/CMSSW_8_0_24_patch1/src/SusyAnalysis/SusyCode/macro/Trig2017.json");
  sea->IncludeAJson("/uscms_data/d3/areinsvo/CMSSW_8_0_24_patch1/src/SusyAnalysis/SusyCode/macro/golden_ReReco.json");
  //  sea->IncludeAJson("/uscms_data/d3/areinsvo/CMSSW_8_0_24_patch1/src/SusyAnalysis/SusyCode/macro/golden.json");
  //  sea->IncludeAJson("/uscms_data/d3/areinsvo/CMSSW_8_0_24_patch1/src/SusyAnalysis/SusyCode/macro/EraC.json");
  sea->Loop();
  //sea->DR03();
  //sea->Pileup();
  //sea->Filter();
  //sea->PhotonId();
  //sea->HiggsAna();

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
