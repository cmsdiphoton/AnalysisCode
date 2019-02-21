// Come back to this directory and do
// $ make
// $ root -b -q -l ana.C
// will produce filename

//Need to first compile SusyEventAnalyzer.cc with any changes with
//root -l
//.L SusyEventAnalyzer.cc+


{
  bool data=false;
  bool powheg=false;
  bool amcatnlo=false;
  bool trig2016 = false;
  bool trig2017 = false;
  bool ZZ = false;
  bool T5Wg = false;
  bool T6Wg = false;
  bool DY = false;
  bool QCD = false;
  bool GJet = true;
  //Will write to output file hist_filename_analysis.root
  TString base = "";
  int num = NNNN;
  TString filename = TString(base) + TString(Form("%d",num));

  TStopwatch ts;
  ts.Start();
  
  R__LOAD_LIBRARY(libSusyEvent.so);
  R__LOAD_LIBRARY(SusyEventAnalyzer_cc.so);

  // chain of inputs
   TChain* chain = new TChain("ggNtuplizer/EventTree");
  
  std::cout<<"Grabbing Files..."<<std::endl;

  TString r = ".root";
  TString s = "";
  TString sB = "";
  TString sC = "";
  TString sD = "";
  TString sE = "";
  TString sF1 = "";
  TString sF2 = "";
  TString sG = "";
  TString sH2 = "";
  TString sH3 = "";

  //extension
  //  if(T5Wg) s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/T5Wg_ntuples_ext2/SMS-T5Wg_mGo2150To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_T5Wg_extv2/170630_145320/0000/ggtree_mc_";
  //main
  if(T5Wg) s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/T5Wg_ntuples/SMS-T5Wg_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_T5Wg/170116_182044/0000/ggtree_mc_";
  //main
  // if(T6Wg) s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/T6Wg_ntuplesv2/SMS-T6Wg_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_T6Wgv2/170630_145433/0000/ggtree_mc_";
  //extension
  if(T6Wg)  s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/T6Wg_ntuples_ext2/SMS-T6Wg_mSq1850To2150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_T6Wg_extv2/170630_145547/0000/ggtree_mc_";

  if(QCD) s = "root://cmsxrootd.fnal.gov//store/group/phys_susy/gpaspala/Summer16_bkg/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/crab_QCD_Pt-40toInf_DoubleEMEnriched/180112_151954/0000/ggtree_mc_";

  if(GJet) s = "root://eoscms.cern.ch//store/group/phys_susy/gpaspala/Summer16_bkg/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/crab_GJet_Pt-40toInf_DoubleEMEnriched/180114_174243/0004/ggtree_mc_";

  if(DY) s = "root://cms-xrd-global.cern.ch//store/group/phys_higgs/itopsisg/lisa/crab/crab3/Spring16_80X/job_spring16_DYJetsToLL_m50_2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_job_spring16_DYJetsToLL_m50/170124_105022/0000/ggtree_mc_";
  if(ZZ) s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/ZZToLLNuNu/Output_Skim_";
  if(amcatnlo) s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_amcatnlo/Output_Skim_";
  if(powheg) s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims/Output_Skim_";

  /*      chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_Leftovers/Output_Skim_0.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_Leftovers/Output_Skim_1.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_Leftovers/Output_Skim_2.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_Leftovers/Output_Skim_3.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_Leftovers/Output_Skim_4.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_Leftovers/Output_Skim_5.root");
    chain->Add("root://cmsxrootd.fnal.gov//store/user/areinsvo/TTSkims_Leftovers/Output_Skim_6.root");*/

  //  if(data) s = "root://cmsxrootd.fnal.gov//store/user/areinsvo/ReMINIAOD_Skims/DoubleEG/Output_Skim_";
  if(data){
    sB  = "root://cmsxrootd.fnal.gov//store/user/lpcggntuples/ggNtuples/13TeV/data/V08_00_26_01/job_DoubleEG_Run2016B_FebReminiAOD/ggtree_data_";
    sC  = "root://cmsxrootd.fnal.gov//store/user/lpcggntuples/ggNtuples/13TeV/data/V08_00_26_01/job_DoubleEG_Run2016C_FebReminiAOD/ggtree_data_";
    sD  = "root://cmsxrootd.fnal.gov//store/user/lpcggntuples/ggNtuples/13TeV/data/V08_00_26_01/job_DoubleEG_Run2016D_FebReminiAOD/ggtree_data_";
    sE  = "root://cmsxrootd.fnal.gov//store/user/lpcggntuples/ggNtuples/13TeV/data/V08_00_26_01/job_DoubleEG_Run2016E_FebReminiAOD/ggtree_data_";
    sF1 = "root://cmsxrootd.fnal.gov//store/user/lpcggntuples/ggNtuples/13TeV/data/V08_00_26_01/job_DoubleEG_Run2016F_FebReminiAOD1/ggtree_data_";
    sF2 = "root://cmsxrootd.fnal.gov//store/user/lpcggntuples/ggNtuples/13TeV/data/V08_00_26_01/job_DoubleEG_Run2016F_FebReminiAOD2/ggtree_data_";
    sG  = "root://cmsxrootd.fnal.gov//store/user/lpcggntuples/ggNtuples/13TeV/data/V08_00_26_01/job_DoubleEG_Run2016G_FebReminiAOD/ggtree_data_";
    sH2 = " root://cmsxrootd.fnal.gov//store/user/lpcggntuples/ggNtuples/13TeV/data/V08_00_26_01/job_DoubleEG_Run2016H_FebReminiAODv2/ggtree_data_";
    sH3 = " root://cmsxrootd.fnal.gov//store/user/lpcggntuples/ggNtuples/13TeV/data/V08_00_26_01/job_DoubleEG_Run2016H_FebReminiAODv3/ggtree_data_";
  }

  int j;
  for(int i = AAAA; i < ZZZZ; i++){
    if(data){
      if(i <= 2109){
	j = i;
	s = sB;
      }
      else if(i <= (2109+696) ){
	s = sC;
	j = i - 2109;
      }
      else if(i <= (2109+696+1167) ){
	s = sD;
	j = i - 2109 - 696;
      }
      else if(i <= (2109+696+1167+991) ){
        s = sE;
	j = i - 2109 - 696 - 1167;
      }
      else if(i <= (2109+696+1167+991+627) ){
        s = sF1;
        j = i - 2109 - 696 - 1167 - 991;
      }
      else if(i <= (2109+696+1167+991+627+97) ){
	s = sF2;
        j = i - 2109 - 696 - 1167 - 991 - 627;
      }
      else if(i <= (2109+696+1167+991+627+97+1708) ){
	s = sG;
        j = i - 2109 - 696 - 1167 - 991 - 627 - 97;
      }
      else if(i <= (2109+696+1167+991+627+97+1708+1849) ){
        s = sH2;
        j = i - 2109 - 696 - 1167 - 991 - 627 - 97 - 1708;
      }
      else if(i <= (2109+696+1167+991+627+97+1708+1849+50) ){
        s = sH3;
        j = i - 2109 - 696 - 1167 - 991 - 627 - 97 - 1708 - 1849;
      }
      else{ cout<< "Woah, out of range!" << endl; continue;}
    }
    else{j = i;};

    if(j < 0){ cout << "Not good. j is negative" << endl; continue;}

    TString infilename = TString(s) + Form("%d",j) +TString(r);
    cout << infilename << endl;
    chain->Add(infilename);
  }
   
  std::cout<<"Files grabbed.  Now to define sea:"<<std::endl;
  SusyEventAnalyzer* sea = new SusyEventAnalyzer(chain);
  std::cout<<"sea defined"<<std::endl;
  // configuration parameters
  // any values given here will replace the default values

  if(powheg)   sea->SetPUFileName("powhegPU.root");
  if(amcatnlo) sea->SetPUFileName("amcatnloPU.root");
  if(ZZ)       sea->SetPUFileName("ZZPU.root");

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
  sea->IncludeAJson("./golden_ReReco.json");
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
