//
//  Susyeventanalyzer.cc
//
//
//  Created by Allie Reinsvold Hall on 8/4/14.
//
// Stripped version with fewer unnecessary histograms
// Includes full analysis code 
// Updated September 2016 
 
#define SusyEventAnalyzer_cxx

#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TGraphAsymmErrors.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <map>
#include <set>
#include <cmath>
#include <utility>
#include <fstream>
#include "lester_mt2_bisect.h"
#include <unordered_map>

#include "SusyEventAnalyzer.h"

using namespace std;

void SusyEventAnalyzer::Loop() {
  bool version74X = false;
  cout<<"Inside Loop(). WOOOOOOO!"<<endl;
  cout<<"Did you remember to compile any changes to SusyEventAnalyzer.cc?"<<endl;
  if (fChain == 0) {
    cout << "fChain ==0"<< endl;
    return;
  }
  Long64_t nentries = fChain->GetEntries();
  cout << "here: " << nentries << endl;
  new ggEventTree(fChain);
  cout << "here" << endl;
  event = new ggEventTree(fChain);
  cout << "total events in files  : " << nentries << endl;
  if(processNEvents <= 0 || processNEvents > nentries) processNEvents = nentries;
    
  cout << "events to be processed : " << processNEvents << endl;
  //array to keep track of cut flow of events
  const int NCNT = 20;
  int nCnt[NCNT];
  for(int i=0; i<NCNT; i++) nCnt[i] = 0;

  fstream myfile;
  myfile.open ("allie.txt");

  fstream signalfile;
  signalfile.open ("ggMET100.txt");
    
  if(printLevel>0) cout<<"Open hist file" << endl;
  TFile* fout = new TFile("hist_"+ds+".root","RECREATE");
  if(printLevel>0) cout<<"Hist file opened" << endl;
    
  fout->cd();

  TTree * tree = new TTree("myTree","Skims");

  bool isGG = false;
  bool isFF = false;
  bool isFF_anyR9 = false;
  float met = 0;
  float HT  = 0;
  float Rho = 0;
  float NVertex = 0;
  double diEMpt = 0;   
  int NJets = 0;
  Long64_t  eventNum = 0;
  float trailPt = 0.0;
  float leadPt = 0.0;
  float leadR9 = 0.0;
  float trailR9 = 0.0;

  tree->Branch("trailR9", &trailR9);
  tree->Branch("leadR9", &leadR9);
  tree->Branch("trailPt", &trailPt);
  tree->Branch("leadPt", &leadPt);
  tree->Branch("event",&eventNum);
  tree->Branch("met", &met);
  tree->Branch("isGG", &isGG);
  tree->Branch("isFF", &isFF);
  tree->Branch("isFF_anyR9", &isFF_anyR9);
  tree->Branch("HT", &HT);
  tree->Branch("Rho", &Rho);
  tree->Branch("NVertex", &NVertex);
  tree->Branch("NJet", &NJets);
  tree->Branch("Diempt", &diEMpt);

    
  //Define histograms here
  //....
  TH1F * eeweights;
  TH1F * gfweights;
  TH1F * ffweights;
  TH1F* puweights;
  
  //  const int numMetBins = 27;
  //Double_t MetBins[numMetBins+1]  = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,75,85,100,115,130,150,185,250,325,425,550,700};

  const int numMetBins = 24;                                                                                                                                                                                                              
  Double_t MetBins[numMetBins+1]  = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,75,85,100,115,130,150,185,250,350}; 

  TH1F* h_puTrue = new TH1F("puTrue","Number of PU interactions; Number; Events",100,0.,100.);
  TH1F* h_puTrue_ee = new TH1F("puTrue_ee","Number of PU interactions in ee events; Number; Events",100,0.,100.);

  TH1I* h_NumMuons = new TH1I("NumMuons","Number of muons that pass quality cuts;;Events",20,0,20);
  //-------------Electron study- how many electrons are we missing by only using the PF Photon collection? This matters for HT
  TH1I* h_NumPFElectrons = new TH1I("NumPFElectrons","Number of PF electrons that don't overlap with another analysis object;;Events",20,0,20);
  TH1F* h_ElectronMissingHT = new TH1F("ElectronMissingHT","Total HT missing from the event by not including PF Electrons;HT (GeV);Events",400,0.,1000);
  TH1F* h_PFElectronPT = new TH1F("PFElectronPT","PT of PF Electrons in event;PT (GeV);Events",200,0.,1000.);
  //-------------End of Electron study------------------------------

  TH1I* h_FakeCuts = new TH1I("FakeCuts","no cuts,#eta,E_{T},hOverE,pixelCut,IsoUpperCut,!(CHiso || sihih),!CHiso,!sihih",10,0,10);

  TH2F* h_signalMass = new TH2F("signalMass","Mass grid for signal sample",6000,1300,2500,500,0,2500);
  TH2F* h_signalEGMass = new TH2F("signalEGMass","Number of eg events per mass point",6000,1300,2500,500,0,2500);

  TH1F* h_ggMetSig = new TH1F("ggMetSig",";Met Significance;Events",500,0.,50.);
  TH2F* h_ggMetSigVsHT = new TH2F("ggMetSigVsHT",";HT;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_ggMetSigVsST = new TH2F("ggMetSigVsST",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_ggMetSigVsMET = new TH2F("ggMetSigVsMET",";MET;Met Significance",800,0.,800.,500,0.,50.);

  TH1F* h_gammafakeMetSig = new TH1F("gammafakeMetSig",";Met Significance;Events",500,0.,50.);
  TH2F* h_gammafakeMetSigVsHT = new TH2F("gammafakeMetSigVsHT",";HT;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_gammafakeMetSigVsST = new TH2F("gammafakeMetSigVsST",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_gammafakeMetSigVsMET = new TH2F("gammafakeMetSigVsMET",";MET;Met Significance",800,0.,800.,500,0.,50.);

  TH2F* h_gammafakeMetSigVsST_0Jets = new TH2F("gammafakeMetSigVsST_0Jets",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_gammafakeMetSigVsST_1Jets = new TH2F("gammafakeMetSigVsST_1Jets",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_gammafakeMetSigVsST_2Jets = new TH2F("gammafakeMetSigVsST_2Jets",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_gammafakeMetSigVsST_3pJets = new TH2F("gammafakeMetSigVsST_3pJets",";ST;Met Significance",500,0.,5000.,500,0.,50.);

  TH1F* h_ffMetSig = new TH1F("ffMetSig",";Met Significance;Events",500,0.,50.);
  TH2F* h_ffMetSigVsMET = new TH2F("ffMetSigVsMET",";MET;Met Significance",800,0.,800.,500,0.,50.);
  TH2F* h_ffMetSigVsHT = new TH2F("ffMetSigVsHT",";HT;Met Significance",500,0.,5000.,500,0.,50.);   
  TH2F* h_ffMetSigVsST = new TH2F("ffMetSigVsST",";ST;Met Significance",500,0.,5000.,500,0.,50.);

  TH1F* h_egMetSig = new TH1F("egMetSig",";Met Significance;Events",500,0.,50.);
  TH2F* h_egMetSigVsMET = new TH2F("egMetSigVsMET",";MET;Met Significance",800,0.,800.,500,0.,50.);
  TH2F* h_egMetSigVsHT = new TH2F("egMetSigVsHT",";HT;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_egMetSigVsST = new TH2F("egMetSigVsST",";ST;Met Significance",500,0.,5000.,500,0.,50.);

  TH2F* h_egMetSigVsST_1Jets = new TH2F("egMetSigVsST_1Jets",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_egMetSigVsST_2Jets = new TH2F("egMetSigVsST_2Jets",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_egMetSigVsST_3pJets = new TH2F("egMetSigVsST_3pJets",";ST;Met Significance",500,0.,5000.,500,0.,50.);

  TH1F* h_eeMetSig_ST150400 = new TH1F("eeMetSig_ST150400",";Met Significance;Events",500,0.,50.);
  TH1F* h_eeMetSig_ST400 = new TH1F("eeMetSig_ST400",";Met Significance;Events",500,0.,50.);

  TH1F* h_eeMetSig = new TH1F("eeMetSig",";Met Significance;Events",500,0.,50.);
  TH2F* h_eeMetSigVsMET = new TH2F("eeMetSigVsMET",";MET;Met Significance",800,0.,800.,500,0.,50.);
  TH2F* h_eeMetSigVsHT = new TH2F("eeMetSigVsHT",";HT;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_eeMetSigVsST = new TH2F("eeMetSigVsST",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  Double_t MetSigBins[22] = {0,0.5,1,1.5,2,2.5,3,3.5,4,5,6,7,8,9,10,12,14,16,20,30,40,50};
  Double_t STBins[32] = {160,165,170,175,180,185,190,195,200,210,220,230,240,250,275,300,325,350,400,450,500,600,700,800,1000,1200,1400,1600,1800,2000,2500,3000};
  //TH2F* h_eeMetSigVsST = new TH2F("eeMetSigVsST",";ST;Met Significance",31,STBins,21,MetSigBins);  

  //---------------ef plots------------------------------
  TH2F* h_efMetSigVsST = new TH2F("efMetSigVsST",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH1F* h_efMet = new TH1F("efMet","Missing transverse energy in events with a fake and an electrons;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* h_efMet_HighR9 = new TH1F("efMet_HighR9","Missing transverse energy in events with a fake and an electrons;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);

  //---------------End ef plots-----------------------

  TH2F* h_eeMetSigVsST_1Jets = new TH2F("eeMetSigVsST_1Jets",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_eeMetSigVsST_2Jets = new TH2F("eeMetSigVsST_2Jets",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_eeMetSigVsST_3pJets = new TH2F("eeMetSigVsST_3pJets",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_eeMetSigVsST_3Jets = new TH2F("eeMetSigVsST_3Jets",";ST;Met/sqrt(HT)",500,0.,5000.,500,0.,50.);
  TH2F* h_eeMetSigVsST_4Jets = new TH2F("eeMetSigVsST_4Jets",";ST;Met/sqrt(HT)",500,0.,5000.,500,0.,50.);
  TH2F* h_eeMetSigVsST_5Jets = new TH2F("eeMetSigVsST_5Jets",";ST;Met/sqrt(HT)",500,0.,5000.,500,0.,50.);
  TH2F* h_eeMetSigVsST_6pJets = new TH2F("eeMetSigVsST_6pJets",";ST;Met/sqrt(HT)",500,0.,5000.,500,0.,50.);

  TH2F* h_ggMetSigVsST_1Jets = new TH2F("ggMetSigVsST_1Jets",";ST;Met/sqrt(HT)",500,0.,5000.,500,0.,50.);
  TH2F* h_ggMetSigVsST_2Jets = new TH2F("ggMetSigVsST_2Jets",";ST;Met/sqrt(HT)",500,0.,5000.,500,0.,50.);
  TH2F* h_ggMetSigVsST_3pJets = new TH2F("ggMetSigVsST_3pJets",";ST;Met/sqrt(HT)",500,0.,5000.,500,0.,50.);

  TH2F* h_ffPVsST = new TH2F("ffPVsST",";ST;Met/sqrt(ST)",500,0.,5000.,500,0.,50.);

  TH2F* h_eeMetSigVsMET_1Jets = new TH2F("eeMetSigVsMET_1Jets",";MET;Met Significance",800,0.,800.,500,0.,50.);
  TH2F* h_eeMetSigVsMET_2Jets = new TH2F("eeMetSigVsMET_2Jets",";MET;Met Significance",800,0.,800.,500,0.,50.);
  TH2F* h_eeMetSigVsMET_3pJets = new TH2F("eeMetSigVsMET_3pJets",";MET;Met Significance",800,0.,800.,500,0.,50.);

  TH2F* h_eeMetSigVsHT_1Jets = new TH2F("eeMetSigVsHT_1Jets",";HT;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_eeMetSigVsHT_2Jets = new TH2F("eeMetSigVsHT_2Jets",";HT;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_eeMetSigVsHT_3pJets = new TH2F("eeMetSigVsHT_3pJets",";HT;Met Significance",500,0.,5000.,500,0.,50.);

  TH2F* h_ffMetSigVsHT_1Jets = new TH2F("ffMetSigVsHT_1Jets",";HT;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_ffMetSigVsHT_2Jets = new TH2F("ffMetSigVsHT_2Jets",";HT;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_ffMetSigVsHT_3pJets = new TH2F("ffMetSigVsHT_3pJets",";HT;Met Significance",500,0.,5000.,500,0.,50.);

  /*  TH2F* h_ffMetSigVsST_lowFracSTEM = new TH2F("ffMetSigVsST_lowFracSTEM",";ST;Met Significance",500,0.,5000.,500,0.,50.);
      TH2F* h_ggMetSigVsST_lowFracSTEM = new TH2F("ggMetSigVsST_lowFracSTEM",";ST;Met Significance",500,0.,5000.,500,0.,50.);
      TH2F* h_eeMetSigVsST_lowFracSTEM = new TH2F("eeMetSigVsST_lowFracSTEM",";ST;Met Significance",500,0.,5000.,500,0.,50.);*/

  TH2F* h_ffMetSigVsST_1Jets = new TH2F("ffMetSigVsST_1Jets",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_ffMetSigVsST_2Jets = new TH2F("ffMetSigVsST_2Jets",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_ffMetSigVsST_3Jets = new TH2F("ffMetSigVsST_3Jets",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_ffMetSigVsST_4Jets = new TH2F("ffMetSigVsST_4Jets",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_ffMetSigVsST_5Jets = new TH2F("ffMetSigVsST_5Jets",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_ffMetSigVsST_6pJets = new TH2F("ffMetSigVsST_6pJets",";ST;Met Significance",500,0.,5000.,500,0.,50.);
  TH2F* h_ffMetSigVsST_3pJets = new TH2F("ffMetSigVsST_3pJets",";ST;Met Significance",500,0.,5000.,500,0.,50.);

  TH1F* h_ggEt_diffCalib = new TH1F("ggEt_diffCalib","Difference between raw and calib transverse energy;Diff (GeV);Events", 150,0.,15.);

  TH1F* h_ggSCE_diffRaw = new TH1F("ggSCE_diffRaw","Difference between raw and corrected SC energy;Diff (GeV);Events", 150,0.,15.);
  TH1F* h_eeSCE_diffRaw = new TH1F("eeSCE_diffRaw","Difference between raw and corrected SC energy;Diff (GeV);Events", 150,0.,15.);
  TH1F* h_ffSCE_diffRaw = new TH1F("ffSCE_diffRaw","Difference between raw and corrected SC energy;Diff (GeV);Events", 150,0.,15.);

  TH1F* h_ggSCE_normalRaw = new TH1F("ggSCE_normalRaw","Correct SC energy to raw SC energy;Ratio;Events",200,0,2);
  TH1F* h_eeSCE_normalRaw = new TH1F("eeSCE_normalRaw","Correct SC energy to raw SC energy;Ratio;Events",200,0,2);  
  TH1F* h_ffSCE_normalRaw = new TH1F("ffSCE_normalRaw","Correct SC energy to raw SC energy;Ratio;Events",200,0,2);

  TH1F* h_ggSCE = new TH1F("ggSCE","Corrected SC energy;Energy (GeV);Events", 500,0.,1000.);
  TH1F* h_eeSCE = new TH1F("eeSCE","Corrected SC energy;Energy (GeV);Events", 500,0.,1000.);
  TH1F* h_ffSCE = new TH1F("ffSCE","Corrected SC energy;Energy (GeV);Events", 500,0.,1000.);

  TH1F* h_ggRawSCE = new TH1F("ggRawSCE","Raw SC energy;Energy (GeV);Events", 500,0.,1000.);
  TH1F* h_eeRawSCE = new TH1F("eeRawSCE","Raw SC energy;Energy (GeV);Events", 500,0.,1000.);
  TH1F* h_ffRawSCE = new TH1F("ffRawSCE","Raw SC energy;Energy (GeV);Events", 500,0.,1000.);

  TH1F* h_ggPt = new TH1F("ggPt","gg Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_ggPtLead = new TH1F("ggPtLead","gg Transverse Momentum of leading object;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_ggPtTrail = new TH1F("ggPtTrail","gg Transverse Momentum of trailing object;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_gammafakePt = new TH1F("gammafakePt","gammafake Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_gammafakePtLead = new TH1F("gammafakePtLead","gammafake Transverse Momentum of leading object;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_gammafakePtTrail = new TH1F("gammafakePtTrail","gammafake Transverse Momentum of trailing object;p_{T} (GeV);Events", 200,0.,1000.);

  TH1F* h_egPt = new TH1F("egPt","eg Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_eePt = new TH1F("eePt","ee Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_ffPt = new TH1F("ffPt","ff Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);

  TH1F* h_gammafakeEta = new TH1F("gammafakeEta","gammafake #eta;#eta;Events", 350,-3.5,3.5);
  TH1F* h_gammafakePhi = new TH1F("gammafakePhi","gammafake #phi;#phi;Events", 64,-3.2,3.2);

  TH1F* h_ggEta = new TH1F("ggEta","gg #eta;#eta;Events", 350,-3.5,3.5);
  TH1F* h_ggPhi = new TH1F("ggPhi","gg #phi;#phi;Events", 64,-3.2,3.2);
  TH1F* h_met = new TH1F("met","missing transverse energy;#slash{E}_{T} (GeV);Events",100,0.,200.);

  TH1F* h_ggSumEt = new TH1F("ggSumEt","Scalar sum of all calorimeter energy in gg events;#sigmaE_{T} (GeV);Events",500,0.0,5000.);
  TH1F* h_eeSumEt = new TH1F("eeSumEt","Scalar sum of all calorimeter energy in ee events;#sigmaE_{T} (GeV);Events",500,0.0,5000.);
  TH1F* h_gammafakeSumEt = new TH1F("gammafakeSumEt","Scalar sum of all calorimeter energy in #gammaf events;#sigmaE_{T} (GeV);Events",500,0.0,5000.);
  TH1F* h_ffSumEt = new TH1F("ffSumEt","Scalar sum of all calorimeter energy in ff events;#sigmaE_{T} (GeV);Events",500,0.0,5000.);

  TH2F* h_ffSumEtVsMET = new TH2F("ffSumEtVsMET","SumEt vs MET in ff events;MET (GeV);SumET (GeV)",400,0.,800.,500,0.,5000.);

  TH2F* h_ggSumEtVsST = new TH2F("ggSumEtVsST","SumEt vs ST in gg events;ST (GeV);SumET (GeV)",500,0.,5000.,500,0.,5000.);
  TH2F* h_eeSumEtVsST = new TH2F("eeSumEtVsST","SumEt vs ST in ee events;ST (GeV);SumET (GeV)",500,0.,5000.,500,0.,5000.);
  TH2F* h_ffSumEtVsST = new TH2F("ffSumEtVsST","SumEt vs ST in ff events;ST (GeV);SumET (GeV)",500,0.,5000.,500,0.,5000.);

  TH1F* h_eeMetMinusDiEM = new TH1F("eeMetMinusDiEM","MET minus two EM objects; MET minus di-EM (GeV);Events",numMetBins,MetBins);
  TH1F* h_ggMetMinusDiEM = new TH1F("ggMetMinusDiEM","MET minus two EM objects; MET minus di-EM (GeV);Events",numMetBins,MetBins);
  TH1F* h_ffMetMinusDiEM = new TH1F("ffMetMinusDiEM","MET minus two EM objects; MET minus di-EM (GeV);Events",numMetBins,MetBins);

  TH1F* h_eeMetMinusDiEM_Reweighted = new TH1F("eeMetMinusDiEM_Reweighted","MET minus two EM objects; MET minus di-EM (GeV);Events",numMetBins,MetBins);
  TH1F* h_ffMetMinusDiEM_Reweighted = new TH1F("ffMetMinusDiEM_Reweighted","MET minus two EM objects; MET minus di-EM (GeV);Events",numMetBins,MetBins);

  Double_t HTBins[43] = {0,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,275,300,325,350,400,450,500,550,600,700,800,1000,1200,1400,1600,1800,2000,2500,3000};

  TH1F* h_eeHT = new TH1F("eeHT","HT in ee events;HT;Events",42,HTBins);
  TH1F* h_eeHTReweighted = new TH1F("eeHTReweighted","di-EM pT reweighted HT in ee events;HT;Events",42,HTBins);

  TH1F* h_ggHT = new TH1F("ggHT","HT in gg events;HT;Events",42,HTBins);
  TH1F* h_ffHT = new TH1F("ffHT","HT in ff events;HT;Events",42,HTBins);
  TH1F* h_ffHTReweighted = new TH1F("ffHTReweighted","di-EM pT reweighted HT in ff events;HT;Events",42,HTBins);
  TH1F* h_gammafakeHT = new TH1F("gammafakeHT","HT in gamma fake events;HT;Events",42,HTBins);
  TH1F* h_gammafakeHTReweighted = new TH1F("gammafakeHTReweighted","di-EM pT reweighted HT in gamma fake events;HT;Events",42,HTBins);

  TH1F* h_eeST = new TH1F("eeST","ST in ee events;ST;Events",800,0.,4000.);
  TH1F* h_ffST = new TH1F("ffST","ST in ff events;ST;Events",800,0.,4000.);
  TH1F* h_ffSTReweighted = new TH1F("ffSTReweighted","di-EM pT reweighted ST in ff events;ST;Events",800,0.,4000.);
  TH1F* h_eeSTReweighted = new TH1F("eeSTReweighted","di-EM pT reweighted ST in ee events;ST;Events",800,0.,4000.);

  TH2F* h_ggMETvsMHT = new TH2F("ggMETvsMHT","Missing HT versus MET for gg events;MET (GeV);Missing HT (GeV)",200,0.,1000.,200,0.,1000.);
  TH2F* h_ffMETvsMHT = new TH2F("ffMETvsMHT","Missing HT versus MET for ff events;MET (GeV);Missing HT (GeV)",200,0.,1000.,200,0.,1000.);
  TH2F* h_egMETvsMHT = new TH2F("egMETvsMHT","Missing HT versus MET for eg events;MET (GeV);Missing HT (GeV)",200,0.,1000.,200,0.,1000.);
  TH2F* h_gammafakeMETvsMHT = new TH2F("gammafakeMETvsMHT","Missing HT versus MET for gammafake events;MET (GeV);Missing HT (GeV)",200,0.,1000.,200,0.,1000.);
  TH2F* h_eeMETvsMHT = new TH2F("eeMETvsMHT","Missing HT versus MET for ee events;MET (GeV);Missing HT (GeV)",200,0.,1000.,200,0.,1000.);

  TH2F* h_ggMETvsMHT_CloseJetMET_ST200 = new TH2F("ggMETvsMHT_CloseJetMET_ST200","Missing HT versus MET for gg events where dPhi(MET,jet)<0.3;MET (GeV);Missing HT (GeV)",200,0.,1000.,200,0.,1000.);
  TH2F* h_ggMETvsMHT_ST200 = new TH2F("ggMETvsMHT_ST200","Missing HT versus MET for gg events;MET (GeV);Missing HT (GeV)",200,0.,1000.,200,0.,1000.);


  TH2F* h_ggMETvsMHT_CloseJetMET = new TH2F("ggMETvsMHT_CloseJetMET","Missing HT versus MET for gg events where dPhi(MET,jet)<0.3;MET (GeV);Missing HT (GeV)",200,0.,1000.,200,0.,1000.);
  TH2F* h_ffMETvsMHT_CloseJetMET = new TH2F("ffMETvsMHT_CloseJetMET","Missing HT versus MET for ff events where dPhi(MET,jet)<0.3;MET (GeV);Missing HT (GeV)",200,0.,1000.,200,0.,1000.);
  TH2F* h_egMETvsMHT_CloseJetMET = new TH2F("egMETvsMHT_CloseJetMET","Missing HT versus MET for eg events where dPhi(MET,jet)<0.3;MET (GeV);Missing HT (GeV)",200,0.,1000.,200,0.,1000.);
  TH2F* h_gammafakeMETvsMHT_CloseJetMET = new TH2F("gammafakeMETvsMHT_CloseJetMET","Missing HT versus MET for gammafake events where dPhi(MET,jet)<0.3;MET (GeV);Missing HT (GeV)",200,0.,1000.,200,0.,1000.);
  TH2F* h_eeMETvsMHT_CloseJetMET = new TH2F("eeMETvsMHT_CloseJetMET","Missing HT versus MET for ee events where dPhi(MET,jet)<0.3;MET (GeV);Missing HT (GeV)",200,0.,1000.,200,0.,1000.);

  TH1F* h_realggMet = new TH1F("realggMet","Missing transverse energy in events with two gen-level photons;#slash{E}_{T} (GeV);Events",200,0.,1000.); //to be compared with the eeMet histogram for DY to ee events to see how well the DY events model the actual MET distribution of digamma events. This histogram is only made for MC.
  TH1F* h_realggHT = new TH1F("realggHT","HT in events with two gen-level photons;HT;Events",800,0.,4000.);
  TH1F* h_ggMet = new TH1F("ggMet","Missing transverse energy in gg events;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  h_ggMet->Sumw2(kFALSE);
  h_ggMet->SetBinErrorOption(TH1::kPoisson);
  TH1F* h_ggMetUnweighted = new TH1F("ggMetUnweighted","Missing transverse energy in gg events;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  h_ggMetUnweighted->Sumw2(kFALSE);
  h_ggMetUnweighted->SetBinErrorOption(TH1::kPoisson);
  TH1F* h_ffMetUnweighted = new TH1F("ffMetUnweighted","Missing transverse energy in ff events;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  h_ffMetUnweighted->Sumw2(kFALSE);
  h_ffMetUnweighted->SetBinErrorOption(TH1::kPoisson);
  TH1F* h_ffMet_HighR9 = new TH1F("ffMet_HighR9","Missing transverse energy in ff events;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  h_ffMet_HighR9->Sumw2(kFALSE);
  h_ffMet_HighR9->SetBinErrorOption(TH1::kPoisson);

  TH1F* h_ffMet_MixR9 = new TH1F("ffMet_MixR9","Missing transverse energy in ff events;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  h_ffMet_MixR9->Sumw2(kFALSE);
  h_ffMet_MixR9->SetBinErrorOption(TH1::kPoisson);

  TH1F* h_ffMet = new TH1F("ffMet","Missing transverse energy in ff events;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  h_ffMet->Sumw2(kFALSE);
  h_ffMet->SetBinErrorOption(TH1::kPoisson);

  TH1F* h_ffMetReweighted = new TH1F("ffMetReweighted","Reweighted missing transverse energy in ff events;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* h_ffMetReweighted_NJet = new TH1F("ffMetReweighted_NJEt","Reweighted missing transverse energy in ff events;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);

  TH1F* h_gammafakeMet_HighR9 = new TH1F("gammafakeMet_HighR9","Missing transverse energy in events with one photon and one fake;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);

  TH1F* h_gammafakeMet = new TH1F("gammafakeMet","Missing transverse energy in events with one photon and one fake;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* h_gammafakeMetReweighted = new TH1F("gammafakeMetReweighted","Missing transverse energy in events with one photon and one fake;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* h_fgMet = new TH1F("fgMet","Missing transverse energy in events with one photon and one fake;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* h_gfMet = new TH1F("gfMet","Missing transverse energy in events with one photon and one fake;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);

  TH1F* h_eeMetReweighted = new TH1F("eeMetReweighted","Missing transverse energy in events with two electrons;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* h_eeMetReweighted_NJet = new TH1F("eeMetReweighted_NJet","Missing transverse energy in events with two electrons;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* h_eeMetReweighted_JESUp = new TH1F("eeMetReweighted_JESUp","Missing transverse energy in events with two electrons;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* h_eeMetReweighted_JESDown = new TH1F("eeMetReweighted_JESDown","Missing transverse energy in events with two electrons;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* h_eeMetReweighted_NJet_JESUp = new TH1F("eeMetReweighted_NJet_JESUp","Missing transverse energy in events with two electrons;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* h_eeMetReweighted_NJet_JESDown = new TH1F("eeMetReweighted_NJet_JESDown","Missing transverse energy in events with two electrons;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);

  TH1F* h_egTransverseMass = new TH1F("egTransverseMass","MT(e,MET) in events with one photon and one electron;Transverse Mass (GeV);Events",numMetBins,MetBins);
  TH2F* h_egTransverseMass_vs_MET = new TH2F("egTransverseMass_vs_MET","MT(e,MET) vs MET;Transverse Mass (GeV);MET",numMetBins,MetBins,numMetBins,MetBins);

  TH1F* h_egMet = new TH1F("egMet","Missing transverse energy in events with one photon and one electron;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* h_eeMet = new TH1F("eeMet","Missing transverse energy in events with two electrons;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* h_eeMet_Funky = new TH1F("eeMet_Funky","Missing transverse energy in events with two electrons;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* h_eeMetUnweighted = new TH1F("eeMetUnweighted","Missing transverse energy in events with two electrons;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);

  //TH1F* h_eeMet = new TH1F("eeMet","Missing transverse energy in events with two electrons;#slash{E}_{T} (GeV);Events",500,0.,1000.); 
  h_egMet->Sumw2(kFALSE);
  h_egMet->SetBinErrorOption(TH1::kPoisson);
  h_eeMet->Sumw2(kFALSE);
  h_eeMet->SetBinErrorOption(TH1::kPoisson);
  h_eeMetUnweighted->Sumw2(kFALSE);
  h_eeMetUnweighted->SetBinErrorOption(TH1::kPoisson);

  TH1F* hMet_Test = new TH1F("Met_Test","Missing transverse energy in events with two electrons;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* hMet_Control = new TH1F("Met_Control","Missing transverse energy in events with two electrons;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);

  TEfficiency * egMTCut = new TEfficiency("egMTCut","Efficiency of applying MT < 100 GeV cut to eg sample",numMetBins,MetBins);

  TEfficiency * crossCheckError = new TEfficiency("crossCheckError","Get errors on gg/ff ratio",numMetBins,MetBins);
  crossCheckError->SetConfidenceLevel(0.683);
  crossCheckError->SetStatisticOption(TEfficiency::kMidP);

  const int numDiemptBins = 25;
  float diemptbins[numDiemptBins+1] = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,70,85,100,115, 130, 150,175,200,300,1000};
  Double_t diemptbins2[numDiemptBins+1] = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,70,85,100,115, 130, 150,175,200,300,1000};
  TH1F* h_realggDiEMPt = new TH1F("realggDiEMPt","real gg DiEMPt;DiEMPt (GeV);Events", 400,0.,1000.);h_realggDiEMPt->Sumw2();
  TH1F* h_ggDiEMPt = new TH1F("ggDiEMPt","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  h_ggDiEMPt->Sumw2(kFALSE);
  h_ggDiEMPt->SetBinErrorOption(TH1::kPoisson);
  TH1F* h_eeDiEMPt = new TH1F("eeDiEMPt","ee DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  h_eeDiEMPt->Sumw2(kFALSE);
  h_eeDiEMPt->SetBinErrorOption(TH1::kPoisson);
  TH1F* h_ffDiEMPt = new TH1F("ffDiEMPt","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  h_ffDiEMPt->Sumw2(kFALSE);
  h_ffDiEMPt->SetBinErrorOption(TH1::kPoisson);
  TH1F* h_gfDiEMPt = new TH1F("gfDiEMPt","gf DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);h_gfDiEMPt->Sumw2();
  TH1F* h_egDiEMPt = new TH1F("egDiEMPt","eg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);h_egDiEMPt->Sumw2();

  TH2F* h_eeDiEMPt_vs_MET = new TH2F("eeDiEMPt_vs_MET","DiEMPt vs MET;Di-EM pT (GeV);MET (GeV)",numDiemptBins,diemptbins2,numMetBins,MetBins);
  TH2F* h_ggDiEMPt_vs_MET = new TH2F("ggDiEMPt_vs_MET","DiEMPt vs MET;Di-EM pT (GeV);MET (GeV)",numDiemptBins,diemptbins2,numMetBins,MetBins);

  TH2F* h_eeDiEMPtWeight_vs_MET = new TH2F("eeDiEMPtWeight_vs_MET","DiEMPt weight vs MET;Di-EM pT weight;MET (GeV)",50,0,5,numMetBins,MetBins);
 
  TH2F* h_eeDiEMPt_vs_HT = new TH2F("eeDiEMPt_vs_HT","DiEMPt vs HT;Di-EM pT (GeV);HT (GeV)",numDiemptBins,diemptbins2,42,HTBins);
  TH2F* h_ggDiEMPt_vs_HT = new TH2F("ggDiEMPt_vs_HT","DiEMPt vs HT;Di-EM pT (GeV);HT (GeV)",numDiemptBins,diemptbins2,42,HTBins);
  TH2F* h_ffDiEMPt_vs_HT = new TH2F("ffDiEMPt_vs_HT","DiEMPt vs HT;Di-EM pT (GeV);HT (GeV)",numDiemptBins,diemptbins2,42,HTBins);

  Double_t njetbins[5] = {-0.5,0.5,1.5,2.5,4.5};
  TH2F* h_eeDiEMPt_vs_NJet = new TH2F("eeDiEMPt_vs_NJet","ee DiEMPt vs NJet;DiEMPt (GeV);NJets", numDiemptBins,diemptbins2,4,njetbins);h_eeDiEMPt_vs_NJet->Sumw2();
  TH2F* h_ggDiEMPt_vs_NJet = new TH2F("ggDiEMPt_vs_NJet","gg DiEMPt vs NJet;DiEMPt (GeV);NJets", numDiemptBins,diemptbins2,4,njetbins);h_ggDiEMPt_vs_NJet->Sumw2();
  TH2F* h_ffDiEMPt_vs_NJet = new TH2F("ffDiEMPt_vs_NJet","ff DiEMPt vs NJet;DiEMPt (GeV);NJets", numDiemptBins,diemptbins2,4,njetbins);h_eeDiEMPt_vs_NJet->Sumw2();

  TH1F * h_eeNJet = new TH1F("eeNJet","Number of jets; N_jets;Events",4,njetbins);
  TH1F * h_ggNJet = new TH1F("ggNJet","Number of jets; N_jets;Events",4,njetbins);
  TH1F * h_ffNJet = new TH1F("ffNJet","Number of jets; N_jets;Events",4,njetbins);

  TH1F * h_eeNJet_Reweighted = new TH1F("eeNJet_Reweighted","Number of jets; N_jets;Events",4,njetbins);
  TH1F * h_ffNJet_Reweighted = new TH1F("ffNJet_Reweighted","Number of jets; N_jets;Events",4,njetbins);

  TH1F* h_realPhosSigIetaIeta = new TH1F("realPhosSigIetaIeta","SigmaIetaIeta of Generator Level Photons Passing All Cuts;#sigma_{i#etai#eta};Events", 150,0.,.03);
  TH1F* h_fakePhosSigIetaIeta = new TH1F("fakePhosSigIetaIeta","SigmaIetaIeta of Fake Photons Passing All Cuts;#sigma_{i#etai#eta};Events", 150,0.,.03);

  TH1F* h_nTotPhotons = new TH1F("nTotPhotons","Total # of Photons Before any Cuts",10,0.,10.);
    
  TH1F* h_ggSigIetaIeta = new TH1F("ggSigIetaIeta","gg SigmaIetaIeta;#sigma_{i#etai#eta};Events", 150,0.,.03);
  TH1F* h_egSigIetaIeta = new TH1F("egSigIetaIeta","eg SigmaIetaIeta;#sigma_{i#etai#eta};Events", 150,0.,.03);
  TH1F* h_eeSigIetaIeta = new TH1F("eeSigIetaIeta","ee SigmaIetaIeta;#sigma_{i#etai#eta};Events", 150,0.,.03);
  TH1F* h_ffSigIetaIeta = new TH1F("ffSigIetaIeta","ff SigmaIetaIeta;#sigma_{i#etai#eta};Events", 150,0.,.03);
  TH1F* h_gfSigIetaIeta = new TH1F("gfSigIetaIeta","gf SigmaIetaIeta;#sigma_{i#etai#eta};Events", 150,0.,.03);

  TH1F* h_ffR9 = new TH1F("ffR9","ff R9;R9;Events", 105,0.0,1.05);
  TH1F* h_eeR9 = new TH1F("eeR9","ee R9;R9;Events", 105,0.0,1.05);
  TH1F* h_ggR9 = new TH1F("ggR9","gg R9;R9;Events", 105,0.0,1.05);

  TH2F* h_phoCands_SigIetaIeta_vs_ChHadIso = new TH2F("phoCands_SigIetaIeta_vs_ChHadIso","Sigma Ieta Ieta vs Charged Hadron Isolation of Photon Candidates passing all other cuts;#sigma_{i#etai#eta};Charged Hadron Isolation (GeV)",30,0.,.03,60,0.,30.);
  TH2F* h_ffSigIetaIeta_vs_ChHadIso = new TH2F("ffSigIetaIeta_vs_ChHadIso","Sigma Ieta Ieta vs Charged Hadron Isolation of Fakes in ff Events;#sigma_{i#etai#eta};Charged Hadron Isolation (GeV)",30,0.,.03,120,0.,30.);
  TH2F* h_ffSigIetaIeta_vs_ChHadIso_Met100 = new TH2F("ffSigIetaIeta_vs_ChHadIso_Met100","Sigma Ieta Ieta vs Charged Hadron Isolation of Fakes in ff Events;#sigma_{i#etai#eta};Charged Hadron Isolation (GeV)",40,0.,.04,100,0.,50.);
  //  TH2F* h_ggSigIetaIeta_vs_ChHadIso = new TH2F("ggSigIetaIeta_vs_ChHadIso","Sigma Ieta Ieta vs Charged Hadron Isolation of Photons in gg Events;#sigma_{i#etai#eta};Charged Hadron Isolation (GeV)",30,0.,.03,120,0.,30.);
  TH2F* h_gfSigIetaIeta_vs_ChHadIso = new TH2F("gfSigIetaIeta_vs_ChHadIso","Sigma Ieta Ieta vs Charged Hadron Isolation of Fakes in gf Events;#sigma_{i#etai#eta};Charged Hadron Isolation (GeV)",30,0.,.03,120,0.,30.);
  TH1F* h_ffChHadIso = new TH1F("ffChHadIso","Charged Hadron Isolation of Fakes in ff Events;Charged Hadron Isolation (GeV);Fakes",120,0.,40.);
  TH1F* h_ffChHadIsoMin = new TH1F("ffChHadIsoMin","Charged Hadron Isolation of Fakes in ff Events;Charged Hadron Isolation (GeV);Fakes",120,0.,40.);
  TH1F* h_ffChHadIsoMax = new TH1F("ffChHadIsoMax","Charged Hadron Isolation of Fakes in ff Events;Charged Hadron Isolation (GeV);Fakes",120,0.,40.);

  TH2F* h_ffChHadIsoMin_vs_minDR = new TH2F("ffChHadIsoMin_vs_minDR","Charged Hadron Isolation of Fakes in ff vs min DR;Charged Hadron Iso;min DR(fake, jet)",150,0,15,150,0,3);
  TH2F* h_ffChHadIsoMin_vs_MET = new TH2F("ffChHadIsoMin_vs_MET","Charged Hadron Isolation of Fakes in ff vs min DR;Charged Hadron Iso;MET (GeV)",150,0,15,numMetBins,MetBins);

  TH2F* h_ggMinDR_vs_MET = new TH2F("ggMinDR_vs_MET","minDR(photon,jet) vs MET;min DR(photon, jet);MET",150,0.,3.0,numMetBins,MetBins);
  TH2F* h_ggMinDR_vs_NJet = new TH2F("ggMinDR_vs_NJet","minDR(photon,jet) vs MET;min DR(photon, jet);NJet",150,0,3,15,0,15);
  TH2F* h_ggNJet_vs_MET = new TH2F("ggNJet_vs_MET","NJet vs MET;NJet;MET",15,0.,15.0, numMetBins,MetBins);

  TH2F* h_ffMinDR_vs_MET = new TH2F("ffMinDR_vs_MET","minDR(fake,jet) vs MET;min DR(fake, jet);MET",150,0.,3.0,numMetBins,MetBins);
  TH2F* h_ffMinDR_vs_NJet = new TH2F("ffMinDR_vs_NJet","minDR(fake,jet) vs MET;min DR(fake, jet);NJet",150,0,3,15,0,15);
  TH2F* h_ffNJet_vs_MET = new TH2F("ffNJet_vs_MET","NJet vs MET;NJet;MET",15,0.,15.0, numMetBins,MetBins);

  TH2F* h_ffNVertex_vs_MET = new TH2F("ffNVertex_vs_MET","NVertex vs MET;NVertex;MET",80,0.,80.0, numMetBins,MetBins);
  TH2F* h_ggNVertex_vs_MET = new TH2F("ggNVertex_vs_MET","NVertex vs MET;NVertex;MET",80,0.,80.0, numMetBins,MetBins);

  TH2F* h_ffRho_vs_MET = new TH2F("ffRho_vs_MET","Rho vs MET;Rho;MET",300,0.,60.0, numMetBins,MetBins);
  TH2F* h_ggRho_vs_MET = new TH2F("ggRho_vs_MET","Rho vs MET;Rho;MET",300,0.,60.0, numMetBins,MetBins);


  TH2F* h_ffSecondMinDR_vs_MET = new TH2F("ffSecondMinDR_vs_MET","minDR(fake,jet) vs MET;min DR(fake, jet);MET",150,0.,3.0,numMetBins,MetBins);
  TH2F* h_ffSecondMinDR_vs_NJet = new TH2F("ffSecondMinDR_vs_NJet","minDR(fake,jet) vs MET;min DR(fake, jet);NJet",150,0,3,15,0,15);

  //  TH1F* h_ggChHadIso = new TH1F("ggChHadIso","Charged Hadron Isolation of Photons in gg Events;Charged Hadron Isolation (GeV);Photons",120,0.,30.);
  TH1F* h_gfChHadIso = new TH1F("gfChHadIso","Charged Hadron Isolation of Fakes in gf Events;Charged Hadron Isolation (GeV);Fakes",120,0.,30.);

  TH2F* h_eeInvMass_vs_Met = new TH2F("eeInvMass_vs_Met","ee Invariant Mass vs MET;MET (Gev);Invariant Mass (GeV)",numMetBins,MetBins,28,74.,130.);

  TH1F* h_ggInvarMass = new TH1F("ggInvarMass","gg Invariant Mass;(GeV);Events", 200,0.,1000.);
  TH1F* h_eeInvarMassFullRange = new TH1F("eeInvarMassFullRange","ee Invariant Mass for all InvarMass;(GeV);Events", 2004,0.,1002.);
    
  //------------------StudyB: Investigating events with high P, low ST---------------------
  TH1F* h_eeMet_StudyB = new TH1F("eeMet_StudyB","Missing transverse energy in ee events with high MetSig, low ST;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* h_eeMHT_StudyB = new TH1F("eeMHT_StudyB","Missing HT in ee events;Missing HT (GeV);Events",500,0,1000);

  TH1F* h_eeST_StudyB = new TH1F("eeST_StudyB","ST in ee events;ST;Events",800,0.,4000.);
  TH1F* h_eeHT_StudyB = new TH1F("eeHT_StudyB","HT in ee events;HT;Events",800,0.,4000.);

  TH1F* h_eeInvarMass_StudyB = new TH1F("eeInvarMass_StudyB","ee Invariant Mass for events with high MetSig, low ST;Invariant Mass (GeV);Events", 2004,0.,1002.);
  TH1F* h_nJets_eeStudyB = new TH1F("nJets_eeStudyB","Number of jets in ee events with high MetSig, low ST",15,0,15);
  TH1F* h_nMuons_eeStudyB = new TH1F("nMuons_eeStudyB","Number of jets in ee sample",15,0,15);
  TH1F* h_nPhotons_eeStudyB = new TH1F("nPhotons_eeStudyB","Number of jets in ee sample",15,0,15);
  TH1F* h_nElectrons_eeStudyB = new TH1F("nElectrons_eeStudyB","Number of jets in ee sample",15,0,15);
  TH1F* h_ee_JetMETdPhi_StudyB = new TH1F("ee_JetMETdPhi_StudyB","#Delta#phi between nearest jet and MET in ee events;#Delta#phi;Events",32,0.,3.2);
  TH1F* h_ee_diEMPtMETdPhi_StudyB = new TH1F("ee_diEMPtMETdPhi_StudyB","#Delta#phi between nearest diempt and MET in ee events;#DeltaR;Events",32,0.,3.2);
  TH1F* h_eeDiEMPt_StudyB = new TH1F("eeDiEMPt_StudyB","ee DiEMPt;DiEMPt (GeV);Events", 200,0.,1000.);

  TH1F* h_eeMet_MetSig3ST400 = new TH1F("eeMet_MetSig3ST400","Missing transverse energy in ee events;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* h_ffMet_MetSig3ST400 = new TH1F("ffMet_MetSig3ST400","Missing transverse energy in ee events;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  //-----------------End of study of high P, low ST events---------------------

  //------Trying to understand correlation - Study C --------------------------
  /*  TH1F* h_ff_fracSTEM_ST80100 = new TH1F("ff_fracSTEM_ST80100","Fraction of ST from 2 photon candidates;Fraction;Events",50,0.,1.0);
      TH1F* h_ff_fracSTEM_ST100120 = new TH1F("ff_fracSTEM_ST100120","Fraction of ST from 2 photon candidates;Fraction;Events",50,0.,1.0);
      TH1F* h_ff_fracSTEM_ST120150 = new TH1F("ff_fracSTEM_ST120150","Fraction of ST from 2 photon candidates;Fraction;Events",50,0.,1.0);
      TH1F* h_ff_fracSTEM_ST150200 = new TH1F("ff_fracSTEM_ST150200","Fraction of ST from 2 photon candidates;Fraction;Events",50,0.,1.0);
      TH1F* h_ff_fracSTEM_ST200 = new TH1F("ff_fracSTEM_ST200","Fraction of ST from 2 photon candidates;Fraction;Events",50,0.,1.0);
      TH2F* h_ffMetSigVsfracSTEM = new TH2F("ffMetSigVsfracSTEM","diPT ST fraction vs MetSig;Fraction;Met Significance",50,0.,1.0,500,0.,50.);
      TH2F* h_ffSTVsfracSTEM = new TH2F("ffSTVsfracSTEM","diPT ST fraction vs ST;Fraction;ST",50,0.,1.,500,0.,5000.);

      TH2F* h_ggMetSigVsfracSTEM = new TH2F("ggMetSigVsfracSTEM","diPT ST fraction vs MetSig;Fraction;Met Significance",50,0.,1.0,500,0.,50.);
      TH2F* h_ggSTVsfracSTEM = new TH2F("ggSTVsfracSTEM","diPT ST fraction vs ST;Fraction;ST",50,0.,1.,500,0.,5000.);

      TH1F* h_ee_fracSTEM_ST80100 = new TH1F("ee_fracSTEM_ST80100","Fraction of ST from 2 photon candidates;Fraction;Events",50,0.,1.0);
      TH1F* h_ee_fracSTEM_ST100120 = new TH1F("ee_fracSTEM_ST100120","Fraction of ST from 2 photon candidates;Fraction;Events",50,0.,1.0);
      TH1F* h_ee_fracSTEM_ST120150 = new TH1F("ee_fracSTEM_ST120150","Fraction of ST from 2 photon candidates;Fraction;Events",50,0.,1.0);
      TH1F* h_ee_fracSTEM_ST150200 = new TH1F("ee_fracSTEM_ST150200","Fraction of ST from 2 photon candidates;Fraction;Events",50,0.,1.0);
      TH1F* h_ee_fracSTEM_ST200 = new TH1F("ee_fracSTEM_ST200","Fraction of ST from 2 photon candidates;Fraction;Events",50,0.,1.0);
      TH2F* h_eeMetSigVsfracSTEM = new TH2F("eeMetSigVsfracSTEM","diPT ST fraction vs MetSig;Fraction;Met Significance",50,0.,1.0,500,0.,50.);
      TH2F* h_eeSTVsfracSTEM = new TH2F("eeSTVsfracSTEM","diPT ST fraction vs ST;Fraction;ST",50,0.,1.,500,0.,5000.);*/
  //------------------End of Study C -----------------------------------------

  //-------------------Study of ee events with high Diempt----------------
  TH1F* h_eeInvarMass_NJet2Diempt10 = new TH1F("eeInvarMass_NJet2Diempt10","ee Invariant Mass;(GeV);Events", 2004,0.,1002.);
  TH1F* h_nJets_eeNJet2Diempt10 = new TH1F("nJets_eeNJet2Diempt10","Number of jets in ee sample",15,0,15);
  TH1F* h_nMuons_eeNJet2Diempt10 = new TH1F("nMuons_eeNJet2Diempt10","Number of muons in ee sample",15,0,15);
  TH1F* h_nPhotons_eeNJet2Diempt10 = new TH1F("nPhotons_eeNJet2Diempt10","Number of photons in ee sample",15,0,15);
  TH1F* h_nElectrons_eeNJet2Diempt10 = new TH1F("nElectrons_eeNJet2Diempt10","Number of electrons in ee sample",15,0,15);
  TH1F* h_eeDiEMPt_NJet2Diempt10 = new TH1F("eeDiEMPt_NJet2Diempt10","ee DiEMPt;DiEMPt (GeV);Events", 200,0.,1000.);
  TH1F* h_eeMet_NJet2Diempt10 = new TH1F("eeMet_NJet2Diempt10","Missing transverse energy in ee events with 2 jets and low diempt;#slash{E}_{T} (GeV);Events",numMetBins,MetBins);
  TH1F* h_ee_JetMETdPhi_NJet2Diempt10 = new TH1F("ee_JetMETdPhi_NJet2Diempt10","#Delta#phi between nearest jet and MET in ee events;#Delta#phi;Events",32,0.,3.2);
  TH1F* h_ee_diEMPtMETdPhi = new TH1F("ee_diEMPtMETdPhi","#Delta#phi between nearest diempt and MET in ee events;#DeltaR;Events",32,0.,3.2);
  TH1F* h_ee_diEMPtMETdPhi_NJet2Diempt10 = new TH1F("ee_diEMPtMETdPhi_NJet2Diempt10","#Delta#phi between nearest diempt and MET in ee events;#DeltaR;Events",32,0.,3.2);

  TH1F* h_nMuons_ee = new TH1F("nMuons_ee","Number of muons in ee sample",15,0,15);
  TH1F* h_nPhotons_ee = new TH1F("nPhotons_ee","Number of photons in ee sample",15,0,15);
  TH1F* h_nElectrons_ee = new TH1F("nElectrons_ee","Number of electrons in ee sample",15,0,15);



  //-------------------End Study-----------------------------


  TH1F* h_nHTJets_gg = new TH1F("nHTJets_gg","Number of HT jets in gg sample",15,0,15);
  TH1F* h_nHTJets_gammafake = new TH1F("nHTJets_gammafake","Number of HT jets in gammafake sample",15,0,15);
  TH1F* h_nHTJets_ee = new TH1F("nHTJets_ee","Number of HT jets in ee sample",15,0,15);
  TH1F* h_nHTJets_ff = new TH1F("nHTJets_ff","Number of HT jets in ff sample",15,0,15);
  TH1F* h_nHTJets_realgg = new TH1F("nHTJets_realgg","Number of HT jets in real gg sample",15,0,15);

  TH1F* h_ggHighestJetEnergy = new TH1F("ggHighestJetEnergy","Pt of most energetic jet in gg events;Energy (GeV);Events",250,0,500.);
  TH1F* h_eeHighestJetEnergy = new TH1F("eeHighestJetEnergy","Pt of most energetic jet in ee events;Energy (GeV);Events",250,0,500.);
  TH1F* h_gammafakeHighestJetEnergy = new TH1F("gammafakeHighestJetEnergy","Pt of most energetic jet in gf events;Energy (GeV);Events",250,0,500.);
  TH1F* h_ffHighestJetEnergy = new TH1F("ffHighestJetEnergy","Pt of most energetic jet in ff events;Energy (GeV);Events",250,0,500.);

  TH1F* h_ffJetEnergy = new TH1F("ffJetEnergy","Pt of jets in ff events;Energy (GeV);Events",150,0,300.);
  TH1F* h_ggJetEnergy = new TH1F("ggJetEnergy","Pt of jets in gg events;Energy (GeV);Events",150,0,300.);
  TH1F* h_eeJetEnergy = new TH1F("eeJetEnergy","Pt of jets in ee events;Energy (GeV);Events",150,0,300.);
  TH1F* h_gammafakeJetEnergy = new TH1F("gammafakeJetEnergy","Pt of jets in gammafake events;Energy (GeV);Events",150,0,300.);

  TH1F* h_ggcomposition = new TH1F("ggcomposition","Truth level composition of gg events; ; Events",3,0,3);    

  TH1F* h_ggMHT = new TH1F("ggMHT","Missing HT in gg events;Missing HT (GeV);Events",500,0,1000);
  TH1F* h_ffMHT = new TH1F("ffMHT","Missing HT in ff events;Missing HT (GeV);Events",500,0,1000);
  TH1F* h_gammafakeMHT = new TH1F("gammafakeMHT","Missing HT in gammafake events;Missing HT (GeV);Events",500,0,1000);
  TH1F* h_eeMHT = new TH1F("eeMHT","Missing HT in ee events;Missing HT (GeV);Events",500,0,1000);
  TH1F* h_egMHT = new TH1F("egMHT","Missing HT in eg events;Missing HT (GeV);Events",500,0,1000);

  TH1F* h_ggMHTReweighted = new TH1F("ggMHTReweighted","DiEM-Pt Reweighted Missing HT in gg events;Missing HT (GeV);Events",500,0,1000);
  TH1F* h_ffMHTReweighted = new TH1F("ffMHTReweighted","DiEm-Pt Reweighted Missing HT in ff events;Missing HT (GeV);Events",500,0,1000);
  TH1F* h_gammafakeMHTReweighted = new TH1F("gammafakeMHTReweighted","DiEm-Pt Reweighted Missing HT in gammafake events;Missing HT (GeV);Events",500,0,1000);
  TH1F* h_eeMHTReweighted = new TH1F("eeMHTReweighted","DiEm-Pt Reweighted Missing HT in ee events;Missing HT (GeV);Events",500,0,1000);
  TH1F* h_egMHTReweighted = new TH1F("egMHTReweighted","DiEm-Pt Reweighted Missing HT in eg events;Missing HT (GeV);Events",500,0,1000);

  TH2F* h_ggMet_vs_diPT = new TH2F("ggMetVsdiPT","Scalar sum of photon PT vs MET in gg events;Missing ET (GeV);Diphoton scalar PT sum (GeV)",500,0.,1000.,1000,0.,4000.);
  TH2F* h_ggHT_vs_diPT = new TH2F("ggHTVsdiPT","Scalar sum of photon PT vs HT in gg events;Missing ET (GeV);Diphoton scalar PT sum (GeV)",800,0.,4000.,1000,0.,4000.);

  TH2F* h_eeMet_vs_diPT = new TH2F("eeMetVsdiPT","Scalar sum of photon PT vs MET in ee events;Missing ET (GeV);Diphoton scalar PT sum (GeV)",500,0.,1000.,1000,0.,4000.);
  TH2F* h_eeHT_vs_diPT = new TH2F("eeHTVsdiPT","Scalar sum of photon PT vs HT in ee events;Missing ET (GeV);Diphoton scalar PT sum (GeV)",800,0.,4000.,1000,0.,4000.);

  TH2F* h_ffMet_vs_diPT = new TH2F("ffMetVsdiPT","Scalar sum of photon PT vs MET in ff events;Missing ET (GeV);Diphoton scalar PT sum (GeV)",500,0.,1000.,1000,0.,4000.);
  TH2F* h_ffHT_vs_diPT = new TH2F("ffHTVsdiPT","Scalar sum of photon PT vs HT in ff events;Missing ET (GeV);Diphoton scalar PT sum (GeV)",800,0.,4000.,1000,0.,4000.);

  TH2F* h_ggMet_vs_ST = new TH2F("ggMetvs_ST","ST vs MET in gg events;Missing ST (GeV);ST (GeV)",500,0.,1000.,800,0.,4000.);
  TH2F* h_ggMet_vs_HT = new TH2F("ggMetvs_HT","HT vs MET in gg events;Missing ST (GeV);HT (GeV)",500,0.,1000.,800,0.,4000.);
  TH2F* h_ffMet_vs_HT = new TH2F("ffMetvs_HT","HT vs MET in ff events;MET (GeV);HT (GeV)",500,0.,1000.,800,0.,4000.);
  TH2F* h_gammafakeMet_vs_HT = new TH2F("gammafakeMetvs_HT","ST vs MET in gammafake events;Missing ST (GeV);ST (GeV)",500,0.,1000.,800,0.,4000.);
  TH2F* h_eeMet_vs_HT = new TH2F("eeMetvs_HT","HT vs MET in ee events;MET (GeV);HT (GeV)",500,0.,1000.,800,0.,4000.);
  TH2F* h_egMet_vs_HT = new TH2F("egMetvs_HT","ST vs MET in eg events;Missing ST (GeV);ST (GeV)",500,0.,1000.,800,0.,4000.);

  TH2F* h_ffMet_vs_ST = new TH2F("ffMetvs_ST","ST vs MET in ff events;MET (GeV);ST (GeV)",500,0.,1000.,800,0.,4000.);
  TH2F* h_eeMet_vs_ST = new TH2F("eeMetvs_ST","ST vs MET in ee events;MET (GeV);ST (GeV)",500,0.,1000.,800,0.,4000.);


  TH2F* h_ggMHT_vs_HT = new TH2F("ggMHT_vs_HT","ST vs missing ST in gg events;Missing ST (GeV);ST (GeV)",500,0.,1000.,800,0.,4000.);
  TH2F* h_ffMHT_vs_HT = new TH2F("ffMHT_vs_HT","ST vs missing ST in ff events;Missing ST (GeV);ST (GeV)",500,0.,1000.,800,0.,4000.);
  TH2F* h_gammafakeMHT_vs_HT = new TH2F("gammafakeMHT_vs_HT","ST vs missing ST in gammafake events;Missing ST (GeV);ST (GeV)",500,0.,1000.,800,0.,4000.);
  TH2F* h_eeMHT_vs_HT = new TH2F("eeMHT_vs_HT","ST vs missing ST in ee events;Missing ST (GeV);ST (GeV)",500,0.,1000.,800,0.,4000.);
  TH2F* h_egMHT_vs_HT = new TH2F("egMHT_vs_HT","ST vs missing ST in eg events;Missing ST (GeV);ST (GeV)",500,0.,1000.,800,0.,4000.);

  TH1F* h_eg_JetMETdPhi = new TH1F("eg_JetMETdPhi","#Delta#phi between nearest jet and MET in eg events;#Delta#phi;Events",32,0.,3.2);
  TH1F* h_gg_JetMETdPhi = new TH1F("gg_JetMETdPhi","#Delta#phi between nearest jet and MET in gg events;#Delta#phi;Events",32,0.,3.2);
  TH1F* h_ff_JetMETdPhi = new TH1F("ff_JetMETdPhi","#Delta#phi between nearest jet and MET in ff events;#Delta#phi;Events",32,0.,3.2);
  TH1F* h_ee_JetMETdPhi = new TH1F("ee_JetMETdPhi","#Delta#phi between nearest jet and MET in ee events;#Delta#phi;Events",32,0.,3.2);
  TH1F* h_gammafake_JetMETdPhi = new TH1F("gammafake_JetMETdPhi","#Delta#phi between nearest jet and MET in gf events;#Delta#phi;Events",32,0.,3.2);

  TH2F* h_gg_JetMETdPhi_vs_Met_ST200 = new TH2F("gg_JetMETdPhi_vs_Met_ST200","#Delta#phi between nearest jet and MET in gg events with ST > 200;E_{T}^{miss};#Delta#phi",200,0.,1000.,32,0.,3.2);
  TH2F* h_ff_JetMETdPhi_vs_Met_ST200 = new TH2F("ff_JetMETdPhi_vs_Met_ST200","#Delta#phi between nearest jet and MET in ff events with ST > 200;E_{T}^{miss};#Delta#phi",200,0.,1000.,32,0.,3.2);
  TH2F* h_ee_JetMETdPhi_vs_Met_ST200 = new TH2F("ee_JetMETdPhi_vs_Met_ST200","#Delta#phi between nearest jet and MET in ee events with ST > 200;E_{T}^{miss};#Delta#phi",200,0.,1000.,32,0.,3.2);
  TH2F* h_gammafake_JetMETdPhi_vs_Met_ST200 = new TH2F("gammafake_JetMETdPhi_vs_Met_ST200","#Delta#phi between nearest jet and MET in gf events with ST > 200;E_{T}^{miss};#Delta#phi",200,0.,1000.,32,0.,3.2);
  TH2F* h_eg_JetMETdPhi_vs_Met_ST200 = new TH2F("eg_JetMETdPhi_vs_Met_ST200","#Delta#phi between nearest jet and MET in eg events with ST > 200;E_{T}^{miss};#Delta#phi",200,0.,1000.,32,0.,3.2);

  TH2F* h_gg_JetMETdPhi_vs_Met = new TH2F("gg_JetMETdPhi_vs_Met","#Delta#phi between nearest jet and MET in gg events;E_{T}^{miss};#Delta#phi",200,0.,1000.,32,0.,3.2);
  TH2F* h_ff_JetMETdPhi_vs_Met = new TH2F("ff_JetMETdPhi_vs_Met","#Delta#phi between nearest jet and MET in ff events;E_{T}^{miss};#Delta#phi",200,0.,1000.,32,0.,3.2);
  TH2F* h_ee_JetMETdPhi_vs_Met = new TH2F("ee_JetMETdPhi_vs_Met","#Delta#phi between nearest jet and MET in ee events;E_{T}^{miss};#Delta#phi",200,0.,1000.,32,0.,3.2);
  TH2F* h_gammafake_JetMETdPhi_vs_Met = new TH2F("gammafake_JetMETdPhi_vs_Met","#Delta#phi between nearest jet and MET in gf events;E_{T}^{miss};#Delta#phi",200,0.,1000.,32,0.,3.2);
  TH2F* h_eg_JetMETdPhi_vs_Met = new TH2F("eg_JetMETdPhi_vs_Met","#Delta#phi between nearest jet and MET in eg events;E_{T}^{miss};#Delta#phi",200,0.,1000.,32,0.,3.2);

  TH1F* h_ggMT2 = new TH1F("ggMT2","MT2 spectrum of gg events;MT2 (GeV);Events",400,0.,2000.);
  TH2F* h_ggMT2_vs_ST = new TH2F("ggMT2_vs_ST","MT2 vs ST for gg events;MT2 (Gev); ST (GeV)",400,0.,2000.,400,0.,2000.);
  TH2F* h_ggMT2_vs_MET =new TH2F("ggMT2_vs_MET","MT2 vs MET for gg events;MT2 (Gev); MET (GeV)",400,0.,2000.,400,0.,1000.);
  TH2F* h_ggMT2_vs_MetSig =new TH2F("ggMT2_vs_MetSig","MT2 vs MetSig for gg events;MT2 (Gev); MetSig",400,0.,2000.,300,0.,30.);

  Double_t MT2bins[37] = {1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64,67,70,75,80,85,90,95,100,105,110,120,130,150,200,250};

  TH1F* h_eeMT2 = new TH1F("eeMT2","MT2 spectrum of ee events;MT2 (GeV);Events",36,MT2bins);
  TH2F* h_eeMT2_vs_ST = new TH2F("eeMT2_vs_ST","MT2 vs ST for ee events;MT2 (Gev); ST (GeV)",36,MT2bins,400,0.,2000.); 
  TH2F* h_eeMT2_vs_HT = new TH2F("eeMT2_vs_HT","MT2 vs HT for ee events;MT2 (Gev); HT (GeV)",36,MT2bins,400,0.,2000.);
  TH2F* h_eeMT2_vs_MET =new TH2F("eeMT2_vs_MET","MT2 vs MET for ee events;MT2 (Gev); MET (GeV)",36,MT2bins,400,0.,1000.);
  TH2F* h_eeMT2_vs_MetSig =new TH2F("eeMT2_vs_MetSig","MT2 vs MetSig for ee events;MT2 (Gev); MetSig",36,MT2bins,300,0.,30.);

  /*  TH1F* h_eeMT2 = new TH1F("eeMT2","MT2 spectrum of ee events;MT2 (GeV);Events",400,0.,2000.);
      TH2F* h_eeMT2_vs_ST = new TH2F("eeMT2_vs_ST","MT2 vs ST for ee events;MT2 (Gev); ST (GeV)",500,0.,500.,400,0.,2000.);
      TH2F* h_eeMT2_vs_HT = new TH2F("eeMT2_vs_HT","MT2 vs HT for ee events;MT2 (Gev); HT (GeV)",500,0.,500.,400,0.,2000.);
      TH2F* h_eeMT2_vs_MET =new TH2F("eeMT2_vs_MET","MT2 vs MET for ee events;MT2 (Gev); MET (GeV)",500,0.,500.,400,0.,1000.);
      TH2F* h_eeMT2_vs_MetSig =new TH2F("eeMT2_vs_MetSig","MT2 vs MetSig for ee events;MT2 (Gev); MetSig",500,0.,500.,300,0.,30.);*/

  TH1F* h_ffMT2 = new TH1F("ffMT2","MT2 spectrum of ff events;MT2 (GeV);Events",400,0.,2000.);
  TH2F* h_ffMT2_vs_ST = new TH2F("ffMT2_vs_ST","MT2 vs ST for ff events;MT2 (Gev); ST (GeV)",400,0.,2000.,400,0.,2000.);
  TH2F* h_ffMT2_vs_MET =new TH2F("ffMT2_vs_MET","MT2 vs MET for ff events;MT2 (Gev); MET (GeV)",400,0.,2000.,400,0.,1000.);
  TH2F* h_ffMT2_vs_MetSig =new TH2F("ffMT2_vs_MetSig","MT2 vs MetSig for ff events;MT2 (Gev); MetSig",400,0.,2000.,300,0.,30.);

  TH1F* h_diffJetPTOverlapPT = new TH1F("diffJetPTOverlapPT","Difference between rejected jet PT and overlapping object PT;Delta PT (GeV);Events",250,0.,1000.);
  TH1F* h_percentDiffJetPTOverlapPT = new TH1F("percentDiffJetPTOverlapPT","Percent Difference between rejected jet PT and overlapping object PT;Percent difference;Events",250,0.,250.);
  TH1F* h_JetPT_overlappingPho = new TH1F("JetPT_overlappingPho","Jet PT for jets within dR < 0.7 of a photon;PT (GeV);Events",250,0.,1000.);




  float xbins2[] = {0,2,5,8,11,14,17,20,25,30,35,40,45,50,57,64,71,80,90, 100,110,120,140,300};
  int binnum2 = sizeof(xbins2)/sizeof(Float_t) - 1;

  float nvtxBins = 50;
  float nvtxMax = 100;

  TH1F *h_All_Signal_MET = new TH1F("h_All_Signal_MET","All Signal MET; Events; MET (GeV)", binnum2, xbins2);
  TH1F *h_All_Signal_Nvtx = new TH1F("h_All_Signal_Nvtx","All Signal Nvtx; Events; Number of Vertices", nvtxBins, 0, nvtxMax );
  TH1F *h_DoublePhoton_Signal_MET = new TH1F("h_DoublePhoton_Signal_MET","Di-photon signal MET; Events; MET (GeV)", binnum2, xbins2);
  TH1F *h_DoublePhoton_Signal_MET_NoTriggerCut = new TH1F("h_DoublePhoton_Signal_MET_NoTriggerCut","Di-photon signal MET; Events; MET (GeV)", binnum2, xbins2);
  TH1F *h_DoublePhoton_Signal_MET_AsymPtCut = new TH1F("h_DoublePhoton_Signal_MET_AsymPtCut","Di-photon signal MET; Events; MET (GeV)", binnum2, xbins2);
  TH1F *h_DoublePhoton_Signal_MET_InvMassCut = new TH1F("h_DoublePhoton_Signal_MET_InvMassCut","Di-photon signal MET; Events; MET (GeV)", binnum2, xbins2);
  TH1F *h_DoublePhoton_Signal_MET_InvMassAndAsymPtCut = new TH1F("h_DoublePhoton_Signal_MET_InvMassAndAsymPtCut","Di-photon signal MET; Events; MET (GeV)", binnum2, xbins2);
  TH1F *h_DoublePhoton_Signal_MET_HighInvMassAndAsymPtCut = new TH1F("h_DoublePhoton_Signal_MET_HighInvMassAndAsymPtCut","Di-photon signal MET; Events; MET (GeV)", binnum2, xbins2);
  TH1F *h_DoublePhotonFromThree_Signal_MET = new TH1F("h_DoublePhotonFromThree_Signal_MET","Di-photon signal MET; Events; MET (GeV)", binnum2, xbins2);
  float xbins[] = {100, 110, 120, 140, 300};
  int binnum = sizeof(xbins)/sizeof(Float_t) - 1;
  TH1F *h_DoublePhoton_Signal_METExpected = new TH1F("h_DoublePhoton_Signal_METExpected","Double Photon MET, MET>100; Events; MET (GeV)", binnum, xbins);
  TH1F *h_DoublePhoton_Signal_METJESUpExpected = new TH1F("h_DoublePhoton_Signal_METJESUpExpected","Double Photon MET, MET>100; Events; MET (GeV)", binnum, xbins);
  TH1F *h_DoublePhoton_Signal_METJESDownExpected = new TH1F("h_DoublePhoton_Signal_METJESDownExpected","Double Photon MET, MET>100; Events; MET (GeV)", binnum, xbins);

  TH1F *h_DoublePhoton_Signal_METExpectedUnweighted = new TH1F("h_DoublePhoton_Signal_METExpectedUnweighted","Double Photon MET, MET>100; Events; MET (GeV)", binnum, xbins);
  TH1F *h_DoublePhoton_Signal_METJESUpExpectedUnweighted = new TH1F("h_DoublePhoton_Signal_METJESUpExpectedUnweighted","Double Photon MET, MET>100; Events; MET (GeV)", binnum, xbins);
  TH1F *h_DoublePhoton_Signal_METJESDownExpectedUnweighted = new TH1F("h_DoublePhoton_Signal_METJESDownExpectedUnweighted","Double Photon MET, MET>100; Events; MET (GeV)", binnum, xbins);

  //  if(!M1M2){
  Double_t xbins3[26]={1275,1325,1375, 1425, 1475, 1525, 1575, 1625, 1675, 1725, 1775, 1825, 1875, 1925, 1975, 2025, 2075, 2125, 2175, 2225, 2275, 2325, 2375, 2425, 2475, 2525};
  int xbinnum3 = 25;
  //   }
  //else{
  //   Double_t xbins3[29]={175, 225, 275, 325, 375, 425, 475, 525, 575, 625, 675, 725, 775, 825, 875, 925, 975, 1025, 1075, 1125, 1175, 1225, 1275, 1325, 1375, 1425, 1475, 1525};
  //int xbinnum3 = 28;
    //}  */

  
  /*  Double_t ybins[135] = {0,17.5,37.5,75,125,175,250,
    350,450,550,650,750,850,950,
    1025.001, 1050.001, 1075.001, 
    1125.001, 1150.001, 1175.001, 
    1225.001, 1250.001,
    1262.5, 1275.001, 1282.5, 1300.001, 1312.5, 1325.001, 1332.5, 1350.001,
    1362.5, 1375.001, 1382.5, 1400.001, 1412.5, 1425.001, 1432.5, 1450.001,
    1462.5, 1475.001, 1482.5, 1500.001, 1512.5, 1525.001, 1532.5, 1550.001,
    1562.5, 1575.001, 1582.5, 1600.001, 1612.5, 1625.001, 1632.5, 1650.001,
    1662.5, 1675.001, 1682.5, 1700.001, 1712.5, 1725.001, 1732.5, 1750.001,
    1762.5, 1775.001, 1782.5, 1800.001, 1812.5, 1825.001, 1832.5, 1850.001,
    1862.5, 1875.001, 1882.5, 1900.001, 1912.5, 1925.001, 1932.5, 1950.001,
    1962.5, 1975.001, 1982.5, 2000.001, 2012.5, 2025.001, 2032.5, 2050.001,
    2062.5, 2075.001, 2082.5, 2100.001, 2112.5, 2125.001, 2132.5, 2150.001,
    2162.5, 2175.001, 2182.5, 2200.001, 2212.5, 2225.001, 2232.5, 2250.001,
    2262.5, 2275.001, 2282.5, 2300.001, 2312.5, 2325.001, 2332.5, 2350.001,
    2362.5, 2375.001, 2382.5, 2400.001, 2412.5, 2425.001, 2432.5, 2450.001, 
    2462.5, 2475.001, 2482.5, 2500.001, 2512.5};
    int ybinnum = 134;   */

  //  if(!M1M2){
  Double_t ybins[123] = {0, 17.5, 37.5, 75, 125, 175, 250, 350, 450, 550, 650, 750, 850, 950, 1025.001, 1050.001, 1075.001, 1125.001, 1150.001, 1175.001, 1225.001, 1250.001, 1262.5, 1275.001, 1282.5, 1300.001, 1312.5, 1325.001, 1332.5, 1350.001, 1362.5, 1375.001, 1382.5, 1400.001, 1412.5, 1425.001, 1432.5, 1450.001, 1462.5, 1475.001, 1482.5, 1500.001, 1512.5, 1525.001, 1532.5, 1550.001, 1562.5, 1575.001, 1582.5, 1600.001, 1612.5, 1625.001, 1632.5, 1650.001, 1662.5, 1675.001, 1682.5, 1700.001, 1712.5, 1725.001, 1732.5, 1750.001, 1762.5, 1775.001, 1782.5, 1800.001, 1812.5, 1825.001, 1832.5, 1850.001, 1862.5, 1875.001, 1882.5, 1900.001, 1912.5, 1925.001, 1932.5, 1950.001, 1962.5, 1975.001, 1982.5, 2000.001, 2012.5, 2025.001, 2032.5, 2050.001, 2062.5, 2075.001, 2082.5, 2100.001, 2112.5, 2125.001, 2132.5, 2150.001, 2162.5, 2175.001, 2182.5, 2200.001, 2212.5, 2225.001, 2232.5, 2250.001,2262.5, 2275.001, 2282.5, 2300.001, 2312.5, 2325.001, 2332.5, 2350.001, 2362.5, 2375.001, 2382.5, 2400.001, 2412.5, 2425.001, 2432.5, 2450.001,2462.5, 2475.001, 2482.5, 2500.001, 2512.5};
      int ybinnum = 122;
    /*  }
  else{
    Double_t ybins[29]={175, 225, 275, 325, 375, 425, 475, 525, 575, 625, 675, 725, 775, 825, 875, 925, 975, 1025, 1075, 1125, 1175, 1225, 1275, 1325, 1375, 1425, 1475, 1525};
    int ybinnum = 28;
    }*/
  TH2D *h_GridAllEvents = new TH2D("h_GridAllEvents","All T5gg Events per mass point; Gluino Mass (GeV); Neutralino Mass (GeV)",xbinnum3, xbins3, ybinnum, ybins);

  TH2D *h_GridAllGGEvents_noInvMass = new TH2D("h_GridAllGGEvents_noInvMass","All gg Events per mass point; Gluino Mass (GeV); Neutralino Mass (GeV)",xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_GridAllGGEvents = new TH2D("h_GridAllGGEvents","All gg Events per mass point; Gluino Mass (GeV); Neutralino Mass (GeV)",xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_GridAllGGEvents2 = new TH2D("h_GridAllGGEvents2","All gg Events per mass point; Gluino Mass (GeV); Neutralino Mass (GeV)",xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_Grid_PhoSF = new TH2D("h_Grid_PhoSF","Photon scale factors per mass point; Gluino Mass (GeV); Neutralino Mass (GeV)",xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_Grid_PhoSF_Err = new TH2D("h_Grid_PhoSF_Err","Photon scale factor uncertainty per mass point; Gluino Mass (GeV); Neutralino Mass (GeV)",xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_Grid_PhoSF_Average = new TH2D("h_Grid_PhoSF_Average","Average photon scale factor per mass point; Gluino Mass (GeV); Neutralino Mass (GeV)",xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_Grid_PhoSF_Err_Average = new TH2D("h_Grid_PhoSF_Err_Average","Average photon SF uncertainty per mass point; Gluino Mass (GeV); Neutralino Mass (GeV)",xbinnum3, xbins3, ybinnum, ybins);

  TH2D *h_MET100to115_GridGG_LowNvtx = new TH2D("h_MET100to115_GridGG_LowNvtx","Events for MET >100 and MET < 115; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET115to130_GridGG_LowNvtx = new TH2D("h_MET115to130_GridGG_LowNvtx","Events for 115 < MET < 130 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET130to150_GridGG_LowNvtx = new TH2D("h_MET130to150_GridGG_LowNvtx","Events for 130 < MET < 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150to185_GridGG_LowNvtx = new TH2D("h_MET150to185_GridGG_LowNvtx","Events for MET >150 and MET < 185 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET185to250_GridGG_LowNvtx = new TH2D("h_MET185to250_GridGG_LowNvtx","Events for 185 < MET < 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET250_GridGG_LowNvtx = new TH2D("h_MET250_GridGG_LowNvtx","Events for MET > 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);

  TH2D *h_MET100to115_GridGG_HighNvtx = new TH2D("h_MET100to115_GridGG_HighNvtx","Events for MET >100 and MET < 115; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET115to130_GridGG_HighNvtx = new TH2D("h_MET115to130_GridGG_HighNvtx","Events for 115 < MET < 130 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET130to150_GridGG_HighNvtx = new TH2D("h_MET130to150_GridGG_HighNvtx","Events for 130 < MET < 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150to185_GridGG_HighNvtx = new TH2D("h_MET150to185_GridGG_HighNvtx","Events for MET >150 and MET < 185 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET185to250_GridGG_HighNvtx = new TH2D("h_MET185to250_GridGG_HighNvtx","Events for 185 < MET < 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET250_GridGG_HighNvtx = new TH2D("h_MET250_GridGG_HighNvtx","Events for MET > 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);

  TH2D *h_MET100to115_GridAll_LowNvtx = new TH2D("h_MET100to115_GridAll_LowNvtx","Events for MET >100 and MET < 115; Gluino Mass (GeV); Neutralino Mass (GeV)", 
						 xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET115to130_GridAll_LowNvtx = new TH2D("h_MET115to130_GridAll_LowNvtx","Events for 115 < MET < 130 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", 
						 xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET130to150_GridAll_LowNvtx = new TH2D("h_MET130to150_GridAll_LowNvtx","Events for 130 < MET < 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", 
						 xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150to185_GridAll_LowNvtx = new TH2D("h_MET150to185_GridAll_LowNvtx","Events for MET >150 and MET < 185 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", 
						 xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET185to250_GridAll_LowNvtx = new TH2D("h_MET185to250_GridAll_LowNvtx","Events for 185 < MET < 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", 
						 xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET250_GridAll_LowNvtx = new TH2D("h_MET250_GridAll_LowNvtx","Events for MET > 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", 
					    xbinnum3, xbins3, ybinnum, ybins);

  TH2D *h_MET100to115_GridAll_HighNvtx = new TH2D("h_MET100to115_GridAll_HighNvtx","Events for MET >100 and MET < 115; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET115to130_GridAll_HighNvtx = new TH2D("h_MET115to130_GridAll_HighNvtx","Events for 115 < MET < 130 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET130to150_GridAll_HighNvtx = new TH2D("h_MET130to150_GridAll_HighNvtx","Events for 130 < MET < 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150to185_GridAll_HighNvtx = new TH2D("h_MET150to185_GridAll_HighNvtx","Events for MET >150 and MET < 185 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET185to250_GridAll_HighNvtx = new TH2D("h_MET185to250_GridAll_HighNvtx","Events for 185 < MET < 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET250_GridAll_HighNvtx = new TH2D("h_MET250_GridAll_HighNvtx","Events for MET > 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);

  TH2D *h_MET100to115_GridAll = new TH2D("h_MET100to115_GridAll","Events for MET >100 and MET < 115; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET115to130_GridAll = new TH2D("h_MET115to130_GridAll","Events for 115 < MET < 130 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET130to150_GridAll = new TH2D("h_MET130to150_GridAll","Events for 130 < MET < 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150to185_GridAll = new TH2D("h_MET150to185_GridAll","Events for MET >150 and MET < 185 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET185to250_GridAll = new TH2D("h_MET185to250_GridAll","Events for 185 < MET < 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET250_GridAll = new TH2D("h_MET250_GridAll","Events for MET > 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150_GridAll = new TH2D("h_MET150_GridAll","Events for MET > 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);

  TH2D *h_MET100to115_Grid = new TH2D("h_MET100to115_Grid","Events for MET >100 and MET < 115; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET115to130_Grid = new TH2D("h_MET115to130_Grid","Events for 115 < MET < 130 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET130to150_Grid = new TH2D("h_MET130to150_Grid","Events for 130 < MET < 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150to185_Grid = new TH2D("h_MET150to185_Grid","Events for MET >150 and MET < 185 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET185to250_Grid = new TH2D("h_MET185to250_Grid","Events for 185 < MET < 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET250_Grid = new TH2D("h_MET250_Grid","Events for MET > 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150_Grid = new TH2D("h_MET150_Grid","Events for MET > 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);

  TH2D *h_MET100to115_GridUnweighted = new TH2D("h_MET100to115_GridUnweighted","Events for MET >100 and MET < 115; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET115to130_GridUnweighted = new TH2D("h_MET115to130_GridUnweighted","Events for 115 < MET < 130 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET130to150_GridUnweighted = new TH2D("h_MET130to150_GridUnweighted","Events for 130 < MET < 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150to185_GridUnweighted = new TH2D("h_MET150to185_GridUnweighted","Events for MET >150 and MET < 185 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET185to250_GridUnweighted = new TH2D("h_MET185to250_GridUnweighted","Events for 185 < MET < 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET250_GridUnweighted = new TH2D("h_MET250_GridUnweighted","Events for MET > 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150_GridUnweighted = new TH2D("h_MET150_GridUnweighted","Events for MET > 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);

  TH2D *h_genMET100to115_GridUnweighted = new TH2D("h_genMET100to115_GridUnweighted","Events for genMET >100 and genMET < 115; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum,ybins);
  TH2D *h_genMET115to130_GridUnweighted = new TH2D("h_genMET115to130_GridUnweighted","Events for 115 < genMET < 130 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_genMET130to150_GridUnweighted = new TH2D("h_genMET130to150_GridUnweighted","Events for 130 < genMET < 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_genMET150to185_GridUnweighted = new TH2D("h_genMET150to185_GridUnweighted","Events for genMET >150 and genMET < 185 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_genMET185to250_GridUnweighted = new TH2D("h_genMET185to250_GridUnweighted","Events for 185 < genMET < 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_genMET250_GridUnweighted = new TH2D("h_genMET250_GridUnweighted","Events for genMET > 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);


  TH2D *h_MET100to115_Grid_JESUp = new TH2D("h_MET100to115_Grid_JESUp","Events for MET >100 and MET < 115; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET115to130_Grid_JESUp = new TH2D("h_MET115to130_Grid_JESUp","Events for 115 < MET < 130 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET130to150_Grid_JESUp = new TH2D("h_MET130to150_Grid_JESUp","Events for 130 < MET < 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150to185_Grid_JESUp = new TH2D("h_MET150to185_Grid_JESUp","Events for MET >150 and MET < 185 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum,ybins);
  TH2D *h_MET185to250_Grid_JESUp = new TH2D("h_MET185to250_Grid_JESUp","Events for 185 < MET < 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET250_Grid_JESUp = new TH2D("h_MET250_Grid_JESUp","Events for MET > 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150_Grid_JESUp = new TH2D("h_MET150_Grid_JESUp","Events for MET > 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);

  TH2D *h_MET100to115_Grid_JESDown = new TH2D("h_MET100to115_Grid_JESDown","Events for MET >100 and MET < 115; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum,ybins);
  TH2D *h_MET115to130_Grid_JESDown = new TH2D("h_MET115to130_Grid_JESDown","Events for 115 < MET < 130 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET130to150_Grid_JESDown = new TH2D("h_MET130to150_Grid_JESDown","Events for 130 < MET < 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150to185_Grid_JESDown = new TH2D("h_MET150to185_Grid_JESDown","Events for MET >150 and MET < 185 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET185to250_Grid_JESDown = new TH2D("h_MET185to250_Grid_JESDown","Events for 185 < MET < 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET250_Grid_JESDown = new TH2D("h_MET250_Grid_JESDown","Events for MET > 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150_Grid_JESDown = new TH2D("h_MET150_Grid_JESDown","Events for MET > 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);

  TH2D *h_MET100to115_GridUnweighted_JESUp = new TH2D("h_MET100to115_GridUnweighted_JESUp","Events for MET >100 and MET < 115; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET115to130_GridUnweighted_JESUp = new TH2D("h_MET115to130_GridUnweighted_JESUp","Events for 115 < MET < 130 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET130to150_GridUnweighted_JESUp = new TH2D("h_MET130to150_GridUnweighted_JESUp","Events for 130 < MET < 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150to185_GridUnweighted_JESUp = new TH2D("h_MET150to185_GridUnweighted_JESUp","Events for MET >150 and MET < 185 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET185to250_GridUnweighted_JESUp = new TH2D("h_MET185to250_GridUnweighted_JESUp","Events for 185 < MET < 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET250_GridUnweighted_JESUp = new TH2D("h_MET250_GridUnweighted_JESUp","Events for MET > 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150_GridUnweighted_JESUp = new TH2D("h_MET150_GridUnweighted_JESUp","Events for MET > 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);

  TH2D *h_MET100to115_GridUnweighted_JESDown = new TH2D("h_MET100to115_GridUnweighted_JESDown","Events for MET >100 and MET < 115; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET115to130_GridUnweighted_JESDown = new TH2D("h_MET115to130_GridUnweighted_JESDown","Events for 115 < MET < 130 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3,xbins3, ybinnum, ybins);
  TH2D *h_MET130to150_GridUnweighted_JESDown = new TH2D("h_MET130to150_GridUnweighted_JESDown","Events for 130 < MET < 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3,xbins3, ybinnum, ybins);
  TH2D *h_MET150to185_GridUnweighted_JESDown = new TH2D("h_MET150to185_GridUnweighted_JESDown","Events for MET >150 and MET < 185 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET185to250_GridUnweighted_JESDown = new TH2D("h_MET185to250_GridUnweighted_JESDown","Events for 185 < MET < 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3,xbins3, ybinnum, ybins);
  TH2D *h_MET250_GridUnweighted_JESDown = new TH2D("h_MET250_GridUnweighted_JESDown","Events for MET > 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150_GridUnweighted_JESDown = new TH2D("h_MET150_GridUnweighted_JESDown","Events for MET > 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);

  TH2D *h_MET100to115_GridStatError = new TH2D("h_MET100to115_GridStatError","Events for MET >100 and MET < 115; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET115to130_GridStatError = new TH2D("h_MET115to130_GridStatError","Events for 115 < MET < 130 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET130to150_GridStatError = new TH2D("h_MET130to150_GridStatError","Events for 130 < MET < 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150to185_GridStatError = new TH2D("h_MET150to185_GridStatError","Events for MET >150 and MET < 185 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET185to250_GridStatError = new TH2D("h_MET185to250_GridStatError","Events for 185 < MET < 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET250_GridStatError = new TH2D("h_MET250_GridStatError","Events for MET > 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150_GridStatError = new TH2D("h_MET150_GridStatError","Events for MET > 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);

  TH2D *h_MET100to115_GridJESUpError = new TH2D("h_MET100to115_GridJESUpError","Events for MET >100 and MET < 115; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET115to130_GridJESUpError = new TH2D("h_MET115to130_GridJESUpError","Events for 115 < MET < 130 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET130to150_GridJESUpError = new TH2D("h_MET130to150_GridJESUpError","Events for 130 < MET < 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150to185_GridJESUpError = new TH2D("h_MET150to185_GridJESUpError","Events for MET >150 and MET < 185 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET185to250_GridJESUpError = new TH2D("h_MET185to250_GridJESUpError","Events for 185 < MET < 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET250_GridJESUpError = new TH2D("h_MET250_GridJESUpError","Events for MET > 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150_GridJESUpError = new TH2D("h_MET150_GridJESUpError","Events for MET > 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);

  TH2D *h_MET100to115_GridJESDownError = new TH2D("h_MET100to115_GridJESDownError","Events for MET >100 and MET < 115; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET115to130_GridJESDownError = new TH2D("h_MET115to130_GridJESDownError","Events for 115 < MET < 130 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET130to150_GridJESDownError = new TH2D("h_MET130to150_GridJESDownError","Events for 130 < MET < 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150to185_GridJESDownError = new TH2D("h_MET150to185_GridJESDownError","Events for MET >150 and MET < 185 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET185to250_GridJESDownError = new TH2D("h_MET185to250_GridJESDownError","Events for 185 < MET < 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET250_GridJESDownError = new TH2D("h_MET250_GridJESDownError","Events for MET > 250 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
  TH2D *h_MET150_GridJESDownError = new TH2D("h_MET150_GridJESDownError","Events for MET > 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);

  TH2D *h_ggDiEMPt400_Grid = new TH2D("ggDiEMPt400_Grid","Events for DiEMPt > 400 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3,  xbins3, ybinnum, ybins);
  TH2D *h_ffDiEMPt400_Grid = new TH2D("ffDiEMPt400_Grid","Events for DiEMPt > 400 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3,  xbins3, ybinnum, ybins);

  TH2D *h_ggDiEMPt300_Grid = new TH2D("ggDiEMPt300_Grid","Events for DiEMPt > 300 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3,  xbins3, ybinnum,ybins);
  TH2D *h_ffDiEMPt300_Grid = new TH2D("ffDiEMPt300_Grid","Events for DiEMPt > 300 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3,  xbins3, ybinnum,ybins);
  TH2D *h_ffRWMET150_Grid = new TH2D("ffRWMET150_Grid","Reweighted ff events with MET > 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3,  xbins3, ybinnum,ybins);
  TH2D *h_ffMET150_Grid = new TH2D("ffMET150_Grid","ff events with MET > 150 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3,  xbins3, ybinnum,ybins);

  TH1F* h_ffDiEMPt_2000_1900 = new TH1F("ffDiEMPt_2000_1900","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  TH1F* h_ffDiEMPt_2000_1500 = new TH1F("ffDiEMPt_2000_1500","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  TH1F* h_ffDiEMPt_2000_1000 = new TH1F("ffDiEMPt_2000_1000","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  TH1F* h_ffDiEMPt_2000_500  = new TH1F("ffDiEMPt_2000_500","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  TH1F* h_ffDiEMPt_2000_100  = new TH1F("ffDiEMPt_2000_100","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);

  TH1F* h_ffDiEMPt_1700_1500 = new TH1F("ffDiEMPt_1700_1500","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  TH1F* h_ffDiEMPt_1700_1000 = new TH1F("ffDiEMPt_1700_1000","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  TH1F* h_ffDiEMPt_1700_500  = new TH1F("ffDiEMPt_1700_500","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  TH1F* h_ffDiEMPt_1700_100  = new TH1F("ffDiEMPt_1700_100","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);

  TH1F* h_ffDiEMPt_1600_1500 = new TH1F("ffDiEMPt_1600_1500","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  TH1F* h_ffDiEMPt_1600_1000 = new TH1F("ffDiEMPt_1600_1000","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  TH1F* h_ffDiEMPt_1600_500  = new TH1F("ffDiEMPt_1600_500","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  TH1F* h_ffDiEMPt_1600_100  = new TH1F("ffDiEMPt_1600_100","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);

  TH1F* h_ffDiEMPt_1500_1400 = new TH1F("ffDiEMPt_1500_1400","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  TH1F* h_ffDiEMPt_1500_1000 = new TH1F("ffDiEMPt_1500_1000","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  TH1F* h_ffDiEMPt_1500_500  = new TH1F("ffDiEMPt_1500_500","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  TH1F* h_ffDiEMPt_1500_100  = new TH1F("ffDiEMPt_1500_100","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);

  TH1F* h_ffDiEMPt_1800_1500 = new TH1F("ffDiEMPt_1800_1500","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  TH1F* h_ffDiEMPt_1800_1000 = new TH1F("ffDiEMPt_1800_1000","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  TH1F* h_ffDiEMPt_1800_500  = new TH1F("ffDiEMPt_1800_500","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);
  TH1F* h_ffDiEMPt_1800_100  = new TH1F("ffDiEMPt_1800_100","ff DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);

  TH1F* h_ggDiEMPt_2000_1900 = new TH1F("ggDiEMPt_2000_1900","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 
  TH1F* h_ggDiEMPt_2000_1500 = new TH1F("ggDiEMPt_2000_1500","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 
  TH1F* h_ggDiEMPt_2000_1000 = new TH1F("ggDiEMPt_2000_1000","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 
  TH1F* h_ggDiEMPt_2000_500  = new TH1F("ggDiEMPt_2000_500","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 
  TH1F* h_ggDiEMPt_2000_100  = new TH1F("ggDiEMPt_2000_100","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 

  TH1F* h_ggDiEMPt_1700_1500 = new TH1F("ggDiEMPt_1700_1500","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 
  TH1F* h_ggDiEMPt_1700_1000 = new TH1F("ggDiEMPt_1700_1000","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 
  TH1F* h_ggDiEMPt_1700_500  = new TH1F("ggDiEMPt_1700_500","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 
  TH1F* h_ggDiEMPt_1700_100  = new TH1F("ggDiEMPt_1700_100","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 

  TH1F* h_ggDiEMPt_1600_1500 = new TH1F("ggDiEMPt_1600_1500","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 
  TH1F* h_ggDiEMPt_1600_1000 = new TH1F("ggDiEMPt_1600_1000","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 
  TH1F* h_ggDiEMPt_1600_500  = new TH1F("ggDiEMPt_1600_500","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 
  TH1F* h_ggDiEMPt_1600_100  = new TH1F("ggDiEMPt_1600_100","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 

  TH1F* h_ggDiEMPt_1500_1400 = new TH1F("ggDiEMPt_1500_1400","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 
  TH1F* h_ggDiEMPt_1500_1000 = new TH1F("ggDiEMPt_1500_1000","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 
  TH1F* h_ggDiEMPt_1500_500  = new TH1F("ggDiEMPt_1500_500","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 
  TH1F* h_ggDiEMPt_1500_100  = new TH1F("ggDiEMPt_1500_100","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 

  TH1F* h_ggDiEMPt_1800_1500 = new TH1F("ggDiEMPt_1800_1500","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 
  TH1F* h_ggDiEMPt_1800_1000 = new TH1F("ggDiEMPt_1800_1000","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 
  TH1F* h_ggDiEMPt_1800_500  = new TH1F("ggDiEMPt_1800_500","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins); 
  TH1F* h_ggDiEMPt_1800_100  = new TH1F("ggDiEMPt_1800_100","gg DiEMPt;DiEMPt (GeV);Events", numDiemptBins,diemptbins);

  TH1F* h_gfMet_2000_1900 = new TH1F("gfMet_2000_1900","gf Met;Met (GeV);Events", numMetBins,MetBins); 
  TH1F* h_gfMet_2000_1500 = new TH1F("gfMet_2000_1500","gf Met;Met (GeV);Events", numMetBins,MetBins); 
  TH1F* h_gfMet_2000_1000 = new TH1F("gfMet_2000_1000","gf Met;Met (GeV);Events", numMetBins,MetBins); 
  TH1F* h_gfMet_2000_500  = new TH1F("gfMet_2000_500","gf Met;Met (GeV);Events", numMetBins,MetBins); 
  TH1F* h_gfMet_2000_100  = new TH1F("gfMet_2000_100","gf Met;Met (GeV);Events", numMetBins,MetBins); 

  TH1F* h_gfMet_1700_1500 = new TH1F("gfMet_1700_1500","gf Met;Met (GeV);Events", numMetBins,MetBins); 
  TH1F* h_gfMet_1700_1000 = new TH1F("gfMet_1700_1000","gf Met;Met (GeV);Events", numMetBins,MetBins); 
  TH1F* h_gfMet_1700_500  = new TH1F("gfMet_1700_500","gf Met;Met (GeV);Events", numMetBins,MetBins); 
  TH1F* h_gfMet_1700_100  = new TH1F("gfMet_1700_100","gf Met;Met (GeV);Events", numMetBins,MetBins); 

  TH1F* h_gfMet_1600_1500 = new TH1F("gfMet_1600_1500","gf Met;Met (GeV);Events", numMetBins,MetBins); 
  TH1F* h_gfMet_1600_1000 = new TH1F("gfMet_1600_1000","gf Met;Met (GeV);Events", numMetBins,MetBins); 
  TH1F* h_gfMet_1600_500  = new TH1F("gfMet_1600_500","gf Met;Met (GeV);Events", numMetBins,MetBins); 
  TH1F* h_gfMet_1600_100  = new TH1F("gfMet_1600_100","gf Met;Met (GeV);Events", numMetBins,MetBins); 

  TH1F* h_gfMet_1500_1400 = new TH1F("gfMet_1500_1400","gf Met;Met (GeV);Events", numMetBins,MetBins); 
  TH1F* h_gfMet_1500_1000 = new TH1F("gfMet_1500_1000","gf Met;Met (GeV);Events", numMetBins,MetBins); 
  TH1F* h_gfMet_1500_500  = new TH1F("gfMet_1500_500","gf Met;Met (GeV);Events", numMetBins,MetBins); 
  TH1F* h_gfMet_1500_100  = new TH1F("gfMet_1500_100","gf Met;Met (GeV);Events", numMetBins,MetBins); 

  TH1F* h_gfMet_1800_1500 = new TH1F("gfMet_1800_1500","gf Met;Met (GeV);Events", numMetBins,MetBins); 
  TH1F* h_gfMet_1800_1000 = new TH1F("gfMet_1800_1000","gf Met;Met (GeV);Events", numMetBins,MetBins); 
  TH1F* h_gfMet_1800_500  = new TH1F("gfMet_1800_500","gf Met;Met (GeV);Events", numMetBins,MetBins); 
  TH1F* h_gfMet_1800_100  = new TH1F("gfMet_1800_100","gf Met;Met (GeV);Events", numMetBins,MetBins); 

  TH1F* h_ffDiffVtx_x = new TH1F("ffDiffVtx_x","Difference between assigned vertex and photon vertex;Distance;Events", 75,0.,0.01);
  TH1F* h_ffDiffVtx_y = new TH1F("ffDiffVtx_y","Difference between assigned vertex and photon vertex;Distance;Events", 75,0.,0.01);
  TH1F* h_ffDiffVtx_z = new TH1F("ffDiffVtx_z","Difference between assigned vertex and photon vertex;Distance;Events", 75,0.,0.01);

  TH1F* h_ggDiffVtx_x = new TH1F("ggDiffVtx_x","Difference between assigned vertex and photon vertex;Distance;Events", 75,0.,0.01);
  TH1F* h_ggDiffVtx_y = new TH1F("ggDiffVtx_y","Difference between assigned vertex and photon vertex;Distance;Events", 75,0.,0.01);
  TH1F* h_ggDiffVtx_z = new TH1F("ggDiffVtx_z","Difference between assigned vertex and photon vertex;Distance;Events", 75,0.,0.01);

  TH1F * phoJetDR  = new TH1F("phoJetDR", "Distance to nearest jet for photons; min dR(photon, jets); photons",150,0.,3.);
  TH1F * phoJetDR2  = new TH1F("phoJetDR2", "Distance to second-closest jet for photons; min dR(photon, jets); photons",150,0.,3.);
  TH1F * fakeJetDR = new TH1F("fakeJetDR","Distance to nearest jet for fakes; min dR(fake, jets); fakes",150,0.,3.);
  TH1F * fakeJetDR2 = new TH1F("fakeJetDR2","Distance to second-closest jet for fakes; min dR(fake, jets); fakes",150,0.,3.);
  TH1F * genPhoJetDR = new TH1F("genPhoJetDR","Distance to nearest jet for all gen-level photons; min dR(gen photon, jets); gen photons",150,0.,3.);
  TEfficiency * phoID_PhoJetDR = new TEfficiency("phoID_PhoJetDR","Eff of photon ID wrt distance to nearest jet; min dR(photon, jets); efficiency",150,0.,3.);

  TH1F* h_GenPhoton_ChargedHadronIsolation= new TH1F("GenPhoton_ChargedHadronIsolation","Charged Hadron Isolation of all gen-level photons from neutralinos;Charged Hadron Isolation (GeV);Photons",300,0.,30.);  

  const int numRhoBins = 37;    

  Double_t RhoBins[numRhoBins+1]= {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,32,34,36,38,40,45,60};
  TH1F* h_ggRho = new TH1F("ggRho",";#rho (Energy/Area);Events/GeV",numRhoBins,RhoBins);
  TH1F* h_ffRho = new TH1F("ffRho",";#rho (Energy/Area);Events/GeV",numRhoBins,RhoBins);

  TH1F* h_rho = new TH1F("rho",";#rho (Energy/Area);Events/GeV",300,0.,30.);

  //insert new histogram definitions here

  fout->cd();
        
  int numPassGenPho = 0, numTotalGenPho = 0;
  int temp = 0, temp2 = 0;
  //Define job wide variables
  bool M1M2 = false;

  int nGoodPho = 0;
  int ggTriggerIndex = 14; //bit with trigger information in the ggNtuples
  int eeTriggerIndex = 15;

  int numPassDefault = 0;

  float lumi = 35883;
  double InvMass;
  TVector2* metvec = new TVector2();
  int nGenPhos;     std::vector<int> genPhosIndices;
  double muMass = 0.1057;  //muon mass, to be used to make the Lorentz vectors needed to find the invariant mass of two objects
  int nTwoCands = 0; //number of events with two photon candidates that pass all the cuts
  bool weird;
  //Define cross sections based on gluino mass                                                                                                                                       
  unordered_map<int,float> xsec_map;
  if(T5Wg){
    xsec_map = {{1300, 0.0460525},
		{1350, 0.0340187},        
		{1400, 0.0252977},
		{1450, 0.0188887},
		{1500, 0.0141903},
		{1550, 0.0107027},
		{1600, 0.00810078},
		{1650, 0.00616072},
		{1700, 0.00470323},
		{1750, 0.00359842},
		{1800, 0.00276133},
		{1850, 0.00212345},
		{1900, 0.00163547},
		{1950, 0.0012642},
		{2000, 0.000981077},
		{2050, 0.000761286},
		{2100, 0.000591918},
		{2150, 0.000460941},
		{2200, 0.000359318},
		{2250, 0.00028065},
		{2300, 0.000219049},
		{2350, 0.000171031},
		{2400, 0.000133965},
		{2450, 0.000104886},
		{2500, 0.0000820068}
    };
  }
  else{
    xsec_map = {{1300, 0.0086557},
		{1350, 0.00637816},
		{1400, 0.00472062},
		{1450, 0.00350484},
		{1500, 0.00261222},
		{1550, 0.0019636},
		{1600, 0.00148553},
		{1650, 0.00112118},
		{1700, 0.000846524},
		{1750, 0.000646271},
		{1800, 0.000494604},
		{1850, 0.000378227},
		{1900, 0.000289986},
		{1950, 0.000223058},
		{2000, 0.000171615},
		{2050, 0.000132986},
		{2100, 0.000102747},
		{2150, 0.000079509}
    };
  }




  Long64_t weirdEvents[] = {151291039};
  unsigned int numWeird = 1;

  float PUweight=1.0;
  int nff_Met100 = 0;
    
  //Indices and momentum of the two photon candidates
  int PhoOne =-1, PhoTwo =-1, PhoThree = -1;
  TLorentzVector PhoOneVec, PhoThreeVec;
  TLorentzVector PhoTwoVec;// = &(new TLorentzVector());
  TLorentzVector PhoOneUnCalibVec, PhoTwoUnCalibVec;
  TLorentzVector thisPho,thisEle,thisMuon,thisJet,thisFake;
  std::vector<TLorentzVector> looseMuons;
  std::vector<int> vetoElesIndices;
  std::vector<TLorentzVector> cleanedPhotons, cleanedEles, cleanedFakes;
  std::vector<TLorentzVector> cleanJets;

  typedef std::multimap<float, int, std::greater<float> > mmid;
  //   The photons will be sorted in the maps by pt so we can select the two with the largest Pt.                                                           
  mmid  pho_Cands;
  mmid  uncleanedPhotons, uncleanedFakes, uncleanedEles;
  mmid  halfcleanedPhotons, halfcleanedFakes, halfcleanedEles;

  float dPhi;
  //counters for the number of events
  int  ngg=0, neg=0, nff=0, nee =0, nee_funky =0;
  int nggMET0to50 = 0, nggMET50to100 =0, nggInvMass110 = 0;
  int nPhoCandsTotal=0,nFakeCandsTotal=0;
  int nFFIso = 0, nFFShape = 0;
  //counters for things I want to keep track of
  int numGammaFake =0;

  int ggMet100to150 = 0, ggMet150to250 =0, ggMet250 =0;
  int numEG_100to150 = 0, numEG_150to250 = 0, numEG_250to350 = 0;
  int numEGSignal[9] = {0,0,0,0,0,0};

  fout->cd();
  //Calculate the PU weights for MC
  if(printLevel>0)cout<<"Read in data PU histogram" << endl;
  TFile *f_dataPU = TFile::Open("MyDataPileupHistogram.root","READ");
  f_dataPU->cd();
  TH1F* h_dataPU = (TH1F *)f_dataPU->Get("pileup");
  h_dataPU->SetName("dataPU");
  h_dataPU->Draw();
  h_dataPU->Scale(1.0 / h_dataPU->Integral() );

  if(printLevel>0)cout<<"Read in MC PU histogram" << endl;
  TFile *f_mcPU = TFile::Open(puFileName,"READ");
  f_mcPU->cd();

  TH1F* h_mcPU = (TH1F *)f_mcPU->Get("puTrue");
  h_mcPU->SetName("mcPU");
  h_mcPU->Draw();
  h_mcPU->Scale(1.0 / h_mcPU->Integral() ) ;
  fout->cd();

  puweights = (TH1F*) f_dataPU->Get("dataPU");
  puweights->Divide( (TH1F*) f_mcPU->Get("mcPU") );

  //Read in the weights histograms located in reweights.root
  if(printLevel>0)cout<<"Read in weights histogram" << endl;
  TFile *f_weights = TFile::Open("reweights.root","READ");
  if(printLevel>0)cout<<"Weights histogram read in" << endl;
  f_weights->cd();
    
  TString eeWeightString = "ee_reweights";
  TH1F* eeweights1 = (TH1F *) f_weights->Get(eeWeightString);
  eeweights1->SetName("eeweights1");
  eeweights1->Draw();
    
  TH1F* ffweights1= (TH1F *) f_weights->Get("ff_reweights");
  ffweights1->SetName("ffweights1");
  ffweights1->Draw();
  
  TH1F* gfweights1 = (TH1F *) f_weights->Get("ff_reweights");
  gfweights1->SetName("gfweights1");
  gfweights1->Draw();
  
  fout->cd();
  ffweights = (TH1F*)f_weights->FindObject("ffweights1");
  ffweights->SetName("ffweights");

  gfweights = (TH1F*)f_weights->FindObject("gfweights1");
  gfweights->SetName("gfweights");
    
  eeweights=(TH1F*)f_weights->FindObject("eeweights1");
  eeweights->SetName("eeweights");
    
    
  if(printLevel>0)cout<<"Line: "<<__LINE__<<endl;
  fout->cd();
  if(printLevel>0)cout<<"Line: "<<__LINE__<<endl;
    
  fout->cd();


  //Get 2D reweights for njet vs diempt systematic
  TFile* f_weights_2D = TFile::Open("NJet_reweights.root","READ");
  f_weights_2D->cd();
  TH1F* ratioee= (TH1F *) f_weights_2D->Get("h_ratio");
  ratioee->SetName("ratioee");
  ratioee->Draw();

  TH1F* ratioff= (TH1F *) f_weights_2D->Get("h_ratio");
  ratioff->SetName("ratioff");
  ratioff->Draw();

  fout->cd();
  TH2F *h2D_weights_ee = (TH2F*) f_weights_2D->FindObject("ratioee");
  h2D_weights_ee->SetName("h2D_weights_ee");

  TH2F *h2D_weights_ff = (TH2F*) f_weights_2D->FindObject("ratioff");
  h2D_weights_ff->SetName("h2D_weights_ff");

  cout<<"Past all of the weight histrogram definitions!"<<endl;

  fout->cd();

  //Generate histograms to be used for diempt systematic uncertainty
  if(printLevel>0) cout<<"Create Toys directory" << endl;
  TDirectory *toys = fout->mkdir("Toys");
  toys->cd();
  if(printLevel>0) cout<<"Toys directory created" << endl;
    
  TH1F *RatioGausEE[1000];
  TH1F *RatioGausFF[1000];
  TH1F *ReweightedEE[1000];
  TH1F *ReweightedFF[1000];
  TH2F *RatioGaus2DEE[1000];
  TH1F *ReweightedEE2D[1000];    

  char *histnameRatio    = new char[27];
  char *histtitleRatio   = new char[50];
    
  TRandom3 r;
  TRandom3 r2;
  TRandom3 r3;    

  if(doToys){
    //Create 1000 toys
    for(int k=0; k<1000; k++){
      if(k%50 == 0 && printLevel >0)printf("processed diEMPt and njets: %d \n", k);
      sprintf(histnameRatio, "DiEMPtee%d",k+1);
      sprintf(histtitleRatio,"random DiEMPt ratio for ee %d",k+1);
      RatioGausEE[k] = (TH1F*)eeweights->Clone(histnameRatio);
      RatioGausEE[k]->Reset();
      RatioGausEE[k]->SetTitle(histtitleRatio);

      sprintf(histnameRatio, "DiEMPtff%d",k+1);
      sprintf(histtitleRatio,"random DiEMPt ratio for ff %d",k+1);
      RatioGausFF[k] = (TH1F*)ffweights->Clone(histnameRatio);
      RatioGausFF[k]-> Reset();
      RatioGausFF[k]-> SetTitle(histtitleRatio);
       
      sprintf(histnameRatio, "eeMetReweighted%d",k+1);
      sprintf(histtitleRatio,"randomly reweighted MET for ee %d",k+1);
      ReweightedEE[k] = (TH1F*)h_eeMetReweighted->Clone(histnameRatio);
      ReweightedEE[k]-> Reset();
      ReweightedEE[k]-> SetTitle(histtitleRatio);
      
      sprintf(histnameRatio, "ffMetReweighted%d",k+1);
      sprintf(histtitleRatio,"randomly reweighted MET for ff %d",k+1);
      ReweightedFF[k] = (TH1F*)h_ffMetReweighted->Clone(histnameRatio);
      ReweightedFF[k]-> Reset();
      ReweightedFF[k]-> SetTitle(histtitleRatio);

      for(int j=0; j<eeweights->GetXaxis()->GetNbins()+1; j++){
	double u = r.Gaus(eeweights->GetBinContent(j+1), eeweights->GetBinError(j+1) );
	RatioGausEE[k]->SetBinContent(j+1,u);
	
	double v = r2.Gaus(ffweights->GetBinContent(j+1), ffweights->GetBinError(j+1));
	RatioGausFF[k]->SetBinContent(j+1,v);
      }

      sprintf(histnameRatio, "2Dreweightsee%d",k+1);
      sprintf(histtitleRatio,"random 2D ratio for ee %d",k+1);
      RatioGaus2DEE[k] = (TH2F*) h2D_weights_ee->Clone(histnameRatio);
      RatioGaus2DEE[k]->Reset();
      RatioGaus2DEE[k]->SetTitle(histtitleRatio);

      sprintf(histnameRatio, "eeMetReweighted_NJet%d",k+1);
      sprintf(histtitleRatio,"randomly reweighted MET for 2D ee %d",k+1);
      ReweightedEE2D[k] = (TH1F*)h_eeMetReweighted->Clone(histnameRatio);
      ReweightedEE2D[k]-> Reset();
      ReweightedEE2D[k]-> SetTitle(histtitleRatio);

      for(int i=0; i< h2D_weights_ee->GetXaxis()->GetNbins()+1; i++){
	for(int j =0;j < h2D_weights_ee->GetYaxis()->GetNbins()+1; j++){
	  double w = r3.Gaus(h2D_weights_ee->GetBinContent(i+1,j+1), h2D_weights_ee->GetBinError(i+1,j+1) );
	  RatioGaus2DEE[k]->SetBinContent(i+1,j+1,w);
	}
      }
    }
  }

  // to check duplicated events
  std::map<int, std::set<int> > allEvents;
  fout->cd();

  // start event looping
  //  Long64_t nbytes = 0, 
  Long64_t nb = 0;
  for (Long64_t jentry=0; jentry < processNEvents; jentry++) {
    

    if(jentry%10000==0){
      cout<<"Processing event number " << jentry << " of " << processNEvents <<" events." <<endl;
    }
    if(printLevel>0)cout<<"-------------\n Processing event number " << jentry << " of " << processNEvents <<" events." <<endl;


    nb = event->GetEntry(jentry);   

    /*    try {nbytes += nb; }
	  catch(std::bad_alloc&){
	  cout <<"nbytes is the problem!" << endl;
	  break;
	  }
	  catch(...){
	  cout <<"nbytes is the problem! (But got the exception wrong)" << endl;
	  break;
	  }*/

    isGG = false;
    isFF = false;
    isFF_anyR9 = false;
    HT  = 0;
    Rho = 0;
    NVertex = 0;
    diEMpt = 0;
    NJets = 0;
    trailPt = 0;
    leadPt = 0;
    trailR9 = 0;
    leadR9 = 0;

    eventNum = event->event_;

    met = event->pfMET_;
    if(met > 350) met =300; //hack to deal with overflow MET events 
      

    if(!event->isData_){
      //Only consider the one quarter of events that are T5gg (vs T5Wg) events 
      int nGenPhotons = 0;
      for(int Part_it = 0; Part_it < event->nMC_; Part_it++){
	float partPdgID = event->mcPID->at(Part_it);
	if(fabs(partPdgID)==22){
	  if(fabs(event->mcMomPID->at(Part_it)) == 1000023){ //10023 = neutralino  
	    nGenPhotons++;
	    //	    if(event->mcPt->at(Part_it) < 1) cout << "Huh. " << event->mcPt->at(Part_it) << endl;
	  }
	}
      }
      if (nGenPhotons < 2 && doSignal) continue;
    }

    int gmass = 0;
    int nmass = 0;
    float weight = 0;
    if(doSignal){
      TString tag = *(event->EventTag_);
      if (tag == ""){
	cout << "TAG is empty!" << endl;
	break;
      }

      TString gmass_str = "";
      TString nmass_str = "";

      string s = string(tag);
      if (s.find("GGM") != std::string::npos){
	M1M2 = true;
        try{
	  int start  = s.find("_",0);
	  int mid    = s.find("_",start+1);
	  int length = (mid - start - 2);

	  gmass_str = tag(start+2,length);
          nmass_str = tag(mid+2,tag.Sizeof() -1 );
        }
        catch(...){
          cout<<"Uh oh! Setting gmass and nmass failed." << endl;
          break;
        }

      }

      else{
	try{
	  gmass_str = tag(5,4);
	  nmass_str = tag(10,tag.Sizeof() -1 );
	}
	catch(...){
	  cout<<"Uh oh! Setting gmass and nmass failed" << endl;
	  break;
	}
      }
      gmass = gmass_str.Atoi();
      nmass = nmass_str.Atoi();

      h_signalMass->Fill(gmass,nmass);

      if(gmass < 1300 && !M1M2) continue;

      if(gmass % 50 != 0 && !M1M2){
	cout << "gmass no good: " << gmass << endl;
	continue;
      }
      //      if(nmass % 100 != 0) {
      //	continue;
      // }

      try{
	weight = xsec_map[gmass];
      }
      catch(...){
	cout<<"Uh oh" << endl;
	break;
      }
      if(weight ==0){ 
	cout << "KILLED BECAUSE OF XSEC MAP! " << gmass << " " << nmass << endl;
	break;
      }
      h_GridAllEvents->Fill( gmass, nmass);
     
      h_All_Signal_MET->Fill(met, weight);

      h_All_Signal_Nvtx->Fill(event->nVtx_);

      //Bin 1 
      if(met >= 100 && met < 115){
	h_MET100to115_GridAll->Fill(gmass, nmass);
	if(event->nVtx_ < 20){
	  h_MET100to115_GridAll_LowNvtx->Fill(gmass, nmass);
	}
	else{
          h_MET100to115_GridAll_HighNvtx->Fill(gmass, nmass);
	}
      }      
      //Bin 2
      else if(met >= 115 && met < 130){
	h_MET115to130_GridAll->Fill(gmass, nmass);
	if(event->nVtx_< 20){
          h_MET115to130_GridAll_LowNvtx->Fill(gmass, nmass);
        }
        else{
          h_MET115to130_GridAll_HighNvtx->Fill(gmass, nmass);
	}
      }
      //Bin 3
      else if(met >= 130 && met < 150){ 
	h_MET130to150_GridAll->Fill(gmass, nmass);
	if(event->nVtx_< 20){
          h_MET130to150_GridAll_LowNvtx->Fill(gmass, nmass);
        }
        else{
          h_MET130to150_GridAll_HighNvtx->Fill(gmass, nmass);
	}
      }
      //Bin 4
      else if(met >= 150 && met < 185){ 
	h_MET150to185_GridAll->Fill(gmass, nmass);
	if(event->nVtx_< 20){
          h_MET150to185_GridAll_LowNvtx->Fill(gmass, nmass);
        }
        else{
          h_MET150to185_GridAll_HighNvtx->Fill(gmass, nmass);
	}
      }
      //Bin 5
      else if(met >= 185 && met < 250){ 
	h_MET185to250_GridAll->Fill(gmass, nmass);
	if(event->nVtx_< 20){
          h_MET185to250_GridAll_LowNvtx->Fill(gmass, nmass);
        }
        else{
          h_MET185to250_GridAll_HighNvtx->Fill(gmass, nmass);
	}
      }	 
      //Bin 6
      else if(met >= 250){              
	h_MET250_GridAll->Fill(gmass, nmass);
	if(event->nVtx_< 20){
          h_MET250_GridAll_LowNvtx->Fill(gmass, nmass);
        }
        else{
          h_MET250_GridAll_HighNvtx->Fill(gmass, nmass);
	}
      }	
    } //End of if(doSignal)
                        

    weird = false;
    for(unsigned int i = 0; i < numWeird; i++){
      if( event->event_== weirdEvents[i]){
        weird = true;
        cout << "================================================" << endl;
	cout << "Weird event " << jentry << " of " << processNEvents <<" events." <<endl;
        cout << "Event to check! Event number "<< event->event_ << " and run number "<<event->run_ << " and lumi " << event->lumis_ << endl;
        cout << "Trigger bits in  HLTPho = " << ( event->HLTPho_>>14&1)<< ", " <<( event->HLTPho_>>15&1)<< ", " <<( event->HLTPho_>>17&1) << endl;
	cout << "Met filters = " << event->metFilters_ << endl;
	cout << "Number of vertices = " << event->nVtx_ << endl;
        cout << "Met = " << event->pfMET_ << endl;
      }
    }


    if(!event->isData_){
      int pu_it = 0; //index of pu 
      for(vector<int>::const_iterator bx = event->puBX_->begin(); bx != event->puBX_->end(); bx++){
	if( *bx == 0){
	  h_puTrue->Fill(event->puTrue_->at(pu_it));
	  if(!doSignal){
	    int bin = puweights->GetXaxis()->FindBin( event->puTrue_->at(pu_it) );
	    PUweight = puweights->GetBinContent(bin);
	  }
	  continue;
	}
	pu_it++;
      }
      //PUweight *= event->genWeight_;
      
      if(printLevel>0) cout << "About to calculate top pT scale factor" << endl;

      if(doTTBar){
	//calculate top pT reweighting factor. 
	float topPTweight = 1.0;
	float SF1 = 1.0, SF2 = 1.0;
	for(int Part_it = 0; Part_it < event->nMC_; Part_it++){
	  float partPdgID = event->mcPID->at(Part_it);
	  if(partPdgID == 6)       SF1 = TMath::Exp(0.0615 - 0.0005*event->mcPt->at(Part_it)) ;
	  else if(partPdgID == -6) SF2 = TMath::Exp(0.0615 - 0.0005*event->mcPt->at(Part_it)) ;
	}
	if(SF1 == 1.0 || SF2 == 1.0){
	  cout << "Could not find a top and an anti-top! Skipping this event for now" << endl;
	  continue;
	}
	topPTweight = TMath::Sqrt(SF1*SF2);
	//	PUweight *= topPTweight;
      }
    }        
    

    if(event->run_>=0){
      nCnt[0]++; // total number of events
            
      if(printLevel>0) cout << "Apply good run (JSON) list." << endl;
      if(!event->isData_) useJSON=false;
      if(useJSON){
	if(weird) cout << "Event rejected because of JSON file? " <<  !isInJson(event->run_,event->lumis_) << " Run number " << event->run_ << " and lumi " << event->lumis_ << endl;
	if(!isInJson(event->run_,event->lumis_)) continue;
      }
    }
        
    nCnt[1]++; // total number of events that pass Json
        
    if(printLevel>0) cout << "Check duplicated events for data only." << endl;
        
    if(event->isData_){
      if(  ! (allEvents[event->run_].insert(event->event_)).second ){
	cout<<"Duplicate Event! Run "<<event->run_<<" Event "<<event->event_<<endl;
	continue;
      }
    }
    
        
    nCnt[2]++;//number of events that pass duplicate check
        
    if(printLevel>0) cout << "Apply trigger selection in the event." << endl;
    bool passHLT = true;
    bool passggTrigger = ((useTrigger && event->isData_ ) ? event->HLTPho_>>ggTriggerIndex&1 : true);
    bool passeeTrigger = ((useTrigger && event->isData_) ? event->HLTPho_>>eeTriggerIndex&1 : true);
    //bool passeeTrigger = ((useTrigger) ? event->HLTPho_>>eeTriggerIndex&1 : true); 
    if(passggTrigger) numPassDefault++;
    if(!version74X){
      passHLT = passggTrigger || passeeTrigger;
    }
    if(weird && !passHLT) cout << "Rejecting event due to trigger requirements" << endl;
    if(!passHLT ) continue;//only accept events which pass our trigger(s)
    nCnt[3]++;// number of events that pass trigger
       
    Rho = event->rho_;
    h_rho->Fill(Rho);
    NVertex=event->nVtx_;
        
    //Require at least one good vertex
    if(weird && (NVertex < 1 || ! (event->isPVGood_) ) ) cout << "Failed vertex requirements! Nvertex = " << NVertex << " and PV is good = " << event->isPVGood_ << endl;
    if(NVertex<1){
      if(printLevel>0){cout<<"No Good Vertex!!!!  Run: "<<event->run_<<"  Event: "<<event->event_<<endl;}
      continue;
    }
    if(! (event->isPVGood_)  ) continue;
    nCnt[4]++;// number of events that pass nVertex>=1
        
    if(weird){
      cout << "Checking MET filters: " ; 
      if(event->isData_){
	cout << (event->metFilters_>>0&1) << " " << (event->metFilters_>>1&1) << " " << 
	  (event->metFilters_>>2&1) << " " << (event->metFilters_>>3&1) << " " << 
	  (event->metFilters_>>4&1) << " " << (event->metFilters_>>5&1) << " " << 
	  (event->metFilters_>>5&1) << " " << (event->metFilters_>>7&1) << " " << 
	  (event->metFilters_>>8&1) << " " << !(event->metFilters_>>9&1) << " " <<
          !(event->metFilters_>>10&1) << " " << (event->metFilters_>>11&1) ;
      }
    }
    //pass MET filters
    if(event->isData_){
      if(event->metFilters_>>0&1) continue;
      if(event->metFilters_>>1&1) continue;
      if(event->metFilters_>>2&1) continue;
      if(event->metFilters_>>3&1) continue;
      if(event->metFilters_>>4&1) continue;
      if(event->metFilters_>>5&1) continue;
      if(event->metFilters_>>6&1) continue;
      if(event->metFilters_>>7&1) continue;
      if(event->metFilters_>>8&1) continue;
      if(!(event->metFilters_>>9&1)) continue;
      if(!(event->metFilters_>>10&1)) continue;
      if(event->metFilters_>>11&1) continue;
    }
    if(weird) cout << "passed." << endl;

    nCnt[5]++;//number of events that pass met filters
  
    metvec->SetMagPhi(met, event->pfMETPhi_);
    h_met->Fill(met);

    if(weird) cout << "Passed met filter, nvertex, and json cuts" << endl;

    //-------------Gen Photons------------------------------    
    //clean up from last round first
    nGenPhos = 0;
    genPhosIndices.clear();
    //now see how many gen photons we have in this event 
    if(!event->isData_){
      for(int Part_it = 0; Part_it < event->nMC_; Part_it++){
	int partPdgID = event->mcPID->at(Part_it);
	//Find all of the real photons in the event, so we can check later to see if the reconstructed photons match any of the gen-level photons
	if(fabs(partPdgID)==22){
	  nGenPhos++;
	}                    
	genPhosIndices.push_back(Part_it);   
      }
    }

    //-------------End of gen Photons------------------------------ 


    //------------Muons--------------------------------------
    //We trust muons more than any other object, so make a list of them first. Any other photons, electrons, or jets that overlap with a muon will be removed
    //    std::vector<TLorentzVector> looseMuons;//we will veto any jets which overlap with these muons when calculating ST  
    looseMuons.clear();
    bool muVeto = false;
    float muPtCut = 25;
    float muEtaCut = 2.4;
    for(int it_Mu = 0; it_Mu < event->nMu_; it_Mu++){
      float pt = event->muPt_->at(it_Mu);
      float eta = event->muEta_->at(it_Mu);

      //      float combIso  = event->muPFChIso_->at(it_Mu) + TMath::Max(0.,event->muPFNeuIso_->at(it_Mu) + event->muPFPhoIso_->at(it_Mu) - 0.5*event->muPFPUIso_->at(it_Mu));
      if(weird){
	cout <<"Muons: pT = " << pt << ", eta = " << eta << ", iso = " << event->muPFMiniIso_->at(it_Mu)  << ", and loose id = " << event->muIDbit_->at(it_Mu) << endl;
      }

      if(pt > muPtCut && 
	 event->muIDbit_->at(it_Mu)>>1&1 && 
	 abs(eta) < muEtaCut && 
	 event->muPFMiniIso_->at(it_Mu) < 0.2  && 
	 abs(event->muD0_->at(it_Mu))   < 0.05 &&
	 abs(event->muDz_->at(it_Mu))  < 0.1) {
	thisMuon = MassLorentzVector(pt, eta, event->muPhi_->at(it_Mu), muMass);
	looseMuons.push_back(thisMuon);
	muVeto = true;
      }
    }

    if(weird){
      cout << "Event number " << event->event_ << endl;
      for(unsigned int imu = 0;imu<looseMuons.size(); imu++){
	cout << "Muon " << imu << ": pT = " << looseMuons.at(imu).Pt() << ", Eta = " << looseMuons.at(imu).Eta() << ", and phi = " << looseMuons.at(imu).Phi() << endl;
      }
      cout <<"About to veto on a muon? " << muVeto << endl;
    }
    if(muVeto) continue;

    h_NumMuons->Fill((int)looseMuons.size());
    //------------End muons------------------------------------
        
    //    cout << "Passed the muons" << endl;
    //---------------Electrons for lepton veto ----------------
    float elePtCut = 25.0;
    float eleEtaBarrel = 1.4442;
    float eleEtaGap = 1.56;
    float eleEtaEndcap = 2.5;
    bool endcapEle = false;
    //    std::vector<int> vetoElesIndices;
    vetoElesIndices.clear();


    for(int it_Ele = 0; it_Ele<event->nEle_; it_Ele++){
      float eleEta = abs(event->eleSCEta_->at(it_Ele) );
      if(weird){
	cout << "Electron index " << it_Ele << " and eta = " << eleEta << endl;
	if(eleEta < eleEtaBarrel){
	  cout << (event->eleCalibPt_->at(it_Ele) < elePtCut)<< endl;
	  cout << (event->eleR9_->at(it_Ele)             < 0.5) << endl;
	  cout << (event->eleSigmaIEtaIEtaFull5x5_->at(it_Ele)  > 0.00998)<< endl;
	  cout << (abs(event->eledEtaAtVtx_->at(it_Ele)) > 0.00311)<< endl;
	  cout << (abs(event->eledPhiAtVtx_->at(it_Ele)) > 0.103)   << endl;
	  cout << (event->eleHoverE_->at(it_Ele)         > 0.253)   << endl;
	  cout << (abs(event->eleEoverPInv_->at(it_Ele)) > 0.134)   << endl;
	  cout << (event->eleMissHits_->at(it_Ele)       > 1)       << endl;
	  cout << (event->eleConvVeto_->at(it_Ele)      == 0)       << endl;
	  cout << (event->elePFMiniIso_->at(it_Ele)      > 0.1)     << endl;
	  cout << event->eledEtaAtVtx_->at(it_Ele) << endl;
	  cout << event->eleEoverPInv_->at(it_Ele) << endl;
	}
      }
      if (event->eleCalibPt_->at(it_Ele) < elePtCut ) continue;
      if( eleEta < eleEtaBarrel){
	if(event->eleR9_->at(it_Ele)             < 0.5)     continue;
	if(event->eleSigmaIEtaIEtaFull5x5_->at(it_Ele)  > 0.00998) continue;
	if(abs(event->eledEtaAtVtx_->at(it_Ele)) > 0.00311) continue;
	if(abs(event->eledPhiAtVtx_->at(it_Ele)) > 0.103)   continue;
	if(event->eleHoverE_->at(it_Ele)         > 0.253)   continue;
	if(abs(event->eleEoverPInv_->at(it_Ele)) > 0.134)   continue;
	if(event->eleMissHits_->at(it_Ele)       > 1)       continue;
	if(event->eleConvVeto_->at(it_Ele)      == 0)       continue;
	if(event->elePFMiniIso_->at(it_Ele)      > 0.1)     continue;
	vetoElesIndices.push_back(it_Ele);
	if(weird) cout << "Adding electron index " << it_Ele << " to veto electrons. " << endl;
      }
      else if(eleEta > eleEtaGap && eleEta < eleEtaEndcap){
	if(event->eleR9_->at(it_Ele)             < 0.8)     continue;
	if(event->eleSigmaIEtaIEtaFull5x5_->at(it_Ele)  > 0.0298)  continue;
	if(abs(event->eledEtaAtVtx_->at(it_Ele)) > 0.00609) continue;
	if(abs(event->eledPhiAtVtx_->at(it_Ele)) > 0.045)   continue;
	if(event->eleHoverE_->at(it_Ele)         > 0.0878)  continue;
	if(abs(event->eleEoverPInv_->at(it_Ele)) > 0.13)    continue;
	if(event->eleMissHits_->at(it_Ele)       > 1)       continue;
	if(event->eleConvVeto_->at(it_Ele)      == 0)       continue;
	if(event->elePFMiniIso_->at(it_Ele)      > 0.1)     continue;
	endcapEle = true;
	if(weird) cout<< "Electron passed all cuts. Index = " << it_Ele << ", eta = " << eleEta << ", and pT = " << event->eleCalibPt_->at(it_Ele) << endl;
      }
    }
    //---------------End veto electrons------------------------
    //    cout << "Passed the electrons" << endl;






    //find photons, sort by energy, define cuts on photons
    //---------------------Photons (and electrons and fakes)------------------------------------------
    //Medium cut based photon id developed during Spring15: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2
    //cut values for Run2
    double maxHoverE = 0.0396;
    double minSietaieta = 0.005;
    double looseSietaieta = 0.01022;
    double maxSietaieta = 0.01022; // 0.03; //0.0106;
    double upperCutSietaieta = 0.0150;
    double looseUpperCutSietaieta = 0.020;
    double looseChIso = 0.441;
    double maxChIso = 0.441;  //this is the maximum rho-corrected isolation.
    double fakeChIso = 1.0;
    double upperCutChIso = 15.0;//40.0;//15.0; // an upper cut on the charged hadron isolation for fakes
    double looseUpperCutChIso = 25.0;
    double maxNeuIso = 2.725; // + 0.0148*pho_pt + 0.000017*pho_pt^2
    double maxPhoIso = 2.571, //+.0047*Photon pt
      maxEta = 2.5;
    double minR9 = 0.5;
        
        
    //find photons, sort by energy, define cuts on photons
    double EtCutValue    = 40.0;
    // maxMETdphi    =2.8;
    double minJetMETdphi = 3.2;
        
    //    typedef std::multimap<float, int, std::greater<float> > mmid;
    //   The photons will be sorted in the maps by pt so we can select the two with the largest Pt.
    pho_Cands.clear();
    uncleanedPhotons.clear();
    uncleanedFakes.clear();
    uncleanedEles.clear();
    halfcleanedPhotons.clear();
    halfcleanedFakes.clear();
    halfcleanedEles.clear();


    int nFakeCands = 0, nEleCands = 0, nPhoCands =0;

    //    std::vector<TLorentzVector> cleanedPhotons, cleanedEles, cleanedFakes;
    cleanedPhotons.clear(), cleanedEles.clear(), cleanedFakes.clear();
    bool funky = false;
    for(int it_Ele = 0; it_Ele<event->nEle_; it_Ele++){  
      if( (std::abs((event->eleSCEta_->at(it_Ele)) - event->eleEta_->at(it_Ele) )> 0.2) && event->eleIDbit_->at(it_Ele)>>2&1){
	//	cout << "NOOOOOOOOO! Electron with large eta difference. Eta = " << event->eleEta_->at(it_Ele) << " but SC Eta = " << event->eleSCEta_->at(it_Ele) << ". Event number " << event->event_ << ". Met = " << met << endl;
	hMet_Test->Fill(met);
	funky = true;
      }
    }
    //For ttbar events, check that we have two prompt electrons
    float numPromptEle =0;
    if(!event->isData_ && doTTBar){
      for(int i = 0; i < event->nMC_; i++) {
	if(std::abs(event->mcPID->at(i)) == 11 && std::abs(event->mcMomPID->at(i)) == 24 && std::abs(event->mcGMomPID->at(i)) == 6){
	  numPromptEle++;
	}
      }
      if(numPromptEle < 2 ) continue;
    }

    //To minimize ttbar contribution, remove any events that have 2 or more b-tagged jets
    //Using the tight cut defined here https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
    float numBJets    = 0;
    for(int i = 0; i < event->nJet_; i++){
      if(event->jetCSV2BJetTags_->at(i) > 0.9535) numBJets++;
    }
    //    if(numBJets >= 2 ) continue;
    
    if(!funky) hMet_Control->Fill(met);
    if(funky && weird){
      cout << "Weird, funky event!! " << endl;
      for(int it_Ele = 0; it_Ele<event->nEle_; it_Ele++){
	cout << "Electron " << it_Ele << ", dEta = " <<std::abs((event->eleSCEta_->at(it_Ele)) - event->eleEta_->at(it_Ele) ) <<  ", eleIDbit = " <<  event->eleIDbit_->at(it_Ele) << endl;
      }
    }
    if(funky) continue;
    //Loop over the photon collection!!!
    nGoodPho = 0;
    int nTotPhotons = event->nPho_;
    h_nTotPhotons->Fill(nTotPhotons);
    //ignore event if there are not at least 2 photon candidates
    if(event->isData_){
      if(nTotPhotons < 2) continue; //Stop analyzing the event if there are not at least 2 pho candidates
    }
  
    bool photonID, fakeID, electronID;
    float chIsoEA=0, neuIsoEA=0, phoIsoEA=0;
  
    for(int it_Pho = 0; it_Pho<nTotPhotons; it_Pho++){
      if(printLevel>0) cout<<"looping over photon collection"<<endl;
      //do things that want all photons with no cuts here: ie fill all histograms and such
      h_FakeCuts->Fill(0);
      double phoEta = std::abs(event->phoSCEta_->at(it_Pho));
    
      if(std::abs((event->phoSCEta_->at(it_Pho)) - event->phoEta_->at(it_Pho))> 0.2) cout << "NOOOOOOOOO! Event with large eta difference. Eta = " << event->phoEta_->at(it_Pho) << " but SC Eta = " << event->phoSCEta_->at(it_Pho) << ". Event number " << event->event_ << endl;

      if(phoEta >= maxEta){ //only consider barrel photons
	continue;
      }
      //bool PhoEtaCut = phoEta < maxEta;
      //--------------Define rho-corrected isolations------  see https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2
            
            
      if(phoEta<1.0){
	chIsoEA  = 0.0360;
	neuIsoEA = 0.0597;
	phoIsoEA = 0.1210;
      }
      else if (phoEta <1.479 && phoEta >=1.0){
	chIsoEA  = 0.0377;
	neuIsoEA = 0.0807;
	phoIsoEA = 0.1107;
      }
      else if (phoEta <2.0 && phoEta >=1.479){
	chIsoEA  = 0.0306;
	neuIsoEA = 0.0629;
	phoIsoEA = 0.0699;
      }
      else if (phoEta <2.2 && phoEta >=2.0){
	chIsoEA  = 0.0283;
	neuIsoEA = 0.0197;
	phoIsoEA = 0.1056;
      }
      else if (phoEta <2.3 && phoEta >=2.2){
	chIsoEA  = 0.0254;
	neuIsoEA = 0.0184;
	phoIsoEA=  0.1457;
      }
      else if (phoEta <2.4 && phoEta >=2.3){
	chIsoEA  = 0.0217;
	neuIsoEA = 0.0284;
	phoIsoEA = 0.1719;
      }
      else if (phoEta >=2.4){
	chIsoEA  = 0.0167;
	neuIsoEA = 0.0591;
	phoIsoEA = 0.1998;
      }
                
      //----------------set up cuts-------------------
      float PhoEt= event->phoCalibEt_->at(it_Pho);
      float PhoHoverE = event->phoHoverE_->at(it_Pho);            
      float PhoSIetaIeta = event->phoSigmaIEtaIEtaFull5x5_->at(it_Pho);
      float PhoR9 = event->phoR9_->at(it_Pho);
      //uncorrected isolations
      float uncorrNeuHadIso = event->phoPFNeuIso_->at(it_Pho);
      float uncorrChHadIso = event->phoPFChIso_->at(it_Pho);
      float uncorrPhoIso = event->phoPFPhoIso_->at(it_Pho);
            
      //These are the pile-up corrected isolations
      float chargedHadronIso= uncorrChHadIso - Rho*chIsoEA > 0. ? uncorrChHadIso - Rho*chIsoEA : 0.00 ;
      float neutralHadronIso= ( uncorrNeuHadIso - Rho*neuIsoEA - 0.0148*PhoEt - 0.000017*PhoEt*PhoEt ) > 0.00 ? uncorrNeuHadIso - Rho*neuIsoEA - 0.0148*PhoEt - 0.000017*PhoEt*PhoEt: 0.00;
      float photonIso = ( uncorrPhoIso - Rho*phoIsoEA - 0.0047*PhoEt ) > 0.00 ? uncorrPhoIso - Rho*phoIsoEA - 0.0047*PhoEt : 0.00;

      float PhoUnCalibEt = event->phoEt_->at(it_Pho);
      float neutralHadronIsoUnCalib= ( uncorrNeuHadIso - Rho*neuIsoEA - 0.0148*PhoUnCalibEt - 0.000017*PhoUnCalibEt*PhoUnCalibEt ) > 0.00 ? uncorrNeuHadIso - Rho*neuIsoEA - 0.0148*PhoUnCalibEt  - 0.000017*PhoUnCalibEt*PhoUnCalibEt: 0.00;
      float photonIsoUnCalib = ( uncorrPhoIso - Rho*phoIsoEA - 0.0047*PhoUnCalibEt ) > 0.00 ? uncorrPhoIso - Rho*phoIsoEA - 0.0047*PhoUnCalibEt : 0.00;



      // Et cuts
      bool EtCut = (PhoEt > EtCutValue);
            
      // H/E
      bool heCut = (PhoHoverE < maxHoverE);
      // sigma_ietaieta
      bool sIetaCut = (PhoSIetaIeta < maxSietaieta);
      bool upperSIetaCut = (PhoSIetaIeta< upperCutSietaieta);
      bool looseUpperSIetaCut = (PhoSIetaIeta< looseUpperCutSietaieta);
      bool looseSIetaCut = (PhoSIetaIeta< looseSietaieta);
      bool minSIetaCut = (PhoSIetaIeta > minSietaieta);
      //R9
      bool R9Cut = (PhoR9 < 1.0 && PhoR9 > minR9);

      //isolation cuts
      bool looseChHadIsoCut = chargedHadronIso < looseChIso;
      bool chHadIsoCut = chargedHadronIso < maxChIso;
      bool fakeChHadIsoCut = chargedHadronIso > fakeChIso;
      bool upperChIsoCut = chargedHadronIso < upperCutChIso; //for fakes
      bool looseUpperChIsoCut = chargedHadronIso < looseUpperCutChIso;
      bool neuHadIsoCut = neutralHadronIso < maxNeuIso;
      bool phoIsoCut = photonIso < maxPhoIso;

      bool pixelCut = !event->phohasPixelSeed_->at(it_Pho);

      if(weird){
	cout << "---------------------------------------------------" << endl;
        cout << "Event number " << event->event_ << endl;
	cout <<"Photon candidate " << it_Pho << " has Pt " << PhoEt<< ", Phi = ,"<<event->phoPhi_->at(it_Pho) << " and Eta " << phoEta <<endl;
	cout << "Pho (non-calib) ET = " << event->phoEt_->at(it_Pho) << endl;
	cout << "Pho iso (non-calib) = " << photonIsoUnCalib << endl;
	cout <<"Neu iso (non-calib) = " << neutralHadronIsoUnCalib <<endl;
	cout << "Neu iso (uncorrected) = " << uncorrNeuHadIso << endl;
	cout << "Rho = " << Rho << endl;
	cout <<"H/E = " << PhoHoverE <<  ", cut = " << heCut<< endl;
	cout <<"Sig ieta ieta = " << PhoSIetaIeta <<  ", cut = " << sIetaCut << endl;
	cout <<"PhoR9 = " << PhoR9 <<  ", cut = " << R9Cut << endl;
	cout <<"Corrected charged hadron isolation = " << chargedHadronIso <<  ", cut = " << chHadIsoCut <<endl;
	cout <<"Corrected photon isolation  = " << photonIso <<  ", cut = " << phoIsoCut <<endl;
	cout <<"Corrected neutral hadron isolation = " << neutralHadronIso << ", cut = " << neuHadIsoCut<< endl;
	cout << "Pixel seed veto = " << event->phohasPixelSeed_->at(it_Pho) << ", cut = " << pixelCut << endl;
	cout << "---------------------------------------------------" << endl;
      }

      if(printLevel>0) cout<<"photon cuts defined"<<endl;
      if(EtCut && heCut && neuHadIsoCut && phoIsoCut) nGoodPho++;

      //photonID = ( EtCut && (event->phoIDbit_->at(it_Pho)>>1&1) && pixelCut && R9Cut); 
      photonID = ( EtCut && heCut && neuHadIsoCut && phoIsoCut && chHadIsoCut && sIetaCut && pixelCut && R9Cut && minSIetaCut);
      electronID = ( EtCut && heCut && neuHadIsoCut && phoIsoCut && chHadIsoCut && sIetaCut && !pixelCut && R9Cut && minSIetaCut);
      // fakeID =   ( EtCut && heCut && neuHadIsoCut && phoIsoCut && R9Cut &&  upperChIsoCut && (fakeChHadIsoCut || !sIetaCut) && upperSIetaCut && pixelCut && !(!chHadIsoCut && !sIetaCut) && minSIetaCut );

      //alternate fake ID - both high and low R9
      //fakeID =   ( EtCut && heCut  && neutralHadronIso < 15.0 && photonIso < 15.0 && PhoR9 > 0.5 /*&& PhoR9 < 0.9 */&&  looseUpperChIsoCut && (!chHadIsoCut || !sIetaCut) && upperSIetaCut && pixelCut /*&& !(fakeChHadIsoCut && !sIetaCut)*/ && minSIetaCut);
      //alternate fake ID
      fakeID =   ( EtCut && heCut  && neutralHadronIso < 15.0 && photonIso < 15.0 && PhoR9 > 0.5 && PhoR9 < 0.9 &&  looseUpperChIsoCut && (!chHadIsoCut || !sIetaCut) && upperSIetaCut && pixelCut /*&& !(fakeChHadIsoCut && !sIetaCut)*/ && minSIetaCut);
      //High R9 version
      //      fakeID =   ( EtCut && heCut  && neutralHadronIso < 15.0 && photonIso < 15.0 && PhoR9 > 0.9 &&  looseUpperChIsoCut && (!chHadIsoCut || !sIetaCut) && upperSIetaCut && pixelCut /*&& !(fakeChHadIsoCut && !sIetaCut)*/ && minSIetaCut);
      //fakeID =   ( EtCut && heCut && neutralHadronIso < 15.0 && photonIso < 15.0 && PhoR9 > 0.5 && PhoR9 < 0.9 &&  looseUpperChIsoCut && (!chHadIsoCut || !sIetaCut) && upperSIetaCut && pixelCut /*&& !(fakeChHadIsoCut && !sIetaCut)*/ && minSIetaCut);
      //fakeID =   ( EtCut && heCut && neuHadIsoCut && phoIsoCut && R9Cut && upperChIsoCut && upperSIetaCut && pixelCut && !sIetaCut );
      //loose fake ID 
      //      fakeID =   ( EtCut && heCut && /* neuHadIsoCut &&  phoIsoCut && */ R9Cut && looseUpperChIsoCut && (fakeChHadIsoCut || !sIetaCut) && looseUpperSIetaCut &&  pixelCut /* && !(!chHadIsoCut && !sIetaCut)*/ );
      //fakeID =   ( EtCut && /*heCut && */ neuHadIsoCut &&  phoIsoCut &&  R9Cut && chHadIsoCut && !sIetaCut && looseUpperSIetaCut &&  pixelCut /* && !(!chHadIsoCut && !sIetaCut)*/ );

      if(weird) cout << "First half " << (( EtCut && heCut  && neutralHadronIso < 15.0 && photonIso < 15.0 && PhoR9 > 0.5 && PhoR9 < 0.9 )) << " second half " << (looseUpperChIsoCut && (!chHadIsoCut || !sIetaCut) && upperSIetaCut && pixelCut && minSIetaCut) << endl;
      if( weird)  cout<< "photonID " << photonID << ", electronID " << electronID <<", fakeID " << fakeID <<endl;


      if(EtCut && heCut && neuHadIsoCut && phoIsoCut && R9Cut){
	h_phoCands_SigIetaIeta_vs_ChHadIso->Fill(PhoSIetaIeta ,chargedHadronIso );
      }
             
      h_FakeCuts->Fill(1);
      if(EtCut){
	h_FakeCuts->Fill(2);
	if(heCut){
	  h_FakeCuts->Fill(3);
	  if(neuHadIsoCut){
	    h_FakeCuts->Fill(4);
	    if(phoIsoCut){
	      h_FakeCuts->Fill(5);
	      if(chHadIsoCut){
		h_FakeCuts->Fill(6);
		if(sIetaCut){
		  h_FakeCuts->Fill(7);
		  if(pixelCut){
		    h_FakeCuts->Fill(8);
		    if(R9Cut){
		      h_FakeCuts->Fill(9);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }

         
      if(photonID){
	nPhoCandsTotal++;
	nPhoCands++;
	uncleanedPhotons.insert(mmid::value_type(PhoEt, it_Pho));
      }
      if(fakeID){
	nFakeCands++;
	nFakeCandsTotal++;
	uncleanedFakes.insert(mmid::value_type(PhoEt, it_Pho));
      }
      if(electronID){
	nEleCands++;
	uncleanedEles.insert(mmid::value_type(PhoEt, it_Pho));
      }
            
      if(printLevel>0) cout<<"End of Photon Loop"<<endl;
    }//end of photon loop
    //    cout<<"end of photon loop"<<endl;
    //    if(nGoodPho < 2) cout << "PHEW! Skimming let this through: " <<event->event_ <<  endl;

    //--------------------------------------------------------------------
    //---------------------Object cleaning--------------------------------
  
    //Now clean eles, photons, and fakes from muons and from each other
    //Remove electrons that overlap within dR < 0.3 of a muon
    for (mmid::const_iterator iele = uncleanedEles.begin(); iele != uncleanedEles.end(); iele++){
      bool overlapMu = false;
      thisEle.SetPtEtaPhiE(iele->first,
			   event->phoEta_->at(iele->second),
			   event->phoPhi_->at(iele->second),
			   event->phoCalibE_->at(iele->second));
      for(unsigned int imu = 0;imu<looseMuons.size(); imu++){
	if(getDR(looseMuons.at(imu),thisEle) < 0.3) overlapMu = true;
      }
      if(!overlapMu){
	halfcleanedEles.insert(mmid::value_type(iele->first, iele->second ) );
      }
      else if(overlapMu && weird) cout << "Electron " << iele->second << " removed due to cleaning!" << endl;
    }
    //Remove photons that overlap within dR < 0.3 of an electron or muon
    for (mmid::const_iterator ipho = uncleanedPhotons.begin(); ipho != uncleanedPhotons.end(); ipho++){
      bool overlap = false;
      thisPho.SetPtEtaPhiE(ipho->first, 
				 event->phoEta_->at(ipho->second), 
				 event->phoPhi_->at(ipho->second),
				 event->phoCalibE_->at(ipho->second));
      for(unsigned int imu = 0;imu<looseMuons.size(); imu++){
        if(getDR(looseMuons.at(imu),thisPho) < 0.3) overlap = true;
      }
      for(mmid::const_iterator  iele = halfcleanedEles.begin(); iele != halfcleanedEles.end(); iele++){
	thisEle.SetPtEtaPhiE(iele->first,
				   event->phoEta_->at(iele->second),
				   event->phoPhi_->at(iele->second),
				   event->phoCalibE_->at(iele->second));
	if(getDR(thisEle, thisPho) < 0.6) overlap = true;
      }
      if(!overlap){
	halfcleanedPhotons.insert(mmid::value_type(ipho->first, ipho->second ) );
      } 
      if(overlap && weird){
	cout << "Photon " << ipho->second << " removed due to cleaning!" << endl; 
      }
    }
  
    //Remove fakes that overlap within dR < 0.4 of an electron, muon or photon (same cleaning as for jets)
    for (mmid::const_iterator ifake = uncleanedFakes.begin(); ifake != uncleanedFakes.end(); ifake++){
      bool overlap = false;
      thisFake.SetPtEtaPhiE(ifake->first, 
				  event->phoEta_->at(ifake->second), 
				  event->phoPhi_->at(ifake->second),
				  event->phoCalibE_->at(ifake->second));
      for(unsigned int imu = 0;imu<looseMuons.size(); imu++){
	if(getDR(looseMuons.at(imu),thisFake) < 0.4) overlap = true;
      }
      for(mmid::const_iterator  iele = halfcleanedEles.begin(); iele != halfcleanedEles.end(); iele++){
	thisEle.SetPtEtaPhiE(iele->first,
				   event->phoEta_->at(iele->second),
				   event->phoPhi_->at(iele->second),
				   event->phoCalibE_->at(iele->second));
	if(getDR(thisEle, thisFake) < 0.6) overlap = true;
      }
      for(mmid::const_iterator ipho = halfcleanedPhotons.begin(); ipho != halfcleanedPhotons.end(); ipho++){
	thisPho.SetPtEtaPhiE(ipho->first,
				   event->phoEta_->at(ipho->second),
				   event->phoPhi_->at(ipho->second),
				   event->phoCalibE_->at(ipho->second));
        if(getDR(thisPho, thisFake) < 0.6) overlap = true;
      }
      if(!overlap){
	halfcleanedFakes.insert(mmid::value_type(ifake->first, ifake->second ) );
      }
      else if (weird){
	cout << "Fake Photon " << ifake->second << " removed due to cleaning!" << endl;
      }	
    }

    //-----------------------------------------------------//
    //-------------Jets------------------------------------//
    int nJet = event->nJet_;
    int numCleanJets = 0;
    //vector<TLorentzVector> cleanJets;        //looseID+eta+pt+gamma/e/mu cuts
    cleanJets.clear();
    double jetPt_cut = 30.0;
    double jetEta_cut = 5.0;
    double closeDR_cut = 0.15;
    double closePt_cut = 0.2;
  
    //First remove jets that overlap within closeDR_cut and closePt_cut of photons, electrons, muons, or fakes
    for(int ijet=0; ijet<nJet; ijet++){            
      if( event->jetPt_->at(ijet) < jetPt_cut) continue;      
      if( event->jetPFLooseId_->at(ijet) ==0 ) continue;            
      if( TMath::Abs(event->jetEta_->at(ijet) ) > jetEta_cut ) continue;

      thisJet.SetPtEtaPhiM(event->jetPt_->at(ijet), event->jetEta_->at(ijet), event->jetPhi_->at(ijet), 0.00);
       
      bool overlapping = false;
      for (mmid::const_iterator ipho = halfcleanedPhotons.begin(); ipho != halfcleanedPhotons.end(); ipho++){
	thisPho.SetPtEtaPhiE(ipho->first,
				   event->phoEta_->at(ipho->second),
				   event->phoPhi_->at(ipho->second),
				   event->phoCalibE_->at(ipho->second));
	if(getDR(thisPho, thisJet) < closeDR_cut){
	  h_diffJetPTOverlapPT->Fill(abs(thisJet.Pt() - thisPho.Pt()));
	  h_percentDiffJetPTOverlapPT->Fill(100*abs(thisJet.Pt() - thisPho.Pt())/thisPho.Pt());
          if(abs(thisJet.Pt() - thisPho.Pt())/thisPho.Pt() < closePt_cut ) {
	    //            h_diffJetPTOverlapPT->Fill(abs(thisJet.Pt() - thisPho.Pt()));
            overlapping = true;
          }
        }//end of if (dR < closeDR_cut)
      }//end of photon loop     
      

      for (mmid::const_iterator iele = halfcleanedEles.begin(); iele != halfcleanedEles.end(); iele++){
	thisEle.SetPtEtaPhiE(iele->first,
				   event->phoEta_->at(iele->second),
				   event->phoPhi_->at(iele->second),
				   event->phoCalibE_->at(iele->second));
	if(getDR(thisEle, thisJet) < closeDR_cut){
	  if(abs(thisJet.Pt() - thisEle.Pt())/thisEle.Pt() < closePt_cut) {
	    h_diffJetPTOverlapPT->Fill(abs(thisJet.Pt() - thisEle.Pt()));
	    overlapping = true;
	  }
	}  
      }

      for(unsigned int imu = 0;imu<looseMuons.size(); imu++){
	if(getDR(looseMuons.at(imu),thisJet) < closeDR_cut){
	  if(abs(thisJet.Pt() - looseMuons.at(imu).Pt())/looseMuons.at(imu).Pt() < closePt_cut){
	    h_diffJetPTOverlapPT->Fill(abs(thisJet.Pt() - looseMuons.at(imu).Pt()));
	    overlapping = true; 
	  } //end of PT check
	} 
      } //end of muon loop
    
      for (mmid::const_iterator ifake = halfcleanedFakes.begin(); ifake != halfcleanedFakes.end(); ifake++){
	thisFake.SetPtEtaPhiE(ifake->first,
				    event->phoEta_->at(ifake->second),
				    event->phoPhi_->at(ifake->second),
				    event->phoCalibE_->at(ifake->second));
        if(getDR(thisFake, thisJet) < closeDR_cut ){
          if( abs(thisJet.Pt() - thisFake.Pt())/thisFake.Pt() < closePt_cut ) {
	    overlapping = true;
	    h_diffJetPTOverlapPT->Fill(abs(thisJet.Pt() - thisFake.Pt()));
	  }
	}
      }

      if(overlapping){
	if(weird) cout << "Throw away jet " << ijet << " for overlapping!" << endl;
	continue;
      }
      if (abs(getDphi(event->pfMETPhi_, thisJet.Phi()) ) < minJetMETdphi){
	minJetMETdphi = abs(getDphi(event->pfMETPhi_, thisJet.Phi()));
      }
      if(thisJet.Pt() < jetPt_cut){
	continue;
      }
      numCleanJets++; 
      cleanJets.push_back(thisJet);
    }
  
    if(weird) cout << "Number of clean jets is " << numCleanJets << endl;
    //--------------------------------------------------------------------
    //Now that we have clean jets, we need to throw away any photons, electrons, or fakes that overlap within rejectDR_cut of a clean jet
    float fakeRejectDR_cut = 0.0;
    float phoRejectDR_cut  = 0.0;
    for (mmid::const_iterator ipho = halfcleanedPhotons.begin(); ipho != halfcleanedPhotons.end(); ipho++){
      bool overlap = false;
      float minDR = 1000;
      float minDR2 = 1000;
      thisPho.SetPtEtaPhiE(ipho->first,
				 event->phoEta_->at(ipho->second),
				 event->phoPhi_->at(ipho->second),
				 event->phoCalibE_->at(ipho->second));

      for(unsigned int jj=0; jj<numCleanJets; jj++){
	float dR = getDR(thisPho, cleanJets[jj]);
	if(dR < minDR) minDR = dR;
	else if (dR < minDR2) minDR2 = dR;
	if(dR < phoRejectDR_cut && dR > 0.2){
	  overlap = true;
	}
      }
      phoJetDR->Fill(minDR);
      phoJetDR2->Fill(minDR2);
      if(!overlap){
	cleanedPhotons.push_back(thisPho);
	if(abs(event->phoSCEta_->at(ipho->second)) < 1.4442){
	  pho_Cands.insert(mmid::value_type(ipho->first, ipho->second ) );
	}
      }
    }

    for (mmid::const_iterator iele = halfcleanedEles.begin(); iele != halfcleanedEles.end(); iele++){
      bool overlap = false;
      thisEle.SetPtEtaPhiE(iele->first,
						event->phoEta_->at(iele->second),
						event->phoPhi_->at(iele->second),
						event->phoCalibE_->at(iele->second));
      for(unsigned int jj=0; jj<numCleanJets; jj++){
	if(getDR(thisEle, cleanJets[jj]) < phoRejectDR_cut && getDR(thisEle, cleanJets[jj]) > 0.2){
	  overlap = true;
	}
      }
      if(!overlap){
	cleanedEles.push_back(thisEle);
	if(abs(event->phoSCEta_->at(iele->second)) < 1.4442){
	  pho_Cands.insert(mmid::value_type(iele->first, iele->second ) );
	}
      }
    }
  
    for (mmid::const_iterator ifake = halfcleanedFakes.begin(); ifake != halfcleanedFakes.end(); ifake++){
      bool overlap = false;
      float minDR = 1000.;
      float minDR2 = 1000.;
      thisFake.SetPtEtaPhiE(ifake->first,
				  event->phoEta_->at(ifake->second),
				  event->phoPhi_->at(ifake->second),
				  event->phoCalibE_->at(ifake->second));
      for(unsigned int jj=0; jj<numCleanJets; jj++){
	float dR = getDR(thisFake, cleanJets[jj]);
	//	minDR =(getDR(thisFake, thisJet) < minDR) ? getDR(thisFake, thisJet) : minDR;
	/*if(getDR(thisFake, thisJet) < fakeRejectDR_cut ){
	  overlap = true;
	  }*/
	if(dR < minDR) minDR = dR;
	else if( dR < minDR2) minDR2 = dR;
      }
      //      if(minDR2 < fakeRejectDR_cut) overlap = true;
      if(minDR < fakeRejectDR_cut) overlap = true;    
      if(overlap && weird) cout << "Removing fake "<< ifake->second << " for overlapping with jet" << endl;
      fakeJetDR->Fill(minDR);
      fakeJetDR2->Fill(minDR2);
      if(!overlap){
	cleanedFakes.push_back(thisFake);
	if(abs(event->phoSCEta_->at(ifake->second)) < 1.4442){
	  pho_Cands.insert(mmid::value_type(ifake->first, ifake->second ) );
	}
      }
    }

    for(int Part_it = 0; Part_it < event->nMC_; Part_it++){
      float partPdgID = event->mcPID->at(Part_it);
      if(fabs(partPdgID)==22){
	if(event->mcPt->at(Part_it) < 5.0) continue;
	float minDR = 1000.;
	TLorentzVector part_it_vec = MassLorentzVector(event->mcPt->at(Part_it),
						       event->mcEta->at(Part_it),
						       event->mcPhi->at(Part_it),
						       event->mcMass->at(Part_it));
	for(unsigned int jj=0; jj<numCleanJets; jj++){
	  minDR =(getDR(part_it_vec, cleanJets[jj]) < minDR) ? getDR(part_it_vec, cleanJets[jj]) : minDR;
	}
	genPhoJetDR->Fill(minDR);
	if(minDR > phoRejectDR_cut && abs(part_it_vec.Eta()) < 1.4442){
	  bool matchToGood = false;
	  for(unsigned int ipho = 0; ipho<cleanedPhotons.size(); ipho++){
	    if(isSameObject(cleanedPhotons[ipho], part_it_vec,0.01) && (part_it_vec.Pt() - cleanedPhotons[ipho].Pt()) / cleanedPhotons[ipho].Pt() < 0.05){
	      phoID_PhoJetDR->Fill(true, minDR);
	      matchToGood = true;
	      numPassGenPho++;
	      break;
	    }
	  }
	  if(!matchToGood) phoID_PhoJetDR->Fill(false, minDR);
	  numTotalGenPho++;
	}
      }
    }
  


    //-----------------------------------------------------------------------------
    //-------------------------------Find photon candidates to use ----------------
     
    if(printLevel>0)cout<<"phoCands size= "<<pho_Cands.size()<<endl;
    if(printLevel>0)cout<<"About to create PhoOne and PhoTwo"<<endl;

    if(weird){
      cout << "Photons included in sorting process: ";
      for(mmid::const_iterator pho_it = pho_Cands.begin();pho_it != pho_Cands.end(); pho_it++){
        cout << pho_it->second << ", ";
      }
      cout << endl;
    }

    bool breakPho=false;
    bool dogg=false,doee=false,doeg=false,doff=false,dogammafake=false,dogf=false,dofg=false,doef=false;
    bool use13, use23;
    bool doffmix = false, doffhigh = false;
    bool TwoPhos=false;
    float phoDR_cut =0.6;
    PhoOne = -1;
    PhoTwo = -1;
    PhoThree = -1;

    //If there are at least 2 photon candidates, check to make sure they are well separated
    if( pho_Cands.size()>=2){
      if(printLevel>0)cout<<"Looping over cleaned photon candidates"<<endl;

      //This makes sure the two pho objects are separated by dR > phoDR_cut 
      breakPho=false;
      mmid::const_iterator phoendminusone = pho_Cands.end();
      phoendminusone--;

      for(mmid::const_iterator pho_it = pho_Cands.begin(); pho_it != (phoendminusone); pho_it++){
        if(breakPho) break;
        //recall PhoOne is the index of the first photon  
        PhoOne= pho_it->second;
        if(printLevel>0)cout<<"PhoOne is "<<PhoOne<<endl;

        PhoOneVec.SetPtEtaPhiE(pho_it->first, 
				     event->phoEta_->at(pho_it->second), 
				     event->phoPhi_->at(pho_it->second),
				     event->phoCalibE_->at(pho_it->second));

        TLorentzVector PhoTwoCandVector, PhoThreeCandVector;
	mmid::const_iterator pho_it_next = pho_it;
        pho_it_next++;
        for(mmid::const_iterator pho_it2 = pho_it_next; pho_it2!=pho_Cands.end(); pho_it2++){
          PhoTwoCandVector.SetPtEtaPhiE(pho_it2->first, 
					      event->phoEta_->at(pho_it2->second), 
					      event->phoPhi_->at(pho_it2->second),
					      event->phoCalibE_->at(pho_it2->second));

          if( !isSameObject(PhoOneVec, PhoTwoCandVector,phoDR_cut)){
            PhoTwo= pho_it2->second;
            PhoTwoVec = PhoTwoCandVector;
            if(printLevel>0)cout<<"PhoTwo is "<<pho_it2->second <<endl;

            TwoPhos=true;
            pho_it_next = pho_it2;
            pho_it_next++;
            breakPho = true;

            for( mmid::const_iterator pho_it3 = pho_it_next; pho_it3 != pho_Cands.end(); pho_it3++){
              if(weird) cout << "Considering photon " << pho_it3->second << " for photon 3."<< endl;
              if(pho_it3->second > nTotPhotons) cout << "PROBLEM with pho_it3!!!" << endl;

              PhoThreeCandVector.SetPtEtaPhiE(pho_it3->first, 
						    event->phoEta_->at(pho_it3->second), 
						    event->phoPhi_->at(pho_it3->second),
						    event->phoCalibE_->at(pho_it3->second));

              if(   !isSameObject(PhoOneVec, PhoThreeCandVector,phoDR_cut) && 
		    !isSameObject(PhoTwoVec, PhoThreeCandVector, phoDR_cut) ){

                PhoThree = pho_it3->second;
                PhoThreeVec = PhoThreeCandVector;
		if(weird) cout << "Assigning pho three to " << PhoThree << endl;
                break;
              }
              else if(weird) cout<<"Skipping photon " << pho_it3->second << " because it overlaps with photon 1 or photon 2 within " << phoDR_cut << endl;
            }
            break;  //Use break to stop looping with pho_it2
          }           //end check that the two photon candidates are not the same object
          else{      //if photons are too close, keep trying to find a second candidate 
            if(printLevel>0)cout <<"Photons fail DR separation requirement"<< "  Run: "<<event->run_<<"  Event: "<<event->event_<<"\n";
            if(weird) cout << "PhosFailDR!"<< "  Run: "<<event->run_<<"  Event: "<<event->event_<< endl;
          }
        }//end pho_it2 loop
      } //end pho_it loop
    }//end of if two photons

    if(!TwoPhos){  //If there aren't two photon candidates passing all the selection criteria, move on to the next event
      if(weird)  cout << " Less than two good photons!" << endl;
      if(printLevel>0)cout<<"Failed to find to phos. Moving to the next event."<<endl;
      continue;
    }
    nTwoCands++;
    //From here on out we are only looking at cases with 2 photon candidates
    if(printLevel>0)cout<<"Found two photons"<<endl;

    //--------------------------------------------------------------------------------
    //------------------------Create event level quantities, HT, MHT, ST---------------

    bool closeJetMET = (minJetMETdphi < 0.3);

    //obtain HT and MHT from cleaned jets
    float ST = 0;
    float mHT = 0;
    TLorentzVector MHTvec;
    int nHTJets = numCleanJets;
    float highestJetEnergy = 0;
    for(unsigned int jj=0; jj<numCleanJets; jj++){
      HT     = HT + cleanJets[jj].Pt();
      MHTvec = MHTvec + cleanJets[jj];
      ST = ST + cleanJets[jj].Pt();
      if(cleanJets[jj].Pt() > highestJetEnergy){
	highestJetEnergy = cleanJets[jj].Pt();
      }
    }
    //Add electron pt to ST
    for(unsigned int iele = 0;iele<cleanedEles.size(); iele++){
      ST = ST + cleanedEles.at(iele).Pt();
      MHTvec = MHTvec + cleanedEles.at(iele);
    }
    //Now add the contribution from muons
    for(unsigned int imu = 0;imu<looseMuons.size(); imu++){
      ST = ST + looseMuons.at(imu).Pt();
      MHTvec = MHTvec + looseMuons.at(imu);
    } 
    //Now add photon contribution
    for(unsigned int ipho = 0; ipho<cleanedPhotons.size(); ipho++){
      ST = ST+ cleanedPhotons.at(ipho).Pt();
      MHTvec = MHTvec + cleanedPhotons.at(ipho);
    }
    //Now add fake contribution 
    for(unsigned int ifake = 0; ifake<cleanedFakes.size(); ifake++){
      ST = ST+ cleanedFakes.at(ifake).Pt();
      MHTvec = MHTvec + cleanedFakes.at(ifake);
    }
    //	cout<<"HT = " << HT<<endl;
    TLorentzVector MHtvec = - MHTvec;
    mHT = MHtvec.Pt();
    float P = -1;
    if(HT > 0 ){
      P = mHT/sqrt(HT);
      if( numCleanJets == 0) cout << "HT > 0 but numCleanJets is zero! HT = " << HT << endl;
    }
    else if (numCleanJets > 0){
      cout << "HT = 0 but numCleanJets is nonzero! numCleanJets = " << numCleanJets << endl;
      P = -1;
    }
    float diPT = PhoOneVec.Pt() + PhoTwoVec.Pt();

    //    float fracSTEM = diPT/ST;


    //-------------Electron study- how many electrons are we missing by only using the PF Photon collection? This matters for HT
    bool overlapping =false;
    float missingElectronHT = 0;
    int numMissingElectrons =0;
    for(int it_ele = 0; it_ele < event->nEle_; it_ele++){
      if(event->eleIDbit_->at(it_ele)>>2&1 && abs(event->eleSCEta_->at(it_ele) )< 2.4){
	thisEle.SetPtEtaPhiE(event->eleCalibPt_->at(it_ele), 
			     event->eleEta_->at(it_ele), 
			     event->elePhi_->at(it_ele),
			     event->eleCalibEn_->at(it_ele));
	overlapping = false;
	for(unsigned int ipho = 0; ipho<cleanedPhotons.size(); ipho++){
	  if(getDR(cleanedPhotons.at(ipho), thisEle) < 0.4) overlapping = true; //remove electrons that overlap with photons
	}
	for(unsigned int iele = 0;iele< cleanedEles.size(); iele++){
	  if(getDR(cleanedEles.at(iele),thisEle) < 0.4) overlapping = true; //remove electrons that overlap with electrons (ie pixel seeded photons) 
	}
	for(unsigned int imu = 0;imu<looseMuons.size(); imu++){
	  if(getDR(looseMuons.at(imu),thisEle) < 0.4) overlapping = true; //remove electrons that overlap with muons   
	}
	for(unsigned int ifake = 0;ifake< cleanedFakes.size(); ifake++){
	  if(getDR(cleanedFakes.at(ifake),thisEle) < 0.4) overlapping = true; //remove electrons that overlap with fakes
	}
	for(unsigned int jj=0; jj<numCleanJets; jj++){
	  if( getDR(cleanJets[jj], thisEle) < 0.4) overlapping = true; //remove electrons that have already been counted as jets
	}
	if(!overlapping){
	  h_PFElectronPT->Fill(event->eleCalibPt_->at(it_ele));
	  missingElectronHT += event->eleCalibPt_->at(it_ele);
	  numMissingElectrons++;
	  //	  if(event->eleCalibPt_->at(it_ele) >30) cout << "Missing PF Electron with PT = " << event->eleCalibPt_->at(it_ele) << endl;
	}
      }
    }
    if(numMissingElectrons > 0 ){
      h_NumPFElectrons->Fill(numMissingElectrons);
      h_ElectronMissingHT->Fill(missingElectronHT);
    }
      
    //-------------End of Electron study------------------------------ 
        
        
    //Sort the events!
    if(weird) cout << "About to categorize the events! Using photons " <<PhoOne <<", " << PhoTwo<<" and " << PhoThree<<endl;
    CategorizeEvents(PhoOne,PhoTwo,PhoThree,Rho,dogg,doee,doeg,doff,dogammafake,dogf,dofg,doef,use13,use23,doffmix,doffhigh);
    if(use13 && use23) cout << "Problem! Set to use both 13 and 23!" << endl;
    if(use13){
      PhoTwoVec = PhoThreeVec;
      PhoTwo = PhoThree;
    }
    if(use23){
      PhoOneVec = PhoTwoVec;
      PhoOne = PhoTwo;
      PhoTwo = PhoThree;
      PhoTwoVec = PhoThreeVec;
    }
    InvMass=InvariantMass(PhoOneVec, PhoTwoVec);        
    diEMpt=GetDiEmPt(PhoOneVec,PhoTwoVec);
    if(diEMpt > 1000){ 
      //      cout << "OVERFLOW diempt " << diEMpt << " set to 1000 GeV" << endl;
      diEMpt = 800.0;
    }
    dPhi = std::fabs(TVector2::Phi_mpi_pi(event->phoPhi_->at(PhoOne) - event->phoPhi_->at(PhoTwo)));    
    TLorentzVector diEMPtVec = PhoOneVec + PhoTwoVec;
    float diEmPtMETDPhi = getDphi(diEMPtVec.Phi(),metvec->Phi());

    float mmdiempt = sqrt((diEMPtVec.Px()-metvec->Px())*(diEMPtVec.Px()-metvec->Px()) + (diEMPtVec.Py()-metvec->Py())*(diEMPtVec.Py()-metvec->Py()) ) ;

    bool eleVeto = false;
    //veto events with extra electrons
    if(weird) cout << "Number of vetoEles = " << vetoElesIndices.size() << endl;
    for(unsigned int i = 0;i<vetoElesIndices.size(); i++){
      int iele = vetoElesIndices[i]; 
      TLorentzVector eleVec = PhoLorentzVector(event->eleCalibPt_->at(iele), 
					       event->eleEta_->at(iele), 
					       event->elePhi_->at(iele),
					       event->eleCalibEn_->at(iele));
      if(!isSameObject(eleVec, PhoOneVec, 0.3) && !isSameObject(eleVec,PhoTwoVec, 0.3) ) eleVeto = true;
      if(weird) cout << "Checking DR between veto ele and photon objects: " << getDR(eleVec, PhoOneVec) << " " << getDR(eleVec, PhoTwoVec) << endl;
    }

    if(weird && (eleVeto || endcapEle) ) cout <<"Reject event because of electron veto! " << endl;

    if(eleVeto || endcapEle) continue;

    bool otherCut = true; //numCleanJets >= 2);//true; //(ST >= 400 && numCleanJets >= 2); 
    if (!otherCut) continue;

    //Calculate MT2. Set mass of photons and gravitinos to 0. 
    double MT2 =  asymm_mt2_lester_bisect::get_mT2(0, PhoOneVec.Px(), PhoOneVec.Py(),
						   0, PhoTwoVec.Px(), PhoTwoVec.Py(),
						   metvec->Px(),metvec->Py(),
						   0,0,0);

   
    if(weird){
      cout << "Event number " << event->event_ << endl;
      cout << "Invariant mass = " << InvMass << endl;
      cout << "PhoOne = " << PhoOne << ", and PhoTwo = "  << PhoTwo << endl;
      cout << "ee event? " << doee << endl;
      cout << "gg event? " << dogg <<endl;
      cout << "ff event? " << doff << endl;
      cout << "gf event? " << dogf << endl;
      cout << "diempt = " << diEMpt<< endl;
      cout << "Number of jets" << numCleanJets << endl;
      for(unsigned int i = 0;i<cleanJets.size(); i++){
	cout << "Jet " << i << " has pT " << cleanJets.at(i).Pt() << " and phi " << cleanJets.at(i).Phi() << " and eta " << cleanJets.at(i).Eta() << endl;; 
      }
    }

    if(dogg && doSignal) {
      h_DoublePhoton_Signal_MET_NoTriggerCut->Fill(met);
      h_GridAllGGEvents_noInvMass->Fill(gmass,nmass);
    }
    if( dogg && InvMass > 105 && passggTrigger){ //gg specific requirements
      //cout << "gg Event! Event Number " << event->event_ << endl;

      isGG = true;
      leadPt = PhoOneVec.Pt();
      trailPt = PhoTwoVec.Pt();
      leadR9 = event->phoR9_->at(PhoOne);
      trailR9 = event->phoR9_->at(PhoTwo);

      if(weird) cout << "Yes, this is a gg event I have. DR between photons is "<< getDR(PhoOneVec,PhoTwoVec) << " and invariant mass is " << InvMass << endl;

      if(weird){
	PhoOneUnCalibVec.SetPtEtaPhiE(event->phoEt_->at(PhoOne),
				      event->phoEta_->at(PhoOne),
				      event->phoPhi_->at(PhoOne),
				      event->phoE_->at(PhoOne));
	PhoTwoUnCalibVec.SetPtEtaPhiE(event->phoEt_->at(PhoTwo),
                                      event->phoEta_->at(PhoTwo),
                                      event->phoPhi_->at(PhoTwo),
				      event->phoE_->at(PhoTwo));
	cout << "Invariant mass (uncalib) = " << InvariantMass(PhoOneUnCalibVec, PhoTwoUnCalibVec) << endl;
      }
      //      if( met > 100) cout << "gg event with met " << met << "! Event number " << event->event_  << " run number " << event->run_ << " lumi = " << event->lumis_ << endl;

      if(printLevel>0)cout<<"Inside gg"<<endl;
      ngg++;
      //      if( met > 100) myfile << event->event_ << "\n";  
      if (met < 50) nggMET0to50++;
      if( met > 50 && met < 100) nggMET50to100++;
      if(InvMass > 110) nggInvMass110++;
      h_ggMT2->Fill(MT2);
      h_ggMT2_vs_ST->Fill(MT2,HT);
      h_ggMT2_vs_MET->Fill(MT2,met);
      if(P >= 0)  h_ggMT2_vs_MetSig->Fill(MT2,P);

      //Try to find the gen level particles that match the photons most closely. The matching particle will have the smallest dR between it and the photon, and the difference between its pt and the photon's pt must be less than 15%
      float drClose1=999., drClose2=999.;
      int PhoOneGen, PhoTwoGen;
      bool goForGen=false,  goForGen1=false, goForGen2=false;
            
      for(int i =0; i < event->nMC_; i++){
	TLorentzVector part_it_vec = MassLorentzVector(event->mcPt->at(i), event->mcEta->at(i), event->mcPhi->at(i), event->mcMass->at(i));
	if(event->mcStatus->at(i)==1){
	  if(getDR(PhoOneVec, part_it_vec)<drClose1 && abs(event->phoCalibEt_->at(PhoOne) - event->mcPt->at(i))/(event->phoCalibEt_->at(PhoOne)) < 0.15){
	    drClose1 = getDR(PhoOneVec, part_it_vec);
	    PhoOneGen=i;
	    goForGen1=true;
	  }
	  if(getDR(PhoTwoVec, part_it_vec)<drClose2 && abs(event->phoCalibEt_->at(PhoTwo) - event->mcPt->at(i))/(event->phoCalibEt_->at(PhoTwo)) < 0.15){
	    drClose2 = getDR(PhoTwoVec, part_it_vec);
	    PhoTwoGen=i;
	    goForGen2=true;
	  }
	}
	if(goForGen1 && goForGen2){
	  goForGen=true; //This flag is set if both photon objects are matched to gen level objects.
	}
      }
            
      bool phoOneMatch = false, phoTwoMatch = false;
      if(goForGen1){
	//If the matching gen particle is a photon, fill histograms
	if(event->mcPID->at(PhoOneGen)==22){
	  phoOneMatch = true;
	  //fill sigma ieta ieta histogram to get the shape for real photons
	  h_realPhosSigIetaIeta->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoOne));
	}
	else{//If the matching gen particle isn't a photon, it is a fake
	  phoOneMatch = false;
	  //fill fake sigma ieta ieta histogram to get the shape for fakes
	  h_fakePhosSigIetaIeta->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoOne));
	}
      }
      else{
	//if the photon object is undefined, fill the fake histogram
	h_fakePhosSigIetaIeta->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoOne));
      }
      if(goForGen2){
	if(event->mcPID->at(PhoTwoGen)==22){
	  phoTwoMatch = true;
	  //fill sigma ieta ieta histogram to get the shape for real photons
	  h_realPhosSigIetaIeta->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoTwo));
	}
	else{
	  phoTwoMatch = false;
	  //fill fake sigma ieta ieta histogram to get the shape for fakes
	  h_fakePhosSigIetaIeta->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoTwo));
	}
      }
      else{
	//if the photon object is undefined, fill the fakes histogram
	h_fakePhosSigIetaIeta->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoTwo));
      }
      if(phoOneMatch && phoTwoMatch){//If both of the photons correspond to actual photons, then we have a real di-gamma event. Plot the MET of these events, to be compared with the MET spectrum of the DY to EE events we use to model the di-gamma distribution.
	h_realggHT->Fill(HT);
	h_realggMet->Fill(met);
	h_realggDiEMPt->Fill(diEMpt);
	h_ggcomposition->Fill(0.5);
	h_nHTJets_realgg->Fill(nHTJets); 
      }
      else if((phoOneMatch && !phoTwoMatch)||(!phoOneMatch && phoTwoMatch) ){//gf events
	h_ggcomposition->Fill(1.5);
      }
      else if(!phoOneMatch&&!phoTwoMatch){//ff events
	h_ggcomposition->Fill(2.5);
      }
      else{//if anything gets placed here then a mistake has been made. All events should eith be gg, gf, or ff.
	h_ggcomposition->Fill(3.5);
      }
      h_gg_JetMETdPhi_vs_Met->Fill(met,minJetMETdphi);
      if(HT > 200){
	h_gg_JetMETdPhi_vs_Met_ST200->Fill(met,minJetMETdphi);
	if(closeJetMET){
	  h_ggMETvsMHT_CloseJetMET_ST200->Fill(met,mHT);
	}
	h_ggMETvsMHT_ST200->Fill(met,mHT);
      }
      
      //      if(met >140) cout<<"ggEvent with met > 140. Run " <<event->run_<<"  Event: "<<event->event_<< endl;

      h_gg_JetMETdPhi->Fill(minJetMETdphi);	
      if(closeJetMET){
	h_ggMETvsMHT_CloseJetMET->Fill(met,mHT);
      }
      h_ggSCE->Fill(event->phoSCE_->at(PhoOne));
      h_ggSCE->Fill(event->phoSCE_->at(PhoTwo));
      h_ggRawSCE->Fill(event->phoSCRawE_->at(PhoOne));
      h_ggRawSCE->Fill(event->phoSCRawE_->at(PhoTwo));
      h_ggSCE_diffRaw->Fill( abs(event->phoSCE_->at(PhoOne) - event->phoSCRawE_->at(PhoOne) ));
      h_ggSCE_diffRaw->Fill( abs(event->phoSCE_->at(PhoTwo) - event->phoSCRawE_->at(PhoTwo) ));
      h_ggSCE_normalRaw->Fill( event->phoSCE_->at(PhoOne) / event->phoSCRawE_->at(PhoOne) );
      h_ggSCE_normalRaw->Fill( event->phoSCE_->at(PhoTwo) / event->phoSCRawE_->at(PhoTwo) );


      h_ggEt_diffCalib->Fill(abs(event->phoEt_->at(PhoOne) - event->phoCalibEt_->at(PhoOne) ) );
      h_ggEt_diffCalib->Fill(abs(event->phoEt_->at(PhoTwo) - event->phoCalibEt_->at(PhoTwo) ) );

      h_ggPt->Fill(PhoOneVec.Pt());
      h_ggPt->Fill(PhoTwoVec.Pt());
      h_ggPtLead->Fill(PhoOneVec.Pt());
      h_ggPtTrail->Fill(PhoTwoVec.Pt());
      h_ggEta->Fill(PhoOneVec.Eta());
      h_ggPhi->Fill(PhoOneVec.Phi());
      h_ggEta->Fill(PhoTwoVec.Eta());
      h_ggPhi->Fill(PhoTwoVec.Phi());

      h_ggSumEt->Fill(event->pfMETsumEt_);
      h_ggSumEtVsST->Fill(ST, event->pfMETsumEt_);


      float minDR1 = 1000., minDR2 = 1000., dR = 0.;
      for(unsigned int i = 0;i<cleanJets.size(); i++){
	dR = getDR(PhoOneVec, cleanJets[i]) < minDR1; 
	if(dR < minDR1) minDR1 = dR;
	if(dR < 0.7) h_JetPT_overlappingPho->Fill(cleanJets[i].Pt());

	dR = getDR(PhoTwoVec, cleanJets[i]) < minDR2;
        if(dR < minDR2) minDR2 = dR;
	if(dR <0.7) h_JetPT_overlappingPho->Fill(cleanJets[i].Pt());
      }


      h_ggMinDR_vs_NJet->Fill(min(minDR1,minDR2),numCleanJets);
      h_ggNJet_vs_MET->Fill(numCleanJets,met);
      h_ggMinDR_vs_MET->Fill(min(minDR1,minDR2),met);
      h_ggNVertex_vs_MET->Fill(NVertex,met);
      h_ggRho_vs_MET->Fill(Rho,met);
      h_ggRho->Fill(Rho);
      bool blinded = false; //(met >= 100 && event->isData_);
      if(met >= 100 && event->isData_) signalfile << "Event: " << event->event_ << ", Run = " << event->run_ << ", Lumi = " << event->lumis_ << ", met = " << met << endl;

      if(!blinded){
	h_ggMet_vs_diPT->Fill(mHT, diPT);
	h_ggHT_vs_diPT->Fill(HT, diPT);
	h_ggMet->Fill(met,PUweight);
	h_ggMetUnweighted->Fill(met);
	crossCheckError->Fill(kTRUE,met);
	h_ggNJet->Fill(numCleanJets);

	h_ggMetMinusDiEM->Fill(mmdiempt);
	if( met > 100 && met < 150) ggMet100to150++;
	else if( met > 150 && met < 250) ggMet150to250++;
	else if( met > 250) ggMet250++;
	h_ggHT->Fill(HT);
	h_ggMHT->Fill(mHT);
	h_ggMETvsMHT->Fill(met,mHT);
	h_ggMHT_vs_HT->Fill(mHT, HT);
        h_ggMet_vs_ST->Fill(met,ST);
	h_ggMet_vs_HT->Fill(met,HT);
	if(P >= 0) {
	  h_ggMetSig->Fill(P);
	  h_ggMetSigVsHT->Fill(HT,P);
	  h_ggMetSigVsST->Fill(ST, P);
	  h_ggMetSigVsMET->Fill(met, P);
	  //	  if(fracSTEM < 0.7) h_ggMetSigVsST_lowFracSTEM->Fill(ST, P);
	}
      }
      h_nHTJets_gg->Fill(nHTJets);
      h_ggHighestJetEnergy->Fill(highestJetEnergy);
      for(int i =0; i < nHTJets; i++){
        h_ggJetEnergy->Fill(cleanJets[i].Pt() );
      }

      if(!event->isData_){
	for(int Part_it = 0; Part_it < event->nMC_; Part_it++){
	  float partPdgID = event->mcPID->at(Part_it);
	  if( fabs(partPdgID)==22 && fabs(event->mcMomPID->at(Part_it)) == 1000023){ //10023 = neutralino
	    TLorentzVector part_it_vec = MassLorentzVector(event->mcPt->at(Part_it), event->mcEta->at(Part_it), event->mcPhi->at(Part_it), event->mcMass->at(Part_it));	    
	    if(isSameObject(part_it_vec,PhoOneVec,0.01) && abs(part_it_vec.Pt() - PhoOneVec.Pt()) / PhoOneVec.Pt() < 0.05){
	      h_ggDiffVtx_y->Fill(abs(event->mcVty->at(Part_it) - event->vty_) );
	      h_ggDiffVtx_z->Fill(abs(event->mcVtz->at(Part_it)- event->vtz_) );
	      h_ggDiffVtx_x->Fill(abs(event->mcVtx->at(Part_it)- event->vtx_) );
	    }
	  }
	}
      }
      


      h_ggSigIetaIeta->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoOne));
      h_ggSigIetaIeta->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoTwo));
      h_ggR9->Fill(event->phoR9_->at(PhoOne) ) ;
      h_ggR9->Fill(event->phoR9_->at(PhoTwo) ) ;
      h_ggDiEMPt->Fill(diEMpt);
      h_ggDiEMPt_vs_MET->Fill(diEMpt,met);
      h_ggDiEMPt_vs_HT->Fill(diEMpt,HT);
      h_ggDiEMPt_vs_NJet->Fill(diEMpt,numCleanJets);
      h_ggInvarMass->Fill(InvMass,PUweight);
      if(P >= 0){
	if(numCleanJets ==1) h_ggMetSigVsST_1Jets->Fill(ST, P);
	if(numCleanJets ==2) h_ggMetSigVsST_2Jets->Fill(ST, P);
	if(numCleanJets >=3) h_ggMetSigVsST_3pJets->Fill(ST, P);
	//	h_ggMetSigVsfracSTEM->Fill(fracSTEM, P);
	///h_ggSTVsfracSTEM->Fill(fracSTEM, ST);
      }


      if(doSignal){
	if(gmass == 2000){
	  if(nmass == 1900){ h_ggDiEMPt_2000_1900->Fill(diEMpt); }
          if(nmass == 1500){ h_ggDiEMPt_2000_1500->Fill(diEMpt); }
          if(nmass == 1000){ h_ggDiEMPt_2000_1000->Fill(diEMpt); }
          if(nmass == 500){  h_ggDiEMPt_2000_500->Fill(diEMpt); }
          if(nmass == 100){  h_ggDiEMPt_2000_100->Fill(diEMpt); }
	}
	else if(gmass == 1800){
          if(nmass == 1500){ h_ggDiEMPt_1800_1500->Fill(diEMpt); }
          if(nmass == 1000){ h_ggDiEMPt_1800_1000->Fill(diEMpt); }
          if(nmass == 500){  h_ggDiEMPt_1800_500->Fill(diEMpt); }
          if(nmass == 100){  h_ggDiEMPt_1800_100->Fill(diEMpt); }
        }
        else if(gmass == 1700){
          if(nmass == 1500){  h_ggDiEMPt_1700_1500->Fill(diEMpt); }
          if(nmass == 1000){  h_ggDiEMPt_1700_1000->Fill(diEMpt); }
          if(nmass == 500){ h_ggDiEMPt_1700_500->Fill(diEMpt); }
          if(nmass == 100){ h_ggDiEMPt_1700_100->Fill(diEMpt); }
        }
        else if(gmass == 1600){
          if(nmass == 1500){  h_ggDiEMPt_1600_1500->Fill(diEMpt); }
          if(nmass == 1000){  h_ggDiEMPt_1600_1000->Fill(diEMpt); }
          if(nmass == 500){   h_ggDiEMPt_1600_500->Fill(diEMpt); }
          if(nmass == 100){   h_ggDiEMPt_1600_100->Fill(diEMpt); }
        }
        else if(gmass == 1500){
          if(nmass == 1400){ h_ggDiEMPt_1500_1400->Fill(diEMpt); }
          if(nmass == 1000){ h_ggDiEMPt_1500_1000->Fill(diEMpt); }
          if(nmass == 500){  h_ggDiEMPt_1500_500->Fill(diEMpt); }
          if(nmass == 100){  h_ggDiEMPt_1500_100->Fill(diEMpt); }
        }
	if(diEMpt > 400) h_ggDiEMPt400_Grid->Fill(gmass,nmass);
        if(diEMpt > 300) h_ggDiEMPt300_Grid->Fill(gmass,nmass);

        h_GridAllGGEvents->Fill(gmass,nmass);
	h_GridAllGGEvents2->Fill(gmass,nmass,2);
	float met_jesup = event->pfMET_T1JESDo_;
	float met_jesdown = event->pfMET_T1JESUp_;


	//Bin 1
	if(met >= 100 && met < 115){
	  h_MET100to115_GridUnweighted->Fill(gmass, nmass);
	  if(event->nVtx_ < 20){
	    h_MET100to115_GridGG_LowNvtx->Fill(gmass, nmass);
	  }
	  else{
	    h_MET100to115_GridGG_HighNvtx->Fill(gmass, nmass);
	  }
	} 
	//Bin 2
	else if(met >= 115 && met < 130){
	  h_MET115to130_GridUnweighted->Fill(gmass, nmass);
	  if(event->nVtx_< 20){
	    h_MET115to130_GridGG_LowNvtx->Fill(gmass, nmass);
	  }
	  else{
	    h_MET115to130_GridGG_HighNvtx->Fill(gmass, nmass);
	  }
	}
	//Bin 3
	else if(met >= 130 && met < 150){
	  h_MET130to150_GridUnweighted->Fill(gmass, nmass);
	  if(event->nVtx_< 20){
	    h_MET130to150_GridGG_LowNvtx->Fill(gmass, nmass);
	  }
	  else{
	    h_MET130to150_GridGG_HighNvtx->Fill(gmass, nmass);
	  }
	}
	//Bin 4	
	else if(met >= 150 && met < 185){
	  h_MET150to185_GridUnweighted->Fill(gmass, nmass);
	  if(event->nVtx_< 20){
	    h_MET150to185_GridGG_LowNvtx->Fill(gmass, nmass);
	  }
	  else{
	    h_MET150to185_GridGG_HighNvtx->Fill(gmass, nmass);
	  }
	}
	//Bin 5
	else if(met >= 185 && met < 250){
	  h_MET185to250_GridUnweighted->Fill(gmass, nmass);
	  if(event->nVtx_< 20){
	    h_MET185to250_GridGG_LowNvtx->Fill(gmass, nmass);
	  }
	  else{
	    h_MET185to250_GridGG_HighNvtx->Fill(gmass, nmass);
	  }
	}
	//Bin 6
	else if(met >= 250){
	  h_MET250_GridUnweighted->Fill(gmass, nmass);
	  if(event->nVtx_< 20){
	    h_MET250_GridGG_LowNvtx->Fill(gmass, nmass);
	  }
	  else{
	    h_MET250_GridGG_HighNvtx->Fill(gmass, nmass);
	  }
	}

	float genMet = event->genMET_;

	if(genMet >= 100 && genMet < 115){
	  h_genMET100to115_GridUnweighted->Fill(gmass,nmass);
	}
	else if(genMet >= 115 && genMet < 130){
	  h_genMET115to130_GridUnweighted->Fill(gmass,nmass);
	}
	else if(genMet >= 130 && genMet < 150 ){
	  h_genMET130to150_GridUnweighted->Fill(gmass,nmass);
	}
	else if(genMet >= 150 && genMet < 185){
	  h_genMET150to185_GridUnweighted->Fill(gmass,nmass);
	}
	else if(genMet >= 185 && genMet < 250){
	  h_genMET185to250_GridUnweighted->Fill(gmass,nmass);
	}
	else if(genMet >= 250){
	  h_genMET250_GridUnweighted->Fill(gmass,nmass);
	}

	if(met_jesup >= 100 && met_jesup < 115){
	  h_MET100to115_GridUnweighted_JESUp->Fill(gmass,nmass);
	}
	else if(met_jesup >= 115 && met_jesup < 130){
	  h_MET115to130_GridUnweighted_JESUp->Fill(gmass,nmass);
	}
	else if(met_jesup >= 130 && met_jesup < 150){
	  h_MET130to150_GridUnweighted_JESUp->Fill(gmass,nmass);
	}
	else if(met_jesup >= 150 && met_jesup < 185){
	  h_MET150to185_GridUnweighted_JESUp->Fill(gmass,nmass);
	}
	else if(met_jesup >= 185 && met_jesup < 250){
	  h_MET185to250_GridUnweighted_JESUp->Fill(gmass,nmass);
	}
	else if(met_jesup >= 250){
	  h_MET250_GridUnweighted_JESUp->Fill(gmass,nmass);
	}
        if(met_jesup >= 150){
          h_MET150_GridUnweighted_JESUp->Fill(gmass,nmass);
        }

	if(met_jesdown >= 100 && met_jesdown < 115){
	  h_MET100to115_GridUnweighted_JESDown->Fill(gmass,nmass);
	}
	else if(met_jesdown >= 115 && met_jesdown < 130){
	  h_MET115to130_GridUnweighted_JESDown->Fill(gmass,nmass);
	}
	else if(met_jesdown >= 130 && met_jesdown < 150){
	  h_MET130to150_GridUnweighted_JESDown->Fill(gmass,nmass);
	}
	else if(met_jesdown >= 150 && met_jesdown < 185){
	  h_MET150to185_GridUnweighted_JESDown->Fill(gmass,nmass);
	}
	else if(met_jesdown >= 185 && met_jesdown < 250){
	  h_MET185to250_GridUnweighted_JESDown->Fill(gmass,nmass);
	}
	else if(met_jesdown >= 250){
	  h_MET250_GridUnweighted_JESDown->Fill(gmass,nmass);
	}
        if(met_jesdown >= 150){
          h_MET150_GridUnweighted_JESDown->Fill(gmass,nmass);
        }
      }

    }//end gg loop
    


    else if(doee && passeeTrigger){ 
      h_eeInvMass_vs_Met->Fill(met,InvMass, PUweight);
      if(InvMass < 105 && InvMass > 75){
	nee++;
	if(funky){
	  nee_funky++;
	  h_eeMet_Funky->Fill(met);
	}

	if(!event->isData_){
	  int pu_it = 0; //index of pu  
	  for(vector<int>::const_iterator bx = event->puBX_->begin(); bx != event->puBX_->end(); bx++){
	    if( *bx == 0)  h_puTrue_ee->Fill(event->puTrue_->at(pu_it), event->genWeight_);
	    pu_it++;
	  }
	}
	//	if( met > 50 && met < 60) myfile << event->event_ << "\n";

	//	if( met <= 115 && met >= 100) myfile << "ee Event with 100<=MET <= 115 GeV ! Event Number " << event->event_ << " run Number " << event->run_ << " lumi sect " << event->lumis_ << " and MET " << met << " GeV""\n";

	//	if( met < 75 && met > 50) myfile << "ee Event with MET between 50 and 75 GeV! Event Number " << event->event_ << "\n";
	//	  float eeinvmas = InvMass(PhoOneVec, PhoTwoVec);
	float diempt_weight=1; //The weights are stored in the histogram eeweights
	int bin = eeweights->FindFixBin(float(diEMpt));
	diempt_weight=eeweights->GetBinContent(bin);

	//	if(met > 115 && met < 130) cout << "Weight = " << diempt_weight << " * " << PUweight << ". Diempt = " << diEMpt << endl;
	float NJet_weight = 1;
	bin = h2D_weights_ee->FindFixBin(float(diEMpt), float(numCleanJets));
	NJet_weight = h2D_weights_ee->GetBinContent(bin);
        h_eeMetReweighted_NJet->Fill(met,NJet_weight*PUweight);
        if(!event->isData_){
          h_eeMetReweighted_NJet_JESDown->Fill(event->pfMET_T1JESDo_,NJet_weight*PUweight);
          h_eeMetReweighted_NJet_JESUp->Fill(event->pfMET_T1JESUp_,NJet_weight*PUweight);
	}
	h_eeMetMinusDiEM->Fill(mmdiempt);
        h_eeMetMinusDiEM_Reweighted->Fill(mmdiempt,diempt_weight);
	//if(met < 185 && met > 150){
	//cout << "Original weight = " <<  diempt_weight << ", and new weight = " << NJet_weight << endl;
	//}

	h_eeInvarMassFullRange->Fill(InvMass);
	h_ee_JetMETdPhi_vs_Met->Fill(met,minJetMETdphi);
	if(HT > 200){
	  h_ee_JetMETdPhi_vs_Met_ST200->Fill(met,minJetMETdphi);
	}
	h_ee_JetMETdPhi->Fill(minJetMETdphi);
	if(closeJetMET){
	  h_eeMETvsMHT_CloseJetMET->Fill(met,mHT);
	}
        h_eeMet_vs_diPT->Fill(mHT, diPT);
        h_eeHT_vs_diPT->Fill(HT, diPT);
	//Study C   
	h_eeMet_vs_diPT->Fill(mHT, diPT);
	h_eeHT_vs_diPT->Fill(HT, diPT);
	/*	if(P >= 0){
	  if(ST>= 80 && ST<100) h_ee_fracSTEM_ST80100->Fill(fracSTEM);
	  if(ST>= 100 && ST<120) h_ee_fracSTEM_ST100120->Fill(fracSTEM);
	  if(ST>= 120 && ST<150) h_ee_fracSTEM_ST120150->Fill(fracSTEM);
	  if(ST>= 150 && ST<200) h_ee_fracSTEM_ST150200->Fill(fracSTEM);
	  if(ST>= 200 ) h_ee_fracSTEM_ST200->Fill(fracSTEM);
	  h_eeMetSigVsfracSTEM->Fill(fracSTEM, P);
	  h_eeSTVsfracSTEM->Fill(fracSTEM, ST);
	  }*/
	//End Study C 
	h_eeMet->Fill(met,PUweight);
	h_eeNJet->Fill(numCleanJets);
	h_eeNJet_Reweighted->Fill(numCleanJets,diempt_weight);
	if(P > 3 && ST > 400) h_eeMet_MetSig3ST400->Fill(met);
	h_eeHT->Fill(HT);
	h_eeST->Fill(ST);
	h_eeMHT->Fill(mHT);
	if(P >= 0){
	  h_eeMetSig->Fill(P);
	  //if(numCleanJets >= 2){
	  if(ST >150 && ST < 400) h_eeMetSig_ST150400->Fill(P);
	  if(ST > 400) h_eeMetSig_ST400->Fill(P);
	    //}
	}
	h_eeMet_vs_HT->Fill(met,HT);
	h_eeMet_vs_ST->Fill(met,ST);
	if(numCleanJets >= 2 && diEMpt < 10){
	  h_eeInvarMass_NJet2Diempt10->Fill(InvMass);
	  h_nJets_eeNJet2Diempt10->Fill(numCleanJets);
	  h_nMuons_eeNJet2Diempt10->Fill(looseMuons.size());
	  h_nPhotons_eeNJet2Diempt10->Fill(cleanedPhotons.size());
	  h_nElectrons_eeNJet2Diempt10->Fill(cleanedEles.size());
	  h_eeDiEMPt_NJet2Diempt10->Fill(diEMpt);
	  h_ee_JetMETdPhi_NJet2Diempt10->Fill(minJetMETdphi);
	  h_ee_diEMPtMETdPhi_NJet2Diempt10->Fill(diEmPtMETDPhi);
	  h_eeMet_NJet2Diempt10->Fill(met);
	}
	h_nMuons_ee->Fill(looseMuons.size());
	h_nPhotons_ee->Fill(cleanedPhotons.size());
	h_nElectrons_ee->Fill(cleanedEles.size());
	h_ee_diEMPtMETdPhi->Fill(diEmPtMETDPhi);

	//	if(met > 140) cout << "ee event with MET > 140! Event number "<< event->event_ << ", lumi "<< event->lumis_ << " and run number "<<event->run_ << endl;

	//if(met>50){
	if(P >= 0) {
	  h_eeMetSigVsMET->Fill(met,P);
	  if(numCleanJets ==1) h_eeMetSigVsMET_1Jets->Fill(met,P);
	  else if(numCleanJets ==2) h_eeMetSigVsMET_2Jets->Fill(met,P);
	  else if(numCleanJets >=3) h_eeMetSigVsMET_3pJets->Fill(met,P);

	  h_eeMetSigVsHT->Fill(HT,P);
	  h_eeMetSigVsST->Fill(ST,P);
	  //          if(fracSTEM <0.7) h_eeMetSigVsST_lowFracSTEM->Fill(ST, P);
	}
	if( P >= 10 && ST < 150){
	  h_eeMet_StudyB->Fill(met);
	  h_eeMHT_StudyB->Fill(mHT);
	  h_eeST_StudyB->Fill(ST);
	  h_eeHT_StudyB->Fill(HT);
	  h_eeInvarMass_StudyB->Fill(InvMass);
	  h_nJets_eeStudyB->Fill(numCleanJets);
          h_nMuons_eeStudyB->Fill(looseMuons.size());
          h_nPhotons_eeStudyB->Fill(cleanedPhotons.size());
          h_nElectrons_eeStudyB->Fill(cleanedEles.size());
          h_eeDiEMPt_StudyB->Fill(diEMpt);
          h_ee_JetMETdPhi_StudyB->Fill(minJetMETdphi);
          h_ee_diEMPtMETdPhi_StudyB->Fill(diEmPtMETDPhi);

	}
        h_eeSumEt->Fill(event->pfMETsumEt_);
	h_eeSumEtVsST->Fill(ST, event->pfMETsumEt_);
	h_eeMT2->Fill(MT2);
	h_eeMT2_vs_ST->Fill(MT2,ST);
        h_eeMT2_vs_HT->Fill(MT2,HT);
	h_eeMT2_vs_MET->Fill(MT2,met);
	if(P >= 0){
	  h_eeMT2_vs_MetSig->Fill(MT2,P);
	  if(numCleanJets == 1){
	    h_eeMetSigVsST_1Jets->Fill(ST,P);
	    h_eeMetSigVsHT_1Jets->Fill(HT,P);
	  }
	  else if(numCleanJets == 2){
	    h_eeMetSigVsST_2Jets->Fill(ST,P);
	    h_eeMetSigVsHT_2Jets->Fill(HT,P);
	  }
	  else if(numCleanJets >= 3){
	    h_eeMetSigVsST_3pJets->Fill(ST,P);
	    h_eeMetSigVsHT_3pJets->Fill(HT,P);
	    if(numCleanJets == 3) h_eeMetSigVsST_3Jets->Fill(ST,P);
	    if(numCleanJets == 4)h_eeMetSigVsST_4Jets->Fill(ST,P);
	    if(numCleanJets == 5)h_eeMetSigVsST_5Jets->Fill(ST,P);
	    if(numCleanJets >= 6)h_eeMetSigVsST_6pJets->Fill(ST,P);
	  }
	}
	h_eeMETvsMHT->Fill(met,mHT);
	h_eeMet_vs_HT->Fill(met,HT);
	h_nHTJets_ee->Fill(nHTJets);
        h_eeHighestJetEnergy->Fill(highestJetEnergy);
	for(int i =0; i < nHTJets; i++){
	  h_eeJetEnergy->Fill(cleanJets[i].Pt() );
	}

	//Fill toy histograms
	if(doToys){
	  toys->cd();
	  float toy_diempt_weight=1; //The weights are stored in the histograms in RatioGausEE
	  bin = eeweights->FindFixBin(float(diEMpt));
	  for(int k=0; k<1000; k++){
	    toy_diempt_weight = RatioGausEE[k]->GetBinContent(bin);
	    ReweightedEE[k]->Fill(met,toy_diempt_weight);
	  }
	  toy_diempt_weight=1;
	  bin =  h2D_weights_ee->FindFixBin(float(diEMpt), float(numCleanJets));
	  for(int k=0; k<1000; k++){
	    toy_diempt_weight = RatioGaus2DEE[k]->GetBinContent(bin);
	    ReweightedEE2D[k]->Fill(met,toy_diempt_weight);
	  }
	fout->cd();
	}
	//Done with toy histograms
	h_eeMetReweighted->Fill(met,diempt_weight*PUweight);
	if(!event->isData_){
	  h_eeMetReweighted_JESDown->Fill(event->pfMET_T1JESDo_,diempt_weight*PUweight);
	  h_eeMetReweighted_JESUp->Fill(event->pfMET_T1JESUp_,diempt_weight*PUweight);
	  h_eeMetUnweighted->Fill(met,event->genWeight_);
	}
	else{ h_eeMetUnweighted->Fill(met); }

	h_eeMHTReweighted->Fill(mHT, diempt_weight);
	h_eeHTReweighted->Fill(HT, diempt_weight);
	//	if(diEMpt < 10 && numCleanJets >= 2) cout<<"ee event with low diempt! Diempt = " << diEMpt << " Event number " << event->event_ << " and run number " << event->run_ << " and number of jets " << numCleanJets << endl;
	h_eeDiEMPt->Fill(diEMpt);
	h_eeDiEMPt_vs_MET->Fill(diEMpt,met);
	h_eeDiEMPtWeight_vs_MET->Fill(diempt_weight,met);
	h_eeDiEMPt_vs_HT->Fill(diEMpt,HT);
	h_eeDiEMPt_vs_NJet->Fill(diEMpt,numCleanJets);
	h_eePt->Fill(PhoOneVec.Pt());
	h_eeSCE->Fill(event->phoSCE_->at(PhoOne));
	h_eeSCE->Fill(event->phoSCE_->at(PhoTwo));
	h_eeRawSCE->Fill(event->phoSCRawE_->at(PhoOne));
	h_eeRawSCE->Fill(event->phoSCRawE_->at(PhoTwo));
	h_eeSCE_diffRaw->Fill( abs(event->phoSCE_->at(PhoOne) - event->phoSCRawE_->at(PhoOne) ));
	h_eeSCE_diffRaw->Fill( abs(event->phoSCE_->at(PhoTwo) - event->phoSCRawE_->at(PhoTwo) ));
	h_eeSCE_normalRaw->Fill( event->phoSCE_->at(PhoOne) / event->phoSCRawE_->at(PhoOne) );
	h_eeSCE_normalRaw->Fill( event->phoSCE_->at(PhoTwo) / event->phoSCRawE_->at(PhoTwo) );

	h_eePt->Fill(PhoTwoVec.Pt());
	h_eeSigIetaIeta->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoOne));
	h_eeSigIetaIeta->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoTwo));
	h_eeR9->Fill(event->phoR9_->at(PhoOne) ) ;
	h_eeR9->Fill(event->phoR9_->at(PhoTwo) ) ;
      }
    }//end of ee
        
        
        
    else if(doeg && passggTrigger){ 
      if( InvMass > 105){
	int EleIt = event->phohasPixelSeed_->at(PhoOne) ? PhoOne : PhoTwo;
	float transverseMass = TMath::Sqrt(2 * event->phoCalibEt_->at(EleIt) * event->pfMET_ * (1 - TMath::Cos( (event->pfMETPhi_) - event->phoPhi_->at(EleIt) ) ) );
	h_egTransverseMass->Fill(transverseMass);
	h_egTransverseMass_vs_MET->Fill(transverseMass,met);
	if(weird){
	  cout << "====================" << endl;
	  cout << "Met = " << event->pfMET_ << endl;
	  cout << "MetPhi = " << event->pfMETPhi_ << endl;
	  cout << "Electron pt = " << event->phoCalibEt_->at(EleIt) << endl;
	  cout << "Electron phi = " << event->phoPhi_->at(EleIt) << endl;
	  cout << "Transverse mass = " << transverseMass << endl;
	  cout << "====================" << endl;
	}
	egMTCut->Fill((transverseMass < 100 ), met);

	if(doSignal && met >250){
	  h_signalEGMass->Fill(gmass,nmass);
	}

	//if(transverseMass > 100 ) continue;

	float eMetdPhi = abs(event->phoPhi_->at(EleIt) - event->pfMETPhi_);
	//if(eMetdPhi < 0.3) continue;

	//	if(met > 250) myfile << "eg Event with MET>250 GeV and MT<100! Event Number " << event->event_ << "  run number  " << event->run_ << "\n";

	if(printLevel>0)cout<<"Inside eg"<<endl;
	neg++;
	if(met >= 100 && met < 115) numEGSignal[0]++;
	if(met >= 115 && met < 130) numEGSignal[1]++;
	if(met >= 130 && met < 150) numEGSignal[2]++;
	if(met >= 150 && met < 185) numEGSignal[3]++;
	if(met >= 185 && met < 250) numEGSignal[4]++;
	if(met >= 250) numEGSignal[5]++;
	//	if(met >= 135 && met < 140) myfile << "EG event with MET between 135 and 140: Event Number "<<event->event_ << "\n"; 
	/*	if(met >= 100 && met < 110) numEGSignal[0]++; 
        if(met >= 110 && met < 120) numEGSignal[1]++; 
        if(met >= 120 && met < 130) numEGSignal[2]++; 
        if(met >= 130 && met < 145) numEGSignal[3]++; 
        if(met >= 145 && met < 170) numEGSignal[4]++;
        if(met >= 170 && met < 200) numEGSignal[5]++;
        if(met >= 200 && met < 240) numEGSignal[6]++;
        if(met >= 240 && met < 290) numEGSignal[7]++;
        if(met >= 290) numEGSignal[8]++;*/

	if(numCleanJets >= 2){
	  if(met >= 100 && met < 150) numEG_100to150++;
	  else if(met >=150 && met < 250) numEG_150to250++;
	  else if(met >=250) numEG_250to350++;
	}

	h_egMHT->Fill(mHT);
	h_egMet->Fill(met);
	h_egDiEMPt->Fill(diEMpt);
	h_eg_JetMETdPhi->Fill(minJetMETdphi);
	h_eg_JetMETdPhi_vs_Met->Fill(met,minJetMETdphi);
	if(HT > 200){
	  h_eg_JetMETdPhi_vs_Met_ST200->Fill(met,minJetMETdphi);
	}
	if(closeJetMET){
	  h_egMETvsMHT_CloseJetMET->Fill(met,mHT);
	}
	if(P >= 0){
	  h_egMetSig->Fill(P);
	  h_egMetSigVsMET->Fill(met,P);

	  h_egMetSigVsHT->Fill(HT,P);
	  h_egMetSigVsST->Fill(ST,P);
	  if(numCleanJets == 1) h_egMetSigVsST_1Jets->Fill(ST,P);
	  if(numCleanJets == 2) h_egMetSigVsST_2Jets->Fill(ST,P);
	  if(numCleanJets >= 3) h_egMetSigVsST_3pJets->Fill(ST,P);
	}

	h_egMETvsMHT->Fill(met,mHT);
	h_egMet_vs_HT->Fill(met,HT);
	h_egMHT_vs_HT->Fill(mHT, HT);
	h_egMet_vs_HT->Fill(met,HT);
	h_egSigIetaIeta->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoOne));
	h_egSigIetaIeta->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoTwo));
	h_egPt->Fill(PhoOneVec.Pt());
	h_egPt->Fill(PhoTwoVec.Pt());
      }
    }//end of eg
        
        
        
    //the dogammafake flag is set to true if there is one photon and one fake, irrespective of which of the two objects has the higher Pt value
    if(dogammafake && InvMass > 105){ 

      int fake = event->phohasPixelSeed_->at(PhoOne) ? PhoTwo : PhoOne;
      bool low = event->phoR9_->at(fake) < 0.9;

      if (!low){
	h_gammafakeMet_HighR9->Fill(met);
      }

      else{
      numGammaFake++;		

      if(doSignal){
        if(gmass == 2000){
          if(nmass == 1900){ h_gfMet_2000_1900->Fill(met); }
          if(nmass == 1500){ h_gfMet_2000_1500->Fill(met);  }
          if(nmass == 1000){ h_gfMet_2000_1000->Fill(met);}
          if(nmass == 500){  h_gfMet_2000_500->Fill(met); }
          if(nmass == 100){  h_gfMet_2000_100->Fill(met); }
        }
        else if(gmass == 1800){
          if(nmass == 1500){ h_gfMet_1800_1500->Fill(met);  }
          if(nmass == 1000){ h_gfMet_1800_1000->Fill(met); }
          if(nmass == 500){ h_gfMet_1800_500->Fill(met); }
          if(nmass == 100){ h_gfMet_1800_100->Fill(met); }
        }
	else if(gmass == 1700){
	  if(nmass == 1500){ h_gfMet_1700_1500->Fill(met);}
	  if(nmass == 1000){ h_gfMet_1700_1000->Fill(met); }
	  if(nmass == 500){ h_gfMet_1700_500->Fill(met); }
	  if(nmass == 100){ h_gfMet_1700_100->Fill(met); }
	}
	else if(gmass == 1600){
	  if(nmass == 1500){ h_gfMet_1600_1500->Fill(met); }
	  if(nmass == 1000){ h_gfMet_1600_1000->Fill(met);}
	  if(nmass == 500){ h_gfMet_1600_500->Fill(met); }
	  if(nmass == 100){ h_gfMet_1600_100->Fill(met);  }
	}
	else if(gmass == 1500){
	  if(nmass == 1400){ h_gfMet_1500_1400->Fill(met);}
	  if(nmass == 1000){ h_gfMet_1500_1000->Fill(met);}
	  if(nmass == 500){ h_gfMet_1500_500->Fill(met); }
	  if(nmass == 100){ h_gfMet_1500_100->Fill(met); }
	}
      }
      float diempt_weight=1; //The weights are stored in the histogram gfweights
      int bin = gfweights->FindFixBin(float(diEMpt));
      diempt_weight=gfweights->GetBinContent(bin);
      h_gammafakeMet->Fill(met);
      h_gfDiEMPt ->Fill(diEMpt);
      h_gammafake_JetMETdPhi->Fill(minJetMETdphi);
      h_gammafake_JetMETdPhi_vs_Met->Fill(met,minJetMETdphi);
      if(HT > 200){
        h_gammafake_JetMETdPhi_vs_Met_ST200->Fill(met,minJetMETdphi);
      }
      if(closeJetMET){
        h_gammafakeMETvsMHT_CloseJetMET->Fill(met,mHT);
      }
      h_gammafakeMetReweighted->Fill(met,diempt_weight);
      h_gammafakeHTReweighted->Fill(HT,diempt_weight);
      h_gammafakeHT->Fill(HT);
      h_gammafakeMHT->Fill(mHT);
      h_gammafakeMETvsMHT->Fill(met,mHT);
      h_gammafakeMHT_vs_HT->Fill(mHT, HT);
      h_gammafakeMet_vs_HT->Fill(met,HT);
      h_gammafakeMetSig->Fill(event->pfMETmEtSig_);
      h_gammafakeMetSigVsHT->Fill(HT,event->pfMETmEtSig_);
      h_gammafakeMetSigVsST->Fill(ST,event->pfMETmEtSig_);
      h_gammafakeMetSigVsMET->Fill(met,event->pfMETmEtSig_);
      if(numCleanJets == 0) h_gammafakeMetSigVsST_0Jets->Fill(ST,event->pfMETmEtSig_);
      if(numCleanJets == 1) h_gammafakeMetSigVsST_1Jets->Fill(ST,event->pfMETmEtSig_);
      if(numCleanJets == 2) h_gammafakeMetSigVsST_2Jets->Fill(ST,event->pfMETmEtSig_);
      if(numCleanJets >= 3) h_gammafakeMetSigVsST_3pJets->Fill(ST,event->pfMETmEtSig_);

      for(int i =0; i < nHTJets; i++){
	h_gammafakeJetEnergy->Fill(cleanJets[i].Pt() );
      }
      h_gammafakeHighestJetEnergy->Fill(highestJetEnergy);
      h_nHTJets_gammafake->Fill(nHTJets);
                
      if(printLevel>0)cout<<"Inside gammafake"<<endl;
      h_gammafakePt->Fill(PhoOneVec.Pt());
      h_gammafakePt->Fill(PhoTwoVec.Pt());
      h_gammafakePtLead->Fill(PhoOneVec.Pt()); 
      h_gammafakePtTrail->Fill(PhoTwoVec.Pt());

      h_gammafakeEta->Fill(PhoOneVec.Eta());
      h_gammafakeEta->Fill(PhoTwoVec.Eta());
      h_gammafakePhi->Fill(PhoOneVec.Phi());
      h_gammafakePhi->Fill(PhoTwoVec.Phi());
                
      float gammafakeEta = abs(event->phoEta_->at(PhoOne));
      float chIsoEA= 0;
      if(gammafakeEta<1.0){
	chIsoEA = .00;
      }
      else if (gammafakeEta <1.479 && gammafakeEta >=1.0){
	chIsoEA = .00;
      }
      float gammafakeChargedHadronIso= event->phoPFChIso_->at(PhoOne) - Rho*chIsoEA > 0. ? event->phoPFChIso_->at(PhoOne) - Rho*chIsoEA : 0.00 ;

      h_gfSigIetaIeta->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoOne));
      h_gfSigIetaIeta_vs_ChHadIso->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoOne), gammafakeChargedHadronIso);
      h_gfChHadIso->Fill(gammafakeChargedHadronIso);

      gammafakeEta = abs(event->phoEta_->at(PhoTwo));
      if(gammafakeEta<1.0){
	chIsoEA = .00;
      }
      else if (gammafakeEta <1.479 && gammafakeEta >=1.0){
	chIsoEA = .00;
      }
      gammafakeChargedHadronIso= event->phoPFChIso_->at(PhoTwo) - Rho*chIsoEA > 0. ? event->phoPFChIso_->at(PhoTwo) - Rho*chIsoEA : 0.00 ;

      h_gfSigIetaIeta->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoTwo));
      h_gfSigIetaIeta_vs_ChHadIso->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoTwo), gammafakeChargedHadronIso);
      h_gfChHadIso->Fill(gammafakeChargedHadronIso);
             
      h_gammafakeSumEt->Fill(event->pfMETsumEt_);
      

      //now do gf case specifically
      if( dogf ){
	int PhoOneGen;
	bool goForGen1 = false;
	float drClose1=999.;
	for(int i =0; i < event->nMC_; i++){
	  TLorentzVector part_it_vec = MassLorentzVector(event->mcPt->at(i), event->mcEta->at(i), event->mcPhi->at(i), event->mcMass->at(i));
	  //Try to find the gen level particles that match the photons most closely. The matching particle will have the smallest dR between it and the photon, and the difference between it's pt and the photon's pt must be less than 15%
	  if(event->mcStatus->at(i)==1){
	    if(getDR(PhoOneVec, part_it_vec)<drClose1 && abs(event->phoCalibEt_->at(PhoOne) - event->mcPt->at(i))/(event->phoCalibEt_->at(PhoOne)) < 0.15){
	      drClose1 = getDR(PhoOneVec, part_it_vec);
	      PhoOneGen=i;
	      goForGen1=true;
	    }
	  }
	}
	
	h_gfMet->Fill(met);

      }//end of if( dogf )
                
      //Now do fg case
      else if(dofg){ //lines 7166 to 7403 of Dave's code
	int PhoTwoGen;
	bool goForGen2 = false;
	float drClose2=999.;
	for(int i =0; i < event->nMC_; i++){
	  TLorentzVector part_it_vec = MassLorentzVector(event->mcPt->at(i), event->mcEta->at(i), event->mcPhi->at(i), event->mcMass->at(i));
	  //Try to find the gen level particles that match the photons most closely.The matching particle will have the smallest dR between it and the photon, and the difference between it's pt and the photon's pt must be less than 15%
	  if(event->mcStatus->at(i)==1){                            
	    if(getDR(PhoTwoVec, part_it_vec)<drClose2 && abs(event->phoCalibEt_->at(PhoTwo) - event->mcPt->at(i))/(event->phoCalibEt_->at(PhoTwo)) < 0.15){
	      drClose2 = getDR(PhoTwoVec, part_it_vec);
	      PhoTwoGen=i;
	      goForGen2=true;
	    }
	  }
	}//end of loop over MC particles
	h_fgMet->Fill(met);
      }//end of if(dofg)
      }          
    }//end dogammafake if-statement
        

    if(doef && InvMass > 75 && InvMass <  105 && passeeTrigger){
      int fake = event->phohasPixelSeed_->at(PhoOne) ? PhoTwo : PhoOne;
      bool lowR9  = event->phoR9_->at(fake) < 0.9;
      if(lowR9) h_efMet->Fill(met);
      else h_efMet_HighR9->Fill(met);
      h_efMetSigVsST->Fill(ST, event->pfMETmEtSig_);


    }

    if(doffmix && InvMass > 105 && passggTrigger){
      h_ffMet_MixR9->Fill(met);
      isFF_anyR9 = true;
      leadPt = PhoOneVec.Pt();
      trailPt = PhoTwoVec.Pt();
      leadR9 = event->phoR9_->at(PhoOne);
      trailR9 = event->phoR9_->at(PhoTwo);
    }

    if(doffhigh && InvMass > 105 && passggTrigger){
      h_ffMet_HighR9->Fill(met);
      isFF_anyR9 = true;
      leadPt = PhoOneVec.Pt();
      trailPt = PhoTwoVec.Pt();
      leadR9 = event->phoR9_->at(PhoOne);
      trailR9 = event->phoR9_->at(PhoTwo);
    }


    if(doff && InvMass > 105 && passggTrigger ){//&& HT < 200){
      //      cout << "ff Event!  Event Number " << event->event_ << ", run number " << event->run_ << ", lumi " << event->lumis_ << endl; 
      if( met < 100 && met > 85)  myfile << event->event_ << "\n";
      h_ffMet->Fill(met, PUweight);
     crossCheckError->Fill(kFALSE,met);
     h_ffNJet->Fill(numCleanJets);
     h_ffMetMinusDiEM->Fill(mmdiempt);
     isFF = true;
     isFF_anyR9 = true;
     leadPt = PhoOneVec.Pt();
     trailPt = PhoTwoVec.Pt();
     leadR9 = event->phoR9_->at(PhoOne);
     trailR9 = event->phoR9_->at(PhoTwo);



     if(!event->isData_){
       for(int Part_it = 0; Part_it < event->nMC_; Part_it++){
	 float partPdgID = event->mcPID->at(Part_it);
	 if( fabs(partPdgID)==22 && fabs(event->mcMomPID->at(Part_it)) == 1000023){ //10023 = neutralino                                                     
	   TLorentzVector part_it_vec = MassLorentzVector(event->mcPt->at(Part_it), event->mcEta->at(Part_it), event->mcPhi->at(Part_it), event->mcMass->at(Part_it)); 
	   if(isSameObject(part_it_vec,PhoOneVec,0.01) && abs(part_it_vec.Pt() -PhoOneVec.Pt()) / PhoOneVec.Pt() < 0.05){
	     h_ffDiffVtx_y->Fill(abs(event->mcVty->at(Part_it) - event->vty_) );
	     h_ffDiffVtx_z->Fill(abs(event->mcVtz->at(Part_it)- event->vtz_) );
	     h_ffDiffVtx_x->Fill(abs(event->mcVtx->at(Part_it)- event->vtx_) );
	   }
	 }
       }
     }

     if(doSignal){
       if(diEMpt > 400) h_ffDiEMPt400_Grid->Fill(gmass,nmass); 
       if(diEMpt > 300) h_ffDiEMPt300_Grid->Fill(gmass,nmass);
       if(met > 150) h_ffMET150_Grid->Fill(gmass,nmass);
       if(gmass == 2000){
	 if(nmass == 1900){ h_ffDiEMPt_2000_1900->Fill(diEMpt); }
	 if(nmass == 1500){ h_ffDiEMPt_2000_1500->Fill(diEMpt); }
	 if(nmass == 1000){ h_ffDiEMPt_2000_1000->Fill(diEMpt); }
	 if(nmass == 500){  h_ffDiEMPt_2000_500->Fill(diEMpt); }
	 if(nmass == 100){  h_ffDiEMPt_2000_100->Fill(diEMpt); }
       }
       else if(gmass == 1800){
	 if(nmass == 1500){ h_ffDiEMPt_1800_1500->Fill(diEMpt); }
	 if(nmass == 1000){ h_ffDiEMPt_1800_1000->Fill(diEMpt); }
	 if(nmass == 500){  h_ffDiEMPt_1800_500->Fill(diEMpt); }
	 if(nmass == 100){  h_ffDiEMPt_1800_100->Fill(diEMpt); }
       }
       else if(gmass == 1700){
	 if(nmass == 1500){  h_ffDiEMPt_1700_1500->Fill(diEMpt); }
	 if(nmass == 1000){  h_ffDiEMPt_1700_1000->Fill(diEMpt); }
	 if(nmass == 500){ h_ffDiEMPt_1700_500->Fill(diEMpt); }
	 if(nmass == 100){ h_ffDiEMPt_1700_100->Fill(diEMpt); }
       }
       else if(gmass == 1600){
	 if(nmass == 1500){  h_ffDiEMPt_1600_1500->Fill(diEMpt); }
	 if(nmass == 1000){  h_ffDiEMPt_1600_1000->Fill(diEMpt); }
	 if(nmass == 500){   h_ffDiEMPt_1600_500->Fill(diEMpt); }
	 if(nmass == 100){   h_ffDiEMPt_1600_100->Fill(diEMpt); }
       }
       else if(gmass == 1500){
	 if(nmass == 1400){ h_ffDiEMPt_1500_1400->Fill(diEMpt); }
	 if(nmass == 1000){ h_ffDiEMPt_1500_1000->Fill(diEMpt); }
	 if(nmass == 500){  h_ffDiEMPt_1500_500->Fill(diEMpt); }
	 if(nmass == 100){  h_ffDiEMPt_1500_100->Fill(diEMpt); }
       }
     }

      //Study C
      if(P >= 0){
	if(P > 3 && ST > 400) h_ffMet_MetSig3ST400->Fill(met);
	h_ffMet_vs_diPT->Fill(mHT, diPT);
	h_ffHT_vs_diPT->Fill(HT, diPT);
	/*	if( P >= 0){
	  if(ST>= 80 && ST<100) h_ff_fracSTEM_ST80100->Fill(fracSTEM);
	  if(ST>= 100 && ST<120) h_ff_fracSTEM_ST100120->Fill(fracSTEM);
	  if(ST>= 120 && ST<150) h_ff_fracSTEM_ST120150->Fill(fracSTEM);
	  if(ST>= 150 && ST<200) h_ff_fracSTEM_ST150200->Fill(fracSTEM);
	  if(ST>= 200 ) h_ff_fracSTEM_ST200->Fill(fracSTEM);
	  h_ffMetSigVsfracSTEM->Fill(fracSTEM, P);
	  h_ffSTVsfracSTEM->Fill(fracSTEM, ST);
	  }*/
      }



      nff++;
      float diempt_weight=1; //The weights are stored in the histogram ffweights
      int bin = ffweights->FindFixBin(float(diEMpt));
      diempt_weight=ffweights->GetBinContent(bin);


      if(met > 100) {
	nff_Met100++;      
      //      if(met > 150) cout << "ff event with MET > 150! Event number "<< event->event_ << ", lumi "<< event->lumis_ << " and run number "<<event->run_ << endl; 
      //Output to help with synchronization and to check specific events
	float chIsoEA=0, neuIsoEA=0, phoIsoEA=0;
	double phoEta = std::abs(event->phoSCEta_->at(PhoOne));
	if(phoEta<1.0){
	  chIsoEA  = 0.0360;
	  neuIsoEA = 0.0597;
	  phoIsoEA = 0.1210;
	}
	else if (phoEta <1.479 && phoEta >=1.0){
	  chIsoEA  = 0.0377;
	  neuIsoEA = 0.0807;
	  phoIsoEA = 0.1107;
	}

	float PhoEt= event->phoCalibEt_->at(PhoOne);
	float uncorrNeuHadIso = event->phoPFNeuIso_->at(PhoOne);
	float uncorrChHadIso = event->phoPFChIso_->at(PhoOne);
	float uncorrPhoIso = event->phoPFPhoIso_->at(PhoOne);

	//These are the pile-up corrected isolations
	float chargedHadronIso1= uncorrChHadIso - Rho*chIsoEA > 0. ? uncorrChHadIso - Rho*chIsoEA : 0.00 ;
	float neutralHadronIso1= ( uncorrNeuHadIso - Rho*neuIsoEA - 0.0148*PhoEt - 0.000017*PhoEt*PhoEt ) > 0.00 ? uncorrNeuHadIso - Rho*neuIsoEA - 0.0148*PhoEt - 0.000017*PhoEt*PhoEt: 0.00;
	float photonIso1 = ( uncorrPhoIso - Rho*phoIsoEA - 0.0047*PhoEt ) > 0.00 ? uncorrPhoIso - Rho*phoIsoEA - 0.0047*PhoEt : 0.00;

	phoEta = std::abs(event->phoSCEta_->at(PhoTwo));
	if(phoEta<1.0){
	  chIsoEA  = 0.0360;
	  neuIsoEA = 0.0597;
	  phoIsoEA = 0.1210;
	}
	else if (phoEta <1.479 && phoEta >=1.0){
	  chIsoEA  = 0.0377;
	  neuIsoEA = 0.0807;
	  phoIsoEA = 0.1107;
	}

        PhoEt= event->phoCalibEt_->at(PhoTwo);
	uncorrNeuHadIso = event->phoPFNeuIso_->at(PhoTwo);
	uncorrChHadIso = event->phoPFChIso_->at(PhoTwo);
	uncorrPhoIso = event->phoPFPhoIso_->at(PhoTwo);

	//These are the pile-up corrected isolations                                                                                  
	float chargedHadronIso2= uncorrChHadIso - Rho*chIsoEA > 0. ? uncorrChHadIso - Rho*chIsoEA : 0.00 ;
	float neutralHadronIso2= ( uncorrNeuHadIso - Rho*neuIsoEA - 0.0148*PhoEt - 0.000017*PhoEt*PhoEt ) > 0.00 ? uncorrNeuHadIso - Rho*neuIsoEA - 0.0148*PhoEt - 0.000017*PhoEt*PhoEt: 0.00;
	float photonIso2 = ( uncorrPhoIso - Rho*phoIsoEA - 0.0047*PhoEt ) > 0.00 ? uncorrPhoIso - Rho*phoIsoEA - 0.0047*PhoEt : 0.00;

	//	if(diEMpt > 300 ) myfile << met << " " << diEMpt << " " << diempt_weight << endl;
	//        myfile << met << " " << diEMpt << " " << diempt_weight << endl;
	//" " << event->phoCalibEt_->at(PhoOne) << " " << event->phoCalibEt_->at(PhoTwo) << " " << event->phoEta_->at(PhoOne) << " " << event->phoEta_->at(PhoTwo) << " " << event->phoPhi_->at(PhoOne) << " " << event->phoPhi_->at(PhoTwo) << " " << event->phoR9_->at(PhoOne) << " " << event->phoR9_->at(PhoTwo) << " " << event->phoHoverE_->at(PhoOne) << " " << event->phoHoverE_->at(PhoTwo) << " " << event-> phoSigmaIEtaIEtaFull5x5_->at(PhoOne) << " " <<  event-> phoSigmaIEtaIEtaFull5x5_->at(PhoTwo) << " " << chargedHadronIso1 << " " << chargedHadronIso2 << " " << neutralHadronIso1 << " " << neutralHadronIso2 << " " << photonIso1 << " " << photonIso2 <<  endl;
	/*	cout << "ff event with MET > 140! Event number "<< event->event_ << ", lumi "<< event->lumis_ << " and run number "<<event->run_ << endl;
	cout << "Pho One Index = " << PhoOne << " and PhoTwo Index= " << PhoTwo << endl;
	cout << "Pho One Pt, Eta, Phi =" << event-> phoCalibEt_->at(PhoOne) << ", " << event->phoEta_->at(PhoOne) << ", " << event->phoPhi_->at(PhoOne) << endl;
	cout <<"Pho Two Pt, Eta, Phi =" << event-> phoCalibEt_->at(PhoTwo) << ", " << event->phoEta_->at(PhoTwo) <<", " <<event->phoPhi_->at(PhoTwo) << endl;
	cout << "Bit 14 of HLT Pho = " << ( event->HLTPho_>>triggerIndex&1) << endl;
	cout << "Met filters = " << event->metFilters_ << endl;
	cout << "Number of vertices = " << event->nVtx_ << endl;
	cout << "Met = " << met << endl;
	cout << " PhoOne sig ieta ieta = " << event->phoSigmaIEtaIEta_->at(PhoOne) << endl;
	cout << " PhoTwo sig ieta ieta = " << event->phoSigmaIEtaIEta_->at(PhoTwo) << endl;
	cout << " PhoOne charged hadron isolation = " << event->phoPFChIso_->at(PhoOne) << endl;
	cout << " PhoTwo charged hadron isolation = " << event->phoPFChIso_->at(PhoTwo) << endl;
	cout << "Invariant mass = " << InvMass << endl;*/
      }

      //      if(met < 100 && met > 85) myfile<< "Event number " << event->event_ << ", diEMpT = " << diEMpt << " and weight << " << diempt_weight << ". MET = " << met << endl;
      //      if(met > 75 && met < 85) myfile<< event->event_ <<  "\n";
      if(met > 250) cout << "FF EVENT WITH MET > 250 GEV! EVENT NUMBER = " << event->event_ << endl;
      h_ffMetUnweighted->Fill(met);
      h_ffDiEMPt->Fill(diEMpt);
      h_ffDiEMPt_vs_NJet->Fill(diEMpt,numCleanJets);
      h_ffDiEMPt_vs_HT->Fill(diEMpt,HT);
      h_ffNJet_Reweighted->Fill(numCleanJets,diempt_weight);
      h_ffMetMinusDiEM_Reweighted->Fill(mmdiempt,diempt_weight);
      float NJet_weight = 1;
      bin = h2D_weights_ff->FindFixBin(float(diEMpt), float(numCleanJets));
      NJet_weight = h2D_weights_ff->GetBinContent(bin);
      h_ffMetReweighted_NJet->Fill(met,NJet_weight);
      if(doSignal && met > 150) h_ffRWMET150_Grid->Fill(gmass,nmass,diempt_weight);

      //fill toy histograms
      if(doToys){
	toys->cd();
	float toy_diempt_weight=1; //The weights are stored in the histograms in RatioGausFF
	bin = ffweights->FindFixBin(float(diEMpt));
	for(int k=0; k<1000; k++){
	  toy_diempt_weight = RatioGausFF[k]->GetBinContent(bin);
	  ReweightedFF[k]->Fill(met,toy_diempt_weight);
	}
	fout->cd();
      }
      //Done with toy histograms    

      h_ffMetReweighted->Fill(met,diempt_weight);
      h_ffHTReweighted->Fill(HT,diempt_weight);
      h_ffMHTReweighted->Fill(mHT, diempt_weight);
      h_ffMet_vs_HT->Fill(met,HT);
      h_ffMet_vs_ST->Fill(met,ST);

      h_ffR9->Fill(event->phoR9_->at(PhoOne) ) ;
      h_ffR9->Fill(event->phoR9_->at(PhoTwo) ) ;

      h_ffHT->Fill(HT);
      h_ffST->Fill(ST);
      h_ffMHT->Fill(mHT);
      h_ff_JetMETdPhi->Fill(minJetMETdphi);
      h_ff_JetMETdPhi_vs_Met->Fill(met,minJetMETdphi);
      if(HT > 200){
        h_ff_JetMETdPhi_vs_Met_ST200->Fill(met,minJetMETdphi);
      }
      if(closeJetMET){
        h_ffMETvsMHT_CloseJetMET->Fill(met,mHT);
      }
      h_ffMETvsMHT->Fill(met,mHT);
      h_ffMHT_vs_HT->Fill(mHT, HT);
      if(P >= 0){
	h_ffMetSig->Fill(P);
	//      if(met>50){
	h_ffMetSigVsMET->Fill(met,P);
	h_ffMetSigVsST->Fill(ST,P);
	//	if(fracSTEM <0.7) h_ffMetSigVsST_lowFracSTEM->Fill(ST, P);
      }
      h_ffSumEtVsMET->Fill(met, event->pfMETsumEt_);
      h_ffSumEtVsST->Fill(ST, event->pfMETsumEt_);
      h_ffSumEt->Fill(event->pfMETsumEt_);

      h_ffMT2->Fill(MT2);
      h_ffMT2_vs_ST->Fill(MT2,HT);
      h_ffMT2_vs_MET->Fill(MT2,met);
      if(P >= 0){   
	h_ffMT2_vs_MetSig->Fill(MT2,P);
	h_ffPVsST->Fill(ST, P);
	if(numCleanJets == 1){
	  h_ffMetSigVsST_1Jets->Fill(ST,P);
	  h_ffMetSigVsHT_1Jets->Fill(HT,P);
	}
	else if(numCleanJets == 2){
	  h_ffMetSigVsST_2Jets->Fill(ST,P);
	  h_ffMetSigVsHT_2Jets->Fill(HT,P);
	}
	else if(numCleanJets >= 3){
	  h_ffMetSigVsST_3pJets->Fill(ST,P);
	  h_ffMetSigVsHT_3pJets->Fill(HT,P);
	  if(numCleanJets == 3) h_ffMetSigVsST_3Jets->Fill(ST,P);
	  if(numCleanJets== 4) h_ffMetSigVsST_4Jets->Fill(ST,P);
	  if(numCleanJets== 5) h_ffMetSigVsST_5Jets->Fill(ST,P);
	  if(numCleanJets>= 6) h_ffMetSigVsST_6pJets->Fill(ST,P);
	}
	h_ffMetSigVsHT->Fill(HT,P);
      }
      for(int i =0; i < nHTJets; i++){
	h_ffJetEnergy->Fill(cleanJets[i].Pt() );
      }
      h_ffHighestJetEnergy->Fill(highestJetEnergy);
      h_nHTJets_ff->Fill(nHTJets);

      bool leadIsoFake, leadShapeFake, trailIsoFake, trailShapeFake;
      float ffEta = abs(event->phoEta_->at(PhoOne));
      float chIsoEA=0.0;
      if(ffEta<1.0){
	chIsoEA = 0.0360;
      }
      else if (ffEta <1.479 && ffEta >=1.0){
        chIsoEA = 0.0377;
      }
      float ffChargedHadronIsoLead = event->phoPFChIso_->at(PhoOne) - Rho*chIsoEA > 0. ? event->phoPFChIso_->at(PhoOne) - Rho*chIsoEA : 0.00 ;

      if(ffChargedHadronIsoLead > 0.441){ leadIsoFake = true; nFFIso++;}
      else{ leadIsoFake = false;}
      if(event->phoSigmaIEtaIEtaFull5x5_->at(PhoOne) > 0.01022){ leadShapeFake = true; nFFShape++;}
      else {leadShapeFake = false;}

      h_ffSigIetaIeta->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoOne));
      h_ffSigIetaIeta_vs_ChHadIso->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoOne), ffChargedHadronIsoLead);
      if(met > 100)    h_ffSigIetaIeta_vs_ChHadIso_Met100->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoOne), ffChargedHadronIsoLead);
      h_ffChHadIso->Fill(ffChargedHadronIsoLead);

      ffEta = abs(event->phoEta_->at(PhoTwo));
      if(ffEta<1.0){
	chIsoEA = .0360;
      }
      else if (ffEta <1.479 && ffEta >=1.0){
	chIsoEA = .0377;
      }

      float ffChargedHadronIsoTrail= event->phoPFChIso_->at(PhoTwo) - Rho*chIsoEA > 0. ? event->phoPFChIso_->at(PhoTwo) - Rho*chIsoEA : 0.00 ;

      if(ffChargedHadronIsoTrail > 0.441){ trailIsoFake= true; nFFIso++;}
      else{ trailIsoFake = false;}
      if(event->phoSigmaIEtaIEtaFull5x5_->at(PhoTwo) > 0.01022){ trailShapeFake =true; nFFShape++;}
      else{ trailShapeFake = false;}

      h_ffSigIetaIeta->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoTwo));
      h_ffSigIetaIeta_vs_ChHadIso->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoTwo), ffChargedHadronIsoTrail);
      if(met > 100)    h_ffSigIetaIeta_vs_ChHadIso_Met100->Fill(event->phoSigmaIEtaIEtaFull5x5_->at(PhoTwo), ffChargedHadronIsoTrail);
      h_ffChHadIso->Fill(ffChargedHadronIsoTrail);

      float minDR1 = 1000., minDR2 = 1000.;
      float sminDR1 = 1000., sminDR2 = 1000.;
      float dR;
      for(unsigned int i = 0;i<cleanJets.size(); i++){
	dR = getDR(PhoOneVec, cleanJets[i]);
	if(dR < minDR1) minDR1 = dR;
	else if(dR  < sminDR1) sminDR1 = dR;

	dR = getDR(PhoTwoVec, cleanJets[i]);
	if(dR <minDR2) minDR2 = dR;
	else if(dR  < sminDR2) sminDR2 = dR;
      }

      h_ffMinDR_vs_NJet->Fill(min(minDR1,minDR2),numCleanJets);
      h_ffNJet_vs_MET->Fill(numCleanJets,met);
      h_ffNVertex_vs_MET->Fill(NVertex,met);
      h_ffRho_vs_MET->Fill(Rho,met);
      h_ffRho->Fill(Rho);

      h_ffMinDR_vs_MET->Fill(min(minDR1,minDR2),met);

      h_ffSecondMinDR_vs_NJet->Fill(min(sminDR1,sminDR2),numCleanJets);
      h_ffSecondMinDR_vs_MET->Fill(min(sminDR1,sminDR2),met);

      if(trailIsoFake && leadIsoFake){
	h_ffChHadIsoMax->Fill(max(ffChargedHadronIsoLead,ffChargedHadronIsoTrail));
	h_ffChHadIsoMin->Fill(min(ffChargedHadronIsoLead,ffChargedHadronIsoTrail));
	h_ffChHadIsoMin_vs_minDR->Fill(min(ffChargedHadronIsoLead,ffChargedHadronIsoTrail),
                                       min(minDR1,minDR2) );
	h_ffChHadIsoMin_vs_MET->Fill(min(ffChargedHadronIsoLead,ffChargedHadronIsoTrail), met);
      }
      else if(trailIsoFake){
        h_ffChHadIsoMax->Fill(ffChargedHadronIsoTrail);
        h_ffChHadIsoMin->Fill(ffChargedHadronIsoTrail);
	h_ffChHadIsoMin_vs_minDR->Fill(ffChargedHadronIsoTrail,
                                       min(minDR1,minDR2) );
	h_ffChHadIsoMin_vs_MET->Fill(ffChargedHadronIsoTrail, met);
      }
      else if(leadIsoFake){
        h_ffChHadIsoMax->Fill(ffChargedHadronIsoLead);
        h_ffChHadIsoMin->Fill(ffChargedHadronIsoLead);
        h_ffChHadIsoMin_vs_minDR->Fill(ffChargedHadronIsoLead,
                                       min(minDR1,minDR2) );
	h_ffChHadIsoMin_vs_MET->Fill(ffChargedHadronIsoLead, met);
      }
      h_ffSCE->Fill(event->phoSCE_->at(PhoOne));
      h_ffSCE->Fill(event->phoSCE_->at(PhoTwo));
      h_ffRawSCE->Fill(event->phoSCRawE_->at(PhoOne));
      h_ffRawSCE->Fill(event->phoSCRawE_->at(PhoTwo));
      h_ffSCE_diffRaw->Fill( abs(event->phoSCE_->at(PhoOne) - event->phoSCRawE_->at(PhoOne) ));
      h_ffSCE_diffRaw->Fill( abs(event->phoSCE_->at(PhoTwo) - event->phoSCRawE_->at(PhoTwo) ));
      h_ffSCE_normalRaw->Fill( event->phoSCE_->at(PhoOne) / event->phoSCRawE_->at(PhoOne) );
      h_ffSCE_normalRaw->Fill( event->phoSCE_->at(PhoTwo) / event->phoSCRawE_->at(PhoTwo) );

      h_ffPt->Fill(PhoOneVec.Pt());
      h_ffPt->Fill(PhoTwoVec.Pt());

    }//end of fake-fake if-statement
        
    //Fill tree and clean up for next round
    NJets = numCleanJets;
    if(isFF || isGG || isFF_anyR9) {
      tree->Fill();
    }

  }//end of looping over the events

  if(doSignal){
    //Finish up acceptance*efficiency histograms now that we are done looping through the events
    for(int i = 1; i <= xbinnum3;i++){
      for(int j = 1; j <= ybinnum;j++){
	float mGlu = xbins3[i] - 25;
	float mNeu = (ybins[j] + ybins[j-1]) / 2.0;
	if(mNeu >= mGlu) continue;
	float xsec = 0;
	try{
	  xsec = xsec_map[mGlu];
	}
	catch(...){
	  cout<<"Uh oh" << endl;
	  break;
	}

	if(h_GridAllGGEvents->GetBinContent(i,j) > 0 ){
	  h_Grid_PhoSF_Err_Average->SetBinContent( i, j, ( h_Grid_PhoSF_Err->GetBinContent(i,j) / h_GridAllGGEvents2->GetBinContent(i,j))   );
	  h_Grid_PhoSF_Average->SetBinContent( i, j, ( h_Grid_PhoSF->GetBinContent(i,j) / h_GridAllGGEvents2->GetBinContent(i,j)) );
	}
	if( h_GridAllEvents->GetBinContent(i,j) > 0 ){
	  float scaleByLumi = lumi*xsec/h_GridAllEvents->GetBinContent(i,j);
	  cout << "scale factor for " << mGlu <<", " << mNeu << " is " << scaleByLumi << endl;
	  h_MET100to115_Grid->SetBinContent(i,j,h_MET100to115_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
	  h_MET115to130_Grid->SetBinContent(i,j,h_MET115to130_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
	  h_MET130to150_Grid->SetBinContent(i,j,h_MET130to150_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
	  h_MET150to185_Grid->SetBinContent(i,j,h_MET150to185_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
	  h_MET185to250_Grid->SetBinContent(i,j,h_MET185to250_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
	  h_MET250_Grid->SetBinContent(i,j,h_MET250_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
          h_MET150_Grid->SetBinContent(i,j,h_MET150_GridUnweighted->GetBinContent(i,j)*scaleByLumi);

	  h_MET100to115_Grid_JESDown->SetBinContent(i,j,h_MET100to115_GridUnweighted_JESDown->GetBinContent(i,j)*scaleByLumi);
	  h_MET115to130_Grid_JESDown->SetBinContent(i,j,h_MET115to130_GridUnweighted_JESDown->GetBinContent(i,j)*scaleByLumi);
	  h_MET130to150_Grid_JESDown->SetBinContent(i,j,h_MET130to150_GridUnweighted_JESDown->GetBinContent(i,j)*scaleByLumi);
	  h_MET150to185_Grid_JESDown->SetBinContent(i,j,h_MET150to185_GridUnweighted_JESDown->GetBinContent(i,j)*scaleByLumi);
	  h_MET185to250_Grid_JESDown->SetBinContent(i,j,h_MET185to250_GridUnweighted_JESDown->GetBinContent(i,j)*scaleByLumi);
	  h_MET250_Grid_JESDown->SetBinContent(i,j,h_MET250_GridUnweighted_JESDown->GetBinContent(i,j)*scaleByLumi);
          h_MET150_Grid_JESDown->SetBinContent(i,j,h_MET150_GridUnweighted_JESDown->GetBinContent(i,j)*scaleByLumi);

	  h_MET100to115_Grid_JESUp->SetBinContent(i,j,h_MET100to115_GridUnweighted_JESUp->GetBinContent(i,j)*scaleByLumi);
	  h_MET115to130_Grid_JESUp->SetBinContent(i,j,h_MET115to130_GridUnweighted_JESUp->GetBinContent(i,j)*scaleByLumi);
	  h_MET130to150_Grid_JESUp->SetBinContent(i,j,h_MET130to150_GridUnweighted_JESUp->GetBinContent(i,j)*scaleByLumi);
	  h_MET150to185_Grid_JESUp->SetBinContent(i,j,h_MET150to185_GridUnweighted_JESUp->GetBinContent(i,j)*scaleByLumi);
	  h_MET185to250_Grid_JESUp->SetBinContent(i,j,h_MET185to250_GridUnweighted_JESUp->GetBinContent(i,j)*scaleByLumi);
	  h_MET250_Grid_JESUp->SetBinContent(i,j,h_MET250_GridUnweighted_JESUp->GetBinContent(i,j)*scaleByLumi);
          h_MET150_Grid_JESUp->SetBinContent(i,j,h_MET150_GridUnweighted_JESUp->GetBinContent(i,j)*scaleByLumi);

	}

	h_MET100to115_GridStatError->Fill(mGlu,mNeu,h_MET100to115_GridUnweighted->GetBinError(i,j)/h_MET100to115_GridUnweighted->GetBinContent(i,j));
	h_MET115to130_GridStatError->Fill(mGlu,mNeu,h_MET115to130_GridUnweighted->GetBinError(i,j)/h_MET115to130_GridUnweighted->GetBinContent(i,j));
	h_MET130to150_GridStatError->Fill(mGlu,mNeu,h_MET130to150_GridUnweighted->GetBinError(i,j)/h_MET130to150_GridUnweighted->GetBinContent(i,j));
	h_MET150to185_GridStatError->Fill(mGlu,mNeu,h_MET150to185_GridUnweighted->GetBinError(i,j)/h_MET150to185_GridUnweighted->GetBinContent(i,j));
	h_MET185to250_GridStatError->Fill(mGlu,mNeu,h_MET185to250_GridUnweighted->GetBinError(i,j)/h_MET185to250_GridUnweighted->GetBinContent(i,j));
	h_MET250_GridStatError->Fill(mGlu,mNeu,h_MET250_GridUnweighted->GetBinError(i,j)/h_MET250_GridUnweighted->GetBinContent(i,j));
        h_MET150_GridStatError->Fill(mGlu,mNeu,h_MET150_GridUnweighted->GetBinError(i,j)/h_MET150_GridUnweighted->GetBinContent(i,j));

	h_MET100to115_GridJESUpError->Fill(mGlu,mNeu, abs( h_MET100to115_Grid_JESUp->GetBinContent(i,j) - h_MET100to115_Grid->GetBinContent(i,j)) / h_MET100to115_Grid->GetBinContent(i,j) );
	h_MET115to130_GridJESUpError->Fill(mGlu,mNeu, abs( h_MET115to130_Grid_JESUp->GetBinContent(i,j) - h_MET115to130_Grid->GetBinContent(i,j)) / h_MET115to130_Grid->GetBinContent(i,j) );
	h_MET130to150_GridJESUpError->Fill(mGlu,mNeu, abs( h_MET130to150_Grid_JESUp->GetBinContent(i,j) - h_MET130to150_Grid->GetBinContent(i,j)) / h_MET130to150_Grid->GetBinContent(i,j) );
	h_MET150to185_GridJESUpError->Fill(mGlu,mNeu, abs( h_MET150to185_Grid_JESUp->GetBinContent(i,j) - h_MET150to185_Grid->GetBinContent(i,j)) / h_MET150to185_Grid->GetBinContent(i,j) );
	h_MET185to250_GridJESUpError->Fill(mGlu,mNeu, abs( h_MET185to250_Grid_JESUp->GetBinContent(i,j) - h_MET185to250_Grid->GetBinContent(i,j)) / h_MET185to250_Grid->GetBinContent(i,j) );
	h_MET250_GridJESUpError->Fill(mGlu,mNeu, abs( h_MET250_Grid_JESUp->GetBinContent(i,j) - h_MET250_Grid->GetBinContent(i,j)) / h_MET250_Grid->GetBinContent(i,j) );
        h_MET150_GridJESUpError->Fill(mGlu,mNeu, abs( h_MET150_Grid_JESUp->GetBinContent(i,j) - h_MET150_Grid->GetBinContent(i,j)) / h_MET150_Grid->GetBinContent(i,j) );

	h_MET100to115_GridJESDownError->Fill(mGlu,mNeu, abs( h_MET100to115_Grid_JESDown->GetBinContent(i,j) - h_MET100to115_Grid->GetBinContent(i,j)) / h_MET100to115_Grid->GetBinContent(i,j) );
	h_MET115to130_GridJESDownError->Fill(mGlu,mNeu, abs( h_MET115to130_Grid_JESDown->GetBinContent(i,j) - h_MET115to130_Grid->GetBinContent(i,j)) / h_MET115to130_Grid->GetBinContent(i,j) );
	h_MET130to150_GridJESDownError->Fill(mGlu,mNeu, abs( h_MET130to150_Grid_JESDown->GetBinContent(i,j) - h_MET130to150_Grid->GetBinContent(i,j)) / h_MET130to150_Grid->GetBinContent(i,j) );
	h_MET150to185_GridJESDownError->Fill(mGlu,mNeu, abs( h_MET150to185_Grid_JESDown->GetBinContent(i,j) - h_MET150to185_Grid->GetBinContent(i,j)) / h_MET150to185_Grid->GetBinContent(i,j) );
	h_MET185to250_GridJESDownError->Fill(mGlu,mNeu, abs( h_MET185to250_Grid_JESDown->GetBinContent(i,j) - h_MET185to250_Grid->GetBinContent(i,j)) / h_MET185to250_Grid->GetBinContent(i,j) );
	h_MET250_GridJESDownError->Fill(mGlu,mNeu, abs( h_MET250_Grid_JESDown->GetBinContent(i,j) - h_MET250_Grid->GetBinContent(i,j)) / h_MET250_Grid->GetBinContent(i,j) );
        h_MET150_GridJESDownError->Fill(mGlu,mNeu, abs( h_MET150_Grid_JESDown->GetBinContent(i,j) - h_MET150_Grid->GetBinContent(i,j)) / h_MET150_Grid->GetBinContent(i,j) );

      }
    }
  }

  if(nCnt[3] == 0) cout << " No events passed trigger!" << endl;
  cout<< nCnt[5] << " events passed JSON file, trigger, and nvertex and met filter cuts." << endl; 
  float eff = 0.0;// temp2/temp;
  cout << "Efficiency of loose ID is " << temp2 << "/" << temp << " = " << eff << endl; 
  //  eff = numPassGenPho/numTotalGenPho;
  cout << "Efficiency with new object cleaning is " << numPassGenPho << "/" << numTotalGenPho << endl;
  cout<< "Total number of photons passing all selection cuts is " << nPhoCandsTotal<<endl;
  cout<< "Total number of fake candidates is " << nFakeCandsTotal << endl;
  cout << "We have " << ngg << " gamma gamma events, " << nff << " fake fake events ("<< nFFIso << " iso and " << nFFShape << " shape), and " << numGammaFake << " gamma fake events." << endl;
  cout<<"Of the "<< ngg << " gg events, there are "<< nggMET0to50 << " events with MET < 50 GeV and "<<nggMET50to100 << " events with MET between 50 and 100 GeV." <<endl;
  cout<< "There are "<< nggInvMass110 << " gg events with invariant mass above 110 GeV." << endl;
  cout << "We have " << nff_Met100 << " ff events with MET > 100" << endl;

  cout << "Considered " << nCnt[2] << " events. "<<nCnt[3] << " passed trigger." << endl;

  cout << "Num gg = " << ngg << endl;  
  cout << "Num ff = " << nff <<endl;
  cout << "Num ee = " << nee <<endl;
  cout << "Funky ee events = " << nee_funky << endl;
  /*
  cout << "For the egamma sample with njets >=2: " << endl;
  cout << "     100 to 150: " << numEG_100to150 << endl;
  cout << "     150 to 250: " << numEG_150to250 << endl;
  cout << "     250 and up: " << numEG_250to350 << endl;

  cout << "For the gg sample: " << endl;
  cout << "     100 to 150: " << ggMet100to150 << endl;
  cout << "     150 to 250: " << ggMet150to250 << endl;
  cout << "     250 and up: " << ggMet250 << endl;*/

  cout << "Number of egamma events in the signal region!" << endl;
  cout << "Bin 1:  " << numEGSignal[0] << endl;
  cout << "Bin 2:  " <<numEGSignal[1] << endl;
  cout << "Bin 3:  " <<numEGSignal[2] << endl;
  cout << "Bin 4:  " <<numEGSignal[3] << endl;
  cout << "Bin 5:  " <<numEGSignal[4] << endl;
  cout << "Bin 6:  " <<numEGSignal[5] << endl;
  cout << "Bin 7:  " <<numEGSignal[6] << endl;
  cout << "Bin 8:  " <<numEGSignal[7] << endl;
  cout << "Bin 9:  " <<numEGSignal[8] << endl;

  cout << "Number of events with two photon candidates is " << nTwoCands <<endl;
  cout<<"Writing analysis root output to: hist_"<<ds<<".root"<<endl;
  // close the output file
    
  fout->cd();
  fout->Write();
  fout->Close();
  delete fout;
    
    
  }//end of Loop() function





//---------------------------------------------------------------




//void SusyEventAnalyzer::CategorizeEvents(int pho1, int pho2, float Rho, bool &gogg, bool &goee, bool &goeg, bool &goff, bool &gogammafake, bool &gogf, bool &gofg){
//return;
//}



void SusyEventAnalyzer::CategorizeEvents(int pho1, int pho2, int pho3, float Rho, bool &gogg, bool &goee, bool &goeg, bool &goff, bool &gogammafake, bool &gogf, bool &gofg, bool &goef, bool &use13, bool&use23, bool& goffmix, bool& goffhigh){
    //All of the photon candidates have had to pass cuts on H over E, photon isolation, and neutral hadron isolation
    //electrons and photons must also pass sigma ieta ieta and charged hadron isolation cuts. Photons must not have a pixel seed.
    //initialize all booleans to false
    bool g1=false, g2=false, f1=false, f2=false, e1=false, e2=false;
    bool f1high = false, f2high = false;
    gogg=false;goee=false;goeg=false;goff=false;gogammafake=false;gogf=false;gofg=false;
    use23 = false, use13 = false, goffmix = false, goffhigh = false; 
    if((pho1>=event->nPho_)||(pho2>=event->nPho_)){
        cout << "Photon index higher than the number of photons in the event!! Error!!"<<endl;
        cout<< "PhoOne= " <<pho1<< "  PhoTwo= "<<pho2<<"  nPho= "<<event->nPho_<<endl;
        return;
    }
    
    float maxSihih = 0.01022; //.03; //.0106;
    float maxChHadIso = 0.441;
    
    //Define the effective area for pho1 (see https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2 )
    float pho1Eta = std::abs(event->phoSCEta_->at(pho1));
    //    float pho1Eta = abs(event->phoEta_->at(pho1) );
    float chIsoEA1=0;
    
    if(pho1Eta<1.0){
        chIsoEA1 = 0.0360;
    }
    else if (pho1Eta <1.479 && pho1Eta >=1.0){
        chIsoEA1 = 0.0377;
    }
    else if (pho1Eta <2.0 && pho1Eta >=1.479){
        chIsoEA1 = 0.0306;
    }
    else if (pho1Eta <2.2 && pho1Eta >=2.0){
        chIsoEA1 = .0283;
    }
    else if (pho1Eta <2.3 && pho1Eta >=2.2){
        chIsoEA1 = .0254;
    }
    else if (pho1Eta <2.4 && pho1Eta >=2.3){
        chIsoEA1 = .0217;
    }
    else if (pho1Eta >=2.4){
        chIsoEA1 = .0167;
    }
    
    //Define EA for pho 2
    float pho2Eta = std::abs(event->phoSCEta_->at(pho2));    
    float chIsoEA2=0.0;
    
    if(pho2Eta<1.0){
        chIsoEA2 = .0360;
    }
    else if (pho2Eta <1.479 && pho2Eta >=1.0){
        chIsoEA2 = .0377;
    }
    else if (pho2Eta <2.0 && pho2Eta >=1.479){
        chIsoEA2 = .0306;
    }
    else if (pho2Eta <2.2 && pho2Eta >=2.0){
        chIsoEA2 = .0283;
    }
    else if (pho2Eta <2.3 && pho2Eta >=2.2){
        chIsoEA2 = .0254;
    }
    else if (pho2Eta <2.4 && pho2Eta >=2.3){
        chIsoEA2 = .0217;
    }
    else if (pho2Eta >=2.4){
        chIsoEA2 = .0167;
    }
    
    //Sort pho1
    float pho1iso  = (event->phoPFChIso_->at(pho1) - Rho*chIsoEA1) > 0.00 ? (event->phoPFChIso_->at(pho1) - Rho*chIsoEA1) : 0.00;
    if(event->phoSigmaIEtaIEtaFull5x5_->at(pho1) < maxSihih && pho1iso < maxChHadIso){
        if(event->phohasPixelSeed_->at(pho1)){
            e1 = true;
        }
        else{
            g1 = true;
        }
    }
    else{  //since all the photon candidates are either electrons, photons, or fakes, if it failed the sigma ieta ieta or the charged hadron isolation cut, it must be a fake
        f1 = true;
	f1high = event->phoR9_->at(pho1) > 0.9;
    }
    
    //Sort pho2
    float pho2iso  = (event->phoPFChIso_->at(pho2) - Rho*chIsoEA2) > 0.00 ? (event->phoPFChIso_->at(pho2) - Rho*chIsoEA2) : 0.00;
    if(event->phoSigmaIEtaIEtaFull5x5_->at(pho2) < maxSihih && pho2iso < maxChHadIso){
        if(event->phohasPixelSeed_->at(pho2)){
            e2 = true;
        }
        else{
            g2 = true;
        }
    }
    else{  //since all the photon candidates are either electrons, photons, or fakes, if it failed the sigma ieta ieta or the charged hadron isolation cut, it must be a fake
        f2 = true;
	f2high = event->phoR9_->at(pho2) > 0.9;
    }
    if( (e1 && f2) || (e2 && f1))   goef=true;
    if(g1 && g2)                    gogg=true; 
    if( (g1 && e2) || (e1 && g2))   goeg=true;
    if(e1 && e2)                    goee=true; 
    if( (f1 && !f1high) && (f2&& !f2high) )   goff=true;
    if(f1 && f2){
      if(f1high){
	if(!f2high)                 goffmix = true; //1 is high and 2 is low
	else                        goffhigh = true; //both are high
      }
      else{
	if(!f2high)                 goff = true; //both are low
	else                        goffmix = true; //1 is low and 2 is high
      }
    }
    if( (g1 && f2) || (f1 && g2) ){
      gogammafake=true;
      if(g1 && f2)    gogf=true;
      else if(f1 && g2)   gofg=true;
    }

    if( ! (gogg || goee || goff || goeg) ) {
      if(pho3 > -1){
	//if the first two photons didn't fall into a useful category, use the third photon object
	float chIsoEA3 = 0.0;
	float pho3Eta = std::abs(event->phoSCEta_->at(pho3));
	if(pho3Eta<1.0){
	  chIsoEA3 = .0360;
	}
	else if (pho3Eta <1.479 && pho3Eta >=1.0){
	  chIsoEA3 = .0377;
	}
	else if (pho3Eta <2.0 && pho3Eta >=1.479){
	  chIsoEA3 = .0306;
	}
	else if (pho3Eta <2.2 && pho3Eta >=2.0){
	  chIsoEA3 = .0283;
	}
	else if (pho3Eta <2.3 && pho3Eta >=2.2){
	  chIsoEA3 = .0254;
	}
	else if (pho3Eta <2.4 && pho3Eta >=2.3){
	  chIsoEA3 = .0217;
	}
	else if (pho3Eta >=2.4){
	  chIsoEA3 = .0167;
	}

 
	float pho3iso  = (event->phoPFChIso_->at(pho3) - Rho*chIsoEA3) > 0.00 ? (event->phoPFChIso_->at(pho3) - Rho*chIsoEA3) : 0.00;
	bool g3 = false, e3 = false, f3 = false;
	if(event->phoSigmaIEtaIEtaFull5x5_->at(pho3) < maxSihih && pho3iso < maxChHadIso){
	  if(event->phohasPixelSeed_->at(pho3)){
            e3 = true;
	  }
	  else{
            g3 = true;
	  }
	}
	else{  //since all the photon candidates are either electrons, photons, or fakes, if it failed the sigma ieta ieta or the charged hadron isolation cut, it must be a fake      
	  f3 = true;
	}

	if(g1 && g3){ gogg = true; use13 = true;}
	else if( (g1 && e3) || (e1 && g3) ){ goeg = true; use13 = true;}
	else if( e1 && e3) {goee = true; use13 = true;}
	else if(f1 && f3) {goff = true; use13 = true;}
	else if( g2&&g3) {gogg = true; use23 = true;}
	else if( (g2 && e3) || (e2 && g3) ){ goeg = true; use23 = true;}
	else if (e2 && e3){ goee = true; use23 = true;}
	else if (f2 && f3 ){ goff = true;  use23 = true;}
      }
    }
    
    if(gogg&&goeg)cout<<"gg AND eg event!!!!! ------- PROBLEM!"<<endl;
    if(gogg&&goee)cout<<"gg AND ee event!!!!! ------- PROBLEM!"<<endl;
    if(gogg&&goff)cout<<"gg AND ff event!!!!! ------- PROBLEM!"<<endl;
    //    if(gogg&&gogammafake)cout<<"gg AND gammafake event!!!!! ------- PROBLEM!"<<endl;
    if(goeg&&goee)cout<<"eg AND ee event!!!!! ------- PROBLEM!"<<endl;
    if(goeg&&goff)cout<<"eg AND ff event!!!!! ------- PROBLEM!"<<endl;
    //    if(goeg&&gogammafake)cout<<"eg AND gammafake event!!!!! ------- PROBLEM!"<<endl;
    if(goee&&goff)cout<<"ee AND ff event!!!!! ------- PROBLEM!"<<endl;
    //    if(goee&&gogammafake)cout<<"ee AND gammafake event!!!!! ------- PROBLEM!"<<endl;
    //    if(goff&&gogammafake)cout<<"ff AND gammafake event!!!!! ------- PROBLEM!"<<endl;
    if(gogf&&gofg)cout<<"gf AND fg event!!!!! ------- PROBLEM!"<<endl;
    if((gogammafake && !gogf) && (gogammafake && !gofg))cout<<"gammafake BUT NOT gf OR fg event!!!!! ------- PROBLEM!"<<endl;
    return;
} //end of Categorize Event function



void SusyEventAnalyzer::MatchPhosToJets(TLorentzVector pOne, TLorentzVector pTwo, std::vector<TLorentzVector*> jetvecs, TLorentzVector &jet1, TLorentzVector &jet2, bool &hasdijetpt, float dR){
    hasdijetpt=false;
    for(std::vector<TLorentzVector*>::iterator jet_it1 = jetvecs.begin(); jet_it1 != jetvecs.end(); jet_it1++){
        if(isSameObject(pOne,*(*jet_it1),dR)){
            for(std::vector<TLorentzVector*>::iterator jet_it2 = jetvecs.begin(); jet_it2 != jetvecs.end(); jet_it2++){
                if( !isSameObject(*(*jet_it1), *(*jet_it2),0.1) ){
                    if(isSameObject(pTwo, *(*jet_it2), dR)){
                        jet1=**jet_it1;
                        jet2=**jet_it2;
                        hasdijetpt=true;
                        break;
                    }//end jet2 match
                }//end jet2!=jet1
            }//end jet2 iterator
            if(hasdijetpt){break;}
        }//end jet1 match
    }//end jet1 iterator
    return;
}//end of MatchPhosToJets function

//  LocalWords:  cout
