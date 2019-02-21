#define T5gg_analyzer_cxx
#include "T5gg_analyzer.h"
#include "lester_mt2_bisect.h"
#include <TH2.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TVector3.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include "TRandom3.h"
#include <TTree.h>
#include <algorithm>
#include <utility>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <tuple>
#include <fstream>
Long64_t weirdEvents[] = {54894719,16790527};
unsigned int numWeird = 2;


Int_t   NbinsMET      = 24;
Float_t XbinsMET[25]  = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,75,85,100,115, 130, 150,185,250,350};


typedef std::map<int,std::map<int,bool> > RunLumiFlagHolder;  //define map that holds json list
RunLumiFlagHolder goodrunlumilist;  // instantiate

//constants

float pi             = acos(-1.);
float Z0Mass_nominal = 91.186;

//Fiducial cuts

float nF1 = 1.4442;
float nF2 = 1.566;
float nF3 = 2.5;

//Counters
int Ientry_   =0; //total events
int IIentry_  =0; //events that pass the trigger
int IIIentry_ =0; //events that pass the METfilters
int IVentry_  =0; //events after 2 good object selection
int Ventry_   =0;//events after GG/EE/FF selection
int numGG = 0;
int numFF = 0;
int numEE = 0;
int numEG = 0;
int numEEGG = 0;
int recoEEGG =0;

//photon counters//

int Ngg_RemovedFromJets = 0;
int Ngg_NotOverlapping  = 0;


fstream myfile;

//Declare the output tree name

TTree *GGtree_;
TTree *EEtree_;
TTree *FFtree_;
TTree *EGtree_;


//histograms to be saved in the output tree
TH1F* h_events;
TH1F* h_NJet;
TH1F* ggDiEmPt;
TH1F* eeDiEmPt;
TH1F* ffDiEmPt;

TH1F* ggNvertex;
TH1F* eeNvertex;
TH1F* ffNvertex;


TH1F* h_ee_MET;
TH1F* h_ff_MET;
TH1F* h_gg_MET;

TH1F* ggNJets;
TH1F* eeNJets;
TH1F* ffNJets;

TH1F* h_E_MET;
TH1F* h_E_METpred;

TH1F* h_EG_MET;

TH1F* gg_events;
TH1F* ee_events;
TH1F* ff_events;
TH1F* eg_events;


TH1F* gg_events_nj;
TH1F* ee_events_nj;
TH1F* ff_events_nj;

TH1F* gg_events_st400;
TH1F* ee_events_st400;
TH1F* ff_events_st400;

TH1F* gg_events_STbin1;
TH1F* ee_events_STbin1;
TH1F* ff_events_STbin1;

TH1F* gg_events_STbin2;
TH1F* ee_events_STbin2;
TH1F* ff_events_STbin2;

TH1F* gg_events_STbin3;
TH1F* ee_events_STbin3;
TH1F* ff_events_STbin3;

TH1F* gg_events_STbin4;
TH1F* ee_events_STbin4;
TH1F* ff_events_STbin4;


TH2F* DiEmPTVSNJets_ee;
TH2F* DiEmPTVSNJets_ff;
TH2F* DiEmPTVSNJets_gg;


TH2F* DiEmPTVSNJetsVarbin_ee;
TH2F* DiEmPTVSNJetsVarbin_ff;
TH2F* DiEmPTVSNJetsVarbin_gg;

TH1F* egMEtMTCut;
TH1F* egMt;
TH1F* egMEtMtDeltaPhiCut;

TH1F* egInvMassMTCut;
TH1F* egInvMassnoMTCut;

TH1F* egInvMassCutDeltaPhiCut;

TH1F* eeggEvents;

//TH1F* egDeltaPhi;


Int_t Nbins = 24;
float Xbins[25] = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,70,85,100,115, 130, 150,175,200,300};


float BinsNJ[6] = {0,1,2,3,5,10};

//TH1F* h_ST;

//TH1F* h_MHT;

//TH1F* h_HT;

///////////////////////////////////////////////////

//branches to be saved in the output tree
//general info

Long64_t	Run_;
Long64_t        Event_;
Long64_t        isData_;
Long64_t        Lumis_;
float 		Rho_;
//PU info
float		PUtrue_;
float		PUweight_;
//met info
int		MET_Filters_;
float	      	MET_;
float	      	MET_Phi_;
float	      	MET_sumEt_;
float	      	MET_mEtSig_;
float	      	MET_Sig_;
float           MT2_;
float           MT_;
//photon info
ULong64_t	NPho_;
ULong64_t	NFake_;
float		Pho_E_[50];
float		Pho_Et_[50];
float		Pho_Eta_[50];
float		Pho_Phi_[50];
int     	Pho_hasPixelSeed_[50];
int     	Pho_EleVeto_[50];
float		Pho_R9_[50];
float		Pho_HoverE_[50];
float		Pho_SigmaIEtaIEta_[50];
float		Pho_SigmaIEtaIPhi_[50];
float		Pho_SigmaIPhiIPhi_[50];
float		Pho_IDMVA_[50];
int	 	Pho_IDbit_[50];
float 		Dipho_Pt_;
float 		Dipho_Mass_;
float 		value[50];

//electron info
ULong64_t	NEle_;
float		Ele_Pt_[50];
float		Ele_Eta_[50];
float		Ele_Phi_[50];
float		Ele_SigmaIEtaIEta_[50];
float		Ele_SigmaIEtaIPhi_[50];
float		Ele_SigmaIPhiIPhi_[50];
float		Ele_SigmaIEtaIEtaFull5x5_[50];
float		Ele_SigmaIPhiIPhiFull5x5_[50];
ULong64_t	NMu_;
ULong64_t       NJet_;
ULong64_t       NJet_n5_;


float gg_sample =0;
float ee_sample =0;
float ff_sample =0;

float gg_sample_nj =0;
float ee_sample_nj =0;
float ff_sample_nj =0;

float gg_sample_st400 =0;
float ee_sample_st400 =0;
float ff_sample_st400 =0;

float gg_sample_bin1=0;
float ee_sample_bin1=0;
float ff_sample_bin1=0;

float gg_sample_bin2=0;
float ee_sample_bin2=0;
float ff_sample_bin2=0;

float gg_sample_bin3=0;
float ee_sample_bin3=0;
float ff_sample_bin3=0;

float gg_sample_bin4=0;
float ee_sample_bin4=0;
float ff_sample_bin4=0;



float HT_;
float ST_;
float MHT_;
float MHTPhi_;



float EffAreas[7][3] = {
  {0.0360, 0.0597, 0.1210},
  {0.0377, 0.0807, 0.1107},
  {0.0306, 0.0629, 0.0699},
  {0.0283, 0.0197, 0.1056},
  {0.0254, 0.0184, 0.1457},
  {0.0217, 0.0284, 0.1719},
  {0.0167, 0.0591, 0.1998}
};

//output root file

TFile fout("./analysis_NewMuonRecoVetoRemoveFinal.root","RECREATE");



void T5gg_analyzer::Begin(TTree * /*tree*/)
{



GGtree_ = new TTree("ggtree","Event data");
EEtree_ = new TTree("eetree","Event data");
FFtree_ = new TTree("fftree","Event data");
EGtree_ = new TTree("egtree","Event data");


////////////////////////////////////////////////

h_events = new TH1F("h_events", "events breakdown",5,0,5);
h_NJet   = new TH1F("h_NJet","h_NJet",5,0,5);
ggDiEmPt = new TH1F("ggDiEmPt","gg",Nbins,Xbins);
eeDiEmPt = new TH1F("eeDiEmPt","ee",Nbins,Xbins);
ffDiEmPt = new TH1F("ffDiEmPt","ff",Nbins,Xbins);


gg_events = new TH1F("gg_events","gg_events",1,0,1);
ee_events = new TH1F("ee_events","ee_events",1,0,1);
ff_events = new TH1F("ff_events","ff_events",1,0,1);

gg_events_st400 = new TH1F("gg_events_st400","gg_events_st400",1,0,1);
ee_events_st400 = new TH1F("ee_events_st400","ee_events_st400",1,0,1);
ff_events_st400 = new TH1F("ff_events_st400","ff_events_st400",1,0,1);


gg_events_nj = new TH1F("gg_events_nj","gg_events_nj",1,0,1);
ee_events_nj = new TH1F("ee_events_nj","ee_events_nj",1,0,1);
ff_events_nj = new TH1F("ff_events_nj","ff_events_nj",1,0,1);

gg_events_STbin1 = new TH1F("gg_events_STbin1","gg_events_STbin1",1,0,1);
ee_events_STbin1 = new TH1F("ee_events_STbin1","ee_events_STbin1",1,0,1);
ff_events_STbin1 = new TH1F("ff_events_STbin1","ff_events_STbin1",1,0,1);

gg_events_STbin2 = new TH1F("gg_events_STbin2","gg_events_STbin2",1,0,1);
ee_events_STbin2 = new TH1F("ee_events_STbin2","ee_events_STbin2",1,0,1);
ff_events_STbin2 = new TH1F("ff_events_STbin2","ff_events_STbin2",1,0,1);

gg_events_STbin3 = new TH1F("gg_events_STbin3","gg_events_STbin3",1,0,1);
ee_events_STbin3 = new TH1F("ee_events_STbin3","ee_events_STbin3",1,0,1);
ff_events_STbin3 = new TH1F("ff_events_STbin3","ff_events_STbin3",1,0,1);

gg_events_STbin4 = new TH1F("gg_events_STbin4","gg_events_STbin4",1,0,1);
ee_events_STbin4 = new TH1F("ee_events_STbin4","ee_events_STbin4",1,0,1);
ff_events_STbin4 = new TH1F("ff_events_STbin4","ff_events_STbin4",1,0,1);


DiEmPTVSNJets_ee =new TH2F("DiEmPTVSNJets_ee","DiEmPTVSNJets_ee",50,0,600,10,0,10);
DiEmPTVSNJets_ff =new TH2F("DiEmPTVSNJets_ff","DiEmPTVSNJets_ff",50,0,600,10,0,10);
DiEmPTVSNJets_gg =new TH2F("DiEmPTVSNJets_gg","DiEmPTVSNJets_gg",50,0,600,10,0,10);

ggNvertex = new TH1F("ggNvertex","gg PU distribution",50,0,100);
ggNJets   = new TH1F("ggNJets","gg Jet Distribution" ,10,0,10);

eeNvertex = new TH1F("eeNvertex","ee PU distribution",50,0,100);
eeNJets   = new TH1F("eeNJets","ee Jet Distribution" ,10,0,10);

ffNvertex = new TH1F("ffNvertex","ff PU distribution",50,0,100);
ffNJets   = new TH1F("ffNJets","ff Jet Distribution" ,10,0,10);

DiEmPTVSNJetsVarbin_ee = new TH2F("DiEmPTVSNJetsVarbin_ee","DiEmPTVSNJetsVarbin_ee", 26, Xbins, 5,BinsNJ);
DiEmPTVSNJetsVarbin_ff = new TH2F("DiEmPTVSNJetsVarbin_ff","DiEmPTVSNJetsVarbin_ff", 26, Xbins, 5,BinsNJ);
DiEmPTVSNJetsVarbin_gg = new TH2F("DiEmPTVSNJetsVarbin_gg","DiEmPTVSNJetsVarbin_gg", 26, Xbins, 5,BinsNJ);


h_E_MET    = new TH1F ("h_E_MET","h_E_MET",30,0,150);

h_E_METpred= new TH1F ("h_E_METpred","h_E_METpred",NbinsMET,XbinsMET);


h_EG_MET = new TH1F ("h_EG_MET","h_EG_MET",NbinsMET,XbinsMET);


h_ee_MET= new TH1F("h_ee_MET","h_ee_MET",NbinsMET,XbinsMET);
h_ff_MET= new TH1F("h_ff_MET","h_ff_MET",NbinsMET,XbinsMET);
h_gg_MET= new TH1F("h_gg_MET","h_gg_MET",NbinsMET,XbinsMET);

eeggEvents = new TH1F("eeggEvents","eeggEvents",1,0,1);


egMEtMTCut         = new TH1F("egMEtMTCut","egMEtMTCut",NbinsMET,XbinsMET);
egMt               = new TH1F("egMt","egMt",NbinsMET,XbinsMET);
egMEtMtDeltaPhiCut = new TH1F("egMEtMtDeltaPhiCut","egMEtMtDeltaPhiCut",NbinsMET,XbinsMET);

egInvMassMTCut            = new TH1F("egInvMassMTCut","egInvMassMTCut",50,75,130);
egInvMassnoMTCut          = new TH1F("egInvMassnoMTCut","egInvMassnoMTCut",50,75,130);

egInvMassCutDeltaPhiCut   = new TH1F("egInvMassCutDeltaPhiCut","eg with MT<100 GeV and",50,75,130);
//egDeltaPhi                = new TH1F();


//Branches



//general event info
 GGtree_->Branch("Run",&Run_,"Run/I");
 GGtree_->Branch("Event",&Event_,"Event/I");
 GGtree_->Branch("Lumis",&Lumis_,"Lumis/I");
 GGtree_->Branch("Rho",&Rho_,"Rho/F");
 GGtree_->Branch("PUtrue",&PUtrue_,"PUtrue/F");
 GGtree_->Branch("PUweight",&PUweight_,"PUweight/F");
 GGtree_->Branch("NJet",&NJet_,"NJet/I");

 //MET - genMET info
 GGtree_->Branch("MET_Filters",&MET_Filters_,"MET_Filters/I");
 GGtree_->Branch("MET",&MET_,"MET/F");
 GGtree_->Branch("MET_Phi",&MET_Phi_,"MET_Phi/F");
 GGtree_->Branch("MET_sumEt",&MET_sumEt_,"MET_sumEt/F");
 GGtree_->Branch("MET_mEtSig",&MET_mEtSig_,"MET_mEtSig/F");
 GGtree_->Branch("MET_Sig",&MET_Sig_,"MET_Sig/F");
 GGtree_->Branch("MT2",&MT2_,"MT2/F");




 //photon info
 GGtree_->Branch("NPho",&NPho_,"NPho/I");
 GGtree_->Branch("Pho_Et",&Pho_Et_,"Pho_Et[NPho]/F");
 GGtree_->Branch("Pho_Eta",&Pho_Eta_,"Pho_Eta[NPho]/F");
 GGtree_->Branch("Pho_Phi",&Pho_Phi_,"Pho_Phi[NPho]/F");
 GGtree_->Branch("Pho_hasPixelSeed",&Pho_hasPixelSeed_,"Pho_hasPixelSeed[NPho]/I");
 GGtree_->Branch("Pho_EleVeto",&Pho_EleVeto_,"Pho_EleVeto[NPho]/I");
 GGtree_->Branch("Pho_R9",&Pho_R9_,"Pho_R9[NPho]/F");
 GGtree_->Branch("Pho_HoverE",&Pho_HoverE_,"Pho_HoverE[NPho]/F");
 GGtree_->Branch("Pho_SigmaIEtaIEta",&Pho_SigmaIEtaIEta_,"Pho_SigmaIEtaIEta[NPho]/F");
 GGtree_->Branch("Pho_SigmaIEtaIPhi",&Pho_SigmaIEtaIPhi_,"Pho_SigmaIEtaIPhi[NPho]/F");
 GGtree_->Branch("Pho_SigmaIPhiIPhi",&Pho_SigmaIPhiIPhi_,"Pho_SigmaIPhiIPhi[NPho]/F");
 GGtree_->Branch("Pho_IDMVA",&Pho_IDMVA_,"Pho_IDMVA[NPho]/F");
 GGtree_->Branch("Pho_IDbit",&Pho_IDbit_,"Pho_IDbit[NPho]/I");

 //ST, MHT, HT info

 GGtree_->Branch("HT",&HT_,"HT/F");
 GGtree_->Branch("ST",&ST_,"ST/F");
 GGtree_->Branch("MHT",&MHT_,"MHT/F");
 GGtree_->Branch("MHTPhi",&MHTPhi_,"MHTPhi/F");


 //other new stuff
 GGtree_->Branch("Dipho_Pt",&Dipho_Pt_,"Dipho_Pt/F");
 GGtree_->Branch("Dipho_Mass",&Dipho_Mass_,"Dipho_Mass/F");


 //general event info
 EEtree_->Branch("Run",&Run_,"Run/I");
 EEtree_->Branch("Event",&Event_,"Event/I");
 EEtree_->Branch("Lumis",&Lumis_,"Lumis/I");
 EEtree_->Branch("Rho",&Rho_,"Rho/F");
 EEtree_->Branch("PUtrue",&PUtrue_,"PUtrue/F");
 EEtree_->Branch("PUweight",&PUweight_,"PUweight/F");
 EEtree_->Branch("NJet",&NJet_,"NJet/I");



 //MET - genMET info

 EEtree_->Branch("MET_Filters",&MET_Filters_,"MET_Filters/I");
 EEtree_->Branch("MET",&MET_,"MET/F");
 EEtree_->Branch("MET_Phi",&MET_Phi_,"MET_Phi/F");
 EEtree_->Branch("MET_sumEt",&MET_sumEt_,"MET_sumEt/F");
 EEtree_->Branch("MET_mEtSig",&MET_mEtSig_,"MET_mEtSig/F");
 EEtree_->Branch("MET_Sig",&MET_Sig_,"MET_Sig/F");
 EEtree_->Branch("MT2",&MT2_,"MT2/F");



 //electron info
 EEtree_->Branch("NEle",&NEle_,"NEle/I");
 EEtree_->Branch("Dipho_Pt",&Dipho_Pt_,"Dipho_Pt/F");
 EEtree_->Branch("Dipho_Mass",&Dipho_Mass_,"Dipho_Mass/F");
 EEtree_->Branch("Ele_Pt",&Ele_Pt_,"Ele_Pt[NEle]/F");
 EEtree_->Branch("Ele_Eta",&Ele_Eta_,"Ele_Eta[NEle]/F");
 EEtree_->Branch("Ele_Phi",&Ele_Phi_,"Ele_Phi[NEle]/F");
 EEtree_->Branch("Ele_SigmaIEtaIEta",&Ele_SigmaIEtaIEta_,"Ele_SigmaIEtaIEta[NEle]/F");
 EEtree_->Branch("Ele_SigmaIEtaIPhi",&Ele_SigmaIEtaIPhi_,"Ele_SigmaIEtaIPhi[NEle]/F");
 EEtree_->Branch("Ele_SigmaIPhiIPhi",&Ele_SigmaIPhiIPhi_,"Ele_SigmaIPhiIPhi[NEle]/F");


 //ST, MHT, HT info

 EEtree_->Branch("HT",&HT_,"HT/F");
 EEtree_->Branch("ST",&ST_,"ST/F");
 EEtree_->Branch("MHT",&MHT_,"MHT/F");
 EEtree_->Branch("MHTPhi",&MHTPhi_,"MHTPhi/F");



 //general event info
 FFtree_->Branch("Run",&Run_,"Run/I");
 FFtree_->Branch("Event",&Event_,"Event/I");
 FFtree_->Branch("Lumis",&Lumis_,"Lumis/I");
 FFtree_->Branch("Rho",&Rho_,"Rho/F");
 FFtree_->Branch("PUtrue",&PUtrue_,"PUtrue/F");
 FFtree_->Branch("PUweight",&PUweight_,"PUweight/F");
 FFtree_->Branch("NJet",&NJet_,"NJet/I");



 //fake info MET
 FFtree_->Branch("MET_Filters",&MET_Filters_,"MET_Filters/I");
 FFtree_->Branch("MET",&MET_,"MET/F");
 FFtree_->Branch("MET_Phi",&MET_Phi_,"MET_Phi/F");
 FFtree_->Branch("MET_sumEt",&MET_sumEt_,"MET_sumEt/F");
 FFtree_->Branch("MET_mEtSig",&MET_mEtSig_,"MET_mEtSig/F");
 FFtree_->Branch("MET_Sig",&MET_Sig_,"MET_Sig/F");
 FFtree_->Branch("MT2",&MT2_,"MT2/F");


 //ST, MHT, HT info

 FFtree_->Branch("HT",&HT_,"HT/F");
 FFtree_->Branch("ST",&ST_,"ST/F");
 FFtree_->Branch("MHT",&MHT_,"MHT/F");
 FFtree_->Branch("MHTPhi",&MHTPhi_,"MHTPhi/F");


 //fake info
 FFtree_->Branch("NFake",&NFake_,"NFake/I");
 FFtree_->Branch("Dipho_Pt",&Dipho_Pt_,"Dipho_Pt/F");
 FFtree_->Branch("Dipho_Mass",&Dipho_Mass_,"Dipho_Mass/F");
 FFtree_->Branch("Pho_E",&Pho_E_,"Pho_E[NFake]/F");
 FFtree_->Branch("Pho_Et",&Pho_Et_,"Pho_Et[NFake]/F");
 FFtree_->Branch("Pho_Eta",&Pho_Eta_,"Pho_Eta[NFake]/F");
 FFtree_->Branch("Pho_Phi",&Pho_Phi_,"Pho_Phi[NFake]/F");
 FFtree_->Branch("Pho_hasPixelSeed",&Pho_hasPixelSeed_,"Pho_hasPixelSeed[NFake]/I");
 FFtree_->Branch("Pho_EleVeto",&Pho_EleVeto_,"Pho_EleVeto[NFake]/I");
 FFtree_->Branch("Pho_R9",&Pho_R9_,"Pho_R9[NFake]/F");
 FFtree_->Branch("Pho_HoverE",&Pho_HoverE_,"Pho_HoverE[NFake]/F");
 FFtree_->Branch("Pho_SigmaIEtaIEta",&Pho_SigmaIEtaIEta_,"Pho_SigmaIEtaIEta[NFake]/F");
 FFtree_->Branch("Pho_SigmaIEtaIPhi",&Pho_SigmaIEtaIPhi_,"Pho_SigmaIEtaIPhi[NFake]/F");
 FFtree_->Branch("Pho_SigmaIPhiIPhi",&Pho_SigmaIPhiIPhi_,"Pho_SigmaIPhiIPhi[NFake]/F");
 FFtree_->Branch("Pho_IDMVA",&Pho_IDMVA_,"Pho_IDMVA[NFake]/F");
 FFtree_->Branch("Pho_IDbit",&Pho_IDbit_,"Pho_IDbit[NFake]/I");




 //general event info
EGtree_->Branch("Run",&Run_,"Run/I");
EGtree_->Branch("Event",&Event_,"Event/I");
EGtree_->Branch("Lumis",&Lumis_,"Lumis/I");
EGtree_->Branch("Rho",&Rho_,"Rho/F");
EGtree_->Branch("PUtrue",&PUtrue_,"PUtrue/F");
EGtree_->Branch("PUweight",&PUweight_,"PUweight/F");
EGtree_->Branch("NJet",&NJet_,"NJet/I");



//MET - genMET info

EGtree_->Branch("MET_Filters",&MET_Filters_,"MET_Filters/I");
EGtree_->Branch("MET",&MET_,"MET/F");
EGtree_->Branch("MET_Phi",&MET_Phi_,"MET_Phi/F");
EGtree_->Branch("MET_sumEt",&MET_sumEt_,"MET_sumEt/F");
EGtree_->Branch("MET_mEtSig",&MET_mEtSig_,"MET_mEtSig/F");
EGtree_->Branch("MET_Sig",&MET_Sig_,"MET_Sig/F");
EGtree_->Branch("MT2",&MT2_,"MT2/F");
EGtree_->Branch("MT2",&MT2_,"MT2/F");
EGtree_->Branch("MT",&MT_,"MT/F");


//electron info
EGtree_->Branch("NEle",&NEle_,"NEle/I");
EGtree_->Branch("Dipho_Pt",&Dipho_Pt_,"Dipho_Pt/F");
EGtree_->Branch("Dipho_Mass",&Dipho_Mass_,"Dipho_Mass/F");
EGtree_->Branch("Ele_Pt",&Ele_Pt_,"Ele_Pt[NEle]/F");
EGtree_->Branch("Ele_Eta",&Ele_Eta_,"Ele_Eta[NEle]/F");
EGtree_->Branch("Ele_Phi",&Ele_Phi_,"Ele_Phi[NEle]/F");
EGtree_->Branch("Ele_SigmaIEtaIEta",&Ele_SigmaIEtaIEta_,"Ele_SigmaIEtaIEta[NEle]/F");
EGtree_->Branch("Ele_SigmaIEtaIPhi",&Ele_SigmaIEtaIPhi_,"Ele_SigmaIEtaIPhi[NEle]/F");
EGtree_->Branch("Ele_SigmaIPhiIPhi",&Ele_SigmaIPhiIPhi_,"Ele_SigmaIPhiIPhi[NEle]/F");


//ST, MHT, HT info

EGtree_->Branch("HT",&HT_,"HT/F");
EGtree_->Branch("ST",&ST_,"ST/F");
EGtree_->Branch("MHT",&MHT_,"MHT/F");
EGtree_->Branch("MHTPhi",&MHTPhi_,"MHTPhi/F");


 IncludeAJson("Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.json");

}




////////////////////////////////////////////////////////
void T5gg_analyzer::SlaveBegin(TTree * /*tree*/)
{}
///////////////////////////////////////////////////////


Bool_t T5gg_analyzer::Process(Long64_t entry)
{




 GetAllEntries(entry);
 if(Ientry_%10000==0) cout<<"Processing #entry: "<<Ientry_<<endl;
 Ientry_++;





  if (!isInJson(run, lumis)) return kFALSE;

  std::map<int, std::set<int> > allEvents;
   if(isData){
     bool duplicateEvent = ! (allEvents[run].insert(event) ).second;
     if(duplicateEvent){
       cout<<"Duplicate Event! Run "<< run <<" Event "<< event<<endl;
       return kFALSE;
     }
   }



  //////////////////////////////////////////////////////
 //////////////////////
 //  information //
 //////////////////////
 if(!(HLTPho>>14&1) && !(HLTPho>>15&1))  return kFALSE;
 IIentry_++;

 //cout << "value of MET filter " << metFilters>>10&0 << endl;
 if(!(metFilters >>9&1) && !(metFilters>>10&1)) return kFALSE;


 for (int i=1; i< 9; i++){

 if (metFilters>>i&1) return kFALSE;


 }


 if ( metFilters >>11&1 ) return kFALSE;

 //IIIentry_++;




   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ////////////////////
 // MC information //
 ////////////////////


    ////////////////////////////
   // Good Vertex Selection  //
  ////////////////////////////

 if(nVtx<1 || !isPVGood==true ) return kFALSE;


  //////////////////
 // Reco Photons //
 /////////////////

 int NDirtyPho = 0;
 std::map<float,int>    phoMap;
 vector<TLorentzVector> dirtyPhotons;
 for(int ipho=0; ipho<nPho; ipho++){

  TLorentzVector thisPhoton;
  thisPhoton = p4vec((*phoCalibEt)[ipho], (*phoPhi)[ipho], (*phoCalibEt)[ipho]*sinh((*phoEta)[ipho]));
  phoMap[thisPhoton.Pt()] = ipho;

  if ((*phoCalibEt)[ipho]<40.)              continue; //Et>40
  if (abs((*phoSCEta)[ipho]) > 1.4442)  continue; //EB only
  if ((*phohasPixelSeed)[ipho] !=0)    continue; //PSV
  if((*phoR9)[ipho] >= 1.0)              continue; //R9<1
  if((*phoR9)[ipho] <= 0.5)             continue; //R9>0.5
  if((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.005) continue; //sinin>0.005
  if(IsGoodPhoton(ipho,1,0) != 1)      continue; //mediumID


  dirtyPhotons.push_back(thisPhoton);
  NDirtyPho++;

  }



  ///////////
 // reco  fakes //
 ///////////
 int NDirtyFake = 0;
 std::map<float,int>    fakeMap;
 vector<TLorentzVector> dirtyFakes;
 for(int ipho=0; ipho<nPho; ipho++){

   TLorentzVector thisFake;
   thisFake = p4vec((*phoCalibEt)[ipho], (*phoPhi)[ipho], (*phoCalibEt)[ipho]*sinh((*phoEta)[ipho]));
   fakeMap[thisFake.Pt()] = ipho;

   if((*phoCalibEt)[ipho] < 40.)              continue; //Et>40
   if(abs((*phoSCEta)[ipho]) > 1.4442)   continue; //EB only
   if((*phohasPixelSeed)[ipho] != 0)     continue; //PSV
   if((*phoR9)[ipho] >= 0.9)              continue; //R9<1
   if((*phoR9)[ipho] <= 0.5)             continue; //R9>0.5
   if((*phoSigmaIEtaIEtaFull5x5)[ipho] < 0.005) continue; //sinin>0.005
   if((*phoSigmaIEtaIEtaFull5x5)[ipho] > 0.015) continue;
   if(IsGoodPhoton(ipho,1,0) != -1)      continue; //fake mediumID


   dirtyFakes.push_back(thisFake);
   NDirtyFake++;
 }



      ////////////////////
     // Reco Electrons //
    ////////////////////


  int NDirtyEle = 0;
  int NEEele =0;

  std::map<float,int> eleMap;
  std::map<float,int> EEeleMap;
  std::map<float,int> EleVetoMap;
  vector<TLorentzVector> dirtyElectrons;
  vector<TLorentzVector> EndCupElectrons;
  vector<TLorentzVector> VetoElectrons;

  for(int iele=0; iele<nPho; iele++){

   TLorentzVector thisElectron;
   thisElectron = p4vec((*phoCalibEt)[iele], (*phoPhi)[iele], (*phoCalibEt)[iele]*sinh((*phoEta)[iele]));
   eleMap[thisElectron.Pt()] = iele;


   //if(abs((*phoSCEta)[iele]) < 2.5 && (abs((*phoSCEta)[iele]) > 1.566)) return kFALSE; //EndCupElectrons.push_back(thisElectron); //continue; //electrons in the barrel


   if((*phoCalibEt)[iele] < 40.)             continue; //Pt>40
   if(abs((*phoSCEta)[iele]) > 1.4442)  continue; //EB only
   //if(abs((*phoSCEta)[iele]) > 2.4)  continue;


   if((*phohasPixelSeed)[iele] == 0)    continue; //PSV reversed
   if((*phoR9)[iele] >= 1.0)             continue; //R9<1
   if((*phoR9)[iele] <= 0.5)             continue; //R9>0.5
   if((*phoSigmaIEtaIEtaFull5x5)[iele] < 0.005) continue; //sinin>0.005
   if(IsGoodPhoton(iele,1,0) != 1)        continue; //mediumID





   dirtyElectrons.push_back(thisElectron);
   NDirtyEle++;



 }

for(int iele=0; iele<nEle; iele++){

 if ((*eleIDbit)[iele]>>2&1) {  //medium ID

     //if (abs((*eleSCEta)[iele]-(*eleEta)[iele])>0.2) return kFALSE;
   if (abs((*eleSCEta)[iele]-(*eleEta)[iele])>0.2) return kFALSE;


   }

}

 for(int iele=0; iele<nPho; iele++){


   TLorentzVector thisEEElectron;
   thisEEElectron = p4vec((*phoCalibEt)[iele], (*phoPhi)[iele], (*phoCalibEt)[iele]*sinh((*phoEta)[iele]));
   EEeleMap[thisEEElectron.Pt()] = iele;

   if((*phoCalibEt)[iele] < 40.)             continue; //Pt>40
   if(abs((*phoSCEta)[iele]) > 2.5) continue;

   if (abs((*phoSCEta)[iele]) < 1.566)  continue;


   if((*phohasPixelSeed)[iele] == 0)    continue; //PSV reversed
   if((*phoR9)[iele] >= 1.0)             continue; //R9<1
   if((*phoR9)[iele] <= 0.5)             continue; //R9>0.5
   if((*phoSigmaIEtaIEtaFull5x5)[iele] < 0.005) continue; //sinin>0.005
   if(IsGoodPhoton(iele,1,0) != 2)        continue; //mediumID


   EndCupElectrons.push_back(thisEEElectron);  //electrons in the barrel

   NEEele++;

 }


   /////////////////////////////////////////////
  ////RECO photon+lepton's analysis electrons ///
 ///////////////////////////////////////////////

  int nVetoEleEB=0;
  int nVetoEleEE=0;
  int More3VetoEle;
  int TotalVetoEle=0;


  for(int iEle=0; iEle<nEle; iEle++){


   TLorentzVector thisVetoElectron;
   thisVetoElectron = p4vec((*eleCalibPt)[iEle], (*elePhi)[iEle], (*eleCalibPt)[iEle]*sinh((*eleEta)[iEle]));
   EleVetoMap[thisVetoElectron.Pt()] = iEle;


   if((*eleCalibPt)[iEle]<25) continue;

   if (fabs((*eleSCEta)[iEle])< 1.4442){
       if ((*eleR9)[iEle]<0.5)                    continue;
       if ((*eleSigmaIEtaIEtaFull5x5)[iEle]>0.00998)      continue;
       if (fabs((*eledEtaAtVtx)[iEle])>0.00311)   continue;
       if (fabs((*eledPhiAtVtx)[iEle])>0.103)     continue;
       if ((*eleHoverE)[iEle]>0.253)         continue;
       if (fabs((*eleEoverPInv)[iEle])>0.134)   continue;
       if ((*eleMissHits)[iEle]>1)            continue;
       if ((*eleConvVeto)[iEle]== 0)           continue;
       if ((*elePFMiniIso)[iEle]>0.1)         continue;

  nVetoEleEB++;
  VetoElectrons.push_back(thisVetoElectron);

  }

  if (fabs((*eleSCEta)[iEle])>1.56 && fabs((*eleSCEta)[iEle])<2.5){

       if ((*eleR9)[iEle]<0.8)                    continue;
       if ((*eleSigmaIEtaIEtaFull5x5)[iEle]>0.0298)      continue;
       if (fabs((*eledEtaAtVtx)[iEle])>0.00609)   continue;
       if (fabs((*eledPhiAtVtx)[iEle])>0.045)     continue;
       if ((*eleHoverE)[iEle]>0.0878)         continue;
       if (fabs((*eleEoverPInv)[iEle])>0.13)   continue;
       if ((*eleMissHits)[iEle]>1)            continue;
       if ((*eleConvVeto)[iEle]== 0)          continue;
       if ((*elePFMiniIso)[iEle]>0.1)         continue;

  VetoElectrons.push_back(thisVetoElectron);
  nVetoEleEE++;

 }

  TotalVetoEle = (nVetoEleEB +nVetoEleEE);

  }

 sort(VetoElectrons.begin(), VetoElectrons.end(), P4SortCondition);

 if (TotalVetoEle >=3) {
   More3VetoEle++;
  // cout <<"3 or more veto electrons " << More3VetoEle << " for event " << event << endl;

 }



  ////////////////
 // RECO muons //
 ////////////////
 NMu_ = 0;
 float MuNumber;
 std::map<float,int>    muMap;
 vector<TLorentzVector> muons;
 for(int imu=0; imu<nMu; imu++){

   TLorentzVector thisMuon;
   thisMuon = p4vec((*muPt)[imu], (*muPhi)[imu], (*muPt)[imu]*sinh((*muEta)[imu]));
   muMap[thisMuon.Pt()] = imu;

   float muIso = ( (*muPFChIso)[imu] + std::max(0., (*muPFNeuIso)[imu]+(*muPFPhoIso)[imu]-0.5*((*muPFPUIso)[imu])) );

  /* old muon reconstruction
    if((*muPt)[imu] < 30.)          continue; //Pt>30
   //if(abs((*muEta)[imu]) > 1.4442) continue; //EB only
   if(abs((*muEta)[imu]) > 2.4) continue;
   //if((*muIsMediumID)[imu] == 0 )  continue; //mediumID
   if (!( ((*muIDbit)[imu])>>0&1)) continue; //muon loose ID
   if(muIso/((*muPt)[imu]) > 0.25 ) continue; //PFbased isolation
   //if(((*muIsoTrk)[imu])/((*muPt)[imu]) > 0.10) continue; //TKbased isolation

   */

   if((*muPt)[imu] < 25.)          continue; //Pt>25
   if(abs((*muEta)[imu]) > 2.40)   continue;
   if (!( ((*muIDbit)[imu])>>1&1)) continue; //muon medium ID
   if (fabs((*muD0)[imu])>0.05 )    continue; //muon D0<005
   if (fabs((*muDz)[imu])>0.1)      continue;
   if ((*muPFMiniIso)[imu]>0.2)    continue;

   muons.push_back(thisMuon);
   NMu_++;
 }
 sort(muons.begin(), muons.end(), P4SortCondition);

 MuNumber = NMu_ ;
 //cout<<"number of muons  " << NMu_  <<"for event "<< event <<endl;


 ///////////////
 // RECO Jets //
 ///////////////
 int NDirtyJet = 0;
 std::map<float,int>    jetMap;
 vector<TLorentzVector> dirtyJets_n5;
 vector<TLorentzVector> dirtyJets;


 bool weird = false;


 for(unsigned int i = 0; i < numWeird; i++){
   if( event == weirdEvents[i]){
     weird = true;
     ///cout << "============================================" << endl;
    /// cout << "Processing event " << event << endl;
    // cout << "Number of raw jets is " << nJet << endl;
   }
 }



 for(int ijet=0; ijet<nJet; ijet++){

   TLorentzVector thisJet;
   thisJet = p4vec((*jetPt)[ijet], (*jetPhi)[ijet], (*jetPt)[ijet]*sinh((*jetEta)[ijet]));
   jetMap[thisJet.Pt()] = ijet;



   if(TMath::Abs((*jetEta)[ijet]) > 5.0 )   continue; //n<5
   if((*jetPt)[ijet] < 30. )                continue; //Pt>30
   if((*jetPFLooseId)[ijet] == 0 )          continue; //looseID

   dirtyJets_n5.push_back(thisJet);

   if(TMath::Abs((*jetEta)[ijet]) > 2.4 )   continue; //n<2.4

   dirtyJets.push_back(thisJet);
   NDirtyJet++;
 }

 if(weird) cout << "Number of jets passing eta, loose ID, and pt " << NDirtyJet << endl;


 ////////////////////////////
 // Object Cleaning    /////
 //////////////////////////

 NPho_   = 0;
 NEle_   = 0;
 NJet_   = 0;
 NJet_n5_ =0;
 NFake_   =0;
 NMu_     = 0;

 float NumberEl;
 float NumElMore3;

 NumElMore3=0;
 vector<TLorentzVector> photons;
 vector<TLorentzVector> electrons;
 vector<TLorentzVector> jets_n5;
 vector<TLorentzVector> jets;
 vector<TLorentzVector> fakes;

 bool isElectron;

 //e cleaning; remove e that overlap with mu with dR<0.3
 for(unsigned int Iele=0; Iele<dirtyElectrons.size(); Iele++){

   TLorentzVector thisElectronToClean = dirtyElectrons[Iele];
   std::pair<float,int> pairEleMu  = minDR_ValueAndIndex(thisElectronToClean, muons);

   if(pairEleMu.first < 0.3) continue;

   electrons.push_back(thisElectronToClean);

   NEle_++;

 }


NumberEl=NEle_;

if ( NumberEl>=3 ) NumElMore3++;
 sort(electrons.begin(), electrons.end(), P4SortCondition);

 float More3_Ele;

 More3_Ele = NumElMore3;

 //gamma cleaning; remove gamma that overlap with e, mu with dR<0.3, reject Menglei's electrons that do NOT overlap within dR<0.3 with photons
 for(unsigned int Ipho=0; Ipho<dirtyPhotons.size(); Ipho++){

   TLorentzVector thisPhotonToClean = dirtyPhotons[Ipho];
   std::pair<float,int> pairPhoEle = minDR_ValueAndIndex(thisPhotonToClean, electrons);
   std::pair<float,int> pairPhoMu     = minDR_ValueAndIndex(thisPhotonToClean, muons);
   std::pair<float,int> pairPhoVetoEle= minDR_ValueAndIndex(thisPhotonToClean, VetoElectrons);

   if(pairPhoEle.first       < 0.3) continue;
   if(pairPhoMu.first        < 0.3) continue;
  // if(pairPhoVetoEle.first   > 0.3) continue;

   photons.push_back(thisPhotonToClean);
   NPho_++;
 }


 sort(photons.begin(), photons.end(), P4SortCondition);

 //fake cleaning; remove jets that overlap with gamma,e,mu with dR<0.4
 for(unsigned int Ifake=0; Ifake<dirtyFakes.size(); Ifake++){

   TLorentzVector thisFakeToClean = dirtyFakes[Ifake];
   std::pair<float,int> pairFakePho = minDR_ValueAndIndex(thisFakeToClean, photons);
   std::pair<float,int> pairFakeEle = minDR_ValueAndIndex(thisFakeToClean, electrons);
   std::pair<float,int> pairFakeMu  = minDR_ValueAndIndex(thisFakeToClean, muons);
   if(pairFakePho.first < 0.4) continue;
   if(pairFakeEle.first < 0.4) continue;
   if(pairFakeMu.first  < 0.4) continue;

   fakes.push_back(thisFakeToClean);
   NFake_++;



 }
 sort(fakes.begin(), fakes.end(), P4SortCondition);


 //jet cleaning; remove jets that overlap with gamma,e,mu with dR<0.4
 for(unsigned int Ijet=0; Ijet<dirtyJets.size(); Ijet++){

   bool overlapping = false;

   TLorentzVector thisJetToClean = dirtyJets[Ijet];
   std::pair<float,int> pairJetPho  = minDR_ValueAndIndex(thisJetToClean, photons);
   std::pair<float,int> pairJetEle  = minDR_ValueAndIndex(thisJetToClean, electrons);
   std::pair<float,int> pairJetMu   = minDR_ValueAndIndex(thisJetToClean, muons);
   std::pair<float,int> pairJetFake = minDR_ValueAndIndex(thisJetToClean, fakes);


   float percent_pho[100];
   float percent_ele[100];
   float percent_mu[100];
   float percent_fake[100];


   for (unsigned int Iele=0; Iele<electrons.size(); Iele++){

   TLorentzVector this_electron = electrons[Iele];

   percent_ele[Iele] = abs((thisJetToClean.Pt()- this_electron.Pt())/(this_electron.Pt()));

     if (thisJetToClean.DeltaR(electrons[Iele])<0.4){

       if(percent_ele[Iele]<0.2){
	 overlapping= true;

       }

     }

   }

   //////////////////////////////

   //loop in all photons to calculate the percentage

   for(unsigned int Ipho=0; Ipho<photons.size(); Ipho++){

     TLorentzVector this_photon = photons[Ipho];
      percent_pho[Ipho] = abs((thisJetToClean.Pt()- this_photon.Pt())/(this_photon.Pt()));

   //loop in all photons and compare with the jet

    if (thisJetToClean.DeltaR(photons[Ipho])<0.4){

      if (percent_pho[Ipho]<0.2){
	overlapping= true;

      }


     }

    }


   /////////////////////////////


   for(unsigned int Imu=0; Imu<muons.size(); Imu++){

   TLorentzVector this_muon = muons[Imu];

   percent_mu[Imu] = abs((thisJetToClean.Pt()- this_muon.Pt())/(this_muon.Pt()));

    if (thisJetToClean.DeltaR(muons[Imu])<0.4 ){

      if (percent_mu[Imu]<0.2){
	overlapping= true;

       }

     }



   }



   for(unsigned int Ifake=0; Ifake<fakes.size(); Ifake++){

    TLorentzVector this_fake = fakes[Ifake];


    percent_fake[Ifake] = abs((thisJetToClean.Pt()- this_fake.Pt())/(this_fake.Pt()));

   if (thisJetToClean.DeltaR(fakes[Ifake])<0.4){
     if (percent_fake[Ifake]<0.2) {
       overlapping= true;

      }

     }


   }
   if(overlapping && weird) cout << "Removing jet " << Ijet << " due to object cleaning!" << endl;

   if(overlapping || thisJetToClean.Pt()<30 ) continue; //ths will make the code move to the next step

   jets.push_back(thisJetToClean);
   NJet_++;


 }
 sort(jets.begin(), jets.end(), P4SortCondition);

 h_NJet->SetBinContent(1,jets.size());

 for(unsigned int IJet=0; IJet<dirtyJets_n5.size(); IJet++){

  bool n5_overlapping = false;

  TLorentzVector this_n5JetToClean = dirtyJets_n5[IJet];
  float percent_pho[100];
  float percent_ele[100];
  float percent_mu[100];
  float percent_fake[100];

   std::pair<float,int> pair_n5JetPho = minDR_ValueAndIndex(this_n5JetToClean, photons);
   std::pair<float,int> pair_n5JetEle = minDR_ValueAndIndex(this_n5JetToClean, electrons);
   std::pair<float,int> pair_n5JetMu  = minDR_ValueAndIndex(this_n5JetToClean, muons);
   std::pair<float,int> pair_n5JetFake= minDR_ValueAndIndex(this_n5JetToClean, fakes);


   for (unsigned int Iele=0; Iele<electrons.size(); Iele++){

   TLorentzVector this_electron = electrons[Iele];

   percent_ele[Iele] = abs((this_n5JetToClean.Pt()- this_electron.Pt())/(this_electron.Pt()));


   if (this_n5JetToClean.DeltaR(electrons[Iele])<0.4){

   if (percent_ele[Iele]<0.2) n5_overlapping= true;

     else

     this_n5JetToClean = this_n5JetToClean - electrons[Iele];

     }

   }



   //////////////////////////////


   //loop in all photons and compare with the jet

   for (unsigned int Ipho=0; Ipho<photons.size(); Ipho++){

    TLorentzVector this_photon = photons[Ipho];

     percent_pho[Ipho] = abs((this_n5JetToClean.Pt()- this_photon.Pt())/(this_photon.Pt()));


    if (this_n5JetToClean.DeltaR(photons[Ipho])<0.4) {

      if( percent_pho[Ipho]<0.2) n5_overlapping= true;


      else

         this_n5JetToClean = this_n5JetToClean - photons[Ipho];

     }

    }



   /////////////////////////////



   for(unsigned int Imu=0; Imu<muons.size(); Imu++){

    TLorentzVector this_muon = muons[Imu];
    percent_mu[Imu] = abs((this_n5JetToClean.Pt()- this_muon.Pt())/(this_muon.Pt()));


    if (this_n5JetToClean.DeltaR(muons[Imu])<0.4){

    if(percent_mu[Imu]<0.2) n5_overlapping= true;

      else

      this_n5JetToClean = this_n5JetToClean - muons[Imu];


     }
   }



   for(unsigned int Ifake=0; Ifake<fakes.size(); Ifake++){

   TLorentzVector this_fake = fakes[Ifake];
   percent_fake[Ifake] = abs((this_n5JetToClean.Pt()- this_fake.Pt())/(this_fake.Pt()));


   if (this_n5JetToClean.DeltaR(fakes[Ifake])<0.4){

    if (percent_fake[Ifake]<0.2) n5_overlapping= true;


      this_n5JetToClean = this_n5JetToClean - fakes[Ifake];

     }


   }



  if(n5_overlapping || this_n5JetToClean.Pt()<30)continue; //ths will make the code move to the next step


  jets_n5.push_back(this_n5JetToClean);

  NJet_n5_++;


 }

  sort(jets_n5.begin(), jets_n5.end(), P4SortCondition);
  //////////////////////

  // loop over the photons and jets  and remove any PHOTONS that overal within dR<0.7//

  // loop over cleaned photons//
  vector<TLorentzVector> photonsCleaned;
  vector<TLorentzVector> photonsRemoved;

  int NPho_Cleaned=0;
  for(unsigned int iPho=0; iPho<photons.size(); iPho++){

   TLorentzVector thisPhotonToClean2 = photons[iPho];
   bool overlap=kFALSE;
   //loop over cleaned jets//
    for (unsigned int i=0; i<jets.size(); i++){

    if(thisPhotonToClean2.DeltaR(jets[i])>0.2 && thisPhotonToClean2.DeltaR(jets[i])<0.7) overlap =kTRUE;

   }

  //if(overlap) continue; //skip this photon and move on to the next one
  photonsCleaned.push_back(thisPhotonToClean2);
  NPho_Cleaned++;

  }

  sort(photonsCleaned.begin(), photonsCleaned.end(), P4SortCondition);


  //do the same for fakes ///

  vector<TLorentzVector> fakesCleaned;
  int NFakes_Cleaned=0;
  for(unsigned int iFake=0; iFake<fakes.size(); iFake++){

   TLorentzVector thisFakeToClean2 = fakes[iFake];
   bool overlap = kFALSE;
   //loop over cleaned jets//
    for (unsigned int i=0; i<jets.size(); i++){

	if(thisFakeToClean2.DeltaR(jets[i])<0.7) overlap = kTRUE;

    }

//  if(overlap) continue; //skip this photon and move onto the next one
  fakesCleaned.push_back(thisFakeToClean2);
  NFakes_Cleaned++;

  }

  sort(fakesCleaned.begin(), fakesCleaned.end(), P4SortCondition);



  /// do the same for electrons ///


  vector<TLorentzVector> electronsCleaned;
  int NElectrons_Cleaned=0;
  for(unsigned int iEle=0; iEle<electrons.size(); iEle++){

   TLorentzVector thisEleToClean2 = electrons[iEle];
   bool overlap = kFALSE;
   //loop over cleaned jets//
     for (unsigned int i=0; i<jets.size(); i++){

	if(thisEleToClean2.DeltaR(jets[i])>0.2 && thisEleToClean2.DeltaR(jets[i])<0.7) overlap = kTRUE;

   }

  //if(overlap) continue; //skip this photon and move onto the next one
  electronsCleaned.push_back(thisEleToClean2);
  NElectrons_Cleaned++;

  }

  sort(electronsCleaned.begin(), electronsCleaned.end(), P4SortCondition);




  //////////////
 // RECO MET //
 //////////////
 float mET  = pfMET;
 float mETx = pfMET*cos(pfMETPhi);
 float mETy = pfMET*sin(pfMETPhi);
 TLorentzVector met;

 met.SetPxPyPzE( mETx, mETy, 0., mET );




//Constract leading object pairs


 vector<TLorentzVector> allObjects;
 vector<TLorentzVector> allObjectsSorted;
 std::map<float,int>    unsortedObjectsMap;

 for(unsigned int i=0; i<photonsCleaned.size(); i++) 	allObjects.push_back(photonsCleaned[i]);
 for(unsigned int i=0; i<electronsCleaned.size(); i++) allObjects.push_back(electronsCleaned[i]);
 for(unsigned int i=0; i<fakesCleaned.size(); i++) 	allObjects.push_back(fakesCleaned[i]);

 if(allObjects.size()<2) return kFALSE; //require at least 2 good objects
 IVentry_++;

 for(unsigned int i=0; i<allObjects.size(); i++) unsortedObjectsMap[allObjects[i].Pt()] = i;

 allObjectsSorted = allObjects;
 sort(allObjectsSorted.begin(), allObjectsSorted.end(), P4SortCondition);

 unsigned int L = 3; //unsigned int L = allObjects.size();
 if(allObjects.size() < L) L = allObjects.size();

 bool  isGG = false, isEE = false, isFF = false, isEG = false, isEEGG = false, isEEGG_reco=false;

 float ggIDs[6]  = {111,112,113,131,311,110};
 float eeIDs[6]  = {222,221,223,232,322,220};
 float ffIDs[8]  = {133,233,333,331,313,332,323,330};
 float egIDs[12] = {121,122,123,132,212,211,213,231,312,321,120,210};
 float eventID = 0;
 float eventID2 = 0;

 for(unsigned int i=0; i<L; i++){
   unsigned int index = unsortedObjectsMap[allObjectsSorted[i].Pt()];
   int particleID = 0;

   if( index  <  photonsCleaned.size() )                                                       particleID = 1;
   if( (index >= photonsCleaned.size()) && (index < (photonsCleaned.size() + electronsCleaned.size())) ) {

   particleID = 2;
   //isElectron = true;
   }
   if( (index >= (photonsCleaned.size() + electronsCleaned.size())) && (index < allObjects.size()) )  particleID = 3;

   float factor = 100./(pow(10.,i));

   eventID += factor * particleID;


 }

 //////////////////////////////////////

 unsigned int L2 = 4;
 for(unsigned int i=0; i<L2; i++){

   if(allObjects.size() < 4) L2 = allObjects.size();

   unsigned int index = unsortedObjectsMap[allObjectsSorted[i].Pt()];
   int particleID = 0;

   if( index  <  photons.size() )                                                       particleID = 1;
   if( (index >= photons.size()) && (index < (photons.size() + electrons.size())) ) {

   particleID = 2;
   //isElectron = true;
   }
   if( (index >= (photons.size() + electrons.size())) && (index < allObjects.size()) )  particleID = 3;


   float factor = 1000./(pow(10.,i));
   eventID2 += factor *particleID ;

 }





/*
 if( std::find(ggIDs, ggIDs+6,  eventID) != ggIDs+6 )  isGG = true;
 if( std::find(eeIDs, eeIDs+6,  eventID) != eeIDs+6 )  isEE = true;
 if( std::find(ffIDs, ffIDs+8,  eventID) != ffIDs+8 )  isFF = true;
 if( std::find(egIDs, egIDs+12, eventID) != egIDs+12 ) isEG = true;

 *///had to changed that 15 October since we veto on Menglei's electrons and not ours

 if (eventID ==111 || eventID ==113 || eventID ==131 || eventID ==311 || eventID==110 || eventID==112) isGG=true;
 //if (eventID ==220 || eventID ==221 || eventID ==223 || eventID ==322 || eventID==122 ) isEE=true;
 if (eventID ==220 || eventID ==221 || eventID ==223 || eventID ==322 || eventID==232 || eventID==222 ) isEE=true;
 if (eventID ==333 || eventID ==133 || eventID ==331 || eventID ==330 || eventID==313 || eventID==332 || eventID==323 || eventID == 233 ) isFF=true;
 if (eventID==121 || eventID==120  || eventID ==210 || eventID ==123 || eventID==213 || eventID==312 || eventID==132 || eventID==321 || eventID == 211 || eventID==231 || eventID==122 || eventID==212) isEG=true;

 //construct eegg pairs //
 if (eventID2 ==2211) isEEGG = true;


 if (eventID2 ==2211 || eventID2==1122 || eventID2==1221  || eventID2==2112 || eventID2 == 2121 ) isEEGG_reco = true;

 if (isEEGG_reco ) recoEEGG++;


 if( !isGG && !isEE && !isFF && !isEG) return kFALSE;
 Ventry_++;


 /////////////////// MHT /////////////////////////

 float Ht = 0;
 float St = 0;
 TLorentzVector mMHTvec;
 TLorentzVector MHTvec;

  for(unsigned int j=0; j<jets_n5.size(); j++) mMHTvec = mMHTvec + jets_n5[j];
  for(unsigned int Ipho=0; Ipho<photons.size(); Ipho++) mMHTvec= mMHTvec +photons[Ipho];
  for(unsigned int Iele=0; Iele<electrons.size(); Iele++) mMHTvec = mMHTvec + electrons[Iele];
  for(unsigned int Imu=0; Imu<muons.size(); Imu++) mMHTvec = mMHTvec + muons[Imu];
  for(unsigned int Ifake=0; Ifake<fakes.size(); Ifake++) mMHTvec = mMHTvec + fakes[Ifake];


  MHTvec=-mMHTvec;

  /////HT////////////

  for (unsigned int j=0; j<jets_n5.size(); j++)  Ht = Ht + jets_n5[j].Pt();

  St= Ht;

  for(unsigned int Ipho=0; Ipho<photons.size(); Ipho++)   St = St + photons[Ipho].Pt();
  for(unsigned int Iele=0; Iele<electrons.size(); Iele++) St = St + electrons[Iele].Pt();
  for(unsigned int Imu=0; Imu<muons.size(); Imu++)        St = St + muons[Imu].Pt();
  for(unsigned int Ifake=0; Ifake<fakes.size(); Ifake++)  St = St + fakes[Ifake].Pt();



//Construct MET  ///


TLorentzVector mMETvec;
TLorentzVector METvec;

for(unsigned int j=0; j<jets_n5.size(); j++) mMETvec = mMETvec + jets_n5[j];
for(unsigned int Ipho=0; Ipho<photons.size(); Ipho++) mMETvec= mMETvec +photons[Ipho];
for(unsigned int Iele=0; Iele<electrons.size(); Iele++) mMETvec = mMETvec + electrons[Iele];
for(unsigned int Imu=0; Imu<muons.size(); Imu++) mMETvec = mMETvec + muons[Imu];
for(unsigned int Ifake=0; Ifake<fakes.size(); Ifake++) mMETvec = mMETvec + fakes[Ifake];


METvec=-mMETvec;


if ( event == 821151459 ) cout << "It survived!!! "<<endl;


 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //////////////////
 // Filling tree //
 //////////////////



if (MuNumber!=0) return kFALSE; //veto on muons after object cleaning

bool GGtrue = false;
bool EEtrue = false;



 if(isGG && (HLTPho>>14&1)  && (metFilters >>9&1) && (metFilters>>10&1)){



   //Take leading and subleading photons and their 4vector
   TLorentzVector leadpho    = photonsCleaned[0];
   TLorentzVector subleadpho = photonsCleaned[1];
   TLorentzVector dipho      = leadpho + subleadpho;


   //////////////////////////////////////////////////////////
  /////////////// MT2 construction ////////////////////////
  ////////////////////////////////////////////////////////
 //Calculate MT2. Set mass of photons and gravitinos to 0
 asymm_mt2_lester_bisect::disableCopyrightMessage();

 double mT2 = asymm_mt2_lester_bisect::get_mT2(0, leadpho.Px(), leadpho.Py(),
                                                   0, subleadpho.Px(), subleadpho.Py(),
                                                   met.Px(),met.Py(),
                                                   0,0,0);


   if((dipho.M()>105.) && (leadpho.DeltaR(subleadpho)>0.6)){


   for (unsigned int i=0; i<VetoElectrons.size(); i++){

        if (leadpho.DeltaR(VetoElectrons[i])>0.3 && subleadpho.DeltaR(VetoElectrons[i])>0.3) return kFALSE;

     }


     numGG++;
     GGtrue = true;



     Run_    	  = run;
     Event_  	  = event;
     isData_      = isData;
     Lumis_  	  = lumis;
     Rho_    	  = rho;
     //PUtrue_      = (*puTrue);
     //PUweight_    = puweight;
     MET_Filters_ = metFilters;
     MET_ 	  = pfMET;
     MET_Phi_ 	  = pfMETPhi;
     MET_sumEt_   = pfMETsumEt;
     MET_mEtSig_  = pfMETmEtSig;
     MET_Sig_ 	  = pfMETSig;

     MT2_=mT2;



     //reject events that  do NOT overlap with MengLei's electrons collection


     for(unsigned int i=0; i<photonsCleaned.size(); i++){

      Pho_E_[i] 	      = (*phoCalibE)[phoMap[photons[i].Pt()]];
      Pho_Et_[i] 	      = (*phoCalibEt)[phoMap[photons[i].Pt()]];
      Pho_Eta_[i] 	      = (*phoEta)[phoMap[photons[i].Pt()]];
      Pho_Phi_[i] 	      = (*phoPhi)[phoMap[photons[i].Pt()]];
      Pho_hasPixelSeed_[i]    = (*phohasPixelSeed)[phoMap[photons[i].Pt()]];
      Pho_EleVeto_[i] 	      = (*phoEleVeto)[phoMap[photons[i].Pt()]];
      Pho_R9_[i] 	      = (*phoR9)[phoMap[photons[i].Pt()]];
      Pho_HoverE_[i] 	      = (*phoHoverE)[phoMap[photons[i].Pt()]];
      Pho_SigmaIEtaIEta_[i]   = (*phoSigmaIEtaIEtaFull5x5)[phoMap[photons[i].Pt()]];
      Pho_SigmaIEtaIPhi_[i]   = (*phoSigmaIEtaIPhiFull5x5)[phoMap[photons[i].Pt()]];
      Pho_SigmaIPhiIPhi_[i]   = (*phoSigmaIPhiIPhiFull5x5)[phoMap[photons[i].Pt()]];
      Pho_IDMVA_[i] 	      = (*phoIDMVA)[phoMap[photons[i].Pt()]];
      Pho_IDbit_[i] 	      = (*phoIDbit)[phoMap[photons[i].Pt()]];
     }

     Dipho_Pt_        = dipho.Pt();
     Dipho_Mass_      = dipho.M();



      ST_         = St;
      MHT_	  = MHTvec.Pt();
      MHTPhi_ 	  = MHTvec.Phi();




    //if (MET_>6 and MET_< 9) myfile <<event <<"\n";

     ggNvertex ->Fill(nVtx);
     ggNJets   ->Fill(nJet);



     ggDiEmPt->Fill(dipho.Pt());


     GGtree_->Fill();

     if (MET_<50){

     gg_sample++;

     }


     if (MET_<50 && NJet_>=2){

     gg_sample_nj++;

     }



     if (MET_<50 && St>=400){

     gg_sample_st400++;

     }



     if (MET_ <50 && NJet_>=2 && St>80 && St<150) {
     gg_sample_bin1++;
    }


    if (MET_ <50 && NJet_>=2 && St>150 && St<250) {
     gg_sample_bin2++;
    }


    if (MET_ <50 && NJet_>=2 && St>250 && St<400) {
     gg_sample_bin3++;
    }



    if (MET_ <50 && NJet_>=2 && St>=400 && St<999999) {
     gg_sample_bin4++;
    }


    if (event == 104524039 || event== 733594020 || event== 691908143) cout <<"event with 1 muon passed " << endl;

    //if (MET_<100) DiEmPTVSNJets_gg ->Fill(Dipho_Pt_, NJet_);

    //if (MET_<100) DiEmPTVSNJetsVarbin_gg ->Fill(Dipho_Pt_, NJet_);



   DiEmPTVSNJets_gg ->Fill(Dipho_Pt_, NJet_);

   DiEmPTVSNJetsVarbin_gg ->Fill(Dipho_Pt_, NJet_);

   h_gg_MET->Fill(MET_);

   }//diphoMass end
 }//gg end  event!= 821151459-->event with 4 electron




 if(isEE && (HLTPho>>15&1)  && (metFilters >>9&1) && (metFilters>>10&1)){


   TLorentzVector leade      = electronsCleaned[0];
   TLorentzVector subleade   = electronsCleaned[1];
   TLorentzVector die        = leade + subleade;

   //////////////////////////////////////////////////////////
  /////////////// MT2 construction ////////////////////////
  ////////////////////////////////////////////////////////
 //Calculate MT2. Set mass of photons and gravitinos to 0
 asymm_mt2_lester_bisect::disableCopyrightMessage();

 double mT2 = asymm_mt2_lester_bisect::get_mT2(0, leade.Px(), leade.Py(),
                                                   0, subleade.Px(), subleade.Py(),
                                                   met.Px(),met.Py(),
                                                   0,0,0);




   if((die.M()>75.) && (die.M()<105.) && (leade.DeltaR(subleade)>0.6)){



     //reject events that  do NOT overlap with MengLei's electrons collection


      for (unsigned int i=0; i<VetoElectrons.size(); i++){

        if (leade.DeltaR(VetoElectrons[i])>0.3 && subleade.DeltaR(VetoElectrons[i])>0.3) return kFALSE;

     }


     numEE++;

     Run_    	  = run;
     Event_  	  = event;
     isData_      = isData;
     Lumis_  	  = lumis;
     Rho_    	  = rho;


     MET_Filters_ 	= metFilters;
     MET_ = pfMET;
     MET_Phi_ = pfMETPhi;
     MET_sumEt_ = pfMETsumEt;
     MET_mEtSig_ = pfMETmEtSig;
     MET_Sig_ = pfMETSig;

     MT2_=mT2;
     //MT_ =mT;


     for(unsigned int i=0; i<electronsCleaned.size(); i++){



       Ele_Pt_[i] 	= (*phoCalibEt)[eleMap[electrons[i].Pt()]];
       Ele_Eta_[i] = (*phoEta)[eleMap[electrons[i].Pt()]];
       Ele_Phi_[i] = (*phoPhi)[eleMap[electrons[i].Pt()]];
       Ele_SigmaIEtaIEta_[i] = (*phoSigmaIEtaIEtaFull5x5)[eleMap[electrons[i].Pt()]];
       Ele_SigmaIEtaIPhi_[i] = (*phoSigmaIEtaIPhiFull5x5)[eleMap[electrons[i].Pt()]];
       Ele_SigmaIPhiIPhi_[i] = (*phoSigmaIPhiIPhiFull5x5)[eleMap[electrons[i].Pt()]];
       Ele_SigmaIEtaIEtaFull5x5_[i] = (*phoSigmaIEtaIEtaFull5x5)[eleMap[electrons[i].Pt()]];
       Ele_SigmaIPhiIPhiFull5x5_[i] = (*phoSigmaIPhiIPhiFull5x5)[eleMap[electrons[i].Pt()]];
     }

     Dipho_Pt_ = die.Pt();
     Dipho_Mass_ = die.M();

     Run_    	  = run;
     Event_  	  = event;
     isData_      = isData;
     Lumis_  	  = lumis;
     Rho_    	  = rho;


     ST_	    = St;
     MHT_	    = MHTvec.Pt();
     MHTPhi_        = MHTvec.Phi();


     eeDiEmPt->Fill(die.Pt());

     eeNvertex ->Fill(nVtx);
     eeNJets   ->Fill(nJet);

     h_E_MET->Fill(MET_);
     h_E_METpred->Fill(MET_);

     //if (MET_<100) DiEmPTVSNJets_ee ->Fill(Dipho_Pt_, NJet_);


	DiEmPTVSNJets_ee ->Fill(Dipho_Pt_, NJet_);

	DiEmPTVSNJetsVarbin_ee ->Fill(Dipho_Pt_, NJet_);

     h_ee_MET->Fill(MET_);

     EEtree_->Fill();

    if (MET_>50 and MET_<60) myfile <<event <<"\n";



     if (MET_<50){

     ee_sample++;

     }


     if (MET_<50 && NJet_>=2){

     ee_sample_nj++;

     }



     if (MET_<50 && St>=400){

     ee_sample_st400++;

     }



     if (MET_ <50 && NJet_>=2 && St>80 && St<150) {
     ee_sample_bin1++;
    }


     if (MET_ <50 && NJet_>=2 && St>150 && St<250) {
     ee_sample_bin2++;
    }


    if (MET_ <50 && NJet_>=2 && St>250 && St<400) {
     ee_sample_bin3++;
    }



    if (MET_ <50 && NJet_>=2 && St>=400 && St<999999) {
     ee_sample_bin4++;
    }


    if (MuNumber>=1) cout << "event with 1 or more muons passed  " << event << " " << MuNumber << endl;

   }//diEMass end


  }//ee end


 if (isEG && (HLTPho>>14&1)  && (metFilters >>9&1) && (metFilters>>10&1)){

 TLorentzVector Electron = electronsCleaned[0];
 TLorentzVector Photon   = photonsCleaned[0];
 TLorentzVector EGamma   = Electron + Photon;


 if (event==3875836596) {

 cout << "Electron Pt " << Electron.Pt() << "Photon Pt " << Photon.Pt() << endl;
 cout << "Number of electrons " << NEle_ << " Number of photons " << NPho_ << "number of fakes " << NFake_ <<endl;
 cout << "diPair mass "<< EGamma.M()<<endl;



 }



   if((EGamma.M()>105.) && (Photon.DeltaR(Electron)>0.6)){


    for (unsigned int i=0; i<VetoElectrons.size(); i++){

        if (Photon.DeltaR(VetoElectrons[i])>0.3 && Electron.DeltaR(VetoElectrons[i])>0.3) return kFALSE;

     }


     numEG++;

//construct MT variable

 float mT=0;

  mT = TMath::Sqrt(2*Electron.Pt()*pfMET*(1-TMath::Cos( Electron.Phi() - pfMETPhi )));


 //if (mT>100) return kFALSE;

 if (event==54835416) cout << "Event Passed with " << EGamma.M()<<endl;


     Run_    	  = run;
     Event_  	  = event;
     isData_      = isData;
     Lumis_  	  = lumis;
     Rho_    	  = rho;



     MET_Filters_ = metFilters;
     MET_ = pfMET;
     MET_Phi_ = pfMETPhi;
     MET_sumEt_ = pfMETsumEt;
     MET_mEtSig_ = pfMETmEtSig;
     MET_Sig_ = pfMETSig;


     ST_	    = St;
     MHT_	    = MHTvec.Pt();
     MHTPhi_        = MHTvec.Phi();

     MT_ = mT;



    if (mT<100){
    egMEtMTCut    ->Fill(MET_);
    egInvMassMTCut->Fill(EGamma.M());

   }



 //DeltaPhi between lepton and MET

 if ( abs(Electron.Phi()-pfMETPhi) > 0.3 && mT<100){

 egInvMassCutDeltaPhiCut->Fill(EGamma.M());
 egMEtMtDeltaPhiCut->Fill(MET_);

 }



egMt     ->Fill(mT);
egInvMassnoMTCut ->Fill(EGamma.M());


 h_EG_MET->Fill(MET_); //no MT or DeltaPhi Cut

 EGtree_->Fill();

 }


 }


 if(isFF && (HLTPho>>14&1) && (metFilters >>9&1) && (metFilters>>10&1)){


   TLorentzVector leadfake    = fakesCleaned[0];
   TLorentzVector subleadfake = fakesCleaned[1];
   TLorentzVector difake      = leadfake + subleadfake;

   //////////////////////////////////////////////////////////
  /////////////// MT2 construction ////////////////////////
  ////////////////////////////////////////////////////////
 //Calculate MT2. Set mass of photons and gravitinos to 0
 asymm_mt2_lester_bisect::disableCopyrightMessage();

 double mT2 = asymm_mt2_lester_bisect::get_mT2(0, leadfake.Px(), leadfake.Py(),
                                                   0, subleadfake.Px(), subleadfake.Py(),
                                                   met.Px(),met.Py(),
                                                   0,0,0);


   if((difake.M()>105.) && (leadfake.DeltaR(subleadfake)>0.6)){


     //reject events that  do NOT overlap with MengLei's electrons collection


      for (unsigned int i=0; i<VetoElectrons.size(); i++){

        if (leadfake.DeltaR(VetoElectrons[i])>0.3 && subleadfake.DeltaR(VetoElectrons[i])>0.3) return kFALSE;

     }

      numFF++;


     Run_    	  = run;
     Event_  	  = event;
     isData_      = isData;
     Lumis_  	  = lumis;
     Rho_    	  = rho;



     MET_Filters_ 	= metFilters;
     MET_ = pfMET;
     MET_Phi_ = pfMETPhi;
     MET_sumEt_ = pfMETsumEt;
     MET_mEtSig_ = pfMETmEtSig;
     MET_Sig_ = pfMETSig;

     MT2_=mT2;



     for(unsigned int i=0; i<fakesCleaned.size(); i++){
      Pho_E_[i] 		    = (*phoCalibE)[fakeMap[fakes[i].Pt()]];
      Pho_Et_[i] 		    = (*phoCalibEt)[fakeMap[fakes[i].Pt()]];
      Pho_Eta_[i] 	    = (*phoEta)[fakeMap[fakes[i].Pt()]];
      Pho_Phi_[i] 	    = (*phoPhi)[fakeMap[fakes[i].Pt()]];
      Pho_hasPixelSeed_[i]    = (*phohasPixelSeed)[fakeMap[fakes[i].Pt()]];
      Pho_EleVeto_[i] 	    = (*phoEleVeto)[fakeMap[fakes[i].Pt()]];
      Pho_R9_[i] 		    = (*phoR9)[fakeMap[fakes[i].Pt()]];
      Pho_HoverE_[i] 	    = (*phoHoverE)[fakeMap[fakes[i].Pt()]];
      Pho_SigmaIEtaIEta_[i]   = (*phoSigmaIEtaIEtaFull5x5)[fakeMap[fakes[i].Pt()]];
      Pho_SigmaIEtaIPhi_[i]   = (*phoSigmaIEtaIPhiFull5x5)[fakeMap[fakes[i].Pt()]];
      Pho_SigmaIPhiIPhi_[i]   = (*phoSigmaIPhiIPhiFull5x5)[fakeMap[fakes[i].Pt()]];




     }
     Dipho_Pt_ = difake.Pt();

     ffDiEmPt->Fill(difake.Pt());

     ffNvertex ->Fill(nVtx);
     ffNJets   ->Fill(nJet);

     ST_	    = St;
     MHT_	    = MHTvec.Pt();
     MHTPhi_        = MHTvec.Phi();




   DiEmPTVSNJets_ff ->Fill(Dipho_Pt_, NJet_);

   DiEmPTVSNJetsVarbin_ff ->Fill(Dipho_Pt_, NJet_);

   h_ff_MET->Fill(MET_);


     FFtree_->Fill();

     if (MET_<50){

     ff_sample++;

     }


     if (MET_<50 && NJet_>=2){

     ff_sample_nj++;

     }



     if (MET_<50 && St>=400){

     ff_sample_st400++;

     }



     if (MET_ <50 && NJet_>=2 && St>80 && St<150) {
     ff_sample_bin1++;
     }

     if (MET_ <50 && NJet_>=2 && St>150 && St<250) {
     ff_sample_bin2++;
    }


    if (MET_ <50 && NJet_>=2 && St>250 && St<400) {
     ff_sample_bin3++;
    }



    if (MET_ <50 && NJet_>=2 && St>=400 && St<999999) {
     ff_sample_bin4++;
    }


   }//diFMass end

}//ff end

  ////////////////////
 ///////////////////
 // clear vectors //
 ///////////////////
 phoMap.clear();

 eleMap.clear();
 muMap.clear();
 jetMap.clear();
 fakeMap.clear();
 unsortedObjectsMap.clear();

 dirtyPhotons.clear();
 dirtyElectrons.clear();
 dirtyJets.clear();
 dirtyJets_n5.clear();
 jets_n5.clear();
 jets.clear();
 muons.clear();
 electrons.clear();
 photons.clear();
 fakes.clear();


 photonsCleaned.clear();
 fakesCleaned.clear();
 electronsCleaned.clear();



 allObjects.clear();
 allObjectsSorted.clear();



   return kTRUE;
}
////////////////////////////////////////////
void T5gg_analyzer::SlaveTerminate()
{}
//////////////////////////////////////////

//////////////////////////////////////////
void T5gg_analyzer::Terminate()
{

   cout<<"\nFinished.\n--Summary: "<<Ientry_<<" total analyzed events.\n"<<endl;
   cout << "Num gg = "<< numGG << endl;
   cout << "Num ff = " << numFF << endl;
   cout << "Num ee = " << numEE << endl;
   cout << "Num eg = " << numEG << endl;
   cout << "Num eegg= " << numEEGG<<endl;
   cout << "Reco EEGG= " << recoEEGG<<endl;

 h_events->SetBinContent(1,Ientry_);

 gg_events->SetBinContent(1,gg_sample);
 ee_events->SetBinContent(1,ee_sample);
 ff_events->SetBinContent(1,ff_sample);
 eeggEvents->SetBinContent(1,numEEGG);

 gg_events_nj->SetBinContent(1,gg_sample_nj);
 ee_events_nj->SetBinContent(1,ee_sample_nj);
 ff_events_nj->SetBinContent(1,ff_sample_nj);

 gg_events_st400->SetBinContent(1,gg_sample_st400);
 ee_events_st400->SetBinContent(1,ee_sample_st400);
 ff_events_st400->SetBinContent(1,ff_sample_st400);

 gg_events_STbin1->SetBinContent(1,gg_sample_bin1);
 ee_events_STbin1->SetBinContent(1,ee_sample_bin1);
 ff_events_STbin1->SetBinContent(1,ff_sample_bin1);

 gg_events_STbin2->SetBinContent(1,gg_sample_bin2);
 ee_events_STbin2->SetBinContent(1,ee_sample_bin2);
 ff_events_STbin2->SetBinContent(1,ff_sample_bin2);

 gg_events_STbin3->SetBinContent(1,gg_sample_bin3);
 ee_events_STbin3->SetBinContent(1,ee_sample_bin3);
 ff_events_STbin3->SetBinContent(1,ff_sample_bin3);

 gg_events_STbin4->SetBinContent(1,gg_sample_bin4);
 ee_events_STbin4->SetBinContent(1,ee_sample_bin4);
 ff_events_STbin4->SetBinContent(1,ff_sample_bin4);


 fout.cd();
 fout.Write();
 fout.Close();

 myfile.close();

}

//////////////////////////////////////////


TLorentzVector T5gg_analyzer::p4vec(float pt, float phi, float pz)
{
 TLorentzVector thisp4vec;
 thisp4vec.SetPxPyPzE( pt*cos(phi), pt*sin(phi), pz, sqrt(pt*pt + pz*pz) );

 return thisp4vec;
}


///////////////////////////////////////////////////////
std::pair<float,int> T5gg_analyzer::minDR_ValueAndIndex(TLorentzVector lv ,vector<TLorentzVector> v)
{
 float minDR      = 1000.;
 int   drMinIndex = 1000;
 float percent[50];
 for(unsigned int i=0; i<v.size(); i++){
   if(lv.DeltaR(v[i]) < minDR){
     minDR = lv.DeltaR(v[i]);
     drMinIndex = i;
     percent[i] = abs((lv.Pt()-v[i].Pt())/(v[i].Pt())) ;

   }
 }


 return std::make_pair(minDR,drMinIndex);

}

/*
std::pair<float int> T5gg_analyzer::DR_ValueAndIndex(TLorentzVector lv ,vector<TLorentzVector> v)
{
  float DRvalue;
  float DRIndex;
  for(unsigned int i=0; i<v.size(); i++){
   if lv.D
   DRValue=lv.DeltaR(v[i]);
   DRIndex=i;


  }


}
*/

///////////////////////////////////////////////////////
int T5gg_analyzer::IsGoodPhoton(int iPho, int iSel, int verbose)
{
 if(verbose!=0) cout<<"--Inside IsGoodPhoton..."<<endl;
 int photonIdentity = 9;

 float  phEt             = (*phoCalibEt)[iPho];
 float  phSCEta          = (*phoSCEta)[iPho];
 float  phSinin          = (*phoSigmaIEtaIEtaFull5x5)[iPho];
 float  phHoE            = (*phoHoverE)[iPho];
 //photon ID cuts EB, EE     Loose; iSel=0			       Medium; iSel=1			        Tight; iSel=2
 float  HoE_EB[3]   	 = { 0.0597,  				       0.0396,				        0.00269					};
 float  Sinin_EB[3] 	 = { 0.01031, 				       0.01022,  			        0.00994 				};
 float  HoE_EE[3]   	 = { 0.0481,  				       0.0219,				        0.0213					};
 float  Sinin_EE[3] 	 = { 0.03013, 				       0.03001,  			        0.03000  				};
 //photon isolation cuts Barrel
 double CHiso_EB[3]   	 = { 1.295,				       0.441,				        0.202					};
 double NHiso_EB[3]   	 = { 10.910 + 0.0148*phEt + 0.000017*phEt*phEt,2.725 + 0.0148*phEt + 0.000017*phEt*phEt,0.264 + 0.0148*phEt + 0.000017*phEt*phEt};
 double PHiso_EB[3]   	 = { 3.630 + 0.0047*phEt,   	     	       2.571 + 0.0047*phEt,		        2.362 + 0.0047*phEt			};
 //photon isolation cuts Endcap
 double CHiso_EE[3]   	 = { 1.011, 			     	       0.442,				        0.034					};
 double NHiso_EE[3]   	 = { 5.931 + 0.0163*phEt + 0.000014*phEt*phEt, 1.715 + 0.0163*phEt + 0.000014*phEt*phEt,0.586 + 0.0163*phEt + 0.000014*phEt*phEt};
 double PHiso_EE[3]   	 = { 6.641 + 0.0034*phEt, 		       3.863 + 0.0034*phEt,		        2.617 + 0.0034*phEt			};
 //Effective areas
 float  chEA,nhEA,phEA,chIsoCor,nhIsoCor,phIsoCor;
 if      (fabs(phSCEta) < 1.0)  { chEA = EffAreas[0][0]; nhEA = EffAreas[0][1]; phEA = EffAreas[0][2]; }
 else if (fabs(phSCEta) < 1.479){ chEA = EffAreas[1][0]; nhEA = EffAreas[1][1]; phEA = EffAreas[1][2]; }
 else if (fabs(phSCEta) < 2.0)  { chEA = EffAreas[2][0]; nhEA = EffAreas[2][1]; phEA = EffAreas[2][2]; }
 else if (fabs(phSCEta) < 2.2)  { chEA = EffAreas[3][0]; nhEA = EffAreas[3][1]; phEA = EffAreas[3][2]; }
 else if (fabs(phSCEta) < 2.3)  { chEA = EffAreas[4][0]; nhEA = EffAreas[4][1]; phEA = EffAreas[4][2]; }
 else if (fabs(phSCEta) < 2.4)  { chEA = EffAreas[5][0]; nhEA = EffAreas[5][1]; phEA = EffAreas[5][2]; }
 else                           { chEA = EffAreas[6][0]; nhEA = EffAreas[6][1]; phEA = EffAreas[6][2]; }
 //isolation variables correction
 chIsoCor = (*phoPFChIso)[iPho]  - rho*chEA;
 nhIsoCor = (*phoPFNeuIso)[iPho] - rho*nhEA;
 phIsoCor = (*phoPFPhoIso)[iPho] - rho*phEA;

 //photon ID for Barrel
 if(fabs(phSCEta) < nF1){
   if(    (phHoE                          < HoE_EB[iSel])
       && (phSinin                        < Sinin_EB[iSel])
       && (std::max(chIsoCor, float(0.0)) < CHiso_EB[iSel])
       && (std::max(nhIsoCor, float(0.0)) < NHiso_EB[iSel])
       && (std::max(phIsoCor, float(0.0)) < PHiso_EB[iSel])
      ) photonIdentity = 1;
   else photonIdentity = 0;
 }
 //photon ID for Endcap
 else if((fabs(phSCEta) > nF2) && (fabs(phSCEta) < nF3)){
   if(    (phHoE                          < HoE_EE[iSel])
       && (phSinin                        < Sinin_EE[iSel])
       && (std::max(chIsoCor, float(0.0)) < CHiso_EE[iSel])
       && (std::max(nhIsoCor, float(0.0)) < NHiso_EE[iSel])
       && (std::max(phIsoCor, float(0.0)) < PHiso_EE[iSel])
      ) photonIdentity = 2;
   else photonIdentity = 0;
 }
 else photonIdentity = 0;
 //fake ID for Barrel
 if(fabs(phSCEta) < nF1){
   if(    (phHoE                          < 0.0396)
       && (std::max(nhIsoCor, float(0.0)) < 15)
       && (std::max(phIsoCor, float(0.0)) < 15)
       && (    ( (phSinin                 > 0.01022) || ( (std::max(chIsoCor,float(0.0)) > 0.441) ))
       && (std::max(chIsoCor,float(0.0))  < 25.0)  )
       && (phSinin                        < 0.015)
       //&& (std::max(chIsoCor, float(0.0)) > CHiso_EB[iSel])


      ) photonIdentity = -1;
 }

 //so if a photon passed photon ID it gives 1(EB), or 2(EE)
 //if it did not pass it gives 0 unless it is EB fake, so gives -1
 return photonIdentity;

}

//////////////////////////////////


 void T5gg_analyzer::IncludeAJson(std::string jsonfile) {

 // Fairly primitive brute force json parser -- opens the json file named in the argument
 // and adds that to the goodrunlumilist map.  Overlapping jsons are merged inclusively.

 char thing;


 ifstream jsonInput;

 std::cout << "Sucking in Json file: " << jsonfile << " which includes: " << std::endl;

 jsonInput.open(jsonfile.c_str());

 if (!jsonInput.good()) {
   std::cout << "Problem reading Json file...  Didn't suck anything in... " << std::endl;
   return;
 }
 jsonInput.get(thing);
 while (jsonInput.good()) {
   if (thing=='{') {  // start of list
     while (thing != '}') {
       int runnum;
       if (thing == '"') {
         std::string srunnum;
         jsonInput.get(thing); // get stuff inside ""

         while (thing != '"') {
           srunnum+=thing; // get stuff inside ""
           jsonInput.get(thing);

         }
         sscanf(srunnum.c_str(),"%i",&runnum);
         std::cout << " runnum: " << runnum << std::endl;
         bool newrun=true;


       } // inside ""
       if (thing == '[') {
         jsonInput.get(thing); // get stuff inside []
         while (thing != ']') {
           if (thing == '[') {
             jsonInput.get(thing); // get stuff inside series []

             std::string lumiseries;
             int firstlumi,lastlumi;
             while (thing !=']') {
               lumiseries+=thing;
               jsonInput.get(thing); // get stuff inside series []
             }
             sscanf(lumiseries.c_str(),"%i,%i",&firstlumi,&lastlumi);
             std::cout << "  lumis  " << firstlumi << " to " << lastlumi << std::endl;

             // At this point have runnum, first lumi, last lumi -- so can fill map here...
             for (int l=firstlumi;l<=lastlumi;l++) {
               goodrunlumilist[runnum][l]=true;
	       //cout<< "filled "<< goodrunlumilist[runnum][l] <<endl;
             }

           } // inside actual series []
           jsonInput.get(thing); // get stuff inside []
         }
       } // inside []
       jsonInput.get(thing); // get another char looking for "
     }
   } // inside {}
   jsonInput.get(thing); // get another char looking for {
 } // EOF
 jsonInput.close();



 //return goodrunlumilist;

 //if (goodrunlumilist[run][lumi]) return true;


}



bool T5gg_analyzer::isInJson(Int_t run_, Int_t lumi){

  if (goodrunlumilist[run_][lumi]) return true;

   return false;

}


/////////////////////////////

void T5gg_analyzer::GetAllEntries(Long64_t entry)
{// List of branches

b_run->GetEntry(entry);
b_event->GetEntry(entry);
b_lumis->GetEntry(entry);
b_isData->GetEntry(entry);
b_nVtx->GetEntry(entry);
b_nGoodVtx->GetEntry(entry);
b_nTrksPV->GetEntry(entry);
b_isPVGood->GetEntry(entry);
b_vtx->GetEntry(entry);
b_vty->GetEntry(entry);
b_vtz->GetEntry(entry);
b_rho->GetEntry(entry);
b_rhoCentral->GetEntry(entry);
b_HLTEleMuX->GetEntry(entry);
b_HLTPho->GetEntry(entry);
b_HLTJet->GetEntry(entry);
b_HLTEleMuXIsPrescaled->GetEntry(entry);
b_HLTPhoIsPrescaled->GetEntry(entry);
b_HLTJetIsPrescaled->GetEntry(entry);
b_phoPrescale->GetEntry(entry);
b_metFilters->GetEntry(entry);
b_pfMET->GetEntry(entry);
b_pfMETPhi->GetEntry(entry);
b_pfMETsumEt->GetEntry(entry);
b_pfMETmEtSig->GetEntry(entry);
b_pfMETSig->GetEntry(entry);
b_pfMET_T1JERUp->GetEntry(entry);
b_pfMET_T1JERDo->GetEntry(entry);
b_pfMET_T1JESUp->GetEntry(entry);
b_pfMET_T1JESDo->GetEntry(entry);
b_pfMET_T1UESUp->GetEntry(entry);
b_pfMET_T1UESDo->GetEntry(entry);
b_pfMETPhi_T1JESUp->GetEntry(entry);
b_pfMETPhi_T1JESDo->GetEntry(entry);
b_pfMETPhi_T1UESUp->GetEntry(entry);
b_pfMETPhi_T1UESDo->GetEntry(entry);
b_nPho->GetEntry(entry);
b_phoE->GetEntry(entry);
b_phoEt->GetEntry(entry);
b_phoEta->GetEntry(entry);
b_phoPhi->GetEntry(entry);
b_phoCalibE->GetEntry(entry);
b_phoCalibEt->GetEntry(entry);
b_phoSCE->GetEntry(entry);
b_phoSCRawE->GetEntry(entry);
b_phoESEn->GetEntry(entry);
b_phoESEnP1->GetEntry(entry);
b_phoESEnP2->GetEntry(entry);
b_phoSCEta->GetEntry(entry);
b_phoSCPhi->GetEntry(entry);
b_phoSCEtaWidth->GetEntry(entry);
b_phoSCPhiWidth->GetEntry(entry);
b_phoSCBrem->GetEntry(entry);
b_phohasPixelSeed->GetEntry(entry);
b_phoEleVeto->GetEntry(entry);
b_phoR9->GetEntry(entry);
b_phoHoverE->GetEntry(entry);
b_phoE1x3->GetEntry(entry);
b_phoE1x5->GetEntry(entry);
b_phoE2x2->GetEntry(entry);
b_phoE2x5Max->GetEntry(entry);
b_phoE5x5->GetEntry(entry);
b_phoESEffSigmaRR->GetEntry(entry);
b_phoSigmaIEtaIEtaFull5x5->GetEntry(entry);
b_phoSigmaIEtaIPhiFull5x5->GetEntry(entry);
b_phoSigmaIPhiIPhiFull5x5->GetEntry(entry);
b_phoE1x3Full5x5->GetEntry(entry);
b_phoE1x5Full5x5->GetEntry(entry);
b_phoE2x2Full5x5->GetEntry(entry);
b_phoE2x5MaxFull5x5->GetEntry(entry);
b_phoE5x5Full5x5->GetEntry(entry);
b_phoR9Full5x5->GetEntry(entry);
b_phoPFChIso->GetEntry(entry);
b_phoPFPhoIso->GetEntry(entry);
b_phoPFNeuIso->GetEntry(entry);
b_phoPFChWorstIso->GetEntry(entry);
b_phoCITKChIso->GetEntry(entry);
b_phoCITKPhoIso->GetEntry(entry);
b_phoCITKNeuIso->GetEntry(entry);
b_phoIDMVA->GetEntry(entry);
b_phoFiredSingleTrgs->GetEntry(entry);
b_phoFiredDoubleTrgs->GetEntry(entry);
b_phoFiredL1Trgs->GetEntry(entry);
b_phoSeedTime->GetEntry(entry);
b_phoSeedEnergy->GetEntry(entry);
b_phoxtalBits->GetEntry(entry);
b_phoIDbit->GetEntry(entry);
b_npfPho->GetEntry(entry);
b_pfphoEt->GetEntry(entry);
b_pfphoEta->GetEntry(entry);
b_pfphoPhi->GetEntry(entry);
b_nEle->GetEntry(entry);
b_eleCharge->GetEntry(entry);
b_eleChargeConsistent->GetEntry(entry);
b_eleEn->GetEntry(entry);
b_eleSCEn->GetEntry(entry);
b_eleESEn->GetEntry(entry);
b_eleESEnP1->GetEntry(entry);
b_eleESEnP2->GetEntry(entry);
b_eleD0->GetEntry(entry);
b_eleDz->GetEntry(entry);
b_eleSIP->GetEntry(entry);
b_elePt->GetEntry(entry);
b_eleEta->GetEntry(entry);
b_elePhi->GetEntry(entry);
b_eleR9->GetEntry(entry);
b_eleCalibPt->GetEntry(entry);
b_eleCalibEn->GetEntry(entry);
b_eleSCEta->GetEntry(entry);
b_eleSCPhi->GetEntry(entry);
b_eleSCRawEn->GetEntry(entry);
b_eleSCEtaWidth->GetEntry(entry);
b_eleSCPhiWidth->GetEntry(entry);
b_eleHoverE->GetEntry(entry);
b_eleEoverP->GetEntry(entry);
b_eleEoverPout->GetEntry(entry);
b_eleEoverPInv->GetEntry(entry);
b_eleBrem->GetEntry(entry);
b_eledEtaAtVtx->GetEntry(entry);
b_eledPhiAtVtx->GetEntry(entry);
b_eledEtaAtCalo->GetEntry(entry);
b_eleSigmaIEtaIEtaFull5x5->GetEntry(entry);
b_eleSigmaIPhiIPhiFull5x5->GetEntry(entry);
b_eleConvVeto->GetEntry(entry);
b_eleMissHits->GetEntry(entry);
b_eleESEffSigmaRR->GetEntry(entry);
b_elePFChIso->GetEntry(entry);
b_elePFPhoIso->GetEntry(entry);
b_elePFNeuIso->GetEntry(entry);
b_elePFPUIso->GetEntry(entry);
b_elePFClusEcalIso->GetEntry(entry);
b_elePFClusHcalIso->GetEntry(entry);
b_elePFMiniIso->GetEntry(entry);
b_eleIDMVA->GetEntry(entry);
b_eleIDMVAHZZ->GetEntry(entry);
b_eledEtaseedAtVtx->GetEntry(entry);
b_eleE1x5->GetEntry(entry);
b_eleE2x5->GetEntry(entry);
b_eleE5x5->GetEntry(entry);
b_eleE1x5Full5x5->GetEntry(entry);
b_eleE2x5Full5x5->GetEntry(entry);
b_eleE5x5Full5x5->GetEntry(entry);
b_eleR9Full5x5->GetEntry(entry);
b_eleEcalDrivenSeed->GetEntry(entry);
b_eleDr03EcalRecHitSumEt->GetEntry(entry);
b_eleDr03HcalDepth1TowerSumEt->GetEntry(entry);
b_eleDr03HcalDepth2TowerSumEt->GetEntry(entry);
b_eleDr03HcalTowerSumEt->GetEntry(entry);
b_eleDr03TkSumPt->GetEntry(entry);
b_elecaloEnergy->GetEntry(entry);
b_eleTrkdxy->GetEntry(entry);
b_eleKFHits->GetEntry(entry);
b_eleKFChi2->GetEntry(entry);
b_eleGSFChi2->GetEntry(entry);
b_eleGSFPt->GetEntry(entry);
b_eleGSFEta->GetEntry(entry);
b_eleGSFPhi->GetEntry(entry);
b_eleGSFCharge->GetEntry(entry);
b_eleGSFHits->GetEntry(entry);
b_eleGSFMissHits->GetEntry(entry);
b_eleGSFNHitsMax->GetEntry(entry);
b_eleGSFVtxProb->GetEntry(entry);
b_eleGSFlxyPV->GetEntry(entry);
b_eleGSFlxyBS->GetEntry(entry);
b_eleBCEn->GetEntry(entry);
b_eleBCEta->GetEntry(entry);
b_eleBCPhi->GetEntry(entry);
b_eleBCS25->GetEntry(entry);
b_eleBCS15->GetEntry(entry);
b_eleBCSieie->GetEntry(entry);
b_eleBCSieip->GetEntry(entry);
b_eleBCSipip->GetEntry(entry);
b_eleFiredSingleTrgs->GetEntry(entry);
b_eleFiredDoubleTrgs->GetEntry(entry);
b_eleFiredL1Trgs->GetEntry(entry);
b_eleIDbit->GetEntry(entry);
b_npfHF->GetEntry(entry);
b_pfHFEn->GetEntry(entry);
b_pfHFECALEn->GetEntry(entry);
b_pfHFHCALEn->GetEntry(entry);
b_pfHFPt->GetEntry(entry);
b_pfHFEta->GetEntry(entry);
b_pfHFPhi->GetEntry(entry);
b_pfHFIso->GetEntry(entry);
b_nMu->GetEntry(entry);
b_muPt->GetEntry(entry);
b_muEn->GetEntry(entry);
b_muEta->GetEntry(entry);
b_muPhi->GetEntry(entry);
b_muCharge->GetEntry(entry);
b_muType->GetEntry(entry);
b_muIDbit->GetEntry(entry);
b_muD0->GetEntry(entry);
b_muDz->GetEntry(entry);
b_muSIP->GetEntry(entry);
b_muChi2NDF->GetEntry(entry);
b_muInnerD0->GetEntry(entry);
b_muInnerDz->GetEntry(entry);
b_muTrkLayers->GetEntry(entry);
b_muPixelLayers->GetEntry(entry);
b_muPixelHits->GetEntry(entry);
b_muMuonHits->GetEntry(entry);
b_muStations->GetEntry(entry);
b_muMatches->GetEntry(entry);
b_muTrkQuality->GetEntry(entry);
b_muIsoTrk->GetEntry(entry);
b_muPFChIso->GetEntry(entry);
b_muPFPhoIso->GetEntry(entry);
b_muPFNeuIso->GetEntry(entry);
b_muPFPUIso->GetEntry(entry);
b_muPFMiniIso->GetEntry(entry);
b_muFiredTrgs->GetEntry(entry);
b_muFiredL1Trgs->GetEntry(entry);
b_muInnervalidFraction->GetEntry(entry);
b_musegmentCompatibility->GetEntry(entry);
b_muchi2LocalPosition->GetEntry(entry);
b_mutrkKink->GetEntry(entry);
b_muBestTrkPtError->GetEntry(entry);
b_muBestTrkPt->GetEntry(entry);
b_nJet->GetEntry(entry);
b_jetPt->GetEntry(entry);
b_jetEn->GetEntry(entry);
b_jetEta->GetEntry(entry);
b_jetPhi->GetEntry(entry);
b_jetRawPt->GetEntry(entry);
b_jetRawEn->GetEntry(entry);
b_jetMt->GetEntry(entry);
b_jetArea->GetEntry(entry);
b_jetLeadTrackPt->GetEntry(entry);
b_jetLeadTrackEta->GetEntry(entry);
b_jetLeadTrackPhi->GetEntry(entry);
b_jetLepTrackPID->GetEntry(entry);
b_jetLepTrackPt->GetEntry(entry);
b_jetLepTrackEta->GetEntry(entry);
b_jetLepTrackPhi->GetEntry(entry);
b_jetCSV2BJetTags->GetEntry(entry);
b_jetJetProbabilityBJetTags->GetEntry(entry);
b_jetpfCombinedMVAV2BJetTags->GetEntry(entry);
b_jetPFLooseId->GetEntry(entry);
b_jetID->GetEntry(entry);
b_jetPUID->GetEntry(entry);
b_jetPUFullID->GetEntry(entry);
b_jetJECUnc->GetEntry(entry);
b_jetFiredTrgs->GetEntry(entry);
b_jetCHF->GetEntry(entry);
b_jetNHF->GetEntry(entry);
b_jetCEF->GetEntry(entry);
b_jetNEF->GetEntry(entry);
b_jetNCH->GetEntry(entry);
b_jetNNP->GetEntry(entry);
b_jetMUF->GetEntry(entry);
b_jetVtxPt->GetEntry(entry);
b_jetVtxMass->GetEntry(entry);
b_jetVtxNtrks->GetEntry(entry);
b_jetVtx3DVal->GetEntry(entry);
b_jetVtx3DSig->GetEntry(entry);
b_nAK8Jet->GetEntry(entry);
b_AK8JetPt->GetEntry(entry);
b_AK8JetEn->GetEntry(entry);
b_AK8JetRawPt->GetEntry(entry);
b_AK8JetRawEn->GetEntry(entry);
b_AK8JetEta->GetEntry(entry);
b_AK8JetPhi->GetEntry(entry);
b_AK8JetMass->GetEntry(entry);
b_AK8Jet_tau1->GetEntry(entry);
b_AK8Jet_tau2->GetEntry(entry);
b_AK8Jet_tau3->GetEntry(entry);
b_AK8JetCHF->GetEntry(entry);
b_AK8JetNHF->GetEntry(entry);
b_AK8JetCEF->GetEntry(entry);
b_AK8JetNEF->GetEntry(entry);
b_AK8JetNCH->GetEntry(entry);
b_AK8JetNNP->GetEntry(entry);
b_AK8JetMUF->GetEntry(entry);
b_AK8Jetnconstituents->GetEntry(entry);
b_AK8JetPFLooseId->GetEntry(entry);
b_AK8JetPFTightLepVetoId->GetEntry(entry);
b_AK8JetSoftDropMass->GetEntry(entry);
b_AK8JetSoftDropMassCorr->GetEntry(entry);
b_AK8JetPrunedMass->GetEntry(entry);
b_AK8JetPrunedMassCorr->GetEntry(entry);
b_AK8JetpfBoostedDSVBTag->GetEntry(entry);
b_AK8JetDSVnewV4->GetEntry(entry);
b_AK8JetCSV->GetEntry(entry);
b_AK8JetJECUnc->GetEntry(entry);
b_AK8JetL2L3corr->GetEntry(entry);
b_AK8puppiPt->GetEntry(entry);
b_AK8puppiMass->GetEntry(entry);
b_AK8puppiEta->GetEntry(entry);
b_AK8puppiPhi->GetEntry(entry);
b_AK8puppiTau1->GetEntry(entry);
b_AK8puppiTau2->GetEntry(entry);
b_AK8puppiTau3->GetEntry(entry);
b_AK8puppiSDL2L3corr->GetEntry(entry);
b_AK8puppiSDMass->GetEntry(entry);
b_AK8puppiSDMassL2L3Corr->GetEntry(entry);
b_nAK8SDSJ->GetEntry(entry);
b_AK8SDSJPt->GetEntry(entry);
b_AK8SDSJEta->GetEntry(entry);
b_AK8SDSJPhi->GetEntry(entry);
b_AK8SDSJMass->GetEntry(entry);
b_AK8SDSJE->GetEntry(entry);
b_AK8SDSJCharge->GetEntry(entry);
b_AK8SDSJFlavour->GetEntry(entry);
b_AK8SDSJCSV->GetEntry(entry);
b_nAK8puppiSDSJ->GetEntry(entry);
b_AK8puppiSDSJPt->GetEntry(entry);
b_AK8puppiSDSJEta->GetEntry(entry);
b_AK8puppiSDSJPhi->GetEntry(entry);
b_AK8puppiSDSJMass->GetEntry(entry);
b_AK8puppiSDSJE->GetEntry(entry);
b_AK8puppiSDSJCharge->GetEntry(entry);
b_AK8puppiSDSJFlavour->GetEntry(entry);
b_AK8puppiSDSJCSV->GetEntry(entry);

}
