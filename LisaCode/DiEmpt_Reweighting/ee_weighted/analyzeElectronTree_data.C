#define analyzeElectronTree_data_cxx
#include "analyzeElectronTree_data.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include "TRandom3.h"
#include <TTree.h>
#include <algorithm>
#include <utility>
/////////////////////////////////////////


Int_t   NbinsMET      = 21;
Float_t XbinsMET[22]  = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,75,85,100,115, 130, 150};


Int_t Nbins = 24;
float Xbins[25] = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,70,85,100,115, 130, 150,175,200,300};


TFile *f_ratio   = new TFile("../EE_DiEmPtRatio.root");
TH1F  *eeweights = (TH1F*)f_ratio->Get("h_ee_DiEmPt_ratio");

//Counters
int entry_ = 0;
/////////////////////////////////////////
TTree *ANtree_;

//--------------------
TH1F* hh_St;
TH1F* hh_diphoPt;
TH1F* hh_MEt;
TH1F* hh_w0MEt;
TH1F* hh_wMEt;

TH1F *h_Electron_Dipho_Pt;
TH1F *h_Electron_Dipho_Pt_w;
TH1F *h_Electron_Dipho_Pt_ew;

TH1F *h_Electron_St;
TH1F *h_Electron_St_w;
TH1F *h_Electron_St_ew;

TH1F *h_Electron_MEt;
TH1F *h_Electron_MEt_w;
TH1F *h_Electron_MEt_ew;




//--------------------
int   run_;
int   event_;
int   totEvents_;
//--------------------

float diphoPt_;
float St_;
float MEt_;
float St_weighted;
float MEt_weighted;





TFile f("./analysisElectronWithWeightsAllMetBin.root","RECREATE");



void analyzeElectronTree_data::Begin(TTree * /*tree*/)
{


 ANtree_ = new TTree("Analysis_tree", "Event data");
 ///////////////////////////////////
 hh_St          = new TH1F("hh_St",          "St plot", 100, 0., 800);
 hh_diphoPt     = new TH1F("hh_diphoPt",     "diphoPt plot", 100, 0., 800);
 hh_MEt         = new TH1F("hh_MEt",         "MEt plot", 100, 0., 300);
 hh_wMEt        = new TH1F("hh_wMEt",        "weighted MEt plot", 100, 0., 300);

 h_Electron_Dipho_Pt	  = new TH1F("h_Electron_Dipho_Pt",   "h_Electron_Dipho_Pt",Nbins,Xbins);
 h_Electron_Dipho_Pt_w	  = new TH1F("h_Electron_Dipho_Pt_w", "h_Electron_Dipho_Pt_w",Nbins,Xbins);
 h_Electron_Dipho_Pt_ew 	= new TH1F("h_Electron_Dipho_Pt_ew","h_Electron_Dipho_Pt_ew",Nbins,Xbins);

 h_Electron_St 		= new TH1F("h_Electron_St","h_Electron_St",30,80,600);
 h_Electron_St_w 		= new TH1F("h_Electron_St_w","h_Electron_St_w",30,80,600);
 h_Electron_St_ew 		= new TH1F("h_Electron_St_ew","h_Electron_St_ew",30,80,600);

 h_Electron_MEt 		= new TH1F("h_Electron_MEt","h_Electron_MEt",NbinsMET,XbinsMET);
 h_Electron_MEt_w 		= new TH1F("h_Electron_MEt_w","h_Electron_MEt_w",NbinsMET,XbinsMET);
 h_Electron_MEt_ew 		= new TH1F("h_Electron_MEt_ew","h_Electron_MEt_ew",NbinsMET,XbinsMET);





 ///////////////////////////////////
 ANtree_->Branch("run",  	&run_,		"run/I");
 ANtree_->Branch("event",	&event_,	"event/I");
 ANtree_->Branch("totEvents",	&totEvents_,	"totEvents/I");
 ///////////////////////////////////
 ANtree_->Branch("diphoPt",     &diphoPt_,	"diphoPt/F");
 ANtree_->Branch("St",		&St_,		"St/F");
 ANtree_->Branch("MEt",		&MEt_,		"MEt/F");


}

////////////////////////////////////////////////////////
void analyzeElectronTree_data::SlaveBegin(TTree * /*tree*/)
{}

//////////////////////////////////////////////////////

Bool_t analyzeElectronTree_data::Process(Long64_t entry)
{

 GetAllEntries(entry);
// if(entry_%1000 == 0) cout<<"#Run = "<<Run<<", #event = "<<Event<<", #entry = "<<entry_<<endl;
 entry_++;


diphoPt_	= Dipho_Pt;
 St_   		= ST;
 MEt_  		= MET;



 hh_St->Fill(St_);
 hh_diphoPt->Fill(diphoPt_);
 hh_MEt->Fill(MEt_);
 h_Electron_Dipho_Pt->Fill(diphoPt_);
 h_Electron_MEt->Fill(MEt_);
 h_Electron_St->Fill(St_);


 float diempt_weight=1;

 int bin = eeweights->FindFixBin(float(Dipho_Pt));
 diempt_weight = eeweights->GetBinContent(bin);

 h_Electron_MEt_w->Fill(MEt_,diempt_weight);

   ANtree_->Fill();

   return kTRUE;
}

void analyzeElectronTree_data::SlaveTerminate() {}

void analyzeElectronTree_data::Terminate()
{

float lastbinMET = h_Electron_MEt_w->Integral(NbinsMET,NbinsMET+1);

h_Electron_MEt_w->SetBinContent(Nbins,lastbinMET);


cout<<"\nTerminating...\n--Summary: "<<entry_<<" total analyzed events.\n"<<endl;
 f.cd();
 f.Write();
 f.Close();


}


void analyzeElectronTree_data::GetAllEntries(Long64_t entry)
{

//List of branches

b_Run->GetEntry(entry);
b_Event->GetEntry(entry);
b_Lumis->GetEntry(entry);
b_Rho->GetEntry(entry);
b_PUtrue->GetEntry(entry);
b_PUweight->GetEntry(entry);
b_MET_Filters->GetEntry(entry);
b_MET->GetEntry(entry);
b_MET_Phi->GetEntry(entry);
b_MET_sumEt->GetEntry(entry);
b_MET_mEtSig->GetEntry(entry);
b_MET_Sig->GetEntry(entry);
b_NJet->GetEntry(entry);
b_MT2->GetEntry(entry);
b_NEle->GetEntry(entry);
b_Dipho_Pt->GetEntry(entry);
b_Dipho_Mass->GetEntry(entry);
b_Ele_Pt->GetEntry(entry);
b_Ele_Eta->GetEntry(entry);
b_Ele_Phi->GetEntry(entry);
b_Ele_SigmaIEtaIEta->GetEntry(entry);
b_Ele_SigmaIEtaIPhi->GetEntry(entry);
b_Ele_SigmaIPhiIPhi->GetEntry(entry);
b_HT->GetEntry(entry);
b_ST->GetEntry(entry);
b_MHT->GetEntry(entry);
b_MHTPhi->GetEntry(entry);





}
