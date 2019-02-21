#define analyzeFakeTree_data_cxx
#include "analyzeFakeTree_data.h"
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
 
//Int_t   Nbins        = 24;
//Float_t Xbins[25]    = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,75,85,100,115, 130, 150, 185,300 ,400};
 
Int_t Nbins = 24;
float Xbins[25] = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,70,85,100,115, 130, 150,175,200,300};
 

TFile *f_ratio   = new TFile("../FF_DiEMPt_ratio.root");
TH1F  *ff_weights = (TH1F*)f_ratio->Get("h_ff_DiEmPt_ratio"); 
 


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

TH1F *h_Fake_Dipho_Pt;
TH1F *h_Fake_Dipho_Pt_w;
TH1F *h_Fake_Dipho_Pt_ew;

TH1F *h_Fake_St; 	
TH1F *h_Fake_St_w; 	
TH1F *h_Fake_St_ew; 	

TH1F *h_Fake_MEt; 	
TH1F *h_Fake_MEt_w; 	
TH1F *h_Fake_MEt_ew; 	




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


TFile f("./analysisFakeWithWeightsNJets2AllMetBin.root","RECREATE");



void analyzeFakeTree_data::Begin(TTree * /*tree*/)
{


 ANtree_ = new TTree("Analysis_tree", "Event data");
 ///////////////////////////////////
 hh_St          = new TH1F("hh_St",          "St plot", 100, 0., 800);
 hh_diphoPt     = new TH1F("hh_diphoPt",     "diphoPt plot", 100, 0., 800);
 hh_MEt         = new TH1F("hh_MEt",         "MEt plot", 100, 0., 300);
 hh_wMEt        = new TH1F("hh_wMEt",        "weighted MEt plot", 100, 0., 300);
 
 h_Fake_Dipho_Pt	  = new TH1F("h_Fake_Dipho_Pt",   "h_Fake_Dipho_Pt",Nbins,Xbins);
 h_Fake_Dipho_Pt_w	  = new TH1F("h_Fake_Dipho_Pt_w", "h_Fake_Dipho_Pt_w",Nbins,Xbins);
 h_Fake_Dipho_Pt_ew 	= new TH1F("h_Fake_Dipho_Pt_ew","h_Fake_Dipho_Pt_ew",Nbins,Xbins);
 
 h_Fake_St 		= new TH1F("h_Fake_St","h_Fake_St",30,80,600);
 h_Fake_St_w 		= new TH1F("h_Fake_St_w","h_Fake_St_w",30,80,600);
 h_Fake_St_ew 		= new TH1F("h_Fake_St_ew","h_Fake_St_ew",30,80,600);
 
 h_Fake_MEt 		= new TH1F("h_Fake_MEt","h_Fake_MEt",NbinsMET,XbinsMET);
 h_Fake_MEt_w 		= new TH1F("h_Fake_MEt_w","h_Fake_MEt_w",NbinsMET,XbinsMET);
 h_Fake_MEt_ew 		= new TH1F("h_Fake_MEt_ew","h_Fake_MEt_ew",NbinsMET,XbinsMET);
 
 
 
 
 
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
void analyzeFakeTree_data::SlaveBegin(TTree * /*tree*/)
{}

//////////////////////////////////////////////////////

Bool_t analyzeFakeTree_data::Process(Long64_t entry)
{
 
 GetAllEntries(entry);
 //if(entry_%1000 == 0) cout<<"#Run = "<<Run<<", #event = "<<Event<<", #entry = "<<entry_<<endl;
 entry_++; 
 
 
diphoPt_	= Dipho_Pt;
 St_   		= ST;
 MEt_  		= MET;
 
 
 
 hh_St->Fill(St_);	  
 hh_diphoPt->Fill(diphoPt_);
 hh_MEt->Fill(MEt_);
 h_Fake_Dipho_Pt->Fill(diphoPt_);
 h_Fake_MEt->Fill(MEt_);
 h_Fake_St->Fill(St_);	
 
   
   
 float diempt_weight=1;
 
 int bin = ff_weights->FindFixBin(float(Dipho_Pt));  
 diempt_weight = ff_weights->GetBinContent(bin);  
 
 h_Fake_MEt_w->Fill(MEt_,diempt_weight);
 

if( MEt_>85 && MEt_<100){

cout <<"Event number " << Event << ", " << "diEMpT = " << Dipho_Pt <<" and weight " << diempt_weight <<"MET = "<< MEt_ <<endl;

}


   ANtree_->Fill();
   

   return kTRUE;
}

void analyzeFakeTree_data::SlaveTerminate() {}

void analyzeFakeTree_data::Terminate()
{

float lastbinMET = h_Fake_MEt_w->Integral(NbinsMET,NbinsMET+1);

h_Fake_MEt_w->SetBinContent(Nbins,lastbinMET);


cout<<"\nTerminating...\n--Summary: "<<entry_<<" total analyzed events.\n"<<endl;
 f.cd();
 f.Write();
 f.Close();

   
}


void analyzeFakeTree_data::GetAllEntries(Long64_t entry)
{

//List of branches


b_Run->GetEntry(entry);
b_Event->GetEntry(entry);
b_Lumis->GetEntry(entry);
b_Rho->GetEntry(entry);
b_PUtrue->GetEntry(entry);
b_PUweight->GetEntry(entry);
b_NJet->GetEntry(entry);
b_MET_Filters->GetEntry(entry);
b_MET->GetEntry(entry);
b_MET_Phi->GetEntry(entry);
b_MET_sumEt->GetEntry(entry);
b_MET_mEtSig->GetEntry(entry);
b_MET_Sig->GetEntry(entry);
b_MT2->GetEntry(entry);
b_HT->GetEntry(entry);
b_ST->GetEntry(entry);
b_MHT->GetEntry(entry);
b_MHTPhi->GetEntry(entry);
b_NFake->GetEntry(entry);
b_Dipho_Pt->GetEntry(entry);
b_Dipho_Mass->GetEntry(entry);
b_Pho_E->GetEntry(entry);
b_Pho_Et->GetEntry(entry);
b_Pho_Eta->GetEntry(entry);
b_Pho_Phi->GetEntry(entry);
b_Pho_hasPixelSeed->GetEntry(entry);
b_Pho_EleVeto->GetEntry(entry);
b_Pho_R9->GetEntry(entry);
b_Pho_HoverE->GetEntry(entry);
b_Pho_SigmaIEtaIEta->GetEntry(entry);
b_Pho_SigmaIEtaIPhi->GetEntry(entry);
b_Pho_SigmaIPhiIPhi->GetEntry(entry);
b_Pho_IDMVA->GetEntry(entry);
b_Pho_IDbit->GetEntry(entry);








}
