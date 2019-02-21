#include <math.h>

void mkDataFakeRatio(){
 gROOT->SetStyle("Plain");
 gStyle->SetOptTitle(0);
 gStyle->SetOptStat(0);


 //Int_t   Nbins     = 24;
 //Float_t Xbins[25] = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,75,85,100,115, 130, 150, 185,250,400};

 Int_t Nbins = 24;
 float Xbins[25] = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,70,85,100,115, 130, 150,175,200,300};


////////////////Background//////////////

 ////////////datatoInf/////////////

TFile *f1  = TFile::Open("../analysis_NewMuonRecoVetoNewMetBinning.root");
TTree *gg1 = (TTree*)gDirectory->Get("ggtree");
TTree *ee1 = (TTree*)gDirectory->Get("eetree");
TTree *ff1 = (TTree*)gDirectory->Get("fftree");


TH1F *h_G_data = new TH1F("h_G_data","h_G_data",Nbins,Xbins);
TH1F *h_F_data = new TH1F("h_F_data","h_F_data",Nbins,Xbins);
TH1F *h_E_data = new TH1F("h_E_data","h_E_data",Nbins,Xbins);
gg1->Draw("Dipho_Pt>>h_G_data");
ff1->Draw("Dipho_Pt>>h_F_data");
ee1->Draw("Dipho_Pt>>h_E_data");
h_G_data->Sumw2();
h_F_data->Sumw2();
h_E_data->Sumw2();


float lastbinG = h_G_data->Integral(Nbins,Nbins+1);
float lastbinF = h_F_data->Integral(Nbins,Nbins+1);
float lastbinE = h_E_data->Integral(Nbins,Nbins+1);

//h_G_data->SetBinContent(Nbins,lastbinG);
//h_F_data->SetBinContent(Nbins,lastbinF);
//h_E_data->SetBinContent(Nbins,lastbinE);

h_G_data->Scale(1./h_G_data->Integral());
h_F_data->Scale(1./h_F_data->Integral());
h_E_data->Scale(1./h_E_data->Integral());


TH1F *h_ff_DiEmPt_ratio = (TH1F*)h_G_data->Clone("h_ff_DiEmPt_ratio");
 h_ff_DiEmPt_ratio->SetLineColor(kBlack);
 h_ff_DiEmPt_ratio->SetMinimum(0.0);
 h_ff_DiEmPt_ratio->SetMaximum(4.5);
 h_ff_DiEmPt_ratio->SetStats(0);
 h_ff_DiEmPt_ratio->Divide(h_F_data);
 //h_ff_DiEmPt_ratio->Divide(h_E_BG);
 h_ff_DiEmPt_ratio->SetMarkerColor(1);
 h_ff_DiEmPt_ratio->SetMarkerSize(1);
 h_ff_DiEmPt_ratio->SetMarkerStyle(20);


 TFile outputfile ("FF_DiEMPt_ratio.root","RECREATE");

 h_ff_DiEmPt_ratio->Write();
 outputfile.Close();


cout<<"\nwithout signal for data electrons;"<<endl;
 for(int i=1; i<28; i++){
   cout<<h_ff_DiEmPt_ratio->GetBinContent(i)<<", ";
 }


  cout<<"\nRatio error for data elecrrons;"<<endl;
 for(int i=1; i<27; i++){
   cout<<(h_ff_DiEmPt_ratio->GetBinError(i))<<", ";
 }



 }
