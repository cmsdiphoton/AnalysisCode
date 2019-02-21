#include "CMS_lumi.C"
#include <iostream>
#include <TH1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TAttFill.h>
#include <TFile.h>
#include <iostream>
#include <TFractionFitter.h>
#include <TCanvas.h>
#include <TPad.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>



TH1F *getOverflow(TH1F *h_Sample){
  int bin = h_Sample->GetNbinsX();
  float lastBinValue = h_Sample->GetBinContent(bin);
  float lastBinError = h_Sample->GetBinError(bin);

  float lastBinOverflowValue = h_Sample->GetBinContent(bin+1);
  float lastBinOverflowError = h_Sample->GetBinError(bin+1);

  float finalValue = lastBinValue + lastBinOverflowValue;
  float finalError = sqrt(lastBinError*lastBinError + lastBinOverflowError*lastBinOverflowError);

  h_Sample->SetBinContent(bin, finalValue);
  h_Sample->SetBinContent(bin+1, 0);
  h_Sample->SetBinError(bin, finalError);
  h_Sample->SetBinError(bin+1, 0);

  return h_Sample;

}

void mkDiemPtReweight(){

gROOT->SetStyle("Plain");
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);



Int_t   NbinsMET      = 24;
Float_t XbinsMET[25]  = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,75,85,100,115, 130, 150, 185,250 ,350};


TFile *f_data = TFile::Open("../newSkimm/analysis_NewMuonRecoVetoRemoveFinal_test.root");


TH1F *h_gg_MEt = new TH1F("h_gg_MEt","Missing Transverse Energy ",NbinsMET,XbinsMET);

TTree *gammatree = (TTree*)f_data->Get("ggtree");

gammatree->Draw("MET>>h_gg_MEt","MET<100");
//gammatree->Draw("MET>>h_gg_MEt");


TFile *f_ff_wgh = TFile::Open("../newSkimm/DiEMpt_Reweighting/ff_weighted/analysisFakeWithWeightsNJets2AllMetBin.root");

TH1F* h_ff_MEt = (TH1F*)f_ff_wgh->Get("h_Fake_MEt_w");


h_ff_MEt->Sumw2();
h_gg_MEt->Sumw2();


h_ff_MEt= getOverflow(h_ff_MEt);
h_gg_MEt= getOverflow(h_gg_MEt);


h_ff_MEt->Scale(1.0,"width");
h_gg_MEt->Scale(1.0,"width");


h_ff_MEt->SetLineColor(kBlue);
h_gg_MEt->SetLineColor(kRed);

h_ff_MEt->SetMarkerStyle(20);
h_gg_MEt->SetMarkerStyle(20);

h_ff_MEt->SetMarkerColor(kBlue);
h_gg_MEt->SetMarkerColor(kRed);

/// compute scale factor for ff ->scale for ggMET<50 GeV //

float denomf = h_ff_MEt->Integral(0, 14);
float numf = h_gg_MEt->Integral(0, 14);
float scalefactorF;

cout << "sum ff "<< denomf<< endl;
cout << "sum gg" << numf << endl;
			if(denomf !=0){
				scalefactorF = numf/denomf;

			}

		else {

			scalefactorF =1;
		}

		h_ff_MEt->Scale(scalefactorF);



cout << "scale factor for ff " << scalefactorF <<endl;



TLine *line = new TLine(0.,1.,350.,1.);
line->SetLineColor(kBlack);
line->SetLineWidth(2);
line->SetLineStyle(2);


///ready to make the plots!!!

//TCanvas *c = new TCanvas("c","canvas",50,50,W_ref,H_ref);

TCanvas *c = new TCanvas("c","canvas",950,950);

c->cd();

TPad *pad1 = new TPad("pad1","pad1",0, 0.3, 1, 1.0);

 pad1->SetBottomMargin(0); // Upper and lower plot are joined
 pad1->SetLogy();
 //pad1->SetGridx();         // Vertical grid
 //pad1->SetGridy();         // Vertical grid
 pad1->Draw();             // Draw the upper pad: pad1
 pad1->cd();               // pad1 becomes the current pad
     // Draw h1


 h_ff_MEt->Draw();
 h_gg_MEt->Draw("sames");

 h_ff_MEt->GetYaxis()->SetTitle("Events/GeV");

 //h_ff_MEt->SetTitle("Missing Transverse Energy");
 //h_ff_MEt->GetXaxis()->SetTitle("Events");
 h_ff_MEt->GetYaxis()->SetTitleSize(0.045);
 h_ff_MEt ->GetYaxis()->SetTitleOffset(0.8);

TLegend *leg1 = new TLegend(0.29, 0.72, 0.88, 0.88);
  leg1->SetFillColor(kWhite);
  leg1->SetBorderSize(0);
  leg1->SetShadowColor(kWhite);

  leg1->SetTextSize(0.04);
  //leg1->SetHeader("E_{T}^{miss} distributions");
  leg1->AddEntry(h_gg_MEt, " Candidate #gamma#gamma  Sample" , "LP");
  leg1->AddEntry(h_ff_MEt, " di-EM p_{T} Reweighted Background Prediction", "LP");

  leg1->Draw();

 CMS_lumi(pad1,4,1);

c->cd();          // Go back to the main canvas before defining pad2
 //TPad *pad2 = new TPad("pad2", "pad2", 0, 0.02, 1, 0.31);
 TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.31);
 pad2->SetTopMargin(0);
 pad2->SetBottomMargin(0.2);
 //pad2->SetGridx(); // vertical grid
 pad2->Draw();
 pad2->cd();       // pad2 becomes the current pad


  h_gg_MEt->GetYaxis()->SetTitleSize(28);
   h_gg_MEt->GetYaxis()->SetTitleFont(43);
 //hCand10b->GetYaxis()->SetTitleOffset(1.55);
 h_gg_MEt ->GetYaxis()->SetTitleOffset(1.0);
 h_gg_MEt ->GetYaxis()->SetTitle("events");

 h_gg_MEt->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
 h_gg_MEt->GetYaxis()->SetLabelSize(24);



 // Define the ratio plot

 TH1F *h3 = (TH1F*)h_gg_MEt->Clone("h3");
 h3->SetMinimum(0.5);  // Define Y ..
 h3->SetMaximum(1.5); // .. range
 //h3->SetMaximum(2.8); // .. range
 h3->Divide(h_ff_MEt);

 h3->SetMarkerStyle(kFullCircle);
 h3->SetMarkerSize(1.2);



//line->Draw();

h3->Sumw2();
 h3->SetLineColor(kBlack);
 h3->SetMarkerColor(kBlack);

 h3->Draw("epsame");
 line->Draw();

 h3->SetTitle("Ratios");

 h3->GetYaxis()->SetTitle("#gamma#gamma/ff");
 h3->GetYaxis()->SetNdivisions(505);
 h3->GetYaxis()->SetTitleSize(22);
 h3->GetYaxis()->SetTitleFont(43);
 h3->GetYaxis()->SetTitleOffset(0.8);
 h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
 h3->GetYaxis()->SetLabelSize(19);

 // X axis ratio plot settings
 h3->GetXaxis()->SetTitleSize(22);
 h3->GetXaxis()->SetTitleFont(43);
 h3->GetXaxis()->SetTitleOffset(3.2);
 h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
 h3->GetXaxis()->SetLabelSize(24);
 h3->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");



}
