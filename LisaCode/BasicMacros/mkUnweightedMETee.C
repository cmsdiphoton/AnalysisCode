
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

void mkUnweightedMETee(){

//gROOT->SetStyle("Plain");
//gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);


    int H_ref = 476;
    int W_ref = 720;
    
    // references for T, B, L, R
    float T = 0.08*H_ref;
    float B = 0.12*H_ref;
    float L = 0.12*W_ref;
    float R = 0.04*W_ref;
    
    TString cmsText     = "CMS";
    float cmsTextFont   = 61;  // default is helvetic-bold
    float cmsTextSize   = 0.5;
    
    TString extraText   = "Preliminary";
    float extraTextFont = 52;  // default is helvetica-italics
    float extraOverCmsTextSize  = 0.7;
    float extraTextSize = extraOverCmsTextSize*cmsTextSize;
    
    TString lumiText = "35.8 fb^{-1}";


Int_t   NbinsMET      = 24;
Float_t XbinsMET[25]  = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,75,85,100,115, 130, 150, 185,250 ,350};
 

TFile *f_data = TFile::Open("../analysis_NewMuonRecoVeto.root");


h_gg_MEt = new TH1F("h_gg_MEt","Missing Transverse Energy ",NbinsMET,XbinsMET);

TTree *gammatree = (TTree*)f_data->Get("ggtree");

gammatree->Draw("MET>>h_gg_MEt","MET<100");


TH1F* h_ee_MEt = (TH1F*)f_data->Get("h_ee_MET");
TH1F* h_ff_MEt = (TH1F*)f_data->Get("h_ff_MET");

h_ee_MEt->Sumw2();
h_gg_MEt->Sumw2();
h_ff_MEt->Sumw2();

h_ee_MEt= getOverflow(h_ee_MEt);
h_gg_MEt= getOverflow(h_gg_MEt);
h_ff_MEt= getOverflow(h_ff_MEt);

h_ee_MEt->Scale(1.0,"width");
h_gg_MEt->Scale(1.0,"width");
h_ff_MEt->Scale(1.0,"width");

h_ee_MEt->SetLineColor(kBlack);
h_gg_MEt->SetLineColor(kRed);
h_ff_MEt->SetLineColor(kBlue);

h_ee_MEt->SetMarkerStyle(20);
h_gg_MEt->SetMarkerStyle(20);
h_ff_MEt->SetMarkerStyle(20);

h_ee_MEt->SetMarkerColor(kBlack);
h_gg_MEt->SetMarkerColor(kRed);
h_ff_MEt->SetMarkerColor(kBlue);

/// compute scale factor for ee->scale for ggMET<50 GeV //

float denome = h_ee_MEt->Integral(0, 14);
float nume = h_gg_MEt->Integral(0, 14);
float scalefactorE;
		
cout << "sum ee "<< denome<< endl;		
cout << "sum gg" << nume << endl;
			if(denome !=0){
				scalefactorE = nume/denome;
			
			}
 		
		else {
		
			scalefactorE =1;
		}
		
		h_ee_MEt->Scale(scalefactorE);
		
		

cout << "scale factor for ee " << scalefactorE <<endl;


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
line->SetLineColor(kRed);
line->SetLineWidth(1);
//line->SetLineStyle(2);


///ready to make the plots!!!

TCanvas *c = new TCanvas("c","canvas",50,50,W_ref,H_ref);


c->cd();
c->SetFillColor(0);
c->SetBorderMode(0);
c->SetFrameFillStyle(0);
c->SetFrameBorderMode(0);
c->SetLeftMargin( L/W_ref );
c->SetRightMargin( R/W_ref );
c->SetTopMargin( T/H_ref );
c->SetBottomMargin( B/H_ref );
c->SetTickx(0);
c->SetTicky(0);
c->SetGrid();
gStyle->SetGridStyle(3);


c->cd();

TPad *pad1 = new TPad("pad1","pad1",0, 0.3, 1, 1.0);

 pad1->SetBottomMargin(0); // Upper and lower plot are joined
 pad1->SetLogy();
 pad1->SetGridx();         // Vertical grid
 pad1->SetGridy();         // Vertical grid
 pad1->Draw();             // Draw the upper pad: pad1
 pad1->cd();               // pad1 becomes the current pad
     // Draw h1
 
 
 h_ee_MEt->Draw();
 h_gg_MEt->Draw("sames");
 h_ff_MEt->Draw("sames");
 h_gg_MEt->Draw("sames");
 
 h_ee_MEt->SetTitle("Missing Transverse Energy");
 h_ee_MEt->GetXaxis()->SetTitle("Events");
 

TLegend *leg1 = new TLegend(0.59, 0.80, 0.90, 0.90);
  leg1->SetFillColor(kWhite);
  leg1->SetTextSize(0.03);
  //leg1->SetHeader("E_{T}^{miss} distributions");
  leg1->AddEntry(h_gg_MEt, "    #gamma#gamma", "LP");
  leg1->AddEntry(h_ee_MEt, "    ee", "LP");
  leg1->AddEntry(h_ff_MEt, "    ff", "LP");
  leg1->Draw();
 

c->cd();          // Go back to the main canvas before defining pad2
 TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
 pad2->SetTopMargin(0);
 pad2->SetBottomMargin(0.2);
 pad2->SetGridx(); // vertical grid
 pad2->Draw();
 pad2->cd();       // pad2 becomes the current pad
 



  h_gg_MEt->GetYaxis()->SetTitleSize(20);
 h_gg_MEt->GetYaxis()->SetTitleFont(43);
 //hCand10b->GetYaxis()->SetTitleOffset(1.55);
 h_gg_MEt ->GetYaxis()->SetTitleOffset(1.3);
 h_gg_MEt ->GetYaxis()->SetTitle("events");
 
 
 
 
 
 
 // Define the ratio plot
 
 TH1F *h3 = (TH1F*)h_gg_MEt->Clone("h3");
 h3->SetMinimum(0.5);  // Define Y ..
 h3->SetMaximum(2.4); // .. range
 h3->Divide(h_ee_MEt);




//line->Draw();

h3->Sumw2();
 h3->SetLineColor(kBlack);
 h3->SetMarkerColor(kBlack);
 h3->SetMarkerStyle(2);
 h3->Draw("epsame");
 line->Draw();
 
 /////////////
 
 
 // Define the ratio plot
 
 TH1F *h4 = (TH1F*)h_gg_MEt->Clone("h4");
 h4->SetMinimum(0.5);  // Define Y ..
 h4->SetMaximum(2.4); // .. range
 h4->Divide(h_ff_MEt);


 h4->Sumw2();
 h4->SetLineColor(kBlack);
 h4->SetMarkerColor(kBlack);
 h4->SetMarkerStyle(2);
 h4->Draw("epsame"); 
 
 /////////////
 h4->SetTitle("Ratios");

 h4->GetYaxis()->SetTitle("gg/ee");
 h4->GetYaxis()->SetNdivisions(505);
 h4->GetYaxis()->SetTitleSize(15);
 h4->GetYaxis()->SetTitleFont(43);
 h4->GetYaxis()->SetTitleOffset(1.0);
 h4->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
 h4->GetYaxis()->SetLabelSize(15);

 // X axis ratio plot settings
 h3->GetXaxis()->SetTitleSize(13);
 h3->GetXaxis()->SetTitleFont(43);
 h3->GetXaxis()->SetTitleOffset(3.0);
 h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
 h3->GetXaxis()->SetLabelSize(15);
 h3->GetXaxis()->SetTitle("E_{T}^{miss} (GeV)");
 


}
