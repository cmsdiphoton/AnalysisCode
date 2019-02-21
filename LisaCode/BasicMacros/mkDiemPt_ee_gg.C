
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

void mkDiemPt_ee_gg(){

//gROOT->SetStyle("Plain");
//gStyle->SetOptTitle(0);
//gStyle->SetOptStat(0);


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
 


 //new DiEmPt bins 
 
 Int_t Nbins = 26;
 float Xbins[27] = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,70,85,100,115, 130, 150,175,200,300,400,600};
 

TFile *f_data = TFile::Open("../analysis_NewMuonRecoVeto.root");


TH1F *h_gg_DiEmPt = new TH1F("h_gg_DiEmPt","DiEM p_{T} Distrubution ",Nbins,Xbins);

TTree *gammatree = (TTree*)f_data->Get("ggtree");

gammatree->Draw("Dipho_Pt>>h_gg_DiEmPt");

TH1F* h_ee_DiEmPt = new TH1F("h_ee_DiEmPt","h_ee_DiEmPt",Nbins,Xbins);

TTree *eletree = (TTree*)f_data->Get("eetree");

eletree->Draw("Dipho_Pt>>h_ee_DiEmPt");


h_ee_DiEmPt->Sumw2();
h_gg_DiEmPt->Sumw2();


h_ee_DiEmPt= getOverflow(h_ee_DiEmPt);
h_gg_DiEmPt= getOverflow(h_gg_DiEmPt);


h_ee_DiEmPt->Scale(1.0,"width");
h_gg_DiEmPt->Scale(1.0,"width");


h_ee_DiEmPt->SetLineColor(kBlack);
h_gg_DiEmPt->SetLineColor(kRed);

h_ee_DiEmPt->SetMarkerStyle(20);
h_gg_DiEmPt->SetMarkerStyle(20);

h_ee_DiEmPt->SetMarkerColor(kBlue);
h_gg_DiEmPt->SetMarkerColor(kRed);

/// compute scale factor for ff ->scale for ggMET<50 GeV //

float denomf = h_ee_DiEmPt->Integral(0, 14);
float numf = h_gg_DiEmPt->Integral(0, 14);
float scalefactorF;
		
cout << "sum ff "<< denomf<< endl;		
cout << "sum gg" << numf << endl;
			if(denomf !=0){
				scalefactorF = numf/denomf;
			
			}
 		
		else {
		
			scalefactorF =1;
		}
		
		h_ee_DiEmPt->Scale(scalefactorF);
		
		

cout << "scale factor for ff " << scalefactorF <<endl;



TLine *line = new TLine(0.,1.,600.,1.);
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
 
 
 h_ee_DiEmPt->Draw("histo");
 h_gg_DiEmPt->Draw("histosames");

 
 h_ee_DiEmPt->SetTitle("DiEM-P_{T} Distributions");
 h_ee_DiEmPt->GetXaxis()->SetTitle("Events");
 

TLegend *leg1 = new TLegend(0.59, 0.80, 0.90, 0.90);
  leg1->SetFillColor(kWhite);
  leg1->SetTextSize(0.03);
  //leg1->SetHeader("E_{T}^{miss} distributions");
  leg1->AddEntry(h_gg_DiEmPt, "    #gamma#gamma", "L");
  leg1->AddEntry(h_ee_DiEmPt, "    ee", "L");
  
  leg1->Draw();
 

c->cd();          // Go back to the main canvas before defining pad2
 TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
 pad2->SetTopMargin(0);
 pad2->SetBottomMargin(0.2);
 pad2->SetGridx(); // vertical grid
 pad2->Draw();
 pad2->cd();       // pad2 becomes the current pad
 

gStyle->SetOptStat(0);


  h_gg_DiEmPt->GetYaxis()->SetTitleSize(20);
 h_gg_DiEmPt->GetYaxis()->SetTitleFont(43);
 //hCand10b->GetYaxis()->SetTitleOffset(1.55);
 h_gg_DiEmPt ->GetYaxis()->SetTitleOffset(1.3);
 h_gg_DiEmPt ->GetYaxis()->SetTitle("events");
 
 
 
 
 
 
 // Define the ratio plot
 
 TH1F *h3 = (TH1F*)h_gg_DiEmPt->Clone("h3");
 h3->SetMinimum(0.0);  // Define Y ..
 h3->SetMaximum(4.5); // .. range
 h3->Divide(h_ee_DiEmPt);




//line->Draw();

h3->Sumw2();
 h3->SetLineColor(kBlack);
 h3->SetMarkerColor(kBlack);
 h3->SetMarkerStyle(2);
 h3->Draw("epsame");
 line->Draw();
 
 h3->SetTitle("Ratios");

 h3->GetYaxis()->SetTitle("gg/ee");
 h3->GetYaxis()->SetNdivisions(505);
 h3->GetYaxis()->SetTitleSize(15);
 h3->GetYaxis()->SetTitleFont(43);
 h3->GetYaxis()->SetTitleOffset(1.0);
 h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
 h3->GetYaxis()->SetLabelSize(15);

 // X axis ratio plot settings
 h3->GetXaxis()->SetTitleSize(13);
 h3->GetXaxis()->SetTitleFont(43);
 h3->GetXaxis()->SetTitleOffset(3.5);
 h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
 h3->GetXaxis()->SetLabelSize(15);
 h3->GetXaxis()->SetTitle("DiEM-p_{T} (GeV)");
 


}
