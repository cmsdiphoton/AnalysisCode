
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


TH1F * scaleWidth(TH1F* h1){
    int numBins = h1->GetXaxis()->GetNbins();
    for(int i = 0; i <= numBins; i++){
        float errorUp  = h1->GetBinErrorUp(i) ;
        float errorLow = h1->GetBinErrorLow(i) ;
        float width = h1->GetBinWidth(i);
        float content = h1->GetBinContent(i);
        float averageErr = (errorUp + errorLow )/(2*width);
        h1->SetBinContent(i,content/width);
        h1->SetBinError(i,averageErr);
    }
    return h1;
}

TH1F * scaleHist(TH1F* h1, float scale){
    int numBins = h1->GetXaxis()->GetNbins();
    for(int i = 0; i <= numBins; i++){
        float errorUp  = h1->GetBinErrorUp(i) ;
        float errorLow = h1->GetBinErrorLow(i) ;
        float content = h1->GetBinContent(i);
        float averageErr = scale*(errorUp + errorLow )/(2);
        h1->SetBinContent(i,scale*content);
        h1->SetBinError(i,averageErr);
    }
    return h1;
}

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



void FinalTexMaker(TLatex *&tex1, TLatex *&tex2, TLatex *&tex3){
    //tex1 = new TLatex(0.88,0.909,"35.9 fb^{-1} (13 TeV)");
    tex1 = new TLatex(0.88,0.92,"35.9 fb^{-1} (13 TeV)");
    tex1->SetNDC();
    tex1->SetTextAlign(31);
    tex1->SetTextFont(42);
//    tex1->SetTextSize(0.0546);
    tex1->SetTextSize(0.07);
    tex1->SetLineWidth(2);
    tex1->Draw();
    
    tex2 = new TLatex(0.12,0.967,"CMS");//(0.15,0.88,"CMS");
    tex2->SetNDC();
    tex2->SetTextAlign(13);
    tex2->SetTextFont(61);
    tex2->SetTextSize(0.075);
    tex2->SetLineWidth(1);
    tex2->Draw();
    
//    tex3 = new TLatex(0.23,0.96,"Preliminary");//(0.15,0.84,"Simulation")
    tex3 = new TLatex(0.22,0.955,"Preliminary");//(0.15,0.84,"Simulation")
    tex3->SetNDC();
    tex3->SetTextAlign(13);
    tex3->SetTextFont(52);
    tex3->SetTextSize(0.072);
    tex3->SetLineWidth(1);
//    tex3->Draw();
    
}


void FinalPlot(){
    
    int xlow(0),xup(350);
    bool blinded = false;
    gStyle->SetLineWidth(3);
    
    int   NbinsMET      = 24;
    Float_t XbinsMET[25]  = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,75,85,100,115, 130, 150, 185,250 ,350};
    
    
    /// Open MC signal //
    //------------------------------------------
    TFile *f_mGlu2000nNeu1000 = new TFile("hist_T5gg_2000_1000.root","READ");
    TH1F  *h_mGlu2000nNeu1000 = (TH1F*)f_mGlu2000nNeu1000->Get("ggMet");

    TFile *f_mGlu1700nNeu1000= new TFile("hist_T5gg_1700_1000.root","READ");
    TH1F  *h_mGlu1700nNeu1000 = (TH1F*)f_mGlu1700nNeu1000->Get("ggMet");
    
    float numEvents = 11942;
    float xsec = 0.00470323*1000;
    float lumi   = 35.9218;
    float scale1700 = lumi*xsec/numEvents;
    h_mGlu1700nNeu1000 = scaleHist(h_mGlu1700nNeu1000,scale1700);
    h_mGlu1700nNeu1000 = scaleWidth(h_mGlu1700nNeu1000);
    
    h_mGlu1700nNeu1000->SetMarkerStyle(kFullTriangleDown);
    h_mGlu1700nNeu1000->SetLineColor(kMagenta+1);
    h_mGlu1700nNeu1000->SetMarkerColor(kMagenta+1);
    h_mGlu1700nNeu1000->SetMarkerSize(1.1);
    
    numEvents = 5273;
    xsec = 0.000981077*1000;
    float scale2000 = lumi*xsec/numEvents;
    h_mGlu2000nNeu1000 = scaleHist(h_mGlu2000nNeu1000,scale2000);
    h_mGlu2000nNeu1000 = scaleWidth(h_mGlu2000nNeu1000);
    
    h_mGlu2000nNeu1000->SetMarkerStyle(kFullDiamond);
    h_mGlu2000nNeu1000->SetLineColor(kViolet-6);
    h_mGlu2000nNeu1000->SetMarkerColor(kViolet-6);
    h_mGlu2000nNeu1000->SetMarkerSize(1.2);
    /*Float_t sigscale = 35.882548708 * 0.981077 / 5273;
    
    Double_t sigx[24], sigy[24], sigxeu[24], sigxel[24], sigyeu[24], sigyel[24];
    for(int i = 0; i < NbinsMET;i++){
        float binwidth = h_mGlu2000nNeu1000->GetBinWidth(i+1);
        sigx[i] = h_mGlu2000nNeu1000->GetBinCenter(i+1);
        sigy[i] = h_mGlu2000nNeu1000->GetBinContent(i+1)*sigscale/binwidth;
        sigxel[i] = binwidth/2;
        sigxeu[i] = binwidth/2;
        sigyel[i] = h_mGlu2000nNeu1000->GetBinErrorUp(i+1)*sigscale/binwidth;
        sigyeu[i] = h_mGlu2000nNeu1000->GetBinErrorLow(i+1)*sigscale/binwidth;
    }
    TGraphAsymmErrors * gr_mGlu2000nNeu1000 = new TGraphAsymmErrors(NbinsMET,sigx,sigy,sigxel,sigxeu,sigyel,sigyeu);
    gr_mGlu2000nNeu1000->SetName("gr_mGlu2000nNeu1000");
    
    gr_mGlu2000nNeu1000 ->SetMarkerStyle(kFullTriangleUp);
    gr_mGlu2000nNeu1000 ->SetLineColor(kCyan);
    gr_mGlu2000nNeu1000 ->SetMarkerColor(kCyan);
    
   
   
    //------------------------------------------
    sigscale = 35.882548708 * 0.981077 / 4956;

    TFile *f_mGlu2000nNeu100 = new TFile("T5Wg_2000_100_unMET.root","READ");
    TH1F  *h_mGlu2000nNeu100 = (TH1F*)f_mGlu2000nNeu100->Get("h_T5Wg_2000_100");
    
    
    h_mGlu2000nNeu100=getOverflow(h_mGlu2000nNeu100);
    
    for(int i = 0; i < NbinsMET;i++){
        float binwidth = h_mGlu2000nNeu100->GetBinWidth(i+1);
        sigx[i] = h_mGlu2000nNeu100->GetBinCenter(i+1);
        sigy[i] = h_mGlu2000nNeu100->GetBinContent(i+1)*sigscale/binwidth;
        sigxel[i] = binwidth/2;
        sigxeu[i] = binwidth/2;
        sigyel[i] = h_mGlu2000nNeu100->GetBinErrorUp(i+1)*sigscale/binwidth;
        sigyeu[i] = h_mGlu2000nNeu100->GetBinErrorLow(i+1)*sigscale/binwidth;
    }
    TGraphAsymmErrors * gr_mGlu2000nNeu100 = new TGraphAsymmErrors(NbinsMET,sigx,sigy,sigxel,sigxeu,sigyel,sigyeu);
    gr_mGlu2000nNeu100->SetName("gr_mGlu2000nNeu100");
    
    gr_mGlu2000nNeu100 ->SetMarkerStyle(kFullDiamond);
    gr_mGlu2000nNeu100 ->SetLineColor(kMagenta);
    gr_mGlu2000nNeu100 ->SetMarkerColor(kMagenta);
    //------------------------------------------
    sigscale = 35.882548708 * 0.981077 / 6060;

    TFile *f_mGlu2000nNeu1900 = new TFile("T5Wg_2000_1900_unMET.root","READ");
    TH1F  *h_mGlu2000nNeu1900 = (TH1F*)f_mGlu2000nNeu1900->Get("h_T5Wg_2000_1900");
    
    h_mGlu2000nNeu1900=getOverflow(h_mGlu2000nNeu1900);
    
    for(int i = 0; i < NbinsMET;i++){
        float binwidth = h_mGlu2000nNeu1900->GetBinWidth(i+1);
        sigx[i] = h_mGlu2000nNeu1900->GetBinCenter(i+1);
        sigy[i] = h_mGlu2000nNeu1900->GetBinContent(i+1)*sigscale/binwidth;
        sigxel[i] = binwidth/2;
        sigxeu[i] = binwidth/2;
        sigyel[i] = h_mGlu2000nNeu1900->GetBinErrorUp(i+1)*sigscale/binwidth;
        sigyeu[i] = h_mGlu2000nNeu1900->GetBinErrorLow(i+1)*sigscale/binwidth;
    }
    TGraphAsymmErrors * gr_mGlu2000nNeu1900 = new TGraphAsymmErrors(NbinsMET,sigx,sigy,sigxel,sigxeu,sigyel,sigyeu);
    gr_mGlu2000nNeu1900->SetName("gr_mGlu2000nNeu1900");
    
    gr_mGlu2000nNeu1900 ->SetMarkerStyle(kFullTriangleDown);
    gr_mGlu2000nNeu1900 ->SetLineColor(kGreen);
    gr_mGlu2000nNeu1900 ->SetMarkerColor(kGreen);
    //------------------------------------------

    */
    
   // TFile *f1   = new TFile("finalQCD.root", "READ");
    TFile *f1   = new TFile("OutputFiles/finalBkg.root", "READ");
    TH1F *h_QCD = (TH1F*)f1->Get("finalQCDMet");
    TH1F *h_EWK = (TH1F*)f1->Get("finalEWKMet");
    TH1F *h_DoublePhoton_MET;
    TH1F  *h_ggZGGToNuNuGG=(TH1F*)f1->Get("hZGGMet");

    TH1F * tempGG = (TH1F*)f1->Get("ggMet");

    if(blinded){
        const int numMetBins = 18;
        Double_t MetBins[numMetBins+1]  = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,75,85,100};
        TH1F * tempGG = (TH1F*)f1->Get("ggMet");
        h_DoublePhoton_MET = new TH1F("h_DoublePhoton_MET","",numMetBins,MetBins);
        for(int i =0; i <= numMetBins;i++){
            float value = tempGG->GetBinContent(i);
            float error = (tempGG->GetBinErrorUp(i) + tempGG->GetBinErrorLow(i))/2;
            h_DoublePhoton_MET->SetBinContent(i,value);
            h_DoublePhoton_MET->SetBinError(i,error);
        }
        
    }
    else{
        h_DoublePhoton_MET =(TH1F*)f1->Get("ggMet");
    }
    TGraphAsymmErrors *gr_met_error = (TGraphAsymmErrors*)f1->Get("finalBkg");
    
   // TFile *fZGG = new TFile("ZGGToNuNuGGV2_analyzed.root");
    //TH1F  *h_ggZGGToNuNuGG=(TH1F*)fZGG->Get("h_gg_MET");
    
    float s_ZGGToNuNuGG = (lumi*1000*0.07477)/994078;
   // h_ggZGGToNuNuGG ->Scale(s_ZGGToNuNuGG);

    
    
    //make Stack Plot
    
    THStack *hs = new THStack("hs",";;Events/GeV");

    int offset = -3;
    h_QCD->SetFillColor(kRed+offset);
    h_QCD->SetLineColor(kRed+offset);
    
    h_EWK->SetLineColor(kBlue+offset);
    h_EWK->SetFillColor(kBlue+offset);
    
    h_ggZGGToNuNuGG -> SetLineColor(kGreen+offset);
    h_ggZGGToNuNuGG -> SetFillColor(kGreen+offset);

    hs->Add(h_ggZGGToNuNuGG);
    hs->Add(h_EWK);
    hs->Add(h_QCD);
    
    
    //hs->Add(h_QCDwithTTJ);
    
    h_DoublePhoton_MET->SetMarkerStyle(kFullCircle);
    h_DoublePhoton_MET->SetMarkerSize(1.0);
    
    h_DoublePhoton_MET->SetLineColor(1);
    h_DoublePhoton_MET->SetMarkerColor(1);
    h_DoublePhoton_MET->GetXaxis()->SetRangeUser(xlow,xup);
    
 //   gr_met_error->SetFillColor(kBlack);
    gr_met_error->SetFillColor(kOrange);
    gr_met_error->SetLineColor(kOrange);
    gr_met_error->SetFillStyle(3144);
    
    
    
 /*   h_QCD	      ->Scale(1.0, "width");
    h_EWK	      ->Scale(1.0, "width");
    h_DoublePhoton_MET->Scale(1.0, "width");
    h_ggZGGToNuNuGG   ->Scale(1.0, "width");*/
    
    int H_ref = 850;
    int W_ref = 900;
    
    // references for T, B, L, R
    float T = 0.08*H_ref;
    float B = 0.12*H_ref;
    float L = 0.12*W_ref;
    float R = 0.15*W_ref;
    
    
    TCanvas *ctemp = new TCanvas("canv1","canv1",50,50,W_ref,H_ref);
    ctemp->cd();
    ctemp->SetFillColor(0);
    ctemp->SetBorderMode(0);
    ctemp->SetFrameFillStyle(0);
    ctemp->SetFrameBorderMode(0);
    ctemp->SetLeftMargin( L/W_ref );
    ctemp->SetRightMargin( R/W_ref );
    ctemp->SetTopMargin( T/H_ref );
    ctemp->SetBottomMargin( B/H_ref );
    ctemp->cd();
    
    //gPad->SetLogy();
    
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    pad1->SetFrameLineWidth(3);

    pad1->SetBottomMargin(0);
    gStyle->SetOptStat(0);
    
    hs->Draw("Histo");
    hs->GetXaxis()->SetTitleSize(0.12);
    hs->GetYaxis()->SetTitleSize(0.061);
    hs->GetYaxis()->SetLabelSize(0.058);
//    hs->GetYaxis()->CenterTitle();
    hs->GetYaxis()->SetTitleOffset(0.87);
    
    TLine *line2 = new TLine(100,0.0002,100,9000);
    line2->SetLineStyle(2);
    line2->Draw();

    gr_met_error->Draw("e2 same");
    h_DoublePhoton_MET->Draw("same ep");
    h_mGlu1700nNeu1000->Draw("pe same");
    h_mGlu2000nNeu1000->Draw("pe same");
    
    cout << "Asymm errs on background: +" << gr_met_error->GetErrorYhigh(23) << ", -" << gr_met_error->GetErrorYlow(23) << endl;
    cout << "Stat errs on data:        +" << h_DoublePhoton_MET->GetBinErrorUp(24) << ", -" << h_DoublePhoton_MET->GetBinErrorLow(24) << endl;
    gr_met_error->Print();
    
//    gr_met_error->Draw("e2 same");
    
    hs->SetMinimum(0.001);
    hs->SetMaximum(4000);
    
   /* gr_mGlu2000nNeu1000->Draw("same PZ");
    gr_mGlu2000nNeu100 ->Draw("same PZ");
    gr_mGlu2000nNeu1900->Draw("same PZ");
    */
    TLegend *leg = new TLegend(0.57,0.46,0.9,0.9); // cms wants 0.5,0.6,0.9,0.9
    leg->SetFillColor(kWhite);
    leg->SetTextFont(42); // cms wants 42
    leg->SetBorderSize(0);
    leg->SetShadowColor(kWhite);
    leg->SetFillStyle(0);
    leg->AddEntry(h_DoublePhoton_MET,"Observed","lep");
    leg->AddEntry(h_QCD,"QCD","f");
    leg->AddEntry(h_EWK,"EWK","f");
    leg->AddEntry(h_ggZGGToNuNuGG,"Z#gamma#gamma","f");
    leg->AddEntry(gr_met_error, "Syst + stat uncert","f");
    leg->AddEntry(h_mGlu1700nNeu1000,"T5gg 1700, 1000","lp");
    leg->AddEntry(h_mGlu2000nNeu1000,"T5gg 2000, 1000","lp");


 /*   leg->AddEntry(gr_mGlu2000nNeu100 ,"T5Wg 2000, 100 ","P");
    leg->AddEntry(gr_mGlu2000nNeu1000,"T5Wg 2000, 1000","P");
    leg->AddEntry(gr_mGlu2000nNeu1900,"T5Wg 2000, 1900","P");*/
    leg->Draw();
    
    
    TLatex *tex1;
    TLatex *tex2;
    TLatex *tex3;
    FinalTexMaker(tex1,tex2,tex3);
    
    TLine *line5 = new TLine(350,0.0001,350,1);
    line5->SetLineWidth(3);
    line5->Draw();
    
    ctemp->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.00,1,0.3);
    
    pad2->Draw();
    pad2->cd();
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.34);
    pad2->SetFrameLineWidth(3);

    TH1F *h2 = (TH1F*)h_DoublePhoton_MET->Clone("h2");
    h2->Reset();
    h2->GetYaxis()->SetTitle("#gamma#gamma/bkg");
    h2->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)");
    h2->GetXaxis()->SetTitleOffset(0.82);
    
    TH1F *h3 = (TH1F*)tempGG->Clone("h3");
    h3->Reset();
    h3->GetYaxis()->SetTitle("#gamma#gamma/bkg");
    h3->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)");
    h3->GetXaxis()->SetTitleOffset(0.82);
    
    const int nLoop = 25;
    Double_t y[nLoop-1],
    yel[nLoop-1],
    yeu[nLoop-1];
    Double_t y2[nLoop-1],
    y2el[nLoop-1],
    y2eu[nLoop-1];
    
    Double_t * x = gr_met_error->GetX();
    Double_t * xel = gr_met_error->GetEXhigh();
    Double_t * xeu = gr_met_error->GetEXlow();
    
    TH1F * hackPoisson = new TH1F("hackPoisson","",110,-0.5,109.5);
    hackPoisson->SetBinErrorOption(TH1::kPoisson);
    
    for(int i = 0; i < 110; i++){
        for (int j = 0; j < i; j++ ) hackPoisson->Fill(i);
    }
    
    for(int i=0;i<=tempGG->GetNbinsX();++i){
        double qcd  = h_QCD->GetBinContent(i+1);
        double ewk  = h_EWK->GetBinContent(i+1);
        double ZGG  = h_ggZGGToNuNuGG->GetBinContent(i+1);
        double central = qcd + ewk + ZGG;
        double sys = gr_met_error->GetErrorY(i)/(central);
        //cout << gr_met_error->GetErrorY(i) << endl;
        double err_bkg = gr_met_error->GetErrorY(i);
        
        if(i==22) cout << h_QCD->GetBinCenter(i) << " " << sys << endl;

        h3->SetBinContent(i+1,1);
        //h3->SetBinError(i+1,sys);
        h3->SetBinError(i+1,0.0);
        
        y[i] = 1.0;
        if(central > 0.0){
            yeu[i] = (gr_met_error->GetErrorYhigh(i)/central);
            yel[i] = (gr_met_error->GetErrorYlow(i)/central);
        }
        if(i <= h_DoublePhoton_MET->GetNbinsX()){
            double yi = h_DoublePhoton_MET->GetBinContent(i+1);
            double erry = h_DoublePhoton_MET->GetBinError(i+1);
            
            double erz(0);
            if(yi!=0 && central!=0) erz= (erry/central);//erz= (yi/central)*sqrt((erry/yi)*(erry/yi))+(err_bkg/central)*(err_bkg/central));
            
            if(central!=0)
            {
                h2->SetBinContent(i+1,yi/central);
                y2[i] = yi/central;
            }
            h2->SetBinError(i+1,0.0);
            //cout << erry/central << endl;
            
            if(i >= 18){
                int unweighted = 0.0;
                if(i == 18) unweighted = 105;
                else if(i == 19) unweighted = 39;
                else if(i == 20) unweighted = 21;
                else if(i == 21) unweighted = 21;
                else if(i == 22) unweighted = 11;
                else if(i == 23) unweighted = 12;
                
                //cout << "int?" << unweighted << " vs " << (yi*h_DoublePhoton_MET->GetBinWidth(i+1)) << endl;
                y2el[i] = hackPoisson->GetBinErrorLow(unweighted)/unweighted*yi/central;
                y2eu[i] = hackPoisson->GetBinErrorUp(unweighted)/unweighted*yi/central;
                cout << "center = " << yi/central << ", poisson error on " << unweighted << " is -" << hackPoisson->GetBinErrorLow(unweighted) << ", + " << hackPoisson->GetBinErrorUp(unweighted) <<" and the total error is then " << y2el[i] << endl;
            }
            else{
                y2el[i] = erz;
                y2eu[i] = erz;
            }
        }
    }
    
    TGraphAsymmErrors * hshaded = new TGraphAsymmErrors(nLoop-1,x,y,xel,xeu,yel,yeu);
    TGraphAsymmErrors * hratio = new TGraphAsymmErrors(nLoop-1,x,y2,xel,xeu,y2el,y2eu);
    hratio->SetLineWidth(2);
    
    cout << "x: " << x[23] << ", +" << xeu[23] << ", -" << xel[23] << endl;
    cout << "y:" << y[23] << ", +" << yeu[23] << ", -" << yel[23] <<endl;

    cout << "x: " << x[21] << ", +" << xeu[21] << ", -" << xel[21] << endl;
    cout << "y:" << y[21] << ", +" << yeu[21] << ", -" << yel[21] <<endl;

    
       cout << "h2:" << h2->GetBinErrorUp(24)  << " " << h2->GetBinErrorLow(24) << " " << h2->GetBinError(24) << " " << h2->GetBinCenter(24) << endl;
    cout << "h3:" << h3->GetBinErrorUp(24)  << " " << h3->GetBinErrorLow(24) << " " << h3->GetBinError(24) << " " << h3->GetBinCenter(24) << endl;
    
    h3->GetXaxis()->SetRangeUser(xlow,xup);
    h3->GetXaxis()->SetLabelSize(0.13);
    h3->GetXaxis()->SetTitleSize(0.14);
    h3->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)");
    h3->GetXaxis()->SetTitleOffset(1.0);

    
    h3->SetTitle("");
    h3->GetYaxis()->SetRangeUser(0.0,3.25);
    h3->GetYaxis()->SetNdivisions(6);
    h3->GetXaxis()->SetTickLength(0.07);

    h3->GetYaxis()->SetTitle("Obs/Bkgd");
    h3->GetYaxis()->SetTitleOffset(0.34);
    h3->GetYaxis()->SetLabelSize(0.13);
    h3->GetYaxis()->SetTitleSize(0.14);
    h3->GetYaxis()->CenterTitle();
    hshaded->SetFillStyle(1001);
    hshaded->SetFillColor(kGray);
    hshaded->SetLineColor(kGray);
    hshaded->SetMarkerColor(kGray);
    

    TLine *line1 = new TLine(xlow,1,xup,1);
    line1->SetLineStyle(2);
    h3->Draw("e2");
    hshaded->Draw("e2 same");
    h2->Draw("ep same");
    hratio->Draw("pZ same");
    //leg2->Draw("same");
    line1->Draw();
    
    TLine *line3 = new TLine(100,0.0,100,3.25);
    line3->SetLineStyle(2);
    line3->Draw();
    
    TLine *line4 = new TLine(350,0.0,350,3.25);
    line4->SetLineWidth(3);
    line4->Draw();
    
    
    //h3->Draw("sameaxis");
    
   /* TLegend *leg2 = new TLegend(0.78,0.77,0.85,0.91);
    leg2->SetText
    
    leg2->SetHeader("data/bkg");
    leg2->Draw();*/
    
    ctemp->SaveAs("plots/stackPlot.pdf","pdf");
   
}
