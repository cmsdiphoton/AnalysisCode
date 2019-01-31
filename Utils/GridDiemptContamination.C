#include <map>
#include <list>
#include <vector>
#include <iostream>


float getXSec(float mGlu);
float getXSecSquark(float mSquark);
TCanvas * SetUpCanvas(TString canvName);
void plotXSec();

void GridDiemptContamination(){
    
    TFile *infile = TFile::Open("hist_T5Wg_THESIS.root");
   // TFile *infile = TFile::Open("/Users/alliereinsvold/MyWorkingDirectory/FinalPlots2016/Diempt/hist_T5Wg_DiemptGrid.root");
    infile->cd();
    
    Double_t ybins[123] = {0, 17.5, 37.5, 75, 125, 175, 250, 350, 450, 550, 650, 750, 850, 950,
        1025.001, 1050.001, 1075.001, 1125.001, 1150.001, 1175.001,
        1225.001, 1250.001, 1262.5, 1275.001, 1282.5, 1300.001, 1312.5,
        1325.001, 1332.5, 1350.001, 1362.5, 1375.001, 1382.5, 1400.001, 1412.5,
        1425.001, 1432.5, 1450.001, 1462.5, 1475.001, 1482.5, 1500.001, 1512.5,
        1525.001, 1532.5, 1550.001, 1562.5, 1575.001, 1582.5, 1600.001, 1612.5,
        1625.001, 1632.5, 1650.001, 1662.5, 1675.001, 1682.5, 1700.001, 1712.5,
        1725.001, 1732.5, 1750.001, 1762.5, 1775.001, 1782.5, 1800.001, 1812.5,
        1825.001, 1832.5, 1850.001, 1862.5, 1875.001, 1882.5, 1900.001, 1912.5,
        1925.001, 1932.5, 1950.001, 1962.5, 1975.001, 1982.5, 2000.001, 2012.5,
        2025.001, 2032.5, 2050.001, 2062.5, 2075.001, 2082.5, 2100.001, 2112.5,
        2125.001, 2132.5, 2150.001, 2162.5, 2175.001, 2182.5, 2200.001, 2212.5,
        2225.001, 2232.5, 2250.001, 2262.5, 2275.001, 2282.5, 2300.001, 2312.5,
        2325.001, 2332.5, 2350.001, 2362.5, 2375.001, 2382.5, 2400.001, 2412.5,
        2425.001, 2432.5, 2450.001,2462.5, 2475.001, 2482.5, 2500.001, 2512.5};
    int ybinnum = 122;
    
    
    Double_t xbins3[26]={1275, 1325, 1375, 1425, 1475, 1525, 1575, 1625, 1675,
        1725, 1775, 1825, 1875, 1925, 1975, 2025, 2075, 2125, 2175, 2225, 2275,
        2325, 2375, 2425, 2475, 2525};
    int xbinnum3 = 25;
    
    Double_t red[5] = {0.5, 0.5, 1.0, 1.0, 1.0};
    Double_t green[5] = {0.5, 1.0, 1.0, 0.6, 0.5};
    Double_t blue[5] = {1.0, 1.0, 0.5, 0.4, 0.5};
    Double_t length[5] =  {0.0, 0.34, 0.61, 0.84, 1.0};
    
    TColor::CreateGradientColorTable(5, length, red, green, blue, 255);
    
    float lumi = 35921.8;
    float ggData = 473;
    float ffData = 164;
    float dataWeight = ggData/ffData;
    
    cout<< "Goal weight = " << dataWeight << endl;
    
    TH2D *ggDiEMPt300_GridScaled = new TH2D("ggDiEMPt300_GridScaled","Scaled events for DiEMPt > 300 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3,  xbins3, ybinnum, ybins);
    TH2D *ffDiEMPt300_GridScaled = new TH2D("ffDiEMPt300_GridScaled","Scaled events for DiEMPt > 300 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3,  xbins3, ybinnum, ybins);
    
    TH2D *DiemptReweight_Grid = new TH2D("DiemptReweight_Grid","Reweighting value for DiEMPt > 300 GeV; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3,  xbins3, ybinnum, ybins);
    
    TH2D *DiemptReweight_GridPercent = new TH2D("DiemptReweight_GridPercent","Percent diff between value with and without SUSY; Gluino Mass (GeV); Neutralino Mass (GeV)", xbinnum3, xbins3, ybinnum, ybins);
    
    float mGlu, mNeu, xsec, scaleByLumi, ffScaled, ggScaled, newWeight, percentDiff, bin;
    
    for(int i = 1; i <= xbinnum3;i++){
        for(int j = 1; j <= ybinnum;j++){
            
            if(ggDiEMPt300_Grid->GetBinContent(i,j) == 0) continue;
            
            mGlu = xbins3[i] - 25;
            mNeu = (ybins[j] + ybins[j-1]) / 2.0;
            
            if(mNeu >= mGlu) continue;
            xsec = getXSec(mGlu);
            if(xsec == 1){
                cout<<"Uh oh" << endl;
                break;
            }
            bin =h_GridAllEvents->FindFixBin(mGlu,mNeu);
            scaleByLumi = lumi*xsec/h_GridAllEvents->GetBinContent(bin);
            
            ffScaled = ffDiEMPt300_Grid->GetBinContent(i,j)*scaleByLumi;
            ffDiEMPt300_GridScaled->SetBinContent(i,j,ffScaled);
            
            ggScaled = ggDiEMPt300_Grid->GetBinContent(i,j)*scaleByLumi;
            ggDiEMPt300_GridScaled->SetBinContent(i,j,ggScaled);
            //newWeight = (ggScaled+ggData)/(ffData);
            newWeight = (ggScaled+ggData)/(ffScaled+ffData);
            DiemptReweight_Grid->SetBinContent(i,j,newWeight);
            
            percentDiff = TMath::Abs(newWeight - dataWeight)/dataWeight*100.0;
            DiemptReweight_GridPercent->SetBinContent(i,j,percentDiff);
            
// cout << mGlu << " " << mNeu << " " << ggScaled << " " << ffScaled << " " << dataWeight << " " << newWeight << " " << scaleByLumi << " " << percentDiff << endl;
            
        }
    }

    TCanvas *c1 = new TCanvas("c1","c1",700,700);
    c1->cd();
    gStyle->SetOptStat(0);
   // DiemptReweight_GridPercent->GetZaxis()->SetRangeUser(0.1,40);
    DiemptReweight_GridPercent->GetXaxis()->SetRangeUser(1400,2112.5);
    DiemptReweight_GridPercent->GetYaxis()->SetTitleOffset(1.5);
    DiemptReweight_GridPercent->Draw("colz");

    TCanvas *c2 = new TCanvas("c2","c2",700,700);
    c2->cd();
    gStyle->SetOptStat(0);
    //DiemptReweight_Grid->GetZaxis()->SetRangeUser(8,13);
    DiemptReweight_Grid->GetXaxis()->SetRangeUser(1400,2112.5);
    DiemptReweight_Grid->GetYaxis()->SetTitleOffset(1.5);
    DiemptReweight_Grid->Draw("colz");
    
    plotXSec();
    
    return;
}

void plotXSec(){
    TCanvas *c3 = SetUpCanvas("c3");
    c3->SetLogy(kTRUE);
    TH1F* xSecHist = new TH1F("",";Gluino/Squark mass (GeV); Cross section (fb)",18,1300,2200);
    for(int i = 1300; i < 2200;i=i+50){
        float value = 1000.*getXSec(i);
        xSecHist->Fill(i,value);
        int bin = xSecHist->FindFixBin(i);
        float error = getXSecError(i)*value/100;
        xSecHist->SetBinError(bin,error);
    }
    gStyle->SetOptStat(0);
    xSecHist->SetMarkerStyle(kFullTriangleDown);
    xSecHist->SetMarkerSize(1.5);
    xSecHist->SetLineColor(38);
    xSecHist->SetMarkerColor(38);
    xSecHist->GetYaxis()->SetRangeUser(0.05,65);
    xSecHist->GetYaxis()->SetTitleOffset(1.4);
    xSecHist->GetYaxis()->CenterTitle();
    xSecHist->GetXaxis()->CenterTitle();

/*    xSecHist->GetYaxis()->SetTitleSize(0.04);
    xSecHist->GetXaxis()->SetTitleSize(0.03);
    xSecHist->GetYaxis()->SetTitleSize(0.5);
    xSecHist->GetXaxis()->SetTitleSize(0.5);*/
    xSecHist->Draw("p");
    
    TH1F* xSecHistSquark = new TH1F("",";Gluino/Squark mass (GeV); Cross section (fb)",18,1300,2200);
    for(int i = 1300; i < 2200;i=i+50){
        float value = 1000.*getXSecSquark(i);
        xSecHistSquark->Fill(i,value);
        int bin = xSecHistSquark->FindFixBin(i);
        float error = getXSecSquarkError(i)*value/100;
        xSecHistSquark->SetBinError(bin,error);
    }
    xSecHistSquark->SetMarkerSize(1.5);
    xSecHistSquark->SetMarkerStyle(kFullTriangleUp);
    xSecHistSquark->SetLineColor(30);
    xSecHistSquark->SetMarkerColor(30);
    xSecHistSquark->Draw("p same");
    
    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(xSecHist,"Gluino pair production","p");
    leg->AddEntry(xSecHistSquark,"Squark pair production","p");
    leg->Draw("same");
    
    int   align_       = 10*1+1;
    float relPosX      = 0.045;
    float relPosY      = -0.01;
    float relExtraDX   = 1.5;
    
    TString cmsText     = "CMS";
    float cmsTextFont   = 61;  // default is helvetic-bold
    float cmsTextSize   = 0.6;
    
    TString extraText   = "Preliminary";
    float extraTextFont = 52;  // default is helvetica-italics
    float extraOverCmsTextSize  = 0.7;
    float extraTextSize = extraOverCmsTextSize*cmsTextSize;
    
    TString lumiText = "35.9 fb^{-1}";
    
    float r = c3->GetRightMargin();
    float l = c3->GetLeftMargin();
    float posX_ =   l + relPosX*(1-l-r);
    
    
    float t = c3->GetTopMargin();
    float b = c3->GetBottomMargin();
    float posY_ = 1-t - relPosY*(1-t-b);
    
    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);
    latex.SetTextFont(cmsTextFont);
    latex.SetTextSize(cmsTextSize*t);
    latex.SetTextAlign(align_);
    latex.DrawLatex(posX_, posY_, cmsText);
    if( true )
    {
        latex.SetTextFont(extraTextFont);
        latex.SetTextAlign(align_);
        latex.SetTextSize(extraTextSize*t);
        latex.DrawLatex(posX_ + relExtraDX*cmsTextSize*l, posY_, extraText);
    }
    latex.SetTextAlign(31);
    latex.DrawLatex(1-r, posY_, lumiText);
    
    gPad->Update();
    
    
}


float getXSec(float mGlu){
    if(mGlu == 1300){ return 0.0460525;}
    if(mGlu == 1350){ return 0.0340187;}
    if(mGlu == 1400){ return 0.0252977;}
    if(mGlu == 1450){ return 0.0188887;}
    if(mGlu == 1500){ return 0.0141903;}
    if(mGlu == 1550){ return 0.0107027;}
    if(mGlu == 1600){ return 0.00810078;}
    if(mGlu == 1650){ return 0.00616072;}
    if(mGlu == 1700){ return 0.00470323;}
    if(mGlu == 1750){ return 0.00359842;}
    if(mGlu == 1800){ return 0.00276133;}
    if(mGlu == 1850){ return 0.00212345;}
    if(mGlu == 1900){ return 0.00163547;}
    if(mGlu == 1950){ return 0.0012642;}
    if(mGlu == 2000){ return 0.000981077;}
    if(mGlu == 2050){ return 0.000761286;}
    if(mGlu == 2100){ return 0.000591918;}
    if(mGlu == 2150){ return 0.000460941;}
    if(mGlu == 2200){ return 0.000359318;}
    if(mGlu == 2250){ return 0.00028065;}
    if(mGlu == 2300){ return 0.000219049;}
    if(mGlu == 2350){ return 0.000171031;}
    if(mGlu == 2400){ return 0.000133965;}
    if(mGlu == 2450){ return 0.000104886;}
    if(mGlu == 2500){ return 0.0000820068;}
    
    return 1;
}

float getXSecSquark(float mSquark){
    if(mSquark ==1300){ return 0.0086557;}
    if(mSquark ==1350){ return 0.00637816;}
    if(mSquark ==1400){ return 0.00472062;}
    if(mSquark ==1450){ return 0.00350484;}
    if(mSquark ==1500){ return 0.00261222;}
    if(mSquark ==1550){ return 0.0019636;}
    if(mSquark ==1600){ return 0.00148553;}
    if(mSquark ==1650){ return 0.00112118;}
    if(mSquark ==1700){ return 0.000846524;}
    if(mSquark ==1750){ return 0.000646271;}
    if(mSquark ==1800){ return 0.000494604;}
    if(mSquark ==1850){ return 0.000378227;}
    if(mSquark ==1900){ return 0.000289986;}
    if(mSquark ==1950){ return 0.000223058;}
    if(mSquark ==2000){ return 0.000171615;}
    if(mSquark ==2050){ return 0.000132986;}
    if(mSquark ==2100){ return 0.000102747;}
    if(mSquark ==2150){ return 0.000079509;}
    
    return 1;
}

TCanvas * SetUpCanvas(TString canvName){
    int H_ref = 850;
    int W_ref = 900;
    
    // references for T, B, L, R
    float T = 0.08*H_ref;
    float B = 0.12*H_ref;
    float L = 0.12*W_ref;
    float R = 0.04*W_ref;
    
    
    ctemp = new TCanvas(canvName,canvName,50,50,W_ref,H_ref);
    ctemp->cd();
    ctemp->SetFillColor(0);
    ctemp->SetBorderMode(0);
    ctemp->SetFrameFillStyle(0);
    ctemp->SetFrameBorderMode(0);
    ctemp->SetLeftMargin( L/W_ref );
    ctemp->SetRightMargin( R/W_ref );
    ctemp->SetTopMargin( T/H_ref );
    ctemp->SetBottomMargin( B/H_ref );
    ctemp->SetTickx(0);
    ctemp->SetTicky(0);
   // ctemp->SetGrid();
    //gStyle->SetGridStyle(3);
    
    
    return ctemp;
    
}

float getXSecSquarkError(float mSquark){
    if(mSquark == 1300) {return 19.8532;}
    if(mSquark == 1350) {return 20.8815;}
    if(mSquark == 1400) {return 21.9859;}
    if(mSquark == 1450) {return 23.0346;}
    if(mSquark == 1500) {return 24.1631;}
    if(mSquark == 1550) {return 25.3454;}
    if(mSquark == 1600) {return 26.6472;}
    if(mSquark == 1650) {return 28.1032;}
    if(mSquark == 1700) {return 29.1726;}
    if(mSquark == 1750) {return 30.2417;}
    if(mSquark == 1800) {return 31.3049;}
    if(mSquark == 1850) {return 32.6422;}
    if(mSquark == 1900) {return 34.0784;}
    if(mSquark == 1950) {return 35.6109;}
    if(mSquark == 2000) {return 36.8769;}
    if(mSquark == 2050) {return 38.4295;}
    if(mSquark == 2100) {return 39.6182;}
    if(mSquark == 2150) {return 40.9212;}
}

float getXSecError(float mGlu){
if(mGlu == 1300){ return 19.64;}
if(mGlu == 1350){ return 20.3088;}
if(mGlu == 1400){ return 20.9163;}
if(mGlu == 1450){ return 21.9548;}
if(mGlu == 1500){ return 22.7296;}
if(mGlu == 1550){ return 23.4971;}
if(mGlu == 1600){ return 24.2679;}
if(mGlu == 1650){ return 25.138;}
if(mGlu == 1700){ return 26.1021;}
if(mGlu == 1750){ return 27.1502;}
if(mGlu == 1800){ return 28.108;}
if(mGlu == 1850){ return 28.9167;}
if(mGlu == 1900){ return 29.9045;}
if(mGlu == 1950){ return 30.4581;}
if(mGlu == 2000){ return 31.8422;}
if(mGlu == 2050){ return 32.9341;}
if(mGlu == 2100){ return 33.9326;}
if(mGlu == 2150){ return 34.9082;}
}

