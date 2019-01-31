{
    
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
    
    //** commented lines = how plot was made in November 2017 to go in AN
    Double_t red[5] = {0.5, 0.5, 1.0, 1.0, 1.0};
    Double_t green[5] = {0.5, 1.0, 1.0, 0.6, 0.5};
    Double_t blue[5] = {1.0, 1.0, 0.5, 0.4, 0.5};
    Double_t length[5] =  {0.0, 0.34, 0.61, 0.84, 1.0};
    
    TColor::CreateGradientColorTable(5, length, red, green, blue, 255);
    
    /*TFile *fileIn = TFile::Open("hist_T5Wg.root");
  h_GridAllGGEvents->Draw("colz");
    h_GridAllGGEvents->SetStats(kFALSE);
    TH1F* acc = h_GridAllGGEvents->Clone("acc");
    acc->Scale(4);
    h_GridAllEvents->Draw();
    acc->Divide(h_GridAllEvents);
    acc->SetTitle("Acceptance x Efficiency;Gluino Mass (GeV);Neutralino Mass (GeV)");
    acc->GetYaxis()->SetTitleOffset(1.2);
    /*acc->Draw("colz");*/
   
    
    TFile *fileIn = TFile::Open("hist_T5Wg_ext_THESIS.root");
    fileIn->cd();
    h_GridAllGGEvents->Draw("colz");
    h_GridAllEvents->Draw();
    
    TFile *fileIn2 = TFile::Open("hist_T5Wg_THESIS.root");
    fileIn2->cd();
    h_GridAllGGEvents->Draw("colz");
    h_GridAllEvents->Draw();
    
    TFile *fileResults = TFile::Open("accXeff.root","RECREATE");
    fileResults->cd();
    
    TH2F* hAll = fileIn->FindObject("h_GridAllEvents");
    TH2F* hAllGG = fileIn->FindObject("h_GridAllGGEvents");
  
    TH2F* hAll_ext   = fileIn2->FindObject("h_GridAllEvents");
    TH2F* hAllGG_ext = fileIn2->FindObject("h_GridAllGGEvents");
    
    hAllGG_ext->SetStats(kFALSE);
    TH2F* acc = hAllGG_ext->Clone("acc");
    acc->Divide(hAll_ext);
    acc->GetZaxis()->SetRangeUser(0.1,0.3);
    acc->Draw("colz");
//**    acc->Scale(4);
    acc->SetTitle(";Gluino Mass (GeV);Neutralino Mass (GeV)");
    acc->GetYaxis()->SetTitleOffset(1.2);
    
    for(int i = 1; i <= hAll->GetXaxis()->GetNbins(); i++){
        for(int j = 1; j <= hAll->GetYaxis()->GetNbins(); j++){
            float tot = hAll->GetBinContent(i,j);
            if(tot == 0) continue;
            float gmass = hAll->GetXaxis()->GetBinCenter(i);
            float nmass = hAll->GetYaxis()->GetBinCenter(j);
            //** float weight = 4*hAllGG->GetBinContent(i,j)/tot;
            float weight = hAllGG->GetBinContent(i,j)/tot;
            acc->Fill(gmass, nmass, weight);
        }
    }
    
    Double_t xbins3[26]={1275,1325,1375, 1425, 1475, 1525, 1575, 1625, 1675, 1725, 1775, 1825, 1875, 1925, 1975, 2025, 2075, 2125, 2175, 2225, 2275, 2325
        , 2375, 2425, 2475, 2525};
    int xbinnum3 = 25;
    
    Double_t ybins[123] = {0,17.5,37.5,75,125,175,250,
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
    int ybinnum = 122;
    
    for(int i =0; i <xbinnum3; i++){
        for(int j =2; j < ybinnum; j++){
            float tot = acc->GetBinContent(i,j);
            if(tot == 0 ){
                float gmass = acc->GetXaxis()->GetBinLowEdge(i) + acc->GetXaxis()->GetBinWidth(i) - 25.0;
                float nmass = acc->GetYaxis()->GetBinLowEdge(j) + acc->GetYaxis()->GetBinWidth(j);
                if(gmass+0.001 < nmass) continue;
                
                float diff = gmass - nmass;
                float newBin_J = -27;
                
                float mod    = fmod(diff, 50);
                float mod100 = fmod(diff, 100);
                
                if(diff == 17.5) newBin_J = j - 1;
                else if (TMath::Abs(nmass - 1175) < 0.002 || TMath::Abs(nmass - 1075) < 0.002 ){
                    if(diff<300) newBin_J = j - 1;
                    else         newBin_J = j + 1;
                    
                }
                else if(diff < 320){
                    if( TMath::Abs(mod- 24.999 ) < 0.002) newBin_J = j - 2;
                    if(mod == 37.5)  newBin_J  = j - 1;
                    if(mod == 17.5)  newBin_J  = j + 1;
                }
                else if( nmass <= 1251 && TMath::Abs(mod - 49.999) < 0.02 ){
                    newBin_J = j - 1;
                }
                else if(fmod(gmass,100) ==0){
                    if      (TMath::Abs(mod100 - 49.999 ) < 0.002) newBin_J = j -4;
                    else if (TMath::Abs(mod100 - 67.5 )   < 0.002) newBin_J = j -3;
                    else if (TMath::Abs(mod100 - 74.999 ) < 0.002) newBin_J = j -2;
                    else if (TMath::Abs(mod100 - 87.5 )   < 0.002) newBin_J = j -1;
                    else if (TMath::Abs(mod100 - 17.5 )   < 0.002) newBin_J = j +1;
                    else if (TMath::Abs(mod100 - 24.999 ) < 0.002) newBin_J = j +2;
                    else if (TMath::Abs(mod100 - 37.5 )   < 0.002) newBin_J = j +3;
                    
                }
                else{
                    if      (TMath::Abs(mod100 - 99.999 ) < 0.002) newBin_J = j -4;
                    else if (TMath::Abs(mod100 - 17.5 )   < 0.002) newBin_J = j -3;
                    else if (TMath::Abs(mod100 - 24.999 ) < 0.002) newBin_J = j -2;
                    else if (TMath::Abs(mod100 - 37.5 )   < 0.002) newBin_J = j -1;
                    else if (TMath::Abs(mod100 - 67.5 )   < 0.002) newBin_J = j +1;
                    else if (TMath::Abs(mod100 - 74.999 ) < 0.002) newBin_J = j +2;
                    else if (TMath::Abs(mod100 - 87.5 )   < 0.002) newBin_J = j +3;
                }
                
                
                
                
                if(newBin_J == -27) cout<< "No new bin found. gmass = " << gmass << ", nmass = " << nmass <<", diff = " << diff << ", mod = " << mod << endl;
                else{
                    acc->SetBinContent(i,j,acc->GetBinContent(i,newBin_J) );
                }
            }
        }
    }
    acc->GetYaxis()->CenterTitle();
    acc->GetXaxis()->CenterTitle();
    acc->GetZaxis()->CenterTitle();

    acc->GetYaxis()->SetTitleSize(0.042);
    acc->GetXaxis()->SetTitleSize(0.042);
    acc->GetZaxis()->SetTitle("A #times #epsilon");
    acc->GetZaxis()->SetTitleSize(0.042);
    acc->GetZaxis()->SetLabelSize(0.039);
    acc->GetYaxis()->SetLabelSize(0.039);
    acc->GetXaxis()->SetLabelSize(0.039);
    
    acc->GetZaxis()->SetTitleOffset(1.1);
    acc->GetYaxis()->SetRangeUser(0.0,2300);
    acc->GetXaxis()->SetRangeUser(1275,2000);
    //acc->Smooth();
    
    
    acc->Draw("colz");
    
    
    TString cmsText     = "CMS";
    float cmsTextFont   = 61;  // default is helvetic-bold
    float cmsTextSize   = 0.8;
    
    TString extraText   = "Preliminary";
    float extraTextFont = 52;  // default is helvetica-italics
    float extraOverCmsTextSize  = 0.7;
    float extraTextSize = extraOverCmsTextSize*cmsTextSize;
    
    TString lumiText = "35.9 fb^{-1}";
    
    float relPosX      = 0.045;
    float relPosY      = +0.06;
    
    float r = ctemp->GetRightMargin();
    float l = ctemp->GetLeftMargin();
    float posX_ =   l + relPosX*(1-l-r);
    
    float t = ctemp->GetTopMargin();
    float b = ctemp->GetBottomMargin();
    float posY_ = 1-t - relPosY*(1-t-b);
    
    int   align_       = 10*1+1;
    float relExtraDX   = 1.4;
    
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
    latex.DrawLatex(1-r*1.1, posY_, lumiText);
}
