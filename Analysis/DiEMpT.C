TH1F * scaleWidth(TH1F* h1);
TH1F * scaleHist(TH1F* h1, float scale);

void DiEMpT(){
    //Root macro for making reweights.root for diempt reweighting procedure
    //TFile *inputFile = TFile::Open("hist_data.root");
    //TFile *inputFile = TFile::Open("unweighted_njet2.root");
    TFile *inputFile = TFile::Open("hist_data_NewSkims_rw.root");
    inputFile->cd();
    eeDiEMPt->Draw();
    ggDiEMPt->Draw();
    ffDiEMPt->Draw();
    gfDiEMPt->Draw();
    
    TFile* fileResults = TFile::Open("reweights.root","RECREATE");
    fileResults->cd();
    
    bool toscale = true;
    bool scalebybin = true;
    
    int H_ref = 850;
    int W_ref = 900;
    
    // references for T, B, L, R
    float T = 0.08*H_ref;
    float B = 0.12*H_ref;
    float L = 0.12*W_ref;
    float R = 0.04*W_ref;
    
    int   align_       = 10*1+1;
    float relPosX      = 0.045;
    float relPosY      = -0.01;
    float relExtraDX   = 1.5;
    
    TString cmsText     = "CMS";
    float cmsTextFont   = 61;  // default is helvetic-bold
    float cmsTextSize   = 0.5;
    
    TString extraText   = "Preliminary";
    float extraTextFont = 52;  // default is helvetica-italics
    float extraOverCmsTextSize  = 0.7;
    float extraTextSize = extraOverCmsTextSize*cmsTextSize;
    
    TString lumiText = "35.9 fb^{-1}";
    
    //const int numBins = 27;
    //float diemptbins[numBins+1] = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,70,85,100,115, 130, 150,175,200,300,400,600,1000};
    
    
    const int numBins = 25;
    float diemptbins[numBins+1] = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,70,85,100,115, 130, 150,175,200,300,1000};
    
    TCanvas* canv1 = new TCanvas("diempt","diempt",50,50,W_ref,H_ref);
    canv1->cd();
    canv1->SetFillColor(0);
    canv1->SetBorderMode(0);
    canv1->SetFrameFillStyle(0);
    canv1->SetFrameBorderMode(0);
    canv1->SetLeftMargin( L/W_ref );
    canv1->SetRightMargin( R/W_ref );
    canv1->SetTopMargin( T/H_ref );
    canv1->SetBottomMargin( B/H_ref );
    canv1->SetTickx(0);
    canv1->SetTicky(0);
    canv1->SetGrid();
    gStyle->SetGridStyle(3);
    
    TPad *p1a = new TPad("p1a","",0.,.3,1.,1.);
    p1a->SetBottomMargin(0);
    p1a->Draw();
    TPad *p1b = new TPad("p1b","",0.,0.,1.,.3);
    p1b->SetTopMargin(0);
    p1b->SetBottomMargin(0.2);
    p1b->Draw();
    p1a->cd();
    p1a->SetGridx();
    p1a->SetGridy();
    gPad->SetLogy();
    gPad->Update();
    
    TH1F* hggDiempt = inputFile->FindObject("ggDiEMPt");
    hggDiempt->SetLineColor(kRed);
    hggDiempt->SetMarkerStyle(kFullCircle);
    hggDiempt->SetMarkerSize(0.6);

    hggDiempt->SetMarkerColor(kRed);
    
    hggDiempt->Draw();
    gPad->Update();
    TPaveStats* sgg = (TPaveStats*) hggDiempt->GetListOfFunctions()->FindObject("stats");
    sgg->SetLineColor(kRed);
    sgg->SetX2NDC(.9);
    sgg->SetX1NDC(.7);
    sgg->SetY1NDC(.53);
    sgg->SetY2NDC(.68);
    
    TH1F* hffDiempt = inputFile->FindObject("ffDiEMPt");
    hffDiempt->SetLineColor(kBlack);
    hffDiempt->SetMarkerStyle(kFullCircle);
    hffDiempt->SetMarkerColor(kBlack);
    hffDiempt->SetMarkerSize(0.6);

    hffDiempt->SetStats(kTRUE);
    TH1F * hff_reweights = new TH1F("ff_reweights"," ;DiEMPt(GeV);gg / ff", numBins,diemptbins);
    hff_reweights->SetStats(0);
    hff_reweights->SetMarkerStyle(kFullCircle);
    hff_reweights->SetMarkerColor(kBlack);
    hff_reweights->SetMarkerSize(0.6);
    
    float num,denom;
    if(toscale){
        float ffscale =1.0 / hffDiempt->Integral();
        hffDiempt = scaleHist(hffDiempt, ffscale);
        float ggscale =1.0 / hggDiempt->Integral();
        hggDiempt = scaleHist(hggDiempt, ggscale);
    }
    
    if(scalebybin){
        hffDiempt = scaleWidth(hffDiempt);
        hggDiempt = scaleWidth(hggDiempt);
    }
    
    float fakemet, ggDiempt, wvalue, binvalue;
    float errorff, errorgg, error;
    for(int i =0; i <28; i++){
        ggDiempt = hggDiempt->GetBinContent(i);
        ffDiempt = hffDiempt->GetBinContent(i);
        errorgg = hggDiempt->GetBinError(i);
        errorff = hffDiempt->GetBinError(i);
        if(ffDiempt == 0 || ggDiempt ==0){
            wvalue = 0;
            error = 0;
        }
        else{
            wvalue = ggDiempt/ffDiempt;
            error = wvalue*sqrt((errorgg/ggDiempt)*(errorgg/ggDiempt)+(errorff/ffDiempt)*(errorff/ffDiempt));
        }
        binvalue= hff_reweights->GetBinCenter(i);
        hff_reweights->Fill(binvalue, wvalue);
        hff_reweights->SetBinError(i, error);
    }
    
    
    hffDiempt->SetTitle("ff DiEMPt;;Events");
    hffDiempt->GetXaxis()->SetRangeUser(0,1000);
    //hffDiempt->GetYaxis()->SetRangeUser(-5,6000);
    hffDiempt->GetYaxis()->SetLabelSize(0.04);
    hffDiempt->Draw();
    
    gPad->Update();
    TPaveStats* s2 =(TPaveStats*) hffDiempt->GetListOfFunctions()->FindObject("stats");
    s2->SetLineColor(kBlack);
    s2->SetX2NDC(.9);
    s2->SetX1NDC(.7);
    s2->SetY1NDC(.36);
    s2->SetY2NDC(.51);
    
    
    
    hffDiempt->Draw("e");
    hggDiempt->Draw("e same");
    
    //Create a Legend for the plot
    TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
    
    leg->SetTextSize(.04);
    leg->AddEntry(hggDiempt, "gg","lp");
    leg->AddEntry(hffDiempt, "ff", "lp");
    leg->Draw();
    gPad->Update();
    
    float r = p1a->GetRightMargin();
    float l = p1a->GetLeftMargin();
    float posX_ =   l + relPosX*(1-l-r);
    
    
    float t = p1a->GetTopMargin();
    float b = p1a->GetBottomMargin();
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
    
    
    p1b->cd();
    
    p1b->SetGridx();
    p1b->SetGridy();
    hff_reweights->GetXaxis()->SetRangeUser(0,1000);
    hff_reweights->SetLineColor(kBlack);
    hff_reweights->GetXaxis()->SetLabelSize(0.08);
    hff_reweights->GetYaxis()->SetLabelSize(0.08);
    hff_reweights->GetXaxis()->SetTitleSize(0.08);
    hff_reweights->GetYaxis()->SetTitleOffset(0.5);
    hff_reweights->GetYaxis()->CenterTitle();
    hff_reweights->GetYaxis()->SetTitleSize(0.08);
    hff_reweights->GetYaxis()->SetRangeUser(0,12.0);
    hff_reweights->Draw();
    
    TLine *line = new TLine(0,1,1000,1);
    line->SetLineColor(kRed);
    line->Draw("same");
    hff_reweights->Draw("same");
    
    if(scalebybin) canv1->SaveAs("plots/ffDiempt.pdf","pdf");
    
    
    //----------------------------------------------------------------------------
    
    
    cout << "Let's do ee now" << endl;
    
    TCanvas* canv2 = new TCanvas("diempt_ee","diempt_ee",50,50,W_ref,H_ref);
    canv2->cd();
    canv2->SetFillColor(0);
    canv2->SetBorderMode(0);
    canv2->SetFrameFillStyle(0);
    canv2->SetFrameBorderMode(0);
    canv2->SetLeftMargin( L/W_ref );
    canv2->SetRightMargin( R/W_ref );
    canv2->SetTopMargin( T/H_ref );
    canv2->SetBottomMargin( B/H_ref );
    canv2->SetTickx(0);
    canv2->SetTicky(0);
    canv2->SetGrid();
    gStyle->SetGridStyle(3);
    
    TPad* p2a = new TPad("p2a","",0.,.3,1.,1.);
    p2a->SetBottomMargin(0);
    p2a->Draw();
    TPad *p2b = new TPad("p2b","",0.,0.,1.,.3);
    p2b->SetTopMargin(0);
    p2b->SetBottomMargin(0.2);
    p2b->Draw();
    p2a->cd();
    p2a->SetGridx();
    p2a->SetGridy();
    gPad->SetLogy();
    gPad->Update();
    
    hggDiempt->Draw();
    gPad->Update();
    
    TH1F* heeDiempt = inputFile->FindObject("eeDiEMPt");
    heeDiempt->SetLineColor(kBlack);
    heeDiempt->SetStats(kTRUE);
    heeDiempt->Sumw2();
    TH1F * hee_reweights = new TH1F("ee_reweights"," ;DiEMPt(GeV);gg / ee", numBins,diemptbins);
    hee_reweights->SetStats(0);
    hee_reweights->Sumw2();
    
    
    
    if(toscale){
        float eescale =1.0 / heeDiempt->Integral();
        heeDiempt = scaleHist(heeDiempt, eescale);    }
    
    
    if(scalebybin){
        heeDiempt = scaleWidth(heeDiempt);
    }
    
    float fakemet, ggDiempt, wvalue, binvalue;
    float erroree, errorgg, error;
    for(int i =0; i <28; i++){
        ggDiempt = hggDiempt->GetBinContent(i);
        eeDiempt = heeDiempt->GetBinContent(i);
        errorgg = hggDiempt->GetBinError(i);
        erroree = heeDiempt->GetBinError(i);
        if(eeDiempt == 0 || ggDiempt ==0){
            wvalue = 0;
            error = 0;
        }
        else{
            wvalue = ggDiempt/eeDiempt;
            error = wvalue*sqrt((errorgg/ggDiempt)*(errorgg/ggDiempt)+(erroree/eeDiempt)*(erroree/eeDiempt));
        }
        if(i==5){
            cout<<"ggDiempt = " << ggDiempt << " pm " << errorgg<<endl;
            cout<<"eeDiempt = " << eeDiempt << " pm " << erroree << endl;
            cout <<"ratio = " << wvalue << "pm" << error << endl;
            
        }
        binvalue= hee_reweights->GetBinCenter(i);
        hee_reweights->Fill(binvalue, wvalue);
        hee_reweights->SetBinError(i, error);
    }
    
    
    
    heeDiempt->SetTitle("ee DiEMPt;;Events");
    heeDiempt->GetXaxis()->SetRangeUser(0,1000);
    //heeDiempt->GetYaxis()->SetRangeUser(-5,6000);
    heeDiempt->GetYaxis()->SetLabelSize(0.04);
    heeDiempt->Draw();
    
    gPad->Update();
    TPaveStats* s2 =(TPaveStats*) heeDiempt->GetListOfFunctions()->FindObject("stats");
    s2->SetLineColor(kBlack);
    s2->SetX2NDC(.9);
    s2->SetX1NDC(.7);
    s2->SetY1NDC(.36);
    s2->SetY2NDC(.51);
    
    TH1F* heeDiempt2 = heeDiempt->Clone();
    heeDiempt2->SetName("heeDiempt2");
    heeDiempt2->Sumw2(kFALSE);
    // heeDiempt2->SetFillColor(kBlack);
    heeDiempt2->SetLineColor(kBlack);
    //  heeDiempt2->SetFillStyle(3354);
    
    
    heeDiempt2->Draw("same");
    heeDiempt->Draw("same");
    hggDiempt->Draw("same");
    
    //Create a Legend for the plot
    TLegend* leg = new TLegend(0.55,0.7,0.9,0.9);
    
    leg->SetTextSize(.03);
    leg->AddEntry(hggDiempt, "gg","lf");
    leg->AddEntry(heeDiempt2, "ee", "lf");
    leg->Draw();
    gPad->Update();
    
    float r = p2a->GetRightMargin();
    float l = p2a->GetLeftMargin();
    float posX_ =   l + relPosX*(1-l-r);
    
    
    float t = p2a->GetTopMargin();
    float b = p2a->GetBottomMargin();
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
    
    
    p2b->cd();
    
    p2b->SetGridx();
    p2b->SetGridy();
    hee_reweights->GetXaxis()->SetRangeUser(0,1000);
    hee_reweights->SetLineColor(kBlack);
    hee_reweights->GetXaxis()->SetLabelSize(0.08);
    hee_reweights->GetYaxis()->SetLabelSize(0.08);
    hee_reweights->GetXaxis()->SetTitleSize(0.08);
    hee_reweights->GetYaxis()->SetTitleOffset(0.5);
    hee_reweights->GetYaxis()->CenterTitle();
    hee_reweights->GetYaxis()->SetTitleSize(0.08);
    hee_reweights->GetYaxis()->SetRangeUser(0,8.0);
    hee_reweights->Draw();
    
    TLine *line = new TLine(0,1,600,1);
    line->SetLineColor(kRed);
    line->Draw("same");
    hee_reweights->Draw("same");
    if(scalebybin) canv2->SaveAs("plots/eeDiempt.pdf","pdf");
    
    //----------------------------------------------------------------------------
    
    cout << "Let's do gf now" << endl;
    
    TCanvas* canv3 = new TCanvas("diempt_gf","diempt_gf",50,50,W_ref,H_ref);
    canv3->cd();
    canv3->SetFillColor(0);
    canv3->SetBorderMode(0);
    canv3->SetFrameFillStyle(0);
    canv3->SetFrameBorderMode(0);
    canv3->SetLeftMargin( L/W_ref );
    canv3->SetRightMargin( R/W_ref );
    canv3->SetTopMargin( T/H_ref );
    canv3->SetBottomMargin( B/H_ref );
    canv3->SetTickx(0);
    canv3->SetTicky(0);
    canv3->SetGrid();
    gStyle->SetGridStyle(3);
    
    TPad *p3a = new TPad("p3a","",0.,.3,1.,1.);
    p3a->SetBottomMargin(0);
    p3a->Draw();
    TPad *p3b = new TPad("p3b","",0.,0.,1.,.3);
    p3b->SetTopMargin(0);
    p3b->SetBottomMargin(0.2);
    p3b->Draw();
    p3a->cd();
    p3a->SetGridx();
    p3a->SetGridy();
    gPad->SetLogy();
    gPad->Update();
    
    hggDiempt->Draw();
    gPad->Update();
    
    TH1F* hgfDiempt = inputFile->FindObject("gfDiEMPt");
    hgfDiempt->SetLineColor(kBlack);
    hgfDiempt->SetStats(kTRUE);
    hgfDiempt->Sumw2();
    TH1F * hgf_reweights = new TH1F("gf_reweights"," ;DiEMPt(GeV);gg / gf", numBins,diemptbins);
    hgf_reweights->SetStats(0);
    hgf_reweights->Sumw2();
    
    
    
    if(toscale){
        float gfscale =1.0 / hgfDiempt->Integral();
        hgfDiempt = scaleHist(hgfDiempt, gfscale);
    }
    
    
    if(scalebybin){
        hgfDiempt = scaleWidth(hgfDiempt);
    }
    
    float fakemet, ggDiempt, wvalue, binvalue;
    float errorgf, errorgg, error;
    for(int i =0; i <28; i++){
        ggDiempt = hggDiempt->GetBinContent(i);
        gfDiempt = hgfDiempt->GetBinContent(i);
        errorgg = hggDiempt->GetBinError(i);
        errorgf = hgfDiempt->GetBinError(i);
        if(gfDiempt == 0 || ggDiempt ==0){
            wvalue = 0;
            error = 0;
        }
        else{
            wvalue = ggDiempt/gfDiempt;
            error = wvalue*sqrt((errorgg/ggDiempt)*(errorgg/ggDiempt)+(errorgf/gfDiempt)*(errorgf/gfDiempt));
        }
        binvalue= hgf_reweights->GetBinCenter(i);
        hgf_reweights->Fill(binvalue, wvalue);
        hgf_reweights->SetBinError(i, error);
    }
    
    
    
    hgfDiempt->SetTitle("gf DiEMPt;;Events");
    hgfDiempt->GetXaxis()->SetRangeUser(0,1000);
    //hgfDiempt->GetYaxis()->SetRangeUser(-5,6000);
    hgfDiempt->GetYaxis()->SetLabelSize(0.04);
    hgfDiempt->Draw();
    
    gPad->Update();
    TPaveStats* s3 =(TPaveStats*) hgfDiempt->GetListOfFunctions()->FindObject("stats");
    s3->SetLineColor(kBlack);
    s3->SetX2NDC(.9);
    s3->SetX1NDC(.7);
    s3->SetY1NDC(.36);
    s3->SetY2NDC(.51);
    
    TH1F* hgfDiempt2 = hgfDiempt->Clone();
    hgfDiempt2->SetName("hgfDiempt2");
    hgfDiempt2->Sumw2(kFALSE);
    // hgfDiempt2->SetFillColor(kBlack);
    hgfDiempt2->SetLineColor(kBlack);
    //  hgfDiempt2->SetFillStyle(3354);
    
    
    hgfDiempt2->Draw("same");
    hgfDiempt->Draw("same");
    hggDiempt->Draw("same");
    
    //Create a Legend for the plot
    TLegend* leg = new TLegend(0.55,0.7,0.9,0.9);
    
    leg->SetTextSize(.03);
    leg->AddEntry(hggDiempt, "gg","lf");
    leg->AddEntry(hgfDiempt2, "gf", "lf");
    leg->Draw();
    gPad->Update();
    
    float r = p3a->GetRightMargin();
    float l = p3a->GetLeftMargin();
    float posX_ =   l + relPosX*(1-l-r);
    
    
    float t = p3a->GetTopMargin();
    float b = p3a->GetBottomMargin();
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
    
    
    p3b->cd();
    
    p3b->SetGridx();
    p3b->SetGridy();
    hgf_reweights->GetXaxis()->SetRangeUser(0,1000);
    hgf_reweights->SetLineColor(kBlack);
    hgf_reweights->GetXaxis()->SetLabelSize(0.08);
    hgf_reweights->GetYaxis()->SetLabelSize(0.08);
    hgf_reweights->GetXaxis()->SetTitleSize(0.08);
    hgf_reweights->GetYaxis()->SetTitleOffset(0.5);
    hgf_reweights->GetYaxis()->CenterTitle();
    hgf_reweights->GetYaxis()->SetTitleSize(0.08);
    hgf_reweights->GetYaxis()->SetRangeUser(0,8.0);
    hgf_reweights->Draw();
    
    TLine *line = new TLine(0,1,600,1);
    line->SetLineColor(kRed);
    line->Draw("same");
    hgf_reweights->Draw("same");
    if(scalebybin) canv3->SaveAs("plots/gfDiempt.pdf","pdf");
    
    
    
    fileResults->Write();
    
    //close files
    //inputFile->Close();
    //  fileResults->Close();
    
}

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
