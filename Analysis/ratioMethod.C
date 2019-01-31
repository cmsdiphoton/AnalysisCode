#include "TMinuit.h"

TCanvas * SetUpCanvas(TString canvName);


void ratioMethod(){
    //----------- Define variables ---------------
    
    const int nLoop = 24; //Last bin number to process
    const int firstSignalBin  = 19;
    const int numSignalBins   =  6;
    float bins[numSignalBins+1] = {100,115,130,150,185,250,350};
    int nBin50 = 14; //Define top bin to include in normalization
    int nContour = 120; //Number of points on 1 sigma contours used to determine fit uncertainties.
    Double_t MetBins[nLoop+1]  = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,75,85,100,115,130,150,185,250,350};
    TH1F* h_highPred = new TH1F("h_highPred","Predicted MET for Mixed R9;#slash{E}_{T} (GeV);Events",nLoop,MetBins);

    bool twice = false;
    bool subtract = false;
    Double_t transfer = 0.0263;
    Double_t transfer_err = 0.30;
    
    //----------- Get histograms ---------------
    
    TFile *fileIn = TFile::Open("hist_data_testFF_R9.root");
    (TH1F*) fileIn->cd();
    
    TFile *fileEF = TFile::Open("hist_data_testFF_R9.root");
    (TH1F*) fileEF->cd();
    
    TFile *fileResults = TFile::Open("OutputFiles/finalBkg.root","RECREATE");
    fileResults->cd();
    
    TH1F * hHigh = (TH1F *) fileIn->Get("ffMet_MixR9");
    TH1F * hLow = (TH1F *) fileIn->Get("ffMet");
    
    TH1F * hgg = (TH1F *) fileIn->Get("ggMet");
    TH1F * hff = (TH1F *) fileIn->Get("ffMet");
    
    //for subtraction
    TH1F * hefLow = (TH1F *) fileEF->Get("efMet");
    TH1F * hefHigh = (TH1F *) fileEF->Get("efMet_HighR9");

    hefLow->Scale(transfer);
    hefHigh->Scale(transfer);

    if(subtract) {
     //   hHigh->Add(hefLow,-1);
        hgg->Add(hefLow,-1);
    }

    cout << hgg->GetNbinsX() << " " << h_highPred->GetNbinsX() << " " << hHigh->GetNbinsX() << endl;

    //------------ gg/ff ratio method ---------------
    //Setup canvas and fits
    TCanvas* canv_xcheck = SetUpCanvas("ratio");
    canv_xcheck->cd();
    
    TGraphAsymmErrors * hrat = new TGraphAsymmErrors(hHigh);
    hrat->Clear();
    hrat->Divide(hHigh,hLow,"pois");
    hrat->SetLineColor(kBlack);
    
    hrat->SetTitle("Ratio of High R9 to Low R9 ff events;#slash{E}_{T} (GeV); High R9 / Low R9");
    hrat->Draw("AP");
    
    gStyle->SetOptFit(111);
    
    TF1* expFunc = new TF1("expFunc","[0]*exp(-[1]*x)",3,350);
    expFunc->SetParameters(2.24,0);
    
    //Perform exponential fit
    TFitResultPtr results = hrat->Fit(expFunc,"RS EX0","",3,350);
    results->Print();
    gMinuit->SetErrorDef(1);
    TGraph *gr0 = (TGraph *)gMinuit->Contour(nContour,0,1);
    double * err0 = gr0->GetX();
    double * err1 = gr0->GetY();
    
   hrat->Draw("AP");
   // gPad->Update();
    
    TF1 * expFunc2 = hrat->GetFunction("expFunc");
    
    Double_t expEst[nLoop];
    Double_t expErrUp[nLoop];
    Double_t expErrLow[nLoop];
    Double_t x[nLoop], xel[nLoop], xeu[nLoop];
    Double_t expStatUp[nLoop];
    Double_t expStatLow[nLoop];
    Double_t expFitUp[nLoop];
    Double_t expFitLow[nLoop];
    
    Double_t ratioValue[nLoop];
    Double_t ratioErrUp[nLoop];
    Double_t ratioErrLow[nLoop];

    float p0,p1,e0,e1;
    p0 = expFunc2->GetParameter(0);
    p1 = expFunc2->GetParameter(1);
    e0 = expFunc2->GetParError(0);
    e1 = expFunc2->GetParError(1);
    
    cout << p0 << " " << p1 << endl;
    
    TF1* funcBlue = new TF1("func","[0]*exp(-[1]*x)",0,350);
    funcBlue->SetParameters(p0,p1);
    funcBlue->SetLineColor(kBlue);
    funcBlue->SetLineWidth(5);
    funcBlue->Draw("same");
    
    //Calculate exponential result and fit error.
    for(int i =1;i<=nLoop;i++){
        //        cout << i << " " << bins[i] << endl;
        float width      = hHigh->GetBinWidth(i);
        int int1 = hHigh->GetBinLowEdge(i);
        int int2 = int1 + width;
        
        if(!twice) expFunc2->SetParameters(p0,p1);
        else expFunc2->SetParameters(p0*p0,p1*2);
            
        float funcInt    = expFunc2->Integral(int1, int2);
        if(i == nLoop) { width = 1; funcInt = expFunc2->Eval(int1); }
        float low         = hff->GetBinContent(i);
       // cout << low << " " << hLow->GetBinContent(i) << " ";
        float prediction = funcInt * low / width;
        float max = 0;
        float min = 1000000;
        cout << i << " " << prediction << " " << width << " " << funcInt << " " << low << " " << expFunc2->Eval(int1) << " " << int1 << endl;
        
        if(i < nLoop){
            for(int j = 0; j < nContour; j++){
                expFunc2->SetParameters(err0[j],err1[j]);
                if(!twice) expFunc2->SetParameters(err0[j],err1[j]);
                else expFunc2->SetParameters(err0[j]*err0[j],err1[j]*2);
                float varied = expFunc2->Integral(int1, int2) * low / width;
                if(varied < 0) continue;
                max = TMath::Max(max,varied);
                min = TMath::Min(min,varied);
            }
        }
        else{
            cout << "Second case" << endl;
            for(int j = 0; j < nContour; j++){
                if(!twice) expFunc2->SetParameters(err0[j],err1[j]);
                else expFunc2->SetParameters(err0[j]*err0[j],err1[j]*2);
                float varied = expFunc2->Eval(int1)*low;
                max = TMath::Max(max,varied);
                min = TMath::Min(min,varied);
            }
        }
        
        if(min > prediction) cout <<" No way!!! " << min << " " << prediction << " " << int1 << " :( " << endl;
           
        expFunc2->SetParameters(p0,p1);
        
        expEst[i-1] = prediction;
        expFitUp[i-1]  = (max - prediction);
        expFitLow[i-1] = (prediction - min);
        
        x[i-1]   = hLow->GetBinCenter(i);
        xel[i-1] = width/2;
        xeu[i-1] = width/2;
        if(i == nLoop) {
            xel[i-1] = 50;
            xeu[i-1] = 50;
        }
        
        expStatUp[i-1] = hff->GetBinErrorUp(i) / low *prediction;
        expStatLow[i-1] = hff->GetBinErrorLow(i) / low *prediction;

        expErrUp[i-1] = TMath::Sqrt(expStatUp[i-1]*expStatUp[i-1]
                                    + expFitUp[i-1]*expFitUp[i-1]);
        expErrLow[i-1] = TMath::Sqrt(expStatLow[i-1]*expStatLow[i-1]
                                     + expFitLow[i-1]*expFitLow[i-1]);
        
        //cout << prediction << " " << expErrUp[i-1] <<endl;
        h_highPred->Fill(x[i-1], prediction);
        h_highPred->SetBinError(i,TMath::Abs(expErrUp[i-1]+expErrLow[i-1])/2 );
        
       // cout << int1 << " " << prediction << " " << expStatUp[i-1] << " " << expStatLow[i-1]  << " " << hHigh->GetBinContent(i) << " " << hHigh->GetBinErrorUp(i) << " " << hHigh->GetBinErrorLow(i) << endl;
        //cout << funcInt << " " << prediction << " " << pp << " " << pm << " " << mm << " " << mp << " " << max << " " << min << endl;
        //cout << "& ${" << prediction << "}^{+"<< (max-prediction) << "}_{-" << (prediction -min) << "}$" << endl;
        
        //cout << int1 << " " <<expFitUp[i-1]  << " " << expFitLow[i-1]  << " "  << prediction << endl ;
        
        float gg = hgg->GetBinContent(i);
        ratioValue[i-1] = gg/prediction;
       // cout << gg << " " << prediction << " " << ratioValue[i-1] << " ";
       // cout << expErrUp[i-1] << " " << hgg->GetBinErrorUp(i) << " ";
        
       ratioErrUp[i-1] = sqrt( (expErrUp[i-1]  * expErrUp[i-1]) / (prediction * prediction) +
                                (hgg->GetBinErrorUp(i) * hgg->GetBinErrorUp(i)) / ( gg*gg));
        ratioErrLow[i-1] = sqrt( (expErrLow[i-1]  * expErrLow[i-1]) / (prediction * prediction) +
                               (hgg->GetBinErrorLow(i) * hgg->GetBinErrorLow(i)) / ( gg*gg));
        
       // cout <<  ratioErrLow[i-1] << endl;
    }
    
    TCanvas* canv_high_pred= SetUpCanvas("high_pred");
    canv_high_pred->cd();
    //h_highPred->Draw();
    //hgg->Draw("same");
    
    TGraphAsymmErrors * g_highPred  = new TGraphAsymmErrors(nLoop,x,expEst,xel,xeu,expErrLow,expErrUp);
    
    g_highPred->SetName("g_highPred");
    
    TGraphAsymmErrors * g_ratio  = new TGraphAsymmErrors(nLoop,x,ratioValue,xel,xeu,ratioErrLow,ratioErrUp);
    g_ratio->SetName("g_ratio");
    g_ratio->Draw("AP");
    
   // highPred->Draw("AP");
    
    /*
    TGraphAsymmErrors * hrat2 = new TGraphAsymmErrors(h_highPred);
    hrat2->Clear();
    hrat2->Divide(ggMet,h_highPred,"pois");
    hrat2->SetLineColor(kBlack);
   hrat2->Draw("AP");
     hrat2->SetTitle("#gamma#gamma / ff prediction from High R9 ff sample;MET (Gev); #gamma#gamma / ff prediction ");
     c1->SaveAs("ggffPredFromHighFF.pdf")
  */
    g_highPred->Write("g_highPred");
    g_ratio->Write("g_ratio");
    h_highPred->Write();
    hgg->Write();
   fileResults->Close();
    
}

TCanvas * SetUpCanvas(TString canvName){
    int H_ref = 850;
    int W_ref = 900;
    
    // references for T, B, L, R
    float T = 0.08*H_ref;
    float B = 0.12*H_ref;
    float L = 0.12*W_ref;
    float R = 0.04*W_ref;
    
    
    TCanvas* ctemp = new TCanvas(canvName,canvName,50,50,W_ref,H_ref);
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
    ctemp->SetGrid();
    gStyle->SetGridStyle(3);
    
    
    return ctemp;
    
}

