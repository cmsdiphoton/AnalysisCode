#include "TMinuit.h"

TCanvas * SetUpCanvas(TString canvName);


void FinalAnalysis_CrossCheckCentral(){
    
    bool useConstant = false;
    
    //----------------Draw histograms and open files--------------
    TFile *fileIn = TFile::Open("hist_data_NewSkims_rw.root");
    (TH1F*) fileIn->cd();
    
    TFile *fileZGG = TFile::Open("hist_ZGG_THESIS.root");
    
    TFile *fileResults = TFile::Open("OutputFiles/finalBkg.root","RECREATE");
    fileResults->cd();
    
    //----------------Get histograms locally--------------
    TH1F * hegMet = (TH1F*) fileIn->Get("egMet");
    TH1F * hggMet = (TH1F*) fileIn->Get("ggMet");
    TH1F * hffMet = (TH1F*) fileIn->Get("ffMet");
    
    TH1F * hZGGMet = (TH1F*) fileZGG->Get("ggMet");
    hZGGMet->SetName("hZGGMet");
    TH1F * hZGGMet_Unweighted= (TH1F*) fileZGG->Get("ggMetUnweighted");
    hZGGMet_Unweighted->SetName("hZGGMetUnweighted");
    
    TH1F * finalEWKMet          = (TH1F*) hegMet->Clone("finalEWKMet");
    
    //----------------Open files with final numbers----------------

    ofstream zggFile;
    zggFile.open("OutputFiles/ZGG.txt");
    
    ofstream xsecFile;
    xsecFile.open("OutputFiles/CrossCheckUncert.txt");
    
    ofstream sysFile;
    sysFile.open("OutputFiles/QCDSystematicUncertainties.txt");
    
    ofstream tabEWKFile;
    tabEWKFile.open("OutputFiles/EWKTable.txt");
    
    ofstream tabQCDFile;
    tabQCDFile.open("OutputFiles/QCDTable.txt");
    
    ofstream tabAllFile;
    tabAllFile.open("OutputFiles/FinalSystematicsTable.txt");
    
    ofstream EWKFile;
    EWKFile.open("OutputFiles/EWKSystematics.txt");
    
    ofstream tabPredictFile;
    tabPredictFile.open("OutputFiles/FinalPredictionTable.txt");
    
    ofstream tabCompareCrossCheckFile;
    tabCompareCrossCheckFile.open("OutputFiles/CompareCrossCheckTable.txt");
    

    
    const int nLoop = 25; //Last bin number to process
    const int firstSignalBin  = 19;
    const int numSignalBins   =  6;
    float bins[numSignalBins+1] = {100,115,130,150,185,250,350};
    int nBin50 = 14; //Define top bin to include in normalization
    int nContour = 120; //Number of points on 1 sigma contours used to determine fit uncertainties.
    Double_t MetBins[nLoop]  = {0,3,6,9,12,15,18,22,26,29,32,36,40,45,50,60,75,85,100,115,130,150,185,250,350};

    TH1F* finalQCDMet = new TH1F("finalQCDMet","QCD Estimate;MET;Events/GeV",nLoop-1,MetBins);
    
    TCanvas* c2 = SetUpCanvas("c2");
    c2->cd();
    
    //----------------Declare arrays--------------
    Double_t QCDEst[nLoop - 1];
    Double_t QCDErrUp[nLoop - 1];
    Double_t QCDErrLow[nLoop - 1];
    
    Double_t QCDEst_WScaled[nLoop - 1];
    Double_t QCDErrUp_WScaled[nLoop - 1];
    Double_t QCDErrLow_WScaled[nLoop - 1];
    
    Double_t QCDStatErrUp[nLoop - 1];
    Double_t QCDStatErrLow[nLoop - 1];
    
    Double_t QCDFitErrUp[nLoop - 1];
    Double_t QCDFitErrLow[nLoop - 1];
    
    Double_t QCDCrossCheckErr[nLoop - 1];
    
    Double_t EWKEst_WScaled[nLoop - 1];
    Double_t EWKErrUp_WScaled[nLoop - 1];
    Double_t EWKErrLow_WScaled[nLoop - 1];
    
    Double_t EWKStatErrUp[nLoop - 1];
    Double_t EWKStatErrLow[nLoop - 1];
    
    Double_t x[nLoop-1],
    xel[nLoop-1],
    xeu[nLoop-1],
    TotalEst_WScaled[nLoop-1],
    TotalErrUp_WScaled[nLoop-1],
    TotalErrLow_WScaled[nLoop-1];
    
    //-----------------------------------------


    //----------------Calculate ZGG--------------
    
    //First scale to get rid of the effects of PU reweighting
    float prePU = hZGGMet_Unweighted->Integral();
    float postPU = hZGGMet->Integral();
    float PUscale = prePU/postPU;
    hZGGMet->Scale(PUscale);
    cout << "PU scale for ZGG = " << PUscale;
    
    //Then scale by cross section
    float lumi   = 35.9218;
    float ZGGxsec= 0.07477*1000;
    float numberZGGevents = 994078;
    float sf_ZGG = lumi * ZGGxsec / numberZGGevents;
    hZGGMet->Scale(sf_ZGG);
    //-----------------------------------------
    
    //----- Values from cross check---------
    Double_t QCDCrossCheckValue[nLoop-1] = {
        1712.23,
        4631.7,
        7500.02,
        9445.86,
        10630.8,
        11258.7,
        15491.2,
        14774.8,
        10116.9,
        8991.92,
        10214.6,
        8314.94,
        7950.52,
        5644.59,
        6437.72,
        3461.11,
        675.998,
        341.839,
        93.1006,
        30.3169,
        17.1326,
        8.86404,
        2.69912,
        0.826031
    };

    Double_t QCDCrossCheckAverageError[nLoop-1] ={
        38.4371,
        67.4875,
        89.5778,
        101.541,
        106.92,
        108.363,
        131.037,
        124.147,
        95.8334,
        88.6789,
        96.5459,
        85.4901,
        84.7162,
        69.4887,
        79.8708,
        56.3431,
        20.5885,
        13.9426,
        6.78962,
        3.74334,
        2.72362,
        1.86102,
        0.949961,
        0.526043
    };
    
    //------------ Fit gg/ff ratio ---------------
    //Setup canvas and fits
    TCanvas* canv_xcheck = SetUpCanvas("ratio");
    canv_xcheck->cd();
    
    TGraphAsymmErrors * hrat = new TGraphAsymmErrors(hggMet);
    hrat->Clear();
    hrat->Divide(hggMet,hffMet,"pois");
    hrat->SetLineColor(kBlack);
    
    hrat->SetTitle("Ratio of number of #gamma#gamma to ff events;#slash{E}_{T} (GeV); #gamma#gamma/ff");
    hrat->Draw("AP");
    
    TF1* expFunc = new TF1("expFunc","[0]*exp(-[1]*x)",3,100);
    expFunc->SetParameters(2.24,0);
    
    //Perform exponential fit
    TFitResultPtr results = hrat->Fit(expFunc,"RS EX0","",3,100);
    results->Print();
    gMinuit->SetErrorDef(1);
    TGraph *gr0 = (TGraph *)gMinuit->Contour(nContour,0,1);
    double * err0 = gr0->GetX();
    double * err1 = gr0->GetY();
    
    hrat->Draw("AP");
    gPad->Update();
    
    TF1 * expFunc2 = hrat->GetFunction("expFunc");

    float p0,p1,e0,e1;
    p0 = expFunc2->GetParameter(0);
    p1 = expFunc2->GetParameter(1);
    e0 = expFunc2->GetParError(0);
    e1 = expFunc2->GetParError(1);
    
    cout << p0 << " " << p1 << endl;
    
    expFunc->SetParameters(p0,p1);
    
    
    //------------ Scale eg sample by fake rate ---------------
    Double_t transfer = 0.0263;
    Double_t transfer_err = 0.30; //percent error on transfer factor
    finalEWKMet->Scale(transfer);
    
    //------------ File headers ---------------

    sysFile << "ffValue FFStatLow FFStatUp DiffUncertLow DiffUncertUp" << endl;
    tabAllFile<<"\ETmiss Bin (GeV) & QCD Estimate &  Stat. Uncert. & Shape Uncert. & Reweighting Uncert. \n";
    EWKFile<<"METBin EWK_Estimate StatLow_Events StatUp_Events FakeRate_Events StatLow_Percent StatUp_Percent FakeRate_Percent \n" ;
    
    //------------ Loop variables ---------------
    float width, funcInt, ffuw, ffvalue, start, endpt;
    
    
    //------------ Calculate final errors ---------------
    for(int i = 1; i < nLoop;i++){
        int j = i-1;
        width = hffMet->GetXaxis()->GetBinWidth(i);
        start = hffMet->GetBinLowEdge(i);
        endpt = hffMet->GetBinLowEdge(i) + width;
        
        funcInt = expFunc->Integral(start, endpt);
        ffuw    = hffMet->GetBinContent(i);
        QCDEst[j] = funcInt * ffuw / width;
        if(start == 250) QCDEst[j] = ffuw*expFunc->Eval(start);
        
        //------------------Calculate statistical errors---------------
        QCDStatErrUp[j] = hffMet->GetBinErrorUp(i) / ffuw * QCDEst[j];
        QCDStatErrLow[j]  = hffMet->GetBinErrorLow(i)  / ffuw * QCDEst[j];
        
        if(i >= firstSignalBin){
            tabAllFile <<"$" << start << "-" << endpt << "$ & " << QCDEst[j] <<" & +"<<QCDStatErrUp[j] << ", -" << QCDStatErrLow[j] << " & pm";
        }
        sysFile << QCDEst[j] << " " << QCDStatErrUp[j] << " " << QCDStatErrLow[j] << " ";
        
        //------------------Calculate fit uncertainties ---------------
        float max = 0;
        float min = 100000;
        for(int k = 0; k < nContour; k++){
            expFunc2->SetParameters(err0[k],err1[k]);
            float varied = expFunc2->Integral(start, endpt) *ffuw / width;
            if(start ==250) varied = expFunc2->Eval(start)*ffuw;
            max = TMath::Max(max,varied);
            min = TMath::Min(min,varied);
        }
        
        QCDFitErrUp[j]  = (max - QCDEst[j]);
        QCDFitErrLow[j] = (QCDEst[j] - min);
        
        //------------------Calculate uncertainties from cross check ---------------
        float statPlusFitUp  = TMath::Sqrt( QCDStatErrUp[j]*QCDStatErrUp[j]   + QCDFitErrUp[j]*QCDFitErrUp[j]  );
        float statPlusFitLow = TMath::Sqrt( QCDStatErrLow[j]*QCDStatErrLow[j] + QCDFitErrLow[j]*QCDFitErrLow[j] );
        float statPlusFitAve = (statPlusFitUp + statPlusFitLow)/2;
        
        float diff = TMath::Abs(QCDCrossCheckValue[j] - QCDEst[j]);
        float diffErr = sqrt( (QCDCrossCheckAverageError[j] * QCDCrossCheckAverageError[j]) +
                             (statPlusFitAve * statPlusFitAve) );
        
        QCDCrossCheckErr[j] = TMath::Max(diff, diffErr);
        
        sysFile << QCDFitErrUp[j] << " " << QCDFitErrLow[j] << " " << QCDCrossCheckErr[j] ;

        
        //------------------ Calculate final QCD values ---------------
        QCDErrUp[j]   = TMath::Sqrt( QCDStatErrUp[j]*QCDStatErrUp[j]   + QCDFitErrUp[j]*QCDFitErrUp[j]   + QCDCrossCheckErr[j]*QCDCrossCheckErr[j]);
        QCDErrLow[j]  = TMath::Sqrt( QCDStatErrLow[j]*QCDStatErrLow[j] + QCDFitErrLow[j]*QCDFitErrLow[j] + QCDCrossCheckErr[j]*QCDCrossCheckErr[j]);

        QCDEst_WScaled[j]    = QCDEst[j]   / width;
        QCDErrUp_WScaled[j]  = QCDErrUp[j] / width;
        QCDErrLow_WScaled[j] = QCDErrLow[j]/ width;

        x[j]   = hffMet->GetBinCenter(i);
        xel[j] = width/2;
        xeu[j] = width/2;
        
        sysFile << " " << QCDErrUp[j] << " " << QCDErrLow[j] << " \n";
        
        //------------------Calculate EGamma estimate-------------------
        EWKEst_WScaled[j]    = hegMet->GetBinContent(i)*transfer / width;
        EWKErrUp_WScaled[j] = transfer*TMath::Sqrt(hegMet->GetBinErrorUp(i)  *  hegMet->GetBinErrorUp(i)  + transfer_err*transfer_err*hegMet->GetBinContent(i)*hegMet->GetBinContent(i)) / width;
        EWKErrLow_WScaled[j] = transfer*TMath::Sqrt(hegMet->GetBinErrorLow(i) *  hegMet->GetBinErrorLow(i) + transfer_err*transfer_err*hegMet->GetBinContent(i)*hegMet->GetBinContent(i)) / width;
        
        float ewk_value = hegMet->GetBinContent(i)*transfer;
        if(i >= firstSignalBin && i < firstSignalBin+numSignalBins) {
            tabEWKFile<< "$" << start << "-" << endpt << "$ & " << ewk_value << "pm " << EWKErrUp_WScaled[i-1]*width << "\\" << "\n" << "hline" << "\n";

            EWKFile << start << "-" << endpt << " " << ewk_value << " ";
            EWKFile << hegMet->GetBinErrorLow(i)/hegMet->GetBinContent(i)* ewk_value << " ";
            EWKFile << hegMet->GetBinErrorUp(i)/hegMet->GetBinContent(i)* ewk_value << " ";
            EWKFile << transfer_err* ewk_value << " ";
            EWKFile << ( 1.0 - hegMet->GetBinErrorLow(i)/hegMet->GetBinContent(i) ) << " ";
            EWKFile << ( 1.0 + hegMet->GetBinErrorUp(i)/hegMet->GetBinContent(i) )  << " 1.3 XXXXX";
        }
        
        //---------------Calculate ZGG contribution---------------------
        double zggValue_WScaled = (hZGGMet->GetBinContent(i)) /width;
        double zggError_WScaled = zggValue_WScaled / 2;
        if(i >= firstSignalBin) {
            zggFile << "$" << start << "-" << endpt << "$ & $" << zggValue_WScaled*width << " pm " << zggError_WScaled*width << "$ \\" << endl;
        }
        //---------------Combine EWK and QCD and ZGG Errors!---------------------
        TotalErrUp_WScaled[j] = TMath::Sqrt(EWKErrUp_WScaled[j] * EWKErrUp_WScaled[j] +
                                            QCDErrUp_WScaled[j]* QCDErrUp_WScaled[j] +
                                            zggError_WScaled*zggError_WScaled);
        TotalErrLow_WScaled[j] = TMath::Sqrt(EWKErrLow_WScaled[j] * EWKErrLow_WScaled[j] +
                                             QCDErrLow_WScaled[j]* QCDErrLow_WScaled[j] +
                                             zggError_WScaled*zggError_WScaled);
        double test2 = EWKEst_WScaled[j] + QCDEst_WScaled[j];
        double test3 = test2 + zggValue_WScaled;
        TotalEst_WScaled[j] = test3;
        if(i >= firstSignalBin) {
            tabPredictFile<< "$"<< start << "-"<<endpt << "$ & ${ " << TotalEst_WScaled[j]*width << " }^{+ ";
            tabPredictFile<< TotalErrUp_WScaled[j]*width << " }_{- " << TotalErrLow_WScaled[j]*width << " }$ & " << hggMet->GetBinContent(i)<< "\\" << "\n";
        }
        
        //hggMet : errors are statistical only
        float widthError = (hggMet->GetBinErrorUp(i) + hggMet->GetBinErrorLow(i))/(2*width);
        float widthValue= hggMet->GetBinContent(i) / width;
        hggMet->SetBinContent(i,widthValue);
        hggMet->SetBinError(i,widthError);
        
        finalQCDMet->SetBinContent(i,QCDEst_WScaled[j]);
        finalQCDMet->SetBinError(i,QCDErrUp_WScaled[j]);
        
    } //end of loop over bins

    //------------ Make datacard ---------------
    /*ofstream t5wgDatacard;
    t5wgDatacard.open("OutputFiles/datacard_T5Wg.txt");

    t5wgDatacard << "imax 6  number of channels \n \n";
    t5wgDatacard << "jmax 3  number of backgrounds \n \n";
    t5wgDatacard << "kmax *  number of nuisance parameters \n\n";
    t5wgDatacard << "\n------------ \n \n";
    t5wgDatacard << "bin            1      2      3      4      5      6\n";
    t5wgDatacard << "observation    ";

    for(int i=firstSignalBin;i<(firstSignalBin+numSignalBins);i++){
        float unblinded = hggMet->GetBinContent(i)*hggMet->GetBinWidth(i);
        t5wgDatacard << unblinded << "  ";
    }
    t5wgDatacard << "\n------------ \n \n";
    t5wgDatacard << "bin            1      1      1      1      2      2      2      2      3      3      3      3      4      4      4      4      5      5      5      5      6      6      6      6\n\n";
    t5wgDatacard << "process        t5gg   qcd    ewk    zgg    t5gg   qcd    ewk    zgg    t5gg   qcd    ewk    zgg    t5gg   qcd    ewk    zgg    t5gg   qcd    ewk    zgg    t5gg   qcd    ewk    zgg\n\n";
    t5wgDatacard << "process        0      1      2      3      0      1      2      3      0      1      2      3      0      1      2      3      0      1      2      3      0      1      2      3 \n\n";
    t5wgDatacard << "rate           ";
    for(int i=firstSignalBin;i<(firstSignalBin+numSignalBins);i++){
        float width = hffMet->GetBinWidth(i);
        t5wgDatacard << "T5GG" << (i+1-firstSignalBin) << "  " << (y[i-1]*width) << "  " << (ewk[i-1]*width) << "  " << hZGGMet->GetBinContent(i) << "  ";
        
    }
    t5wgDatacard << "\n------------ \n \n";
    //-----------MC stats-----------
    for(int i=1;i<numSignalBins+1;i++){
        t5wgDatacard << "mcStats_"<< (i-1) << "  gmN  MCN"<< (i) << "    ";
        for(int j = 1; j <= i-1;j++){
            t5wgDatacard << "-      -      -      -      ";
        }
        t5wgDatacard << "MCA" << i << "   -      -      -      ";
        for(int j = 1; j <= numSignalBins-i;j++){
            t5wgDatacard << "-      -      -      -      ";
        }
        t5wgDatacard << "\n \n";
    }
    //------------------------------
    //------------QCD stats------------------
    for(int i=1;i<numSignalBins+1;i++){
        int bin = firstSignalBin+i-1;
        t5wgDatacard << "qcdStats_"<< (i-1) << "  gmN   " << hffMet->GetBinContent(bin) << "    ";
        float width = hffMet->GetBinWidth(bin);
        float weight = y[bin-1]*width/hffMet->GetBinContent(bin);
        for(int j = 1; j <= i-1;j++){
            t5wgDatacard << "-      -      -      -      ";
        }
        t5wgDatacard << "-      " << weight<< "   -      -      ";
        for(int j = 1; j <= numSignalBins-i;j++){
            t5wgDatacard << "-      -      -      -      ";
        }
        t5wgDatacard << "\n \n";
    }
    //------------------------------
    //------------EWK stats------------------
    for(int i=1;i<numSignalBins+1;i++){
        int bin = firstSignalBin+i-1;
        t5wgDatacard << "ewkStats_"<< (i-1) << "  gmN   " << hegMet->GetBinContent(bin) << "    ";
        for(int j = 1; j <= i-1;j++){
            t5wgDatacard << "-      -      -      -      ";
        }
        t5wgDatacard << "-      -      " <<transfer << "   -      ";
        for(int j = 1; j <= numSignalBins-i;j++){
            t5wgDatacard << "-      -      -      -      ";
        }
        t5wgDatacard << "\n \n";
    }
    //------------------------------
    //-----------Other uncertainties-----------
    t5wgDatacard << "lumi     lnN   1.023  -      -      -      1.023  -      -      -      1.023  -      -      -      1.023  -      -      -      1.023  -      -      -      1.023  -      -      - \n \n";
    t5wgDatacard << "phoSf    lnN   1.025  -      -      -      1.025  -      -      -      1.025  -      -      -      1.025  -      -      -      1.025  -      -      -      1.025  -      -      -\n \n";
    t5wgDatacard << "jes      lnN   ";
    for(int i=1;i<numSignalBins+1;i++){
        t5wgDatacard << "JES" << i << "   -      -      -      ";
    }
    //------------------------------
    //-----------QCD uncertainties-----------

    t5wgDatacard<<"\n\n";
    t5wgDatacard << "diempt   lnN   ";
    for(int i=firstSignalBin;i<(firstSignalBin+numSignalBins);i++){
        t5wgDatacard << "-      " << 1+diempt_err[i-1] << "  -      -      ";
    }
    
    t5wgDatacard<<"\n\n";
    t5wgDatacard << "qcdShape lnN   ";
    for(int i=0;i<numSignalBins;i++){
        t5wgDatacard << "-      " << percentShape[i] << "  -      -      ";
    }
    //------------------------------
    //------------EWK uncertainty------------------
    t5wgDatacard<<"\n\n";
    t5wgDatacard << "ewkFake  lnN   -      -      1.300  -      -      -      1.300  -      -      -      1.300  -      -      -      1.300  -      -      -      1.300  -      -      -      1.300  -";
    //------------------------------
    //------------ZGG uncertainty------------------
    t5wgDatacard<<"\n\n";
    t5wgDatacard << "zggErr   lnN   -      -      -      1.500  -      -      -      1.500  -      -      -      1.500  -      -      -      1.500  -      -      -      1.500   -      -      -      1.500 \n \n";
    */
    //------------ Save histograms to output root file ---------------

    TGraphAsymmErrors * finalBkg  = new TGraphAsymmErrors(nLoop-1,x,TotalEst_WScaled,xel,xeu,TotalErrLow_WScaled,TotalErrUp_WScaled);
    finalBkg->SetName("finalBkg");
    finalBkg->Write();
    
    TGraphAsymmErrors * finalEWK  = new TGraphAsymmErrors(nLoop-1,x,EWKEst_WScaled,xel,xeu,EWKErrLow_WScaled,EWKErrUp_WScaled);
    finalEWK->SetName("finalEWK");
    finalEWK->Write();
    
    TGraphAsymmErrors * finalQCD  = new TGraphAsymmErrors(nLoop-1,x,QCDEst_WScaled,xel,xeu,QCDErrLow_WScaled,QCDErrUp_WScaled);
    finalQCD->SetName("finalQCD");
    finalQCD->Draw("APZ");
    finalQCD->Write();

    finalQCDMet->Write();
    
    finalEWKMet->Scale(1,"width");
    finalEWKMet->Write();
    
    hggMet->Write();
    
    hZGGMet->Scale(1,"width");
    hZGGMet->Write();
    
   // gr_ee->Write();
    //gr_ff->Write();
    
   /* TCanvas* canv1 = SetUpCanvas("unweighted");
    DrawPlot(hggMet, finaleeMetUnweighted, finalffMetUnweighted, canv1);
    canv1->SaveAs("plots/unweighted.pdf","pdf");
    
    TCanvas* canv2 = SetUpCanvas("weighted");
    DrawPlot(hggMet, finaleeMet, finalffMetReweighted, canv2);
    if(useFF_diempt || useDiempt) canv2->SaveAs("plots/reweighted.pdf","pdf");
    
    TCanvas* canv3 = SetUpCanvas("2Dweighted");
    DrawPlot(hggMet, finaleeMet_NJet, finalffMetReweighted, canv3);
    if(useFF_2D || use2D) canv3->SaveAs("plots/2Dweighted.pdf","pdf");
    
    
    TCanvas* canv4 = SetUpCanvas("eemets");
    DrawPlot(finaleeMet, finaleeMet_NJet, canv4);
    canv4->SaveAs("plots/eeComparison.pdf","pdf");*/
    
}
//******************************************

//******************************************

//******************************************

//******************************************

//******************************************

//******************************************

//******************************************


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


void DrawPlot(TH1F* h1, TH1F* h2, TH1F* h3, TCanvas* ctemp){
    
    ctemp->cd();
    TPad *p1a = new TPad("p1a","",0.,.3,1.,1.);
    p1a->SetBottomMargin(0);
    p1a->Draw();
    TPad *p1b = new TPad("p1b","",0.,0.,1.,.3);
    p1b->SetTopMargin(0);
    p1b->SetBottomMargin(0.2);
    p1b->Draw();
    p1a->cd();
    p1a->SetGridx();
    p1a->SetLogy();
    //  p1a->SetGridy();
    gPad->Update();
    
    TString cmsText     = "CMS";
    float cmsTextFont   = 61;  // default is helvetic-bold
    float cmsTextSize   = 0.5;
    
    TString extraText   = "Preliminary";
    float extraTextFont = 52;  // default is helvetica-italics
    float extraOverCmsTextSize  = 0.7;
    float extraTextSize = extraOverCmsTextSize*cmsTextSize;
    
    TString lumiText = "35.9 fb^{-1}";
    
    float relPosX      = 0.045;
    float relPosY      = -0.01;
    
    float r = p1a->GetRightMargin();
    float l = p1a->GetLeftMargin();
    float posX_ =   l + relPosX*(1-l-r);
    
    float t = p1a->GetTopMargin();
    float b = p1a->GetBottomMargin();
    float posY_ = 1-t - relPosY*(1-t-b);
    
    int   align_       = 10*1+1;
    float relExtraDX   = 1.5;
    
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
    
    h1->SetLineColor(kRed);
    h1->SetMarkerStyle(kFullTriangleUp);
    h1->SetMarkerSize(0.6);
    h1->SetMarkerColor(kRed);
    h1->SetStats(kFALSE);
    
    h2->SetLineColor(kBlack);
    h2->SetMarkerStyle(kFullSquare);
    h2->SetMarkerSize(0.6);
    h2->SetMarkerColor(kBlack);
    h2->SetStats(kFALSE);
    
    h3->SetLineColor(kBlue);
    h3->SetMarkerStyle(kFullCircle);
    h3->SetMarkerSize(0.6);
    h3->SetMarkerColor(kBlue);
    h3->SetStats(kFALSE);
    
    h3->SetTitle("Missing Transverse Energy");
    
    h3->Draw("e");
    h2->Draw("e same");
    h1->Draw("e same");
    
    //Create a Legend for the plot
    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    
    leg->SetTextSize(.05);
    leg->AddEntry(h1, "gg","lp");
    leg->AddEntry(h2, "ee", "lp");
    leg->AddEntry(h3, "ff","lp");
    leg->Draw();
    gPad->Update();
    
    
    p1b->cd();
    
    TH1F * h2_ratio = (TH1F*) h2->Clone("h2_ratio");
    h2_ratio->Reset();
    h2_ratio->SetTitle(" ");
    h2_ratio->SetStats(0);
    h2_ratio->Sumw2();
    
    TH1F * h3_ratio = (TH1F*) h3->Clone("h3_ratio");
    h3_ratio->Reset();
    h3_ratio->SetTitle(" ");
    h3_ratio->SetStats(0);
    h3_ratio->Sumw2();
    
    float ffMet, ggMet, eeMet, wvalue, binvalue;
    float erroree, errorgg, errorff, error;
    for(int i =0; i < h1->GetNbinsX()+1; i++){
        ggMet = h1->GetBinContent(i);
        eeMet = h2->GetBinContent(i);
        errorgg = h1->GetBinError(i);
        erroree = h2->GetBinError(i);
        if(eeMet == 0 || ggMet ==0){
            wvalue = 0;
            error = 0;
        }
        else{
            wvalue = ggMet/eeMet;
            error = wvalue*sqrt((errorgg/ggMet)*(errorgg/ggMet)+(erroree/eeMet)*(erroree/eeMet));
        }
        binvalue= h2_ratio->GetBinCenter(i);
        h2_ratio->Fill(binvalue, wvalue);
        h2_ratio->SetBinError(i, error);
        
        ffMet = h3->GetBinContent(i);
        errorff = h3->GetBinError(i);
        if(ffMet == 0 || ggMet ==0){
            wvalue = 0;
            error = 0;
        }
        else{
            wvalue = ggMet/ffMet;
            error = wvalue*sqrt((errorgg/ggMet)*(errorgg/ggMet)+(errorff/ffMet)*(errorff/ffMet));
        }
        binvalue= h3_ratio->GetBinCenter(i);
        h3_ratio->Fill(binvalue, wvalue);
        h3_ratio->SetBinError(i, error);
        
    }
    
    
    p1b->SetGridx();
    p1b->SetGridy();
    h2_ratio->GetXaxis()->SetRangeUser(0,350);
    h2_ratio->SetLineColor(kBlack);
    h2_ratio->GetXaxis()->SetLabelSize(0.08);
    h2_ratio->GetYaxis()->SetLabelSize(0.08);
    h2_ratio->GetXaxis()->SetTitleSize(0.08);
    h2_ratio->GetYaxis()->SetTitleOffset(0.5);
    h2_ratio->GetYaxis()->CenterTitle();
    h2_ratio->GetYaxis()->SetTitleSize(0.08);
    h2_ratio->GetYaxis()->SetRangeUser(0.5,2.5);
    h2_ratio->SetMarkerStyle(kFullSquare);
    h2_ratio->SetMarkerColor(kBlack);
    h2_ratio->SetMarkerSize(0.4);
    
    h2_ratio->Draw();
    
    h3_ratio->SetMarkerStyle(kFullCircle);
    h3_ratio->SetMarkerColor(kBlue);
    h3_ratio->SetMarkerSize(0.4);
    h3_ratio->SetLineColor(kBlue);
    
    TLine *line = new TLine(0,1,350,1);
    line->SetLineColor(kBlack);
    line->Draw("same");
    h2_ratio->Draw("same");
    h3_ratio->Draw("same");
    
    return;
}


void DrawPlot(TH1F* h1, TH1F* h2, TCanvas* ctemp){
    
    ctemp->cd();
    TPad *p1a = new TPad("p1a","",0.,.3,1.,1.);
    p1a->SetBottomMargin(0);
    p1a->Draw();
    TPad *p1b = new TPad("p1b","",0.,0.,1.,.3);
    p1b->SetTopMargin(0);
    p1b->SetBottomMargin(0.2);
    p1b->Draw();
    p1a->cd();
    p1a->SetGridx();
    p1a->SetLogy();
    //  p1a->SetGridy();
    gPad->Update();
    
    TString cmsText     = "CMS";
    float cmsTextFont   = 61;  // default is helvetic-bold
    float cmsTextSize   = 0.5;
    
    TString extraText   = "Preliminary";
    float extraTextFont = 52;  // default is helvetica-italics
    float extraOverCmsTextSize  = 0.7;
    float extraTextSize = extraOverCmsTextSize*cmsTextSize;
    
    TString lumiText = "35.9 fb^{-1}";
    
    float relPosX      = 0.045;
    float relPosY      = -0.01;
    
    float r = p1a->GetRightMargin();
    float l = p1a->GetLeftMargin();
    float posX_ =   l + relPosX*(1-l-r);
    
    float t = p1a->GetTopMargin();
    float b = p1a->GetBottomMargin();
    float posY_ = 1-t - relPosY*(1-t-b);
    
    int   align_       = 10*1+1;
    float relExtraDX   = 1.5;
    
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
    
    h1->SetLineColor(kRed);
    h1->SetMarkerStyle(kFullTriangleUp);
    h1->SetMarkerSize(0.6);
    h1->SetMarkerColor(kRed);
    h1->SetStats(kFALSE);
    
    h2->SetLineColor(kBlack);
    h2->SetMarkerStyle(kFullSquare);
    h2->SetMarkerSize(0.6);
    h2->SetMarkerColor(kBlack);
    h2->SetStats(kFALSE);

    
    h2->SetTitle("Missing Transverse Energy");
    
    h2->Draw("e");
    h1->Draw("e same");
    
    //Create a Legend for the plot
    TLegend* leg = new TLegend(0.4,0.74,0.9,0.9);
    
    leg->SetTextSize(.04);
    leg->AddEntry(h1, "Di-EM pT reweighted ee","lp");
    leg->AddEntry(h2, "Di-EM pT vs NJet reweighted ee", "lp");
    leg->Draw();
    gPad->Update();
    
    
    p1b->cd();
    
    TH1F * h2_ratio = (TH1F*) h2->Clone("h2_ratio");
    h2_ratio->Reset();
    h2_ratio->SetTitle(" ");
    h2_ratio->SetStats(0);
    h2_ratio->Sumw2();
    
    
    float ffMet, ggMet, eeMet, wvalue, binvalue;
    float erroree, errorgg, errorff, error;
    for(int i =0; i < h1->GetNbinsX()+1; i++){
        ggMet = h1->GetBinContent(i);
        eeMet = h2->GetBinContent(i);
        errorgg = h1->GetBinError(i);
        erroree = h2->GetBinError(i);
        if(eeMet == 0 || ggMet ==0){
            wvalue = 0;
            error = 0;
        }
        else{
            wvalue = ggMet/eeMet;
            error = wvalue*sqrt((errorgg/ggMet)*(errorgg/ggMet)+(erroree/eeMet)*(erroree/eeMet));
        }
        binvalue= h2_ratio->GetBinCenter(i);
        h2_ratio->Fill(binvalue, wvalue);
        h2_ratio->SetBinError(i, error);
        
    }
    
    
    p1b->SetGridx();
    p1b->SetGridy();
    h2_ratio->GetXaxis()->SetRangeUser(0,350);
    h2_ratio->SetLineColor(kBlack);
    h2_ratio->GetXaxis()->SetLabelSize(0.08);
    h2_ratio->GetYaxis()->SetLabelSize(0.08);
    h2_ratio->GetXaxis()->SetTitleSize(0.08);
    h2_ratio->GetYaxis()->SetTitleOffset(0.5);
    h2_ratio->GetYaxis()->CenterTitle();
    h2_ratio->GetYaxis()->SetTitleSize(0.08);
    h2_ratio->GetYaxis()->SetRangeUser(0.5,2.5);
    h2_ratio->SetMarkerStyle(kFullSquare);
    h2_ratio->SetMarkerColor(kBlack);
    h2_ratio->SetMarkerSize(0.4);
    
    h2_ratio->Draw();
    

    
    TLine *line = new TLine(0,1,350,1);
    line->SetLineColor(kBlack);
    line->Draw("same");
    h2_ratio->Draw("same");
    
    return;
}
