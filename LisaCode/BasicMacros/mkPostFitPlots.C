void mkPostFitPlots(){

Double_t total_bkg[6]={64.6,27.68,19.59,15.86,8.71,8.55};
Double_t EWK_bkg[6]={7.51,5.05,4.38,3.61,3.23,1.95};
Double_t QCD_bkg[6]={55.85,21.54,14.14,10.99,4.26,5.52};
Double_t ZGG_bkg[6]={1.24,1.09,1.07,1.26,1.22,1.09};

Double_t tot_err[6]     = {5.95,4.37,3.43,3.24,2.19,2.87};
Double_t tot_err_EWK[6] = {21.38,1.21,1.01,0.94,0.56};
Double_t tot_err_QCD[6] = {6.35,4.59,3.12,2.23,2.93};
Double_t tot_err_ZGG[6] = {0.6,0.53,0.52,0.61,0.59,5.29};

Double_t   channel[6]={1,2,3,4,5,6};
Double_t  err_x[6] = {0.5,0.5,0.5,0.5,0.5,0.5};

/////////////PreFit///////


Double_t Pre_total_bkg[6]={78.7,37.53,31.88,25.76,13.54,8.38};
Double_t Pre_EWK_bkg[6]={8.17,5.5,4.78,3.95,3.52,2.11};
Double_t Pre_QCD_bkg[6]={69.23,30.89,25.98,20.49,8.74,5.13};
Double_t Pre_ZGG_bkg[6]={1.3,1.14,1.12,1.32,1.28,1.14};


Double_t Pre_tot_err[6]     = {21.09,17.05,15.28,11.51,15.07, 18.69};
Double_t Pre_tot_err_EWK[6] = {2.31,1.51,1.32,1.11,1,0.63};
Double_t Pre_tot_err_QCD[6] = {6.35,4.59,3.42,3.12,2.23,2.93};
Double_t Pre_tot_err_ZGG[6] = {0.65,0.57,0.56,0.66,0.64,0.57};


TCanvas *c1 = new TCanvas("c1","c1",1000,1000);

c1->Divide(2,2);
c1->cd(1);
TGraphErrors *gr      = new TGraphErrors(6,channel,total_bkg,err_x,tot_err);
TGraphErrors *Pregr   = new TGraphErrors(6,channel,Pre_total_bkg,err_x,Pre_tot_err);
TMultiGraph *mg = new TMultiGraph();
mg->Add(gr,"AP");
mg->Add(Pregr,"AP");
mg->Draw("a");
   //gr->Draw("AP");
   mg->SetTitle("Total Background");
   gr->SetMarkerColor(kRed);
   gr->SetLineColor(kRed);
   Pregr->SetMarkerColor(kBlue);
   Pregr->SetLineColor(kBlue);
   //mg->SetMarkerStyle(22);
   //mg->SetMarkerSize(0.9);
   mg->GetXaxis()->SetTitle("Channel");

   mg->GetYaxis()->SetTitleOffset(0.4);


 TLegend *leg = new TLegend(0.58, 0.80, 0.899, 0.899);
      leg->SetFillColor(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.03);
      leg->AddEntry( Pregr, "Pre Fit Total background", "LP" );
      leg->AddEntry( gr, "Post Fit Total background", "LP" );
      leg->Draw();

c1->cd(2);
TGraphErrors * grEWK= new TGraphErrors(6,channel,EWK_bkg,err_x,tot_err_EWK);
TGraphErrors *PregrEWK= new TGraphErrors(6,channel,Pre_EWK_bkg,err_x,Pre_tot_err_EWK);

TMultiGraph *mgEWK = new TMultiGraph();

mgEWK->Add(PregrEWK,"AP");
mgEWK->Add(grEWK,"AP");
mgEWK->Draw("a");


   //grEWK->Draw("AP");
   mgEWK->SetTitle(" EWK Background");
   grEWK->SetMarkerColor(kRed);
   grEWK->SetLineColor(kRed);
   PregrEWK->SetMarkerColor(kBlue);
   PregrEWK->SetLineColor(kBlue);
   //grEWK->SetMarkerStyle(21);
   //grEWK->SetMarkerSize(0.9);
   mgEWK->GetXaxis()->SetTitle("Channel");

   grEWK->GetYaxis()->SetTitleOffset(0.4);

   TLegend *leg1 = new TLegend(0.58, 0.80, 0.899, 0.899);
      leg1->SetFillColor(0);
      leg1->SetTextFont(42);
      leg1->SetTextSize(0.03);
      leg1->AddEntry( PregrEWK, "Pre Fit EWK background", "LP" );
      leg1->AddEntry( grEWK, "Post Fit EWK background", "LP" );
      leg1->Draw();



c1->cd(3);


TGraphErrors *grQCD= new TGraphErrors(6,channel,QCD_bkg,err_x,tot_err_QCD);
TGraphErrors *PregrQCD= new TGraphErrors(6,channel,Pre_QCD_bkg,err_x,Pre_tot_err_QCD);

TMultiGraph *mgQCD = new TMultiGraph();
mgQCD->Add(grQCD,"AP");
mgQCD->Add(PregrQCD,"AP");

mgQCD->Draw("a");


  // grQCD->Draw("AP");
   mgQCD->SetTitle("QCD Background");
   grQCD->SetMarkerColor(kRed);
   grQCD->SetLineColor(kRed);
   PregrQCD->SetMarkerColor(kBlue);
   PregrQCD->SetLineColor(kBlue);

   //grQCD->SetMarkerStyle(20);
   //grQCD->SetMarkerSize(0.9);
   mgQCD->GetXaxis()->SetTitle("Channel");
   grQCD->GetYaxis()->SetTitleOffset(0.4);



   TLegend *leg2 = new TLegend(0.58, 0.80, 0.899, 0.899);
      leg2->SetFillColor(0);
      leg2->SetTextFont(42);
      leg2->SetTextSize(0.03);
      leg2->AddEntry( PregrQCD, "Pre Fit QCD background", "LP" );
      leg2->AddEntry( grQCD, "Post Fit QCD background", "LP" );
      leg2->Draw();



c1->cd(4);

   TGraphErrors *grZGG   = new TGraphErrors(6,channel,ZGG_bkg,err_x,tot_err_ZGG);
   TGraphErrors *PregrZGG= new TGraphErrors(6,channel,Pre_ZGG_bkg,err_x,Pre_tot_err_ZGG);

   TMultiGraph *mgZGG = new TMultiGraph();
   mgZGG->Add(grZGG,"AP");
   mgZGG->Add(PregrZGG,"AP");

   mgZGG->Draw("a");

   mgZGG->Draw("AP");
   mgZGG->SetTitle("ZGG Background");
   grZGG->SetMarkerColor(kRed);
   grZGG->SetLineColor(kRed);
   PregrZGG->SetMarkerColor(kBlue);
   PregrZGG->SetLineColor(kBlue);

   //grZGG->SetMarkerStyle(23);
   //grZGG->SetMarkerSize(0.9);
   mgZGG->GetXaxis()->SetTitle("Channel");
   grZGG->GetYaxis()->SetTitleOffset(0.4);
   //grZGG->Draw("AP");


   TLegend *leg3 = new TLegend(0.58, 0.80, 0.899, 0.899);
      leg3->SetFillColor(0);
      leg3->SetTextFont(42);
      leg3->SetTextSize(0.03);
      leg3->AddEntry( PregrZGG, "Pre Fit ZGG background", "LP" );
      leg3->AddEntry( grZGG, "Post Fit ZGG background", "LP" );
      leg3->Draw();



}
