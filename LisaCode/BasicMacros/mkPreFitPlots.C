void mkPreFitPlots(){

Double_t total_bkg[6]={78.7,37.53,31.88,25.76,13.54,8.38};
Double_t EWK_bkg[6]={8.17,5.5,4.78,3.95,3.52,2.11};
Double_t QCD_bkg[6]={69.23,30.89,25.98,20.49,8.74,5.13};
Double_t ZGG_bkg[6]={1.3,1.14,1.12,1.32,1.28,1.14};


Double_t tot_err[6]     = {21.09,17.05,15.28,11.51,15.07, 18.69};
Double_t tot_err_EWK[6] = {2.31,1.51,1.32,1.11,1,0.63}; 
Double_t tot_err_QCD[6] = {6.35,4.59,3.42,3.12,2.23,2.93};
Double_t tot_err_ZGG[6] = {0.65,0.57,0.56,0.66,0.64,0.57};

Double_t   channel[6]={1,2,3,4,5,6};
Double_t  err_x[6] = {0.5,0.5,0.5,0.5,0.5,0.5};

TGraph* gr = new TGraph(6,channel,total_bkg);

TCanvas *c1 = new TCanvas("c1","c1",1000,1000);

c1->Divide(2,2);
c1->cd(1);
gr   = new TGraphErrors(6,channel,total_bkg,err_x,tot_err);


   gr->Draw("AP");
   gr->SetTitle("Pre Fit Total Background");
   gr->SetMarkerColor(kRed);
   gr->SetLineColor(kRed);
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(0.9);
   gr->GetXaxis()->SetTitle("Channel");
  
   gr->GetYaxis()->SetTitleOffset(0.4);


 TLegend *leg = new TLegend(0.58, 0.80, 0.899, 0.899);
      leg->SetFillColor(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.03);
      leg->AddEntry( gr, "Total background", "P" );
     // leg->Draw();

c1->cd(2);

grEWK= new TGraphErrors(6,channel,EWK_bkg,err_x,tot_err_EWK);


   grEWK->Draw("AP");
   grEWK->SetTitle("Pre Fit EWK Background");
   grEWK->SetMarkerColor(kBlue);
   grEWK->SetLineColor(kBlue);
   grEWK->SetMarkerStyle(21);
   grEWK->SetMarkerSize(0.9);
   grEWK->GetXaxis()->SetTitle("Channel");
  
   grEWK->GetYaxis()->SetTitleOffset(0.4);



c1->cd(3);


grQCD= new TGraphErrors(6,channel,QCD_bkg,err_x,tot_err_QCD);
   grQCD->Draw("AP");
   grQCD->SetTitle("Pre Fit QCD Background");
   grQCD->SetMarkerColor(kGreen);
   grQCD->SetLineColor(kGreen);
   grQCD->SetMarkerStyle(20);
   grQCD->SetMarkerSize(0.9);
   grQCD->GetXaxis()->SetTitle("Channel");
   grQCD->GetYaxis()->SetTitleOffset(0.4);


c1->cd(4);

   grZGG= new TGraphErrors(6,channel,ZGG_bkg,err_x,tot_err_ZGG);
   grZGG->Draw("AP");
   grZGG->SetTitle("Pre Fit ZGG Background");
   grZGG->SetMarkerColor(kViolet);
   grZGG->SetLineColor(kViolet);
   grZGG->SetMarkerStyle(23);
   grZGG->SetMarkerSize(0.9);
   grZGG->GetXaxis()->SetTitle("Channel");
   grZGG->GetYaxis()->SetTitleOffset(0.4);


   grZGG->Draw("AP");

}
