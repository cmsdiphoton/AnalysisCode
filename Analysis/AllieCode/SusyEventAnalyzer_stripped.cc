//
//  Susyeventanalyzer.cc
//
//
//  Created by Allie Reinsvold Hall on 8/4/14.
//
// Stripped version with fewer unnecessary histograms
// Includes full analysis code 
// Updated September 2016 
 
#define SusyEventAnalyzer_cxx

#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TGraphAsymmErrors.h>
#include <TRandom.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <map>
#include <set>
#include <cmath>
#include <utility>
#include <fstream>
#include "lester_mt2_bisect.h"

#include "SusyEventAnalyzer.h"

using namespace std;

void SusyEventAnalyzer::Loop() {
  bool version74X = false;
  cout<<"Inside Loop(). WOOOOOOO!"<<endl;
  cout<<"Did you remember to compile any changes to SusyEventAnalyzer.cc?"<<endl;
  if (fChain == 0) {
    cout << "fChain ==0"<< endl;
    return;
  }
  Long64_t nentries = fChain->GetEntries();
  new ggEventTree(fChain);
  event = new ggEventTree(fChain);
  cout << "total events in files  : " << nentries << endl;
  if(processNEvents <= 0 || processNEvents > nentries) processNEvents = nentries;
    
  cout << "events to be processed : " << processNEvents << endl;
  //array to keep track of cut flow of events
  const int NCNT = 20;
  int nCnt[NCNT];
  for(int i=0; i<NCNT; i++) nCnt[i] = 0;

    
  if(printLevel>0) cout<<"Open hist file" << endl;
  TFile* fout = new TFile("trig_"+ds+".root","RECREATE");
  if(printLevel>0) cout<<"Hist file opened" << endl;
    
  fout->cd();
    
  TH1F* h_phoR9 = new TH1F("phoR9","Photon R9;R9;Events", 105,0.0,1.05);
  TH1F* h_phoR9_PT40 = new TH1F("phoR9_PT40","Photon R9 if pT > 40;R9;Events", 105,0.0,1.05);
  TH1F* h_phoPT = new TH1F("phoPT","Photon Transverse Momentum;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_GenPhoPT = new TH1F("GenPhoPT","Truth-level Photon pT;p_{T} (GeV);Events", 200,0.,1000.);
  TH1F* h_puTrue = new TH1F("puTrue","Number of PU interactions in events with two photons; Number; Events",100,0.,100.);
  TH1F* h_ggNVertex = new TH1F("ggNVertex","Number of vertices in events with two photons;N vertex;Events",75,0,150);

  // start event looping
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < processNEvents; jentry++) {
    

    if(jentry%10000==0){
      cout<<"Processing event number " << jentry << " of " << processNEvents <<" events." <<endl;
    }

    nb = event->GetEntry(jentry);   nbytes += nb;
        
    if(event->run_>=0){
      nCnt[0]++; // total number of events
            
    }
        
    int nGenPhotons =0;
    for(int Part_it = 0; Part_it < event->nMC_; Part_it++){
      float partPdgID = event->mcPID->at(Part_it);
      if(fabs(partPdgID)==22){
	if(fabs(event->mcMomPID->at(Part_it)) == 1000023){ //10023 = neutralino 
	  nGenPhotons++;
	  h_GenPhoPT->Fill(event->mcPt->at(Part_it));

	  TLorentzVector part_it_vec = MassLorentzVector(event->mcPt->at(Part_it), event->mcEta->at(Part_it), event->mcPhi->at(Part_it), event->mcMass->at(Part_it));

          for(int it_pho = 0; it_pho<event->nPho_; it_pho++){
            TLorentzVector candVector = PhoLorentzVector(event->phoEt_->at(it_pho),
                                                         event->phoEta_->at(it_pho),
                                                         event->phoPhi_->at(it_pho),
                                                         event->phoCalibE_->at(it_pho));

	    if(getDR(candVector,part_it_vec) < 0.15){
	      h_phoPT->Fill(event->phoEt_->at(it_pho));
	      h_phoR9->Fill( event->phoR9_->at(it_pho) );
	      if(event->phoEt_->at(it_pho) > 40) h_phoR9_PT40->Fill( event->phoR9_->at(it_pho) );
	    }
	  }//end of reco photon loop

	} // end of if mom is a neutralino
      }
    } // end of mc particles loop

    if (nGenPhotons >= 2){
      h_ggNVertex->Fill(event->nVtx_);
  
      int pu_it = 0; 
      for(vector<int>::const_iterator bx = event->puBX_->begin(); bx != event->puBX_->end(); bx++){
	if( *bx == 0){
          h_puTrue->Fill(event->puTrue_->at(pu_it));
        }
        pu_it++;
      }

    }


  }//end of looping over the events

  cout<<"Writing analysis root output to: hist_"<<ds<<".root"<<endl;
  // close the output file
    
  fout->cd();
  fout->Write();
  fout->Close();
  delete fout;
    
    
}//end of Loop() function





//---------------------------------------------------------------




//void SusyEventAnalyzer::CategorizeEvents(int pho1, int pho2, float Rho, bool &gogg, bool &goee, bool &goeg, bool &goff, bool &gogammafake, bool &gogf, bool &gofg){
//return;
//}



void SusyEventAnalyzer::CategorizeEvents(int pho1, int pho2, int pho3, float Rho, bool &gogg, bool &goee, bool &goeg, bool &goff, bool &gogammafake, bool &gogf, bool &gofg, bool &goef, bool &use13, bool&use23){
    //All of the photon candidates have had to pass cuts on H over E, photon isolation, and neutral hadron isolation
    //electrons and photons must also pass sigma ieta ieta and charged hadron isolation cuts. Photons must not have a pixel seed.
    //initialize all booleans to false
    bool g1=false, g2=false, f1=false, f2=false, e1=false, e2=false;
    gogg=false;goee=false;goeg=false;goff=false;gogammafake=false;gogf=false;gofg=false;
    use23 = false, use13 = false; 
    if((pho1>=event->nPho_)||(pho2>=event->nPho_)){
        cout << "Photon index higher than the number of photons in the event!! Error!!"<<endl;
        cout<< "PhoOne= " <<pho1<< "  PhoTwo= "<<pho2<<"  nPho= "<<event->nPho_<<endl;
        return;
    }
    
    float maxSihih = 0.01022; //.03; //.0106;
    float maxChHadIso = 0.441;
    
    //Define the effective area for pho1 (see https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2 )
    
    float pho1Eta = abs(event->phoEta_->at(pho1) );
    float chIsoEA1=0;
    
    if(pho1Eta<1.0){
        chIsoEA1 = 0.0360;
    }
    else if (pho1Eta <1.479 && pho1Eta >=1.0){
        chIsoEA1 = 0.0377;
    }
    else if (pho1Eta <2.0 && pho1Eta >=1.479){
        chIsoEA1 = 0.0306;
    }
    else if (pho1Eta <2.2 && pho1Eta >=2.0){
        chIsoEA1 = .0283;
    }
    else if (pho1Eta <2.3 && pho1Eta >=2.2){
        chIsoEA1 = .0254;
    }
    else if (pho1Eta <2.4 && pho1Eta >=2.3){
        chIsoEA1 = .0217;
    }
    else if (pho1Eta >=2.4){
        chIsoEA1 = .0167;
    }
    
    //Define EA for pho 2
    
    float pho2Eta = abs(event->phoEta_->at(pho2) );
    float chIsoEA2=0.0;
    
    if(pho2Eta<1.0){
        chIsoEA2 = .0360;
    }
    else if (pho2Eta <1.479 && pho2Eta >=1.0){
        chIsoEA2 = .0377;
    }
    else if (pho2Eta <2.0 && pho2Eta >=1.479){
        chIsoEA2 = .0306;
    }
    else if (pho2Eta <2.2 && pho2Eta >=2.0){
        chIsoEA2 = .0283;
    }
    else if (pho2Eta <2.3 && pho2Eta >=2.2){
        chIsoEA2 = .0254;
    }
    else if (pho2Eta <2.4 && pho2Eta >=2.3){
        chIsoEA2 = .0217;
    }
    else if (pho2Eta >=2.4){
        chIsoEA2 = .0167;
    }
    
    //Sort pho1
    float pho1iso  = (event->phoPFChIso_->at(pho1) - Rho*chIsoEA1) > 0.00 ? (event->phoPFChIso_->at(pho1) - Rho*chIsoEA1) : 0.00;
    if(event->phoSigmaIEtaIEtaFull5x5_->at(pho1) < maxSihih && pho1iso < maxChHadIso){
        if(event->phohasPixelSeed_->at(pho1)){
            e1 = true;
        }
        else{
            g1 = true;
        }
    }
    else{  //since all the photon candidates are either electrons, photons, or fakes, if it failed the sigma ieta ieta or the charged hadron isolation cut, it must be a fake
        f1 = true;
    }
    
    //Sort pho2
    float pho2iso  = (event->phoPFChIso_->at(pho2) - Rho*chIsoEA2) > 0.00 ? (event->phoPFChIso_->at(pho2) - Rho*chIsoEA2) : 0.00;
    if(event->phoSigmaIEtaIEtaFull5x5_->at(pho2) < maxSihih && pho2iso < maxChHadIso){
        if(event->phohasPixelSeed_->at(pho2)){
            e2 = true;
        }
        else{
            g2 = true;
        }
    }
    else{  //since all the photon candidates are either electrons, photons, or fakes, if it failed the sigma ieta ieta or the charged hadron isolation cut, it must be a fake
        f2 = true;
    }
    if( (e1 && f2) || (e2 && f1))   goef=true;
    if(g1 && g2)                    gogg=true; 
    if( (g1 && e2) || (e1 && g2))   goeg=true;
    if(e1 && e2)                    goee=true; 
    if(f1 && f2)                    goff=true;
    if( (g1 && f2) || (f1 && g2) ){
      gogammafake=true;
      if(g1 && f2)    gogf=true;
      else if(f1 && g2)   gofg=true;
    }

    if( ! (gogg || goee || goff || goeg) ) {
      if(pho3 > -1){
	//if the first two photons didn't fall into a useful category, use the third photon object
	float chIsoEA3 = 0.0; 
	float pho3iso  = (event->phoPFChIso_->at(pho3) - Rho*chIsoEA3) > 0.00 ? (event->phoPFChIso_->at(pho3) - Rho*chIsoEA3) : 0.00;
	bool g3 = false, e3 = false, f3 = false;
	if(event->phoSigmaIEtaIEtaFull5x5_->at(pho3) < maxSihih && pho3iso < maxChHadIso){
	  if(event->phohasPixelSeed_->at(pho3)){
            e3 = true;
	  }
	  else{
            g3 = true;
	  }
	}
	else{  //since all the photon candidates are either electrons, photons, or fakes, if it failed the sigma ieta ieta or the charged hadron isolation cut, it must be a fake      
	  f3 = true;
	}

	if(g1 && g3){ gogg = true; use13 = true;}
	else if( (g1 && e3) || (e1 && g3) ){ goeg = true; use13 = true;}
	else if( e1 && e3) {goee = true; use13 = true;}
	else if(f1 && f3) {goff = true; use13 = true;}
	else if( g2&&g3) {gogg = true; use23 = true;}
	else if( (g2 && e3) || (e2 && g3) ){ goeg = true; use23 = true;}
	else if (e2 && e3){ goee = true; use23 = true;}
	else if (f2 && f3 ){ goff = true;  use23 = true;}
      }
    }
    
    if(gogg&&goeg)cout<<"gg AND eg event!!!!! ------- PROBLEM!"<<endl;
    if(gogg&&goee)cout<<"gg AND ee event!!!!! ------- PROBLEM!"<<endl;
    if(gogg&&goff)cout<<"gg AND ff event!!!!! ------- PROBLEM!"<<endl;
    //    if(gogg&&gogammafake)cout<<"gg AND gammafake event!!!!! ------- PROBLEM!"<<endl;
    if(goeg&&goee)cout<<"eg AND ee event!!!!! ------- PROBLEM!"<<endl;
    if(goeg&&goff)cout<<"eg AND ff event!!!!! ------- PROBLEM!"<<endl;
    //    if(goeg&&gogammafake)cout<<"eg AND gammafake event!!!!! ------- PROBLEM!"<<endl;
    if(goee&&goff)cout<<"ee AND ff event!!!!! ------- PROBLEM!"<<endl;
    //    if(goee&&gogammafake)cout<<"ee AND gammafake event!!!!! ------- PROBLEM!"<<endl;
    //    if(goff&&gogammafake)cout<<"ff AND gammafake event!!!!! ------- PROBLEM!"<<endl;
    if(gogf&&gofg)cout<<"gf AND fg event!!!!! ------- PROBLEM!"<<endl;
    if((gogammafake && !gogf) && (gogammafake && !gofg))cout<<"gammafake BUT NOT gf OR fg event!!!!! ------- PROBLEM!"<<endl;
    return;
} //end of Categorize Event function



void SusyEventAnalyzer::MatchPhosToJets(TLorentzVector pOne, TLorentzVector pTwo, std::vector<TLorentzVector*> jetvecs, TLorentzVector &jet1, TLorentzVector &jet2, bool &hasdijetpt, float dR){
    hasdijetpt=false;
    for(std::vector<TLorentzVector*>::iterator jet_it1 = jetvecs.begin(); jet_it1 != jetvecs.end(); jet_it1++){
        if(isSameObject(pOne,*(*jet_it1),dR)){
            for(std::vector<TLorentzVector*>::iterator jet_it2 = jetvecs.begin(); jet_it2 != jetvecs.end(); jet_it2++){
                if( !isSameObject(*(*jet_it1), *(*jet_it2),0.1) ){
                    if(isSameObject(pTwo, *(*jet_it2), dR)){
                        jet1=**jet_it1;
                        jet2=**jet_it2;
                        hasdijetpt=true;
                        break;
                    }//end jet2 match
                }//end jet2!=jet1
            }//end jet2 iterator
            if(hasdijetpt){break;}
        }//end jet1 match
    }//end jet1 iterator
    return;
}//end of MatchPhosToJets function
