//
//  SusyEventAnalyzer.h
//  
//
//  Created by Allie Reinsvold on 8/4/14.
//
//

#ifndef SusyEventAnalyzer_h
#define SusyEventAnalyzer_h

#include <iostream>
#include <fstream>
#include <string>
#include<map>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TVector2.h>
#include <TLorentzVector.h>

#include<vector>
#include<algorithm>

//#include "../../../TMVA-v4.1.2/TMVA/Tools.h"
//#include "../../../TMVA-v4.1.2/TMVA/Reader.h"
//#include "../../../TMVA-v4.1.2/TMVA/MethodCuts.h"
#include "../src/ggEventTree.h"

class SusyEventAnalyzer{
    public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain
    
    
    ggEventTree* event;
    
    SusyEventAnalyzer(TTree *tree=0);
    virtual ~SusyEventAnalyzer();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    
    virtual void     Loop();                          // event loop for main analysis
    
     void Initialize();         // global variables needed to be initialized just once
    
//    virtual void     DR03();
//    virtual void     Pileup();
//    virtual void     Filter();
//    virtual void     PhotonId();
    
    //... missing functions:
//    float d0correction(TVector3& beamSpot, susy::Track& track) const;
    float GetDiJetPt(TLorentzVector Jet1, TLorentzVector Jet2);//calculates and returns diJetPt
//    float GetTriEmPt(susy::Photon* Pho1, susy::Photon* Pho2, susy::Photon* Pho3);//calculates and returns TriEmPt

//    float GetRazrMr(susy::Photon* Pho1, susy::Photon* Pho2);
//    float GetRazrR2(susy::Photon* Pho1, susy::Photon* Pho2, susy::MET* Met);
//    std::vector<float> GetMetAndInvMassProbs(float invmass, float met, TH1F* SMmet, TH1F* SMinvmass, TH1F* BGmet, TH1F* BGinvmass, TH1F* NewHiggsMet);
    
    //utility functions
    int   key(int i, int j);    
    float GetDiEmPt(TLorentzVector Pho1, TLorentzVector Pho2);
    void IncludeAJson(std::string jsonfile);  // Call to pull in a json file
    bool isInJson(Int_t run,Int_t lumi);      // JSON based good run list cut...
    float getDphi(float p1, float p2);
    void MatchPhosToJets(TLorentzVector pOne, TLorentzVector pTwo, std::vector<TLorentzVector*> jetvecs, TLorentzVector &jet1, TLorentzVector &jet2, bool &hasdijetpt, float dR);
    float InvariantMass(TLorentzVector P1, TLorentzVector P2); //returns the invariant mas of the two object system
    float InvariantMass(Float_t ggInvMass, TLorentzVector part, TVector2 met);//calculates and returns Invariant Mass
    float InvariantMass(TLorentzVector P1, TLorentzVector P2, TLorentzVector P3);//calculates and returns Invariant Mass of three object system
    TLorentzVector PhoLorentzVector(const Double_t pt, const Double_t eta, const Double_t phi, const Double_t e);
    TLorentzVector MassLorentzVector(Double_t pt, Double_t eta, Double_t phi, Double_t m);
    bool isSameObject(TLorentzVector& p1, TLorentzVector& p2, float dR_Cut);
    bool tooClosePhi(TLorentzVector& p1, TLorentzVector& p2, float phi_Cut);
    float getDR(TLorentzVector& p1, TLorentzVector& p2);
    float TransverseMass(TLorentzVector part, TVector2 met);//calculates and returns Transverse Mass
    float TransverseMassSquare3body(TLorentzVector lep1, TLorentzVector lep2, TVector2 met);//calculates and returns Transverse Mass Squared for 3 object system
    float MT2(TLorentzVector P1, TLorentzVector P2, TLorentzVector P3);//calculates and returns mt2 for 3 body system, e.g. (e+)(e-)(gamma)
    float GetPhotonLessHt(float Ht, TLorentzVector pOne, TLorentzVector pTwo);
    
    
    void CategorizeEvents(int pho1, int pho2, int pho3, float Rho, bool &gogg, bool &goee, bool &goeg, bool &goff, bool &gogammafake, bool &gogf, bool &gofg, bool &goef, bool &use13, bool &use23, bool& goffmix, bool& goffhigh);
    void CategorizeEvents(int pho1, int pho2, float Rho, bool &gogg, bool &goee, bool &goeg,bool &goff, bool &gogammafake, bool &gogf, bool &gofg);
    void CategorizeEvents(int pho1, int pho2, int pho3, float Rho, bool &gogg, bool &goee, bool &goeg, bool &goff, bool &gogammafake, bool &gogf, bool &gofg, bool &goef, bool &use13, bool&use23);


    //parameter configuration functions
    void SetDataset(TString& v) {          ds = v; }
    void SetPrintInterval(int v) {         printInterval = v; }
    void SetPrintLevel(int v) {            printLevel = v; }
    void SetProcessNEvents(int v) {        processNEvents = v; }
    void SetDoToys(bool v) {               doToys = v; }
    void SetDoTTBar(bool v) {              doTTBar = v;}
    void SetUseTrigger(bool v) {           useTrigger = v; }
    void SetUseJSON(bool v) {              useJSON = v; }
    void AddHltName(TString v) {           hltNames.push_back(v); }
    void SetFilter(bool v) {               enableFilter = v; }
    void SetFilteredFileName(TString v) {  filtered_file_name = v; }
    void SetOutputEventNumbers(bool v) {   outputEventNumbers = v; }
    void DoRhoCorrection(bool v) {         doRhoCorrection = v; }
    void DoNvertexCorrection(bool v) {     doNVertexCorrection = v; }
    void SetDR03Rho25Corr(float ecal, float hcal, float track){ PUCorr_ECAL=ecal;PUCorr_HCAL=hcal;PUCorr_TRACK = track; }
    void SetPFisoRho25Corr(float ch, float nh, float ph){ PUCorr_chargedHadron=ch;PUCorr_neutralHadron=nh;PUCorr_photon=ph; }
    void isFastSim(bool isFast){FastSim=isFast;}
    void isFullSim(bool isFull){FullSim=isFull;}
    void SetPUFileName(TString filename) {puFileName = filename; }
    void SetDoSignal(bool sig){ doSignal = sig;}
    void SetDoT5Wg(bool t5wg){T5Wg = t5wg; }
    //Missing functions:
//    float FastSimSmear(susy::Photon* pho, TRandom* rand);
//    float FullSimSmear(susy::Photon* pho, TRandom* rand);
//    susy::PFJet* JECup  (susy::PFJet* jet);
//    susy::PFJet* JECdown(susy::PFJet* jet);
//    TVector2 CalcMet(std::vector<susy::PFJet*> jets, susy::Photon* pho1, susy::Photon* pho2, std::vector<susy::Electron*> eles, std::vector<susy::Muon*> mus);
//    TVector2 CalcMetFromPFandJets(std::vector<susy::PFJet*> jets, TVector2 pfMet);
//    float GetPhoScaleFactor(susy::Photon* pho);
//    float GetEleScaleFactor(susy::Electron* ele);
//    float GetMuScaleFactor(susy::Muon* mu);


    private:
    TString ds;               // dataset name to be used for output hitfile name
    
  
    
    // printLevel
    // 0 : default - no printout
    // 1 : print functional step in every event
    // 2 : print values in collections
    int printLevel;           // print level for event content: defined in Event.h
    bool outputEventNumbers;  // print run and event numbers for gg, eg, ee, ff to txt
    int printInterval;        // print frequency
    int processNEvents;       // number of events to be processed
    
    TString puFileName = "powhegPU.root";
    
    bool T5Wg = true;
    bool doSignal = false;
    bool doTTBar = false;
    bool doToys;
    bool useTrigger;          // flag for using trigger bit selection.
    bool useJSON;             // flag for using JSON selection
    std::vector<TString> hltNames;          // HLT trigger path names
    bool enableFilter;        // filter events of interest
    TString filtered_file_name; // filtered output file name
    bool doRhoCorrection;
    bool doNVertexCorrection;
    float PUCorr_ECAL;
    float PUCorr_HCAL;
    float PUCorr_TRACK;
    float PUCorr_chargedHadron;
    float PUCorr_neutralHadron;
    float PUCorr_photon;
    bool FastSim;
    bool FullSim;
    
    typedef std::map<int,std::map<int,bool> > RunLumiFlagHolder;  //define map that holds json list
    RunLumiFlagHolder goodrunlumilist;  // instantiate it
    
};

#endif

#ifdef SusyEventAnalyzer_cxx
SusyEventAnalyzer::SusyEventAnalyzer(TTree *tree)
{
    if (tree == 0) {
        std::cout << "Error!!! There is no file containing a tree." << std::endl;
    }
    Init(tree);
    Initialize();
}

SusyEventAnalyzer::~SusyEventAnalyzer()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

template<typename T> bool EtGreater(const T* p1, const T* p2) {
  return (p1 < p2);
}


void SusyEventAnalyzer::Init(TTree *tree)
{
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    //fChain->SetMakeClass(1);
    
    //event = new susy::Event;
    
    //event->setInput(*fChain);

}


Long64_t SusyEventAnalyzer::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (!fChain->InheritsFrom(TChain::Class()))  return centry;
    TChain *chain = (TChain*)fChain;
    if (chain->GetTreeNumber() != fCurrent) {
        fCurrent = chain->GetTreeNumber();
    }
    return centry;
}


Int_t SusyEventAnalyzer::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry,0);
}

void SusyEventAnalyzer::Initialize() {
    
    ds = "test";
    printLevel = 0;
    printInterval = 1000;
    processNEvents = -1;
    
}

float SusyEventAnalyzer::getDphi(float p1, float p2)
{
    float dPhi = std::fabs(TVector2::Phi_mpi_pi(p1 - p2));
    return dPhi;
}
void SusyEventAnalyzer::IncludeAJson(std::string jsonfile) {


  // Fairly primitive brute force json parser -- opens the json file named in the argument                
  // and adds that to the goodrunlumilist map.  Overlapping jsons are merged inclusively.                 

  char thing;

  ifstream jsonInput;

   std::cout << "Sucking in Json file: " << jsonfile << " which includes: " << std::endl;

  jsonInput.open(jsonfile.c_str());

  if (!jsonInput.good()) {
    std::cout << "Problem reading Json file...  Didn't suck anything in... " << std::endl;
    return;
  }
  jsonInput.get(thing);

  while (jsonInput.good()) {
    if (thing=='{') {  // start of list                                                                  
      while (thing != '}') {
	int runnum;
	if (thing == '"') {
	  std::string srunnum;
	  jsonInput.get(thing); // get stuff inside ""                                                   

	  while (thing != '"') {
	    srunnum+=thing; // get stuff inside ""                                                       
	    jsonInput.get(thing);

	  }
	  sscanf(srunnum.c_str(),"%i",&runnum);
	  //	  std::cout << " runnum: " << runnum << std::endl;
	  bool newrun=true;

	} // inside ""                                                                                   
	if (thing == '[') {
          jsonInput.get(thing); // get stuff inside []                                                  
	  while (thing != ']') {
	    if (thing == '[') {
	      jsonInput.get(thing); // get stuff inside series []                                        

	      std::string lumiseries;
	      int firstlumi,lastlumi;
	      while (thing !=']') {
		lumiseries+=thing;
                jsonInput.get(thing); // get stuff inside series []                                     
	      }
	      sscanf(lumiseries.c_str(),"%i,%i",&firstlumi,&lastlumi);
	      //  std::cout << "  lumis  " << firstlumi << " to " << lastlumi << std::endl;

	      // At this point have runnum, first lumi, last lumi -- so can fill map here...             
	      for (int l=firstlumi;l<=lastlumi;l++) {
		goodrunlumilist[runnum][l]=true;
	      }

	    } // inside actual series []                                                                 
	    jsonInput.get(thing); // get stuff inside []                                               
	  }
	} // inside []                                         
	jsonInput.get(thing); // get another char looking for "
      }
    } // inside {}
    jsonInput.get(thing); // get another char looking for {                    
  } // EOF                                                                     
  jsonInput.close();
}

bool SusyEventAnalyzer::isInJson(Int_t run, Int_t lumi){
       
    if (goodrunlumilist[run][lumi]) return true;
    
    return false;
    
}

float SusyEventAnalyzer::InvariantMass(TLorentzVector P1, TLorentzVector P2){
    float InvarMass = (P1 + P2).M();
    return InvarMass;
}

float SusyEventAnalyzer::InvariantMass(TLorentzVector P1, TLorentzVector P2, TLorentzVector P3){
    float InvarMass = (P1 + P2 + P3).M();
    return InvarMass;
}

TLorentzVector SusyEventAnalyzer::PhoLorentzVector(const Double_t pt, const Double_t eta, const Double_t phi, const Double_t e){
  TLorentzVector* vec = new TLorentzVector();
    vec->SetPtEtaPhiE(pt, eta, phi, e);
    return *vec;
}

TLorentzVector SusyEventAnalyzer::MassLorentzVector(Double_t pt, Double_t eta, Double_t phi, Double_t m){
  TLorentzVector* vec = new TLorentzVector();
    vec->SetPtEtaPhiM(pt, eta, phi, m);
    return *vec;
}

bool SusyEventAnalyzer::isSameObject(TLorentzVector& p1, TLorentzVector& p2, float dR_Cut) {
    float dEta = p1.Eta() - p2.Eta();
    float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
    float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
    if(dR < dR_Cut) return true;
    return false;
}

bool SusyEventAnalyzer::tooClosePhi(TLorentzVector& p1, TLorentzVector& p2, float phi_Cut) {
    float dPhi = std::fabs(TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi()));
    if(dPhi<phi_Cut){
        //std::cout<<"p1Phi= "<<p1.Phi()<<"  p2Phi= "<<p2.Phi()<<"  dPhi= "<<dPhi<<endl;
        return true;
    }
    return false;
}

float SusyEventAnalyzer::getDR(TLorentzVector& p1, TLorentzVector& p2) {
    
    float dEta = p1.Eta() - p2.Eta();
    float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
    float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
    return dR;
}

float SusyEventAnalyzer::GetDiEmPt(TLorentzVector Pho1, TLorentzVector Pho2){
    float diempt = (Pho1 + Pho2).Pt();
    return diempt;
}

float SusyEventAnalyzer::GetDiJetPt(TLorentzVector Jet1, TLorentzVector Jet2){
  float dijetpt = (Jet1 + Jet2).Pt();
  return dijetpt;
}


float SusyEventAnalyzer::TransverseMass(TLorentzVector part, TVector2 met){
    float dPhi = std::fabs(TVector2::Phi_mpi_pi(part.Phi()-met.Phi()));
    float Mt2 = 2*part.Pt()*met.Mod()*(1-cos(dPhi));
    float Mt=sqrt(Mt2);
    return Mt;
}


float SusyEventAnalyzer::TransverseMassSquare3body(TLorentzVector lep1, TLorentzVector lep2, TVector2 met){
    TLorentzVector temp = lep1+lep2;
    float dPhi = std::fabs(TVector2::Phi_mpi_pi(temp.Phi()-met.Phi()));
    float Mt2 = 2*temp.Pt()*met.Mod()*(1-cos(dPhi));
    //float Mt=sqrt(Mt2);
    return Mt2;
}

float SusyEventAnalyzer::MT2(TLorentzVector P1, TLorentzVector P2, TLorentzVector P3){
    float mt2 = (P1 + P2 + P3).Mt2();
    return mt2;
}

float SusyEventAnalyzer::GetPhotonLessHt(float Ht, TLorentzVector pOne, TLorentzVector pTwo)
{
    float pLessHt = Ht - pOne.Et() - pTwo.Et();
    return pLessHt;
}



float SusyEventAnalyzer::InvariantMass(Float_t ggInvMass, TLorentzVector part, TVector2 met){
    
    float dPhi = std::fabs(TVector2::Phi_mpi_pi(part.Phi()-met.Phi()));
    float Mt2 = ggInvMass*ggInvMass + 2*part.Pt()*met.Mod()*(1-cos(dPhi));
    float InvMass=sqrt(Mt2);
    return InvMass;
    
}


#endif /* defined(____SusyEventAnalyzer__) */
