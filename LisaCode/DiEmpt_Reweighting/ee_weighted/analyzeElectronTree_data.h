//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 17 17:11:53 2016 by ROOT version 5.34/36
// from TTree eetree/Event data
// found on file: ../../../data2016_HLT15_analysed.root
//////////////////////////////////////////////////////////

#ifndef analyzeElectronTree_data_h
#define analyzeElectronTree_data_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class analyzeElectronTree_data : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           Run;
   Int_t           Event;
   Int_t           Lumis;
   Float_t         Rho;
   Float_t         PUtrue;
   Float_t         PUweight;
   Int_t           NJet;
   Int_t           MET_Filters;
   Float_t         MET;
   Float_t         MET_Phi;
   Float_t         MET_sumEt;
   Float_t         MET_mEtSig;
   Float_t         MET_Sig;
   Float_t         MT2;
   Int_t           NEle;
   Float_t         Dipho_Pt;
   Float_t         Dipho_Mass;
   Float_t         Ele_Pt[4];   //[NEle]
   Float_t         Ele_Eta[4];   //[NEle]
   Float_t         Ele_Phi[4];   //[NEle]
   Float_t         Ele_SigmaIEtaIEta[4];   //[NEle]
   Float_t         Ele_SigmaIEtaIPhi[4];   //[NEle]
   Float_t         Ele_SigmaIPhiIPhi[4];   //[NEle]
   Float_t         HT;
   Float_t         ST;
   Float_t         MHT;
   Float_t         MHTPhi;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_Lumis;   //!
   TBranch        *b_Rho;   //!
   TBranch        *b_PUtrue;   //!
   TBranch        *b_PUweight;   //!
   TBranch        *b_NJet;   //!
   TBranch        *b_MET_Filters;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_MET_Phi;   //!
   TBranch        *b_MET_sumEt;   //!
   TBranch        *b_MET_mEtSig;   //!
   TBranch        *b_MET_Sig;   //!
   TBranch        *b_MT2;   //!
   TBranch        *b_NEle;   //!
   TBranch        *b_Dipho_Pt;   //!
   TBranch        *b_Dipho_Mass;   //!
   TBranch        *b_Ele_Pt;   //!
   TBranch        *b_Ele_Eta;   //!
   TBranch        *b_Ele_Phi;   //!
   TBranch        *b_Ele_SigmaIEtaIEta;   //!
   TBranch        *b_Ele_SigmaIEtaIPhi;   //!
   TBranch        *b_Ele_SigmaIPhiIPhi;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_ST;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_MHTPhi;   //!

   analyzeElectronTree_data(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~analyzeElectronTree_data() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
   virtual void    GetAllEntries(Long64_t entry);
   


   ClassDef(analyzeElectronTree_data,0);
};

#endif

#ifdef analyzeElectronTree_data_cxx
void analyzeElectronTree_data::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("Lumis", &Lumis, &b_Lumis);
   fChain->SetBranchAddress("Rho", &Rho, &b_Rho);
   fChain->SetBranchAddress("PUtrue", &PUtrue, &b_PUtrue);
   fChain->SetBranchAddress("PUweight", &PUweight, &b_PUweight);
   fChain->SetBranchAddress("NJet", &NJet, &b_NJet);
   fChain->SetBranchAddress("MET_Filters", &MET_Filters, &b_MET_Filters);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("MET_Phi", &MET_Phi, &b_MET_Phi);
   fChain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   fChain->SetBranchAddress("MET_mEtSig", &MET_mEtSig, &b_MET_mEtSig);
   fChain->SetBranchAddress("MET_Sig", &MET_Sig, &b_MET_Sig);
   fChain->SetBranchAddress("MT2", &MT2, &b_MT2);
   fChain->SetBranchAddress("NEle", &NEle, &b_NEle);
   fChain->SetBranchAddress("Dipho_Pt", &Dipho_Pt, &b_Dipho_Pt);
   fChain->SetBranchAddress("Dipho_Mass", &Dipho_Mass, &b_Dipho_Mass);
   fChain->SetBranchAddress("Ele_Pt", Ele_Pt, &b_Ele_Pt);
   fChain->SetBranchAddress("Ele_Eta", Ele_Eta, &b_Ele_Eta);
   fChain->SetBranchAddress("Ele_Phi", Ele_Phi, &b_Ele_Phi);
   fChain->SetBranchAddress("Ele_SigmaIEtaIEta", Ele_SigmaIEtaIEta, &b_Ele_SigmaIEtaIEta);
   fChain->SetBranchAddress("Ele_SigmaIEtaIPhi", Ele_SigmaIEtaIPhi, &b_Ele_SigmaIEtaIPhi);
   fChain->SetBranchAddress("Ele_SigmaIPhiIPhi", Ele_SigmaIPhiIPhi, &b_Ele_SigmaIPhiIPhi);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("ST", &ST, &b_ST);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("MHTPhi", &MHTPhi, &b_MHTPhi);
}

Bool_t analyzeElectronTree_data::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef analyzeElectronTree_data_cxx
