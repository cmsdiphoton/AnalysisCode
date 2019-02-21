//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 13 18:02:08 2016 by ROOT version 5.34/36
// from TTree fftree/Event data
// found on file: ../../data2016_HLT14_analysed.root
//////////////////////////////////////////////////////////

#ifndef analyzeFakeTree_data_h
#define analyzeFakeTree_data_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class analyzeFakeTree_data : public TSelector {
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
   Float_t         HT;
   Float_t         ST;
   Float_t         MHT;
   Float_t         MHTPhi;
   Int_t           NFake;
   Float_t         Dipho_Pt;
   Float_t         Dipho_Mass;
   Float_t         Pho_E[3];   //[NFake]
   Float_t         Pho_Et[3];   //[NFake]
   Float_t         Pho_Eta[3];   //[NFake]
   Float_t         Pho_Phi[3];   //[NFake]
   Int_t           Pho_hasPixelSeed[3];   //[NFake]
   Int_t           Pho_EleVeto[3];   //[NFake]
   Float_t         Pho_R9[3];   //[NFake]
   Float_t         Pho_HoverE[3];   //[NFake]
   Float_t         Pho_SigmaIEtaIEta[3];   //[NFake]
   Float_t         Pho_SigmaIEtaIPhi[3];   //[NFake]
   Float_t         Pho_SigmaIPhiIPhi[3];   //[NFake]
   Float_t         Pho_IDMVA[3];   //[NFake]
   Int_t           Pho_IDbit[3];   //[NFake]

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
   TBranch        *b_HT;   //!
   TBranch        *b_ST;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_MHTPhi;   //!
   TBranch        *b_NFake;   //!
   TBranch        *b_Dipho_Pt;   //!
   TBranch        *b_Dipho_Mass;   //!
   TBranch        *b_Pho_E;   //!
   TBranch        *b_Pho_Et;   //!
   TBranch        *b_Pho_Eta;   //!
   TBranch        *b_Pho_Phi;   //!
   TBranch        *b_Pho_hasPixelSeed;   //!
   TBranch        *b_Pho_EleVeto;   //!
   TBranch        *b_Pho_R9;   //!
   TBranch        *b_Pho_HoverE;   //!
   TBranch        *b_Pho_SigmaIEtaIEta;   //!
   TBranch        *b_Pho_SigmaIEtaIPhi;   //!
   TBranch        *b_Pho_SigmaIPhiIPhi;   //!
   TBranch        *b_Pho_IDMVA;   //!
   TBranch        *b_Pho_IDbit;   //!

   analyzeFakeTree_data(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~analyzeFakeTree_data() { }
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
   
   
   
   ClassDef(analyzeFakeTree_data,0);
};

#endif

#ifdef analyzeFakeTree_data_cxx
void analyzeFakeTree_data::Init(TTree *tree)
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
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("ST", &ST, &b_ST);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("MHTPhi", &MHTPhi, &b_MHTPhi);
   fChain->SetBranchAddress("NFake", &NFake, &b_NFake);
   fChain->SetBranchAddress("Dipho_Pt", &Dipho_Pt, &b_Dipho_Pt);
   fChain->SetBranchAddress("Dipho_Mass", &Dipho_Mass, &b_Dipho_Mass);
   fChain->SetBranchAddress("Pho_E", Pho_E, &b_Pho_E);
   fChain->SetBranchAddress("Pho_Et", Pho_Et, &b_Pho_Et);
   fChain->SetBranchAddress("Pho_Eta", Pho_Eta, &b_Pho_Eta);
   fChain->SetBranchAddress("Pho_Phi", Pho_Phi, &b_Pho_Phi);
   fChain->SetBranchAddress("Pho_hasPixelSeed", Pho_hasPixelSeed, &b_Pho_hasPixelSeed);
   fChain->SetBranchAddress("Pho_EleVeto", Pho_EleVeto, &b_Pho_EleVeto);
   fChain->SetBranchAddress("Pho_R9", Pho_R9, &b_Pho_R9);
   fChain->SetBranchAddress("Pho_HoverE", Pho_HoverE, &b_Pho_HoverE);
   fChain->SetBranchAddress("Pho_SigmaIEtaIEta", Pho_SigmaIEtaIEta, &b_Pho_SigmaIEtaIEta);
   fChain->SetBranchAddress("Pho_SigmaIEtaIPhi", Pho_SigmaIEtaIPhi, &b_Pho_SigmaIEtaIPhi);
   fChain->SetBranchAddress("Pho_SigmaIPhiIPhi", Pho_SigmaIPhiIPhi, &b_Pho_SigmaIPhiIPhi);
   fChain->SetBranchAddress("Pho_IDMVA", Pho_IDMVA, &b_Pho_IDMVA);
   fChain->SetBranchAddress("Pho_IDbit", Pho_IDbit, &b_Pho_IDbit);
}

Bool_t analyzeFakeTree_data::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef analyzeFakeTree_data_cxx
