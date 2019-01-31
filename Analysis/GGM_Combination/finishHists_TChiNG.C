#include <unordered_map>
#include <sstream>
#include <fstream>


void finishHists_TChiNG(bool both = false){
    TFile *fileIn = TFile::Open("hist_TChiNg.root");
    TFile *fileOut;
    if(both){
      fileOut = TFile::Open("output_TChiNg_N2C1_C1C1.root","RECREATE");
    }
    else{
      fileOut = TFile::Open("output_TChiNg_C1C1.root","RECREATE");
    }
    fileOut->cd();

    float totLA = 0, totLG = 0, totHA = 0, totHG = 0;

    //Take the unweighted histograms and set 
    //values appropriately based on the total normalization 
    //and cross section in each bin. 
    //Output of this script is used in 

    //All events, before gg selection
    TH1F* h_GridAllEvents = (TH1F*) fileIn->Get("h_GridAllEvents");

    //----------------------------------  
    //Final prediction (average of pfMET and genMET)
    TH1F* h_MET100to115_Grid = (TH1F*) fileIn->Get("h_MET100to115_Grid");
    TH1F* h_MET115to130_Grid = (TH1F*) fileIn->Get("h_MET115to130_Grid");
    TH1F* h_MET130to150_Grid = (TH1F*) fileIn->Get("h_MET130to150_Grid");
    TH1F* h_MET150to185_Grid = (TH1F*) fileIn->Get("h_MET150to185_Grid");
    TH1F* h_MET185to250_Grid = (TH1F*) fileIn->Get("h_MET185to250_Grid");
    TH1F* h_MET250_Grid = (TH1F*) fileIn->Get("h_MET250_Grid");

    //Unweighted yields using pfMET
    TH1F* h_MET100to115_GridUnweighted = (TH1F*) fileIn->Get("h_MET100to115_GridUnweighted");
    TH1F* h_MET115to130_GridUnweighted = (TH1F*) fileIn->Get("h_MET115to130_GridUnweighted");
    TH1F* h_MET130to150_GridUnweighted = (TH1F*) fileIn->Get("h_MET130to150_GridUnweighted");
    TH1F* h_MET150to185_GridUnweighted = (TH1F*) fileIn->Get("h_MET150to185_GridUnweighted");
    TH1F* h_MET185to250_GridUnweighted = (TH1F*) fileIn->Get("h_MET185to250_GridUnweighted");
    TH1F* h_MET250_GridUnweighted = (TH1F*) fileIn->Get("h_MET250_GridUnweighted");

    //Unweighted yields using genMET
    TH1F* h_genMET100to115_GridUnweighted = (TH1F*) fileIn->Get("h_genMET100to115_GridUnweighted");
    TH1F* h_genMET115to130_GridUnweighted = (TH1F*) fileIn->Get("h_genMET115to130_GridUnweighted");
    TH1F* h_genMET130to150_GridUnweighted = (TH1F*) fileIn->Get("h_genMET130to150_GridUnweighted");
    TH1F* h_genMET150to185_GridUnweighted = (TH1F*) fileIn->Get("h_genMET150to185_GridUnweighted");
    TH1F* h_genMET185to250_GridUnweighted = (TH1F*) fileIn->Get("h_genMET185to250_GridUnweighted");
    TH1F* h_genMET250_GridUnweighted = (TH1F*) fileIn->Get("h_genMET250_GridUnweighted");

    //----------------------------------
    //All events, before gg selection, ISR weighted
    TH1F* h_GridAllEvents_ISRWeighted = (TH1F*) fileIn->Get("h_GridAllEvents_ISRWeighted");

    //ISR-weighted yields using genMET
    TH1F* h_genMET100to115_GridISRWeighted = (TH1F*) fileIn->Get("h_genMET100to115_GridISRWeighted");
    TH1F* h_genMET115to130_GridISRWeighted = (TH1F*) fileIn->Get("h_genMET115to130_GridISRWeighted");
    TH1F* h_genMET130to150_GridISRWeighted = (TH1F*) fileIn->Get("h_genMET130to150_GridISRWeighted");
    TH1F* h_genMET150to185_GridISRWeighted = (TH1F*) fileIn->Get("h_genMET150to185_GridISRWeighted");
    TH1F* h_genMET185to250_GridISRWeighted = (TH1F*) fileIn->Get("h_genMET185to250_GridISRWeighted");
    TH1F* h_genMET250_GridISRWeighted = (TH1F*) fileIn->Get("h_genMET250_GridISRWeighted");

    //ISR-weighted yields using pf MET
    TH1F* h_MET100to115_GridISRWeighted = (TH1F*) fileIn->Get("h_MET100to115_GridISRWeighted");
    TH1F* h_MET115to130_GridISRWeighted = (TH1F*) fileIn->Get("h_MET115to130_GridISRWeighted");
    TH1F* h_MET130to150_GridISRWeighted = (TH1F*) fileIn->Get("h_MET130to150_GridISRWeighted");
    TH1F* h_MET150to185_GridISRWeighted = (TH1F*) fileIn->Get("h_MET150to185_GridISRWeighted");
    TH1F* h_MET185to250_GridISRWeighted = (TH1F*) fileIn->Get("h_MET185to250_GridISRWeighted");
    TH1F* h_MET250_GridISRWeighted = (TH1F*) fileIn->Get("h_MET250_GridISRWeighted");

    //----------------------------------
    //Weighted yields using genMET
    TH1F* h_genMET100to115_Grid = (TH1F*) h_genMET100to115_GridUnweighted->Clone("h_genMET100to115_Grid");
    TH1F* h_genMET115to130_Grid = (TH1F*) h_genMET100to115_GridUnweighted->Clone("h_genMET115to130_Grid");
    TH1F* h_genMET130to150_Grid = (TH1F*) h_genMET100to115_GridUnweighted->Clone("h_genMET130to150_Grid");
    TH1F* h_genMET150to185_Grid = (TH1F*) h_genMET100to115_GridUnweighted->Clone("h_genMET150to185_Grid");
    TH1F* h_genMET185to250_Grid = (TH1F*) h_genMET100to115_GridUnweighted->Clone("h_genMET185to250_Grid");
    TH1F* h_genMET250_Grid = (TH1F*) h_genMET100to115_GridUnweighted->Clone("h_genMET250_Grid");

    h_genMET100to115_Grid->Reset();
    h_genMET115to130_Grid->Reset();
    h_genMET130to150_Grid->Reset();
    h_genMET150to185_Grid->Reset();
    h_genMET185to250_Grid->Reset();
    h_genMET250_Grid->Reset();

    //Weighted yields using genMET, after ISR weighting
    TH1F* h_genMET100to115_GridISR = (TH1F*) h_genMET100to115_GridUnweighted->Clone("h_genMET100to115_GridISR");
    TH1F* h_genMET115to130_GridISR = (TH1F*) h_genMET100to115_GridUnweighted->Clone("h_genMET115to130_GridISR");
    TH1F* h_genMET130to150_GridISR = (TH1F*) h_genMET100to115_GridUnweighted->Clone("h_genMET130to150_GridISR");
    TH1F* h_genMET150to185_GridISR = (TH1F*) h_genMET100to115_GridUnweighted->Clone("h_genMET150to185_GridISR");
    TH1F* h_genMET185to250_GridISR = (TH1F*) h_genMET100to115_GridUnweighted->Clone("h_genMET185to250_GridISR");
    TH1F* h_genMET250_GridISR = (TH1F*) h_genMET100to115_GridUnweighted->Clone("h_genMET250_GridISR");

    h_genMET100to115_GridISR->Reset();
    h_genMET115to130_GridISR->Reset();
    h_genMET130to150_GridISR->Reset();
    h_genMET150to185_GridISR->Reset();
    h_genMET185to250_GridISR->Reset();
    h_genMET250_GridISR->Reset();

    //---------------------------------- 
    //Weighted yields using pfMET
    TH1F* h_pfMET100to115_Grid = (TH1F*)  h_MET100to115_Grid->Clone("h_pfMET100to115_Grid");
    TH1F* h_pfMET115to130_Grid = (TH1F*)  h_MET100to115_Grid->Clone("h_pfMET115to130_Grid");
    TH1F* h_pfMET130to150_Grid = (TH1F*)  h_MET100to115_Grid->Clone("h_pfMET130to150_Grid");
    TH1F* h_pfMET150to185_Grid = (TH1F*)  h_MET100to115_Grid->Clone("h_pfMET150to185_Grid");
    TH1F* h_pfMET185to250_Grid = (TH1F*)  h_MET100to115_Grid->Clone("h_pfMET185to250_Grid");
    TH1F* h_pfMET250_Grid = (TH1F*)  h_MET100to115_Grid->Clone("h_pfMET250_Grid");

    h_pfMET100to115_Grid->Reset();
    h_pfMET115to130_Grid->Reset();
    h_pfMET130to150_Grid->Reset();
    h_pfMET150to185_Grid->Reset();
    h_pfMET185to250_Grid->Reset();
    h_pfMET250_Grid->Reset();

    //Weighted yields using pfMET, after ISR weighting 
    TH1F* h_pfMET100to115_GridISR = (TH1F*)  h_MET100to115_Grid->Clone("h_pfMET100to115_GridISR");
    TH1F* h_pfMET115to130_GridISR = (TH1F*)  h_MET100to115_Grid->Clone("h_pfMET115to130_GridISR");
    TH1F* h_pfMET130to150_GridISR = (TH1F*)  h_MET100to115_Grid->Clone("h_pfMET130to150_GridISR");
    TH1F* h_pfMET150to185_GridISR = (TH1F*)  h_MET100to115_Grid->Clone("h_pfMET150to185_GridISR");
    TH1F* h_pfMET185to250_GridISR = (TH1F*)  h_MET100to115_Grid->Clone("h_pfMET185to250_GridISR");
    TH1F* h_pfMET250_GridISR = (TH1F*)  h_MET100to115_Grid->Clone("h_pfMET250_GridISR");

    h_pfMET100to115_GridISR->Reset();
    h_pfMET115to130_GridISR->Reset();
    h_pfMET130to150_GridISR->Reset();
    h_pfMET150to185_GridISR->Reset();
    h_pfMET185to250_GridISR->Reset();
    h_pfMET250_GridISR->Reset();

    //----------------------------------  
    //Weighted yields using pfMET_JESUp and Down 
    TH1F* h_MET100to115_Grid_JESDown  = (TH1F*) fileIn->Get("h_MET100to115_Grid_JESDown");
    TH1F* h_MET115to130_Grid_JESDown  = (TH1F*) fileIn->Get("h_MET115to130_Grid_JESDown");
    TH1F* h_MET130to150_Grid_JESDown  = (TH1F*) fileIn->Get("h_MET130to150_Grid_JESDown");
    TH1F* h_MET150to185_Grid_JESDown  = (TH1F*) fileIn->Get("h_MET150to185_Grid_JESDown");
    TH1F* h_MET185to250_Grid_JESDown  = (TH1F*) fileIn->Get("h_MET185to250_Grid_JESDown");
    TH1F* h_MET250_Grid_JESDown  = (TH1F*) fileIn->Get("h_MET250_Grid_JESDown");


    TH1F* h_MET100to115_Grid_JESUp  = (TH1F*) fileIn->Get("h_MET100to115_Grid_JESUp");
    TH1F* h_MET115to130_Grid_JESUp  = (TH1F*) fileIn->Get("h_MET115to130_Grid_JESUp");
    TH1F* h_MET130to150_Grid_JESUp  = (TH1F*) fileIn->Get("h_MET130to150_Grid_JESUp");
    TH1F* h_MET150to185_Grid_JESUp  = (TH1F*) fileIn->Get("h_MET150to185_Grid_JESUp");
    TH1F* h_MET185to250_Grid_JESUp  = (TH1F*) fileIn->Get("h_MET185to250_Grid_JESUp");
    TH1F* h_MET250_Grid_JESUp  = (TH1F*) fileIn->Get("h_MET250_Grid_JESUp");

    //----------------------------------
    //Unweighted yields using pfMET_JESUp and Down
    TH1F* h_MET100to115_GridUnweighted_JESDown = (TH1F*) fileIn->Get("h_MET100to115_GridUnweighted_JESDown");
    TH1F* h_MET115to130_GridUnweighted_JESDown = (TH1F*) fileIn->Get("h_MET115to130_GridUnweighted_JESDown");
    TH1F* h_MET130to150_GridUnweighted_JESDown = (TH1F*) fileIn->Get("h_MET130to150_GridUnweighted_JESDown");
    TH1F* h_MET150to185_GridUnweighted_JESDown = (TH1F*) fileIn->Get("h_MET150to185_GridUnweighted_JESDown");
    TH1F* h_MET185to250_GridUnweighted_JESDown = (TH1F*) fileIn->Get("h_MET185to250_GridUnweighted_JESDown");
    TH1F* h_MET250_GridUnweighted_JESDown = (TH1F*) fileIn->Get("h_MET250_GridUnweighted_JESDown");

    TH1F* h_MET100to115_GridUnweighted_JESUp = (TH1F*) fileIn->Get("h_MET100to115_GridUnweighted_JESUp");
    TH1F* h_MET115to130_GridUnweighted_JESUp = (TH1F*) fileIn->Get("h_MET115to130_GridUnweighted_JESUp");
    TH1F* h_MET130to150_GridUnweighted_JESUp = (TH1F*) fileIn->Get("h_MET130to150_GridUnweighted_JESUp");
    TH1F* h_MET150to185_GridUnweighted_JESUp = (TH1F*) fileIn->Get("h_MET150to185_GridUnweighted_JESUp");
    TH1F* h_MET185to250_GridUnweighted_JESUp = (TH1F*) fileIn->Get("h_MET185to250_GridUnweighted_JESUp");
    TH1F* h_MET250_GridUnweighted_JESUp = (TH1F*) fileIn->Get("h_MET250_GridUnweighted_JESUp");

    //----------------------------------
    //Histograms to hold the JES uncertainties from pfMET_JESUp and Down
    TH1F* h_MET100to115_GridJESDownError =(TH1F*) fileIn->Get("h_MET100to115_GridJESDownError");
    TH1F* h_MET115to130_GridJESDownError =(TH1F*) fileIn->Get("h_MET115to130_GridJESDownError");
    TH1F* h_MET130to150_GridJESDownError =(TH1F*) fileIn->Get("h_MET130to150_GridJESDownError");
    TH1F* h_MET150to185_GridJESDownError =(TH1F*) fileIn->Get("h_MET150to185_GridJESDownError");
    TH1F* h_MET185to250_GridJESDownError =(TH1F*) fileIn->Get("h_MET185to250_GridJESDownError");
    TH1F* h_MET250_GridJESDownError =(TH1F*) fileIn->Get("h_MET250_GridJESDownError");

    TH1F* h_MET100to115_GridJESUpError =(TH1F*) fileIn->Get("h_MET100to115_GridJESUpError");
    TH1F* h_MET115to130_GridJESUpError =(TH1F*) fileIn->Get("h_MET115to130_GridJESUpError");
    TH1F* h_MET130to150_GridJESUpError =(TH1F*) fileIn->Get("h_MET130to150_GridJESUpError");
    TH1F* h_MET150to185_GridJESUpError =(TH1F*) fileIn->Get("h_MET150to185_GridJESUpError");
    TH1F* h_MET185to250_GridJESUpError =(TH1F*) fileIn->Get("h_MET185to250_GridJESUpError");
    TH1F* h_MET250_GridJESUpError =(TH1F*) fileIn->Get("h_MET250_GridJESUpError");

    h_MET100to115_GridJESDownError->Reset();
    h_MET115to130_GridJESDownError->Reset();
    h_MET130to150_GridJESDownError->Reset();
    h_MET150to185_GridJESDownError->Reset();
    h_MET185to250_GridJESDownError->Reset();
    h_MET250_GridJESDownError->Reset();

    h_MET100to115_GridJESUpError->Reset();
    h_MET115to130_GridJESUpError->Reset();
    h_MET130to150_GridJESUpError->Reset();
    h_MET150to185_GridJESUpError->Reset();
    h_MET185to250_GridJESUpError->Reset();
    h_MET250_GridJESUpError->Reset();

    //-------------------------------------
    //Histograms to hold the gen MET uncertainty
    TH1F* h_MET100to115_GridGenMetError =(TH1F*) h_MET100to115_Grid->Clone("h_MET100to115_GridGenMetError");
    TH1F* h_MET115to130_GridGenMetError =(TH1F*) h_MET100to115_Grid->Clone("h_MET115to130_GridGenMetError");
    TH1F* h_MET130to150_GridGenMetError =(TH1F*) h_MET100to115_Grid->Clone("h_MET130to150_GridGenMetError");
    TH1F* h_MET150to185_GridGenMetError =(TH1F*) h_MET100to115_Grid->Clone("h_MET150to185_GridGenMetError");
    TH1F* h_MET185to250_GridGenMetError =(TH1F*) h_MET100to115_Grid->Clone("h_MET185to250_GridGenMetError");
    TH1F* h_MET250_GridGenMetError =(TH1F*) h_MET100to115_Grid->Clone("h_MET250_GridGenMetError");

    h_MET100to115_GridGenMetError->Reset();
    h_MET115to130_GridGenMetError->Reset();
    h_MET130to150_GridGenMetError->Reset();
    h_MET150to185_GridGenMetError->Reset();
    h_MET185to250_GridGenMetError->Reset();
    h_MET250_GridGenMetError->Reset();

    //-------------------------------------
    //Histograms to hold the ISR uncertainty
    TH1F* h_MET100to115_GridISRError =(TH1F*) h_MET100to115_Grid->Clone("h_MET100to115_GridISRError");
    TH1F* h_MET115to130_GridISRError =(TH1F*) h_MET100to115_Grid->Clone("h_MET115to130_GridISRError");
    TH1F* h_MET130to150_GridISRError =(TH1F*) h_MET100to115_Grid->Clone("h_MET130to150_GridISRError");
    TH1F* h_MET150to185_GridISRError =(TH1F*) h_MET100to115_Grid->Clone("h_MET150to185_GridISRError");
    TH1F* h_MET185to250_GridISRError =(TH1F*) h_MET100to115_Grid->Clone("h_MET185to250_GridISRError");
    TH1F* h_MET250_GridISRError =(TH1F*) h_MET100to115_Grid->Clone("h_MET250_GridISRError");

    h_MET100to115_GridISRError->Reset();
    h_MET115to130_GridISRError->Reset();
    h_MET130to150_GridISRError->Reset();
    h_MET150to185_GridISRError->Reset();
    h_MET185to250_GridISRError->Reset();
    h_MET250_GridISRError->Reset();

    //------------------------------------- 
    //Unweighted number of events for low and high NVertex
    TH1F* h1_lowAll = (TH1F*) fileIn->Get("h_MET100to115_GridAll_LowNvtx");
    TH1F* h2_lowAll = (TH1F*) fileIn->Get("h_MET115to130_GridAll_LowNvtx");
    TH1F* h3_lowAll = (TH1F*) fileIn->Get("h_MET130to150_GridAll_LowNvtx");
    TH1F* h4_lowAll = (TH1F*) fileIn->Get("h_MET150to185_GridAll_LowNvtx");
    TH1F* h5_lowAll = (TH1F*) fileIn->Get("h_MET185to250_GridAll_LowNvtx");
    TH1F* h6_lowAll = (TH1F*) fileIn->Get("h_MET250_GridAll_LowNvtx");

    TH1F* h1_highAll = (TH1F*) fileIn->Get("h_MET100to115_GridAll_HighNvtx");
    TH1F* h2_highAll = (TH1F*) fileIn->Get("h_MET115to130_GridAll_HighNvtx");
    TH1F* h3_highAll = (TH1F*) fileIn->Get("h_MET130to150_GridAll_HighNvtx");
    TH1F* h4_highAll = (TH1F*) fileIn->Get("h_MET150to185_GridAll_HighNvtx");
    TH1F* h5_highAll = (TH1F*) fileIn->Get("h_MET185to250_GridAll_HighNvtx");
    TH1F* h6_highAll = (TH1F*) fileIn->Get("h_MET250_GridAll_HighNvtx");

    TH1F* h1_lowGG = (TH1F*) fileIn->Get("h_MET100to115_GridGG_LowNvtx");
    TH1F* h2_lowGG = (TH1F*) fileIn->Get("h_MET115to130_GridGG_LowNvtx");
    TH1F* h3_lowGG = (TH1F*) fileIn->Get("h_MET130to150_GridGG_LowNvtx");
    TH1F* h4_lowGG = (TH1F*) fileIn->Get("h_MET150to185_GridGG_LowNvtx");
    TH1F* h5_lowGG = (TH1F*) fileIn->Get("h_MET185to250_GridGG_LowNvtx");
    TH1F* h6_lowGG = (TH1F*) fileIn->Get("h_MET250_GridGG_LowNvtx");

    TH1F* h1_highGG = (TH1F*) fileIn->Get("h_MET100to115_GridGG_HighNvtx");
    TH1F* h2_highGG = (TH1F*) fileIn->Get("h_MET115to130_GridGG_HighNvtx");
    TH1F* h3_highGG = (TH1F*) fileIn->Get("h_MET130to150_GridGG_HighNvtx");
    TH1F* h4_highGG = (TH1F*) fileIn->Get("h_MET150to185_GridGG_HighNvtx");
    TH1F* h5_highGG = (TH1F*) fileIn->Get("h_MET185to250_GridGG_HighNvtx");
    TH1F* h6_highGG = (TH1F*) fileIn->Get("h_MET250_GridGG_HighNvtx");
    //----------------------------------
    //Histograms to hold R for the full grid
    TH1F* h_MET100to115_NVtxEffRatio =(TH1F*) h_MET100to115_Grid->Clone("h_MET100to115_NVtxEffRatio");
    TH1F* h_MET115to130_NVtxEffRatio =(TH1F*) h_MET100to115_Grid->Clone("h_MET115to130_NVtxEffRatio");
    TH1F* h_MET130to150_NVtxEffRatio =(TH1F*) h_MET100to115_Grid->Clone("h_MET130to150_NVtxEffRatio");
    TH1F* h_MET150to185_NVtxEffRatio =(TH1F*) h_MET100to115_Grid->Clone("h_MET150to185_NVtxEffRatio");
    TH1F* h_MET185to250_NVtxEffRatio =(TH1F*) h_MET100to115_Grid->Clone("h_MET185to250_NVtxEffRatio");
    TH1F* h_MET250_NVtxEffRatio =(TH1F*) h_MET100to115_Grid->Clone("h_MET250_NVtxEffRatio");

    h_MET100to115_NVtxEffRatio->Reset();
    h_MET115to130_NVtxEffRatio->Reset();
    h_MET130to150_NVtxEffRatio->Reset();
    h_MET150to185_NVtxEffRatio->Reset();
    h_MET185to250_NVtxEffRatio->Reset();
    h_MET250_NVtxEffRatio->Reset();

    TH1F* h_MET100to115_NVtxEffRatioErr =(TH1F*) h_MET100to115_Grid->Clone("h_MET100to115_NVtxEffRatioErr");
    TH1F* h_MET115to130_NVtxEffRatioErr =(TH1F*) h_MET100to115_Grid->Clone("h_MET115to130_NVtxEffRatioErr");
    TH1F* h_MET130to150_NVtxEffRatioErr =(TH1F*) h_MET100to115_Grid->Clone("h_MET130to150_NVtxEffRatioErr");
    TH1F* h_MET150to185_NVtxEffRatioErr =(TH1F*) h_MET100to115_Grid->Clone("h_MET150to185_NVtxEffRatioErr");
    TH1F* h_MET185to250_NVtxEffRatioErr =(TH1F*) h_MET100to115_Grid->Clone("h_MET185to250_NVtxEffRatioErr");
    TH1F* h_MET250_NVtxEffRatioErr =(TH1F*) h_MET100to115_Grid->Clone("h_MET250_NVtxEffRatioErr");

    h_MET100to115_NVtxEffRatioErr->Reset();
    h_MET115to130_NVtxEffRatioErr->Reset();
    h_MET130to150_NVtxEffRatioErr->Reset();
    h_MET150to185_NVtxEffRatioErr->Reset();
    h_MET185to250_NVtxEffRatioErr->Reset();
    h_MET250_NVtxEffRatioErr->Reset();

    //---------------------------------- 
    float lumi = 35.9218;

    ifstream file_xsec ("C1C1_xsec.txt");
    ifstream file_xsec2 ("N2C1_xsec.txt");

    unordered_map<int, float> xsec_map;
    unordered_map<int, float> xsec_map2;

    string line;
    string tempstr;
    stringstream ss;
    int mass;
    float xsec, uncert;
    if(file_xsec.is_open()){
      cout << "Opened file" << endl;
      while( file_xsec >> mass >> xsec){
	xsec_map[mass] = xsec;
	//      cout << "Filled " << mass << " with " << xsec << endl; 
      }
    }
    else{
      cout <<"File not open" << endl;
    }

    if(file_xsec2.is_open()){
      cout << "Opened file 2" << endl;
      while( file_xsec2 >> mass >> xsec){
        xsec_map2[mass] = xsec;
        cout << "Filled " << mass << " with " << xsec << endl;
      } 
    }
    else{
      cout <<"File 2 not open" << endl;
    }

    for(int mass =300; mass <= 1300; mass = mass + 25){
      cout << "Processing " << mass << ". " ;
      float xsec = 0;
      int i = h_MET100to115_GridUnweighted->FindFixBin(mass);
      try{
        xsec = xsec_map[mass];
      }
      catch(...){
        cout<<"Uh oh" << endl;
        break;
      }

      if(both){
        try{
          xsec += xsec_map2[mass];
        }
        catch(...){
          cout<<"Uh oh 2" << endl;
          break;
        }
      }

      float scaleByLumi = 0;
      if( h_GridAllEvents->GetBinContent(i) > 0 ){
	scaleByLumi = lumi*xsec/h_GridAllEvents->GetBinContent(i);
	cout << "Scale factor = " << scaleByLumi << ", xsec = " << xsec << endl;
      }
      else{
	cout << "No events for " << mass << " :(" << endl;
      }	  
      //Normalization factor C for ISR uncertainty
      float scaleByISRC = 0;
      if( h_GridAllEvents_ISRWeighted->GetBinContent(i) > 0 ){
        scaleByISRC = h_GridAllEvents->GetBinContent(i)/h_GridAllEvents_ISRWeighted->GetBinContent(i);
      }

      //Scale pf, gen, and JES MET distributions by lumi
      h_pfMET100to115_Grid->SetBinContent(i,h_MET100to115_GridUnweighted->GetBinContent(i)*scaleByLumi);
      h_pfMET115to130_Grid->SetBinContent(i,h_MET115to130_GridUnweighted->GetBinContent(i)*scaleByLumi);
      h_pfMET130to150_Grid->SetBinContent(i,h_MET130to150_GridUnweighted->GetBinContent(i)*scaleByLumi);
      h_pfMET150to185_Grid->SetBinContent(i,h_MET150to185_GridUnweighted->GetBinContent(i)*scaleByLumi);
      h_pfMET185to250_Grid->SetBinContent(i,h_MET185to250_GridUnweighted->GetBinContent(i)*scaleByLumi);
      h_pfMET250_Grid->SetBinContent(i,h_MET250_GridUnweighted->GetBinContent(i)*scaleByLumi);
            
      h_pfMET100to115_GridISR->SetBinContent(i,h_MET100to115_GridISRWeighted->GetBinContent(i)*scaleByLumi*scaleByISRC);
      h_pfMET115to130_GridISR->SetBinContent(i,h_MET115to130_GridISRWeighted->GetBinContent(i)*scaleByLumi*scaleByISRC);
      h_pfMET130to150_GridISR->SetBinContent(i,h_MET130to150_GridISRWeighted->GetBinContent(i)*scaleByLumi*scaleByISRC);
      h_pfMET150to185_GridISR->SetBinContent(i,h_MET150to185_GridISRWeighted->GetBinContent(i)*scaleByLumi*scaleByISRC);
      h_pfMET185to250_GridISR->SetBinContent(i,h_MET185to250_GridISRWeighted->GetBinContent(i)*scaleByLumi*scaleByISRC);
      h_pfMET250_GridISR->SetBinContent(i,h_MET250_GridISRWeighted->GetBinContent(i)*scaleByLumi*scaleByISRC);

      h_MET100to115_Grid_JESDown->SetBinContent(i,h_MET100to115_GridUnweighted_JESDown->GetBinContent(i)*scaleByLumi);
      h_MET115to130_Grid_JESDown->SetBinContent(i,h_MET115to130_GridUnweighted_JESDown->GetBinContent(i)*scaleByLumi);
      h_MET130to150_Grid_JESDown->SetBinContent(i,h_MET130to150_GridUnweighted_JESDown->GetBinContent(i)*scaleByLumi);
      h_MET150to185_Grid_JESDown->SetBinContent(i,h_MET150to185_GridUnweighted_JESDown->GetBinContent(i)*scaleByLumi);
      h_MET185to250_Grid_JESDown->SetBinContent(i,h_MET185to250_GridUnweighted_JESDown->GetBinContent(i)*scaleByLumi);
      h_MET250_Grid_JESDown->SetBinContent(i,h_MET250_GridUnweighted_JESDown->GetBinContent(i)*scaleByLumi);
      
      h_MET100to115_Grid_JESUp->SetBinContent(i,h_MET100to115_GridUnweighted_JESUp->GetBinContent(i)*scaleByLumi);
      h_MET115to130_Grid_JESUp->SetBinContent(i,h_MET115to130_GridUnweighted_JESUp->GetBinContent(i)*scaleByLumi);
      h_MET130to150_Grid_JESUp->SetBinContent(i,h_MET130to150_GridUnweighted_JESUp->GetBinContent(i)*scaleByLumi);
      h_MET150to185_Grid_JESUp->SetBinContent(i,h_MET150to185_GridUnweighted_JESUp->GetBinContent(i)*scaleByLumi);
      h_MET185to250_Grid_JESUp->SetBinContent(i,h_MET185to250_GridUnweighted_JESUp->GetBinContent(i)*scaleByLumi);
      h_MET250_Grid_JESUp->SetBinContent(i,h_MET250_GridUnweighted_JESUp->GetBinContent(i)*scaleByLumi);

      h_genMET100to115_Grid->SetBinContent(i,h_genMET100to115_GridUnweighted->GetBinContent(i)*scaleByLumi);
      h_genMET115to130_Grid->SetBinContent(i,h_genMET115to130_GridUnweighted->GetBinContent(i)*scaleByLumi);
      h_genMET130to150_Grid->SetBinContent(i,h_genMET130to150_GridUnweighted->GetBinContent(i)*scaleByLumi);
      h_genMET150to185_Grid->SetBinContent(i,h_genMET150to185_GridUnweighted->GetBinContent(i)*scaleByLumi);
      h_genMET185to250_Grid->SetBinContent(i,h_genMET185to250_GridUnweighted->GetBinContent(i)*scaleByLumi);
      h_genMET250_Grid->SetBinContent(i,h_genMET250_GridUnweighted->GetBinContent(i)*scaleByLumi);

      h_genMET100to115_GridISR->SetBinContent(i,h_genMET100to115_GridISRWeighted->GetBinContent(i)*scaleByLumi*scaleByISRC);
      h_genMET115to130_GridISR->SetBinContent(i,h_genMET115to130_GridISRWeighted->GetBinContent(i)*scaleByLumi*scaleByISRC);
      h_genMET130to150_GridISR->SetBinContent(i,h_genMET130to150_GridISRWeighted->GetBinContent(i)*scaleByLumi*scaleByISRC);
      h_genMET150to185_GridISR->SetBinContent(i,h_genMET150to185_GridISRWeighted->GetBinContent(i)*scaleByLumi*scaleByISRC);
      h_genMET185to250_GridISR->SetBinContent(i,h_genMET185to250_GridISRWeighted->GetBinContent(i)*scaleByLumi*scaleByISRC);
      h_genMET250_GridISR->SetBinContent(i,h_genMET250_GridISRWeighted->GetBinContent(i)*scaleByLumi*scaleByISRC);

      //Actual prediction is average of pf and gen MET after ISR weighting
      h_MET100to115_Grid->SetBinContent(i, ( h_pfMET100to115_GridISR->GetBinContent(i) + h_genMET100to115_GridISR->GetBinContent(i) ) / 2) ;
      h_MET115to130_Grid->SetBinContent(i, ( h_pfMET115to130_GridISR->GetBinContent(i) + h_genMET115to130_GridISR->GetBinContent(i) ) / 2) ;
      h_MET130to150_Grid->SetBinContent(i, ( h_pfMET130to150_GridISR->GetBinContent(i) + h_genMET130to150_GridISR->GetBinContent(i) ) / 2) ;
      h_MET150to185_Grid->SetBinContent(i, ( h_pfMET150to185_GridISR->GetBinContent(i) + h_genMET150to185_GridISR->GetBinContent(i) ) / 2) ;
      h_MET185to250_Grid->SetBinContent(i, ( h_pfMET185to250_GridISR->GetBinContent(i) + h_genMET185to250_GridISR->GetBinContent(i) ) / 2) ;
      h_MET250_Grid->SetBinContent(i, ( h_pfMET250_GridISR->GetBinContent(i) + h_genMET250_GridISR->GetBinContent(i) ) / 2) ;

      //ISR uncertainty is half the difference between weighted and unweighted prediction 
      float noISRPred, yesISRPred, ISRError; 
      noISRPred = (h_genMET100to115_Grid->GetBinContent(i) + h_pfMET100to115_Grid->GetBinContent(i) )/2;
      yesISRPred = h_MET100to115_Grid->GetBinContent(i); 
      if(yesISRPred > 0) ISRError = TMath::Abs(noISRPred-yesISRPred)/(yesISRPred);
      else ISRError = 0.0;
      h_MET100to115_GridISRError->SetBinContent(i,ISRError);

      noISRPred = (h_genMET115to130_Grid->GetBinContent(i) + h_pfMET115to130_Grid->GetBinContent(i) )/2;
      yesISRPred = h_MET115to130_Grid->GetBinContent(i);
      if(yesISRPred > 0) ISRError = TMath::Abs(noISRPred-yesISRPred)/(yesISRPred);
      else ISRError = 0.0;
      h_MET115to130_GridISRError->SetBinContent(i,ISRError);

      noISRPred = (h_genMET130to150_Grid->GetBinContent(i) + h_pfMET130to150_Grid->GetBinContent(i) )/2;
      yesISRPred = h_MET130to150_Grid->GetBinContent(i);
      if(yesISRPred > 0) ISRError = TMath::Abs(noISRPred-yesISRPred)/(yesISRPred);
      else ISRError = 0.0;
      h_MET130to150_GridISRError->SetBinContent(i,ISRError);

      noISRPred = (h_genMET150to185_Grid->GetBinContent(i) + h_pfMET150to185_Grid->GetBinContent(i) )/2;
      yesISRPred = h_MET150to185_Grid->GetBinContent(i);
      if(yesISRPred > 0) ISRError = TMath::Abs(noISRPred-yesISRPred)/(yesISRPred);
      else ISRError = 0.0;
      h_MET150to185_GridISRError->SetBinContent(i,ISRError);

      noISRPred = (h_genMET185to250_Grid->GetBinContent(i) + h_pfMET185to250_Grid->GetBinContent(i) )/2;
      yesISRPred = h_MET185to250_Grid->GetBinContent(i);
      if(yesISRPred > 0) ISRError = TMath::Abs(noISRPred-yesISRPred)/(yesISRPred);
      else ISRError = 0.0;
      h_MET185to250_GridISRError->SetBinContent(i,ISRError);

      noISRPred = (h_genMET250_Grid->GetBinContent(i) + h_pfMET250_Grid->GetBinContent(i) )/2;
      yesISRPred = h_MET250_Grid->GetBinContent(i);
      if(yesISRPred > 0) ISRError = TMath::Abs(noISRPred-yesISRPred)/(yesISRPred);
      else ISRError = 0.0;
      h_MET250_GridISRError->SetBinContent(i,ISRError);

      //Gen MET uncertainty is half the difference between pf and gen MET
      //So the percent uncertainty is [abs(gen-pf) / 2 ] / [ (gen+pf) / 2 ] = (gen-pf) / (gen+pf)
      if( h_MET100to115_Grid->GetBinContent(i) > 0){
	h_MET100to115_GridGenMetError->
	  SetBinContent(i,  abs(h_genMET100to115_Grid->GetBinContent(i) - h_pfMET100to115_Grid->GetBinContent(i)) / 
			(h_genMET100to115_Grid->GetBinContent(i) + h_pfMET100to115_Grid->GetBinContent(i)) );
      }
      else{
	h_MET100to115_GridGenMetError->SetBinContent(i,0.0);
      }
      if(h_MET115to130_Grid->GetBinContent(i) > 0){
	h_MET115to130_GridGenMetError->
	  SetBinContent(i,  abs(h_genMET115to130_Grid->GetBinContent(i) - h_pfMET115to130_Grid->GetBinContent(i)) /
			(h_genMET115to130_Grid->GetBinContent(i) + h_pfMET115to130_Grid->GetBinContent(i)) );
      }
      else{
	h_MET115to130_GridGenMetError->SetBinContent(i,0.0);
      }
      if(h_MET130to150_Grid->GetBinContent(i) > 0){
	h_MET130to150_GridGenMetError->
	  SetBinContent(i,  abs(h_genMET130to150_Grid->GetBinContent(i) - h_pfMET130to150_Grid->GetBinContent(i)) /
			(h_genMET130to150_Grid->GetBinContent(i) + h_pfMET130to150_Grid->GetBinContent(i)) );
      }
      else{
	h_MET130to150_GridGenMetError->SetBinContent(i,0.0);
      }
      if(h_MET150to185_Grid->GetBinContent(i) > 0){
	h_MET150to185_GridGenMetError->
	  SetBinContent(i,  abs(h_genMET150to185_Grid->GetBinContent(i) - h_pfMET150to185_Grid->GetBinContent(i)) /
			(h_genMET150to185_Grid->GetBinContent(i) + h_pfMET150to185_Grid->GetBinContent(i)) );
      }
      else{
	h_MET150to185_GridGenMetError->SetBinContent(i,0.0);
      }
      if(h_MET185to250_Grid->GetBinContent(i) > 0){
	h_MET185to250_GridGenMetError->
	  SetBinContent(i,  abs(h_genMET185to250_Grid->GetBinContent(i) - h_pfMET185to250_Grid->GetBinContent(i)) /
			(h_genMET185to250_Grid->GetBinContent(i) + h_pfMET185to250_Grid->GetBinContent(i)) );
      }
      else{
	h_MET185to250_GridGenMetError->SetBinContent(i,0.0);
      }
      if(h_MET250_Grid->GetBinContent(i) > 0){
	h_MET250_GridGenMetError->
	  SetBinContent(i,  abs(h_genMET250_Grid->GetBinContent(i) - h_pfMET250_Grid->GetBinContent(i)) /
			(h_genMET250_Grid->GetBinContent(i) + h_pfMET250_Grid->GetBinContent(i)) );
      }
      else{
	h_MET250_GridGenMetError->SetBinContent(i,0.0);
      }
      
      //JES uncertainty is the difference between the pf uncertainty and the JES Met distributions
      //These histograms store the percent difference
      if(h_pfMET100to115_Grid->GetBinContent(i) > 0) {	    
	h_MET100to115_GridJESUpError->Fill(mass, 
					   abs( h_MET100to115_Grid_JESUp->GetBinContent(i) - 
						h_pfMET100to115_Grid->GetBinContent(i)) / 
					   h_pfMET100to115_Grid->GetBinContent(i) );
	h_MET100to115_GridJESDownError->Fill(mass, 
					     abs( h_MET100to115_Grid_JESDown->GetBinContent(i) - 
						  h_pfMET100to115_Grid->GetBinContent(i)) / 
					     h_pfMET100to115_Grid->GetBinContent(i) );
      }
      else{
	h_MET100to115_GridJESUpError->Fill(mass,0.0);
	h_MET100to115_GridJESDownError->Fill(mass,0.0);
      }
      
      if(h_pfMET115to130_Grid->GetBinContent(i) > 0) {
	h_MET115to130_GridJESUpError->Fill(mass,
					   abs( h_MET115to130_Grid_JESUp->GetBinContent(i) -
						h_pfMET115to130_Grid->GetBinContent(i)) /
					   h_pfMET115to130_Grid->GetBinContent(i) );
	h_MET115to130_GridJESDownError->Fill(mass,
					     abs( h_MET115to130_Grid_JESDown->GetBinContent(i) -
						  h_pfMET115to130_Grid->GetBinContent(i)) /
					     h_pfMET115to130_Grid->GetBinContent(i) );
      }
      else{
	h_MET115to130_GridJESUpError->Fill(mass,0.0);
	h_MET115to130_GridJESDownError->Fill(mass,0.0);
      }
      
      
      if(h_pfMET130to150_Grid->GetBinContent(i) > 0) {
	h_MET130to150_GridJESUpError->Fill(mass,
					   abs( h_MET130to150_Grid_JESUp->GetBinContent(i) -
						h_pfMET130to150_Grid->GetBinContent(i)) /
					   h_pfMET130to150_Grid->GetBinContent(i) );
	h_MET130to150_GridJESDownError->Fill(mass,
					     abs( h_MET130to150_Grid_JESDown->GetBinContent(i) -
						  h_pfMET130to150_Grid->GetBinContent(i)) /
					     h_pfMET130to150_Grid->GetBinContent(i) );
      }
      else{
	h_MET130to150_GridJESUpError->Fill(mass,0.0);
	h_MET130to150_GridJESDownError->Fill(mass,0.0);
      }
      
      
      if(h_pfMET150to185_Grid->GetBinContent(i) > 0) {
	h_MET150to185_GridJESUpError->Fill(mass,
					   abs( h_MET150to185_Grid_JESUp->GetBinContent(i) -
						h_pfMET150to185_Grid->GetBinContent(i)) /
					   h_pfMET150to185_Grid->GetBinContent(i) );
	h_MET150to185_GridJESDownError->Fill(mass,
					     abs( h_MET150to185_Grid_JESDown->GetBinContent(i) -
						  h_pfMET150to185_Grid->GetBinContent(i)) /
					     h_pfMET150to185_Grid->GetBinContent(i) );
      }
      else{
	h_MET150to185_GridJESUpError->Fill(mass,0.0);
	h_MET150to185_GridJESDownError->Fill(mass,0.0);
      }
      
      if(h_pfMET185to250_Grid->GetBinContent(i) > 0) {
	h_MET185to250_GridJESUpError->Fill(mass,
					   abs( h_MET185to250_Grid_JESUp->GetBinContent(i) -
						h_pfMET185to250_Grid->GetBinContent(i)) /
					   h_pfMET185to250_Grid->GetBinContent(i) );
	h_MET185to250_GridJESDownError->Fill(mass,
					     abs( h_MET185to250_Grid_JESDown->GetBinContent(i) -
						  h_pfMET185to250_Grid->GetBinContent(i)) /
					     h_pfMET185to250_Grid->GetBinContent(i) );
      }
      else{
	h_MET185to250_GridJESUpError->Fill(mass,0.0);
	h_MET185to250_GridJESDownError->Fill(mass,0.0);
      }
      
      if(h_pfMET250_Grid->GetBinContent(i) > 0) {
	h_MET250_GridJESUpError->Fill(mass,
				      abs( h_MET250_Grid_JESUp->GetBinContent(i) -
					   h_pfMET250_Grid->GetBinContent(i)) /
				      h_pfMET250_Grid->GetBinContent(i) );
	h_MET250_GridJESDownError->Fill(mass,
					abs( h_MET250_Grid_JESDown->GetBinContent(i) -
					     h_pfMET250_Grid->GetBinContent(i)) /
					h_pfMET250_Grid->GetBinContent(i) );
      }
      else{
	h_MET250_GridJESUpError->Fill(mass,0.0);
	h_MET250_GridJESDownError->Fill(mass,0.0);
      }

      //Calculate r = acceptance for high nvtx / acceptance for low nvtx
      //Bin 1
      float lA =  h1_lowAll->GetBinContent(i);
      float lG =  h1_lowGG->GetBinContent(i);
      float hA = h1_highAll->GetBinContent(i);
      float hG = h1_highGG->GetBinContent(i);

      float r = -1;
      float r_err = -1;
      if( lA > 0 && lG > 0 && hA > 0 && hG > 0){
	r = (hG/hA) / (lG/lA);
	r_err = r*TMath::Sqrt( (1/lA) + (1/lG) + (1/hA) + (1/hG) );
      }

      h_MET100to115_NVtxEffRatio->SetBinContent(i, r);
      h_MET100to115_NVtxEffRatioErr->SetBinContent(i, r_err);


      //Bin 2
      lA = h2_lowAll->GetBinContent(i);
      lG = h2_lowGG->GetBinContent(i);
      hA = h2_highAll->GetBinContent(i);
      hG = h2_highGG->GetBinContent(i);

      r = -1;
      r_err = -1;
      if( lA > 0 && lG > 0 && hA > 0 && hG > 0){
	r = (hG/hA) / (lG/lA);
	r_err = r*TMath::Sqrt( (1/lA) + (1/lG) + (1/hA) + (1/hG) );
      }

      h_MET115to130_NVtxEffRatio->SetBinContent(i, r);
      h_MET115to130_NVtxEffRatioErr->SetBinContent(i, r_err);

      totLA += lA;
      totLG += lG;
      totHA += hA;
      totHG += hG;

      //Bin 3
      lA = h3_lowAll->GetBinContent(i);
      lG = h3_lowGG->GetBinContent(i);
      hA = h3_highAll->GetBinContent(i);
      hG = h3_highGG->GetBinContent(i);

      r = -1;
      r_err = -1;
      if( lA > 0 && lG > 0 && hA > 0 && hG > 0){
	r = (hG/hA) / (lG/lA);
	r_err = r*TMath::Sqrt( (1/lA) + (1/lG) + (1/hA) + (1/hG) );
      }

      h_MET130to150_NVtxEffRatio->SetBinContent(i, r);
      h_MET130to150_NVtxEffRatioErr->SetBinContent(i, r_err);

      totLA += lA;
      totLG += lG;
      totHA += hA;
      totHG += hG;

      //Bin 4
      lA = h4_lowAll->GetBinContent(i);
      lG = h4_lowGG->GetBinContent(i);
      hA = h4_highAll->GetBinContent(i);
      hG = h4_highGG->GetBinContent(i);

      r = -1;
      r_err = -1;
      if( lA > 0 && lG > 0 && hA > 0 && hG > 0){
	r = (hG/hA) / (lG/lA);
	r_err = r*TMath::Sqrt( (1/lA) + (1/lG) + (1/hA) + (1/hG) );
      }

      h_MET150to185_NVtxEffRatio->SetBinContent(i, r);
      h_MET150to185_NVtxEffRatioErr->SetBinContent(i, r_err);

      totLA += lA;
      totLG += lG;
      totHA += hA;
      totHG += hG;

      //Bin 5
      lA = h5_lowAll->GetBinContent(i);
      lG = h5_lowGG->GetBinContent(i);
      hA = h5_highAll->GetBinContent(i);
      hG = h5_highGG->GetBinContent(i);

      r = -1;
      r_err = -1;
      if( lA > 0 && lG > 0 && hA > 0 && hG > 0){
	r = (hG/hA) / (lG/lA);
	r_err = r*TMath::Sqrt( (1/lA) + (1/lG) + (1/hA) + (1/hG) );
      }

      h_MET185to250_NVtxEffRatio->SetBinContent(i, r);
      h_MET185to250_NVtxEffRatioErr->SetBinContent(i, r_err);

      totLA += lA;
      totLG += lG;
      totHA += hA;
      totHG += hG;

      //Bin 6
      lA = h6_lowAll->GetBinContent(i);
      lG = h6_lowGG->GetBinContent(i);
      hA = h6_highAll->GetBinContent(i);
      hG = h6_highGG->GetBinContent(i);

      r = -1;
      r_err = -1;
      if( lA > 0 && lG > 0 && hA > 0 && hG > 0){
	r = (hG/hA) / (lG/lA);
	r_err = r*TMath::Sqrt( (1/lA) + (1/lG) + (1/hA) + (1/hG) );
      }

      h_MET250_NVtxEffRatio->SetBinContent(i, r);
      h_MET250_NVtxEffRatioErr->SetBinContent(i, r_err);

      totLA += lA;
      totLG += lG;
      totHA += hA;
      totHG += hG;

    }


    cout << "totLA = " << totLA << endl;
    cout << "totLG = " << totLG << endl;
    cout << "totHA = " << totHA << endl;
    cout << "totHG = " << totHG << endl;
    float rTot = (totHG/totHA)/(totLG/totLA);
    float rErr = rTot*TMath::Sqrt( (1/totHA) + (1/totHG) + (1/totLA) + (1/totLG) );

    cout << "rTot = " << rTot << endl;
    cout << "rErr = " << rErr << endl;

    h_MET100to115_Grid->Write("h_MET100to115_Grid");
    h_MET115to130_Grid->Write("h_MET115to130_Grid");
    h_MET130to150_Grid->Write("h_MET130to150_Grid");
    h_MET150to185_Grid->Write("h_MET150to185_Grid");
    h_MET185to250_Grid->Write("h_MET185to250_Grid");
    h_MET250_Grid->Write("h_MET250_Grid");

    h_MET100to115_NVtxEffRatio->Write("h_MET100to115_NVtxEffRatio");
    h_MET115to130_NVtxEffRatio->Write("h_MET115to130_NVtxEffRatio");
    h_MET130to150_NVtxEffRatio->Write("h_MET130to150_NVtxEffRatio");
    h_MET150to185_NVtxEffRatio->Write("h_MET150to185_NVtxEffRatio");
    h_MET185to250_NVtxEffRatio->Write("h_MET185to250_NVtxEffRatio");
    h_MET250_NVtxEffRatio->Write("h_MET250_NVtxEffRatio");

    h_MET100to115_NVtxEffRatioErr->Write("h_MET100to115_NVtxEffRatioErr");
    h_MET115to130_NVtxEffRatioErr->Write("h_MET115to130_NVtxEffRatioErr");
    h_MET130to150_NVtxEffRatioErr->Write("h_MET130to150_NVtxEffRatioErr");
    h_MET150to185_NVtxEffRatioErr->Write("h_MET150to185_NVtxEffRatioErr");
    h_MET185to250_NVtxEffRatioErr->Write("h_MET185to250_NVtxEffRatioErr");
    h_MET250_NVtxEffRatioErr->Write("h_MET250_NVtxEffRatioErr");

    h_MET100to115_GridUnweighted->Write("h_MET100to115_GridUnweighted");
    h_MET115to130_GridUnweighted->Write("h_MET115to130_GridUnweighted");
    h_MET130to150_GridUnweighted->Write("h_MET130to150_GridUnweighted");
    h_MET150to185_GridUnweighted->Write("h_MET150to185_GridUnweighted");
    h_MET185to250_GridUnweighted->Write("h_MET185to250_GridUnweighted");
    h_MET250_GridUnweighted->Write("h_MET250_GridUnweighted");

    h_MET100to115_GridGenMetError->Write("h_MET100to115_GridGenMetError");
    h_MET115to130_GridGenMetError->Write("h_MET115to130_GridGenMetError");
    h_MET130to150_GridGenMetError->Write("h_MET130to150_GridGenMetError");
    h_MET150to185_GridGenMetError->Write("h_MET150to185_GridGenMetError");
    h_MET185to250_GridGenMetError->Write("h_MET185to250_GridGenMetError");
    h_MET250_GridGenMetError->Write("h_MET250_GridGenMetError");

    h_MET100to115_GridISRError->Write("h_MET100to115_GridISRError");
    h_MET115to130_GridISRError->Write("h_MET115to130_GridISRError");
    h_MET130to150_GridISRError->Write("h_MET130to150_GridISRError");
    h_MET150to185_GridISRError->Write("h_MET150to185_GridISRError");
    h_MET185to250_GridISRError->Write("h_MET185to250_GridISRError");
    h_MET250_GridISRError->Write("h_MET250_GridISRError");

    h_MET100to115_GridJESDownError->Write("h_MET100to115_GridJESDownError");
    h_MET115to130_GridJESDownError->Write("h_MET115to130_GridJESDownError");
    h_MET130to150_GridJESDownError->Write("h_MET130to150_GridJESDownError");
    h_MET150to185_GridJESDownError->Write("h_MET150to185_GridJESDownError");
    h_MET185to250_GridJESDownError->Write("h_MET185to250_GridJESDownError");
    h_MET250_GridJESDownError->Write("h_MET250_GridJESDownError");

    h_MET100to115_GridJESUpError->Write("h_MET100to115_GridJESUpError");
    h_MET115to130_GridJESUpError->Write("h_MET115to130_GridJESUpError");
    h_MET130to150_GridJESUpError->Write("h_MET130to150_GridJESUpError");
    h_MET150to185_GridJESUpError->Write("h_MET150to185_GridJESUpError");
    h_MET185to250_GridJESUpError->Write("h_MET185to250_GridJESUpError");
    h_MET250_GridJESUpError->Write("h_MET250_GridJESUpError");

    fileOut->Write();
    
    return;
}


int key(int i, int j) {
    return (int) i*10000+j;
}
