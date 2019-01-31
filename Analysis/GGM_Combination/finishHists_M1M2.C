#include <unordered_map>
#include <sstream>
#include <fstream>

int key(int i, int j);

void finishHists_M1M2(){
    TFile *fileIn = TFile::Open("hist_GGM_M1M2.root");

    TFile *fileOut = TFile::Open("output_M1M2.root","RECREATE");
    fileOut->cd();

    float totLA = 0, totLG = 0, totHA = 0, totHG = 0;

    //Take the unweighted histograms and set 
    //values appropriately based on the total normalization 
    //and cross section in each bin. 
    //Output of this script is used in 

    //All events, before gg selection
    TH2F* h_GridAllEvents = (TH2F*) fileIn->Get("h_GridAllEvents");

    //----------------------------------  
    //Final prediction (average of pfMET and genMET)
    TH2F* h_MET100to115_Grid = (TH2F*) fileIn->Get("h_MET100to115_Grid");
    TH2F* h_MET115to130_Grid = (TH2F*) fileIn->Get("h_MET115to130_Grid");
    TH2F* h_MET130to150_Grid = (TH2F*) fileIn->Get("h_MET130to150_Grid");
    TH2F* h_MET150to185_Grid = (TH2F*) fileIn->Get("h_MET150to185_Grid");
    TH2F* h_MET185to250_Grid = (TH2F*) fileIn->Get("h_MET185to250_Grid");
    TH2F* h_MET250_Grid = (TH2F*) fileIn->Get("h_MET250_Grid");

    //Unweighted yields using pfMET
    TH2F* h_MET100to115_GridUnweighted = (TH2F*) fileIn->Get("h_MET100to115_GridUnweighted");
    TH2F* h_MET115to130_GridUnweighted = (TH2F*) fileIn->Get("h_MET115to130_GridUnweighted");
    TH2F* h_MET130to150_GridUnweighted = (TH2F*) fileIn->Get("h_MET130to150_GridUnweighted");
    TH2F* h_MET150to185_GridUnweighted = (TH2F*) fileIn->Get("h_MET150to185_GridUnweighted");
    TH2F* h_MET185to250_GridUnweighted = (TH2F*) fileIn->Get("h_MET185to250_GridUnweighted");
    TH2F* h_MET250_GridUnweighted = (TH2F*) fileIn->Get("h_MET250_GridUnweighted");

    //Unweighted yields using genMET
    TH2F* h_genMET100to115_GridUnweighted = (TH2F*) fileIn->Get("h_genMET100to115_GridUnweighted");
    TH2F* h_genMET115to130_GridUnweighted = (TH2F*) fileIn->Get("h_genMET115to130_GridUnweighted");
    TH2F* h_genMET130to150_GridUnweighted = (TH2F*) fileIn->Get("h_genMET130to150_GridUnweighted");
    TH2F* h_genMET150to185_GridUnweighted = (TH2F*) fileIn->Get("h_genMET150to185_GridUnweighted");
    TH2F* h_genMET185to250_GridUnweighted = (TH2F*) fileIn->Get("h_genMET185to250_GridUnweighted");
    TH2F* h_genMET250_GridUnweighted = (TH2F*) fileIn->Get("h_genMET250_GridUnweighted");

    //----------------------------------
    //Weighted yields using genMET
    TH2F* h_genMET100to115_Grid = (TH2F*) h_genMET100to115_GridUnweighted->Clone("h_genMET100to115_Grid");
    TH2F* h_genMET115to130_Grid = (TH2F*) h_genMET100to115_GridUnweighted->Clone("h_genMET115to130_Grid");
    TH2F* h_genMET130to150_Grid = (TH2F*) h_genMET100to115_GridUnweighted->Clone("h_genMET130to150_Grid");
    TH2F* h_genMET150to185_Grid = (TH2F*) h_genMET100to115_GridUnweighted->Clone("h_genMET150to185_Grid");
    TH2F* h_genMET185to250_Grid = (TH2F*) h_genMET100to115_GridUnweighted->Clone("h_genMET185to250_Grid");
    TH2F* h_genMET250_Grid = (TH2F*) h_genMET100to115_GridUnweighted->Clone("h_genMET250_Grid");

    h_genMET100to115_Grid->Reset();
    h_genMET115to130_Grid->Reset();
    h_genMET130to150_Grid->Reset();
    h_genMET150to185_Grid->Reset();
    h_genMET185to250_Grid->Reset();
    h_genMET250_Grid->Reset();

    //---------------------------------- 
    //Weighted yields using pfMET
    TH2F* h_pfMET100to115_Grid = (TH2F*)  h_MET100to115_Grid->Clone("h_pfMET100to115_Grid");
    TH2F* h_pfMET115to130_Grid = (TH2F*)  h_MET100to115_Grid->Clone("h_pfMET115to130_Grid");
    TH2F* h_pfMET130to150_Grid = (TH2F*)  h_MET100to115_Grid->Clone("h_pfMET130to150_Grid");
    TH2F* h_pfMET150to185_Grid = (TH2F*)  h_MET100to115_Grid->Clone("h_pfMET150to185_Grid");
    TH2F* h_pfMET185to250_Grid = (TH2F*)  h_MET100to115_Grid->Clone("h_pfMET185to250_Grid");
    TH2F* h_pfMET250_Grid = (TH2F*)  h_MET100to115_Grid->Clone("h_pfMET250_Grid");

    h_pfMET100to115_Grid->Reset();
    h_pfMET115to130_Grid->Reset();
    h_pfMET130to150_Grid->Reset();
    h_pfMET150to185_Grid->Reset();
    h_pfMET185to250_Grid->Reset();
    h_pfMET250_Grid->Reset();

    //----------------------------------  
    //Weighted yields using pfMET_JESUp and Down 
    TH2F* h_MET100to115_Grid_JESDown  = (TH2F*) fileIn->Get("h_MET100to115_Grid_JESDown");
    TH2F* h_MET115to130_Grid_JESDown  = (TH2F*) fileIn->Get("h_MET115to130_Grid_JESDown");
    TH2F* h_MET130to150_Grid_JESDown  = (TH2F*) fileIn->Get("h_MET130to150_Grid_JESDown");
    TH2F* h_MET150to185_Grid_JESDown  = (TH2F*) fileIn->Get("h_MET150to185_Grid_JESDown");
    TH2F* h_MET185to250_Grid_JESDown  = (TH2F*) fileIn->Get("h_MET185to250_Grid_JESDown");
    TH2F* h_MET250_Grid_JESDown  = (TH2F*) fileIn->Get("h_MET250_Grid_JESDown");


    TH2F* h_MET100to115_Grid_JESUp  = (TH2F*) fileIn->Get("h_MET100to115_Grid_JESUp");
    TH2F* h_MET115to130_Grid_JESUp  = (TH2F*) fileIn->Get("h_MET115to130_Grid_JESUp");
    TH2F* h_MET130to150_Grid_JESUp  = (TH2F*) fileIn->Get("h_MET130to150_Grid_JESUp");
    TH2F* h_MET150to185_Grid_JESUp  = (TH2F*) fileIn->Get("h_MET150to185_Grid_JESUp");
    TH2F* h_MET185to250_Grid_JESUp  = (TH2F*) fileIn->Get("h_MET185to250_Grid_JESUp");
    TH2F* h_MET250_Grid_JESUp  = (TH2F*) fileIn->Get("h_MET250_Grid_JESUp");

    //----------------------------------
    //Unweighted yields using pfMET_JESUp and Down
    TH2F* h_MET100to115_GridUnweighted_JESDown = (TH2F*) fileIn->Get("h_MET100to115_GridUnweighted_JESDown");
    TH2F* h_MET115to130_GridUnweighted_JESDown = (TH2F*) fileIn->Get("h_MET115to130_GridUnweighted_JESDown");
    TH2F* h_MET130to150_GridUnweighted_JESDown = (TH2F*) fileIn->Get("h_MET130to150_GridUnweighted_JESDown");
    TH2F* h_MET150to185_GridUnweighted_JESDown = (TH2F*) fileIn->Get("h_MET150to185_GridUnweighted_JESDown");
    TH2F* h_MET185to250_GridUnweighted_JESDown = (TH2F*) fileIn->Get("h_MET185to250_GridUnweighted_JESDown");
    TH2F* h_MET250_GridUnweighted_JESDown = (TH2F*) fileIn->Get("h_MET250_GridUnweighted_JESDown");

    TH2F* h_MET100to115_GridUnweighted_JESUp = (TH2F*) fileIn->Get("h_MET100to115_GridUnweighted_JESUp");
    TH2F* h_MET115to130_GridUnweighted_JESUp = (TH2F*) fileIn->Get("h_MET115to130_GridUnweighted_JESUp");
    TH2F* h_MET130to150_GridUnweighted_JESUp = (TH2F*) fileIn->Get("h_MET130to150_GridUnweighted_JESUp");
    TH2F* h_MET150to185_GridUnweighted_JESUp = (TH2F*) fileIn->Get("h_MET150to185_GridUnweighted_JESUp");
    TH2F* h_MET185to250_GridUnweighted_JESUp = (TH2F*) fileIn->Get("h_MET185to250_GridUnweighted_JESUp");
    TH2F* h_MET250_GridUnweighted_JESUp = (TH2F*) fileIn->Get("h_MET250_GridUnweighted_JESUp");

    //----------------------------------
    //Histograms to hold the JES uncertainties from pfMET_JESUp and Down

    TH2F* h_MET100to115_GridJESDownError =(TH2F*) fileIn->Get("h_MET100to115_GridJESDownError");
    TH2F* h_MET115to130_GridJESDownError =(TH2F*) fileIn->Get("h_MET115to130_GridJESDownError");
    TH2F* h_MET130to150_GridJESDownError =(TH2F*) fileIn->Get("h_MET130to150_GridJESDownError");
    TH2F* h_MET150to185_GridJESDownError =(TH2F*) fileIn->Get("h_MET150to185_GridJESDownError");
    TH2F* h_MET185to250_GridJESDownError =(TH2F*) fileIn->Get("h_MET185to250_GridJESDownError");
    TH2F* h_MET250_GridJESDownError =(TH2F*) fileIn->Get("h_MET250_GridJESDownError");

    TH2F* h_MET100to115_GridJESUpError =(TH2F*) fileIn->Get("h_MET100to115_GridJESUpError");
    TH2F* h_MET115to130_GridJESUpError =(TH2F*) fileIn->Get("h_MET115to130_GridJESUpError");
    TH2F* h_MET130to150_GridJESUpError =(TH2F*) fileIn->Get("h_MET130to150_GridJESUpError");
    TH2F* h_MET150to185_GridJESUpError =(TH2F*) fileIn->Get("h_MET150to185_GridJESUpError");
    TH2F* h_MET185to250_GridJESUpError =(TH2F*) fileIn->Get("h_MET185to250_GridJESUpError");
    TH2F* h_MET250_GridJESUpError =(TH2F*) fileIn->Get("h_MET250_GridJESUpError");

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
    TH2F* h_MET100to115_GridGenMetError =(TH2F*) h_MET100to115_Grid->Clone("h_MET100to115_GridGenMetError");
    TH2F* h_MET115to130_GridGenMetError =(TH2F*) h_MET100to115_Grid->Clone("h_MET115to130_GridGenMetError");
    TH2F* h_MET130to150_GridGenMetError =(TH2F*) h_MET100to115_Grid->Clone("h_MET130to150_GridGenMetError");
    TH2F* h_MET150to185_GridGenMetError =(TH2F*) h_MET100to115_Grid->Clone("h_MET150to185_GridGenMetError");
    TH2F* h_MET185to250_GridGenMetError =(TH2F*) h_MET100to115_Grid->Clone("h_MET185to250_GridGenMetError");
    TH2F* h_MET250_GridGenMetError =(TH2F*) h_MET100to115_Grid->Clone("h_MET250_GridGenMetError");

    h_MET100to115_GridGenMetError->Reset();
    h_MET115to130_GridGenMetError->Reset();
    h_MET130to150_GridGenMetError->Reset();
    h_MET150to185_GridGenMetError->Reset();
    h_MET185to250_GridGenMetError->Reset();
    h_MET250_GridGenMetError->Reset();

    //------------------------------------- 
    //Unweighted number of events for low and high NVertex
    TH2F* h1_lowAll = (TH2F*) fileIn->Get("h_MET100to115_GridAll_LowNvtx");
    TH2F* h2_lowAll = (TH2F*) fileIn->Get("h_MET115to130_GridAll_LowNvtx");
    TH2F* h3_lowAll = (TH2F*) fileIn->Get("h_MET130to150_GridAll_LowNvtx");
    TH2F* h4_lowAll = (TH2F*) fileIn->Get("h_MET150to185_GridAll_LowNvtx");
    TH2F* h5_lowAll = (TH2F*) fileIn->Get("h_MET185to250_GridAll_LowNvtx");
    TH2F* h6_lowAll = (TH2F*) fileIn->Get("h_MET250_GridAll_LowNvtx");

    TH2F* h1_highAll = (TH2F*) fileIn->Get("h_MET100to115_GridAll_HighNvtx");
    TH2F* h2_highAll = (TH2F*) fileIn->Get("h_MET115to130_GridAll_HighNvtx");
    TH2F* h3_highAll = (TH2F*) fileIn->Get("h_MET130to150_GridAll_HighNvtx");
    TH2F* h4_highAll = (TH2F*) fileIn->Get("h_MET150to185_GridAll_HighNvtx");
    TH2F* h5_highAll = (TH2F*) fileIn->Get("h_MET185to250_GridAll_HighNvtx");
    TH2F* h6_highAll = (TH2F*) fileIn->Get("h_MET250_GridAll_HighNvtx");

    TH2F* h1_lowGG = (TH2F*) fileIn->Get("h_MET100to115_GridGG_LowNvtx");
    TH2F* h2_lowGG = (TH2F*) fileIn->Get("h_MET115to130_GridGG_LowNvtx");
    TH2F* h3_lowGG = (TH2F*) fileIn->Get("h_MET130to150_GridGG_LowNvtx");
    TH2F* h4_lowGG = (TH2F*) fileIn->Get("h_MET150to185_GridGG_LowNvtx");
    TH2F* h5_lowGG = (TH2F*) fileIn->Get("h_MET185to250_GridGG_LowNvtx");
    TH2F* h6_lowGG = (TH2F*) fileIn->Get("h_MET250_GridGG_LowNvtx");

    TH2F* h1_highGG = (TH2F*) fileIn->Get("h_MET100to115_GridGG_HighNvtx");
    TH2F* h2_highGG = (TH2F*) fileIn->Get("h_MET115to130_GridGG_HighNvtx");
    TH2F* h3_highGG = (TH2F*) fileIn->Get("h_MET130to150_GridGG_HighNvtx");
    TH2F* h4_highGG = (TH2F*) fileIn->Get("h_MET150to185_GridGG_HighNvtx");
    TH2F* h5_highGG = (TH2F*) fileIn->Get("h_MET185to250_GridGG_HighNvtx");
    TH2F* h6_highGG = (TH2F*) fileIn->Get("h_MET250_GridGG_HighNvtx");
    //----------------------------------
    //Histograms to hold R for the full grid
    TH2F* h_MET100to115_NVtxEffRatio =(TH2F*) h_MET100to115_Grid->Clone("h_MET100to115_NVtxEffRatio");
    TH2F* h_MET115to130_NVtxEffRatio =(TH2F*) h_MET100to115_Grid->Clone("h_MET115to130_NVtxEffRatio");
    TH2F* h_MET130to150_NVtxEffRatio =(TH2F*) h_MET100to115_Grid->Clone("h_MET130to150_NVtxEffRatio");
    TH2F* h_MET150to185_NVtxEffRatio =(TH2F*) h_MET100to115_Grid->Clone("h_MET150to185_NVtxEffRatio");
    TH2F* h_MET185to250_NVtxEffRatio =(TH2F*) h_MET100to115_Grid->Clone("h_MET185to250_NVtxEffRatio");
    TH2F* h_MET250_NVtxEffRatio =(TH2F*) h_MET100to115_Grid->Clone("h_MET250_NVtxEffRatio");

    h_MET100to115_NVtxEffRatio->Reset();
    h_MET115to130_NVtxEffRatio->Reset();
    h_MET130to150_NVtxEffRatio->Reset();
    h_MET150to185_NVtxEffRatio->Reset();
    h_MET185to250_NVtxEffRatio->Reset();
    h_MET250_NVtxEffRatio->Reset();

    TH2F* h_MET100to115_NVtxEffRatioErr =(TH2F*) h_MET100to115_Grid->Clone("h_MET100to115_NVtxEffRatioErr");
    TH2F* h_MET115to130_NVtxEffRatioErr =(TH2F*) h_MET100to115_Grid->Clone("h_MET115to130_NVtxEffRatioErr");
    TH2F* h_MET130to150_NVtxEffRatioErr =(TH2F*) h_MET100to115_Grid->Clone("h_MET130to150_NVtxEffRatioErr");
    TH2F* h_MET150to185_NVtxEffRatioErr =(TH2F*) h_MET100to115_Grid->Clone("h_MET150to185_NVtxEffRatioErr");
    TH2F* h_MET185to250_NVtxEffRatioErr =(TH2F*) h_MET100to115_Grid->Clone("h_MET185to250_NVtxEffRatioErr");
    TH2F* h_MET250_NVtxEffRatioErr =(TH2F*) h_MET100to115_Grid->Clone("h_MET250_NVtxEffRatioErr");

    h_MET100to115_NVtxEffRatioErr->Reset();
    h_MET115to130_NVtxEffRatioErr->Reset();
    h_MET130to150_NVtxEffRatioErr->Reset();
    h_MET150to185_NVtxEffRatioErr->Reset();
    h_MET185to250_NVtxEffRatioErr->Reset();
    h_MET250_NVtxEffRatioErr->Reset();

    //---------------------------------- 
    float lumi = 35.9218*1000;
    /*  Double_t xbins3[31]={25,75,125,175, 225, 275, 325, 375, 425, 475, 525, 575, 625, 675, 725, 775, 825, 875, 925, 975, 1025, 1075, 1125, 1175, 1225, 1275, 1325, 1375, 1425, 1475, 1525};
    int xbinnum3 = 30;
        
        Double_t ybins[32]={975, 1025, 1075, 1125, 1175, 1225, 1275, 1325, 1375, 1425, 1475, 1525,
                        1575, 1625, 1675, 1725, 1775, 1825, 1875, 1925, 1975,
                        2025, 2075, 2125, 2175, 2225, 2275, 2325, 2375, 2425, 2475, 2525,
    };
    int ybinnum = 31;
    */

    Double_t xbins3[27]={225, 275, 325, 375, 425, 475, 525, 575, 625, 675, 725, 775, 825, 875, 925, 975, 1025, 1075, 1125, 1175, 1225, 1275, 1325, 1375, 1425, 1475, 1525};
    int xbinnum3 = 26;
    Double_t ybins[27]={225, 275, 325, 375, 425, 475, 525, 575, 625, 675, 725, 775, 825, 875, 925, 975, 1025, 1075, 1125, 1175, 1225, 1275, 1325, 1375, 1425, 1475, 1525};
    int ybinnum = 26;


    ifstream M1M2_xsec ("M1M2_xsec.txt");

    unordered_map<int, float> xsec_map_M1M2;
    
    string line;
    string tempstr;
    stringstream ss;
    int M1, M2;
    float xsec, uncert;
    if(M1M2_xsec.is_open()){
      cout << "Opened file" << endl;
        while( M1M2_xsec >> M1 >> M2 >> xsec >> uncert){
            xsec_map_M1M2[key(M1,M2)] = xsec;
	    // cout << "Filled " << M1 << ", " << M2 << " with " << xsec << endl;
        }
    }
    else{
      cout <<"File not open" << endl;
    }    
    for(int i = 1; i <= xbinnum3;i++){
        for(int j = 1; j <= ybinnum;j++){
	  //Notice that the bin histogram starts counting at 1 but the array starts at 0. So the first bin (i=1) is given by xbins3[1] -25 = 275-25 = 250
            float mGlu = xbins3[i] - 25;
            float mNeu = ybins[j] - 25;
            float xsec = 0;
            try{
                xsec = xsec_map_M1M2[key(mGlu,mNeu)];
            }
            catch(...){
                cout<<"Uh oh" << endl;
                break;
            }
	    float scaleByLumi = 0;
            if( h_GridAllEvents->GetBinContent(i,j) > 0 ){
                scaleByLumi = lumi*xsec/h_GridAllEvents->GetBinContent(i,j);
	    }
            else{
              cout << "No events for " << mGlu << ", " <<  mNeu << " :(" << endl;
            }	      
	    //Scale pf, gen, and JES MET distributions by lumi
	    h_pfMET100to115_Grid->SetBinContent(i,j,h_MET100to115_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
	    h_pfMET115to130_Grid->SetBinContent(i,j,h_MET115to130_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
	    h_pfMET130to150_Grid->SetBinContent(i,j,h_MET130to150_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
	    h_pfMET150to185_Grid->SetBinContent(i,j,h_MET150to185_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
	    h_pfMET185to250_Grid->SetBinContent(i,j,h_MET185to250_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
	    h_pfMET250_Grid->SetBinContent(i,j,h_MET250_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
            
	    h_MET100to115_Grid_JESDown->SetBinContent(i,j,h_MET100to115_GridUnweighted_JESDown->GetBinContent(i,j)*scaleByLumi);
	    h_MET115to130_Grid_JESDown->SetBinContent(i,j,h_MET115to130_GridUnweighted_JESDown->GetBinContent(i,j)*scaleByLumi);
	    h_MET130to150_Grid_JESDown->SetBinContent(i,j,h_MET130to150_GridUnweighted_JESDown->GetBinContent(i,j)*scaleByLumi);
	    h_MET150to185_Grid_JESDown->SetBinContent(i,j,h_MET150to185_GridUnweighted_JESDown->GetBinContent(i,j)*scaleByLumi);
	    h_MET185to250_Grid_JESDown->SetBinContent(i,j,h_MET185to250_GridUnweighted_JESDown->GetBinContent(i,j)*scaleByLumi);
	    h_MET250_Grid_JESDown->SetBinContent(i,j,h_MET250_GridUnweighted_JESDown->GetBinContent(i,j)*scaleByLumi);
                
	    h_MET100to115_Grid_JESUp->SetBinContent(i,j,h_MET100to115_GridUnweighted_JESUp->GetBinContent(i,j)*scaleByLumi);
	    h_MET115to130_Grid_JESUp->SetBinContent(i,j,h_MET115to130_GridUnweighted_JESUp->GetBinContent(i,j)*scaleByLumi);
	    h_MET130to150_Grid_JESUp->SetBinContent(i,j,h_MET130to150_GridUnweighted_JESUp->GetBinContent(i,j)*scaleByLumi);
	    h_MET150to185_Grid_JESUp->SetBinContent(i,j,h_MET150to185_GridUnweighted_JESUp->GetBinContent(i,j)*scaleByLumi);
	    h_MET185to250_Grid_JESUp->SetBinContent(i,j,h_MET185to250_GridUnweighted_JESUp->GetBinContent(i,j)*scaleByLumi);
	    h_MET250_Grid_JESUp->SetBinContent(i,j,h_MET250_GridUnweighted_JESUp->GetBinContent(i,j)*scaleByLumi);

	    h_genMET100to115_Grid->SetBinContent(i,j,h_genMET100to115_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
            h_genMET115to130_Grid->SetBinContent(i,j,h_genMET115to130_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
            h_genMET130to150_Grid->SetBinContent(i,j,h_genMET130to150_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
            h_genMET150to185_Grid->SetBinContent(i,j,h_genMET150to185_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
            h_genMET185to250_Grid->SetBinContent(i,j,h_genMET185to250_GridUnweighted->GetBinContent(i,j)*scaleByLumi);
            h_genMET250_Grid->SetBinContent(i,j,h_genMET250_GridUnweighted->GetBinContent(i,j)*scaleByLumi);

	    //Actual prediction is average of pf and gen MET
	    h_MET100to115_Grid->SetBinContent(i,j, ( h_pfMET100to115_Grid->GetBinContent(i,j) + h_genMET100to115_Grid->GetBinContent(i,j) ) / 2) ;
	    h_MET115to130_Grid->SetBinContent(i,j, ( h_pfMET115to130_Grid->GetBinContent(i,j) + h_genMET115to130_Grid->GetBinContent(i,j) ) / 2) ;
            h_MET130to150_Grid->SetBinContent(i,j, ( h_pfMET130to150_Grid->GetBinContent(i,j) + h_genMET130to150_Grid->GetBinContent(i,j) ) / 2) ;
            h_MET150to185_Grid->SetBinContent(i,j, ( h_pfMET150to185_Grid->GetBinContent(i,j) + h_genMET150to185_Grid->GetBinContent(i,j) ) / 2) ;
            h_MET185to250_Grid->SetBinContent(i,j, ( h_pfMET185to250_Grid->GetBinContent(i,j) + h_genMET185to250_Grid->GetBinContent(i,j) ) / 2) ;
            h_MET250_Grid->SetBinContent(i,j, ( h_pfMET250_Grid->GetBinContent(i,j) + h_genMET250_Grid->GetBinContent(i,j) ) / 2) ;

	    //Gen MET uncertainty is half the difference between pf and gen MET
	    //So the percent uncertainty is [abs(gen-pf) / 2 ] / [ (gen+pf) / 2 ] = (gen-pf) / (gen+pf)
	    if( h_MET100to115_Grid->GetBinContent(i,j) > 0){
	      h_MET100to115_GridGenMetError->
		SetBinContent(i,j,  abs(h_genMET100to115_Grid->GetBinContent(i,j) - h_pfMET100to115_Grid->GetBinContent(i,j)) / 
			      (h_genMET100to115_Grid->GetBinContent(i,j) + h_pfMET100to115_Grid->GetBinContent(i,j)) );
	    }
	    else{
	      h_MET100to115_GridGenMetError->SetBinContent(i,j,0.0);
	    }
	    if(h_MET115to130_Grid->GetBinContent(i,j) > 0){
              h_MET115to130_GridGenMetError->
                SetBinContent(i,j,  abs(h_genMET115to130_Grid->GetBinContent(i,j) - h_pfMET115to130_Grid->GetBinContent(i,j)) /
			      (h_genMET115to130_Grid->GetBinContent(i,j) + h_pfMET115to130_Grid->GetBinContent(i,j)) );
            }
            else{
              h_MET115to130_GridGenMetError->SetBinContent(i,j,0.0);
            }
	    if(h_MET130to150_Grid->GetBinContent(i,j) > 0){
              h_MET130to150_GridGenMetError->
                SetBinContent(i,j,  abs(h_genMET130to150_Grid->GetBinContent(i,j) - h_pfMET130to150_Grid->GetBinContent(i,j)) /
			      (h_genMET130to150_Grid->GetBinContent(i,j) + h_pfMET130to150_Grid->GetBinContent(i,j)) );
            }
            else{
              h_MET130to150_GridGenMetError->SetBinContent(i,j,0.0);
            }
	    if(h_MET150to185_Grid->GetBinContent(i,j) > 0){
              h_MET150to185_GridGenMetError->
                SetBinContent(i,j,  abs(h_genMET150to185_Grid->GetBinContent(i,j) - h_pfMET150to185_Grid->GetBinContent(i,j)) /
			      (h_genMET150to185_Grid->GetBinContent(i,j) + h_pfMET150to185_Grid->GetBinContent(i,j)) );
            }
            else{
              h_MET150to185_GridGenMetError->SetBinContent(i,j,0.0);
            }
	    if(h_MET185to250_Grid->GetBinContent(i,j) > 0){
              h_MET185to250_GridGenMetError->
                SetBinContent(i,j,  abs(h_genMET185to250_Grid->GetBinContent(i,j) - h_pfMET185to250_Grid->GetBinContent(i,j)) /
			      (h_genMET185to250_Grid->GetBinContent(i,j) + h_pfMET185to250_Grid->GetBinContent(i,j)) );
            }
            else{
              h_MET185to250_GridGenMetError->SetBinContent(i,j,0.0);
            }
	    if(h_MET250_Grid->GetBinContent(i,j) > 0){
              h_MET250_GridGenMetError->
                SetBinContent(i,j,  abs(h_genMET250_Grid->GetBinContent(i,j) - h_pfMET250_Grid->GetBinContent(i,j)) /
			      (h_genMET250_Grid->GetBinContent(i,j) + h_pfMET250_Grid->GetBinContent(i,j)) );
            }
            else{
              h_MET250_GridGenMetError->SetBinContent(i,j,0.0);
            }

            //JES uncertainty is the difference between the pf uncertainty and the JES Met distributions
	    //These histograms store the percent difference
	    if(h_pfMET100to115_Grid->GetBinContent(i,j) > 0) {	    
	      h_MET100to115_GridJESUpError->Fill(mGlu,mNeu, 
						 abs( h_MET100to115_Grid_JESUp->GetBinContent(i,j) - 
						      h_pfMET100to115_Grid->GetBinContent(i,j)) / 
						 h_pfMET100to115_Grid->GetBinContent(i,j) );
              h_MET100to115_GridJESDownError->Fill(mGlu,mNeu, 
						   abs( h_MET100to115_Grid_JESDown->GetBinContent(i,j) - 
							h_pfMET100to115_Grid->GetBinContent(i,j)) / 
						   h_pfMET100to115_Grid->GetBinContent(i,j) );
	    }
	    else{
	      h_MET100to115_GridJESUpError->Fill(mGlu,mNeu,0.0);
	      h_MET100to115_GridJESDownError->Fill(mGlu,mNeu,0.0);
	    }

            if(h_pfMET115to130_Grid->GetBinContent(i,j) > 0) {
              h_MET115to130_GridJESUpError->Fill(mGlu,mNeu,
                                                 abs( h_MET115to130_Grid_JESUp->GetBinContent(i,j) -
                                                      h_pfMET115to130_Grid->GetBinContent(i,j)) /
                                                 h_pfMET115to130_Grid->GetBinContent(i,j) );
              h_MET115to130_GridJESDownError->Fill(mGlu,mNeu,
                                                   abs( h_MET115to130_Grid_JESDown->GetBinContent(i,j) -
                                                        h_pfMET115to130_Grid->GetBinContent(i,j)) /
                                                   h_pfMET115to130_Grid->GetBinContent(i,j) );
            }
            else{
              h_MET115to130_GridJESUpError->Fill(mGlu,mNeu,0.0);
              h_MET115to130_GridJESDownError->Fill(mGlu,mNeu,0.0);
            }


            if(h_pfMET130to150_Grid->GetBinContent(i,j) > 0) {
              h_MET130to150_GridJESUpError->Fill(mGlu,mNeu,
                                                 abs( h_MET130to150_Grid_JESUp->GetBinContent(i,j) -
                                                      h_pfMET130to150_Grid->GetBinContent(i,j)) /
                                                 h_pfMET130to150_Grid->GetBinContent(i,j) );
              h_MET130to150_GridJESDownError->Fill(mGlu,mNeu,
                                                   abs( h_MET130to150_Grid_JESDown->GetBinContent(i,j) -
                                                        h_pfMET130to150_Grid->GetBinContent(i,j)) /
                                                   h_pfMET130to150_Grid->GetBinContent(i,j) );
            }
            else{
              h_MET130to150_GridJESUpError->Fill(mGlu,mNeu,0.0);
              h_MET130to150_GridJESDownError->Fill(mGlu,mNeu,0.0);
            }


	    if(h_pfMET150to185_Grid->GetBinContent(i,j) > 0) {
              h_MET150to185_GridJESUpError->Fill(mGlu,mNeu,
                                                 abs( h_MET150to185_Grid_JESUp->GetBinContent(i,j) -
                                                      h_pfMET150to185_Grid->GetBinContent(i,j)) /
                                                 h_pfMET150to185_Grid->GetBinContent(i,j) );
              h_MET150to185_GridJESDownError->Fill(mGlu,mNeu,
                                                   abs( h_MET150to185_Grid_JESDown->GetBinContent(i,j) -
                                                        h_pfMET150to185_Grid->GetBinContent(i,j)) /
                                                   h_pfMET150to185_Grid->GetBinContent(i,j) );
            }
            else{
              h_MET150to185_GridJESUpError->Fill(mGlu,mNeu,0.0);
              h_MET150to185_GridJESDownError->Fill(mGlu,mNeu,0.0);
            }

            if(h_pfMET185to250_Grid->GetBinContent(i,j) > 0) {
              h_MET185to250_GridJESUpError->Fill(mGlu,mNeu,
                                                 abs( h_MET185to250_Grid_JESUp->GetBinContent(i,j) -
                                                      h_pfMET185to250_Grid->GetBinContent(i,j)) /
                                                 h_pfMET185to250_Grid->GetBinContent(i,j) );
              h_MET185to250_GridJESDownError->Fill(mGlu,mNeu,
                                                   abs( h_MET185to250_Grid_JESDown->GetBinContent(i,j) -
                                                        h_pfMET185to250_Grid->GetBinContent(i,j)) /
                                                   h_pfMET185to250_Grid->GetBinContent(i,j) );
            }
            else{
              h_MET185to250_GridJESUpError->Fill(mGlu,mNeu,0.0);
              h_MET185to250_GridJESDownError->Fill(mGlu,mNeu,0.0);
            }

            if(h_pfMET250_Grid->GetBinContent(i,j) > 0) {
              h_MET250_GridJESUpError->Fill(mGlu,mNeu,
                                                 abs( h_MET250_Grid_JESUp->GetBinContent(i,j) -
                                                      h_pfMET250_Grid->GetBinContent(i,j)) /
                                                 h_pfMET250_Grid->GetBinContent(i,j) );
              h_MET250_GridJESDownError->Fill(mGlu,mNeu,
                                                   abs( h_MET250_Grid_JESDown->GetBinContent(i,j) -
                                                        h_pfMET250_Grid->GetBinContent(i,j)) /
                                                   h_pfMET250_Grid->GetBinContent(i,j) );
            }
            else{
              h_MET250_GridJESUpError->Fill(mGlu,mNeu,0.0);
              h_MET250_GridJESDownError->Fill(mGlu,mNeu,0.0);
            }

	    //Calculate r = acceptance for high nvtx / acceptance for low nvtx
	    //Bin 1
	    float lA =  h1_lowAll->GetBinContent(i,j);
	    float lG =  h1_lowGG->GetBinContent(i,j);
	    float hA = h1_highAll->GetBinContent(i,j);
	    float hG = h1_highGG->GetBinContent(i,j);

	    float r = -1;
	    float r_err = -1;
	    if( lA > 0 && lG > 0 && hA > 0 && hG > 0){
	      r = (hG/hA) / (lG/lA);
	      r_err = r*TMath::Sqrt( (1/lA) + (1/lG) + (1/hA) + (1/hG) );
	    }

	    h_MET100to115_NVtxEffRatio->SetBinContent(i,j, r);
            h_MET100to115_NVtxEffRatioErr->SetBinContent(i,j, r_err);


	    //Bin 2
	    lA = h2_lowAll->GetBinContent(i,j);
            lG = h2_lowGG->GetBinContent(i,j);
            hA = h2_highAll->GetBinContent(i,j);
            hG = h2_highGG->GetBinContent(i,j);

            r = -1;
            r_err = -1;
            if( lA > 0 && lG > 0 && hA > 0 && hG > 0){
              r = (hG/hA) / (lG/lA);
              r_err = r*TMath::Sqrt( (1/lA) + (1/lG) + (1/hA) + (1/hG) );
            }

            h_MET115to130_NVtxEffRatio->SetBinContent(i,j, r);
            h_MET115to130_NVtxEffRatioErr->SetBinContent(i,j, r_err);

	    totLA += lA;
            totLG += lG;
            totHA += hA;
            totHG += hG;

	    //Bin 3
            lA = h3_lowAll->GetBinContent(i,j);
            lG = h3_lowGG->GetBinContent(i,j);
            hA = h3_highAll->GetBinContent(i,j);
            hG = h3_highGG->GetBinContent(i,j);

            r = -1;
            r_err = -1;
            if( lA > 0 && lG > 0 && hA > 0 && hG > 0){
              r = (hG/hA) / (lG/lA);
              r_err = r*TMath::Sqrt( (1/lA) + (1/lG) + (1/hA) + (1/hG) );
            }

            h_MET130to150_NVtxEffRatio->SetBinContent(i,j, r);
            h_MET130to150_NVtxEffRatioErr->SetBinContent(i,j, r_err);

	    totLA += lA;
            totLG += lG;
            totHA += hA;
            totHG += hG;

	    //Bin 4
            lA = h4_lowAll->GetBinContent(i,j);
            lG = h4_lowGG->GetBinContent(i,j);
            hA = h4_highAll->GetBinContent(i,j);
            hG = h4_highGG->GetBinContent(i,j);

            r = -1;
            r_err = -1;
            if( lA > 0 && lG > 0 && hA > 0 && hG > 0){
              r = (hG/hA) / (lG/lA);
              r_err = r*TMath::Sqrt( (1/lA) + (1/lG) + (1/hA) + (1/hG) );
            }

            h_MET150to185_NVtxEffRatio->SetBinContent(i,j, r);
            h_MET150to185_NVtxEffRatioErr->SetBinContent(i,j, r_err);

	    totLA += lA;
            totLG += lG;
            totHA += hA;
            totHG += hG;

	    //Bin 5
            lA = h5_lowAll->GetBinContent(i,j);
            lG = h5_lowGG->GetBinContent(i,j);
            hA = h5_highAll->GetBinContent(i,j);
            hG = h5_highGG->GetBinContent(i,j);

            r = -1;
            r_err = -1;
            if( lA > 0 && lG > 0 && hA > 0 && hG > 0){
              r = (hG/hA) / (lG/lA);
              r_err = r*TMath::Sqrt( (1/lA) + (1/lG) + (1/hA) + (1/hG) );
            }

            h_MET185to250_NVtxEffRatio->SetBinContent(i,j, r);
            h_MET185to250_NVtxEffRatioErr->SetBinContent(i,j, r_err);

	    totLA += lA;
            totLG += lG;
            totHA += hA;
            totHG += hG;

	    //Bin 6
            lA = h6_lowAll->GetBinContent(i,j);
            lG = h6_lowGG->GetBinContent(i,j);
            hA = h6_highAll->GetBinContent(i,j);
            hG = h6_highGG->GetBinContent(i,j);

            r = -1;
            r_err = -1;
            if( lA > 0 && lG > 0 && hA > 0 && hG > 0){
              r = (hG/hA) / (lG/lA);
              r_err = r*TMath::Sqrt( (1/lA) + (1/lG) + (1/hA) + (1/hG) );
            }

            h_MET250_NVtxEffRatio->SetBinContent(i,j, r);
            h_MET250_NVtxEffRatioErr->SetBinContent(i,j, r_err);

	    totLA += lA;
            totLG += lG;
            totHA += hA;
            totHG += hG;

        }
    }

    cout << "totLA = " << totLA << endl;
    cout << "totLG = " << totLG << endl;
    cout << "totHA = " << totHA << endl;
    cout << "totHG = " << totHG << endl;
    rTot = (totHG/totHA)/(totLG/totLA);
    rErr = rTot*TMath::Sqrt( (1/totHA) + (1/totHG) + (1/totLA) + (1/totLG) );

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
