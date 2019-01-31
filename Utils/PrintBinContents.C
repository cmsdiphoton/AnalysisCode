
void PrintEventContents(TString filename, TString histname, bool verbose){
  TFile *_file0 = TFile::Open(filename);
  TH1F * hist = (TH1F*) _file0->Get(histname);
  for (int i = 0;i < hist->GetNbinsX()+2;i++){
    if(verbose){
      cout << hist->GetBinLowEdge(i) << " to " << hist->GetBinLowEdge(i)+hist->GetBinWidth(i) << " : " << hist->GetBinContent(i) << " events." << endl;
    }
    else{
      cout << hist->GetBinContent(i) << endl;
    }
  }
  return;
}


void PrintErrorUp(TString filename, TString histname, bool verbose){
    TFile *_file0 = TFile::Open(filename);
    TH1F * hist = (TH1F*) _file0->Get(histname);
    for (int i = 0;i < hist->GetNbinsX()+2;i++){
        if(verbose){
            cout << hist->GetBinLowEdge(i) << " to " << hist->GetBinLowEdge(i)+hist->GetBinWidth(i) << " : " << hist->GetBinErrorUp(i) << " events." << endl;
        }
        else{
            cout << hist->GetBinErrorUp(i) << endl;
        }
    }
    return;
}

void PrintErrorDown(TString filename, TString histname, bool verbose){
    TFile *_file0 = TFile::Open(filename);
    TH1F * hist = (TH1F*) _file0->Get(histname);
    for (int i = 0;i < hist->GetNbinsX()+2;i++){
        if(verbose){
            cout << hist->GetBinLowEdge(i) << " to " << hist->GetBinLowEdge(i)+hist->GetBinWidth(i) << " : " << hist->GetBinErrorLow(i) << " events." << endl;
        }
        else{
            cout << hist->GetBinErrorLow(i) << endl;
        }
    }
    return;
}
