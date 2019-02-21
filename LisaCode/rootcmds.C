{
  gROOT->ProcessLine(".x yourchain.C");
  chain->Process("T5gg_analyzer.C++");
}
