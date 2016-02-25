#include <iostream>
#include <fstream>
#include <vector>
#include "TF1.h"
#include "TString.h"

const int nNcoll = 1400;
const int nNpart = 450;
const int nBinNcoll = 1400;
const int nBinNpart = 450;

void test1(Int_t nf = 1, Float_t f0 = 0.0, Float_t f1 = 1.0, Int_t nsigma = 10, Int_t nEvents = 10000)
{
   gSystem->Load("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/build/libCentrality");
   
   
   TH1F *h1d_Mult_TPC_Ref = new TH1F("h1d_Mult_TPC_Ref","M_{TPC}^{ref}",125,0,250);

//     gSystem->LoadMacro("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Glauber/GlauberFitter.cxx");

//     TCanvas *c1 = new TCanvas("c1", "canvas", 1500, 900);

//     TFile *fData = new TFile ("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/root_files/NA61_run_23001.root");
//     TTree* fDataTree = (TTree*) fData->Get("na61_data");

//     fDataTree->Draw("CentralityEventContainer.fDetectorEvents.fWeights >> hData (125, 0, 250)", "CentralityEventContainer.fDetectorEvents.fDetId == 3");
//     TH1F* hData = (TH1F*)gPad->GetPrimitive("hData");
//     hData->GetYaxis()->SetRangeUser(0.1, 2e5);
//     hData->Draw();
//     fData->Close();
    
    GlauberFitter *test = new GlauberFitter;
    
    test->SetNpartMax (450);
    test->SetNcollMax (1400);
    test->SetMultMax (250);
    test->SetBinSize (2);
    test->SetFitMultMin (100);
    
    
    test->SetSimHistos ("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Glauber/root_files/glau_pbpb_ntuple.root");
    
    
//     test->SetGlauberFitHisto (0.5, 0.01, 1, 1000000, false);
//     hData = test->GetGlauberFitHisto ();
    TFile *f = new TFile("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Glauber/root_files/in_histo_1.root");    
//     h1d_Mult_TPC_Ref->Write();
//     f->Close();
//     hData = (TH1F*)f->Get("glaub_0.50_0.0100_1.00_1000000");
    h1d_Mult_TPC_Ref = (TH1F*)f->Get("h1d_Mult_TPC_Ref");
    
    test->SetInputHisto (h1d_Mult_TPC_Ref);
//     test->SetInputHisto (hData);
    test->TestFunc(nf, f0, f1, nsigma, nEvents);
//     test->FitGlauber (20000);
//     test->DrawHistos (true, true, true/*, true*/);

//     test->MinuitFit(0.9, 0.3, 1);
    
//     f->Close();

//     TFile *file1 = new TFile("root_files/out_histo_2.root", "recreate");    
// 
// 
//     test->SetGlauberFitHisto (0.25, 0.170897, 0.783594, 1000000);
//     hData = test->GetGlauberFitHisto ();
//     hData->Write();
// 
//     h1d_Mult_TPC_Ref->SetLineColor(2);
//     hData->Draw();
//     h1d_Mult_TPC_Ref->Draw("same");
//     h1d_Mult_TPC_Ref->Write();
// 
// 
//     test->SetGlauberFitHisto (0.8, 0.304777, 0.751315, 1000000);
//     hData->SetLineColor(3);
//     hData = test->GetGlauberFitHisto ();
//     hData->Write();
// 
//     hData->Draw("same");
// 
//     test->SetGlauberFitHisto (0.5, 0.2, 1, 10000, false);
//     hData = test->GetGlauberFitHisto ();
//     hData->Draw("same");
// 
//     file1->Close();
}