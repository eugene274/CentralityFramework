

void RunGlauberFitter (Int_t nf = 1, Float_t f0 = 1, Float_t f1 = 1, Int_t nsigma = 3, Float_t StepSigma = 0.15, Int_t nEvents = 1000000, Int_t MultMin = 50)
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/git/CentralityFramework/";
//     gStyle->SetOptStat(0000);    
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");    
    
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //     
    TH1F *h1d_Mult_TPC_Ref/* = new TH1F("h1d_Mult_TPC_Ref","M_{TPC}^{ref}", 125, 0, 250)*/;
    TFile *f = new TFile("/lustre/nyx/cbm/users/klochkov/git/CentralityFramework/Macro/hTPC_ref.root");    
    h1d_Mult_TPC_Ref = (TH1F*)f->Get("h1Corr");
    TH1F *hData1 = h1d_Mult_TPC_Ref;
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
    
//     TString DataFileName = CentralityFrameworkDir + "containers/ana_dima_merged.root";
//     TString DataFileName = "/lustre/nyx/cbm/users/klochkov/soft/PidFramework/input/DCM_1M.root";
//     TString DataFileName = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/root_files/na61_cont/na61_container_merged.root";
    
    
//     TFile *DataFile = new TFile ( DataFileName, "read" );    
//     CentralityEventContainer *container = new CentralityEventContainer;
//     TTree *ContTree = (TTree*)DataFile->Get("cbm_data");
// //     TTree *ContTree = (TTree*)DataFile->Get("na61_data");
//     ContTree->SetBranchAddress("CentralityEventContainer", &container);
//     ContTree->Draw("CentralityEventContainer.GetDetectorWeight(0) >> h1(500, 0, 500)", "", "");
//     TH1F *hData1 = (TH1F*)gPad->GetPrimitive("h1");
    
    GlauberFitter *test = new GlauberFitter;

    test->SetNpartMax (450);
    test->SetNcollMax (1100);
    test->SetMultMax (500);
    test->SetBinSize (1);
    
    TString OutDir = "MinMult_";
    OutDir += MultMin;
    
    test->SetFitMultMin (MultMin);
    test->SetNormMultMin (MultMin);
    test->SetOutDirName (OutDir);
    
    test->SetSimHistos ("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Glauber/examples/glau_pbpb_ntuple_signn_31.0_7.6AGeV_CM_30AGeV_LC.root", nEvents);
    test->SetInputHisto (h1d_Mult_TPC_Ref);
       
    test->FitGlauber(nf, f0, f1, nsigma, StepSigma, nEvents);
    test->DrawHistos(true, true, true, false);

/*    test->SetGlauberFitHisto (1,0.879211, 145, 500000, true);
    TH1F *h1 = test->GetGlauberFitHisto ();
    h1->SetLineColor(kRed);
    
    std::cout << "Integral = "  << h1->Integral(0,450) << std::endl;
    
    h1->Draw();
    hData1->Draw("same");
    
    TLegend* legData = new TLegend(0.6,0.75,0.75,0.83);
    legData->AddEntry(h1 ,"Fit", "l");    
    legData->AddEntry(hData1 ,"M_{track}", "l");    
    legData->Draw("same"); */              
}
