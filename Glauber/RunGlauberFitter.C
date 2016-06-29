

void RunGlauberFitter (Int_t nf = 2, Float_t f0 = 0.7, Float_t f1 = 0.8, Int_t nsigma = 3, Float_t StepSigma = 0.15, Int_t nEvents = 100000, Int_t MultMin = 20)
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";
//     gStyle->SetOptStat(0000);    
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");    
    
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //     
//     TH1F *h1d_Mult_TPC_Ref/* = new TH1F("h1d_Mult_TPC_Ref","M_{TPC}^{ref}", 125, 0, 250)*/;
//     TFile *f = new TFile("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Macro/hTPC_ref.root");    
//     h1d_Mult_TPC_Ref = (TH1F*)f->Get("h1Corr");
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
    
    TString DataFileName = CentralityFrameworkDir + "root_files/DCM-QGSM/" + "sis100_electron_SC_ON_SL_OFF.root";
    TFile *DataFile = new TFile ( DataFileName, "read" );    
    CentralityEventContainer *container = new CentralityEventContainer;
    TTree *ContTree = (TTree*)DataFile->Get("cbm_data");
    ContTree->SetBranchAddress("CentralityEventContainer", &container);
    ContTree->Draw("CentralityEventContainer.GetDetectorWeight(3) >> h1(400, 0, 400)");
    TH1F *hData1 = (TH1F*)gPad->GetPrimitive("h1");
    
    
    
    GlauberFitter *test = new GlauberFitter;

    test->SetNpartMax (450);
    test->SetNcollMax (1100);
    test->SetMultMax (400);
    test->SetBinSize (1);
    
    TString OutDir = "MinMult_";
    OutDir += MultMin;
    OutDir = "";
    
    test->SetFitMultMin (MultMin);
    test->SetNormMultMin (MultMin);
    test->SetOutDirName ("");
    
    test->SetSimHistos ("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Glauber/examples/glau_pbpb_ntuple_signn_31.0_7.6AGeV_CM_30AGeV_LC.root", nEvents);
    test->SetInputHisto (hData1);
       
//     f = 0.783333    mu = 0.669471    k = 11.6558
// f = 0.775    mu = 0.669559    k = 3.81179

    test->SetGlauberFitHisto (0.775, 0.669559, 3.81179, 500000, true);
    TH1F *h1 = test->GetGlauberFitHisto ();
    h1->Draw();
    hData1->Draw("same");
    
    TLegend* legData = new TLegend(0.6,0.75,0.75,0.83);
    legData->AddEntry(h1 ,"Fit", "l");    
    legData->AddEntry(hData1 ,"M_{STS}", "l");    
    legData->Draw("same");           

// f = 0.75    mu = 0.404777    k = 66.8897    Chi2Min = 4.00233
// f = 0.783333    mu = 0.692843    k = 11.6558    Chi2Min = 12.7801    
    
    
//     test->FitGlauber(nf, f0, f1, nsigma, StepSigma, nEvents);
//     test->DrawHistos(true, true, false, false);
}
























//     TString DataFileName = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/root_files/DCM-QGSM/sis100_electron_SC_ON_SL_OFF.root";
//     
//     TFile *DataFile = new TFile ( DataFileName, "read" );    
// 
//     CentralityEventContainer *container = new CentralityEventContainer;
//     TTree *ContTree = (TTree*)DataFile->Get("cbm_data");
//     ContTree->SetBranchAddress("CentralityEventContainer", &container);
//     
//     Int_t n = ContTree->GetEntries();
//     std::cout << "n = " << n  << std::endl;
// 
//     TCut psd1_mult = "";//" CentralityEventContainer.GetDetectorWeight(3) > 250 - 250.0/40 * CentralityEventContainer.GetDetectorWeight(0) ";
//     
//     TCanvas *c1 = new TCanvas("c1", "canvas", 1000, 800);
// //     c1->Divide(3, 3);
//     
// //     c1->cd(1);
//     ContTree->Draw("CentralityEventContainer.GetDetectorWeight(3) >> h1(100, 0, 400)", psd1_mult);
//     TH1F *hData1 = (TH1F*)gPad->GetPrimitive("h1");
//     hData1->GetXaxis()->SetTitle( "M_{STS}" );
//     hData1->GetYaxis()->SetTitle("Counts");   
// 
//     hData1->Draw();