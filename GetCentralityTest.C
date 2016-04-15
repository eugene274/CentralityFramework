
int GetCentralityTest(TString DataFileName = "au10au_cbm_test.root")
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";
    TString CentralityFileName = "Slices_Mult_PSD1_0.root";
    
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");

    float c = -1;
    CentralityManager *manager = new CentralityManager;
//     manager->SetDirectory(CentralityFrameworkDir);
//     manager->LoadCentalityDataFile( CentralityFrameworkDir + "root_files/" + CentralityFileName);
//     
//     c = manager->GetCentrality (250, 25);
//     std::cout << "Test Centrality 250, 25 = " << c << std::endl;


    
    
    TFile *DataFile = new TFile ( CentralityFrameworkDir + "root_files/" + DataFileName, "read" );    

    CentralityEventContainer *container = new CentralityEventContainer;
    TTree *ContTree = (TTree*)DataFile->Get("cbm_data");
    ContTree->SetBranchAddress("CentralityEventContainer", &container);
    
    TH1F *hCentr = new TH1F ("hCentr", "", 101, -1, 100);
    TH2F *h2DCor = new TH2F ("h2DCor", "", 100, 0, 500, 100, 0, 50);

    Float_t PSD1, PSD2, PSD3, Msts;
    Float_t Centrality;
    
    Int_t n = ContTree->GetEntries();
    std::cout << "n = " << n  << std::endl;

    TCanvas *c1 = new TCanvas("c1", "canvas", 1500, 800);
//     c1->Divide(2,1);
    
//     c1->cd(1);
    ContTree->Draw("CentralityEventContainer.GetDetectorWeight(1) : CentralityEventContainer.GetDetectorWeight(3) >> h1(100, 0, 400, 100, 0, 70)", "", "colz");
    TH2F *hData = (TH2F*)gPad->GetPrimitive("h1");
    hData->GetXaxis()->SetTitle( "M_{STS}" );
    hData->GetYaxis()->SetTitle("E_{PSD}^{1}, GeV");    
    
//     for (Int_t i=0; i<n; i++)
//     {
//         ContTree->GetEntry(i);
//         PSD1 = container->GetDetectorWeight(0);
//         PSD2 = container->GetDetectorWeight(1);
//         PSD3 = container->GetDetectorWeight(2);
//         Msts = container->GetDetectorWeight(3);
//         
//         Centrality = manager->GetCentrality (Msts, PSD1);
//         
//         if (Centrality > 00 && Centrality < 80)
//         {
// //             std::cout << "PSD1 = " << PSD1 << "  Msts = " << Msts << std::endl;
//             h2DCor->Fill (Msts, PSD1);
//             
//         }
//         hCentr->Fill(Centrality);
//         
// //         std::cout << "PSD1 = " << container->GetDetectorWeight(0) << std::endl;
//     }
    
//     hCentr->Draw();
//     c1->cd(2);
//     h2DCor->Draw("colz");
    
    
    return 0;     
}
