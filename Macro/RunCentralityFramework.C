
void RunCentralityFramework()
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");  
    
    TString dir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";    
    
    TString det1 = "PSD1";
    TString det2 = "PSD2";
    TString det3 = "PSD3";
    TString det4 = "TPC";
    
//     TString ContainerFile = dir + "root_files/" + "na61_container_merged_nocuts.root";
    TString ContainerFile = dir + "root_files/" + "na61_container_merged.root";
    
    TCut cuts = "";//"det1 > 0.65 - det2";
    
    CentralityManager *manager = new CentralityManager;
    manager->SetDirectory(dir);
     
//  For CentralityFinder 
//  ************************************   
    manager->AddDetector("TPC");  
    manager->AddDetector("PSD1");
    manager->AddDetector("PSD2");
    manager->AddDetector("PSD3");
    manager->AddDetector("PSD4");
    manager->AddDetector("PSD12");
    manager->AddDetector("PSD1234");


    manager->SetContainerFileName ( ContainerFile );  
//     manager->CopyNa61ExpDataToContainer(na61datadir);
    manager->SetNormalization (733440);
    manager->IsSimData(false);
    manager->Det1IsInt(true);
//     manager->Do1DAnalisys(false);
    manager->SetDetectorsForCentralityAnalisys ("TPC");
    manager->SetCentralityMax(80);
    manager->SetDirectionCentralEvents(1);
    manager->SetSliceStep (5);
    manager->SetCuts (cuts);
    
//     manager->RunSliceFinder();
//     manager->WriteCentralityFile();
    
//  
//     
//     manager->Det1IsInt(false);
//     manager->Do1DAnalisys(false);
//     manager->SetDetectorsForCentralityAnalisys (det1, det2);
//     manager->SetCentralityMax(100);
//     manager->SetDirectionCentralEvents(0);
//     manager->SetSliceStep (5);
//     manager->SetCuts (cuts);
//     
//     manager->RunSliceFinder();
//     manager->WriteCentralityFile();



/*    
//  For CentralityGetter
//  ************************************   
    
    det4 = "TPC";
    det1 = "PSD12";
    
    float c = -1;
    manager->LoadCentalityDataFile( dir + "root_files/" + Form("Slices_%s_%s.root", det4.Data(), det1.Data()) );
    
    
    TString DataFileName = ContainerFile;
    TFile *DataFile = new TFile ( ContainerFile, "read" );    

    CentralityEventContainer *container = new CentralityEventContainer;
    TTree *ContTree = (TTree*)DataFile->Get("na61_data");
    ContTree->SetBranchAddress("CentralityEventContainer", &container);
    
    TH1F *hCentr = new TH1F ("hCentr", "", 101, -1, 100);

    Float_t PSD1, PSD2, PSD3, PSD12, M;
    Float_t Centrality;
    
    Int_t n = ContTree->GetEntries();
    std::cout << "n = " << n  << std::endl;

    TCanvas *c1 = new TCanvas("c1", "canvas", 1500, 800);

    Int_t RunId = 0;
    for (Int_t i=0; i<n; i++)
    {
        ContTree->GetEntry(i);
        if (RunId != container->GetRunId())
        {
            RunId = container->GetRunId();
            manager->SetGetterRunId (RunId);
        }
        
        M = container->GetDetectorWeight(0);
        PSD1 = container->GetDetectorWeight(1);
        PSD2 = container->GetDetectorWeight(2);
        PSD3 = container->GetDetectorWeight(3);
        PSD12 = container->GetDetectorWeight(5);
//         std::cout << "M = " << M << "   PSD12 = " << PSD12  << std::endl;
        Centrality = manager->GetCentrality (M, PSD12);
        hCentr->Fill(Centrality);
        
    }
    
    hCentr->Draw();
        */
}






