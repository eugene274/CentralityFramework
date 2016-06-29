
void RunCentralityFramework()
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");  
    
    TString dir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";    
    
    TString det1 = "E_{PSD}^{1}";
    TString det2 = "E_{PSD}^{2}";
    TString det3 = "E_{PSD}^{3}";
    TString det4 = "M_{STS}";
    
//     TString ContainerFile = dir + "root_files/" + "na61_container_merged_nocuts.root";
//     TString ContainerFile = dir + "root_files/" + "na61_container_merged.root";
//     TString ContainerFile = dir + "root_files/" + "na61_container_22912.root";
    
    
    
    TString ContainerFile = dir + "root_files/DCM-QGSM/" + "sis100_STS_PSD_SC_OFF_SL_OFF_2016_04_11.root";
//     TString ContainerFile = dir + "root_files/DCM-QGSM/" + "sis100_electron_SC_OFF_SL_OFF_urqmd.root";
    
    TCut cuts = "det1 > 0.65 - det2";
    
    CentralityManager *manager = new CentralityManager;
    manager->SetDirectory(dir);
     
//  For CentralityFinder 
//  ************************************   
    manager->AddDetector(det1);
    manager->AddDetector(det2);
    manager->AddDetector(det3);
    manager->AddDetector(det4);  

//     manager->AddDetector("PSD4");
//     manager->AddDetector("PSD12");
//     manager->AddDetector("PSD1234");


    manager->SetContainerFileName ( ContainerFile );  
//     manager->CopyNa61ExpDataToContainer(na61datadir);
//     manager->SetNormalization (733440);
    manager->IsSimData(true);
    manager->Det1IsInt(false);
    manager->Do1DAnalisys(true);
    manager->SetDetectorsForCentralityAnalisys (det1, det4);
    manager->SetCentralityMax(60);
    manager->SetDirectionCentralEvents(0);
    manager->SetSliceStep (5);
    manager->SetCuts (cuts);
    
    manager->RunSliceFinder();
    manager->WriteCentralityFile();
    manager->QA();


    
//     manager->IsSimData(true);
//     manager->Det1IsInt(true);
// //     manager->Do1DAnalisys(false);
//     manager->SetDetectorsForCentralityAnalisys (det4);
//     manager->SetCentralityMax(100);
//     manager->SetDirectionCentralEvents(1);
//     manager->SetSliceStep (5);
//     manager->SetCuts (cuts);
//     
//     manager->RunSliceFinder();
//     manager->WriteCentralityFile();
//     manager->QA();


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






