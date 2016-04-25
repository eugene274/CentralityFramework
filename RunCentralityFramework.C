
void RunCentralityFramework(Int_t RunId = 23005)
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");  
    
    TString na61datadir = "/lustre/nyx/cbm/users/vblinov/NA61/QAtree/QA_15_12_09/";
    TString dir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";    
    
    TString det1 = "PSD1";
    TString det2 = "PSD2";
    TString det3 = "PSD3";
    TString det4 = "TPC";
    
    TString ContainerFile = dir + "root_files/" + Form ("na61_container_%d.root", RunId);
    
    TCut cuts = "";//"det1 > 0.65 - det2";
    
    CentralityManager *manager = new CentralityManager;
    manager->SetDirectory(dir);
    manager->SetRunId (RunId);
     
//  For CentralityFinder 
//  ************************************   

    manager->AddDetector("PSD1");
    manager->AddDetector("PSD2");
    manager->AddDetector("PSD3");
    manager->AddDetector("TPC");  
    manager->SetContainerFileName ( ContainerFile );  
    manager->CopyNa61ExpDataToContainer(na61datadir);

    manager->IsSimData(false);
    manager->Det1IsInt(true);
    manager->Do1DAnalisys(false);
    manager->SetDetectorsForCentralityAnalisys (det4, det1);
    manager->SetCentralityMax(100);
    manager->SetDirectionCentralEvents(0);
    manager->SetSliceStep (5);
    manager->SetCuts (cuts);
    
    manager->RunSliceFinder();
    manager->WriteCentralityFile();
    
//  
//     
    manager->Det1IsInt(false);
    manager->Do1DAnalisys(false);
    manager->SetDetectorsForCentralityAnalisys (det1, det2);
    manager->SetCentralityMax(100);
    manager->SetDirectionCentralEvents(1);
    manager->SetSliceStep (5);
    manager->SetCuts (cuts);
    
    manager->RunSliceFinder();
    manager->WriteCentralityFile();



    
//  For CentralityGetter
//  ************************************   
    
    float c = -1;
    manager->LoadCentalityDataFile( dir + "root_files/" + Form("Slices_%s_%s_%d.root", det4.Data(), det1.Data(), RunId) );
    c = manager->GetCentrality (100, 2000);
    std::cout << "Centrality = " << c << std::endl;
        
        
}






