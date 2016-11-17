
void RunCentralityFramework()
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/git/CentralityFramework/";
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");  
    
    TString dir = "/lustre/nyx/cbm/users/klochkov/git/CentralityFramework/";    
    
    TString sts_det = "M_{STS}";
    TString psd1_det = "E_{PSD}^{1}";
    TString psd2_det = "E_{PSD}^{2}";
    TString psd3_det = "E_{PSD}^{3}";
    
//     TString ContainerFile = dir + "root_files/" + "na61_container_merged_nocuts.root";
//     TString ContainerFile = dir + "root_files/" + "na61_container_merged.root";
//     TString ContainerFile = dir + "root_files/na61_cont/" + "na61_container_merged.root";
    
//     TString ContainerFile = dir + "containers/" + "ana_dima_merged.root";
    TString ContainerFile = "/lustre/nyx/cbm/users/dblau/CentralityFramework.new/containers/cbm_urqmd_10AGeV_1.root";
    TCut cuts = "";//"det1 > 0.6 - 0.6/0.75*det2";
    
    CentralityManager *manager = new CentralityManager;
    manager->SetDirectory(dir);
     
//  For CentralityFinder 
//  ************************************   
    manager->AddDetector(sts_det );  
    manager->AddDetector(psd1_det);
    manager->AddDetector(psd2_det);
    manager->AddDetector(psd3_det);

    manager->SetContainerFileName ( ContainerFile );  
//     manager->SetNormalization (979057);                  // set full integral (calculated with Glauber model, for mc simulation we don't neeed it?)

    manager->IsSimData(true);                               // write impact parameter value - true for CBM (mc simulations)
    manager->Det1IsInt(true);                               // if multiplicity used as det1 - true; if as det2 - manager->Det2IsInt(true);
//     manager->Do1DAnalisys(true);                            // 1D slicing
    manager->SetDetectorsForCentralityAnalisys (sts_det, psd1_det);   // for 2D slicing give 2 arguments for function (for example sts_det, psd1_det)
    manager->SetCentralityMax(50);                         // set maximum value for centrality, more peripheral events will not be analysed
    manager->SetDirectionCentralEvents(1);                  // 1 - if impact parameter and detector signal is correlated = 1 (PSD1), if anticorelated (STS) = 0   
    manager->SetSliceStep (5);                              // slice step (width) in percents
    manager->SetCuts (cuts);                                // set additional cuts (obsolete probably, since we have config file for containers filling)
    
    manager->RunSliceFinder();
    manager->WriteCentralityFile();
    manager->QA();


    
//     manager->IsSimData(true);
//     manager->Det1IsInt(true);
//     manager->Do1DAnalisys(false);
//     manager->SetDetectorsForCentralityAnalisys (sts_det, psd1_det);
//     manager->SetCentralityMax(100);
//     manager->SetDirectionCentralEvents(1);
//     manager->SetSliceStep (5);
//     manager->SetCuts (cuts);
//     
//     manager->RunSliceFinder();
//     manager->WriteCentralityFile();
//     manager->QA();


}






