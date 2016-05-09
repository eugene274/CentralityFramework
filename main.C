#include "CentralityGetter.h"
#include "CentralitySlicesFinder.h"
#include "NA61DataEvent.h"
#include "CentralityManager.h"

int main()
{
    Int_t RunId = 23006;
    TString na61datadir = "/lustre/nyx/cbm/users/vblinov/NA61/QAtree/QA_15_12_09/";
    TString dir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";    
    TString ContainerFile = dir + "root_files/na61_test.root";
    
    TCut cuts = "";//"det1 > 0.65 - det2";
    
    CentralityManager *manager = new CentralityManager;
    manager->SetDirectory(dir);
     
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
    manager->SetDetectorsForCentralityAnalisys ("TPC", "PSD1");
    manager->SetCentralityMax(60);
    manager->SetDirectionCentralEvents(1);
    manager->SetSliceStep (5);
    manager->SetCuts (cuts);
    
    manager->RunSliceFinder();

    manager->WriteCentralityFile();
//     
//     
// //  For CentralityGetter
// //  ************************************   
//     
//     float c = -1;
//     manager->LoadCentalityDataFile( dir + "root_files/Slices_TPC_PSD1_23006.root");
//     c = manager->GetCentrality (250, 25);
//     std::cout << "Centrality = " << c << std::endl;
        
    
    return 0;     
}





