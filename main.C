#include "CentralityGetter.h"
#include "CentralitySlicesFinder.h"
#include "NA61DataEvent.h"
#include "CentralityManager.h"

int main()
{

    Int_t RunId = 0;
   
//     TString na61datadir = "/lustre/nyx/cbm/users/vblinov/NA61/QAtree/QA_15_12_09/";
//     TString InFileName = "";
    TString dir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";    
    
    TString ContainerFile = dir + "root_files/cbm_shield_1e4ev_au10au.root";

    
    CentralityManager *manager = new CentralityManager;
    manager->SetDirectory(dir);
    
    manager->AddDetector("PSD1");
    manager->AddDetector("PSD2");
    manager->AddDetector("PSD3");
    manager->AddDetector("Mult");

    manager->SetRunId (RunId);
    manager->SetContainerFileName ( ContainerFile );   //output for NA61 data 

//     manager->CopyNa61ExpDataToContainer(na61datadir);

    manager->SetDetectorsForCentralityAnalisys ("Mult", "PSD1");
    manager->SetCentralityMax(70);
    manager->SetDirectionCentralEvents(1);
//     manager->Det1IsInt(true);
    manager->Do1DAnalisys(false);
    manager->SetSliceStep (5);
    manager->RunSliceFinder();

    manager->WriteCentralityFile();
    
    
//     manager->SetDetectorsForCentralityAnalisys ("PSD1", "PSD2");
//     manager->Do1DAnalisys(false);
//     manager->SetDirectionCentralEvents(0);
//     manager->RunSliceFinder();
// 
//     manager->WriteCentralityFile();

    
    
    float c = -1;
    manager->LoadCentalityDataFile( dir + "root_files/Slices_Mult_PSD1_0.root");
    c = manager->GetCentrality (250, 25);
    std::cout << "Centrality = " << c << std::endl;
// 
    delete manager;    

    return 0;     
}





