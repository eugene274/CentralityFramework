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
//     TString ContainerFile = dir + "root_files/SC_ON_SL_OFF.root";
//     TString ContainerFile = dir + "root_files/urqmd_sis100_electron_SC_OFF_SL_OFF.root";
    TString ContainerFile = dir + "root_files/DCM-QGSM/sis100_STS_PSD_SC_OFF_SL_OFF_2016_04_11.root";

    TCut cuts = "det1 > 0.65 - det2";
    
    
    CentralityManager *manager = new CentralityManager;
    manager->SetDirectory(dir);
    
    manager->AddDetector("PSD1");
    manager->AddDetector("PSD2");
    manager->AddDetector("PSD3");
    manager->AddDetector("Mult");

    manager->SetRunId (RunId);
    manager->SetContainerFileName ( ContainerFile );   //output for NA61 data 

//     manager->CopyNa61ExpDataToContainer(na61datadir);
    manager->IsSimData(true);
    manager->Det1IsInt(true);
    manager->Do1DAnalisys(false);
    manager->SetDetectorsForCentralityAnalisys ("Mult", "PSD1");
    manager->SetCentralityMax(60);
    manager->SetDirectionCentralEvents(1);
    manager->SetSliceStep (5);
    manager->SetCuts (cuts);
    
    manager->RunSliceFinder();

    manager->WriteCentralityFile();
    
    
//     manager->SetDetectorsForCentralityAnalisys ("PSD1", "PSD2");
//     manager->Do1DAnalisys(false);
//     manager->SetDirectionCentralEvents(0);
//     manager->RunSliceFinder();
// 
//     manager->WriteCentralityFile();

    
    
//     float c = -1;
//     manager->LoadCentalityDataFile( dir + "root_files/Slices_Mult_PSD1_0.root");
//     c = manager->GetCentrality (250, 25);
//     std::cout << "Centrality = " << c << std::endl;
// // 
//     delete manager;    
// 
    return 0;     
}





