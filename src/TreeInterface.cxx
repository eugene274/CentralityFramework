#include "TreeInterface.h"

ClassImp(TreeInterface)

// -----   Default constructor   -------------------------------------------
TreeInterface::TreeInterface() 
  : TNamed(),
    fEvent (new DataTreeEvent)
{
}

void TreeInterface::Test ()
{
    
    
    fContainer = new CentralityEventContainer;
    
    CentralityDetectorEvent psd; //= new CentralityDetectorEvent;
    CentralityDetectorEvent tpc ; //= new CentralityDetectorEvent;
    
    if (true) { fContainer->AddDetector(psd); }
    if (true) { fContainer->AddDetector(psd); }
    if (true) { fContainer->AddDetector(psd); }
    if (true) { fContainer->AddDetector(tpc);  }
    
    fOutTree = new TTree ( "cbm_data", "cbm_data" );
    fOutTree->Branch("CentralityEventContainer", "CentralityEventContainer", &fContainer);
    
    fOutTree->SetDirectory(0);  

    
    
    
    
    TString fInfileName;
    fInfileName = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/InputTree/test_tree.root";
    
    TFile *fInFile;
    fInFile = new TFile (fInfileName);
    
    TTree *fInTree;
    fInTree = (TTree*) fInFile->Get("fTreeQA");
    
    fInTree->SetBranchAddress("DTEvent", &fEvent);
    
    Int_t nEntries = fInTree->GetEntries();
    
    DataTreePSDModule *PdsMod = new DataTreePSDModule;
    
    for (Int_t iEntry = 0; iEntry<1/*nEntries*/; iEntry++)
    {
        fInTree->GetEntry (iEntry);
        Double_t B = fEvent->GetImpactParameter();
                
        std::vector <float> psd1 = SetPsdVector(1);
        std::vector <float> psd2 = SetPsdVector(2);
        std::vector <float> psd3 = SetPsdVector(3);
        std::vector <float> sts (1,1);

        fContainer->AddDetectorEvent (0, psd1);
        fContainer->AddDetectorEvent (1, psd2);
        fContainer->AddDetectorEvent (2, psd3);        
        fContainer->AddDetectorEvent (3, sts);
        fContainer->SetRunId (0);
        fContainer->SetB (B);

        fOutTree->Fill();


        for (Int_t i=0; i<psd1.size(); i++)
            std::cout << psd1.at(i) << ", ";
        std::cout << std::endl;

        for (Int_t i=0; i<psd2.size(); i++)
            std::cout << psd2.at(i) << ", ";
        std::cout << std::endl;

        for (Int_t i=0; i<psd3.size(); i++)
            std::cout << psd3.at(i) << ", ";
        std::cout << std::endl;

//         for (Int_t iMod = 1; iMod<=44; iMod++){
//             PdsMod = fEvent->GetPSDModule(iMod-1);
//             Double_t energy = PdsMod->GetEnergy();
//             
//             std::cout << "ModuleId = " << iMod << "    Energy = " << energy << std::endl;
//             
//         }        
        
        

        std::cout << "EventId = " << iEntry << "    Impact Parameter = " << B << std::endl;
        
        
    }

    
}

std::vector <Float_t> TreeInterface::SetPsdVector(Int_t subgroup)
{

    std::vector <Float_t> psdEnergies;
    DataTreePSDModule *PdsMod = new DataTreePSDModule;
    
    if  (true/*fPsdGeomConfig == "CBM_44"*/){   //NOTE later can be used for NA61 ?
        
        std::vector <Int_t> PsdPos;
        if ( subgroup == 1 )
            PsdPos = {18, 19, 26, 27};
        else if ( subgroup == 2 )
            PsdPos = {9, 10, 11, 12, 17, 20, 33, 34, 35, 36, 25, 28};
        else if ( subgroup == 3 )
            PsdPos = {1,2,3,4,5,6,7,8,13,14,15,16,21,22,23,24,29,30,31,32,37,38,39,40,41,42,43,44};            
 
        for (Int_t i=0; i<PsdPos.size(); i++)
        {
            PdsMod = fEvent->GetPSDModule(PsdPos.at(i)-1);
            Float_t edep = 0;
            if ( PdsMod ) edep = PdsMod->GetEnergy();
            psdEnergies.push_back(edep);
        }
    }
    
//     else if (fPsdGeomConfig == "NA61"){
//         
//         std::vector <Int_t> PsdPos;
//         if ( fDetName == "psd1" || fDetName == "PSD1" )
//             PsdPos = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,45};
//         else if ( fDetName == "psd2" || fDetName == "PSD2" )
//             PsdPos = {22,23,24,25,28,29,32,33,36,37,38,39};
//         else if ( fDetName == "psd3" || fDetName == "PSD3" )
//             PsdPos = {17,18,19,20,21,30,31,34,35,40,41,42,43,44};            
// 
//         for (Int_t i=0; i<PsdPos.size(); i++)
//         {
//             hit = (CbmPsdHit*) flistPSDhit->At(PsdPos.at(i)-1);
//             Float_t edep = 0;
//             if ( hit ) edep = hit->GetEdep();
//             psdEnergies.push_back(edep);
//         }
//     }    
//     std::cout << "psd1Energy = " << psd1Energy << std::endl;    
    
    return psdEnergies;
}
