#include "TreeInterface.h"
#include "../config/ContConfig.C"

ClassImp(TreeInterface)

using std::cout;
using std::endl;
using std::flush;


// -----   Default constructor   -------------------------------------------
TreeInterface::TreeInterface() 
  : TNamed(),
    fEvent (new DataTreeEvent),
    fPsdGeomConfig ("CBM_44"),
    nEntries (0)
{
}

void TreeInterface::WriteCentralityContainer ()
{
    
    TFile *fInFile = new TFile (fInFileName);
    TTree *fInTree = (TTree*) fInFile->Get("fDataTree"); 
    fInTree->SetBranchAddress("DTEvent", &fEvent);

    fContainer = new CentralityEventContainer;
    
    CentralityDetectorEvent psd; 
    CentralityDetectorEvent tpc ; 
    
    if (fPsdGeomConfig == "CBM_44"){
    
        if (true) { fContainer->AddDetector(psd); }
        if (true) { fContainer->AddDetector(psd); }
        if (true) { fContainer->AddDetector(psd); }
        if (true) { fContainer->AddDetector(tpc);  }
    }
    
    else if (fPsdGeomConfig == "NA61"){
    
        if (true) { fContainer->AddDetector(psd); }
        if (true) { fContainer->AddDetector(psd); }
        if (true) { fContainer->AddDetector(psd); }
        if (true) { fContainer->AddDetector(psd); }
        if (true) { fContainer->AddDetector(psd); }
        if (true) { fContainer->AddDetector(psd); }
        if (true) { fContainer->AddDetector(tpc);  }
    } 
            
    fOutFile = new TFile (fOutFileName.Data(), "recreate");    
    fOutTree = new TTree ( "container", "container" );
    fOutTree->Branch("CentralityEventContainer", "CentralityEventContainer", &fContainer);
    fOutTree->SetDirectory(0);  
    
    if (!nEntries) nEntries = fInTree->GetEntries();
    
    cout << "Reading " << nEntries << " events from chain..." << endl;
    for (Int_t iEntry = 0; iEntry<nEntries; iEntry++)
    {
        
        if ((iEntry+1) % 100 == 0) cout << "Event # " << iEntry+1 << "... \r" << flush; 
        
        fInTree->GetEntry (iEntry);
        
        if (!isSelectedEvent (fEvent)) continue;
        
        int nTracks_Ref = GetRefMultiplicity(fEvent);
                
        std::vector <float> psd1 = SetPsdVector(1);
        std::vector <float> psd2 = SetPsdVector(2);
        std::vector <float> psd3 = SetPsdVector(3);
        std::vector <float> sts (1, nTracks_Ref);

        fContainer->AddDetectorEvent (0, sts);
        fContainer->AddDetectorEvent (1, psd1);
        fContainer->AddDetectorEvent (2, psd2);        
        fContainer->AddDetectorEvent (3, psd3);        

        fContainer->SetRunId (fEvent->GetRunId());
        fContainer->SetEventId(fEvent->GetEventId());
        fContainer->SetB (fEvent->GetImpactParameter());

        fOutTree->Fill();        
        
    }

    cout << "Writing centrility container..." << endl;
    fOutTree->Write();
    fOutFile->Close();
    cout << "Done! Centrality container " <<  fOutFileName <<  " created!" << endl;
    
}

std::vector <Float_t> TreeInterface::SetPsdVector(Int_t subgroup)
{

    std::vector <Float_t> psdEnergies;
    DataTreePSDModule *PdsMod = new DataTreePSDModule;
    
    if  (fPsdGeomConfig == "CBM_44"){   //NOTE later can be used for NA61 ?
        
        std::vector <Int_t> PsdPos;
        if ( subgroup == 1 )
            PsdPos = {18, 19, 26, 27};
        else if ( subgroup == 2 )
            PsdPos = {9, 10, 11, 12, 17, 20, 33, 34, 35, 36, 25, 28};
        else if ( subgroup == 3 )
            PsdPos = {1,2,3,4,5,6,7,8,13,14,15,16,21,22,23,24,29,30,31,32,37,38,39,40,41,42,43,44};            
        else 
            std::cout << "TreeInterface::SetPsdVector:  Wrong PSD subgroup identificator!" << std::endl;
        
        
        for (UInt_t i=0; i<PsdPos.size(); i++)
        {
            PdsMod = fEvent->GetPSDModule(PsdPos.at(i)-1);
            Float_t edep = 0;
            if ( PdsMod ) edep = PdsMod->GetEnergy();
            psdEnergies.push_back(edep);
        }
    }
    
    else if (fPsdGeomConfig == "NA61"){
        
        std::vector <Int_t> PsdPos;
        if ( subgroup == 1 )
            PsdPos = {6,7,10,11};
        else if ( subgroup == 2 )
            PsdPos = {1,2,3,4,5,8,9,12,13,14,15,16};
        else if ( subgroup == 12 )
            PsdPos = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,45};
        else if ( subgroup == 3 )
            PsdPos = {22,23,24,25,28,29,32,33,36,37,38,39};
        else if ( subgroup == 4 )
            PsdPos = {17,18,19,20,21,30,31,34,35,40,41,42,43,44};            
        else if ( subgroup == 1234 )
            PsdPos = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45};            
        else 
            std::cout << "TreeInterface::SetPsdVector:  Wrong PSD subgroup identificator!" << std::endl;
        
        for (UInt_t i=0; i<PsdPos.size(); i++)
        {
            PdsMod = fEvent->GetPSDModule(PsdPos.at(i)-1);
            Float_t edep = 0;
            if ( PdsMod ) edep = PdsMod->GetEnergy();
            psdEnergies.push_back(edep);
        }
    }    
    
    return psdEnergies;
}



