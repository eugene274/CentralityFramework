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
    fPsdGeomConfig ("CBM_44")
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
    
//     if (fPsdGeomConfig == "CBM_44"){
    
        if (true) { fContainer->AddDetector(psd); }
        if (true) { fContainer->AddDetector(psd); }
        if (true) { fContainer->AddDetector(psd); }
//         if (true) { fContainer->AddDetector(psd); }
//         if (true) { fContainer->AddDetector(psd); }
//         if (true) { fContainer->AddDetector(psd); }
        if (true) { fContainer->AddDetector(tpc);  }
//     }
    
    fOutFile = new TFile (fOutFileName.Data(), "recreate");    
    fOutTree = new TTree ( "container", "container" );
    fOutTree->Branch("CentralityEventContainer", "CentralityEventContainer", &fContainer);
    fOutTree->SetDirectory(0);  
    
    Int_t nEntries = fInTree->GetEntries();
    
    cout << "Reading " << nEntries << " events from chain..." << endl;
    for (Int_t iEntry = 0; iEntry<nEntries; iEntry++)
    {
        
        if ((iEntry+1) % 50 == 0) cout << "Event # " << iEntry+1 << "... \r" << flush; 
        
        fInTree->GetEntry (iEntry);
        
        if (!isSelectedEvent (fEvent)) continue;
        
        int nTracks_Ref = ContConfig(fEvent);
                
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
    cout << "Done! Centrality container " <<  fOutFileName <<  " created!" << endl;}

std::vector <Float_t> TreeInterface::SetPsdVector(Int_t subgroup)
{

    std::vector <Float_t> psdEnergies;
    DataTreePSDModule *PdsMod = new DataTreePSDModule;
    
    if  (fPsdGeomConfig == "CBM_44"){   //NOTE later can be used for CBM ?
        
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


bool TreeInterface::isRefMultTrack(int iTrk)
{
    Float_t Beam_eta = 2.0792;
    DataTreeTrack* track = fEvent->GetTrack(iTrk);
    if (track->GetEta(0) < Beam_eta - 0.5 || track->GetEta(0) > Beam_eta + 0.5) return false;
    if (track->GetNofHits(0,1) > 72 || track->GetNofHits(0,2) > 72 || track->GetNofHits(0,3) > 89) return false;
    if (track->GetDCAComponent(0,0) > 2 || track->GetDCAComponent(0,1) > 2) return false;
    if (TMath::Sqrt(TMath::Power(track->GetDCAComponent(0,0),2)+TMath::Power(track->GetDCAComponent(0,1),2)) > 2) return false;
    if (track->GetNofHits(0,1) < 6 && track->GetNofHits(0,2) < 10 && track->GetNofHits(0,3) < 6) return false;
    return true;
}

bool TreeInterface::isGoodEvent(Int_t nTPC_Tracks_Ref)
{
    Float_t fPSD_Energy_Total = fEvent->GetPSDEnergy();
    if (fEvent->GetVertexPositionComponent(2) > -584.5 || fEvent->GetVertexPositionComponent(2) < -589) return false;
    if (nTPC_Tracks_Ref < 0) return false;
//     if (fPSD_Energy_Total < 1200 || fPSD_Energy_Total > 6000) return false;   //TODO check
    if (fEvent->GetTrigger(7)->GetSignal() != 1) return false;
    if (fPSD_Energy_Total  < 4500 - nTPC_Tracks_Ref*(4500.-1200.)/180.) return false;  //TODO check
    return true;
}


