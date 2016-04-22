#define NA61DataEvent_cxx

#include <iostream>
#include <fstream>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMath.h>

#include "CentralityEventContainer.h"
#include "CentralityDetectorEvent.h"

#include "NA61DataEvent.h"

using std::cout;
using std::endl;
using std::flush;


void NA61DataEvent::Loop()
{

    if (fChain == 0) return;

//     Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
   
    TCut T4 = "T4==1";
    TCut FitVertexZ = "FitVertexZ > - 589  && FitVertexZ < - 584.5";
//     TCut FitVertexQ = "FitVertexQ == 0";
//     TCut psd1_tpc = "PSD_1_Energy > 7000 - 7/1.4*TPC_Multiplicity";
    
    TCut allCuts = T4 && FitVertexZ  ;
    
    fChain->Draw(">>elist", allCuts);
    TEventList *cutList = (TEventList*)gDirectory->Get("elist");
    int iCut = cutList->GetN();
    
    TFile *OutFile = new TFile (fOutFile.Data(), "recreate");    
    
    CentralityEventContainer *container = new CentralityEventContainer;
    
    CentralityDetectorEvent psd1; //= new CentralityDetectorEvent;
    CentralityDetectorEvent tpc ; //= new CentralityDetectorEvent;
    
    if (isPSD1) { container->AddDetector(psd1); }
    if (isPSD2) { container->AddDetector(psd1); }
    if (isPSD3) { container->AddDetector(psd1); }
    if (isTPC) { container->AddDetector(tpc);  }
    
    TTree *OutTree = new TTree ( "na61_data", "na61_data" );
    OutTree->Branch("CentralityEventContainer", "CentralityEventContainer", &container);
    
    cout << "Reading events from chain..." << endl;
    for (Long64_t jentry=0; jentry<iCut;jentry++) {
        
        if ((jentry+1) % 50 == 0) cout << "Event # " << jentry+1 << "... \r" << flush; 
        
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(cutList->GetEntry(jentry));   nbytes += nb;

//         cout << "Event_Id = " << Event_Id << endl; 

        int j=0;
        int nGoodTracks = 0;
        while (TPC_track_eta[j] != -999 && j<2000)
        {
            if (isRefMultTrack(j))
                nGoodTracks++ ;
            j++;
        }

        float PSD_Mult_coor = PSD_Energy-5717.11+23.8262*nGoodTracks; 
        if (fabs(PSD_Mult_coor)>1200) continue;    
        
        std::vector <float> psd1_sig;
        std::vector <float> psd2_sig;
        std::vector <float> psd3_sig;
        std::vector <float> tpc_sig;
        
        psd1_sig.push_back(PSD_1_Energy);
        psd2_sig.push_back(PSD_2_Energy);
        psd3_sig.push_back(PSD_3_Energy);
        tpc_sig.push_back(nGoodTracks);

//         cout << "nGoodTracks = " << nGoodTracks << endl; 
        
        container->AddDetectorEvent (PSD1id, psd1_sig);
        container->AddDetectorEvent (PSD2id, psd2_sig);
        container->AddDetectorEvent (PSD3id, psd3_sig);        
        container->AddDetectorEvent (TPCid, tpc_sig);
        container->SetRunId (Run_Id);
        
        OutTree->Fill();

    }
    cout << "Writing centrility container..." << endl;
    OutTree->Write();
    OutFile->Close();

    cout << "Done! Centrality container " << fOutFile <<  " created!" << endl;

    
    
}

bool NA61DataEvent::isGoodTrack(int iTrk)
{
    if (TPC_track_eta[iTrk] >4.7) return false;
    if (TPC_track_eta[iTrk] <1.) return false;
    if (TPC_track_nClusters_TPCV1[iTrk]>72) return false;
    if (TPC_track_nClusters_TPCV2[iTrk]>72) return false;
    if (fabs(TPC_track_nClusters_TPCVmain[iTrk])>89) return false;
    
    if (fabs(TPC_track_DCAtoVertex_X[iTrk])>3.) return false;
    if (fabs(TPC_track_DCAtoVertex_Y[iTrk])>3.) return false;
    if (sqrt(TMath::Power(TPC_track_DCAtoVertex_X[iTrk],2)+TMath::Power(TPC_track_DCAtoVertex_Y[iTrk],2))>3.) return false;
    //                         if (fabs(TPC_track_chi2[iTrk])>1000000000000) return false;
    if (TPC_track_nClustersPotential_Total[iTrk] ==0 ) return false;
    if (TPC_track_nClusters_Total[iTrk]/(TPC_track_nClustersPotential_Total[iTrk]+0.00001) > 1.) return false;
    if (TPC_track_nClusters_Total[iTrk]/(TPC_track_nClustersPotential_Total[iTrk]+0.00001) < 0.7) return false;
//     if (TPC_track_nClusters_Total[iTrk]==TPC_track_nClustersPotential_Total[iTrk]) return false;
    if
        (
            TPC_track_nClusters_TPCV1[iTrk] == 0 &&
            TPC_track_nClusters_TPCV2[iTrk] == 0 &&
            TPC_track_nClusters_TPCVmain[iTrk] == 0
        )  return false;
    bool reject = true;
    if(TPC_track_nClusters_TPCVmain[iTrk] > 10 )  reject=false;
    if(TPC_track_nClusters_TPCV1[iTrk] > 6 )  reject=false;
    if(TPC_track_nClusters_TPCV2[iTrk] > 10 )  reject=false;
    if (reject) return false;
    return true;
}

bool NA61DataEvent::isRefMultTrack(int iTrk)
{
    if (!isGoodTrack(iTrk)) return false;
    if (fabs(TPC_track_DCAtoVertex_X[iTrk])>2) return false;
    if (fabs(TPC_track_DCAtoVertex_Y[iTrk])>2) return false;
    if (sqrt(TMath::Power(TPC_track_DCAtoVertex_X[iTrk],2)+TMath::Power(TPC_track_DCAtoVertex_Y[iTrk],2))>2.) return false;
    if (TPC_track_nClusters_Total[iTrk]/(TPC_track_nClustersPotential_Total[iTrk]+0.00001) < 0.7) return false;
//     if (abs(TPC_track_eta[iTrk]-2.0792)>0.5) return false;
    return true;
}



