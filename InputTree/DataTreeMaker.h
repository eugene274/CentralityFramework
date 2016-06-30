#ifndef DataTreeMaker_H
#define DataTreeMaker_H 1

#include "CbmMCEventData.h"
#include "CbmPsdEventData.h"
#include "CbmStsEventData.h"
// #include "CbmMCEventHeader.h"
//#include "CbmAnaParticleData.h"

#include "KFMCParticle.h"
#include "KFParticleMatch.h"
#include "CbmKFPartEfficiencies.h"
#include "CbmKFParticleFinder.h"

#include "FairTask.h"
#include "CbmVertex.h"
#include <vector>
#include "TLorentzVector.h"
#include <map>
#include <cstring>

#include "UEvent.h"
#include "TGraphErrors.h"
#include <iostream>
// #include "QdetVector.h"
#include "TClonesArray.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "FairMCEventHeader.h"
#include "UEvent.h"

// #include "CbmL1PFFitter.h"
#include "CbmKFVertex.h"
#include "CbmTrackMatch.h"
// #include "L1FieldRegion.h"

#include "DataTreeEvent.h"
#include "DataTreeTrack.h"

// #include "QdetVectorGetter.h"
// #include "CbmL1.h"

// class CbmMCEventHeader;
// class CbmKFParticlesFinder;
class TClonesArray;
class CbmVertex;
class TDirectory;
class TH1F;
class TProfile;
class TH2F;
// class QdetVector;


class DataTreeMaker : public FairTask
{
    
const double SpeedOfLight = 29.9792458;

public:
  
    DataTreeMaker();
    ~DataTreeMaker();
    
    virtual InitStatus Init();
    virtual void Exec(Option_t* opt);
    virtual void Finish();
    
    void SetPSDCoordinatesFileName(TString fileName) { sPSDCoordinatesFileName = fileName; } // File containing PSD module (x, y) coordinates in LAB
    void SetOutputFile(TString filename) { sOutputFileName = filename; }
    void SetInputFile(TString filename) { sInputFileName = filename; }
    void SetCbmKFParticleFinder_MC(CbmKFParticleFinder* kf) {fCbmKFParticleFinder_MC=kf;}
    void SetCbmKFParticleFinder_TOF(CbmKFParticleFinder* kf) {fCbmKFParticleFinder_TOF=kf;}
    void SetGenerator(int fValue){Generator = fValue;}
    
    
private:
    
    DataTreeEvent* DTEvent;
  
    CbmKFParticleFinder* fCbmKFParticleFinder_TOF;
    CbmKFParticleFinder* fCbmKFParticleFinder_MC;
    
    int Generator = 0; //0 -- URQMD, 1 -- DCM_QGSM
    
    int fCurEvent;
//     CbmMCEventHeader* fHeader;
    CbmVertex* fPrimVtx;
    FairMCEventHeader* fHeader;
    TClonesArray* flistPSDhit;
    TClonesArray* flistPSDdigit;
    TClonesArray* flistMCtrack;
    TClonesArray* flistSTSRECOtrack;
    TClonesArray* flistSTStrackMATCH;
    TClonesArray *fGlobalTrackArray; //input reco tracks
    TClonesArray *fTofHitArray; //input reco tracks
        
    void OutputTree_Init();
    void DataTreeEvent_Init();
    void PSD_Init();
    void Clear_Event();
    void Read_Event();
    void Read_PSD();
    void Read_STS();
    void Read_TOF_Hits();
    void Read_MC_Tracks();
    void Read_Generated_Tracks();
    void Read_V0_Candidate(int UseMCpid);
    
    const int nInfinity = 10000000;
    const int nUndefinedValue = -999;
    static const int nV0Daughters = 2;
    
    static const int nPSD_Modules = 44;
    static const int nTPC_Tracks = 500;
    
    TFile* fTreeFile;
    TTree* fDataTree;
    TChain* fTreeEvents;//temp for bad phi data
    UEvent* uEvent;//temp for bad phi data
    TString sInputFileName;//temp for bad phi data
    
    TString sPSDCoordinatesFileName;
    TString sOutputFileName;
    
    int Run_Id;
    int Event_Id;
    float RPAngle;
    float ImpactParameter;
    
    float FitVertexX;
    float FitVertexY;
    float FitVertexZ;
    
    float PSD_module_Energy[nPSD_Modules];
    float PSD_module_X[nPSD_Modules];
    float PSD_module_Y[nPSD_Modules];
    float PSD_module_Z[nPSD_Modules];
    
    float TPC_track_pT[nTPC_Tracks];
    float TPC_track_eta[nTPC_Tracks];
    float TPC_track_phi[nTPC_Tracks];
    
    int TPC_track_NofHits[nTPC_Tracks];
    int TPC_track_PidHypo[nTPC_Tracks];
    int TPC_track_Flag[nTPC_Tracks];
    double TPC_track_ChiSq[nTPC_Tracks];
    int TPC_track_NDF[nTPC_Tracks];
    
    float Sim_TPC_track_pT[nTPC_Tracks];
    float Sim_TPC_track_eta[nTPC_Tracks];
    float Sim_TPC_track_phi[nTPC_Tracks];
    
    ClassDefNV(DataTreeMaker, 1)
//     ClassDef(DataTreeMaker,1);
};

#endif