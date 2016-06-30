#ifndef DataTreeConverterNA61_H
#define DataTreeConverterNA61_H 1

#include <iostream>

#include "TClonesArray.h"

#include "TFile.h"

#include "TString.h"
#include "TChain.h"
#include "TTree.h"
#include "DataTreeEvent.h"
#include "DataTreeTrack.h"

class TClonesArray;

class TDirectory;



class DataTreeConverterNA61
{
    
public:
  
    DataTreeConverterNA61();
    ~DataTreeConverterNA61();
    
    void SetInputFile(TString filename) { sInputFileName = filename; }
    void SetOutputFile(TString filename) { sOutputFileName = filename; }
    void SetNEvents(int fValue) { nEvents = fValue; }
    void Body();
      virtual void Init();
	void Chain_Init();
	void OutputTree_Init();
	void DataTreeEvent_Init();
      virtual void Exec(int event);
	void Clear_Event();
	void Read_Event();
      virtual void Finish();
      
private:
    
    DataTreeEvent* DTEvent;
  
    int fCurEvent;
    int nEvents;
    int nRealEvents;
    
    static const int nInfinity = 10000000;
    static const int nUndefinedValue = -999;
    
    static const int nRawReco = 2;
    static const int nTriggers_Simple = 5;
    static const int nTriggers = 6;
    static const int nBPD = 3;
    static const int nBPDcomponents = 3;
    static const int nPSD_Modules = 45;
    static const int nPSD_Sections = 10;
    static const int nTPC_Tracks = 2000;
    
    TString sInputFileName;
    TChain* fChain;
    
    TString sOutputFileName;
    TFile* fTreeFile;
    TTree* fTreeQA;
    
    int Run_Id;
    int Event_Id;
    int Event_Timestamp;
    float RPAngle;
    float ImpactParameter;
    
    bool T1;
    bool T2;
    bool T4;
    
    float triggersADC[nRawReco][nTriggers];
    float BPD_Position[nRawReco][nBPD][nBPDcomponents];
    
    float FitVertexX;
    float FitVertexY;
    float FitVertexZ;
    int FitVertexQ;
    
    float PSD_module_Energy[nPSD_Modules];
    float PSD_section_Energy[nPSD_Modules][nPSD_Sections];
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
    
    float TPC_track_DCAtoVertex_X[nTPC_Tracks];
    float TPC_track_DCAtoVertex_Y[nTPC_Tracks];
    
    int TPC_track_nClusters_Total[nTPC_Tracks];
    int TPC_track_nClusters_TPCV1[nTPC_Tracks];
    int TPC_track_nClusters_TPCV2[nTPC_Tracks];
    int TPC_track_nClusters_TPCVmain[nTPC_Tracks];
    int TPC_track_nClusters_TPCVgap[nTPC_Tracks];
    
    int TPC_track_nClustersPotential_Total[nTPC_Tracks];
    int TPC_track_nClustersPotential_TPCV1[nTPC_Tracks];
    int TPC_track_nClustersPotential_TPCV2[nTPC_Tracks];
    int TPC_track_nClustersPotential_TPCVmain[nTPC_Tracks];
    int TPC_track_nClustersPotential_TPCVgap[nTPC_Tracks];
    
    ClassDefNV(DataTreeConverterNA61, 1)
//     ClassDef(DataTreeConverterNA61,1);
};

#endif