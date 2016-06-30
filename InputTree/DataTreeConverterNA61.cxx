#include "DataTreeConverterNA61.h"
#include <iostream>
#include <fstream>
#include <vector>
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;

#include "TDirectory.h"
#include "DataTreeEvent.h"
#include "DataTreeTrack.h"
#include "TString.h"
#include "TChain.h"
#include "TTree.h"

//=================================================================> MAIN <===============================================================
DataTreeConverterNA61::DataTreeConverterNA61()
:
  nEvents (1),
  fCurEvent (0)
{
      DTEvent = new DataTreeEvent();
}
//--------------------------------------------------------------------------------------------------
DataTreeConverterNA61::~DataTreeConverterNA61()
{
    
}
//--------------------------------------------------------------------------------------------------
void DataTreeConverterNA61::Body()
{
//     cout << "Init" << endl;
    Init();
    nRealEvents = fChain->GetEntries();
//     cout << nRealEvents << " " << nEvents << endl;//TEST
    if (fChain->GetEntries() > nEvents) nRealEvents = nEvents;
//     cout << nRealEvents << " " << nEvents << endl;//TEST
//     cout << "Exec" << endl;
    for (int i=0;i<nRealEvents;i++)
    {
	Exec(i);
    }
//     cout << "Finish" << endl;
    Finish();
}
//=================================================================> INIT <===============================================================
void DataTreeConverterNA61::Init()
{
//     cout << "Chain_Init" << endl;
    Chain_Init();
//     cout << "OutputTree_Init" << endl;
    OutputTree_Init();
//     cout << "DataTreeEvent_Init" << endl;
    DataTreeEvent_Init();
}
//--------------------------------------------------------------------------------------------------
void DataTreeConverterNA61::DataTreeEvent_Init()
{
    for (int i=0;i<nPSD_Modules;i++)
    {
	DTEvent -> AddPSDModule(10);
    }
    TString sTriggerNames[8] = {"S1","S2","S3","V1","PSD","T1","T2","T4"};	
    for (int i=0;i<8;i++)
    {
	DTEvent -> AddTrigger(sTriggerNames[i]);
    }
    for (int i=0;i<nBPD;i++)
    {
	DTEvent -> AddBPD();
    }
}
//--------------------------------------------------------------------------------------------------
void DataTreeConverterNA61::OutputTree_Init()
{
    fTreeFile = new TFile(sOutputFileName, "RECREATE");
    fTreeFile -> cd();
    fTreeQA = new TTree("fTreeQA","fTreeQA");
    fTreeQA -> SetMaxTreeSize(90000000);
    
    fTreeQA -> Branch("DTEvent", &DTEvent, 256000, 3);    
}
//--------------------------------------------------------------------------------------------------
void DataTreeConverterNA61::Chain_Init()
{
    fChain = new TChain("fTreeQA","fTreeQA");
    fChain -> Add(sInputFileName);
    
    fChain->SetBranchAddress("Run_Id",&Run_Id);
    fChain->SetBranchAddress("Event_Id",&Event_Id);
//     fChain->SetBranchAddress("RPAngle",&RPAngle);
//     fChain->SetBranchAddress("ImpactParameter",&ImpactParameter);
    
    fChain->SetBranchAddress("Event_Timestamp",&Event_Timestamp);
    
    fChain->SetBranchAddress("T1",&T1);
    fChain->SetBranchAddress("T2",&T2);
    fChain->SetBranchAddress("T4",&T4);
    
    fChain->SetBranchAddress("triggersADC",&triggersADC);
    
    fChain->SetBranchAddress("BPD_Position",&BPD_Position);
    
    fChain->SetBranchAddress("FitVertexX",&FitVertexX);
    fChain->SetBranchAddress("FitVertexY",&FitVertexY);
    fChain->SetBranchAddress("FitVertexZ",&FitVertexZ);
    fChain->SetBranchAddress("FitVertexQ",&FitVertexQ);
    
    fChain->SetBranchAddress("PSD_module_Energy",&PSD_module_Energy);
    fChain->SetBranchAddress("PSD_section_Energy",&PSD_section_Energy);
    fChain->SetBranchAddress("PSD_module_X",&PSD_module_X);
    fChain->SetBranchAddress("PSD_module_Y",&PSD_module_Y);
    fChain->SetBranchAddress("PSD_module_Z",&PSD_module_Z);
    
    fChain->SetBranchAddress("TPC_track_pT",&TPC_track_pT);
    fChain->SetBranchAddress("TPC_track_eta",&TPC_track_eta);
    fChain->SetBranchAddress("TPC_track_phi",&TPC_track_phi);

    fChain->SetBranchAddress("TPC_track_DCAtoVertex_X",&TPC_track_DCAtoVertex_X);
    fChain->SetBranchAddress("TPC_track_DCAtoVertex_Y",&TPC_track_DCAtoVertex_Y);
    
    fChain->SetBranchAddress("TPC_track_nClusters_Total",&TPC_track_nClusters_Total);
    fChain->SetBranchAddress("TPC_track_nClusters_TPCV1",&TPC_track_nClusters_TPCV1);
    fChain->SetBranchAddress("TPC_track_nClusters_TPCV2",&TPC_track_nClusters_TPCV2);
    fChain->SetBranchAddress("TPC_track_nClusters_TPCVmain",&TPC_track_nClusters_TPCVmain);
    fChain->SetBranchAddress("TPC_track_nClusters_TPCVgap",&TPC_track_nClusters_TPCVgap);
    
    fChain->SetBranchAddress("TPC_track_nClustersPotential_Total",&TPC_track_nClustersPotential_Total);
    fChain->SetBranchAddress("TPC_track_nClustersPotential_TPCV1",&TPC_track_nClustersPotential_TPCV1);
    fChain->SetBranchAddress("TPC_track_nClustersPotential_TPCV2",&TPC_track_nClustersPotential_TPCV2);
    fChain->SetBranchAddress("TPC_track_nClustersPotential_TPCVmain",&TPC_track_nClustersPotential_TPCVmain);
    fChain->SetBranchAddress("TPC_track_nClustersPotential_TPCVgap",&TPC_track_nClustersPotential_TPCVgap);
}

//=================================================================> EXEC <===============================================================
void DataTreeConverterNA61::Exec(int event)
{
    fChain -> GetEntry(event);
    Read_Event();
    DTEvent -> Process();
    fTreeQA -> Fill();
    Clear_Event();
}
//--------------------------------------------------------------------------------------------------
void DataTreeConverterNA61::Clear_Event()
{
    DTEvent -> ClearEvent();
}
//--------------------------------------------------------------------------------------------------
void DataTreeConverterNA61::Read_Event()
{
    DTEvent -> SetRunId(Run_Id);
    DTEvent -> SetEventId(Event_Id);
    DTEvent -> SetRPAngle(RPAngle);
    DTEvent -> SetImpactParameter(ImpactParameter);
    DTEvent -> SetEventTimestamp(Event_Timestamp);
    
    DTEvent -> GetTrigger(0) -> SetSignal(triggersADC[0][0]);	//S1
    DTEvent -> GetTrigger(1) -> SetSignal(triggersADC[0][1]);	//S2
    DTEvent -> GetTrigger(2) -> SetSignal(triggersADC[0][4]);	//V1p -> S3
    DTEvent -> GetTrigger(3) -> SetSignal(triggersADC[0][3]);	//V1
    DTEvent -> GetTrigger(4) -> SetSignal(triggersADC[0][5]);	//PSD
    DTEvent -> GetTrigger(5) -> SetSignal(T1);			//T1
    DTEvent -> GetTrigger(6) -> SetSignal(T2);			//T2
    DTEvent -> GetTrigger(7) -> SetSignal(T4);			//T4
    
    for (int i=0;i<nBPD;i++)
    {
	DTEvent -> GetBPD(i) -> SetPosition(BPD_Position[1][i][0],BPD_Position[1][i][1],BPD_Position[1][i][2]);
    }
    
    DTEvent -> SetVertexPosition(FitVertexX,FitVertexY,FitVertexZ);
    DTEvent -> SetVertexQuality(FitVertexQ);
    
    for (int i=0;i<nPSD_Modules;i++)
    {
	DataTreePSDModule* mod = DTEvent -> GetPSDModule(i);
	mod -> SetPosition(PSD_module_X[i],PSD_module_Y[i],PSD_module_Z[i]);
	for (int j=0;j<nPSD_Sections;j++)
	{
	    DataTreePSDSection* sec = mod -> GetSection(j);
	    sec -> AddEnergy(PSD_section_Energy[i][j]);
	}
    }
    
    for (int i=0;i<nTPC_Tracks;i++)
    {
	if (TPC_track_pT[i] <= nUndefinedValue) break;
	DTEvent -> AddTrack();
	DataTreeTrack* track = DTEvent -> GetTrack(i);
	
	track -> SetPt(0,TPC_track_pT[i]);
	track -> SetEta(0,1/2.*TPC_track_eta[i]);//NA61 code specifics (bug fixing)
	track -> SetPhi(0,TPC_track_phi[i]);

	track -> SetDCA(0,TPC_track_DCAtoVertex_X[i],TPC_track_DCAtoVertex_Y[i],nUndefinedValue);
	
	track-> SetNofHits(0,TPC_track_nClusters_Total[i],0);
	track-> SetNofHits(0,TPC_track_nClusters_TPCV1[i],1);
	track-> SetNofHits(0,TPC_track_nClusters_TPCV2[i],2);
	track-> SetNofHits(0,TPC_track_nClusters_TPCVmain[i],3);
	track-> SetNofHits(0,TPC_track_nClusters_TPCVgap[i],4);
	
	track-> SetNofHitsPotential(0,TPC_track_nClustersPotential_Total[i],0);
	track-> SetNofHitsPotential(0,TPC_track_nClustersPotential_TPCV1[i],1);
	track-> SetNofHitsPotential(0,TPC_track_nClustersPotential_TPCV2[i],2);
	track-> SetNofHitsPotential(0,TPC_track_nClustersPotential_TPCVmain[i],3);
	track-> SetNofHitsPotential(0,TPC_track_nClustersPotential_TPCVgap[i],4);
    }
}
//================================================================> FINISH <==============================================================
void DataTreeConverterNA61::Finish()
{
    fTreeQA -> Write();
    fTreeFile -> Write();
    fTreeFile -> Close();
}

ClassImp(DataTreeConverterNA61)