#ifndef DataTreeEvent_H
#define DataTreeEvent_H 1

#include <vector>
#include <iostream>
#include "TClonesArray.h"
#include "TObject.h"
#include "DataTreeTrack.h"
#include "DataTreePSDModule.h"
#include "DataTreeTOFSegment.h"
#include "DataTreeV0Candidate.h"
#include "DataTreeMCTrack.h"
#include "DataTreeTrigger.h"
#include "DataTreeBPD.h"

const int nV0Types = 3;

class DataTreeEvent : public TObject
{
    
public:
  
    DataTreeEvent();
    ~DataTreeEvent();
    
    void ClearEvent()
    {
	for(int i=0;i<nPSDModules;i++)
	{
	    GetPSDModule(i)->ClearEvent();
	}
	PSDEnergy = 0.;
	nFiredPSDModules = 0;
	nFiredPSDSections = 0;
	
	arrTracks->Clear();
	nTracks = 0;
	ProcessFlag = false;
    }
    void Process()
    {
	PSDEnergy = 0.;
	nFiredPSDModules = 0;
	nFiredPSDSections = 0;
	MCPSDEnergy = 0.;
	for(int i=0;i<nPSDModules;i++)
	{
	    PSDEnergy+=GetPSDModule(i)->GetEnergy();
	    nFiredPSDSections+=GetPSDModule(i)->GetNFiredSections();
	    if(GetPSDModule(i)->GetEnergy() > 0)nFiredPSDModules++;
	}
	for (int i=0;i<nMCTracks;i++)
	{
	    MCPSDEnergy+=GetMCTrack(i)->GetPSDEnergy();
	}
	ProcessFlag = true;
    }
    
//     bool GetProcessFlag(){return ProcessFlag;}

    int GetRunId(){return RunId;}
    int GetEventId(){return EventId;}
    double GetEventTimestamp(){return EventTimestamp;}
    double GetRPAngle(){return RPAngle;}
    double GetImpactParameter(){return ImpactParameter;}
    double GetVertexPositionComponent(int idx){return VertexPosition[idx];}
    double GetVertexQuality(){return VertexQuality;}
    int GetNFiredPSDModules(){Process(); return nFiredPSDModules;}
    int GetNFiredPSDSections(){Process(); return nFiredPSDSections;}
    double GetPSDEnergy(){Process(); return PSDEnergy;}
    double GetMCPSDEnergy(){Process(); return MCPSDEnergy;}
    double GetMCVertexPositionComponent(int idx){return MCVertexPosition[idx];}
    
    void SetRunId(int fValue){RunId = fValue;}
    void SetEventId(int fValue){EventId = fValue;}
    void SetEventTimestamp(double fValue){EventTimestamp = fValue;}
    void SetRPAngle(double fValue){RPAngle = fValue;}
    void SetImpactParameter(double fValue){ImpactParameter = fValue;}
    void SetVertexPosition(double fX, double fY, double fZ){VertexPosition[0]=fX; VertexPosition[1]=fY; VertexPosition[2]=fZ;}
    void SetVertexPositionComponent(int idx, double fValue){VertexPosition[idx] = fValue;}
    void SetVertexQuality(double fValue){VertexQuality = fValue;}
    void SetMCVertexPosition(double fX, double fY, double fZ){MCVertexPosition[0]=fX; MCVertexPosition[1]=fY; MCVertexPosition[2]=fZ;}
    void SetMCVertexPositionComponent(int idx, double fValue){MCVertexPosition[idx] = fValue;}
        
    int GetNTracks(){return nTracks;}
    DataTreeTrack* GetTrack(int idx){return (DataTreeTrack*)arrTracks->At(idx);}
    DataTreeTrack* GetLastTrack(){return (DataTreeTrack*)arrTracks->At(nTracks-1);}
    void AddTrack(){TClonesArray &arr = *arrTracks; new(arr[nTracks]) DataTreeTrack(nTracks); nTracks++;}

    int GetNPSDModules(){return nPSDModules;}
    DataTreePSDModule* GetPSDModule(int idx){return (DataTreePSDModule*)arrPSDModules->At(idx);}
    DataTreePSDModule* GetLastPSDModule(){return (DataTreePSDModule*)arrPSDModules->At(nPSDModules-1);}
    void AddPSDModule(){TClonesArray &arr = *arrPSDModules; new(arr[nPSDModules]) DataTreePSDModule(nPSDModules); nPSDModules++;}
    void AddPSDModule(int inSections){TClonesArray &arr = *arrPSDModules; new(arr[nPSDModules]) DataTreePSDModule(nPSDModules,inSections); nPSDModules++;}
    
    int GetNTOFSegments(){return nTOFSegments;}
    DataTreeTOFSegment* GetTOFSegment(int idx){return (DataTreeTOFSegment*)arrTOFSegments->At(idx);}
    DataTreeTOFSegment* GetLastTOFSegment(){return (DataTreeTOFSegment*)arrTOFSegments->At(nTOFSegments-1);}
    void AddTOFSegment(){TClonesArray &arr = *arrTOFSegments; new(arr[nTOFSegments]) DataTreeTOFSegment(nTOFSegments); nTOFSegments++;}
    
    int GetNV0Candidates(){return nV0Candidates;}
    DataTreeV0Candidate* GetV0Candidate(int idx){return (DataTreeV0Candidate*)arrV0Candidates->At(idx);}
    DataTreeV0Candidate* GetLastV0Candidate(){return (DataTreeV0Candidate*)arrV0Candidates->At(nV0Candidates-1);}
    void AddV0Candidate(){TClonesArray &arr = *arrV0Candidates; new(arr[nV0Candidates]) DataTreeV0Candidate(nV0Candidates); nV0Candidates++;}

    int GetNV0CandidatesMCpid(){return nV0CandidatesMCpid;}
    DataTreeV0Candidate* GetV0CandidateMCpid(int idx){return (DataTreeV0Candidate*)arrV0CandidatesMCpid->At(idx);}
    DataTreeV0Candidate* GetLastV0CandidateMCpid(){return (DataTreeV0Candidate*)arrV0CandidatesMCpid->At(nV0CandidatesMCpid-1);}
    void AddV0CandidateMCpid(){TClonesArray &arr = *arrV0CandidatesMCpid; new(arr[nV0CandidatesMCpid]) DataTreeV0Candidate(nV0CandidatesMCpid); nV0CandidatesMCpid++;}
    
    int GetNMCTracks(){return nMCTracks;}
    DataTreeMCTrack* GetMCTrack(int idx){return (DataTreeMCTrack*)arrMCTracks->At(idx);}
    DataTreeMCTrack* GetLastMCTrack(){return (DataTreeMCTrack*)arrMCTracks->At(nMCTracks-1);}
    void AddMCTrack(){TClonesArray &arr = *arrMCTracks; new(arr[nMCTracks]) DataTreeMCTrack(nMCTracks); nMCTracks++;}
    
    int GetNTriggers(){return nTriggers;}
    DataTreeTrigger* GetTrigger(int idx){return (DataTreeTrigger*)arrTriggers->At(idx);}
    DataTreeTrigger* GetLastTrigger(){return (DataTreeTrigger*)arrTriggers->At(nTriggers-1);}
    void AddTrigger(){TClonesArray &arr = *arrTriggers; new(arr[nTriggers]) DataTreeTrigger(nTriggers); nTriggers++;} 
    void AddTrigger(TString label){TClonesArray &arr = *arrTriggers; new(arr[nTriggers]) DataTreeTrigger(nTriggers, label); nTriggers++;}     
    
    int GetNBPDs(){return nBPDs;}
    DataTreeBPD* GetBPD(int idx){return (DataTreeBPD*)arrBPDs->At(idx);}
    DataTreeBPD* GetLastBPD(){return (DataTreeBPD*)arrBPDs->At(nBPDs-1);}
    void AddBPD(){TClonesArray &arr = *arrBPDs; new(arr[nBPDs]) DataTreeBPD(nBPDs); nBPDs++;}
    
private:

    bool ProcessFlag;

    int RunId;
    int EventId;
    double EventTimestamp;
    double VertexPosition[3];		//Position of the vertex
    double VertexQuality;		//Quality of vertex fit
    
    double MCVertexPosition[3];		//Position of the vertex in MC
    double RPAngle;			//Reaction plane angle
    double ImpactParameter;		//Impact parameter
    double MCPSDEnergy;			//sum of PSDEnergies for MC tracks
    
    int nTracks;			//Number of tracks
    TClonesArray* arrTracks;		//tracks
    
    int nFiredPSDModules;		//number of modules in PSD having some energy deposit
    int nFiredPSDSections;		//number of sections in PSD having some energy deposit
    double PSDEnergy;			//total PSD energy
    int nPSDModules;			//number of PSD modules
    TClonesArray* arrPSDModules;	//PSD modules
    
    int nTOFSegments;			//number of segments in TOF
    TClonesArray* arrTOFSegments;	//TOF segments
    
    int nV0Candidates;				//number of V0 candidates
    int nV0SpecificCandidates[nV0Types];	//number of V0 candidates of specific type
    TClonesArray* arrV0Candidates;		//V0 candidates
    
    int nV0CandidatesMCpid;			//number of V0 candidates with MC PiD
    int nV0SpecificCandidatesMCpid[nV0Types];	//number of V0 candidates wit MC PiD and of specific type
    TClonesArray* arrV0CandidatesMCpid;		//V0 candidates with MC PiD
    
    int nMCTracks;			//Multiplicity of MC tracks
    TClonesArray* arrMCTracks;		//MC tracks
    
    //BEGIN NEW

    int nTriggers;			//Number of nTriggers
    TClonesArray* arrTriggers;		//Triggers
    
    int nBPDs;				//Number of BPDs
    TClonesArray* arrBPDs;		//BPDs
    
    //END NEW
  
    ClassDefNV(DataTreeEvent, 2)
};

#endif