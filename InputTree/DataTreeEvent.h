#ifndef DataTreeEvent_H
#define DataTreeEvent_H 1

#include <vector>
#include <iostream>
#include "TClonesArray.h"
#include "TObject.h"

#include "DataTreeTrack.h"
#include "DataTreePSDModule.h"
#include "DataTreeTOFHit.h"
#include "DataTreeV0Candidate.h"
#include "DataTreeMCTrack.h"
#include "DataTreeTrigger.h"
#include "DataTreeBPD.h"

const int nV0Types = 3;
const int V0Pdg[nV0Types] = {3122,-3122,-310};

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
	arrMCTracks->Clear();
	nMCTracks = 0;
	arrV0CandidatesTOFpid->Clear();
	nV0CandidatesTOFpid = 0;
	arrV0CandidatesMCpid->Clear();
	nV0CandidatesMCpid = 0;
	for (int i=0;i<nV0Types;i++)
	{
	    nV0SpecificCandidatesTOFpid[i] = 0;
	    nV0SpecificCandidatesMCpid[i] = 0;
	}
	arrTOFHits->Clear();
	nTOFHits = 0;
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
    const int GetNV0Types(){return nV0Types;}
    const int GetNV0Pdg(int idx){return V0Pdg[idx];}
    
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
        
    void SetNV0SpecificCandidatesTOFpid(int idx, int fValue){ nV0SpecificCandidatesTOFpid[idx] = fValue; }
    void SetNV0CandidatesTOFpid(int fValue){ nV0CandidatesTOFpid = fValue; }
    void SetNV0SpecificCandidatesMCpid(int idx, int fValue){ nV0SpecificCandidatesMCpid[idx] = fValue; }
    void SetNV0CandidatesMCpid(int fValue){ nV0CandidatesMCpid = fValue; }
    
    int GetNTracks(){return nTracks;}
    DataTreeTrack* GetTrack(int idx){return (DataTreeTrack*)arrTracks->At(idx);}
    DataTreeTrack* GetLastTrack(){return (DataTreeTrack*)arrTracks->At(nTracks-1);}
    void AddTrack(){TClonesArray &arr = *arrTracks; new(arr[nTracks]) DataTreeTrack(nTracks); nTracks++; ClearTrack(nTracks-1);}
    void ClearTrack(int idx){ GetTrack(idx)->SetUndefinedValues(); }
    
    int GetNPSDModules(){return nPSDModules;}
    DataTreePSDModule* GetPSDModule(int idx){return (DataTreePSDModule*)arrPSDModules->At(idx);}
    DataTreePSDModule* GetLastPSDModule(){return (DataTreePSDModule*)arrPSDModules->At(nPSDModules-1);}
    void AddPSDModule(){TClonesArray &arr = *arrPSDModules; new(arr[nPSDModules]) DataTreePSDModule(nPSDModules); nPSDModules++;}
    void AddPSDModule(int inSections){TClonesArray &arr = *arrPSDModules; new(arr[nPSDModules]) DataTreePSDModule(nPSDModules,inSections); nPSDModules++;}
    
    int GetNTOFHits(){return nTOFHits;}
    DataTreeTOFHit* GetTOFHit(int idx){return (DataTreeTOFHit*)arrTOFHits->At(idx);}
    DataTreeTOFHit* GetLastTOFHit(){return (DataTreeTOFHit*)arrTOFHits->At(nTOFHits-1);}
//     void AddTOFHit(){TClonesArray &arr = *arrTOFHits; new(arr[nTOFHits]) DataTreeTOFHit(nTOFHits); nTOFHits++;}
    DataTreeTOFHit* AddTOFHit()
    {
	TClonesArray &arr = *arrTOFHits;
	new(arr[nTOFHits]) DataTreeTOFHit(nTOFHits);
	nTOFHits++;
	return (DataTreeTOFHit*)arrTOFHits->At(nTOFHits-1);
    }
    
    int GetNV0SpecificCandidatesTOFpid(int idx){ return nV0SpecificCandidatesTOFpid[idx]; }
    int GetNV0CandidatesTOFpid(){return nV0CandidatesTOFpid;}
    DataTreeV0Candidate* GetV0CandidateTOFpid(int idx){return (DataTreeV0Candidate*)arrV0CandidatesTOFpid->At(idx);}
    DataTreeV0Candidate* GetLastV0CandidateTOFpid(){return (DataTreeV0Candidate*)arrV0CandidatesTOFpid->At(nV0CandidatesTOFpid-1);}
    DataTreeV0Candidate* AddV0CandidateTOFpid()
    {
	TClonesArray &arr = *arrV0CandidatesTOFpid;
	new(arr[nV0CandidatesTOFpid]) DataTreeV0Candidate(nV0CandidatesTOFpid);
	nV0CandidatesTOFpid++;
	return (DataTreeV0Candidate*)arrV0CandidatesTOFpid->At(nV0CandidatesTOFpid-1);
    }

    int GetNV0SpecificCandidatesMCpid(int idx){ return nV0SpecificCandidatesMCpid[idx]; }
    int GetNV0CandidatesMCpid(){return nV0CandidatesMCpid;}
    DataTreeV0Candidate* GetV0CandidateMCpid(int idx){return (DataTreeV0Candidate*)arrV0CandidatesMCpid->At(idx);}
    DataTreeV0Candidate* GetLastV0CandidateMCpid(){return (DataTreeV0Candidate*)arrV0CandidatesMCpid->At(nV0CandidatesMCpid-1);}
    DataTreeV0Candidate* AddV0CandidateMCpid()
    {
	TClonesArray &arr = *arrV0CandidatesMCpid;
	new(arr[nV0CandidatesMCpid]) DataTreeV0Candidate(nV0CandidatesMCpid);
	nV0CandidatesMCpid++;
	return (DataTreeV0Candidate*)arrV0CandidatesMCpid->At(nV0CandidatesMCpid-1);
    }
    
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
    
    int nTOFHits;			//number of segments in TOF
    TClonesArray* arrTOFHits;	//TOF segments
    
    int nV0CandidatesTOFpid;				//number of V0 candidates
    int nV0SpecificCandidatesTOFpid[nV0Types];	//number of V0 candidates of specific type
    TClonesArray* arrV0CandidatesTOFpid;		//V0 candidates
    
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