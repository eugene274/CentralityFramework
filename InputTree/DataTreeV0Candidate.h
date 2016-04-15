#ifndef DataTreeV0Candidate_H
#define DataTreeV0Candidate_H 1

#include <vector>
#include <iostream>
#include "TClonesArray.h"
#include "TObject.h"

class DataTreeV0Candidate : public TObject
{
    
public:
  
    DataTreeV0Candidate(int idx = 0);
    ~DataTreeV0Candidate();
    
    int GetId(){return id;}
    
    double GetPt(){return pT;}
    double GetPhi(){return phi;}
    double GetEta(){return eta;}
    int GetPdgId(){return PdgId;}
    int GetMass(){return Mass;}
    double GetChiSq(){return ChiSq;}

    int GetTrackId(){return TrackId;}

    void SetPt(double fPt){pT = fPt;}
    void SetPhi(double fPhi){phi = fPhi;}
    void SetEta(double fEta){eta = fEta;}
    void SetPdgId(double fValue){PdgId = fValue;}
    void SetMass(double fValue){Mass = fValue;}
    void SetChiSq(double fChiSq){ChiSq = fChiSq;}
    
    void SetTrackId(double idx){TrackId = idx;}
    
    int GetNDaughters(){return nDaughters;}
    DataTreeV0Candidate* GetDaughter(int idx){return (DataTreeV0Candidate*)arrDaughters->At(idx);}
    DataTreeV0Candidate* GetLastDaughter(){return (DataTreeV0Candidate*)arrDaughters->At(nDaughters-1);}
    void AddDaughter(){TClonesArray &arr = *arrDaughters; new(arr[nDaughters]) DataTreeV0Candidate(nDaughters); nDaughters++;}
    
private:
    
    void SetId(int idx){id = idx;}
    
    int id;
    double pT;			//pT
    double phi;			//phi
    double eta;			//eta
    double PdgId;		//pdg
    double Mass;		//mass
    double ChiSq;		//Chi squared
    
    int TrackId;		//id of the corresponding track

    int nDaughters;		//Number of daughter tracks
    TClonesArray* arrDaughters;	//Daughter tracks
    
    ClassDefNV(DataTreeV0Candidate, 1)
};

#endif