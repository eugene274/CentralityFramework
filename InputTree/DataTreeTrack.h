#ifndef DataTreeTrack_H
#define DataTreeTrack_H 1

#include <vector>
#include <iostream>
#include <fstream>
#include "TClonesArray.h"
#include "TObject.h"
#include "TMath.h"

// #include "DataTreeConstants.h"

class DataTreeTrack : public TObject
{
    
public:
    static const int nZPositions = 4;
    static const int nMaxStations = 8;
    static const int nSubDetectors = 5;
  
    DataTreeTrack(int idx = 0);
    ~DataTreeTrack();
    
    int GetId(){return id;}

//     void SetUndefinedValues(double fValue);
    void SetUndefinedValues();
    
    static int GetnZPositions(){return nZPositions;}
    static int GetnMaxStations(){return nMaxStations;}
    static int GetnSubDetectors(){return nSubDetectors;}
    
    double GetPt(int idx){return pT[idx];}
//     double GetPhi(int idx){if (phi[idx]<0) {return (2*TMath::Pi()+phi[idx]);} else {return phi[idx];}}
    double GetPhi(int idx){return phi[idx];}
    double GetEta(int idx){return eta[idx];}
    double GetPx(int idx){return px[idx];}
    double GetPy(int idx){return py[idx];}
    double GetPz(int idx){return pz[idx];}
    double GetP(int idx){return p[idx];}
    
    int GetNofHits(int idx, int subidx = 0){return NofHits[idx][subidx];}
    int GetNofHitsPotential(int idx, int subidx = 0){return NofHitsPotential[idx][subidx];}
    int GetFlag(int idx){return Flag[idx];}
    double GetChiSq(int idx){return ChiSq[idx];}
    double GetVtxChiSq(int idx){return VtxChiSq[idx];}
    int GetNDF(int idx){return NDF[idx];}
    double GetDCAComponent(int idx, int jdx){return DCA[idx][jdx];}
    bool GetStation(int idx, int jdx){return Stations[idx][jdx];}
    int GetNStations(int idx){return nStations[idx];}
    int GetSTSHitsPossible(int idx){return nSTSHitsPossible[idx];}
    double GetLengthInSTS(int idx){return LengthInSTS[idx];}
    double GetCharge(int idx){return Charge[idx];}
    
    int GetPSDModuleId(){return PSDModuleId;}
    int GetTOFSegmentId(){return TOFSegmentId;}
    int GetMCTrackId(){return MCTrackId;}    


    void SetPt(int idx, double fPt){if(fPt < 0) std::cout << "pT < 0!!! " << fPt << std::endl; pT[idx] = fPt;}
    void SetPhi(int idx, double fPhi){phi[idx] = fPhi;}
    void SetEta(int idx, double fEta){eta[idx] = fEta;}
    void SetPx(int idx, double fValue){px[idx] = fValue;}
    void SetPy(int idx, double fValue){py[idx] = fValue;}
    void SetPz(int idx, double fValue){pz[idx] = fValue;}
    void SetP(int idx, double fValue){p[idx] = fValue;}
    
    void SetNofHits(int idx, int fNofHits, int subidx = 0){NofHits[idx][subidx] = fNofHits;}
    void SetNofHitsPotential(int idx, int fNofHitsPotential, int subidx = 0){NofHitsPotential[idx][subidx] = fNofHitsPotential;}
    void SetFlag(int idx, int fFlag){Flag[idx] = fFlag;}
    void SetChiSq(int idx, double fChiSq){ChiSq[idx] = fChiSq;}
    void SetVtxChiSq(int idx, double fValue){VtxChiSq[idx] = fValue;}
    void SetNDF(int idx, int fNDF){NDF[idx] = fNDF;}
    void SetDCA(int idx, double fX, double fY, double fZ){DCA[idx][0]=fX; DCA[idx][1]=fY; DCA[idx][2]=fZ;}
    void SetDCAComponent(int idx, int jdx, double fValue){DCA[idx][jdx] = fValue;}
    void SetStation(int idx, int jdx, bool fValue){Stations[idx][jdx] = fValue;}
    void SetNStations(int idx, int fValue){nStations[idx] = fValue;}
    void SetSTSHitsPossible(int idx, int fValue){nSTSHitsPossible[idx] = fValue;}
    void SetLengthInSTS(int idx, double fValue){LengthInSTS[idx] = fValue;}
    void SetCharge(int idx, double fValue){Charge[idx] = fValue;}
    
    void SetPSDModuleId(double idx){PSDModuleId = idx;}
    void SetTOFSegmentId(double idx){TOFSegmentId = idx;}
    void SetMCTrackId(double idx){MCTrackId = idx;}
    
private:
//     const int nUndefinedValue = -999;
    
    void SetId(int idx){id = idx;}
  
    int id;
    
    //kinematics at different z-positions of track: 0--first STS hit, 1--vertex, 2--TOF, 3--PSD
    double pT[nZPositions];
    double phi[nZPositions];
    double eta[nZPositions];
    double px[nZPositions];
    double py[nZPositions];
    double pz[nZPositions];
    double p[nZPositions];
    
    int NofHits[nZPositions][nSubDetectors];
    int NofHitsPotential[nZPositions][nSubDetectors];
    int Flag[nZPositions];
    double ChiSq[nZPositions];
    double VtxChiSq[nZPositions];
    int NDF[nZPositions];
    double DCA[nZPositions][3];
    bool Stations[nZPositions][nMaxStations];
    int nStations[nZPositions];
    int nSTSHitsPossible[nZPositions];
    double LengthInSTS[nZPositions];
    double Charge[nZPositions];
       
    int PSDModuleId;		//id of the corresponding PSD module
    int TOFSegmentId;		//id of the corresponding TOF segment
    int MCTrackId;		//id of the best matched MC track

    

    
    
    ClassDefNV(DataTreeTrack, 1)
};

#endif