#ifndef DataTreeTrack_H
#define DataTreeTrack_H 1

#include <vector>
#include <iostream>
#include <fstream>
#include "TClonesArray.h"
#include "TObject.h"

// #include "DataTreeConstants.h"

class DataTreeTrack : public TObject
{
    
public:
  
    DataTreeTrack(int idx = 0);
    ~DataTreeTrack();
    
    int GetId(){return id;}

    double GetPt(int idx){return pT[idx];}
    double GetPhi(int idx){return phi[idx];}
    double GetEta(int idx){return eta[idx];}
    
    int GetNofHits(int idx){return NofHits[idx];}
    int GetFlag(int idx){return Flag[idx];}
    double GetChiSq(int idx){return ChiSq[idx];}
    int GetNDF(int idx){return NDF[idx];}
    double GetDCAComponent(int idx, int jdx){return DCA[idx][jdx];}
    bool GetStation(int idx, int jdx){return Stations[idx][jdx];}
    int GetNStations(int idx){return nStations[idx];}
    int GetSTSHitsPossible(int idx){return nSTSHitsPossible[idx];}
    double GetLengthInSTS(int idx){return LengthInSTS[idx];}
    
    int GetPSDModuleId(){return PSDModuleId;}
    int GetTOFSegmentId(){return TOFSegmentId;}
    int GetMCTrackId(){return MCTrackId;}
    
//     double GetVtxPt(){return Vtx_pT;}
//     double GetVtxPhi(){return Vtx_phi;}
//     double GetVtxEta(){return Vtx_eta;}
//     double GetPt(){return pT;}
//     double GetPhi(){return phi;}
//     double GetEta(){return eta;}
//     double GetTOFPt(){return TOF_pT;}
//     double GetTOFPhi(){return TOF_phi;}
//     double GetTOFEta(){return TOF_eta;}
//     double GetPSDPt(){return PSD_pT;}
//     double GetPSDPhi(){return PSD_phi;}
//     double GetPSDEta(){return PSD_eta;}
//     
//     int GetNofHits(){return NofHits;}
//     int GetFlag(){return Flag;}
//     double GetChiSq(){return ChiSq;}
//     int GetNDF(){return NDF;}
//     double GetDCAComponent(int idx){return DCA[idx];}
//     bool GetStation(int idx){return Stations[idx];}
//     int GetNStations(){return nStations;}
//     int GetSTSHitsPossible(){return nSTSHitsPossible;}
//     double GetLengthInSTS(){return LengthInSTS;}
    


    void SetPt(int idx, double fPt){if(fPt < 0) std::cout << "pT < 0!!! " << fPt << std::endl; pT[idx] = fPt;}
    void SetPhi(int idx, double fPhi){phi[idx] = fPhi;}
    void SetEta(int idx, double fEta){eta[idx] = fEta;}
    
    void SetNofHits(int idx, int fNofHits){NofHits[idx] = fNofHits;}
    void SetFlag(int idx, int fFlag){Flag[idx] = fFlag;}
    void SetChiSq(int idx, double fChiSq){ChiSq[idx] = fChiSq;}
    void SetNDF(int idx, int fNDF){NDF[idx] = fNDF;}
    void SetDCA(int idx, double fX, double fY, double fZ){DCA[idx][0]=fX; DCA[idx][1]=fY; DCA[idx][2]=fZ;}
    void SetDCAComponent(int idx, int jdx, double fValue){DCA[idx][jdx] = fValue;}
    void SetStation(int idx, int jdx, bool fValue){Stations[idx][jdx] = fValue;}
    void SetNStations(int idx, int fValue){nStations[idx] = fValue;}
    void SetSTSHitsPossible(int idx, int fValue){nSTSHitsPossible[idx] = fValue;}
    void SetLengthInSTS(int idx, double fValue){LengthInSTS[idx] = fValue;}
    
    void SetPSDModuleId(double idx){PSDModuleId = idx;}
    void SetTOFSegmentId(double idx){TOFSegmentId = idx;}
    void SetMCTrackId(double idx){MCTrackId = idx;}
    
//     void SetVtxPt(double fVtxPt){Vtx_pT = fVtxPt;}
//     void SetVtxPhi(double fVtxPhi){Vtx_phi = fVtxPhi;}
//     void SetVtxEta(double fVtxEta){Vtx_eta = fVtxEta;}
//     void SetPt(double fPt){pT = fPt;}
//     void SetPhi(double fPhi){phi = fPhi;}
//     void SetEta(double fEta){eta = fEta;}
//     void SetTOFPt(double fTOFPt){TOF_pT = fTOFPt;}
//     void SetTOFPhi(double fTOFPhi){TOF_phi = fTOFPhi;}
//     void SetTOFEta(double fTOFEta){TOF_eta = fTOFEta;}
//     void SetPSDPt(double fPSDPt){PSD_pT = fPSDPt;}
//     void SetPSDPhi(double fPSDPhi){PSD_phi = fPSDPhi;}
//     void SetPSDEta(double fPSDEta){PSD_eta = fPSDEta;}
//     
//     void SetNofHits(int fNofHits){NofHits = fNofHits;}
//     void SetFlag(int fFlag){Flag = fFlag;}
//     void SetChiSq(double fChiSq){ChiSq = fChiSq;}
//     void SetNDF(int fNDF){NDF = fNDF;}
//     void SetDCA(double fX, double fY, double fZ){DCA[0]=fX; DCA[1]=fY; DCA[2]=fZ;}
//     void SetDCAComponent(int idx, double fValue){DCA[idx] = fValue;}
//     void SetStation(int idx, bool fValue){Stations[idx] = fValue;}
//     void SetNStations(int fValue){nStations = fValue;}
//     void SetSTSHitsPossible(int fValue){nSTSHitsPossible = fValue;}
//     void SetLengthInSTS(double fValue){LengthInSTS = fValue;}
    
private:
    
    void SetId(int idx){id = idx;}
  
    int id;
    
    double pT[4];
    double phi[4];
    double eta[4];
    int NofHits[4];
    int Flag[4];
    double ChiSq[4];
    int NDF[4];
    double DCA[4][3];
    bool Stations[4][8];
    int nStations[4];
    int nSTSHitsPossible[4];
    double LengthInSTS[4];
    
    int PSDModuleId;		//id of the corresponding PSD module
    int TOFSegmentId;		//id of the corresponding TOF segment
    int MCTrackId;		//id of the best matched MC track

//     double Vtx_pT;		//pT for track extrapolated to vertex
//     double Vtx_phi;		//phi for track extrapolated to vertex
//     double Vtx_eta;		//eta for track extrapolated to vertex
//     double pT;			//pT for track in the first hit point
//     double phi;			//phi for track in the first hit point
//     double eta;			//eta for track in the first hit point
//     double TOF_pT;		//pT for track extrapolated to TOF
//     double TOF_phi;		//phi for track extrapolated to TOF
//     double TOF_eta;		//eta for track extrapolated to TOF
//     double PSD_pT;		//pT for track extrapolated to PSD
//     double PSD_phi;		//phi for track extrapolated to PSD
//     double PSD_eta;		//eta for track extrapolated to PSD
// 
//     int NofHits;			//number of hits
//     int Flag;// (?definition)		//flag
//     double ChiSq;			//Chi Squared
//     int NDF;				//Number of Degrees of Freedom
//     double DCA[3]; //(calculate)	//Distance of Closest Approach
//     bool Stations[8];	//Fired stations in STS
//     int nStations;			//Number of fired stations in STS
//     int nSTSHitsPossible;		//Number of possible hits in STS
//     double LengthInSTS;			//Length of the track in STS
    

    
    
    ClassDefNV(DataTreeTrack, 1)
};

#endif