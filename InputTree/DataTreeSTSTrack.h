#ifndef DataTreeSTSTrack_H
#define DataTreeSTSTrack_H 1

#include <vector>
#include <iostream>
#include "TClonesArray.h"
#include "TObject.h"

class DataTreeSTSTrack : public TObject
{
    
public:
  
    DataTreeSTSTrack(int idx);
    ~DataTreeSTSTrack();
    
    int GetId(){return id;}
    double GetPt(){return pT;}
    double GetPhi(){return phi;}
    double GetEta(){return eta;}
    int GetNofHits(){return NofHits;}
    int GetFlag(){return Flag;}
    double GetChiSq(){return ChiSq;}
    int GetNDF(){return NDF;}
    double GetDCAComponent(int idx){return DCA[idx];}
    int GetMCTrackId(){return MCTrack_id;}
    int GetTOFHitId(){return TOFHit_id;}
    
    void SetPt(double fPt){pT = fPt;}
    void SetPhi(double fPhi){phi = fPhi;}
    void SetEta(double fEta){eta = fEta;}
    void SetNofHits(int fNofHits){NofHits = fNofHits;}
    void SetFlag(int fFlag){Flag = fFlag;}
    void SetChiSq(double fChiSq){ChiSq = fChiSq;}
    void SetNDF(int fNDF){NDF = fNDF;}
    void SetDCA(double fX, double fY, double fZ){DCA[0]=fX; DCA[1]=fY; DCA[2]=fZ;}
    void SetDCAComponent(int idx, double fValue){DCA[idx]=fValue;}
    void SetMCTrackId(double idx){MCTrack_id = idx;}
    void SetTOFHitId(double idx){TOFHit_id = idx;}
private:
    
    void SetId(int idx){id = idx;}
  
    int id;
    double pT;
    double phi;
    double eta;
    int NofHits;
    int Flag;// (?definition)
    double ChiSq;
    int NDF;
    double DCA[3]; //(calculate)
    int MCTrack_id;
    int TOFHit_id;

    ClassDefNV(DataTreeSTSTrack, 1)
};

#endif