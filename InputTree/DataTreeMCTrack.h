#ifndef DataTreeMCTrack_H
#define DataTreeMCTrack_H 1

#include <vector>
#include <iostream>
#include "TClonesArray.h"
#include "TObject.h"

class DataTreeMCTrack : public TObject
{
    
public:
  
    DataTreeMCTrack(int idx = 0);
    ~DataTreeMCTrack();
    
    int GetId(){return id;}
    
    bool GetIfPrimary(){return ifPrimary;}
    
    double GetPx(){return px;}
    double GetPy(){return py;}
    double GetPz(){return pz;}
    double GetP(){return p;}
    double GetPt(){return pT;}
    double GetPhi(){return phi;}
    double GetEta(){return eta;}
    double GetPdgId(){return PdgId;}
    double GetMass(){return Mass;}
    double GetCharge(){return Charge;}
    int GetMotherId(){return MotherId;}
    
    double GetTOFPositionComponent(int idx){return TOFPosition[idx];}
    double GetTOFPt(){return TOF_pT;}
    double GetTOFPhi(){return TOF_phi;}
    double GetTOFEta(){return TOF_eta;}
    
    double GetPSDPositionComponent(int idx){return PSDPosition[idx];}
    double GetPSDPt(){return PSD_pT;}
    double GetPSDPhi(){return PSD_phi;}
    double GetPSDEta(){return PSD_eta;}
    
    double GetPSDEnergy(){return PSDEnergy;}
    
    int GetPSDSectionId(){return PSDSectionId;}
    int GetTOFSegmentId(){return TOFSegmentId;}
    
    void SetIfPrimary(bool fValue){ifPrimary = fValue;}

    void SetPt(double fPt){pT = fPt;}
    void SetPhi(double fPhi){phi = fPhi;}
    void SetEta(double fEta){eta = fEta;}
    void SetPx(double fValue){px = fValue;}
    void SetPy(double fValue){py = fValue;}
    void SetPz(double fValue){pz = fValue;}
    void SetP(double fValue){p = fValue;}
    void SetPdgId(double fValue){PdgId = fValue;}
    void SetCharge(double fValue){Charge = fValue;}
    void SetMass(double fValue){Mass = fValue;}
    void SetMotherId(int fValue){MotherId = fValue;}

    void SetTOFPosition(double fX, double fY, double fZ){TOFPosition[0]=fX; TOFPosition[1]=fY; TOFPosition[2]=fZ;}
    void SetTOFPositionComponent(int idx, double fValue){TOFPosition[idx] = fValue;}
    void SetTOFPt(double fTOFPt){TOF_pT = fTOFPt;}
    void SetTOFPhi(double fTOFPhi){TOF_phi = fTOFPhi;}
    void SetTOFEta(double fTOFEta){TOF_eta = fTOFEta;}
    
    void SetPSDPosition(double fX, double fY, double fZ){PSDPosition[0]=fX; PSDPosition[1]=fY; PSDPosition[2]=fZ;}
    void SetPSDPositionComponent(int idx, double fValue){PSDPosition[idx] = fValue;}
    void SetPSDPt(double fPSDPt){PSD_pT = fPSDPt;}
    void SetPSDPhi(double fPSDPhi){PSD_phi = fPSDPhi;}
    void SetPSDEta(double fPSDEta){PSD_eta = fPSDEta;}

    void SetPSDEnergy(double fValue){PSDEnergy = fValue;}
    
    void SetPSDSectionId(double idx){PSDSectionId = idx;}
    void SetTOFSegmentId(double idx){TOFSegmentId = idx;}

private:
    
    void SetId(int idx){id = idx;}
  
    int id;

    bool ifPrimary;		//flag if it is a primary particle  
    
    double px;
    double py;
    double pz;
    double p;
    
    double pT;			//pT for track in the first hit point
    double phi;			//phi for track in the first hit point
    double eta;			//eta for track in the first hit point
    double PdgId;		//pdg code
    double Mass;		//mass
    double Charge;		//charge
    int MotherId;		//mother id (-1 == primary)
    
    double TOFPosition[3];	//Position at TOF
    double TOF_pT;		//pT at TOF
    double TOF_phi;		//phi at TOF
    double TOF_eta;		//eta at TOF

    double PSDPosition[3];	//Position at PSD
    double PSD_pT;		//pT for track extrapolated to PSD
    double PSD_phi;		//phi for track extrapolated to PSD
    double PSD_eta;		//eta for track extrapolated to PSD

    double PSDEnergy;		//sum over PsdPoints for a given MC track
    
    int PSDSectionId;		//id of the corresponding PSD section
    int TOFSegmentId;		//id of the corresponding TOF segment
    
    ClassDefNV(DataTreeMCTrack, 1)
};

#endif