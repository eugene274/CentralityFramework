#ifndef DataTreePSDModule_H
#define DataTreePSDModule_H 1

#include <vector>
#include <iostream>
#include "TClonesArray.h"
#include "TObject.h"
#include "DataTreePSDSection.h"
#include "TMath.h"

class DataTreePSDModule : public TObject
{
    
public:
  
    DataTreePSDModule(int idx = 0);
    DataTreePSDModule(int idx, int inSections);
    ~DataTreePSDModule();
    
    void ClearEvent()
    {
        for(int i=0;i<nSections;i++)
        {
            GetSection(i)->ClearEvent();
        }
        Energy = 0.;
        nFiredSections = 0;
        ProcessFlag = false;
    }
    void Process()
    {
        Energy = 0.;
        nFiredSections = 0;
        for(int i=0;i<SectionNumberCut;i++)
        {
            Energy+=GetSection(i)->GetEnergy();
            if(GetSection(i)->GetEnergy() > 0){nFiredSections++;}
        }
        ProcessFlag = true;
    }
    
//     bool GetProcessFlag(){return ProcessFlag;}
    int GetId(){return id;}
    int GetSectionNumberCut(){return SectionNumberCut;}
    double GetPositionComponent(int idx){return Position[idx];}
//     double GetPhi(){if (TMath::ATan2(Position[1],Position[0]) < 0){return 2*TMath::Pi()+TMath::ATan2(Position[1],Position[0]);}else {return TMath::ATan2(Position[1],Position[0]);}}
    double GetPhi(){return TMath::ATan2(Position[1],Position[0]);}
    double GetEnergy(){ Process(); return Energy; }
    int GetNFiredSections(){ Process(); return nFiredSections; }
    
    void SetSectionNumberCut(int fValue){SectionNumberCut = fValue;}
    void SetPosition(double fX, double fY, double fZ){Position[0]=fX; Position[1]=fY; Position[2]=fZ;}
    void SetPositionComponent(int idx, double fValue){Position[idx] = fValue;}
//     void SetEnergy(double fValue){Energy = fValue;}
//     void SetNFiredSections(int fValue){nFiredSections=fValue;}
    
    int GetNSections(){return nSections;}
    DataTreePSDSection* GetSection(int idx){return (DataTreePSDSection*)arrSections->At(idx);}
    DataTreePSDSection* GetLastSection(){return (DataTreePSDSection*)arrSections->At(nSections-1);}
    void AddSection(){TClonesArray &arr = *arrSections; new(arr[nSections]) DataTreePSDSection(nSections); nSections++;}
    
private:    
    void SetId(int idx){id = idx;}
  
    bool ProcessFlag;
    
    int nSections;			//number of section in the module
    int SectionNumberCut;               //max depth in sections to calculate energy
    
    int id;
    double Position[3];			//Position of the module in lab frame
    double Energy;			//energy deposit in the module
    int nFiredSections;			//number of sections where the energy was deposited
    TClonesArray* arrSections;		//sections in the module
    
    ClassDefNV(DataTreePSDModule, 1)
};

#endif