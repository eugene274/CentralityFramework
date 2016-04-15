#ifndef DataTreeTOFSegment_H
#define DataTreeTOFSegment_H 1

#include <vector>
#include <iostream>
#include "TClonesArray.h"
#include "TObject.h"

class DataTreeTOFSegment : public TObject
{
    
public:
  
    DataTreeTOFSegment(int idx = 0);
    ~DataTreeTOFSegment();
    
    void ClearEvent(){Time = 0.; MassSq = 0.;}
    
    int GetId(){return id;}
    
    double GetPositionComponent(int idx){return Position[idx];}
    double GetTime(){return Time;}
    double GetMassSq(){return MassSq;}
    
    void SetPositionComponent(int idx, double fValue){Position[idx]=fValue;}
    void SetPosition(double fX, double fY, double fZ){Position[0]=fX; Position[1]=fY; Position[2]=fZ;}
    void SetTime(double fValue){Time = fValue;}
    void SetMassSq(double fValue){MassSq = fValue;}
    
private:    
    void SetId(int idx){id = idx;}
    
    int id;
    double Position[3];
    double Time;
    double MassSq;
    
    ClassDefNV(DataTreeTOFSegment, 1)
};

#endif