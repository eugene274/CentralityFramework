#ifndef DataTreeBPD_H
#define DataTreeBPD_H 1

#include <vector>
#include <iostream>
#include "TClonesArray.h"
#include "TObject.h"

class DataTreeBPD : public TObject
{
    
public:
  
    DataTreeBPD(int idx = 0);
    ~DataTreeBPD();
    
    int GetId(){return id;}
    
    double GetPositionComponent(int idx){return Position[idx];}

    void SetPosition(double fX, double fY, double fZ){Position[0]=fX; Position[1]=fY; Position[2]=fZ;}
    void SetPositionComponent(int idx, double fValue){Position[idx] = fValue;}
    
    
private:    
    void SetId(int idx){id = idx;}
    
    int id;
    double Position[3];
    
    ClassDefNV(DataTreeBPD, 1)
};

#endif