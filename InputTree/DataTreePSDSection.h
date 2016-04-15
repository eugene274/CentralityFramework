#ifndef DataTreePSDSection_H
#define DataTreePSDSection_H 1

#include <vector>
#include <iostream>
#include "TClonesArray.h"
#include "TObject.h"

class DataTreePSDSection : public TObject
{
    
public:
  
    DataTreePSDSection(int idx = 0);
    ~DataTreePSDSection();
    
    void ClearEvent(){Energy = 0.;}
    
    int GetId(){return id;}
    double GetEnergy(){return Energy;}

    void SetEnergy(double fValue){Energy = fValue;}
    void AddEnergy(double fValue){Energy += fValue;}
    
private:    
    void SetId(int idx){id = idx;}
    
    int id;
    double Energy;
    
    ClassDefNV(DataTreePSDSection, 1)
};

#endif