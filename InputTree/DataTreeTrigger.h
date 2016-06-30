#ifndef DataTreeTrigger_H
#define DataTreeTrigger_H 1

#include <vector>
#include <iostream>
#include "TClonesArray.h"
#include "TObject.h"
#include "TString.h"

class DataTreeTrigger : public TObject
{
    
public:
  
    DataTreeTrigger(int idx = 0);
    DataTreeTrigger(int idx, TString label);
    ~DataTreeTrigger();
    
    int GetId(){return id;}
    
    double GetSignal(){return Signal;}
    TString GetLabel(){return Label;}
    
    void SetSignal(double fValue){Signal = fValue;}
    void SetLabel(TString fLabel){Label = fLabel;}
    
private:    
    void SetId(int idx){id = idx;}
    
    int id;
    double Signal;
    TString Label;
    
    ClassDefNV(DataTreeTrigger, 1)
};

#endif