#include "DataTreeTrigger.h"
#include <iostream>
#include <vector>
#include "TObject.h"
#include "TString.h"

DataTreeTrigger::DataTreeTrigger(int idx) : TObject()
{
    SetId(idx);
}
DataTreeTrigger::DataTreeTrigger(int idx, TString label) : TObject()
{
    SetId(idx);
    SetLabel(label);
}
DataTreeTrigger::~DataTreeTrigger()
{
    
}

ClassImp(DataTreeTrigger)

