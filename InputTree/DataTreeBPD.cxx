#include "DataTreeBPD.h"
#include <iostream>
#include <vector>
#include "TObject.h"

DataTreeBPD::DataTreeBPD(int idx) : TObject()
{
    SetId(idx);
}
DataTreeBPD::~DataTreeBPD()
{
    
}

ClassImp(DataTreeBPD)

