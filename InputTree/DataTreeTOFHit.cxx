#include "DataTreeTOFHit.h"
#include <iostream>
#include <vector>
#include "TObject.h"

DataTreeTOFHit::DataTreeTOFHit(int idx) : TObject()
{
    SetId(idx);
}
DataTreeTOFHit::~DataTreeTOFHit()
{
    
}

ClassImp(DataTreeTOFHit)

