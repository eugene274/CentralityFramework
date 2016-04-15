#include "DataTreeTOFSegment.h"
#include <iostream>
#include <vector>
#include "TObject.h"

DataTreeTOFSegment::DataTreeTOFSegment(int idx) : TObject()
{
    SetId(idx);
}
DataTreeTOFSegment::~DataTreeTOFSegment()
{
    
}

ClassImp(DataTreeTOFSegment)

