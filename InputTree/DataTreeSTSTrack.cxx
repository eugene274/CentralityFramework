#include "DataTreeTrack.h"
#include <iostream>
#include <vector>
#include "TObject.h"

DataTreeTrack::DataTreeTrack(int idx) : TObject()
{
    SetId(idx);
}
DataTreeTrack::~DataTreeTrack()
{
    
}

ClassImp(DataTreeTrack)

