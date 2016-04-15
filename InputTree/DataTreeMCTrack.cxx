#include "DataTreeMCTrack.h"
#include <iostream>
#include <vector>
#include "TObject.h"

DataTreeMCTrack::DataTreeMCTrack(int idx) : TObject()
{
    SetId(idx);
}
DataTreeMCTrack::~DataTreeMCTrack()
{
    
}

ClassImp(DataTreeMCTrack)

