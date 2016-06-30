#include "DataTreeSTSTrack.h"
#include <iostream>
#include <vector>
#include "TObject.h"

DataTreeSTSTrack::DataTreeSTSTrack(int idx) : TObject()
{
    SetId(idx);
}
DataTreeSTSTrack::~DataTreeSTSTrack()
{
    
}

ClassImp(DataTreeSTSTrack)

