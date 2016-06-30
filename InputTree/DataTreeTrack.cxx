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

void DataTreeTrack::SetUndefinedValues()
{
    double fValue = -999;//TEST
    for (int i=0;i<nZPositions;i++)
    {
	pT[i] = fValue;
	phi[i] = fValue;
	eta[i] = fValue;
	px[i] = fValue;
	py[i] = fValue;
	pz[i] = fValue;
	p[i] = fValue;

	VtxChiSq[i] = fValue;
	
	NofHits[i][nSubDetectors] = fValue;
	NofHitsPotential[i][nSubDetectors] = fValue;
	Flag[i] = fValue;
	ChiSq[i] = fValue;
	NDF[i] = fValue;
	for (int j=0;j<3;j++){ DCA[i][j] = fValue; }
	Stations[i][nMaxStations] = false;
	nStations[i] = fValue;
	nSTSHitsPossible[i] = fValue;
	LengthInSTS[i] = fValue;
	Charge[i] = fValue;
    }
    PSDModuleId = fValue;
    TOFSegmentId = fValue;
    MCTrackId = fValue;	
}