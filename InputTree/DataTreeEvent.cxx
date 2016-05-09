#include "DataTreeEvent.h"
#include <iostream>
#include <vector>
#include "TObject.h"
#include "DataTreeTrack.h"
#include "DataTreePSDModule.h"
#include "DataTreeTOFSegment.h"
#include "DataTreeV0Candidate.h"
#include "DataTreeMCTrack.h"
#include "DataTreeTrigger.h"
#include "DataTreeBPD.h"

DataTreeEvent::DataTreeEvent() : TObject(),
ProcessFlag(false),
arrTracks (new TClonesArray("DataTreeTrack")),
nTracks(0),
arrPSDModules (new TClonesArray("DataTreePSDModule")),
nPSDModules(0),
arrTOFSegments (new TClonesArray("DataTreeTOFSegment")),
nTOFSegments(0),
arrV0Candidates (new TClonesArray("DataTreeV0Candidate")),
nV0Candidates(0),
arrV0CandidatesMCpid (new TClonesArray("DataTreeV0Candidate")),
nV0CandidatesMCpid(0),
arrMCTracks (new TClonesArray("DataTreeMCTrack")),
nMCTracks(0),
arrTriggers (new TClonesArray("DataTreeTrigger")),
nTriggers(0),
arrBPDs (new TClonesArray("DataTreeBPD")),
nBPDs(0)
{
//     std::cout << "Constructor" << std::endl;
}
DataTreeEvent::~DataTreeEvent()
{

}

ClassImp(DataTreeEvent)

