#include "DataTreeEvent.h"
#include <iostream>
#include <vector>
#include "TObject.h"
#include "DataTreeTrack.h"
#include "DataTreePSDModule.h"
#include "DataTreeTOFHit.h"
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
arrTOFHits (new TClonesArray("DataTreeTOFHit")),
nTOFHits(0),
arrV0CandidatesTOFpid (new TClonesArray("DataTreeV0Candidate")),
nV0CandidatesTOFpid(0),
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

