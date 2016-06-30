#include "DataTreeV0Candidate.h"
#include <iostream>
#include <vector>
#include "TObject.h"

DataTreeV0Candidate::DataTreeV0Candidate(int idx) : TObject(),
    arrDaughters (new TClonesArray("DataTreeV0Candidate")),
    nDaughters(0)
{
    SetId(idx);
}
DataTreeV0Candidate::~DataTreeV0Candidate()
{
    
}

ClassImp(DataTreeV0Candidate)

