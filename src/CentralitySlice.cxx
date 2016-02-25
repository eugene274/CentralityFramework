#include "CentralitySlice.h"

ClassImp(CentralitySlice)

// -----   Default constructor   -------------------------------------------
CentralitySlice::CentralitySlice() 
  : TNamed(),
    fRunId(0),
    fFitFunction (new TF1)
{
}

void CentralitySlice::AddSlice (Float_t A, Float_t B, Float_t X)
{
    fAvec.push_back (A);
    fBvec.push_back (B);
    fXvec.push_back (X);
}

void CentralitySlice::ClearData()
{
    fAvec.clear();
    fBvec.clear();
    fXvec.clear();
    MeanX.clear();
    MeanY.clear();
    MeanXY.clear();
    SigmaX.clear();
    SigmaY.clear();
    SigmaXY.clear();
    MeanXY3.clear();

    
}