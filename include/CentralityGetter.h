#ifndef CentralityGetter_H
#define CentralityGetter_H 1


#include <vector>
#include "TCut.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "CentralitySlice.h"

// 

class CentralityGetter 
{
    public:
        
	/** Default constructor **/
	CentralityGetter();
//         CentralityGetter(const char* name);

	/** Destructor **/
	virtual ~CentralityGetter();

        void LoadCentalityDataFile (TString fileName);
        Float_t GetCentrality (Double_t det1);        
        Float_t GetCentrality (Double_t det1, Double_t det2);

        void Det1IsInt (Bool_t is) { isDet1Int = is; }
	void Det2IsInt (Bool_t is) { isDet2Int = is; }
        
        void SetNSlices (Int_t NSlices) { fNSlices = NSlices;}   
        Int_t GetNSlices () { return fNSlices; }

        void SetRunId (Int_t RunId);
        
    private:

        int fNSlices;
        Bool_t isDet1Int, isDet2Int;
        Int_t fRunId;

        // borders of centrality slices y = kx + b
        
        CentralitySlice *fSlice;
        
        TTree  *fCentrTree;
        TFile *fCentrFile;
        Float_t fDet1Norm, fDet2Norm;

        CentralityGetter(const CentralityGetter&);
        CentralityGetter& operator=(const CentralityGetter&);

	ClassDef(CentralityGetter,2);
};


#endif