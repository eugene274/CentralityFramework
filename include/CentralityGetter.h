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

        void SetNSlices (Int_t NSlices) { fNSlices = NSlices;}   
        Int_t GetNSlices () { return fNSlices; }

    private:

        int fNSlices;

        // borders of centrality slices y = kx + b
        
        CentralitySlice *fSlice;
        
        TTree  *fCentrTree;
        TFile *fCentrFile;
        
        CentralityGetter(const CentralityGetter&);
        CentralityGetter& operator=(const CentralityGetter&);

	ClassDef(CentralityGetter,2);
};


#endif