#ifndef CentralitySlicesFinder_H
#define CentralitySlicesFinder_H 1


#include <vector>
#include "TCut.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TNamed.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

#include "CentralityEventContainer.h"
#include "CentralitySlice.h"

// 

class CentralitySlicesFinder  : public TNamed
{
    public:
        
	/** Default constructor **/
	CentralitySlicesFinder();
//         CentralitySlicesFinder(const char* name, Int_t verbose = 1, Double_t En = 0.);

	/** Destructor **/
	virtual ~CentralitySlicesFinder();

        
        //====== Centrality 
        void LoadInputData (Int_t Det1Id, Int_t Det2Id = -1);

        void RunSliceFinder();
        void FindCentralitySlices ();        
        void Fit2DCorrelation ();
        void FitCorrection ( int n_points );
        void FindSlices ();
        void FindMeanSignalsInSlice ( void );
        void WriteOutputData ();
        
        void QA();
        
        void find_norm (Double_t par[], Double_t x, Double_t kb[], int N);
        
        /** Setters **/        
        void SetDet1Name               (TString Name)          { fDet1Name =               Name; }
        void SetDet2Name               (TString Name)          { fDet2Name =               Name; }
        void SetCuts                   (TCut Cuts)             { fCuts =                   Cuts; }
	void SetCentralityMax          (Float_t CentralityMax) { fCentralityMax =          CentralityMax; }
        void SetPrecision              (Float_t Precision)     { fPrecision =              Precision; }
	void SetDirectionCentralEvents (Int_t Direction)       { fDirectionCentralEvents = Direction; }
	void SetSliceStep              (Float_t SliceStep)     { fSliceStep =              SliceStep; }
	void SetNormalization          (Float_t Norm)          { fNormalization =          Norm; }
	void SetInFileName             (TString FileName)      { fInFileName =             FileName; }        
	void Do1DAnalisys              (bool is)               { is1DAnalisys =            is; }
	void Det1IsInt                 (Bool_t is)             { isDet1Int =               is; }
	void Det2IsInt                 (Bool_t is)             { isDet2Int =               is; }
        void SetNumberOfSlices         (Int_t num)             { nIntervals =              num; }
        void SetDir                    (TString dir)           { CFdir =                   dir;}
        void IsSimData                 (Bool_t is = true)      { fIsSimData =              is; }
        /** Getters **/
        TString GetDet1Name  () { return fDet1Name; }
        TString GetDet2Name  () { return fDet2Name; }
        Float_t GetSliceStep () { return fSliceStep; }
 
        
        Double_t polN (Double_t par[], Double_t x, int N)
        {
            Double_t res = par[0];
            Double_t xn = 1;
            for (int i=1; i<N;  i++){
                xn *= x;
                res += par[i]*xn;
            }
            return res;
        }       

    private:
        
        // Calculation parameters
        TString   fDet1Name, fDet2Name;
        TCut      fCuts;
        Float_t   fCentralityMax;
        Float_t   fPrecision;       
        Int_t     fDirectionCentralEvents;
        TString   fInFileName;
        Bool_t    is1DAnalisys;
        Bool_t    isDet1Int, isDet2Int;
        Float_t   fSliceStep;
        int       nIntervals;
        Bool_t    fIsSimData;
        Float_t   fNormalization;
        TString   CFdir;

        // Input
        TFile *fInFile;
        TTree *fInTree;
        CentralityEventContainer *fContainer;
        Float_t fDet1Id, fDet2Id;  

        // Run-by-run QA
        std::vector <Float_t> fDet1NormVec;
        std::vector <Float_t> fDet2NormVec;
        std::vector <Int_t>   fRunIdVec;
        std::vector <Int_t>   nEventsInRunVec;
        
        //Output
        CentralitySlice *fSlice;
        TTree           *fNormTree;
        TTree           *fCentrTree;
        TFile           *fSlicesFile;
        TF1             *fFitFunction;
        Float_t         det1max, det2max;
        TGraphErrors    *fCorrProfile;
        TCanvas         *fCanvas;
        Float_t         det1, det2, fB;

        static const int n_par = 4;      //TODO make params
        static const int nFinalPar = 6;

        
        CentralitySlicesFinder(const CentralitySlicesFinder&);
        CentralitySlicesFinder& operator = (const CentralitySlicesFinder&);

	ClassDef(CentralitySlicesFinder, 2);
};


#endif