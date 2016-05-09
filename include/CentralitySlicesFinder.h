#ifndef CentralitySlicesFinder_H
#define CentralitySlicesFinder_H 1


#include <vector>
#include "TCut.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TNamed.h"
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
        void GetNormalization (Int_t Det1Id, Int_t Det2Id );
        void GetNormalization (Int_t Det1Id);


        void RunSliceFinder();
        void FindCentralitySlices ();        
        void Fit2DCorrelation ();
        void FitCorrection ( int n_points );
        void FindSlices ();
        void FindMeanSignalsInSlice ( void );
        void WriteOutputData ();
        
        void QA();
        
        void find_norm (double par[], double x, double kb[], int N);
        
        /** Setters **/        
        void SetDet1Name (TString Name)   { fDet1Name = Name; }
        void SetDet2Name (TString Name)   { fDet2Name = Name; }
        void SetCuts (TCut Cuts) { fBaseCuts = Cuts; }
	void SetPrecision (Float_t Precision) { fPrecision = Precision; }
	void SetCentralityMax (Float_t CentralityMax) { fCentralityMax = CentralityMax; }
	void SetDirectionCentralEvents (Int_t Direction) { fDirectionCentralEvents = Direction; }
	void SetSliceStep (Float_t SliceStep) { fSliceStep = SliceStep; }
	void SetNormalization (Float_t Norm) { fNormalization = Norm; }
	

        
	void SetInFileName (TString FileName)  {  fInFileName = FileName;   }        
              
	
	void Do1DAnalisys (bool is) { is1DAnalisys = is; }
	void Det1IsInt (Bool_t is) { isDet1Int = is; }
	void Det2IsInt (Bool_t is) { isDet2Int = is; }
        void SetNumberOfSlices (Int_t num) { nIntervals = num; }
        void SetDir (TString dir) {CFdir = dir;}
        void IsSimData (Bool_t is = true) { fIsSimData = is; }
        /** Getters **/
        TString GetDet1Name () { return fDet1Name; }
        TString GetDet2Name () { return fDet2Name; }
        Float_t GetSliceStep () { return fSliceStep; }
 
        
        double polN (double par[], double x, int N)
        {
            double res = par[0];
            for (int i=1; i<N;  i++)
                res += par[i]*TMath::Power(x, i);
            return res;
        }       
        
//         Float_t GetCentrality (Double_t par1, Double_t par2 = 0.);


    private:
        
        //Calculation parameters
        TString fDet1Name, fDet2Name;
        TCut fCuts, fBaseCuts;
        Float_t fCentralityMax;
        Float_t fPrecision;       
        Int_t fDirectionCentralEvents;
        TString fInFileName;
        Bool_t is1DAnalisys;
        Bool_t isDet1Int, isDet2Int;
        Float_t fSliceStep;
        int nIntervals;
        Bool_t fIsSimData;
        Float_t fNormalization;
        TString CFdir;
        
        std::vector <Float_t> fDet1NormVec;
        std::vector <Float_t> fDet2NormVec;
        std::vector <Int_t>   fRunIdVec;
        std::vector <Int_t>   nEventsInRunVec;
        TFile *fInFile;
        TTree *fInTree;
        
        CentralityEventContainer *fContainer;
        //output class
        CentralitySlice *fSlice;
        static const int n_par = 4;
        static const int nFinalPar = 6;

        TTree *fNormTree, *fCentrTree;
        TFile *fSlicesFile;
        Float_t det1, det2, det1max, det2max, fB;
        TF1 *fFitFunction;
        
        CentralitySlicesFinder(const CentralitySlicesFinder&);
        CentralitySlicesFinder& operator = (const CentralitySlicesFinder&);

	ClassDef(CentralitySlicesFinder, 2);
};


#endif