
#ifndef GlauberFitter_H
#define GlauberFitter_H 1

#include <vector>
#include "TString.h"
#include "TNamed.h"
#include "TH1F.h"
#include "TTree.h"
#include "TMinuit.h"

class GlauberFitter : public TNamed
{

    public:

	/**   Default constructor   **/
	GlauberFitter() ;
        void TestFunc(Int_t nf, Float_t f0, Float_t f1, Int_t nsigma, Int_t nEvents);
        
        void SetSimHistos (TString filename, Int_t nEntries = -1);
        void SetGlauberFitHisto (Float_t f, Float_t mu, Float_t k, Int_t n = 10000, Bool_t Norm2Data = true);
        void NormalizeGlauberFit (Int_t StartBin = 30);
        void DrawHistos (Bool_t isSim = true, Bool_t isData = true, Bool_t isGlauber = false, Bool_t isNBD = false);
	
        void FitGlauber (Int_t nf, Float_t f0, Float_t f1, Int_t nsigma, Float_t SigmaStep, Int_t nEvents);
        void FindMuGoldenSection (Float_t *mu, Float_t *chi2, Float_t mu_min, Float_t mu_max, Float_t f, Float_t k, Int_t nEvents = 10000, Int_t nIter = 5);
        
        Double_t GetChi2 (void);
        
        Double_t NBD(Double_t n, Double_t mu, Double_t k) const;
        void SetNBDhist(Double_t mu, Double_t k);
        /**   Destructor   **/
	virtual ~GlauberFitter() {};
        
        
        void SetInputHisto (TH1F *h) { hData = h; }
        
        void SetNpartMax (Int_t max)   { fNpartMax =   max; }
        void SetNcollMax (Int_t max)   { fNcollMax =   max; }
        void SetMultMax  (Int_t max)   { fMultMax = max; }
        void SetFitMultMin  (Int_t min)   { fFitMultMin = min; }
        void SetNormMultMin  (Int_t min)   { fNormMultMin = min; }

        void SetBinSize  (Int_t size)   { fBinSize = size; }
        void SetOutDirName (TString name) { fOutDirName = name; }
        TH1F *GetGlauberFitHisto () { return hGlaub; }
        TH1F *GetDataHisto ()       { return hData;  }
        TH1F *GetNBDHisto ()        { return hNBD;   }
        
        void SetMode (TString mode) { fMode = mode; }
        
    private:

	/**   Data members  **/
        TH1F *hNpart, *hNcoll, *hData, *hNBD, *hGlaub, *hBestFit;
        TTree* fSimTree;
        
        Float_t fNpart, fNcoll;
        Int_t fNpartMax, fNcollMax;
        Int_t fMultMax;
        Int_t fBinSize;
        
        Int_t fFitMultMin;
        Int_t fNormMultMin;

        TString fMode;
        
        TString fOutDirName;
        ClassDef(GlauberFitter, 2);

};


#endif
