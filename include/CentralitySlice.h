
#ifndef CentralitySlice_H
#define CentralitySlice_H 1

#include <vector>
#include "TString.h"
#include "TNamed.h"
#include "TF1.h"

class CentralitySlice : public TNamed
{

    public:

	/**   Default constructor   **/
	CentralitySlice() ;

	/**   Destructor   **/
	virtual ~CentralitySlice() {};
        void ClearData();
	/**   Setters  **/
        void SetRunId (Int_t RunId) { fRunId = RunId;}   
        void AddSlice (Float_t A, Float_t B, Float_t X);
        void SetDetMax (Float_t det1, Float_t det2 = -1) { fDet1Max = det1;  fDet2Max = det2; }
        void AddXPar (Float_t mean, Float_t sigma) { MeanX.push_back(mean);  SigmaX.push_back(sigma);  }
        void AddYPar (Float_t mean, Float_t sigma) { MeanY.push_back(mean);  SigmaY.push_back(sigma);  }
        void AddXYPar (Float_t mean, Float_t sigma) { MeanXY.push_back(mean);  SigmaXY.push_back(sigma);  }

        void AddXY3 (Float_t x3) { MeanXY3.push_back(x3);  }
        void AddB (Float_t mean, Float_t sigma, Float_t dmean = 0, Float_t dsigma = 0) { MeanB.push_back(mean);  SigmaB.push_back(sigma); dB.push_back(dmean); dSigmaB.push_back(dsigma);  }
        
        void SetNSlices (Int_t NSlices) { fNSlices = NSlices;}   
        void SetCentralityMax (Float_t CentralityMax) { fCentralityMax = CentralityMax;}   
        void SetFitFunction (TF1 *FitFunction) { fFitFunction = FitFunction; }
        void SetSlicesStep (Float_t step) { fStep = step;}   
        void SetDirectionCentralEvents (Int_t d) {fDirectionCentralEvents = d;}

        
        void SetDet1NormVec (std::vector <Float_t> v) { fDet1NormVec = v; }        
        void SetDet2NormVec (std::vector <Float_t> v) { fDet2NormVec = v; }        
        void SetRunIdVec    (std::vector <Int_t> v) { fRunIdVec    = v; }  
        
        /**   Getters  **/
        Int_t GetRunId () { return fRunId; }
        Float_t GetAi (Int_t i) { return  fAvec.at(i); }
        Float_t GetBi (Int_t i) { return  fBvec.at(i); }
        Float_t GetXi (Int_t i) { return  fXvec.at(i); }

        Float_t GetMeanXi (Int_t i) { return  MeanX.at(i); }
        Float_t GetMeanYi (Int_t i) { return  MeanY.at(i); }
        Float_t GetMeanXYi (Int_t i) { return  MeanXY.at(i); }

        Float_t GetSigmaXi (Int_t i) { return  SigmaX.at(i); }
        Float_t GetSigmaYi (Int_t i) { return  SigmaY.at(i); }
        Float_t GetSigmaXYi (Int_t i) { return  SigmaXY.at(i); }
        
        Float_t GetDet1Max () { return fDet1Max; }
        Float_t GetDet2Max () { return fDet2Max; }
        
        UInt_t GetNSlices () { return fAvec.size(); }
        Float_t GetCentralityMax () { return fCentralityMax; }
        Float_t GetSlicesStep () { return fStep; }
        
        Int_t GetDirectionCentralEvents() {return fDirectionCentralEvents;}
        
        TF1* GetFitFunction () { return fFitFunction; }

        std::vector <Float_t> GetA () { return fAvec; }
        std::vector <Float_t> GetB () { return fBvec; }
        std::vector <Float_t> GetX () { return fXvec; }

        std::vector <Float_t> GetMeanX () { return MeanX; }
        std::vector <Float_t> GetMeanY () { return MeanY; }
        std::vector <Float_t> GetMeanXY () { return MeanXY; }

        std::vector <Float_t> GetSigmaX () { return SigmaX; }
        std::vector <Float_t> GetSigmaY () { return SigmaY; }
        std::vector <Float_t> GetSigmaXY () { return SigmaXY; }
        
        std::vector <Float_t> GetMeanB () { return MeanB; }
        std::vector <Float_t> GetSigmaB () { return SigmaB; }
        std::vector <Float_t> GetdB ()     { return dB; }
        std::vector <Float_t> GetdSigmaB () { return dSigmaB; }
        std::vector <Float_t> GetMeanXY3 () { return MeanXY3; }
        
        std::vector <Float_t> GetDet1NormVec () { return fDet1NormVec; }
        std::vector <Float_t> GetDet2NormVec () { return fDet2NormVec; }
        std::vector <Int_t> GetRunIdVec    () { return fRunIdVec; }
        
        
        
    private:

	/**   Data members  **/
        Int_t fRunId;

        Float_t fDet1Max;
        Float_t fDet2Max;

        Int_t fNSlices;
        Float_t fCentralityMax;
        Float_t fStep;
        Int_t fDirectionCentralEvents;

        std::vector <Float_t> fAvec;
        std::vector <Float_t> fBvec;
        std::vector <Float_t> fXvec;

        std::vector <Float_t> MeanX;
        std::vector <Float_t> MeanY;
        std::vector <Float_t> MeanXY;

        std::vector <Float_t> SigmaX;
        std::vector <Float_t> SigmaY;
        std::vector <Float_t> SigmaXY;

        std::vector <Float_t> MeanXY3;

        std::vector <Float_t> MeanB;
        std::vector <Float_t> SigmaB;
        std::vector <Float_t> dB;
        std::vector <Float_t> dSigmaB;
        
        std::vector <Float_t> fDet1NormVec;
        std::vector <Float_t> fDet2NormVec;
        std::vector <Int_t>   fRunIdVec;


        TF1 *fFitFunction;
        
        ClassDef(CentralitySlice, 2);

};


#endif
