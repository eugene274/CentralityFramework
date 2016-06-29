#ifndef CentralityContainerNormalizer_H
#define CentralityContainerNormalizer_H 1


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

class CentralityContainerNormalizer  : public TNamed
{
    public:
        
	/** Default constructor **/
	CentralityContainerNormalizer();
//         CentralityContainerNormalizer(const char* name, Int_t verbose = 1, Double_t En = 0.);

	/** Destructor **/
	virtual ~CentralityContainerNormalizer();

        
        //====== Centrality 
        void LoadInputData (Int_t Det1Id, Int_t Det2Id = -1);
        void GetNormalization (Int_t Det1Id, Int_t Det2Id );
        void GetNormalization (Int_t Det1Id);
	void RunByRunCorrection (Int_t Det1Id, Int_t Det2Id);
        void GetMaximum (Int_t Det1Id, Int_t Det2Id);
                        
        /** Setters **/        
        void SetDet1Name (TString Name)   { fDet1Name = Name; }
        void SetDet2Name (TString Name)   { fDet2Name = Name; }
	void SetNormalization (Float_t Norm) { fNormalization = Norm; }
	void SetInFileName (TString FileName)  {  fInFileName = FileName;   }        
	void Do1DAnalisys (bool is) { is1DAnalisys = is; }
	void Det1IsInt (Bool_t is) { isDet1Int = is; }
	void Det2IsInt (Bool_t is) { isDet2Int = is; }
        void IsSimData (Bool_t is = true)  { fIsSimData = is; }

	void SetDir (TString dir) {CFdir = dir;}

	/** Getters **/
        TString GetDet1Name () { return fDet1Name; }
        TString GetDet2Name () { return fDet2Name; }
        Float_t GetDet1Max () { return det1max; }
        Float_t GetDet2Max () { return det2max; }

        
        TTree *GetNormTree (void) { return fNormTree; }

        std::vector <Float_t> GetDet1NormVec    () {  return  fDet1NormVec; }
        std::vector <Float_t> GetDet2NormVec    () {  return  fDet2NormVec; }
        std::vector <Int_t>   GetRunIdVec       () {  return  fRunIdVec; }
        std::vector <Int_t>   GetEventsInRun2Vec () {  return  nEventsInRun2Vec; }
        std::vector <Int_t>   GetEventsInRun1Vec () {  return  nEventsInRun1Vec; }
        
        
//         Float_t GetCentrality (Double_t par1, Double_t par2 = 0.);


    private:
        
        //Calculation parameters
        TString fDet1Name, fDet2Name;
        Int_t fRunId;
        Bool_t is1DAnalisys;
        TString fInFileName;
        Bool_t isDet1Int, isDet2Int;
        Bool_t fIsSimData;
        Float_t fNormalization;
        TString CFdir;
        
        std::vector <Float_t> fDet1NormVec;
        std::vector <Float_t> fDet2NormVec;
        std::vector <Int_t>   fRunIdVec;
        std::vector <Int_t>   nEventsInRun2Vec;
        std::vector <Int_t>   nEventsInRun1Vec;

        TFile *fInFile;
        TTree *fInTree;
        
        CentralityEventContainer *fContainer;
        //output class

        TTree *fNormTree;
        Float_t det1, det2, det1max, det2max, fB;
        
        CentralityContainerNormalizer(const CentralityContainerNormalizer&);
        CentralityContainerNormalizer& operator = (const CentralityContainerNormalizer&);

	ClassDef(CentralityContainerNormalizer, 2);
};


#endif