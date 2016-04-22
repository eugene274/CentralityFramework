
#ifndef CentralityManager_H
#define CentralityManager_H 1

#include <vector>
#include "TString.h"
#include "TNamed.h"
#include "TCut.h"

#include "CentralitySlicesFinder.h"
#include "CentralityGetter.h"
#include "NA61DataEvent.h"

class CentralityManager : public TNamed
{

    public:

	/**   Default constructor   **/
	CentralityManager() ;

	/**   Destructor   **/
	virtual ~CentralityManager();

	/**   Setters  **/
        void SetContainerFileName (TString ContainerFileName) { fContainerFileName = ContainerFileName; }
        void CopyNa61ExpDataToContainer (TString dir);
        void LoadDataFromContainer (Int_t Det1Id, Int_t Det2Id = -1); 
        void RunSliceFinder() { fSlicesFinder->RunSliceFinder(fRunId); }
        void SetRunId (Int_t RunId) { fRunId = RunId;}   
        void SetNumberOfSlices (Int_t num) { fSlicesFinder->SetNumberOfSlices(num); }
        void SetSliceStep (Float_t step) { fSlicesFinder->SetSliceStep(step); }
        
        void Det1IsInt (Bool_t is = true) { fSlicesFinder->Det1IsInt (is); }
	void Det2IsInt (Bool_t is = true) { fSlicesFinder->Det2IsInt (is); }
	void IsSimData (Bool_t is = true) { fSlicesFinder->IsSimData (is); }
        void SetNormalization (Float_t norm) { fSlicesFinder->SetNormalization(norm); }
        void SetCuts (TCut cut)           { fSlicesFinder->SetCuts(cut); }

        void AddDetector (TString DetName);
        void SetDetectorsForCentralityAnalisys (TString Det1Name, TString Det2Name = "");
        void SetCentralityMax (Float_t max) { fSlicesFinder->SetCentralityMax(max); }
        void SetDirectionCentralEvents (Int_t d) { d>0 ? fSlicesFinder->SetDirectionCentralEvents(1) : fSlicesFinder->SetDirectionCentralEvents(0); }
        void Do1DAnalisys (bool is = true) { fSlicesFinder->Do1DAnalisys(is); }
        void WriteCentralityFile () { fSlicesFinder->WriteOutputData(); }
        
        
        void LoadCentalityDataFile (TString file) { fCentralityGetter->LoadCentalityDataFile(file); }
        Float_t GetCentrality (Double_t det1) { return fCentralityGetter->GetCentrality(det1); }        
        Float_t GetCentrality (Double_t det1, Double_t det2) { return fCentralityGetter->GetCentrality(det1, det2); }
        
        void SetDirectory (TString dir) { fSlicesFinder->SetDir(dir); }
	/**   Getters  **/
        Int_t GetRunId () { return fRunId; }
        
    private:

	/**   Data members  **/
        TString fContainerFileName;
        
        CentralitySlicesFinder *fSlicesFinder;
        CentralityGetter *fCentralityGetter;
        NA61DataEvent *fNA61DataEvent;

        Int_t fRunId;
        std::vector <TString> fDetNames;
        
	ClassDef(CentralityManager, 2);

};


#endif
