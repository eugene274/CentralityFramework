
#ifndef CentralityEventContainer_H
#define CentralityEventContainer_H 1

#include <vector>
#include "TString.h"
#include "TNamed.h"

#include "CentralityDetectorEvent.h"


class CentralityEventContainer : public TNamed
{

    public:

	/**   Default constructor   **/
	CentralityEventContainer() ;

	/**   Destructor   **/
	virtual ~CentralityEventContainer() {};

	/**   Setters  **/
        void SetRunId (Int_t RunId) { fRunId = RunId;}   
        void SetB     (Float_t B)     { fB = B;}           
        void SetEventId (Int_t EventId) { fEventId = EventId;}   
        void AddDetector (CentralityDetectorEvent Det) { fDetectorEvents.push_back(Det); }
        void AddDetectorEntry (Int_t i, Float_t Entry) { fDetectorEvents.at(i).AddEntry(Entry); }
        void AddDetectorEvent (Int_t i, std::vector <Float_t> weights) { fDetectorEvents.at(i).SetEntries (weights) ;   fDetectorEvents.at(i).SetDetId(i); }

        /**   Getters  **/
        Int_t GetRunId () { return fRunId; }
        Int_t GetNumberOfDetectors () { return fDetectorEvents.size(); }
        CentralityDetectorEvent GetDetectorData (Int_t i) { return  fDetectorEvents.at(i); }
        Float_t GetB () { return fB; }
        Int_t GetEventId () { return fEventId; }
        Float_t GetDetectorWeight (Int_t i) { return (fDetectorEvents.at(i)).GetTotalSignal (); }
        
    private:

	/**   Data members  **/
        Int_t fRunId;
        Float_t fB;
        Int_t fEventId;
//         std::vector <TString> fDetectorNames;        
        std::vector <CentralityDetectorEvent> fDetectorEvents;
        
	ClassDef(CentralityEventContainer, 2);

};


#endif
