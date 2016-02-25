
#ifndef CentralityDetectorEvent_H
#define CentralityDetectorEvent_H 1

#include <vector>
#include "TString.h"
#include "TNamed.h"

class CentralityDetectorEvent : public TNamed
{

    public:

	/**   Default constructor   **/
	CentralityDetectorEvent() ;

	/**   Destructor   **/
	virtual ~CentralityDetectorEvent() {};

	/**   Setters  **/
        void SetDetId (Int_t DetId) { fDetId = DetId;}   
        void SetEntries (std::vector <Float_t> weights) { fWeights = weights; }
        void AddEntry (Float_t weight) { fWeights.push_back(weight); }
        
        
	/**   Getters  **/
        std::vector <Float_t> GetWeightVector () { return fWeights; }
        Int_t GetDetId () { return fDetId; }
        Int_t GetNumberOfChannels () { return fWeights.size(); }
        
        Float_t GetTotalSignal () 
        { 
            Float_t signal = 0;
            for (unsigned int i=0; i<fWeights.size(); i++)
                signal += fWeights.at(i);
            return signal;
        }
                
    private:

	/**   Data members  **/
        Int_t fDetId;
        std::vector <Float_t> fWeights;
        
	ClassDef(CentralityDetectorEvent, 2);

};


#endif
