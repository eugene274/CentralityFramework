#include "CentralityManager.h"

ClassImp(CentralityManager)

// -----   Default constructor   -------------------------------------------
CentralityManager::CentralityManager() 
  : TNamed(),
    fContainerFileName(""),
    fSlicesFinder(new CentralitySlicesFinder),
    fCentralityGetter(new CentralityGetter),
    fNA61DataEvent(new NA61DataEvent)
{
}

CentralityManager::~CentralityManager() 
{
    if ( fSlicesFinder ) {
	delete fSlicesFinder;
    }

    if ( fCentralityGetter ) {
	delete fCentralityGetter;
    }

    if ( fNA61DataEvent ) {
	delete fNA61DataEvent;
    }
}

void CentralityManager::CopyNa61ExpDataToContainer (TString dir)
{
    for (unsigned int i = 0; i<fDetNames.size(); i++)
    {
        std::cout <<  fDetNames.at(i) << std::endl;
        if  (fDetNames.at(i) == "PSD1")       { fNA61DataEvent->isPSD1 = true;  fNA61DataEvent->PSD1id = i; }
        else if  (fDetNames.at(i) == "PSD2")  { fNA61DataEvent->isPSD2 = true;  fNA61DataEvent->PSD2id = i; }
        else if  (fDetNames.at(i) == "PSD3")  { fNA61DataEvent->isPSD3 = true;  fNA61DataEvent->PSD3id = i; }
        else if  (fDetNames.at(i) == "TPC")   { fNA61DataEvent->isTPC = true;   fNA61DataEvent->TPCid = i; }
        else {
            std::cout << "*** CopyNa61ExpDataToContainer *** : Detector name " << fDetNames.at(i) << " is not supported for NA61. Please choose one of: PSD1, PSD2, PSD3, TPC" << std::endl;
        }
    }

    fNA61DataEvent->SetData ( dir, 0 );
    fNA61DataEvent->fOutFile = fContainerFileName;
    fNA61DataEvent->Loop();
}


void CentralityManager::LoadDataFromContainer (Int_t Det1Id, Int_t Det2Id ) 
{ 
    fSlicesFinder->SetInFileName(fContainerFileName);
    fSlicesFinder->LoadInputData(Det1Id, Det2Id); 
    if (Det2Id == -1) fSlicesFinder->Do1DAnalisys(true);
} 

void  CentralityManager:: AddDetector (TString DetName)
{
    fDetNames.push_back (DetName);
    std::cout << "Detector " << DetName << " added with DetId = " << fDetNames.size()-1 << std::endl;
    
}

void CentralityManager::SetDetectorsForCentralityAnalisys (TString Det1Name, TString Det2Name)
{

    fSlicesFinder->SetDet1Name(Det1Name);
    fSlicesFinder->SetDet2Name(Det2Name);
    
    Int_t id1 = -1, id2 = -1;
    for (unsigned int i = 0; i<fDetNames.size(); i++)
    {
        if (fDetNames.at(i) == Det1Name) id1 = i;
        if (fDetNames.at(i) == Det2Name) id2 = i;        
    }    
    
    if (id1 == -1)
        std::cout << "*** SetDetectorsForCentralityAnalisys *** : Detector name " << Det1Name << " is not supported." << std::endl;
        
    if (id2 == -1 && Det2Name != "")
        std::cout << "*** SetDetectorsForCentralityAnalisys *** : Detector name " << Det1Name << " is not supported." << std::endl;
    
    LoadDataFromContainer(id1, id2);
}