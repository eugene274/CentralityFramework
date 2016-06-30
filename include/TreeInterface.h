
#ifndef TreeInterface_H
#define TreeInterface_H 1

#include <vector>
#include "TNamed.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"



#include "DataTreeEvent.h"
#include "CentralityEventContainer.h"

class TreeInterface : public TNamed
{

    public:
    TreeInterface ();
    std::vector <Float_t> SetPsdVector(Int_t subgroup);
    
    void WriteCentralityContainer();

    bool isGoodEvent(Int_t nTPC_Tracks_Ref);
    bool isRefMultTrack(int iTrk);

    void SetOutFileName (TString OutFileName) {fOutFileName = OutFileName;}
    void SetInFileName  (TString InFileName)  {fInFileName = InFileName;}
    
private:
    Float_t Beam_Eta;
    TString fOutFileName;
    TString fInFileName;
    
/**   Data members  **/
    DataTreeEvent *fEvent;
    
    TTree *fOutTree;
    TFile *fOutFile;
    CentralityEventContainer *fContainer;
    
    TString fPsdGeomConfig;
    
    ClassDef(TreeInterface, 2);

};


#endif
