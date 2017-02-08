
#ifndef TreeInterface_H
#define TreeInterface_H 1

#include <vector>
#include "TNamed.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include "TChain.h"

#include "DataTreeEvent.h"
#include "CentralityEventContainer.h"

class TreeInterface : public TNamed
{

    public:
    TreeInterface ();
    std::vector <Float_t> SetPsdVector(Int_t subgroup);
    
    void WriteCentralityContainer();

    void SetOutFileName (TString OutFileName) {fOutFileName = OutFileName;}
    void SetInFileName  (TString InFileName)  {fInFileName = InFileName;}
    void SetNEntries (int n) { nEntries = n; }

    void SetInTChain  (TChain* chain)  { fInTree = (TTree*)chain; }

private:
    Float_t Beam_Eta;
    TString fOutFileName;
    TString fInFileName;
    
/**   Data members  **/
    DataTreeEvent *fEvent;
    
    TFile *fInFile ;
    TTree *fInTree ;    
    
    TTree *fOutTree;
    TFile *fOutFile;
    CentralityEventContainer *fContainer;
    
    TString fPsdGeomConfig;
    
    Int_t nEntries;
    
    ClassDef(TreeInterface, 2);

};


#endif
