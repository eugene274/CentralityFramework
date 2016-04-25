
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
    
    void Test();
    private:

	/**   Data members  **/
        DataTreeEvent *fEvent;
	
        TTree *fOutTree;
        TFile *fOutFile;
        CentralityEventContainer *fContainer;
        
        
        ClassDef(TreeInterface, 2);

};


#endif
