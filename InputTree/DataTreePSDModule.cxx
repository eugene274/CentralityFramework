#include "DataTreePSDModule.h"
#include <iostream>
#include <vector>
#include "TObject.h"

DataTreePSDModule::DataTreePSDModule(int idx) : TObject(),
ProcessFlag(false),
arrSections (new TClonesArray("DataTreePSDSection")),
nSections(0)
{
    SetId(idx);
}

DataTreePSDModule::DataTreePSDModule(int idx, int inSections) : TObject(),
ProcessFlag(false),
arrSections (new TClonesArray("DataTreePSDSection")),
nSections(0)
{
    SetId(idx);
    for (int i=0;i<inSections;i++)
    {
	AddSection();
    }
    SectionNumberCut = nSections;
}

DataTreePSDModule::~DataTreePSDModule()
{
    
}

ClassImp(DataTreePSDModule)

