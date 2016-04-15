#include "DataTreePSDSection.h"
#include <iostream>
#include <vector>
#include "TObject.h"

DataTreePSDSection::DataTreePSDSection(int idx) : TObject()
{
    SetId(idx);
}
DataTreePSDSection::~DataTreePSDSection()
{
    
}

ClassImp(DataTreePSDSection)

