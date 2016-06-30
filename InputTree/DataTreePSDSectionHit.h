#ifndef DataTreePSDSectionHit_H
#define DataTreePSDSectionHit_H 1

#include <vector>
#include <iostream>
class TClonesArray;

class DataTreePSDSectionHit : public TObject
{
    
public:
  
    DataTreePSDSectionHit();
    ~DataTreePSDSectionHit();
    
    
    
    
private:

  PSD hit
    id
    energy
    module id

  
    ClassDef(DataTreePSDSectionHit, 1)
};

#endif