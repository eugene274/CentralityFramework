#ifndef DataTreePSDModuleHit_H
#define DataTreePSDModuleHit_H 1

#include <vector>
#include <iostream>
class TClonesArray;

class DataTreePSDModuleHit : public TObject
{
    
public:
  
    DataTreePSDModuleHit();
    ~DataTreePSDModuleHit();
    
    
    
    
private:

  PSD hit
    id
    energy
    module id

  
    ClassDef(DataTreePSDModuleHit, 1)
};

#endif