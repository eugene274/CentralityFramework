#ifndef DataTreePSDHit_H
#define DataTreePSDHit_H 1

#include <vector>
#include <iostream>
class TClonesArray;

class DataTreePSDHit : public TObject
{
    
public:
  
    DataTreePSDHit();
    ~DataTreePSDHit();
    
    
    
    
private:

  PSD hit
    id
    energy
    module id

  
    ClassDef(DataTreePSDHit, 1)
};

#endif