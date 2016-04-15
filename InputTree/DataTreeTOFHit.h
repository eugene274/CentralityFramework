#ifndef DataTreeTOFHit_H
#define DataTreeTOFHit_H 1

#include <vector>
#include <iostream>
class TClonesArray;

class DataTreeTOFHit : public TObject
{
    
public:
  
    DataTreeTOFHit();
    ~DataTreeTOFHit();
    
    
    
    
private:
      


  TOF hit
    id
    time
    m2
    hit position: x, y, z
    px, py, pz extrapolated from STS
    MC position: x, y, z
    MC px, py, pz
    MC pdg_code
    

  
    ClassDefNV(DataTreeTOFHit, 1)
};

#endif