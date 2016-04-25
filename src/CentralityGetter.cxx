
#include <iostream>
#include <fstream>

#include "CentralityGetter.h"

#include "TCanvas.h"
#include "TRandom.h"

using std::vector;
using std::cout;
using std::endl;
using std::ifstream;

ClassImp(CentralityGetter)

CentralityGetter::CentralityGetter() :
fNSlices (20),
isDet1Int(false),
isDet2Int(false),
fSlice (new CentralitySlice)
{
}

CentralityGetter::~CentralityGetter()
{
}


void CentralityGetter::LoadCentalityDataFile (TString fileName)
{

    fCentrFile = new TFile(fileName.Data(), "READ");
    fCentrTree = (TTree*) fCentrFile->Get("CentrTree");    
    fCentrTree->SetBranchAddress("CentralitySlice", &fSlice);   
}

Float_t CentralityGetter::GetCentrality (Double_t det1)
{    
    float centrality = -1;
    fCentrTree->GetEntry(0);
    TRandom* random = new TRandom;  
    
    if (isDet1Int){ 
        Float_t rand1 = random->Rndm()/* - 0.5*/;  /*cout << rand1 << endl;*/  
        cout << "det1 = " << det1 << endl;
        det1 += rand1; 
        cout << "  det1 = " << det1 << endl;

    }
    
    det1 /= fSlice->GetDet1Max();
    bool isOK = false;

    Float_t step = fSlice->GetSlicesStep ();
    Int_t NSlices = fSlice->GetNSlices();
    
    cout << "  det1 = " << det1 << endl;
//     cout << "NSlices = " << NSlices << endl;
//     cout << "step = " << step << endl;
    
    Int_t direction = fSlice->GetDirectionCentralEvents();
    
    if ( (direction == 1 && det1 > fSlice->GetXi(0)) || (direction == 0 && det1 < fSlice->GetXi(0)) )  //most central events
    {
        centrality = step*(0.5);
        return centrality;
    }

    if ( (direction == 1 && det1 < fSlice->GetXi(NSlices-1)) || (direction == 0 && det1 > fSlice->GetXi(NSlices-1)) )  //most peripheral events
    {
        centrality = step*(NSlices+0.5);
        return centrality;
    }    
    
    
    for (UInt_t i=1; i<=NSlices; i++)
    {
//         cout << "x = " << fSlice->GetXi(i) << endl;
        if ( (det1-fSlice->GetXi(i))*(det1-fSlice->GetXi(i-1)) <= 0 )
        {
            centrality = step*(i+0.5);
            isOK = true;
            break;
        }
    }

//     if (!isOK)
//         cout << "*** Warning *** GetCentrality: centrality is not found!" << endl;

    
    return centrality;
}

Float_t CentralityGetter::GetCentrality (Double_t det1, Double_t det2)
{    
    float centrality = -1;

    fCentrTree->GetEntry(0);
    Float_t step = fSlice->GetSlicesStep ();
    Int_t NSlices = fSlice->GetNSlices();
    TRandom* random = new TRandom;  
    
    if (isDet1Int){ 
        Float_t rand1 = random->Rndm()/* - 0.5*/;  /*cout << rand1 << endl;*/  
        det1 += rand1; 
    }
    if (isDet2Int){ 
        Float_t rand2 = random->Rndm()/* - 0.5*/;  /*cout << rand1 << endl;*/  
        det2 += rand2; 
    }
    det1 /= fSlice->GetDet1Max();
    det2 /= fSlice->GetDet2Max();    
    bool isOK = false;

    Float_t x0 = (det2-fSlice->GetAi(0))/fSlice->GetBi(0);
    Float_t xn = (det2-fSlice->GetAi(NSlices-1))/fSlice->GetBi(NSlices-1);
    Int_t direction = fSlice->GetDirectionCentralEvents();
    
    if ( (direction == 1 && det1 > x0) || (direction == 0 && det1 < x0) )  //most central events
    {
        centrality = step*(0.5);
        return centrality;
    }
    if ( (direction == 1 && det1 < x0) || (direction == 0 && det1 > x0) )  //most peripheral events
    {
        centrality = step*(NSlices+0.5);
        return centrality;
    }    
    
    for (UInt_t i=1; i<NSlices; i++)
    {
        
        Float_t xi = (det2-fSlice->GetAi(i))/fSlice->GetBi(i);
        Float_t xi_1 = (det2-fSlice->GetAi(i-1))/fSlice->GetBi(i-1);

//         cout << "xi = " << xi << "  xi+1 = " << xi_1  << "   det1 = " << det1 << "   det2 = " << det2  << endl;

        if ( (det1-xi)*(det1-xi_1) <= 0 )
        {
            centrality = step*(i+0.5);
            isOK = true;
            break;
        }
    }        
        
//     if (!isOK)
//         cout << "*** Warning *** GetCentrality: centrality is not found!" << endl;
    
    return centrality;
}

