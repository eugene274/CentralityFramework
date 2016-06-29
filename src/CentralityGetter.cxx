
#include <iostream>
#include <fstream>

#include "CentralityGetter.h"

#include "TCanvas.h"
#include "TRandom.h"
#include "TMath.h"
#include "TGraphErrors.h"



using std::vector;
using std::cout;
using std::endl;
using std::ifstream;

ClassImp(CentralityGetter)

CentralityGetter::CentralityGetter() :
fNSlices (20),
isDet1Int(false),
isDet2Int(false),
fSlice (new CentralitySlice),
fDet1Norm (-1),
fDet2Norm (-1)
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

void CentralityGetter::SetRunId (Int_t RunId)
{
    fCentrTree->GetEntry(0);    
    UInt_t nRuns = fSlice->GetRunIdVec().size();
    Bool_t isOK = false;
    
//     std::cout << nRuns << std::endl;
    
    for (UInt_t i=0; i<nRuns; i++)
    {
        if (fSlice->GetRunIdVec().at(i) == RunId)
        {
//             std::cout << fSlice->GetRunIdVec().at(i) << std::endl;
            isOK = true;
            fDet1Norm = fSlice->GetDet1NormVec().at(i);
            if ( fSlice->GetDet2NormVec().size() > i)   fDet2Norm = fSlice->GetDet2NormVec().at(i);
            fRunId = RunId;
            break;
        }
    }
    if (!isOK)   std::cout << "*** Warning! *** Corrections for run # " << RunId << " was not found! Uncorrected value will be used!" << std::endl;
    if (isOK)    std::cout << "   Correction factor(s) for run # " << RunId << " = " << fDet1Norm ;
    if (isOK && fDet2Norm != -1)    std::cout << "  and  " << fDet2Norm ;
    if (isOK)  std::cout << std::endl;
    
}



Float_t CentralityGetter::GetCentrality (Double_t det1)
{    
    float centrality = -1;
    fCentrTree->GetEntry(0);
    TRandom* random = new TRandom;  
    random->SetSeed();
    if (isDet1Int){ 
        Float_t rand1 = random->Rndm()/* - 0.5*/;  /*cout << rand1 << endl;*/  
//         cout << "det1 = " << det1 << endl;
        det1 += rand1; 
//         cout << "  det1 = " << det1 << endl;

    }
    
    if (fDet1Norm != -1)  det1 *= fDet1Norm;
    det1 /= fSlice->GetDet1Max();
    bool isOK = false;

    Float_t step = fSlice->GetSlicesStep ();
    Int_t NSlices = fSlice->GetNSlices();
    
//     cout << "  det1 = " << det1 << endl;
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
    
    
    for (Int_t i=1; i<=NSlices; i++)
    {
//         cout << "x = " << fSlice->GetXi(i) << endl;
        if ( (det1-fSlice->GetXi(i))*(det1-fSlice->GetXi(i-1)) <= 0 )
        {
            centrality = step*(i+0.5);
            isOK = true;
            break;
        }
    }

    if (!isOK)
        cout << "*** Warning *** GetCentrality: centrality is not found!" << endl;

    
    return centrality;
}

Float_t CentralityGetter::GetCentrality (Double_t det1, Double_t det2)
{    
    float centrality = -1;

    fCentrTree->GetEntry(0);
    Float_t step = fSlice->GetSlicesStep ();
    Int_t NSlices = fSlice->GetNSlices();
    TRandom* random = new TRandom;  
    random->SetSeed();
    if (isDet1Int){ 
        Float_t rand1 = random->Rndm()/* - 0.5*/;  /*cout << rand1 << endl;*/  
        det1 += rand1; 
    }
    if (isDet2Int){ 
        Float_t rand2 = random->Rndm()/* - 0.5*/;  /*cout << rand1 << endl;*/  
        det2 += rand2; 
    }
    
    if (fDet1Norm != -1)  det1 *= fDet1Norm;
    if (fDet2Norm != -1)  det2 *= fDet2Norm;
    
    det1 /= fSlice->GetDet1Max();
    det2 /= fSlice->GetDet2Max();    
    bool isOK = false;

//     cout  << "   det1 = " << det1 << "   det2 = " << det2  << endl;
        
    Float_t x0 = (det2-fSlice->GetAi(0))/fSlice->GetBi(0);
    Float_t xn = (det2-fSlice->GetAi(NSlices-1))/fSlice->GetBi(NSlices-1);
    Int_t direction = fSlice->GetDirectionCentralEvents();

//     cout  << "   x0 = " << x0 << "   xn = " << xn  << endl;

    
    if ( (direction == 1 && det1 > x0) || (direction == 0 && det1 < x0) )  //most central events
    {
        centrality = step*(0.5);
        return centrality;
    }
    if ( (direction == 1 && det1 < xn) || (direction == 0 && det1 > xn) )  //most peripheral events
    {
        centrality = step*(NSlices+0.5);
        return centrality;
    }    
    
    for (Int_t i=1; i<NSlices; i++)
    {
        
        Float_t xi = (det2-fSlice->GetAi(i))/fSlice->GetBi(i);
        Float_t xi_1 = (det2-fSlice->GetAi(i-1))/fSlice->GetBi(i-1);

        if ( (det1-xi)*(det1-xi_1) <= 0 )
        {
            centrality = step*(i+0.5);
            isOK = true;
            break;
        }
    }        
        
    if (!isOK)
        cout << "*** Warning *** GetCentrality: centrality is not found!" << endl;
    
    return centrality;
}


void CentralityGetter::GetGlauberB ()
{
    
    TCanvas *c1 = new TCanvas("c1", "canvas", 1500, 900);

    
    GlauberParGetter *gg = new GlauberParGetter;
    gg->SetSimTree ("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Glauber/examples/glau_pbpb_ntuple_signn_31.0_7.6AGeV_CM_30AGeV_LC.root");
    
    fCentrTree->GetEntry(0);
    
    Float_t step = fSlice->GetSlicesStep ();
    Int_t NSlices = fSlice->GetNSlices();
    Float_t norm = fSlice->GetDet1Max();  

    std::vector <Float_t> centrality;    
    for (UInt_t i=0; i<NSlices; i++)
        centrality.push_back( (i+0.5)*step );
    
    std::vector <Float_t> GlaubSigmaB_B;
    std::vector <Float_t> GlaubSigmaB;
    std::vector <Float_t> GlaubB;
    c1->Divide(2,2);
    
    
    TH1F *hTotal = new TH1F ("hTotal", "hTotal", 200, 0, 20);
    
    
    
    c1->cd(1);
    gPad->SetLogy();
    for (Int_t i=0; i<NSlices; i++)
    {
        
        Float_t MultMax = 1e5, MultMin;
        if (i>0) MultMax = fSlice->GetXi(i-1) * norm;
        MultMin = fSlice->GetXi(i) * norm;
        
        TH1F *h1 = new TH1F ("h1", "", 200, 0, 20);
        h1->SetName(Form("hB_%f_%f", MultMin, MultMax));
        h1->SetLineColor(i+1);
        
        fSlice->GetXi(i);
        gg->GetBHisto(MultMin, MultMax, h1, 1e5);
        
        h1->GetYaxis()->SetRangeUser(1, 500);

        if (i==0)  h1->Draw();
        else       h1->Draw("same");
        
        h1->Fit("gaus", "Q");
        TF1 *fFit = h1->GetFunction("gaus");
        fFit->SetLineColor(i+1);
        Float_t meanB = fFit->GetParameter(1); 
        Float_t sigmaB = fFit->GetParameter(2);        
        
        std::cout << "mean = " << meanB << "    sigma = " << sigmaB << "    area = " << h1->Integral() << std::endl;
        
        GlaubSigmaB_B.push_back (sigmaB/meanB);
        GlaubSigmaB.push_back (sigmaB);
        GlaubB.push_back (meanB);
        
        hTotal->Add(h1);
        
        gPad->Update();
    }
    hTotal->Draw("same");
    
    c1->cd(3);
    
    gPad->SetLogy();
    std::vector <Float_t> TempVecB, TempdB;
    for (UInt_t j=0; j<NSlices; j++){
        
        float temp_dB = fSlice->GetdB().at(j);
        float temp_B = fSlice->GetMeanB().at(j);
        float temp_sB = fSlice->GetSigmaB().at(j);
        float temp_dsB = fSlice->GetdSigmaB().at(j);
        
        TempVecB.push_back(temp_sB/temp_B);
        
        float ddd = TMath::Sqrt ( (temp_dsB/temp_B)*(temp_dsB/temp_B) + (temp_dsB/temp_B/temp_B*temp_dB)*(temp_dsB/temp_B/temp_B*temp_dB) );
        
//         cout << "d_sigma = " << (temp_dsB/temp_B) << endl;
//         cout << "d_b = " << (temp_dsB/temp_B/temp_B*temp_dB) << endl;
        
        TempdB.push_back( ddd );
    }

    TGraphErrors *grSigmaB = new TGraphErrors (NSlices, &(centrality[0]), &(TempVecB[0]), 0, &(TempdB[0]));
    grSigmaB->SetMarkerStyle(23);    
    grSigmaB->GetXaxis()->SetTitle( "Centrality" );
    grSigmaB->GetYaxis()->SetTitle( "sigma_{B}/<B>" );    
    grSigmaB->GetYaxis()->SetRangeUser(0.01, 0.5);
      
    grSigmaB->Draw("APL");

    TGraphErrors *grSigmaB1 = new TGraphErrors (NSlices, &(centrality[0]), &(GlaubSigmaB_B[0]), 0, 0);
    grSigmaB1->SetMarkerStyle(22);    
    grSigmaB1->SetMarkerColor(kRed);    
    grSigmaB1->SetLineColor(kRed);    

    grSigmaB1->Draw("PLsame");


    c1->cd(4);
    
    TGraphErrors *grB = new TGraphErrors (NSlices, &(centrality[0]), &(fSlice->GetMeanB()[0]), 0, &(fSlice->GetSigmaB()[0]));
    grB->SetMarkerStyle(23);    
    grB->GetXaxis()->SetTitle( "Centrality" );
    grB->GetYaxis()->SetTitle( "<B>" );    
    grB->GetYaxis()->SetRangeUser(0.0, 16.0);
      
    grB->Draw("APL");

    TGraphErrors *grB1 = new TGraphErrors (NSlices, &(centrality[0]), &(GlaubB[0]), 0, &(GlaubSigmaB[0]));
    grB1->SetMarkerStyle(22);    
    grB1->SetMarkerColor(kRed);    
    grB1->SetLineColor(kRed);    

    grB1->Draw("PLsame");
    


    
    

    
}






