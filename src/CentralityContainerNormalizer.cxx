
// check in the code:
// - TO DO
// - WARNING (to know for users)

#include <iostream>
#include <fstream>

 
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TProfile.h"
#include "TFile.h"
#include "TCut.h"
#include "TStyle.h"
#include "TEventList.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TPad.h"
#include "TImage.h"
#include "TLine.h"
#include "TRandom.h"

#include "CentralityContainerNormalizer.h"



using std::vector;
using std::cout;
using std::endl;
using std::flush;
using std::ifstream;

ClassImp(CentralityContainerNormalizer)


// -----   Default constructor   -------------------------------------------
CentralityContainerNormalizer::CentralityContainerNormalizer() 
  : TNamed(),
    fDet1Name (""),
    fDet2Name (""),
    fRunId (0),
    is1DAnalisys(false),
    isDet1Int(false),
    isDet2Int(false),
    fIsSimData (false),
    fNormalization (-1)
{

}

// -------------------------------------------------------------------------


// -----   Destructor   ----------------------------------------------------
CentralityContainerNormalizer::~CentralityContainerNormalizer()
{
    
}


void CentralityContainerNormalizer::LoadInputData (Int_t Det1Id, Int_t Det2Id)
{
    TCanvas *c2 = new TCanvas("c2", "canvas", 800, 800);
    
    fInFile = new TFile(fInFileName.Data());
    if (!fInFile)
    {
        cout << "Cannot open input file! Please check file name and path." << endl;
        exit (-1);
    }
    
    if ( fInFile->IsOpen() ) cout << "*** CentralityContainerNormalizer::LoadInputData ***  File opened successfully" << endl;
        
    fInTree = (TTree*) fInFile->Get("na61_data");    //TODO set Tree name as parameter
    fContainer = new CentralityEventContainer;
    fInTree->SetBranchAddress("CentralityEventContainer", &fContainer);   
    
    GetMaximum (Det1Id, Det2Id);
    RunByRunCorrection (Det1Id, Det2Id);
            
    if (Det2Id != -1) GetNormalization (Det1Id, Det2Id);
        else          GetNormalization (Det1Id);
    
    TRandom* random = new TRandom;
    random->SetSeed();
    Float_t rand1 = 0, rand2 = 0;

    Int_t nTotalEvents = fInTree->GetEntries();

    fInTree->GetEntry(0);
    Int_t RunId = fContainer->GetRunId();
    Float_t NormFactor = 1.;
    Int_t IterRunId = 0;
    Float_t NormFactor1 = fDet1NormVec.at(IterRunId);
    Float_t NormFactor2;
    if (Det2Id != -1) NormFactor2 = fDet2NormVec.at(IterRunId);

    std::cout << "Normalizing and applying run-by-run corrections..." << std::endl;
    
    TH1F *h1NotCorr = new TH1F ("h1NotCorr", "h1NotCorr", 300, 0, 300);
    TH1F *h1Corr    = new TH1F ("h1Corr", "h1Corr", 300, 0, 300);

    fNormTree = new TTree ("NormTree", "Norm Tree");
    fNormTree->SetDirectory(0);     
    fNormTree->Branch("det1", &det1, "det1/F");
    if (Det2Id != -1) fNormTree->Branch("det2", &det2, "det2/F");
    if (fIsSimData) fNormTree->Branch("B", &fB, "fB/F");
    fNormTree->Branch("RunId", &fRunId, "RunId/I");

    
    for (Int_t i=0; i<nTotalEvents; i++)
    {
        if ((i+1) % 50 == 0) cout << "Event # " << i+1 << "... \r" << flush;         
        fInTree->GetEntry(i);
        if (fIsSimData)  fB = fContainer->GetB();
        fRunId = fContainer->GetRunId();
        if (RunId != fContainer->GetRunId()/* || i == nTotalEvents-1 */)
        {
            IterRunId++;
            std::cout << "NormFactor for run # " << RunId <<  ", det1 = " << NormFactor1;
            if (Det2Id != -1)  std::cout << " , det2 = " << NormFactor2;
            std::cout << std::endl;
            NormFactor1 = fDet1NormVec.at(IterRunId);
            if (Det2Id != -1)  NormFactor2 = fDet2NormVec.at(IterRunId);
            RunId = fContainer->GetRunId();
        }
        det1 = NormFactor1*fContainer->GetDetectorWeight(Det1Id)/det1max;
        if (Det2Id != -1)  det2 = NormFactor2*fContainer->GetDetectorWeight(Det2Id)/det2max;
        
        if (isDet1Int)  { rand1 = random->Rndm()/* - 0.5*/;  /*cout << rand1 << endl;*/  det1 += rand1/det1max; }
        if (isDet2Int && Det2Id != -1)  { rand2 = random->Rndm()/* - 0.5*/;  det2 += rand2/det2max; }
        
        h1NotCorr->Fill(fContainer->GetDetectorWeight(Det1Id));
        h1Corr   ->Fill(NormFactor1*fContainer->GetDetectorWeight(Det1Id));
                
        fNormTree->Fill();        
    }
    
    
    
    fNormTree->Draw("det1 >> h1 (1000, 0, 1.1)");
    h1Corr->SetLineColor(kRed);

    h1NotCorr->Draw();
    h1Corr   ->Draw("same");
    h1Corr   ->SaveAs("hTPC_ref.root");
    
    c2->Print( CFdir + "QA/norm.pdf");
    c2->Print( CFdir + "QA/norm.root");
    
    c2->Delete();
    
    std::cout << "Normalization and correction done!" << std::endl;
    
    
}

void CentralityContainerNormalizer::RunByRunCorrection (Int_t Det1Id, Int_t Det2Id)
{
    std::cout << "Calculating run-by-run corrections..." << std::endl;

    Int_t nTotalEvents = fInTree->GetEntries();

    fInTree->GetEntry(0);
    Int_t RunId = fContainer->GetRunId();
//     fRunIdVec.push_back(RunId);
    std::cout << "Run # " << RunId << " " << std::endl;

    Float_t Det1Temp = 0, Det2Temp = 0;
    Float_t MeanDet1 = 0, MeanDet2 = 0;
    Float_t sig1, sig2;
    Int_t nEventsInRun1 = 0, nEventsInRun2 = 0;
    
    TProfile *profile1 = new TProfile ("profile1", "profile1", 100, 22914, 23013, 0, det1max );
    TProfile *profile2 = new TProfile ("profile2", "profile2", 100, 22914, 23013, 0, det2max );
    
    
    for (Int_t i=0; i<nTotalEvents; i++)
    {
        fInTree->GetEntry(i);
        
        sig1 = fContainer->GetDetectorWeight(Det1Id);
        if (Det2Id != -1)  sig2 = fContainer->GetDetectorWeight(Det2Id);
        
        if (RunId != fContainer->GetRunId() || i == nTotalEvents-1)
        {
            std::cout << "Run # " << RunId << " " << std::endl;

            fRunIdVec.push_back(RunId);
            nEventsInRun1Vec.push_back(nEventsInRun1);
            nEventsInRun2Vec.push_back(nEventsInRun2);
            fDet1NormVec.push_back( Det1Temp/nEventsInRun1 );
            if (Det2Id != -1)  fDet2NormVec.push_back( Det2Temp/nEventsInRun2 );
            
            RunId = fContainer->GetRunId();
            nEventsInRun1 = nEventsInRun2 = 0;
            Det1Temp = Det2Temp = 0.0;

            continue;
        }
        if ( /*sig1 < 0.7*det1max &&*/ sig1 > 0.3*det1max ){
            Det1Temp += sig1;
            MeanDet1 += sig1;
            profile1->Fill(RunId, sig1, 1);
            nEventsInRun1++;
        }
        if (Det2Id != -1 /*&& sig2 < 0.7*det2max*/ && sig2 > 0.3*det2max) { 
            Det2Temp += sig2;
            MeanDet2 += sig2;
            profile2->Fill(RunId, sig2, 1);
            nEventsInRun2++;
        }
        
    }
        
    Int_t nTotalEvents1 = 0, nTotalEvents2 = 0;
    for (UInt_t i=0; i<nEventsInRun1Vec.size(); i++)
        nTotalEvents1 += nEventsInRun1Vec.at(i);
        
    MeanDet1 /= nTotalEvents1;
    if (Det2Id != -1)  {
        for (UInt_t i=0; i<nEventsInRun2Vec.size(); i++)
            nTotalEvents2 += nEventsInRun2Vec.at(i);
        MeanDet2 /= nTotalEvents2;
    }
    
    for (UInt_t i=0; i<fDet1NormVec.size(); i++)
    {
        fDet1NormVec.at(i) = MeanDet1/fDet1NormVec.at(i);
        if (Det2Id != -1)  fDet2NormVec.at(i) = MeanDet2/fDet2NormVec.at(i);
    }
    std::cout << "Corrections calculated!" << std::endl;
/*
    TCanvas *c3 = new TCanvas ("c3", "c3", 800, 600);
    profile1->Scale(1.0/MeanDet1);
    profile2->Scale(1.0/MeanDet2);
    profile2->SetLineColor(kRed);
    profile1->SetMinimum(0.95);
    profile1->SetMaximum(1.05);
    profile1->Draw();
    profile2->Draw("same");

    gPad->Update();
    
    TString fn = "CutsMidHigh";
    
    c3->SaveAs( fn + ".pdf");
    c3->SaveAs( fn + ".root");
    */
    
    
    
    
}



void CentralityContainerNormalizer::GetMaximum (Int_t Det1Id, Int_t Det2Id)
{
    det1max = 0;
    det2max = 0;
    
    Int_t nTotalEvents = fInTree->GetEntries();
    Float_t sig1, sig2;

    for (Int_t i=0; i<nTotalEvents; i++)
    {
        fInTree->GetEntry(i);
        if ((i+1) % 50 == 0) cout << "Event # " << i+1 << "... \r" << flush; 
        
        sig1 = fContainer->GetDetectorWeight(Det1Id);
        if (sig1 > det1max)    det1max = sig1;

        if (Det2Id != -1){
            sig2 = fContainer->GetDetectorWeight(Det2Id);
            if (sig2 > det2max)    det2max = sig2;        
        }
    }
}



void CentralityContainerNormalizer::GetNormalization (Int_t Det1Id, Int_t Det2Id)
{
    std::cout << "Calculating normalization value..." << std::endl;

    
    Int_t n1 = 100, n2 = 100;
    if (isDet1Int)  n1 = det1max;
    if (isDet2Int)  n2 = det2max;
    
    
    TString DrawPar1 = Form ( "CentralityEventContainer.GetDetectorWeight(%d) >> h1(%d, 0., %f)", Det1Id, n1, det1max) ;
    fInTree->Draw( DrawPar1.Data() );        
    TH1F *h1 = (TH1F*)gPad->GetPrimitive("h1");

    TString DrawPar2 = Form ( "CentralityEventContainer.GetDetectorWeight(%d) >> h2(%d, 0., %f)", Det2Id, n2, det2max) ;
    fInTree->Draw( DrawPar2.Data() );        
    TH1F *h2 = (TH1F*)gPad->GetPrimitive("h2");

    
    for (Int_t i=n1/2; i<n1; i++)
    {
//         std::cout << "h1 = " << h1->GetBinContent(i+1) << std::endl; 
        if (h1->GetBinContent(i+1) < 2)
        {
            det1max *= i/float(n1);
            break;
        }
    }
    
    for (Int_t i=n2/2; i<n2; i++)
    {
//         std::cout << "h1 = " << h1->GetBinContent(i+1) << std::endl; 
        if (h2->GetBinContent(i+1) < 2)
        {
            det2max *= i/float(n2);
            break;
        }
    }    
    std::cout << "   det1max = " << det1max << std::endl; 
    std::cout << "   det2max = " << det2max << std::endl; 
    
    std::cout << "Normalization calculated!" << std::endl;
    
}

void CentralityContainerNormalizer::GetNormalization (Int_t Det1Id)
{
    std::cout << "Calculating normalization value..." << std::endl;
    
    Int_t n1 = 100;
    if (isDet1Int)  n1 = det1max;
    
    TString DrawPar1 = Form ( "CentralityEventContainer.GetDetectorWeight(%d) >> h1(%d, 0., %f)", Det1Id, n1, det1max) ;
    fInTree->Draw( DrawPar1.Data() );        
    TH1F *h1 = (TH1F*)gPad->GetPrimitive("h1");
    
    for (Int_t i=n1/2; i<n1; i++)
    {
//         std::cout << "h1 = " << h1->GetBinContent(i+1) << std::endl; 
        if (h1->GetBinContent(i+1) < 1)
        {
            det1max *= i/(float)n1;
            break;
        }
    }
    
//     c2->Print( CFdir + "QA/norm.pdf");
//     c2->Print( CFdir + "QA/norm.root");
    
    std::cout << "   det1max = " << det1max << std::endl; 
    std::cout << "Normalization calculated!" << std::endl;
}
