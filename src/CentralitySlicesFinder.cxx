

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

#include "CentralitySlicesFinder.h"



using std::vector;
using std::cout;
using std::endl;
using std::ifstream;

ClassImp(CentralitySlicesFinder)


// -----   Default constructor   -------------------------------------------
CentralitySlicesFinder::CentralitySlicesFinder() 
  : TNamed(),
    fDet1Name (""),
    fDet2Name (""),
    fCuts(""),
    fBaseCuts(""),
    fCentralityMax(100),
    fPrecision(0.01),
    fDirectionCentralEvents (1),
    fRunId (0),
    is1DAnalisys(false),
    isDet1Int(false),
    isDet2Int(false),
    fSliceStep (5.),
    nIntervals (20),
    fIsSimData (false),
    fNormalization (-1),
    fSlice (new CentralitySlice),
    fCentrTree (new TTree ("CentrTree", "Tree with slices parameters")),   
    fFitFunction (new TF1)
{
    fCentrTree->Branch("CentralitySlice", "CentralitySlice", &fSlice);

}

// -------------------------------------------------------------------------


// -----   Destructor   ----------------------------------------------------
CentralitySlicesFinder::~CentralitySlicesFinder()
{
    
    if ( fCentrTree ) 
	fCentrTree->Delete();      
    
    if ( fSlicesFile ) {
        fSlicesFile->Delete();
    }
}

void CentralitySlicesFinder::WriteOutputData ()
{

    if ( fSlicesFile->IsOpen() ) {
        fSlicesFile->Close();
        cout << "File is closed" << endl;
    }    
}


void CentralitySlicesFinder::LoadInputData (Int_t Det1Id, Int_t Det2Id)
{
    TCanvas *c2 = new TCanvas("c2", "canvas", 800, 800);
    
    fInFile = new TFile(fInFileName.Data());
    if (!fInFile)
    {
        cout << "Cannot open input file! Please check file name and path." << endl;
        exit (-1);
    }
    
    if ( fInFile->IsOpen() ) cout << "*** CentralitySlicesFinder::LoadInputData ***  File opened successfully" << endl;
        
    fInTree = (TTree*) fInFile->Get("cbm_data");
    
    fContainer = new CentralityEventContainer;
    fInTree->SetBranchAddress("CentralityEventContainer", &fContainer);   
    
    fNormTree = new TTree ("NormTree", "Norm Tree");
    fNormTree->SetDirectory(0); 
    
    fNormTree->Branch("det1", &det1, "det1/F");
    if (Det2Id != -1) fNormTree->Branch("det2", &det2, "det2/F");
    if (fIsSimData) fNormTree->Branch("B", &fB, "fB/F");
    fNormTree->Branch("RunId", &fRunId, "RunId/I");
    
    if (Det2Id != -1) GetNormalization (Det1Id, Det2Id);
        else          GetNormalization (Det1Id);
    
    TRandom* random = new TRandom;
    
    Float_t rand1 = 0, rand2 = 0;

    Int_t nTotalEvents = fInTree->GetEntries();
    
    for (Int_t i=0; i<nTotalEvents; i++)
    {
        fInTree->GetEntry(i);
                
        det1 = (fContainer->GetDetectorWeight(Det1Id))/det1max;
        if (Det2Id != -1)  det2 = fContainer->GetDetectorWeight(Det2Id)/det2max;
        if (fIsSimData)  fB = fContainer->GetB();
        fRunId = fContainer->GetRunId();
        
//         std::cout << "b = " << fContainer->GetB() << std::endl;
        
        if (isDet1Int)  { rand1 = random->Rndm()/* - 0.5*/;  /*cout << rand1 << endl;*/  det1 += rand1/det1max; }
        if (isDet2Int && Det2Id != -1)  { rand2 = random->Rndm() - 0.5;  det2 += rand2/det2max; }
        
        fNormTree->Fill();        
    }
    
    fNormTree->Draw("det1 >> h1 (1000, 0, 1.1)");
    c2->Print( CFdir + "QA/norm.pdf");
    c2->Print( CFdir + "QA/norm.root");
}


void CentralitySlicesFinder::GetNormalization (Int_t Det1Id, Int_t Det2Id)
{
    
    det1max = 0;
    det2max = 0;
    
    Int_t nTotalEvents = fInTree->GetEntries();


    for (Int_t i=0; i<nTotalEvents; i++)
    {
        fInTree->GetEntry(i);

        Float_t sig1 = fContainer->GetDetectorWeight(Det1Id);
        Float_t sig2 = fContainer->GetDetectorWeight(Det2Id);
        
        if (sig1 > det1max)    det1max = sig1;
        if (sig2 > det2max)    det2max = sig2;        
    }
    
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
    std::cout << "det1max = " << det1max << std::endl; 
    std::cout << "det2max = " << det2max << std::endl; 
}

void CentralitySlicesFinder::GetNormalization (Int_t Det1Id)
{
    det1max = 0;
    
    Int_t nTotalEvents = fInTree->GetEntries();


    for (Int_t i=0; i<nTotalEvents; i++)
    {
        fInTree->GetEntry(i);

        Float_t sig1 = fContainer->GetDetectorWeight(Det1Id);
        if (sig1 > det1max)    det1max = sig1;
    }
    
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
    
    std::cout << "det1max = " << det1max << std::endl; 
    
}

void CentralitySlicesFinder::RunSliceFinder (Int_t RunId)
{
    fCentrTree = new TTree ("CentrTree", "Tree with slices parameters");
    fCentrTree->Branch("CentralitySlice", "CentralitySlice", &fSlice);

    FindCentralitySlices( RunId );
    
    TString filename = Form(CFdir + "root_files/Slices_%s_%s_%d.root", fDet1Name.Data(), fDet2Name.Data(), fRunId );
    
    fSlicesFile = new TFile (filename, "RECREATE");
    std::cout << "Output file name is " << filename << std::endl;
    
    fCentrTree->Write(); 
}

void CentralitySlicesFinder::Fit2DCorrelation ()
{

//     TH2F *h2D = new TH2F ("h2D", "det12", 100, 0., 1., 100, 0., 1.);
    TProfile *p2D;
//     TString DrawPar = distrName + " >> h1(100, 0., 1., 100, 0., 1.)";
    fNormTree->Draw( "det2 : det1  >> h1(120, 0., 1.2, 120, 0., 1.2)", fCuts, "colz");
    
    TH2F *h2D = (TH2F*)gPad->GetPrimitive("h1");
    h2D->SetMaximum(40);
    h2D->SetTitle ("2D correlation fit");
    h2D->GetXaxis()->SetTitle( Form ( "%s/%s_{max}", fDet1Name.Data(), fDet1Name.Data() ) );
    h2D->GetYaxis()->SetTitle( Form ( "%s/%s_{max}", fDet2Name.Data(), fDet2Name.Data() ) );   

//     DrawPar = distrName + " >> prof (100, 0., 1., 100, 0., 1.)";
    
    fNormTree->Draw("det2 : det1 >> prof (100, 0., 1., 100, 0., 1.)", fCuts, "sameprof");
    p2D = (TProfile*)gPad->GetPrimitive("prof");
    
    p2D->SetLabelSize(10);
    p2D->SetBarWidth(10);
    p2D->SetLineColor(1);    
    p2D->SetLineWidth(3);
    p2D->Draw("same");
    gPad->Update();

//     c1->Print( CFdir + "QA/profile.pdf");    
//     c1->Print( CFdir + "QA/profile.root");    



    p2D->Fit("pol3", "", "", 0.0, 1.0);
    gPad->Update();

//     c1->Print( CFdir + "QA/profile_fit.pdf");    
//     c1->Print( CFdir + "QA/profile_fit.root");    

    fFitFunction = (TF1*) p2D->GetFunction("pol3");
    
//     TLegend* leg1 = new TLegend(0.7,0.75,0.89,0.89);
//     leg1->AddEntry(fFitFunction ,"start fit","l");    
//     leg1->Draw("same");
    
    gPad->Update();
    
}

void CentralitySlicesFinder::FitCorrection ( int n_points )
{
//     TCanvas *c1 = new TCanvas("c1", "canvas", 1000, 600);

    std::vector <double> X4fit;
    std::vector <double> Y4fit;
    
    std::vector <double> ErrorX;
    std::vector <double> ErrorY;
    
    double par[n_par];
    for (int i=0; i<n_par; i++)
        par[i] = fFitFunction->GetParameter(i);
    
    double width = 0.02;
    double kb[2];
    fNormTree->Draw( "det2 : det1  >> h1(120, 0., 1.2, 120, 0., 1.2)", fCuts, "colz");

    for (int i=0; i<n_points; i++)
    {
        double x = 1.0*(i+0.5)/n_points;
        double y = polN (par, x, n_par);
      
        find_norm (par, x, kb, n_par); 
        
        double tempk = -1/kb[0];                     // fit line 
        double tempb = kb[1] + x*(kb[0] - tempk);    
        
        
        double deltabNorm = kb[0]/TMath::Abs(kb[0]) * (width/2) / TMath::Cos (TMath::ATan (kb[0]));
            
        TCut isRight, isLeft;
            
//         cout << "kb[0] = " << kb[0] <<  "      kb[1] = " << kb[1] << "       deltabNorm = " << deltabNorm << endl;
    
        isRight = Form ( "det1 > abs(det2 - %f)/ abs(%f)",    kb[1] + deltabNorm,  kb[0]);
        isLeft =  Form ( "det1 < abs(det2 - %f)/ abs(%f)",    kb[1] - deltabNorm,  kb[0]);

        TCut isInside = isRight && isLeft;
        
        fNormTree->Draw(">>elist", fCuts && isInside);
        TEventList *insideList = (TEventList*)gDirectory->Get("elist");
        
        int nInside = insideList->GetN();
        int nn = 0;
        double distance;
        
        double mean = 0;
        double sigmaX = 0;
        double sigmaY = 0;
        
        for (int k=0; k<nInside; k++){
            fNormTree->GetEntry(insideList->GetEntry(k));
            
            double dy = det1*tempk - det2 + tempb;    
            double dx = dy*tempk;
                        
            distance = dy / TMath::Sqrt( 1 + tempk*tempk);  // distance from point to fit function
            if (TMath::Abs(distance) < 0.5){
                mean += distance;     
                
                sigmaY += dy*dy;
                sigmaX += dx*dx;
                
                nn++;
            }
        }
        
//         cout << "nn = " << nn << "       nInside = " << nInside << endl;
//         cout << "mean = " << mean << "       tempk = " << tempk << endl;

        if (nn > 1){
            mean /= nn;
            sigmaX = TMath::Sqrt(sigmaX/(nn-1));
            sigmaY = TMath::Sqrt(sigmaY/(nn-1));
        }
        else {
            cout << "ERROR!" << endl;
            continue;
        }
            
        
//         cout << "x = " << x << "       y = " << y << endl;
        X4fit.push_back(x + mean/TMath::Sqrt( 1 + tempk*tempk)*tempk) ;
        Y4fit.push_back(y - mean/TMath::Sqrt( 1 + tempk*tempk)) ;
        ErrorX.push_back (sigmaX);
        ErrorY.push_back (sigmaY);
    }

//     for (int i=0; i<X4fit.size(); i++ )
//         cout << "x = " << X4fit.at(i) << "       y = " << Y4fit.at(i) << endl;
    
    TGraph *gr = new TGraphErrors( X4fit.size(), &(X4fit[0]), &(Y4fit[0]), &(ErrorX[0]), &(ErrorY[0]));
    gr->GetXaxis()->SetLimits(0.0, 1.0);
    gr->GetYaxis()->SetRangeUser(0.0, 1.0);
    gr->SetMarkerStyle(22);  
    gr->SetMarkerSize(1);
    gr->SetMarkerColor(1);
    gr->Draw("Psame");
    
//     c1->Print( CFdir + "QA/corr.pdf");    
//     c1->Print( CFdir + "QA/corr.root");    
    
    
    
    
    
    gr->Fit("pol5", "", "", 0.00, 1.0);
    gr->GetFunction("pol5")->SetLineWidth(2);
    gr->GetFunction("pol5")->SetLineColor(1);
    fFitFunction = (TF1*) gr->GetFunction("pol5");
    
/*    TLegend* leg1 = new TLegend(0.6,0.8,0.75,0.89);
    leg1->AddEntry(fFitFunction ,"corrected fit","l");    
    leg1->Draw("same");  */  

//     c1->Print( CFdir + "QA/fit.pdf");    
//     c1->Print( CFdir + "QA/fit.root");    
    
    
}


void CentralitySlicesFinder::FindCentralitySlices (Int_t RunId )
{
    gStyle->SetOptStat(0000);
    fSlice->ClearData ();
    
    nIntervals = Int_t (100/fSliceStep);
    
    fCuts = fBaseCuts && Form( "RunId == %d", RunId );
    fRunId = RunId;
    std::cout << " fRunId = " << fRunId << endl;     
    
    TString distrName;
    if (!is1DAnalisys)  distrName = Form("%s_%s_run_%d_centrality_QA", fDet1Name.Data(), fDet2Name.Data(), fRunId );
    else                distrName = Form("%s_run_%d_centrality_QA", fDet1Name.Data(), fRunId );
        
    TCanvas *c1 = new TCanvas("c1", "canvas", 1500, 900);

    if (!is1DAnalisys) Fit2DCorrelation ();    
    if (!is1DAnalisys) FitCorrection (50);   

    if (is1DAnalisys)
    {
        fNormTree->Draw( "det1 >> h1(120, 0., 1.1)", fCuts);        
        c1->SetLogy();
    }
    
    fRunId = RunId;
    FindSlices ();    

    /*if (!is1DAnalisys)  */FindMeanSignalsInSlice();
    
    c1->Print( CFdir + "QA/" + distrName + ".pdf");    
    c1->Print( CFdir + "QA/" + distrName + ".root");    
    
    
//     fSlicesFile = new TFile ( Form(CFdir + "root_files/Slices_%s_%s_%d.root", fDet1Name.Data(), fDet2Name.Data(), RunId ) , "RECREATE");

    fSlice->SetDetMax (det1max, det2max);
    fSlice->SetRunId (RunId);
    fSlice->SetFitFunction (fFitFunction);
    fSlice->SetSlicesStep (fSliceStep);
    fSlice->SetDirectionCentralEvents (fDirectionCentralEvents);
    
    fCentrTree->Fill();  
    
//     QA ();
//     fSlicesFile->Write();    
    
}

void CentralitySlicesFinder::FindSlices ( )
{
    
    double finalPar[6] = {0, 0, 0, 0, 0, 0};

    if (!is1DAnalisys)
        for (int i=0; i<nFinalPar; i++)
            finalPar[i] = fFitFunction->GetParameter(i);
    
    fNormTree->Draw(">>elist", Form("RunId == %d", fRunId) );
    TEventList* elist = (TEventList*)gDirectory->Get("elist");
    int nTotalEvents = fNormalization > 0 ? fNormalization : elist->GetN();
   
    std::cout << " fRunId = " << fRunId << endl;     
    std::cout << " nTotalEvents = " << nTotalEvents << endl; 

    fNormTree->Draw(">>elist1", fCuts );
    TEventList* elist1 = (TEventList*)gDirectory->Get("elist1");
    int nTotalEvents1 = elist1->GetN();

    std::cout << "Events with cuts = " << nTotalEvents1 << "   fraction of total = " << nTotalEvents1*100.0/nTotalEvents << " %"  << endl; 

    int nInside = (fSliceStep/100) * nTotalEvents;
    int nIntervalsCut = nIntervals*fCentralityMax/100; //need in case cuts usage
    if (fCentralityMax == 100) nIntervalsCut--;
    
    double x = (!fDirectionCentralEvents) ?  0.0 :  1.0;
    
    for (int i=0; i<nIntervalsCut; i++)
    {
        nInside = (i+1)*1.0/nIntervals * nTotalEvents;
        int nInsideCurrent = 0;
        double step = 0.5;//(!fDirectionCentralEvents) ?  (1-abs(x))/2.0 :  abs(x)/2.0;  //0.5;  
        double kb[2] = {0, 0};

//         std::cout << " i = " << i << endl; 
//         std::cout << " nInside = " << nInside << endl; 
        
        TCut RightOfNorm = "", LeftOfNorm = "";

        while ( TMath::Abs(nInside - nInsideCurrent) > nInside*fPrecision/(i+1) && step > 1e-11 )
        {
            if (!is1DAnalisys){
                find_norm(finalPar, x, kb, nFinalPar);

                if (fDirectionCentralEvents == 0)
                    LeftOfNorm = Form ( "det1 < (det2 - %f)/%f",  kb[1],   kb[0] );

                if (fDirectionCentralEvents == 1)
                    RightOfNorm  = Form ( "det1 >= (det2 - %f)/%f",  kb[1],   kb[0] );

            }
            else {
                if (fDirectionCentralEvents == 0)
                    LeftOfNorm = Form ( "det1 < %12.11f",  x );

                if (fDirectionCentralEvents == 1)
                    RightOfNorm  = Form ( "det1 >= %12.11f",  x );
          
                }
                
//             std::cout << " x = " << x << endl; 
            TCut isInside = RightOfNorm && LeftOfNorm;        
            fNormTree->Draw(">>elist", fCuts && isInside );
            elist = (TEventList*)gDirectory->Get("elist");
            nInsideCurrent =  elist->GetN();  // number of events to pass cuts
//             std::cout << " nInsideCurrent = " << nInsideCurrent << endl; 
            
            if (TMath::Abs(nInside - nInsideCurrent) > nInside*fPrecision/(i+1)){
                if (fDirectionCentralEvents == 0)               x += nInsideCurrent > nInside ? -step : step;
                if (fDirectionCentralEvents == 1)               x += nInsideCurrent > nInside ? step : -step;
            }
            
            if (x<0 && !is1DAnalisys) x = 0.;
            if (x>1 && !is1DAnalisys) x = 1.;
            
            step /= 2;
     
        }
        std::cout << " % inside = " << 100.*nInsideCurrent/nTotalEvents << endl; 
        
        fSlice->AddSlice(kb[1], kb[0], x);
                 
//         std::cout << " X  = " << fSlice->GetXi(i) << endl; 
//         std::cout << " B  = " << fSlice->GetBi(i) << endl; 
//         std::cout << " A = " << fSlice->GetAi(i) << endl; 

        
        double xNormInterval = 0.15 * TMath::Abs ( TMath::Cos (TMath::ATan (kb[0]) ) );
        if (is1DAnalisys)  xNormInterval = 0.15 / 1e5 ;
                                                            
        if (!is1DAnalisys){
            TF1 *norm = new TF1("norm","[0]*x + [1]", x-xNormInterval, x+xNormInterval); 
            norm->SetParameters (kb[0], kb[1]);
            norm->Draw("same");  
        }
        else{
            gPad->Update();
            
            Double_t xr = gPad->GetX2()-gPad->GetX1();
            double x1 = (x-gPad->GetX1())/ xr;
            TLine l2;
            l2.SetLineColor(kRed);
            l2.DrawLineNDC(x1,0.15,x1,0.9);             
        } 
    }
    
}


void CentralitySlicesFinder::FindMeanSignalsInSlice ( void )
{
    
    int nIntervalsCut = nIntervals*fCentralityMax/100; //need in case cuts usage
    int nIntervalsCut1 = nIntervalsCut;
    if (fCentralityMax != 100) nIntervalsCut++;
    
    // TODO check ranges for loop for 100% centrality (last cut)
   
    for (int i=0; i<nIntervalsCut1; i++)
    {
        TCut RightOfNorm, LeftOfNorm;
        
        
        if (i>0 && (fDirectionCentralEvents == 0)){
            if (is1DAnalisys)             
                RightOfNorm = Form ( "det1 >= %f",  fSlice->GetXi(i-1) );
            else
                RightOfNorm = Form ( "det1 >= (det2 - %f)/%f", fSlice->GetAi(i-1), fSlice->GetBi(i-1) );
        }
        if (i>0 && (fDirectionCentralEvents == 1)){
            if (is1DAnalisys)             
                LeftOfNorm = Form ( "det1 < %f", fSlice->GetXi(i-1) );
            else
                LeftOfNorm = Form ( "det1 < (det2 - %f)/%f", fSlice->GetAi(i-1), fSlice->GetBi(i-1) );
        }

        if (i<nIntervalsCut-1 && (fDirectionCentralEvents == 0)){
            if (is1DAnalisys)             
                LeftOfNorm = Form ( "det1 < %f", fSlice->GetXi(i) );
            else
                LeftOfNorm = Form ( "det1 < (det2 - %f)/%f", fSlice->GetAi(i), fSlice->GetBi(i) );
        }
        
        if (i<nIntervalsCut-1 && (fDirectionCentralEvents == 1)){
            if (is1DAnalisys)             
                RightOfNorm = Form ( "det1 >= %f", fSlice->GetXi(i) );
            else
                RightOfNorm = Form ( "det1 >= (det2 - %f)/%f", fSlice->GetAi(i), fSlice->GetBi(i) );
        }                
        
//         cout << "x_i = " << fSlice->GetXi(i) << "     x_i-1 = " << fSlice->GetXi(i-1) << endl;
         
        double tempk, tempb;
        
        if (i<nIntervalsCut-1){
            tempk = -1/fSlice->GetBi(i);
            tempb = fSlice->GetAi(i) + fSlice->GetXi(i)*(fSlice->GetBi(i) - tempk);  
        }
        if (i==nIntervalsCut-1){
            tempk = -1/fSlice->GetBi(i-1);
            tempb = fSlice->GetAi(i-1) + fSlice->GetXi(i-1)*(fSlice->GetBi(i-1) - tempk);  
        }
//         float tempD = TMath::Sqrt( 1 + tempk*tempk);
//         TCut Dist1 = Form ( " abs((%s*%f - %s + %f)/%f) < %f", par1.Data(), tempk, par2.Data(), tempb, tempD, cutWidth );
        
        TCut isInside = RightOfNorm && LeftOfNorm;        
                
        fNormTree->Draw(">>elist", fCuts && isInside);
        TEventList *insideList = (TEventList*)gDirectory->Get("elist");
        
        int nInside = insideList->GetN();
        
        std::vector <double> X(nInside);
        std::vector <double> Y(nInside);    
        std::vector <double> XY(nInside);                   
        std::vector <double> B(nInside);

        double meanB = 0, meanParX = 0, meanParY = 0, meanParXY = 0;
        int nn = 0;
        
        for (int k=0; k<nInside; k++)
        {
            fNormTree->GetEntry(insideList->GetEntry(k));

            X.at(k) = det1;
            Y.at(k) = det2;
            B.at(k) = fB; 

//         std::cout << "fill b = " << fB << std::endl;

            
            XY.at(k) = ( X.at(k)*tempk - Y.at(k) + tempb) / TMath::Sqrt( 1 + tempk*tempk);  

            meanB     += B.at(k);
            meanParX  += X.at(k);
            meanParY  += Y.at(k);
            meanParXY += XY.at(k); 
            nn++;
            
        }
        cout << "nInside = " << nInside << endl;
       
        if (nn > 0){
            meanB     /= nn;
            meanParX  /= nn;
            meanParY  /= nn;
            meanParXY /= nn;        
        }
        else continue;
        
        double sigmaB = 0, sigmaParX = 0, sigmaParY = 0, sigmaParXY = 0, meanParXY3 = 0;
        
        for (int k=0; k<nInside; k++)
        {
            sigmaB     += TMath::Power(B.at(k) - meanB, 2);
            sigmaParX  += TMath::Power(X.at(k) - meanParX, 2);
            sigmaParY  += TMath::Power(Y.at(k) - meanParY, 2);
            sigmaParXY += TMath::Power(XY.at(k) - meanParXY, 2);
            meanParXY3 += TMath::Power(XY.at(k) - meanParXY, 3);
            
        }

        sigmaB = TMath::Power(sigmaB/nn, 0.5);
        sigmaParX = TMath::Power(sigmaParX/nn, 0.5);
        sigmaParY = TMath::Power(sigmaParY/nn, 0.5);
        sigmaParXY = TMath::Power(sigmaParXY/nn, 0.5);
        meanParXY3 = meanParXY3/nn/TMath::Power(sigmaParXY, 3);

        
        fSlice->AddXPar (meanParX, sigmaParX);
        fSlice->AddYPar (meanParY, sigmaParY);
        fSlice->AddXYPar (meanParXY, sigmaParXY);
        fSlice->AddXY3 (meanParXY3);
        fSlice->AddB (meanB, sigmaB);
    }
}



void CentralitySlicesFinder::find_norm (double par[], double x, double kb[], int N)
{
    double dx = 0.001;
    
    double y1, y2;
    double cx, cy, c;   
    
    y1 = polN(par, x - dx, N);
    y2 = polN(par, x + dx, N);
    
    // cx*x + cy*y + c == 0
    
    cx = 1/(y2 - y1);
    cy = 0.5/dx;
        
    c = -cx*x - cy*polN(par, x, N);
    
    kb[0] = -cx/cy;
    kb[1] = -c/cy;

}


void CentralitySlicesFinder::QA ()
{
    
    TCanvas *c1 = new TCanvas("c1", "canvas", 1500, 900);
    c1->Divide (2, 2);

    std::vector <Float_t> centrality;
    Float_t step = fSlice->GetSlicesStep ();
    UInt_t NSlices = fSlice->GetNSlices();
    
    for (UInt_t i=0; i<NSlices; i++)
        centrality.push_back( (i+0.5)*step );
        
    c1->cd(1);
    
    fNormTree->Draw( "det2 : det1  >> h1(120, 0., 1.2, 120, 0., 1.2)", fCuts, "colz");
//     TH2F *h1 = (TH2F*)gPad->GetPrimitive("h1");
    fFitFunction->Draw("same");

    for (UInt_t i=0; i<NSlices; i++)
    {
            Float_t x = fSlice->GetXi(i);
            Float_t xNormInterval = 0.15 * TMath::Abs ( TMath::Cos (TMath::ATan (fSlice->GetBi(i)) ) );
            TF1 *norm = new TF1("norm","[0]*x + [1]", x-xNormInterval, x+xNormInterval); 
            norm->SetParameters (fSlice->GetBi(i), fSlice->GetAi(i));
            norm->Draw("same");  
    }
    
    c1->cd(2);
    
    TGraphErrors *grX = new TGraphErrors(NSlices, &(centrality[0]), &(fSlice->GetMeanX()[0]), 0, &(fSlice->GetSigmaX()[0]));
    grX->SetTitle ("2D fit: X, Y vs centrality");
    grX->SetMarkerStyle(23);    
    grX->GetXaxis()->SetTitle( "Centrality, %" );
    grX->GetYaxis()->SetTitle("X, Y");    

    grX->GetXaxis()->SetLimits(0.0, NSlices*step);
    grX->GetYaxis()->SetRangeUser(0.0, 1);
    grX->Draw("APL");

    TGraphErrors *grY = new TGraphErrors(NSlices, &(centrality[0]), &(fSlice->GetMeanY()[0]), 0, &(fSlice->GetSigmaY()[0]));
    grY->SetMarkerStyle(22);    
    grY->SetMarkerColor(3);    
    grY->Draw("PLsame");


    c1->cd(3);
    
    TGraphErrors *grXY = new TGraphErrors(NSlices, &(centrality[0]), &(fSlice->GetMeanXY()[0]), 0, &(fSlice->GetSigmaXY()[0]));
    grXY->SetTitle ("2D fit: XY vs centrality");
    grXY->SetMarkerStyle(23);    
    grXY->GetXaxis()->SetTitle( "Centrality, %" );
    grXY->GetYaxis()->SetTitle("XY");    
    grXY->GetXaxis()->SetLimits(0.0, NSlices*step);
    grXY->GetYaxis()->SetRangeUser(-0.2, 0.2);
    grXY->Draw("APL");
    
    c1->cd(4);
    
    TGraphErrors *grXY3 = new TGraphErrors(NSlices, &(centrality[0]), &(fSlice->GetMeanXY3()[0]));
    grXY3->SetTitle ("2D fit: XY3 vs centrality");
    grXY3->SetMarkerStyle(23);    
    grXY3->GetXaxis()->SetTitle( "Centrality, %" );
    grXY3->GetYaxis()->SetTitle("XY3");    
    grXY3->GetXaxis()->SetLimits(0.0, NSlices*step);
    grXY3->GetYaxis()->SetRangeUser(-5, 5);
    grXY3->Draw("APL");
    
    
    c1->Print( CFdir + "QA/testQA.pdf");    
    
//     TGraph g(a.size(), &(a[0]), &(b[0]));
    
    
}

