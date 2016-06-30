

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
#include "CentralityContainerNormalizer.h"



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
    fCentralityMax(100),
    fPrecision(0.01),
    fDirectionCentralEvents (1),
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
        cout << "Output file is closed" << endl;
    }    
}


void CentralitySlicesFinder::LoadInputData (Int_t Det1Id, Int_t Det2Id)
{

    fDet1Id = Det1Id;
    fDet2Id = Det2Id;
    
    CentralityContainerNormalizer *fCCNorm = new CentralityContainerNormalizer;

    fCCNorm->SetDet1Name (fDet1Name) ;  
    fCCNorm->SetDet2Name (fDet2Name) ;  
    fCCNorm->SetInFileName (fInFileName) ;     
    fCCNorm->Do1DAnalisys (is1DAnalisys) ;
    fCCNorm->IsSimData (fIsSimData) ;
    fCCNorm->Det1IsInt (isDet1Int) ;
    fCCNorm->Det2IsInt (isDet2Int) ;
    fCCNorm->SetDir (CFdir) ;
    fCCNorm->LoadInputData (Det1Id, Det2Id);
    
    fNormTree = fCCNorm->GetNormTree ();
    fNormTree->SetBranchAddress("det1", &det1);
    if (Det2Id != -1) fNormTree->SetBranchAddress("det2", &det2);
    if (fIsSimData)   fNormTree->SetBranchAddress("B", &fB);
    
    fDet1NormVec    = fCCNorm->GetDet1NormVec    () ;
    fDet2NormVec    = fCCNorm->GetDet2NormVec    () ;
    fRunIdVec       = fCCNorm->GetRunIdVec       () ;
    nEventsInRunVec = fCCNorm->GetEventsInRun1Vec () ; //TODO
    det1max = fCCNorm->GetDet1Max ();
    det2max = fCCNorm->GetDet2Max ();
}


void CentralitySlicesFinder::Fit2DCorrelation ()
{
    std::cout << "Fitting 2D Correlation..." << std::endl;
    fNormTree->Draw( "det2 : det1  >> h1(120, 0., 1.2, 120, 0., 1.2)", fCuts, "colz");    
    fNormTree->Draw("det2 : det1 >> prof (100, 0., 1., 100, 0., 1.)", fCuts, "sameprof");    
    TProfile *p2D = (TProfile*)gPad->GetPrimitive("prof");
    p2D->Fit("pol3", "Q", "", 0.0, 1.0);
    fFitFunction = (TF1*) p2D->GetFunction("pol3");    
    std::cout << "Fitting done!" << std::endl;
}

void CentralitySlicesFinder::FitCorrection ( int n_points )
{
//     TCanvas *c1 = new TCanvas("c1", "canvas", 1000, 600);
    std::cout << "Starting fit correction..." << std::endl;

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
    
    fCorrProfile = new TGraphErrors( X4fit.size(), &(X4fit[0]), &(Y4fit[0]), &(ErrorX[0]), &(ErrorY[0]));
    fCorrProfile->Fit("pol5", "Q", "", 0.00, 1.0);
    fCorrProfile->GetFunction("pol5")->SetLineWidth(2);
    fCorrProfile->GetFunction("pol5")->SetLineColor(1);
    fFitFunction = (TF1*) fCorrProfile->GetFunction("pol5");
    
    std::cout << "Fit correction done!" << std::endl;
}


void CentralitySlicesFinder::FindCentralitySlices ( )
{
//     std::cout << "Start ..." << std::endl;
    gStyle->SetOptStat(0000);
    fSlice->ClearData ();
    
    nIntervals = Int_t (100/fSliceStep);
        
    TString distrName;
    if (!is1DAnalisys)  distrName = Form("%s_%s_centrality_QA", fDet1Name.Data(), fDet2Name.Data() );
    else                distrName = Form("%s_centrality_QA", fDet1Name.Data() );
        
    TCanvas *c1 = new TCanvas("c1", "canvas", 1500, 900);

    if (!is1DAnalisys) Fit2DCorrelation ();    
    if (!is1DAnalisys) FitCorrection (50);   

    if (is1DAnalisys)
    {
        fNormTree->Draw( "det1 >> h1(120, 0., 1.1)", fCuts);        
        c1->SetLogy();
    }
    
    FindSlices ();    

    /*if (!is1DAnalisys)  */FindMeanSignalsInSlice();
    
//     c1->Print( CFdir + "QA/" + distrName + ".pdf");    
//     c1->Print( CFdir + "QA/" + distrName + ".root");    
    
    
//     fSlicesFile = new TFile ( Form(CFdir + "root_files/Slices_%s_%s_%d.root", fDet1Name.Data(), fDet2Name.Data(), RunId ) , "RECREATE");

    fSlice->SetDetMax (det1max, det2max);
    fSlice->SetFitFunction (fFitFunction);
    fSlice->SetSlicesStep (fSliceStep);
    fSlice->SetDirectionCentralEvents (fDirectionCentralEvents);
    
    fSlice->SetDet1NormVec (fDet1NormVec);
    fSlice->SetDet2NormVec (fDet2NormVec);
    fSlice->SetRunIdVec    (fRunIdVec   );
    
    
    fCentrTree->Fill();  
    
//     std::cout << "Slicing done!" << std::endl;
//     QA ();
//     fSlicesFile->Write();    
    
}

void CentralitySlicesFinder::FindSlices ( )
{
    std::cout << "Start slicing..." << std::endl;
    double finalPar[6] = {0, 0, 0, 0, 0, 0};

    if (!is1DAnalisys)
        for (int i=0; i<nFinalPar; i++)
            finalPar[i] = fFitFunction->GetParameter(i);
    
    fNormTree->Draw(">>elist" );
    TEventList* elist = (TEventList*)gDirectory->Get("elist");
    int nTotalEvents = fNormalization > 0 ? fNormalization : elist->GetN();
   
    std::cout << "   TotalEvents = " << nTotalEvents << endl; 

    fNormTree->Draw(">>elist1", fCuts );
    TEventList* elist1 = (TEventList*)gDirectory->Get("elist1");
    int nTotalEvents1 = elist1->GetN();

    std::cout << "   Events with cuts = " << nTotalEvents1 << "   fraction of total = " << nTotalEvents1*100.0/nTotalEvents << " %"  << endl; 

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
        std::cout << "   % of events in slice # " << i << " = " << 100.*nInsideCurrent/nTotalEvents << endl; 
        
        fSlice->AddSlice(kb[1], kb[0], x);
                 
//         std::cout << " X  = " << fSlice->GetXi(i) << endl; 
//         std::cout << " B  = " << fSlice->GetBi(i) << endl; 
//         std::cout << " A = " << fSlice->GetAi(i) << endl; 

    }
    std::cout << "Slicing done!" << std::endl;
}


void CentralitySlicesFinder::FindMeanSignalsInSlice ( void )
{
    std::cout << "Start calculating mean parameters in each slice for QA..." << std::endl;
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
        
//         std::vector <double> X(nInside);
//         std::vector <double> Y(nInside);    
//         std::vector <double> XY(nInside);                   
//         std::vector <double> B(nInside);

        TH1F *hX  = new TH1F ("hX",  "hX",  100, 0, 1);
        TH1F *hY  = new TH1F ("hY",  "hY",  100, 0, 1);
        TH1F *hXY = new TH1F ("hXY", "hXY", 100, -0.5, 0.5);
        TH1F *hB  = new TH1F ("hB",  "hB",  200, 0, 20);
        
        double meanB = 0, meanParX = 0, meanParY = 0, meanParXY = 0;
        double sigmaB = 0, sigmaParX = 0, sigmaParY = 0, sigmaParXY = 0, meanParXY3 = 0, dsigmaB = 0, dB = 0;
        
//         int nn = 0;
        
        for (int k=0; k<nInside; k++)
        {
            fNormTree->GetEntry(insideList->GetEntry(k));

//             X.at(k) = det1;
//             Y.at(k) = det2;
//             B.at(k) = fB; 
//             XY.at(k) = ( X.at(k)*tempk - Y.at(k) + tempb) / TMath::Sqrt( 1 + tempk*tempk);  
//         std::cout << "fill b = " << fB << std::endl;

            hB->Fill(fB);
            hX->Fill(det1);
            if (!is1DAnalisys)  hY->Fill(det2);
            if (!is1DAnalisys)  hXY->Fill(( det1*tempk - det2 + tempb) / TMath::Sqrt( 1 + tempk*tempk));

//             meanB     += B.at(k);
//             meanParX  += X.at(k);
//             meanParY  += Y.at(k);
//             meanParXY += XY.at(k); 
//             nn++;
            
        }
        
//         TCanvas *cTemp = new TCanvas ("cTemp", "cTemp", 800, 600);
//         cTemp->Divide(2,2);
        
//         cTemp->cd(1);
//         hB->Draw();
        
//         cTemp->cd(2);
//         hX->Draw();

//         cTemp->cd(3);
//         hY->Draw();

//         cTemp->cd(4);
//         hXY->Draw();
        hB->Fit("gaus", "Q");
        hX->Fit("gaus", "Q");
        if (!is1DAnalisys)  hY->Fit("gaus", "Q");
        if (!is1DAnalisys)  hXY->Fit("gaus", "Q");
        
//         gPad->Update();
        
//         std::cin >> meanB;
        TF1 *fYfit;
        TF1 *fBfit = hB->GetFunction("gaus");
        TF1 *fXfit = hX->GetFunction("gaus");
        if (!is1DAnalisys)   fYfit = hY->GetFunction("gaus");

        meanB    = fBfit->GetParameter(1);      sigmaB    = fBfit->GetParameter(2);
        meanParX = fXfit->GetParameter(1);      sigmaParX = fXfit->GetParameter(2);
        if (!is1DAnalisys)   { meanParY = fYfit->GetParameter(1);      sigmaParY = fYfit->GetParameter(2); }
        
        dsigmaB = fBfit->GetParError(2);
        dB = fBfit->GetParError(1);
        
//         cout << fBfit->GetParameter(0) << endl;
//         cout << fBfit->GetParameter(1) << endl;        
//         cout << fBfit->GetParameter(2) << endl;        
        
                
        cout << "   Events inside slice = " << nInside << endl;
       
//         if (nn > 0){
//             meanB     /= nn;
//             meanParX  /= nn;
//             meanParY  /= nn;
//             meanParXY /= nn;        
//         }
//         else continue;
        
        
//         for (int k=0; k<nInside; k++)
//         {
//             sigmaB     += TMath::Power(B.at(k) - meanB, 2);
//             sigmaParX  += TMath::Power(X.at(k) - meanParX, 2);
//             sigmaParY  += TMath::Power(Y.at(k) - meanParY, 2);
//             sigmaParXY += TMath::Power(XY.at(k) - meanParXY, 2);
//             meanParXY3 += TMath::Power(XY.at(k) - meanParXY, 3);
//             
//         }

//         sigmaB = TMath::Power(sigmaB/nn, 0.5);
//         sigmaParX = TMath::Power(sigmaParX/nn, 0.5);
//         sigmaParY = TMath::Power(sigmaParY/nn, 0.5);
//         sigmaParXY = TMath::Power(sigmaParXY/nn, 0.5);
//         meanParXY3 = meanParXY3/nn/TMath::Power(sigmaParXY, 3);

        
        fSlice->AddXPar (meanParX, sigmaParX);
        fSlice->AddYPar (meanParY, sigmaParY);
        fSlice->AddXYPar (meanParXY, sigmaParXY);
        fSlice->AddXY3 (meanParXY3);
        fSlice->AddB (meanB, sigmaB, dB, dsigmaB);
        
        hX  ->Delete();
        hY  ->Delete();
        hXY ->Delete();
        hB  ->Delete();
    }
    std::cout << "Mean parameters calculated!" << std::endl;
    
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

void CentralitySlicesFinder::RunSliceFinder ()
{
    std::cout << "Starting SliceFinder" << std::endl;
    fCentrTree = new TTree ("CentrTree", "Tree with slices parameters");
    fCentrTree->Branch("CentralitySlice", "CentralitySlice", &fSlice);

    FindCentralitySlices( );
    
    TString filename = Form(CFdir + "root_files/Slices_%s_%s.root", fDet1Name.Data(), fDet2Name.Data() );
    
    fSlicesFile = new TFile (filename, "RECREATE");
    std::cout << "   Output file name is " << filename << std::endl;
    
    fCentrTree->Write();
    std::cout << "SliceFinder - done!" << std::endl;
    
}

void CentralitySlicesFinder::QA ()
{
    std::cout << "Start QA..." << std::endl;
    TCanvas *c1 = new TCanvas("c1", "canvas", 1500, 900);
    c1->Divide (2, 2);



// ***************************** Run-by-run QA *****************************
    c1->cd(1);
    std::vector <Float_t> TempVec;
    for (UInt_t j=0; j<fRunIdVec.size(); j++)
        TempVec.push_back(fRunIdVec.at(j));

    TGraph *gr1 = new TGraph (Int_t(fRunIdVec.size()), &(TempVec[0]), &(fDet1NormVec[0]) );
    gr1->SetMarkerStyle(23);    
    gr1->GetXaxis()->SetTitle( "Run" );
    gr1->GetYaxis()->SetTitle(  Form ( "%s^{run}/<%s>", fDet1Name.Data(), fDet1Name.Data() ) );    
    gr1->GetYaxis()->SetRangeUser(0.9, 1.1);
 
    gr1->Draw("APL");

    if (fDet2Id != -1){ 
        TGraph *gr2 = new TGraph (Int_t(fRunIdVec.size()), &(TempVec[0]), &(fDet2NormVec[0]) );
        gr1->GetYaxis()->SetTitle(  Form ( "signal^{run}/<signal>" ) );    
        
        gr2->SetMarkerStyle(23);    
        gr2->SetMarkerColor(kRed);    
        gr2->SetLineColor(kRed);    
        gr2->Draw("PLsame");
    
        TLegend* leg1 = new TLegend(0.7,0.75,0.89,0.89);
        leg1->AddEntry(gr1, fDet1Name, "p");    
        leg1->AddEntry(gr2, fDet2Name, "p");    
        leg1->Draw("same");
    }    
// **************************************************************************
    
// ******************************* Slices QA ********************************
    c1->cd(2);    
    std::vector <Float_t> centrality;
    Float_t step = fSlice->GetSlicesStep ();
    UInt_t NSlices = fSlice->GetNSlices();
    
    for (UInt_t i=0; i<NSlices; i++)
        centrality.push_back( (i+0.5)*step );
            
    if (!is1DAnalisys)
    {
        fNormTree->Draw( "det2 : det1  >> h1(1200, 0., 1.2, 1200, 0., 1.2)", fCuts, "colz");
        TH2F *h1 = (TH2F*)gPad->GetPrimitive("h1");
        
        fFitFunction->SetLineColor(kRed);

        h1->SetTitle ("2D correlation");
        h1->GetXaxis()->SetTitle( Form ( "%s/%s_{max}", fDet1Name.Data(), fDet1Name.Data() ) );
        h1->GetYaxis()->SetTitle( Form ( "%s/%s_{max}", fDet2Name.Data(), fDet2Name.Data() ) );   

        fNormTree->Draw("det2 : det1 >> prof (100, 0., 1., 100, 0., 1.)", fCuts, "sameprof");
        TProfile *p2D = (TProfile*)gPad->GetPrimitive("prof");
        p2D->SetMarkerStyle(23);  
        p2D->SetMarkerSize(1);
        p2D->SetLineColor(1);    

        for (UInt_t i=0; i<NSlices; i++)
        {
            Float_t x = fSlice->GetXi(i);
            Float_t xNormInterval = 0.15 * TMath::Abs ( TMath::Cos (TMath::ATan (fSlice->GetBi(i)) ) );
            TF1 *norm = new TF1("norm","[0]*x + [1]", x-xNormInterval, x+xNormInterval); 
            norm->SetParameters (fSlice->GetBi(i), fSlice->GetAi(i));
            norm->Draw("same"); 
        }

        fCorrProfile->SetMarkerStyle(22);  
        fCorrProfile->SetMarkerSize(1);
        fCorrProfile->SetMarkerColor(kRed);
    
        TLegend* leg1 = new TLegend(0.6,0.6,0.9,0.9);
        leg1->AddEntry(p2D,          "profile", "p");            
        leg1->AddEntry(fCorrProfile, "corrected profile", "p");    
        leg1->AddEntry(fFitFunction, "corrected fit", "l");            
    
        fFitFunction->Draw("same");
        p2D->Draw("same");
        fCorrProfile->Draw("Psame");
        leg1->Draw("same");    
        gPad->Update();
        
    }
    else
    {
        fNormTree->Draw( "det1  >> h1(120, 0., 1.2)", fCuts);
        TH1F *h1 = (TH1F*)gPad->GetPrimitive("h1");
        h1->GetXaxis()->SetTitle( Form ("<%s>/%s^{max}", fDet1Name.Data(), fDet1Name.Data()) );
        h1->GetYaxis()->SetTitle( "Counts" );    
        gPad->SetLogy();
        
        for (UInt_t i=0; i<NSlices; i++)
        {
            gPad->Update();
            Float_t x = fSlice->GetXi(i);
            Double_t xr = gPad->GetX2()-gPad->GetX1();
            double x1 = (x-gPad->GetX1()) / xr;
            TLine l2;
            l2.SetLineColor(kRed);
            l2.DrawLineNDC(x1, 0.10, x1, 0.9);             
        }
    }
// **************************************************************************

    c1->cd(3);
    TGraphErrors *grX = new TGraphErrors(NSlices, &(centrality[0]), &(fSlice->GetMeanX()[0]), 0, &(fSlice->GetSigmaX()[0]));
    grX->SetTitle ("X, Y vs centrality");
    grX->SetMarkerStyle(23);    
    grX->GetXaxis()->SetTitle( "Centrality, %" );
    grX->GetYaxis()->SetTitle( Form ("<%s>/%s^{max}", fDet1Name.Data(), fDet1Name.Data()) );    

    grX->GetXaxis()->SetLimits(0.0, NSlices*step);
    grX->GetYaxis()->SetRangeUser(0.0, 1.0);
    grX->Draw("APL");

    TGraphErrors *grY;
    if (!is1DAnalisys)
    {
        grY = new TGraphErrors(NSlices, &(centrality[0]), &(fSlice->GetMeanY()[0]), 0, &(fSlice->GetSigmaY()[0]));
        grY->SetMarkerStyle(22);    
        grY->SetMarkerColor(kRed);    
        grY->SetLineColor(kRed);    
        grY->Draw("PLsame");

        grX->GetYaxis()->SetTitle( "<Signal>/Signal_{max}" );    
        
        TLegend* leg2 = new TLegend(0.7,0.8,0.9,0.9);
        leg2->AddEntry(grX, fDet1Name.Data(), "p");            
        leg2->AddEntry(grY, fDet2Name.Data(), "p");    
        leg2->Draw("same");     
    }
    
    
    
    c1->cd(4);
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

    TGraphErrors *grB = new TGraphErrors (NSlices, &(centrality[0]), &(TempVecB[0]), 0, &(TempdB[0]));
    grB->SetMarkerStyle(23);    
    grB->GetXaxis()->SetTitle( "Centrality" );
    grB->GetYaxis()->SetTitle( "sigma_{B}/<B>" );    
    grB->GetYaxis()->SetRangeUser(0.05, 0.5);
      
    grB->Draw("APL");
    
//     c1->cd(3);
//     
//     TGraphErrors *grXY = new TGraphErrors(NSlices, &(centrality[0]), &(fSlice->GetMeanXY()[0]), 0, &(fSlice->GetSigmaXY()[0]));
//     grXY->SetTitle ("2D fit: XY vs centrality");
//     grXY->SetMarkerStyle(23);    
//     grXY->GetXaxis()->SetTitle( "Centrality, %" );
//     grXY->GetYaxis()->SetTitle("XY");    
//     grXY->GetXaxis()->SetLimits(0.0, NSlices*step);
//     grXY->GetYaxis()->SetRangeUser(-0.2, 0.2);
//     grXY->Draw("APL");
    
//     c1->cd(4);
    
//     TGraphErrors *grXY3 = new TGraphErrors(NSlices, &(centrality[0]), &(fSlice->GetMeanXY3()[0]));
//     grXY3->SetTitle ("2D fit: XY3 vs centrality");
//     grXY3->SetMarkerStyle(23);    
//     grXY3->GetXaxis()->SetTitle( "Centrality, %" );
//     grXY3->GetYaxis()->SetTitle("XY3");    
//     grXY3->GetXaxis()->SetLimits(0.0, NSlices*step);
//     grXY3->GetYaxis()->SetRangeUser(-5, 5);
//     grXY3->Draw("APL");
    
    TString filename = Form(CFdir + "QA/QA_%s_%s", fDet1Name.Data(), fDet2Name.Data() );
    c1->Print( filename + ".pdf" );    
    c1->Print( filename + ".root" );    
    
    std::cout << "QA saved!" << std::endl;
//     TGraph g(a.size(), &(a[0]), &(b[0]));
    
    
}

