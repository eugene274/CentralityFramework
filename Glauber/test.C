#include <iostream>
#include <fstream>
#include <vector>

#include "TF1.h"
#include "TString.h"


const int nNcoll = 1400;
const int nNpart = 450;

const int nBinNcoll = 1400;
const int nBinNpart = 450;

Double_t NcollArr[nBinNcoll];
Double_t NpartArr[nBinNpart];

Double_t fNcoll (const Double_t *x, const Double_t *par)
{
    Double_t ret = 0;
    if ( par[0]*x[0] < nBinNcoll-1 && par[0]*x[0] >= 0)
    {   
        Float_t n = par[0]*x[0];        
        ret = NcollArr[(int)n];
//         ret += (n-(int)n)*NcollArr[(int)(n+1)];
    }
    return ret;    
}

Double_t fNpart (const Double_t *x, const Double_t *par)
{
    Double_t ret = 0;
    if ( par[0]*x[0] < nBinNpart && par[0]*x[0] >= 0)
    {
        Float_t n = par[0]*x[0];        
        ret = NpartArr[(int)n];
    }
    return ret;
}

Double_t fSimpleFit (const Double_t *x, const Double_t *par)
{
    return par[3]*( par[0]*fNpart( x, &(par[1]) ) + (1-par[0])*fNcoll( x, &(par[2]) ) );
}




Double_t fFit (const Double_t *x, const Double_t *par)
{

    Double_t ret = 0;
//     if (x0 < 0)   return 0;
    

    TF1 * funcNBD = new TF1("funcNBD","ROOT::Math::negative_binomial_pdf(x,[0], [1])", 0, 300);
//     TF1 * funcNBD = new TF1("funcNBD","ROOT::Math::gaussian_pdf(x, [0], [1])", 0, 300);
    
//     funcNBD->SetParameter(0, par[3]);// normalization constant
    funcNBD->SetParameter(0, par[4]); // k parameter      (for gauss - sigma)
    funcNBD->SetParameter(1, par[5]); // mean multiplicity        (for gauss - mean)
    
    
    Double_t mean = (1-par[4])/par[4]*par[5];
    Double_t sigma = TMath::Sqrt(mean/par[4] );
 
//     std::cout << "mean = " << mean << "      sigma = " << sigma << std::endl;    
    
    Int_t x0 = int(x[0] - mean);
    Int_t dx = int( 3*sigma );

    Int_t i=x0+dx;
    
    while (i>=0 && i >= x0-dx)
    {
//         std::cout << "i in = " << i << std::endl;
        Double_t func = dx > 0 ? funcNBD->Eval(i-x0) : 1.0;
        Double_t xx = i;
        ret += func * ( par[0]*fNpart( &xx, &(par[1]) ) + (1-par[0])*fNcoll( &xx, &(par[2]) ) );
        
//         std::cout << "ret = " << ret << std::endl;

        i--;
    }
    return par[3]*ret;
}

void test(void)
{

    TCanvas *c1 = new TCanvas("c1", "canvas", 1500, 900);
    c1->Divide (2, 2);
    

    TFile *fData = new TFile ("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/root_files/NA61_run_23001.root");
    TTree* fDataTree = (TTree*) fData->Get("na61_data");
    
    TFile *fSim = new TFile ("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Glauber/glau_pbpb_ntuple.root");
    TTree* fSimTree = (TTree*) fSim->Get("nt_Pb_Pb");

    TPad *c1_1 = (TPad*) c1->GetListOfPrimitives()->FindObject("c1_1");
    c1_1->SetLogy(1);
    TPad *c1_2 = (TPad*) c1->GetListOfPrimitives()->FindObject("c1_2");
    c1_2->SetLogy(1);
    TPad *c1_4 = (TPad*) c1->GetListOfPrimitives()->FindObject("c1_4");
    c1_4->SetLogy(1);

    TH1F *hNpart, *hNcoll, *hData;
    
    c1->cd(1);
    
    TString DrawPar1 = Form("Ncoll >> hNcoll(%d, 0, %d)", nBinNcoll, nNcoll  );
    TString DrawPar2 = Form("Npart >> hNpart(%d, 0, %d)", nBinNpart, nNpart  );
    
    fSimTree->Draw( DrawPar1.Data() );
    hNcoll = (TH1F*)gPad->GetPrimitive("hNcoll");    
    
    fSimTree->Draw( DrawPar2.Data()  );
    hNpart = (TH1F*)gPad->GetPrimitive("hNpart");
    hNcoll->SetLineColor(2);

    hNcoll->Draw();
    hNpart->Draw("same");
    
    TLegend* legSim = new TLegend(0.6,0.75,0.75,0.83);
    legSim->AddEntry(hNpart ,"Npart", "l");    
    legSim->AddEntry(hNcoll ,"Ncoll", "l");    
    legSim->Draw("same");
    
    c1->cd(2);
    fDataTree->Draw("CentralityEventContainer.fDetectorEvents.fWeights >> hData (120, 0, 240)", "CentralityEventContainer.fDetectorEvents.fDetId == 3");
    hData = (TH1F*)gPad->GetPrimitive("hData");
    hData->GetYaxis()->SetRangeUser(0.1, 1e4);
   
    Double_t xmin = 0, xmax = 250;

    for (Int_t i=0; i<nBinNpart; i++)
        NpartArr[i] = hNpart->GetBinContent(i+1);
 
    for (Int_t i=0; i<nBinNcoll; i++)
        NcollArr[i] = hNcoll->GetBinContent(i+1);
    
    
//     for (Int_t i=0; i<nBinNpart; i++)
//         std::cout << "NcollArr = " << NcollArr[i] << "         NpartArr = " << NpartArr[i] << std::endl;
    
    TF1 *f = new TF1("f", /*fFit*/fSimpleFit, xmin, xmax, 4);     
    
    TF1 *fPart = new TF1("fPart", fNpart, xmin, xmax, 1);     
    TF1 *fColl = new TF1("fColl", fNcoll, xmin, xmax, 1);     
    
    Double_t SimMax = 200.;
    Double_t NcollNorm = nBinNcoll/SimMax;
    Double_t NpartNorm = nBinNpart/SimMax;
    Double_t NormY = 0.2;
    
    
    f->SetParameter(0, 0.1);// f (<1)    
    f->SetParameter(1, NpartNorm);   // normalization constant Npart   
    f->SetParameter(2, NcollNorm);   // normalization constant Ncoll   
    f->SetParameter(3, NormY);       // normalization constant Y   
    f->SetParameter(4, 0.7);           // k (<1)     sigma
    f->SetParameter(5, 10);          //mean

//     f->FixParameter (4, 0.5);
//     f->FixParameter (5, 1);
    f->SetParNames ("f", "NpartNorm", "NcollNorm", "Norm", "k", "r");
    
    f->SetParLimits(0, 0, 1);    
    f->SetParLimits(1, 0.5*NpartNorm, 1.5*NpartNorm);    
    f->SetParLimits(2, 0.5*NcollNorm, 1.5*NcollNorm);    
    f->SetParLimits(3, 0, 100);    
    f->SetParLimits(4, 0, 1);    
    f->SetParLimits(5, 0, 50);    

    f->Draw("same");

//     f->SetLineColor(3);
    
    hData->Fit("f","V","", 20, 210);
//     f->Draw("same");
//     fPart->Draw("same");
//     fColl->Draw("same");
//     hData->Fit("NBD");
    
    c1->cd(3);
    TF1 * funcNBD = new TF1("funcNBD","[0]*ROOT::Math::negative_binomial_pdf(x,[1],[2])", 0, 100);
    funcNBD->SetParameter(0, 1);// normalization constant
    funcNBD->SetParameter(1, 0.5); // k parameter
    funcNBD->SetParameter(2, 10); // mean multiplicity        
    funcNBD->Draw();
    gPad->Update();
        
    c1->cd(4);
    fPart->SetLineColor(3);
    fPart->SetParameter(0, NpartNorm);
    fColl->SetParameter(0, NcollNorm);

    fPart->Draw();
    fColl->Draw("same");
    TLegend* legFunc = new TLegend(0.6,0.75,0.75,0.83);
    legFunc->AddEntry(fPart ,"Npart", "l");    
    legFunc->AddEntry(fColl ,"Ncoll", "l");    
    legFunc->Draw("same");
    
}