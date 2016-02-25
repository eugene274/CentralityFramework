#include "GlauberFitter.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom.h"

ClassImp(GlauberFitter)


// -----   Default constructor   -------------------------------------------
GlauberFitter::GlauberFitter() 
  : TNamed()

{

}

void GlauberFitter::TestFunc(Int_t nf, Float_t f0, Float_t f1, Int_t nsigma, Int_t nEvents)
{

//     Float_t mu, k, sigma;
//     
//     mu = 0.2;
//     sigma = 2;
//     k = mu/(sigma/mu - 1);
// 
//     SetNBDhist(mu,  k);
// 
//     Float_t Na = 50000;
//     TH1F *test = new TH1F ("test", "", 50, 0, 50);
//     Float_t nHits = 0;
//     for (Int_t j=0; j<Na; j++){
//         Float_t temp = (int)hNBD->GetRandom();
//         nHits += temp;
//         test->Fill(temp);
// //         std::cout << "GetRandom = " << temp << std::endl;
//     }
//     test->SetLineColor(2);
//     
//     hNBD->DrawNormalized();
//     test->DrawNormalized("same");
    
    
    TString filename = Form ( "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Glauber/root_files/test_%d_%f_%f.root", nf, f0, f1 );
    TFile *file = new TFile(filename, "recreate");    
    TTree *tree = new TTree("test_tree", "tree" );
    
    TH1F *h1 = new TH1F ("h1", "", 125, 0, 250);
           
    Float_t f, mu, k;
    Float_t chi2 = 1e10;
    Float_t sigma;

    tree->Branch("histo", "TH1F", &h1, 32000, 0);
    tree->Branch("f",    &f,    "f/F");   
    tree->Branch("mu",   &mu,   "mu/F");   
    tree->Branch("k",    &k,    "k/F");   
    tree->Branch("chi2", &chi2, "chi2/F");   
    tree->Branch("sigma",&sigma,"sigma/F");   


    for (Int_t i=0; i<nf; i++)
    {
        Float_t temp_f = f0 + (f1-f0)/nf*(i+0.5);
        
        for (Int_t j=0; j<nsigma; j++)
        {
            Float_t temp_mu = fMultMax/(temp_f*(fNpartMax - 75) + (1-temp_f)*(fNcollMax - 150));
            sigma = temp_mu + (j+0.5)/nsigma;
            Float_t temp_k = temp_mu/(sigma/temp_mu - 1);;
            Float_t temp_chi2;
           
            std::cout << "f = " << temp_f << "   mu0 = " << temp_mu << "    k = " << temp_k << "    sigma = " << sigma << std::endl;    
            
            Float_t mu_min = 0.4*temp_mu, mu_max = 0.9*temp_mu;

            FindMuGoldenSection (&temp_mu, &temp_chi2, mu_min, mu_max, temp_f, temp_k, nEvents, 10);

            f = temp_f;
            k = temp_k;
            mu = temp_mu;
            chi2 = temp_chi2;
            h1 = hGlaub;
            
            tree->Fill();
            std::cout << "f = " << temp_f << "   mu1 = " << temp_mu << "    k = " << temp_k << "                Chi2Min = " << temp_chi2 << std::endl;    
        } 
    }
    
    tree->Write ();

    file->Close();
    
}


void GlauberFitter::SetSimHistos(TString filename, Int_t nEntries)
{    
    TFile *fSim = new TFile (filename.Data());
    fSimTree = (TTree*) fSim->Get("nt_Pb_Pb");
    
    if (!fSimTree) {
        std::cout << "SetSimHistos: *** Error - " << std::endl;
        exit(EXIT_FAILURE);
    }
    
    fSimTree->SetBranchAddress("Npart", &fNpart);
    fSimTree->SetBranchAddress("Ncoll", &fNcoll);
    
    if (nEntries < 0){
        std::cout << "SetSimHistos: *** Warning - Number of Entries < 0 or not set. Default value will be used" << std::endl;
        nEntries = fSimTree->GetEntries();
    }
    
    hNpart = new TH1F ("hNpart", "", fNpartMax/fBinSize, 0, fNpartMax );
    hNcoll = new TH1F ("hNcoll", "", fNcollMax/fBinSize, 0, fNcollMax );
    
    for (Int_t i=0; i<nEntries; i++)
    {
        fSimTree->GetEntry(i);
        hNcoll->Fill(fNcoll);
        hNpart->Fill(fNpart);
    }
//     fSim->Close();
}

void GlauberFitter::DrawHistos (Bool_t isSim, Bool_t isData, Bool_t isGlauber, Bool_t isNBD )
{
    TCanvas *c1 = new TCanvas("c1", "canvas", 1500, 900);
    c1->Divide(2,2);
    TPad *c1_1 = (TPad*) c1->GetListOfPrimitives()->FindObject("c1_1");
    c1_1->SetLogy(1);
    TPad *c1_2 = (TPad*) c1->GetListOfPrimitives()->FindObject("c1_2");
    c1_2->SetLogy(1);
    TPad *c1_4 = (TPad*) c1->GetListOfPrimitives()->FindObject("c1_4");
    c1_4->SetLogy(1);
   
    if (isSim){
        c1->cd(1);
        hNcoll->SetLineColor(2);
    
        hNcoll->Draw();
        hNpart->Draw("same");
        
        TLegend* legSim = new TLegend(0.6,0.75,0.75,0.83); 
        legSim->AddEntry(hNpart ,"Npart", "l");    
        legSim->AddEntry(hNcoll ,"hNcoll", "l");    
        legSim->Draw("same");          
    }
  
    if (isData){
        c1->cd(2);
        hData->GetYaxis()->SetRangeUser(0.1, 2e5);    
        hData->Draw();
        
        if (isGlauber){
            hBestFit->SetLineColor (kRed);
            hBestFit->Draw("same");
            
            TLegend* legData = new TLegend(0.6,0.75,0.75,0.83);
            legData->AddEntry(hBestFit ,"Fit", "l");    
            legData->AddEntry(hData ,"Data", "l");    
            legData->Draw("same");           
        }
    }
 
    if (isNBD){
        c1->cd(3);
        hNBD->Draw();    
    }
    
    if (isGlauber){
        c1->cd(4);
        hBestFit->Draw();
    }
}

void GlauberFitter::SetGlauberFitHisto (Float_t f, Float_t mu, Float_t k, Int_t n, Bool_t Norm2Data)
{
    Float_t Na;    
    SetNBDhist(mu,  k);

    hGlaub = new TH1F("htemp", "", fMultMax/fBinSize, 0, fMultMax);
    hGlaub->SetName(Form("glaub_%4.2f_%6.4f_%4.2f_%d", f, mu, k, n ));

    for (Int_t i=0; i<n; i++)
    {
        fSimTree->GetEntry(i);
        Na = f*fNpart + (1-f)*fNcoll;
//         Na = TMath::Power(fNpart, f);    //TODO Try later
//         Na = TMath::Power(fNcoll, f);    //TODO Try later
        
        Float_t nHits = 0;
        for (Int_t j=0; j<Na; j++){
            nHits += (int)hNBD->GetRandom();
//             std::cout << hNBD->GetRandom() << std::endl;
        }
        hGlaub->Fill(nHits);
    }
    
    if (Norm2Data)
        NormalizeGlauberFit(fNormMultMin/fBinSize);
}


void GlauberFitter::NormalizeGlauberFit (Int_t StartBin)
{
    
    Int_t hGlaubInt = 0, hDataInt = 0;
    for (Int_t i=StartBin; i<fMultMax/fBinSize; i++)
    {
        hGlaubInt += hGlaub->GetBinContent(i+1);
        hDataInt += hData->GetBinContent(i+1);
    }

    Float_t ScaleFactor = (float)hDataInt/hGlaubInt;
    
//     std::cout << "Scale = " << Scale << std::endl;
    hGlaub->Scale(ScaleFactor);    
}

void GlauberFitter::FindMuGoldenSection (Float_t *mu, Float_t *chi2, Float_t mu_min, Float_t mu_max, Float_t f, Float_t k, Int_t nEvents, Int_t nIter)
{
    Float_t phi = (1+TMath::Sqrt(5))/2;
   
    Float_t mu_1 = mu_max - (mu_max-mu_min)/phi; 
    Float_t mu_2 = mu_min + (mu_max-mu_min)/phi;
    
    SetGlauberFitHisto (f, mu_1, k, nEvents);
    Double_t chi2_mu1 = GetChi2 ();
    
    SetGlauberFitHisto (f, mu_2, k, nEvents);
    Double_t chi2_mu2 = GetChi2 ();
    
    for (Int_t j=0; j<nIter; j++)
    {
//         std::cout << "mu_min = " << mu_min << "    mu1 = " << mu_1 << "    mu2 = " << mu_2 << "    mu_max = " << mu_max <<  std::endl;
        
        if (chi2_mu1 > chi2_mu2)
        {
            mu_min = mu_1;
            mu_1 = mu_2;
            mu_2 = mu_min + (mu_max-mu_min)/phi;
            chi2_mu1 = chi2_mu2;
            SetGlauberFitHisto (f, mu_2, k, nEvents);
            chi2_mu2 = GetChi2 ();
        }
        else
        {
            mu_max = mu_2;
            mu_2 = mu_1;
            mu_1 = mu_max - (mu_max-mu_min)/phi; 
            chi2_mu2 = chi2_mu1;
            SetGlauberFitHisto (f, mu_1, k, nEvents);
            chi2_mu1 = GetChi2 ();            
        }
        
        std::cout << "mu1 = " << mu_1 << "    mu2 = " << mu_2 << "    chi2_mu1 = " << chi2_mu1  << "    chi2_mu2 = " << chi2_mu2 << std::endl;
    }

    *mu = (chi2_mu1 < chi2_mu2) ? mu_1 : mu_2;
    *chi2 = (chi2_mu1 < chi2_mu2) ? chi2_mu1 : chi2_mu2;

//     std::cout << "mu = " << *mu << "    f = " << f << "    k = " << k  << "    chi2 = " << *chi2 << std::endl;
}


void GlauberFitter::FitGlauber (Int_t nf, Float_t f0, Float_t f1, Int_t nsigma, Float_t SigmaStep, Int_t nEvents)
{
    Float_t f_fit, mu_fit, k_fit;
    Float_t Chi2Min = 1e10;
    TString dirname = fOutDirName;

    TString filename = Form ( "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Glauber/root_files/%s/test_%d_%4.2f_%4.2f.root", dirname.Data(), nf, f0, f1 );
    TFile *file = new TFile(filename, "recreate");    
    TTree *tree = new TTree("test_tree", "tree" );
    
    TH1F *h1 = new TH1F ("h1", "", 125, 0, 250);
           
    Float_t f, mu, k, chi2, sigma;

    tree->Branch("histo", "TH1F", &h1, 32000, 0);
    tree->Branch("f",    &f,    "f/F");   
    tree->Branch("mu",   &mu,   "mu/F");   
    tree->Branch("k",    &k,    "k/F");   
    tree->Branch("chi2", &chi2, "chi2/F");   
    tree->Branch("sigma",&sigma,"sigma/F");   


    for (Int_t i=0; i<nf; i++)
    {
        Float_t temp_f = f0 + (f1-f0)/nf*(i+0.5);
        for (Int_t j=0; j<nsigma; j++)
        {
            Float_t temp_mu = fMultMax/(temp_f*(fNpartMax - 75) + (1-temp_f)*(fNcollMax - 150));
            sigma = temp_mu + (j+0.1)*SigmaStep;
            Float_t temp_k = temp_mu/(sigma/temp_mu - 1);;
            Float_t temp_chi2;
           
            std::cout << "f = " << temp_f << "   mu0 = " << temp_mu << "    k = " << temp_k << "    sigma = " << sigma << std::endl;    
            
            Float_t mu_min = 0.4*temp_mu, mu_max = 0.9*temp_mu;

            FindMuGoldenSection (&temp_mu, &temp_chi2, mu_min, mu_max, temp_f, temp_k, nEvents, 15);

            f = temp_f;
            k = temp_k;
            mu = temp_mu;
            chi2 = temp_chi2;
            h1 = hGlaub;
            
            tree->Fill();
            std::cout << "f = " << temp_f << "   mu1 = " << temp_mu << "    k = " << temp_k << "                chi2 = " << temp_chi2 << std::endl;    
            
            if (chi2 < Chi2Min)
            {
                f_fit = temp_f;
                mu_fit = temp_mu;
                k_fit = temp_k;
                Chi2Min = chi2;
                hBestFit = hGlaub;
            }            
            
            
        } 
    }
    
    tree->Write ();
    file->Close();

    std::cout << "f = " << f_fit << "    mu = " << mu_fit << "    k = " << k_fit << "    Chi2Min = " << Chi2Min << std::endl;    
}

Double_t GlauberFitter::GetChi2 ()
{
    Double_t chi2 = 0.0;
    
    Int_t lowchibin = fFitMultMin/fBinSize, highchibin = fMultMax/fBinSize;
    
    for (Int_t i=lowchibin; i<=highchibin; i++) 
    {
        if (hData->GetBinContent(i) < 1.0) continue;
        Double_t diff = TMath::Power((hGlaub->GetBinContent(i) - hData->GetBinContent(i)),2);
        diff = diff / (TMath::Power(hData->GetBinError(i),2) + TMath::Power(hGlaub->GetBinError(i),2)); 
        chi2 += diff;
    }
    chi2 = chi2 / (highchibin - lowchibin + 1);
    return chi2;
}

void GlauberFitter::SetNBDhist(Double_t mu, Double_t k)
{
    // Interface for TH1F.
    Int_t nBins = 50;
    
    
    hNBD = new TH1F("hNBD","",nBins, 0, nBins);
    hNBD->SetName(Form("nbd_%f_%f",mu,k));
    hNBD->SetDirectory(0);
    
    for (Int_t i=0; i<nBins; ++i) 
    {
        Double_t val = NBD(i, mu, k);
        if (val>1e-20) hNBD->SetBinContent(i+1, val);
//         std::cout << "val " << val << std::endl;    
    }

}


Double_t GlauberFitter::NBD(Double_t n, Double_t mu, Double_t k) const
{
    // Compute NBD.
    Double_t F;
    Double_t f;

    if (n+k > 100.0) 
    {
        // log method for handling large numbers
        F  = TMath::LnGamma(n + k)- TMath::LnGamma(n + 1.)- TMath::LnGamma(k);
        f  = n * TMath::Log(mu/k) - (n + k) * TMath::Log(1.0 + mu/k);
        F = F+f;
        F = TMath::Exp(F);
    } 
    else 
    {
        F  = TMath::Gamma(n + k) / ( TMath::Gamma(n + 1.) * TMath::Gamma(k) );
        f  = n * TMath::Log(mu/k) - (n + k) * TMath::Log(1.0 + mu/k);
        f  = TMath::Exp(f);
        F *= f;
    }

    return F;
}

/*
void GlauberFitter::MinuitFitFCN (Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // FCN for minuit.

    Double_t ff = par[0];
    Double_t mu    = par[1];
    Double_t k     = par[2];
  
    GlauberFitter * obj = (GlauberFitter *) gMinuit->GetObjectFit();
    
    obj->SetGlauberFitHisto (ff, mu, k);
    Double_t chi2 = obj->GetChi2 ();

    f = chi2;
    
    Printf("Minuit step: chi2=%f, f=%f,  mu=%f, k=%f\n", f, ff, mu, k);
}*/

/*
void GlauberFitter::MinuitFit(Double_t f, Double_t mu, Double_t k)
{
    
    SetSimHistos();
    
    Float_t mu_min = 0.05;
    Float_t mu_max = 0.5;
    Float_t mu_step = 0.001;
    
    Float_t f_min = 0.0;
    Float_t f_max = 1.0;
    Float_t f_step = 0.01;
    
    Float_t k_min = 0.1;
    Float_t k_max = 5.0;
    Float_t k_step = 0.1;

  
    Double_t arglist[10];
    Int_t ierflg = 0;
    arglist[0] = 1;
    
//     if(gMinuit) delete gMinuit;
    TMinuit * gMinuit = new TMinuit(3);  //initialize TMinuit with a maximum of 3 params
    gMinuit->mncler();
    gMinuit->SetFCN(GlauberFitter::MinuitFitFCN);
    gMinuit->SetObjectFit(this); 
    
    gMinuit->mnexcm("SET ERR", arglist ,1, ierflg);

    gMinuit->mnparm(0,"f",   f,   f_step,   f_min,   f_max, ierflg);
    gMinuit->mnparm(1,"mu", mu,  mu_step,  mu_min,  mu_max, ierflg);
    gMinuit->mnparm(2,"k",   k,   k_step,   k_min,   k_max, ierflg);
    
    arglist[0] = 1000; // max calls
    arglist[1] = 0.1;  // tolerance    
    gMinuit->mnexcm("MIGrad",arglist, 2, ierflg);
    
    // ______________________ Get chi2 and Fit Status __________
    Double_t amin,edm,errdef;
    Int_t    nvpar, nparx, icstat;
    
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    Double_t chi2min = amin;
    std::cout << "Fit status " << icstat << std::endl;

    Double_t f_fit, mu_fit, k_fit;
    Double_t alpha_mine, mu_mine, k_mine;
    gMinuit->GetParameter(0, f_fit , alpha_mine );
    gMinuit->GetParameter(1, mu_fit    , mu_mine    );
    gMinuit->GetParameter(2, k_fit     , k_mine     );

    // print some infos
    std::cout << "chi2 min is " << chi2min << ", " << f_fit << ", "<< mu_fit<< ", "  << k_fit << std::endl;    


    SetGlauberFitHisto ( f_fit, mu_fit, k_fit);
    DrawHistos ();          
              
}*/
