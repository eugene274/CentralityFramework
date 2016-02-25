

void RunGlauberFitter (Int_t nf = 5, Float_t f0 = 0.0, Float_t f1 = 1.0, Int_t nsigma = 3, Float_t StepSigma = 0.05, Int_t nEvents = 1000000)
{
    gSystem->Load("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/build/libCentrality");
 
    TH1F *h1d_Mult_TPC_Ref = new TH1F("h1d_Mult_TPC_Ref","M_{TPC}^{ref}", 125, 0, 250);

    TFile *f = new TFile("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Glauber/root_files/in_histo_1.root");    
    h1d_Mult_TPC_Ref = (TH1F*)f->Get("h1d_Mult_TPC_Ref");
    
    h1d_Mult_TPC_Ref->Draw();
    GlauberFitter *test = new GlauberFitter;

    test->SetNpartMax (450);
    test->SetNcollMax (1400);
    test->SetMultMax (250);
    test->SetBinSize (2);
    
    test->SetFitMultMin (150);
    test->SetNormMultMin (100);
    test->SetOutDirName ("TestMin_150");
    
    
    test->SetSimHistos ("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Glauber/root_files/glau_pbpb_ntuple.root", nEvents);
    test->SetInputHisto (h1d_Mult_TPC_Ref);
   
    test->SetGlauberFitHisto (0.7, 0.266134, 0.559225, 1000000, true);
    TH1F *h1 = test->GetGlauberFitHisto ();
    h1->Draw("same");
    
    TLegend* legData = new TLegend(0.6,0.75,0.75,0.83);
    legData->AddEntry(h1 ,"Fit", "l");    
    legData->AddEntry(h1d_Mult_TPC_Ref ,"Data", "l");    
    legData->Draw("same");           


//     test->DrawHistos(true, true, false, false);
//     test->FitGlauber(nf, f0, f1, nsigma, StepSigma, nEvents);

}