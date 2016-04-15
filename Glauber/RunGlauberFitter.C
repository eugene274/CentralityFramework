

void RunGlauberFitter (Int_t nf = 1, Float_t f0 = 0.95, Float_t f1 = 0.95, Int_t nsigma = 1, Float_t StepSigma = 0.05, Int_t nEvents = 100000)
{
//     gSystem->Load("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/build/libCentrality");
 
//     TH1F *h1d_Mult_TPC_Ref = new TH1F("h1d_Mult_TPC_Ref","M_{TPC}^{ref}", 125, 0, 250);

//     TFile *f = new TFile("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Glauber/root_files/in_histo_1.root");    
//     h1d_Mult_TPC_Ref = (TH1F*)f->Get("h1d_Mult_TPC_Ref");
    
//     TH1F *hIn = new TH1F("hIn","", 275, 0, 1100);
//     TFile *f = new TFile("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Glauber/Multiplicity.root");    
//     hIn = (TH1F*)f->Get("h1");
    
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";
    gStyle->SetOptStat(0000);    
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");    
        
    TString DataFileName = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/root_files/DCM-QGSM/sis100_electron_SC_ON_SL_OFF.root";
    
    TFile *DataFile = new TFile ( DataFileName, "read" );    

    CentralityEventContainer *container = new CentralityEventContainer;
    TTree *ContTree = (TTree*)DataFile->Get("cbm_data");
    ContTree->SetBranchAddress("CentralityEventContainer", &container);
    
    Int_t n = ContTree->GetEntries();
    std::cout << "n = " << n  << std::endl;

    TCut psd1_mult = "";//" CentralityEventContainer.GetDetectorWeight(3) > 250 - 250.0/40 * CentralityEventContainer.GetDetectorWeight(0) ";
    
    TCanvas *c1 = new TCanvas("c1", "canvas", 1000, 800);
//     c1->Divide(3, 3);
    
//     c1->cd(1);
    ContTree->Draw("CentralityEventContainer.GetDetectorWeight(3) >> h1(100, 0, 400)", psd1_mult);
    TH1F *hData1 = (TH1F*)gPad->GetPrimitive("h1");
    hData1->GetXaxis()->SetTitle( "M_{STS}" );
    hData1->GetYaxis()->SetTitle("Counts");   

    hData1->Draw();
    GlauberFitter *test = new GlauberFitter;

    test->SetNpartMax (450);
    test->SetNcollMax (1400);
    test->SetMultMax (400);
    test->SetBinSize (4);
    
    test->SetFitMultMin (20);
    test->SetNormMultMin (30);
    test->SetOutDirName ("");
    
    test->SetSimHistos ("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Glauber/root_files/glau_pbpb_ntuple_signn_40_1e5ev.root", nEvents);
    test->SetInputHisto (hData1);
       
    test->SetGlauberFitHisto (0.95, 0.778935, 182.491, 100000, true);
    TH1F *h1 = test->GetGlauberFitHisto ();
    h1->Draw("same");
    
    TLegend* legData = new TLegend(0.6,0.75,0.75,0.83);
    legData->AddEntry(h1 ,"Fit", "l");    
    legData->AddEntry(hData1 ,"M_{STS}", "l");    
    legData->Draw("same");           


//     test->FitGlauber(nf, f0, f1, nsigma, StepSigma, nEvents);
//     test->DrawHistos(true, true, false, false);
}