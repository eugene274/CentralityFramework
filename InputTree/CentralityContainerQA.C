
int CentralityContainerQA (TString DataFileName, Int_t max = 30)
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";
    
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");    
        
    TFile *DataFile = new TFile ( DataFileName, "read" );    

    CentralityEventContainer *container = new CentralityEventContainer;
    TTree *ContTree = (TTree*)DataFile->Get("cbm_data");
    ContTree->SetBranchAddress("CentralityEventContainer", &container);
    
    Float_t PSD1, PSD2, PSD3, Msts;
    Float_t Centrality;
    
    Int_t n = ContTree->GetEntries();
    std::cout << "n = " << n  << std::endl;

    TCut psd1_mult = "";//" CentralityEventContainer.GetDetectorWeight(3) > 250 - 250.0/40 * CentralityEventContainer.GetDetectorWeight(0) ";
    
    TCanvas *c1 = new TCanvas("c1", "canvas", 1500, 800);
    c1->Divide(2,2);
    
    c1->cd(1);
    ContTree->Draw("CentralityEventContainer.GetDetectorWeight(0) : CentralityEventContainer.GetDetectorWeight(3) >> h1(100, 0, 400, 100, 0, 70)", psd1_mult, "colz");
    TH2F *hData1 = (TH2F*)gPad->GetPrimitive("h1");
    hData1->GetXaxis()->SetTitle( "M_{STS}" );
    hData1->GetYaxis()->SetTitle("E_{PSD}^{1}, GeV");    
    hData1->SetMaximum(max);

    c1->cd(2);
    ContTree->Draw("CentralityEventContainer.GetDetectorWeight(1) : CentralityEventContainer.GetDetectorWeight(3) >> h2(100, 0, 400, 100, 0, 30)", "", "colz");
    TH2F *hData2 = (TH2F*)gPad->GetPrimitive("h2");
    hData2->GetXaxis()->SetTitle( "M_{STS}" );
    hData2->GetYaxis()->SetTitle("E_{PSD}^{2}, GeV");    
    hData2->SetMaximum(max);

    c1->cd(3);
    ContTree->Draw("CentralityEventContainer.GetDetectorWeight(2) : CentralityEventContainer.GetDetectorWeight(3) >> h3(100, 0, 400, 100, 0, 30)", "", "colz");
    TH2F *hData3 = (TH2F*)gPad->GetPrimitive("h3");
    hData3->GetXaxis()->SetTitle( "M_{STS}" );
    hData3->GetYaxis()->SetTitle("E_{PSD}^{3}, GeV");    
    hData3->SetMaximum(max);

    c1->cd(4);
    ContTree->Draw("CentralityEventContainer.GetDetectorWeight(1) : CentralityEventContainer.GetDetectorWeight(0) >> h4(100, 0, 70, 100, 0, 30)", "", "colz");
    TH2F *hData4 = (TH2F*)gPad->GetPrimitive("h4");
    hData4->GetXaxis()->SetTitle( "E_{PSD}^{1}" );
    hData4->GetYaxis()->SetTitle("E_{PSD}^{2}, GeV");    
    hData4->SetMaximum(max);
    
    
    
    return 0;     
}
