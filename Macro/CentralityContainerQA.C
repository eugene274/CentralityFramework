
int CentralityContainerQA (TString DataFileName)
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";
    gStyle->SetOptStat(0000);    
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");    
        
    TFile *DataFile = new TFile ( DataFileName, "read" );    

    CentralityEventContainer *container = new CentralityEventContainer;
//     TTree *ContTree = (TTree*)DataFile->Get("cbm_data");
    TTree *ContTree = (TTree*)DataFile->Get("container");
//     TTree *ContTree = (TTree*)DataFile->Get("na61_data");
    ContTree->SetBranchAddress("CentralityEventContainer", &container);
    
    Int_t n = ContTree->GetEntries();
    std::cout << "n = " << n  << std::endl;

    TCut psd1_mult = "";//" CentralityEventContainer.GetDetectorWeight(3) > 250 - 250.0/40 * CentralityEventContainer.GetDetectorWeight(0) ";
    
    TCanvas *c1 = new TCanvas("c1", "canvas", 1500, 800);
        
    c1->Divide(3,3);
    
    std::vector <TString> SHisto;
    
    SHisto.push_back ("CentralityEventContainer.GetDetectorWeight(1) : CentralityEventContainer.GetDetectorWeight(0) >> h1(450, 0, 450, 500, 0, 50)");
    SHisto.push_back ("CentralityEventContainer.GetDetectorWeight(2) : CentralityEventContainer.GetDetectorWeight(0) >> h2(450, 0, 450, 500, 0, 30)");
    SHisto.push_back ("CentralityEventContainer.GetDetectorWeight(3) : CentralityEventContainer.GetDetectorWeight(0) >> h3(450, 0, 450, 500, 0, 20)");

    SHisto.push_back ("CentralityEventContainer.GetDetectorWeight(0) + CentralityEventContainer.GetDetectorWeight(1) + CentralityEventContainer.GetDetectorWeight(2) : CentralityEventContainer.GetDetectorWeight(0) >> h4(450, 0, 450, 500, 0, 90)");

    SHisto.push_back ("CentralityEventContainer.GetB() : CentralityEventContainer.GetDetectorWeight(0) >> h5(450, 0, 450, 500, 0, 20");
    SHisto.push_back ("CentralityEventContainer.GetB() : CentralityEventContainer.GetDetectorWeight(0) >> h6(500, 0, 50, 500, 0, 20)");
    SHisto.push_back ("CentralityEventContainer.GetB() : CentralityEventContainer.GetDetectorWeight(1) >> h7(500, 0, 30, 500, 0, 20)");
    SHisto.push_back ("CentralityEventContainer.GetB() : CentralityEventContainer.GetDetectorWeight(2) >> h8(500, 0, 20, 500, 0, 20)");
    
    
    UInt_t nHisto = SHisto.size();

    for (Int_t i=0; i<nHisto; i++)
    {
        c1->cd(i+1);
        ContTree->Draw(SHisto.at(i), psd1_mult, "colz");
        gPad->SetLogz();
    }
    
//     TCanvas *c2 = new TCanvas("c2", "canvas2", 1500, 800);
//     ContTree->Draw ( "CentralityEventContainer.GetB() : CentralityEventContainer.GetDetectorWeight(3) >> h5(450, 0, 450, 500, 0, 18", "", "colz");
//     gPad->SetLogz();
    
    
    return 0;     
}
