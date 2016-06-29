
void RunGetter()
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");  
    gStyle->SetOptStat(0000);    
    TString dir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";    
    
    TString det1 = "E_{PSD}^{1}";
    TString det2 = "E_{PSD}^{2}";
    TString det3 = "E_{PSD}^{3}";
    TString det4 = "M_{STS}";
    
//     TString ContainerFile = dir + "root_files/" + "na61_container_merged_nocuts.root";
//     TString ContainerFile = dir + "root_files/" + "na61_container_merged.root";
//     TString ContainerFile = dir + "root_files/" + "na61_container_22912.root";
//     TString ContainerFile = dir + "root_files/DCM-QGSM/" + "sis100_electron_SC_OFF_SL_OFF_urqmd.root";
    TString ContainerFile = dir + "root_files/DCM-QGSM/" + "sis100_STS_PSD_SC_OFF_SL_OFF_2016_04_11.root";
    
    TCut cuts =  "det1 > 0.65 - det2";
    
    CentralityManager *manager = new CentralityManager;
    manager->SetDirectory(dir);
     
//  For CentralityFinder 
//  ************************************   
    manager->AddDetector(det1);
    manager->AddDetector(det2);
    manager->AddDetector(det3);
    manager->AddDetector(det4);  
        
    float c = -1;
    
//     dir += "root_files/urqmd/";
//     TString SlicesFileName = dir + "Slices_M_{STS}_urqmd.root";

    dir += "root_files/dcm_qgsm_OFF_OFF/";
    TString SlicesFileName = (dir + "Slices_E_{PSD}^{1}_dcmqgsm.root");
    
    
    
    
    manager->LoadCentalityDataFile( SlicesFileName );
    
    TString DataFileName = ContainerFile;
    TFile *DataFile = new TFile ( ContainerFile, "read" );    

    CentralityEventContainer *container = new CentralityEventContainer;
    TTree *ContTree = (TTree*)DataFile->Get("cbm_data");
    ContTree->SetBranchAddress("CentralityEventContainer", &container);
    
    TH1F *hCentr = new TH1F ("hCentr", "", 101, -1, 100);

    Float_t PSD1, PSD2, PSD3, PSD12, M;
    Float_t Centrality;
    
    Int_t n = ContTree->GetEntries();
    std::cout << "n = " << n  << std::endl;

    TCanvas *c1 = new TCanvas("c1", "canvas", 900, 900);

    Int_t RunId = 0;
//     manager->SetGetterRunId (RunId);
    
    gPad->SetLogz();
    
    
    ContTree->Draw("CentralityEventContainer.GetB() : CentralityEventContainer.GetDetectorWeight(0) >> h7(200, 0, 60, 200, 0, 20)", "", "colz");
    TH2F *hData7 = (TH2F*)gPad->GetPrimitive("h7");
//     hData7->GetXaxis()->SetTitle( "E_{PSD}^{1}" );
    hData7->GetXaxis()->SetTitle( "E_{PSD}^{1}, GeV" );
    hData7->GetYaxis()->SetTitle("b, fm");    
    hData7->SetMaximum(200);
    hData7->GetXaxis()->SetTitleSize(0.05);
    hData7->GetYaxis()->SetTitleSize(0.05);   
    hData7->GetXaxis()->SetLabelSize(0.04);
    hData7->GetYaxis()->SetLabelSize(0.04);   
    
    CentralityGetter *getter = manager->GetCentralityGetter();
    
    CentralitySlice *fSlice = new CentralitySlice;
    TFile *fCentrFile = new TFile(SlicesFileName.Data(), "READ");
    TTree *fCentrTree = (TTree*) fCentrFile->Get("CentrTree");    
    fCentrTree->SetBranchAddress("CentralitySlice", &fSlice);   
    fCentrTree->GetEntry(0);        

    Int_t NSlices = fSlice->GetNSlices();
    
    for (Int_t i=1; i<NSlices; i++)
    {
        std::cout << fSlice->GetXi(i)*fSlice->GetDet1Max() << std::endl;
        gPad->Update();
        Float_t x = fSlice->GetXi(i)*fSlice->GetDet1Max();
        Double_t xr = gPad->GetX2()-gPad->GetX1();
        double x1 = (x-gPad->GetX1()) / xr;
        TLine l2;
        l2.SetLineColor(11);
        l2.DrawLineNDC(x1, 0.10, x1, 0.9);             
    }
    
    
//     manager->GetCentralityGetter()->GetGlauberB();
    
//     for (Int_t i=0; i<n; i++)
//     {
//         ContTree->GetEntry(i);
//         if (RunId != container->GetRunId())
//         {
//             RunId = container->GetRunId();
//             manager->SetGetterRunId (RunId);
//         }
//         
//         PSD1 = container->GetDetectorWeight(0);
//         PSD2 = container->GetDetectorWeight(1);
//         PSD3 = container->GetDetectorWeight(2);
//         M = container->GetDetectorWeight(3);
// 
//         Centrality = manager->GetCentrality (M, PSD1);
//         hCentr->Fill(Centrality);
//         
//     }
//     
//     hCentr->Draw();

    
}






