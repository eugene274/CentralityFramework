// #include "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/include/CentralityManager.h"
TString sCentralityCodePath = "/lustre/nyx/cbm/users/dblau/CentralityFramework/";

void TestGetter( TString InFileName = "/lustre/nyx/cbm/users/dblau/cbm/mc/UrQMD/AuAu/10AGeV/sis100_electron/SC_ON/2016_09_01/tree/11111.root", 
                 TString det1 = "M_{STS}", TString det2 = "" ){

    TString CentralityFrameworkDir = sCentralityCodePath;    
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");


    CentralityManager* fCentralityManager;
    fCentralityManager = new CentralityManager();
    fCentralityManager -> SetDirectory(sCentralityCodePath);

    TString dir = sCentralityCodePath + "root_files/" + Form("Slices_%s_%s.root", det1, det2);

    //    fCentralityManager -> LoadCentalityDataFile( sCentralityCodePath +
    //"root_files/" + Form("Slices_%s_%s.root", "M_{STS}", "") );
    fCentralityManager -> LoadCentalityDataFile(dir);
    fCentralityManager -> Det1IsInt(true);                     

    DataTreeEvent *fEvent = new DataTreeEvent;
    
    TFile *fInFile = new TFile (InFileName);
    TTree *fInTree = (TTree*) fInFile->Get("fDataTree"); 
    fInTree->SetBranchAddress("DTEvent", &fEvent);
    Int_t nEntries = fInTree->GetEntries();
    
    TH1F *hTest = new TH1F("hTest", "hTest", 100, 0, 100 );
    
    for(Int_t i=0; i</*nEntries*/10000; i++)
    {
        fInTree->GetEntry(i);
        if ((i+1) % 50 == 0) cout << "Event # " << i+1 << "... \r" << flush; 
        Int_t RefMult = fEvent->GetNTracks();
        double c = fCentralityManager->GetCentrality(RefMult);
//         cout << "Centrality: " << c << endl;
        hTest->Fill(c);
    }
    
    hTest->Draw("E");
    
    
//     cout << dir << endl;
//     double c=fCentralityManager->GetCentrality(100);
//     cout<<"Centrality: "<<c<<endl;

}
