
int RunTreeInterface (Int_t FileNumFirst=0, Int_t FileNumLast=9)
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/git/CentralityFramework/";    
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");

    TString OutFileName = CentralityFrameworkDir + "containers/" + Form("cbm_urqmd_10AGeV_%d_%d.root", FileNumFirst, FileNumLast);
//     TString InFileName = CentralityFrameworkDir + "input/" + "Merged.root" //Form("%d.root", FileNum);
//     TString InFileName =  "/lustre/nyx/cbm/users/dblau/cbm/mc/UrQMD/AuAu/10AGeV/sis100_electron/SC_ON/2016_09_01/tree/11111.root";
    
    TString InFileDir =  "/lustre/nyx/cbm/users/klochkov/mc/urqmd_auau_10gev_mbias_sis100_electron/";
    
    TChain *t = new TChain("fDataTree");
    
    for (Int_t i=FileNumFirst; i<=FileNumLast; ++i)
        t->Add( InFileDir + Form( "%d/tree.root", i ) );
    
    
    TreeInterface *ti = new TreeInterface();
    ti->SetOutFileName (OutFileName) ;
//     ti->SetInFileName (InFileName);
    ti->SetInTChain (t);
    
    ti->SetNEntries (0);
    
    ti->WriteCentralityContainer();
    
}
