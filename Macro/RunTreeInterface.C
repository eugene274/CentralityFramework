
int RunTreeInterface (Int_t FileNum=1)
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";    
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");

    TString OutFileName = CentralityFrameworkDir + "containers/" + Form("cbm_urqmd_10AGeV_%d.root", FileNum);
//     TString InFileName = CentralityFrameworkDir + "input/" + "Merged.root" //Form("%d.root", FileNum);
    TString InFileName =  "/lustre/nyx/cbm/users/dblau/cbm/mc/UrQMD/AuAu/10AGeV/sis100_electron/SC_ON/2016_09_01/tree/11111.root";
    
    TreeInterface *ti = new TreeInterface();
    ti->SetOutFileName (OutFileName) ;
    ti->SetInFileName (InFileName);
    
    ti->SetNEntries (0);
    
    ti->WriteCentralityContainer();
    
}