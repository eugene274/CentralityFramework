
int RunTreeInterface (Int_t FileNum=0)
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";    
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");

    TString OutFileName = CentralityFrameworkDir + "containers/" + Form("root_files/na61_container_%d.root", FileNum);
    TString InFileName = CentralityFrameworkDir + "input/" + "Merged_short.root";
    
    TreeInterface *ti = new TreeInterface();
    ti->SetOutFileName (OutFileName) ;
//     ti->SetInFileName ( Form("/lustre/nyx/cbm/users/vblinov/NA61/InputTree/%d.root", FileNum));
    ti->SetInFileName (InFileName);
    ti->WriteCentralityContainer();
    
}