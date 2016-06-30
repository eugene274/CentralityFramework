
int RunTreeInterface (Int_t FileNum=1)
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";    
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");

    TString OutFileName = CentralityFrameworkDir + "containers/" + Form("cbm_urqmd_CC_%d.root", FileNum);
    TString InFileName = CentralityFrameworkDir + "input/" + Form("%d.root", FileNum);
    
    TreeInterface *ti = new TreeInterface();
    ti->SetOutFileName (OutFileName) ;
//     ti->SetInFileName ( Form("/lustre/nyx/cbm/users/vblinov/NA61/InputTree/%d.root", FileNum));
    ti->SetInFileName (InFileName);
    ti->WriteCentralityContainer();
    
}