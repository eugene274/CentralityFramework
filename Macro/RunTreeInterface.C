
int RunTreeInterface (Int_t FileNum=1)
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";    
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");

    TString OutFileName = CentralityFrameworkDir + "containers/" + Form("cbm_urqmd_CC_%d.root", FileNum);
//     TString InFileName = CentralityFrameworkDir + "input/" + "Merged.root" //Form("%d.root", FileNum);
    TString InFileName =  "/lustre/nyx/cbm/users/vblinov/CBM/mc/URQMD/Au10Au/sis100_electron/SC_ON_SL_30K/2016_06_28/datatree/Merged/Merged.root";
    
    TreeInterface *ti = new TreeInterface();
    ti->SetOutFileName (OutFileName) ;
//     ti->SetInFileName ( Form("/lustre/nyx/cbm/users/vblinov/NA61/InputTree/%d.root", FileNum));
    ti->SetInFileName (InFileName);
    ti->WriteCentralityContainer();
    
}