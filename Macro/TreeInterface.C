
int TreeInterface (Int_t RunId=23005)
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";    
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");

    
    NA61DataEventNew *ti = new NA61DataEventNew();
    ti->SetOutFileName ( CentralityFrameworkDir + Form("root_files/na61_container_%d.root", RunId) );
    ti->SetInFileName ( Form("/lustre/nyx/cbm/users/vblinov/NA61/InputTree/%d.root", RunId));
    ti->WriteCentralityContainer();
    
}