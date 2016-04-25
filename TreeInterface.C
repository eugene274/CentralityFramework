
int TreeInterface (TString DataFileName = "au10au_cbm_test.root")
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";    
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");

    
    TreeInterface *ti = new TreeInterface();
    
    ti->Test();
    
}