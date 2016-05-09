void DrawProfile()
{
    gStyle->SetOptStat(0000);    
    TFile *f1 = new TFile("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Macro/CutsMid.root");
    TFile *f2 = new TFile("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Macro/CutsMidLow.root");
    TFile *f3 = new TFile("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Macro/CutsMidHigh.root");
    TFile *f4 = new TFile("/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Macro/NoCutsMid.root");
    
    TCanvas *c1 = (TCanvas*)f1->Get("c3");
    TCanvas *c2 = (TCanvas*)f2->Get("c3");
    TCanvas *c3 = (TCanvas*)f3->Get("c3");
    TCanvas *c4 = (TCanvas*)f4->Get("c3");

    TProfile *pr1 = (TProfile*) c1->GetListOfPrimitives()->FindObject("profile1");
    TProfile *pr2 = (TProfile*) c1->GetListOfPrimitives()->FindObject("profile2");
    TProfile *pr3 = (TProfile*) c2->GetListOfPrimitives()->FindObject("profile1");
    TProfile *pr4 = (TProfile*) c2->GetListOfPrimitives()->FindObject("profile2");
    TProfile *pr5 = (TProfile*) c3->GetListOfPrimitives()->FindObject("profile1");
    TProfile *pr6 = (TProfile*) c3->GetListOfPrimitives()->FindObject("profile2");
    TProfile *pr7 = (TProfile*) c4->GetListOfPrimitives()->FindObject("profile1");
    gStyle->SetOptStat(0000);    
    pr1->SetName("1");
    pr2->SetName("2");
    pr3->SetName("3");
    pr4->SetName("4");
    pr5->SetName("5");
    pr6->SetName("6");
    pr7->SetName("7");
    
    pr1->SetLineWidth( 3 );
    pr1->SetLineColor( kBlack );
    pr2->SetLineColor( kBlack );
    pr3->SetLineColor( kRed );
    pr4->SetLineColor( kRed );
    pr5->SetLineColor( kBlue );
    pr6->SetLineColor( kBlue );
    pr7->SetLineColor( kGreen+3 );

    gStyle->SetOptStat(0000);    
    TCanvas *cTPC = new TCanvas ("cTPC", "cTPC", 1000, 800);
    pr1->SetTitle ("TPC");
    pr1->GetXaxis()->SetTitle( "Run" );
    pr1->GetYaxis()->SetTitle("<M_{TPC}/<M_{TPC}^{run}>>");   
    gStyle->SetOptStat(0000);    
    pr1->Draw();
    pr3->Draw("same");
    pr5->Draw("same");
    pr7->Draw("same");
    
    TLegend* legTPC = new TLegend(0.6,0.8,0.75,0.89);    
    legTPC->AddEntry(pr1, "30-70%", "l"); 
    legTPC->AddEntry(pr3, "0-70%", "l");
    legTPC->AddEntry(pr5, "30-100%", "l");
    legTPC->AddEntry(pr7, "30-70%, no track cuts", "l");
    legTPC->Draw("same");        
    
    gStyle->SetOptStat(0000);    
    TCanvas *cPSD = new TCanvas ("cPSD", "cPSD", 1000, 800);
    pr2->SetTitle ("PSD");
    pr2->GetXaxis()->SetTitle( "Run" );
    pr2->GetYaxis()->SetTitle("<E_{PSD}>/<E_{PSD}^{run}>");   

    pr2->SetMinimum(0.95);    
    pr2->SetMaximum(1.05);    
    pr2->Draw("");
    pr4->Draw("same");
    pr6->Draw("same");

    TLegend* legPSD = new TLegend(0.6,0.8,0.75,0.89);    
    legPSD->AddEntry(pr2, "30-70%", "l"); 
    legPSD->AddEntry(pr4, "0-70%", "l");
    legPSD->AddEntry(pr6, "30-100%", "l");
    legPSD->Draw("same");     



    
}