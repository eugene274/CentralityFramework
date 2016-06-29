#include <vector> 

#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TGraphErrors.h"


void SlicesQA()
{
    TString CentralityFrameworkDir = "/lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/";
    gStyle->SetOptStat(0000);    
    gSystem->Load( CentralityFrameworkDir + "build/libCentrality");    
    
    TString FileDir = CentralityFrameworkDir + "root_files/";
    
    
    TString dir = "urqmd/";
    std::vector <TString> DataFileNameVec;
    DataFileNameVec.push_back (dir + "Slices_M_{STS}_urqmd.root");
    DataFileNameVec.push_back (dir + "Slices_M_{STS}_E_{PSD}^{1}_urqmd.root");
    DataFileNameVec.push_back (dir + "Slices_E_{PSD}^{1}_urqmd.root");

    dir = "dcm_qgsm_OFF_OFF/";
    DataFileNameVec.push_back (dir + "Slices_M_{STS}_dcmqgsm.root");    
    DataFileNameVec.push_back (dir + "Slices_M_{STS}_E_{PSD}^{1}_dcmqgsm.root");
    DataFileNameVec.push_back (dir + "Slices_E_{PSD}^{1}_dcmqgsm.root");

    Int_t colors[] = {1,2,4,1,2,4};
    Int_t markers[] = {24,25,26,20,21,22};
    
    
    UInt_t NFiles = DataFileNameVec.size();

    std::vector <TFile*> DataFileVec (NFiles, NULL);
    std::vector <TTree*> DataTreeVec (NFiles, NULL);
    std::vector <CentralitySlice*> SlicesVec (NFiles, NULL);
    std::vector <TGraphErrors*> GraphVec (NFiles, NULL);
    std::vector <TGraphErrors*> GraphVecSb (NFiles, NULL);

    for (UInt_t iFile=0; iFile<NFiles; iFile++)
    {
        DataFileVec.at(iFile) = ( new TFile ( FileDir + DataFileNameVec.at(iFile), "read" ));
        SlicesVec.at(iFile) = (new CentralitySlice);
        DataTreeVec.at(iFile) = (TTree*)DataFileVec.at(iFile)->Get("CentrTree");
        DataTreeVec.at(iFile) -> SetBranchAddress("CentralitySlice", &(SlicesVec.at(iFile)) );
        DataTreeVec.at(iFile) -> GetEntry(0);
    }
    

    
    std::vector <Float_t> centrality;
    std::vector <Float_t> sigmaB;
    std::vector <Float_t> dsigmaB;
    
    
    Float_t step = 5;//SlicesVec.at(0)->GetSlicesStep ();
    UInt_t NSlices = 20;//SlicesVec.at(0)->GetNSlices();
    
//     NSlices --;
    for (UInt_t i=0; i<NSlices; i++)
        centrality.push_back( (i+0.5)*step );

    TCanvas *c1 = new TCanvas("c1", "canvas", 900, 900);
//     c1->Divide(2,1);
    
//     c1->cd(1);
    gPad->SetLogy();
    for (UInt_t iFile=0; iFile<NFiles; iFile++){
        
        sigmaB.clear();
        dsigmaB.clear();
        
        NSlices = SlicesVec.at(iFile)->GetNSlices() - 1;
        if (NSlices==18) NSlices++;
        
        for (UInt_t j=0; j<=NSlices; j++){
            float temp_dB = SlicesVec.at(iFile)->GetdB().at(j);
            float temp_B = SlicesVec.at(iFile)->GetMeanB().at(j);
            float temp_sB = SlicesVec.at(iFile)->GetSigmaB().at(j);
            float temp_dsB = SlicesVec.at(iFile)->GetdSigmaB().at(j);
            
            sigmaB.push_back(temp_sB/temp_B);
            
            float ddd = TMath::Sqrt ( (temp_dsB/temp_B)*(temp_dsB/temp_B) + (temp_dsB/temp_B/temp_B*temp_dB)*(temp_dsB/temp_B/temp_B*temp_dB) );
            dsigmaB.push_back(ddd);
        }

        GraphVecSb.at(iFile) = new TGraphErrors(NSlices+1, &(centrality[0]), &(sigmaB[0]), 0, &(dsigmaB[0]) );
        
//         GraphVecSb.at(iFile) -> SetMarkerSize(2);    
        GraphVecSb.at(iFile) -> SetMarkerStyle(markers[iFile]);    
        GraphVecSb.at(iFile) -> SetMarkerColor(colors[iFile]);    
        GraphVecSb.at(iFile) -> SetLineColor(colors[iFile]);    
        
        
        if (iFile == 0){
            GraphVecSb.at(iFile) -> SetTitle ("Impact parameter vs centrality");
            GraphVecSb.at(iFile) -> GetXaxis()->SetTitle( "Centrality, %" );
            GraphVecSb.at(iFile) -> GetYaxis()->SetTitle("#sigma_{b}/<b>");    

            GraphVecSb.at(iFile) -> GetXaxis()->SetLimits(0.0, (NSlices+1.5)*step);
            GraphVecSb.at(iFile) -> GetYaxis()->SetRangeUser(0.03, 0.6);
            
            GraphVecSb.at(iFile) -> GetXaxis()->SetTitleSize(0.05);
            GraphVecSb.at(iFile) -> GetYaxis()->SetTitleSize(0.05);   
            GraphVecSb.at(iFile) -> GetXaxis()->SetLabelSize(0.04);
            GraphVecSb.at(iFile) -> GetYaxis()->SetLabelSize(0.04);   
            
            GraphVecSb.at(iFile) -> Draw("AL");      
        }
        else{
            if (iFile > 2)
                GraphVecSb.at(iFile) -> Draw("P");      
            else 
                GraphVecSb.at(iFile) -> Draw("L");      

            
        }
    }
    
    
    
    TLegend* leg2 = new TLegend(0.6,0.6,0.85,0.85);
    
    leg2->SetNColumns(2);
    
    leg2->AddEntry(GraphVecSb.at(0), "    ",  "l");            
    leg2->AddEntry(GraphVecSb.at(3), "  M_{STS}", "p");            
    leg2->AddEntry(GraphVecSb.at(2), "    ",  "l");        
    leg2->AddEntry(GraphVecSb.at(5), "  E_{PSD}^{1}", "p");              
    leg2->AddEntry(GraphVecSb.at(1), "    ", "l");                                                 
    leg2->AddEntry(GraphVecSb.at(4), "  M_{STS}+E_{PSD}^{1}", "p");   
    
    
    
    leg2->Draw("same");     
    
    
//     c1->cd(2);
//     for (UInt_t iFile=0; iFile<DataFileNameVec.size(); iFile++){
//         
//         NSlices = SlicesVec.at(iFile)->GetNSlices() - 1;
//         if (NSlices==18) NSlices++;
// 
//         GraphVec.at(iFile) =  ( new TGraphErrors(NSlices+1, &(centrality[0]), &(SlicesVec.at(iFile)->GetMeanB()[0]), 0, &(SlicesVec.at(iFile)->GetSigmaB()[0])) );
//         
// //         GraphVec.at(iFile) -> SetMarkerSize(2);    
//         GraphVec.at(iFile) -> SetMarkerStyle(markers[iFile]);    
//         GraphVec.at(iFile) -> SetMarkerColor(colors[iFile]);    
//         GraphVec.at(iFile) -> SetLineColor(colors[iFile]);    
//         
//         
//         if (iFile == 0){
//             GraphVec.at(iFile) -> SetTitle ("Impact parameter vs centrality");
//             GraphVec.at(iFile) -> GetXaxis()->SetTitle( "Centrality, %" );
//             GraphVec.at(iFile) -> GetYaxis()->SetTitle("b, fm");    
// 
//             GraphVec.at(iFile) -> GetXaxis()->SetLimits(0.0, (NSlices+1.5)*step);
//             GraphVec.at(iFile) -> GetYaxis()->SetRangeUser(0.0, 17);
//             
//             GraphVec.at(iFile) -> GetXaxis()->SetTitleSize(0.05);
//             GraphVec.at(iFile) -> GetYaxis()->SetTitleSize(0.05);   
//             GraphVec.at(iFile) -> GetXaxis()->SetLabelSize(0.05);
//             GraphVec.at(iFile) -> GetYaxis()->SetLabelSize(0.05);   
//             
//             GraphVec.at(iFile) -> Draw("APL");      
//         }
//         else{
//             GraphVec.at(iFile) -> Draw("PL");      
//         }
//     }
    
    
}  
    
/*
    TString DataFileName = "Slices_Mult_urqmd_0.root";     //"Slices_Mult_urqmd.root";
    TFile *DataFile = new TFile ( DataFileName, "read" );    
    CentralitySlice *slices = new CentralitySlice;
    TTree *Tree = (TTree*)DataFile->Get("CentrTree");
    Tree->SetBranchAddress("CentralitySlice", &slices);
    Tree->GetEntry(0);

    TString DataFileName1 = "Slices_PSD1_urqmd_0.root";    //"Slices_PSD1_urqmd.root";
    TFile *DataFile1 = new TFile ( DataFileName1, "read" );    
    CentralitySlice *slices1 = new CentralitySlice;
    TTree *Tree1 = (TTree*)DataFile1->Get("CentrTree");
    Tree1->SetBranchAddress("CentralitySlice", &slices1);
    Tree1->GetEntry(0);
    
    TString DataFileName2 = "Slices_Mult_urqmd.root";     //"Slices_PSD1_urqmd.root";
    TFile *DataFile2 = new TFile ( DataFileName2, "read" );    
    CentralitySlice *slices2 = new CentralitySlice;
    TTree *Tree2 = (TTree*)DataFile2->Get("CentrTree");
    Tree2->SetBranchAddress("CentralitySlice", &slices2);
    Tree2->GetEntry(0);

    TString DataFileName3 = "Slices_PSD1_urqmd.root";     //"Slices_PSD1_urqmd.root";
    TFile *DataFile3 = new TFile ( DataFileName3, "read" );    
    CentralitySlice *slices3 = new CentralitySlice;
    TTree *Tree3 = (TTree*)DataFile3->Get("CentrTree");
    Tree3->SetBranchAddress("CentralitySlice", &slices3);
    Tree3->GetEntry(0);



    
    std::vector <Float_t> centrality;
    Float_t step = slices->GetSlicesStep ();
    UInt_t NSlices = slices1->GetNSlices();
    
    NSlices --;
    for (UInt_t i=0; i<=NSlices; i++)
        centrality.push_back( (i+0.5)*step );

    
    if (true)
    {
        TCanvas *c1 = new TCanvas("c1", "canvas", 1000, 800);

        TGraphErrors *grB = new TGraphErrors(NSlices+1, &(centrality[0]), &(slices->GetMeanB()[0]), 0, &(slices->GetSigmaB()[0]));
        grB->SetTitle ("Impact parameter vs centrality");
        grB->SetMarkerStyle(23);    
        grB->SetMarkerSize(2);    
        grB->GetXaxis()->SetTitle( "Centrality, %" );
        grB->GetYaxis()->SetTitle("b, fm");    

        grB->GetXaxis()->SetLimits(0.0, (NSlices+1.5)*step);
        grB->GetYaxis()->SetRangeUser(0.0, 17);
        grB->Draw("APL");
        
        grB->GetXaxis()->SetTitleSize(0.05);
        grB->GetYaxis()->SetTitleSize(0.05);   
        grB->GetXaxis()->SetLabelSize(0.05);
        grB->GetYaxis()->SetLabelSize(0.05);   
    }
    
    
    if (true)
    {
        TCanvas *c2 = new TCanvas("c2", "canvas", 1000, 800);
        c2->SetLogy();
        
        std::vector <Float_t> SigmaB;
        for (UInt_t i=0; i<=NSlices; i++)
            SigmaB.push_back( slices->GetSigmaB()[i]/slices->GetMeanB()[i] );
        
        TGraphErrors *grSigmaB = new TGraph (NSlices+1, &(centrality[0]), &(SigmaB[0]));
        grSigmaB->SetTitle ("#sigma_{b}/b vs centrality");
        grSigmaB->SetMarkerStyle(23);    
        grSigmaB->SetMarkerSize(2);    
        grSigmaB->SetMarkerColor(kBlack);    
        grSigmaB->GetXaxis()->SetTitle( "Centrality, %" );
        grSigmaB->GetYaxis()->SetTitle("#sigma_{b}/b");    

        grSigmaB->GetXaxis()->SetLimits(0.0, (NSlices+1.5)*step);
        grSigmaB->GetYaxis()->SetRangeUser(0.03, 1);
        grSigmaB->Draw("APL");
        
        grSigmaB->GetXaxis()->SetTitleSize(0.05);
        grSigmaB->GetYaxis()->SetTitleSize(0.05);   
        grSigmaB->GetXaxis()->SetLabelSize(0.05);
        grSigmaB->GetYaxis()->SetLabelSize(0.05);  
        
        std::vector <Float_t> SigmaB1;
        for (UInt_t i=0; i<=NSlices; i++)
            SigmaB1.push_back( slices1->GetSigmaB()[i]/slices1->GetMeanB()[i] );
        
        TGraphErrors *grSigmaB1 = new TGraph (NSlices+1, &(centrality[0]), &(SigmaB1[0]));
        grSigmaB1->SetMarkerStyle(22);    
        grSigmaB1->SetMarkerSize(2);    
        grSigmaB1->SetMarkerColor(kRed);    
        grSigmaB1->SetLineColor(kRed);    
        grSigmaB1->Draw("PLsame");


        std::vector <Float_t> SigmaB2;
        for (UInt_t i=0; i<=NSlices; i++)
            SigmaB2.push_back( slices2->GetSigmaB()[i]/slices2->GetMeanB()[i] );
        
        TGraphErrors *grSigmaB2 = new TGraph (NSlices+1, &(centrality[0]), &(SigmaB2[0]));
        grSigmaB2->SetMarkerStyle(23);    
        grSigmaB2->SetMarkerSize(2);    
        grSigmaB2->SetMarkerColor(kBlue);    
        grSigmaB2->SetLineColor(kBlue);    
        grSigmaB2->Draw("PLsame");

        std::vector <Float_t> SigmaB3;
        for (UInt_t i=0; i<=NSlices; i++)
            SigmaB3.push_back( slices3->GetSigmaB()[i]/slices3->GetMeanB()[i] );
        
        TGraphErrors *grSigmaB3 = new TGraph (NSlices+1, &(centrality[0]), &(SigmaB3[0]));
        grSigmaB3->SetMarkerStyle(22);    
        grSigmaB3->SetMarkerSize(2);    
        grSigmaB3->SetMarkerColor(kMagenta);    
        grSigmaB3->SetLineColor(kMagenta);    
        grSigmaB3->Draw("PLsame");
        
        TLegend* leg1 = new TLegend(0.7,0.75,0.89,0.89);
        leg1->AddEntry(grSigmaB ,"M_{STS}, target out","p");    
        leg1->AddEntry(grSigmaB1 ,"E_{PSD}^{1}, target out","p");    
        leg1->AddEntry(grSigmaB2 ,"M_{STS}, target in","p");    
        leg1->AddEntry(grSigmaB3 ,"E_{PSD}^{1}, target in","p");    
        
        
        
//         leg1->AddEntry(grSigmaB2 ,"M_{STS}+E_{PSD}^{1}","p");    
        leg1->Draw("same");
    }   */
