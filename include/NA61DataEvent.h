//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec 21 20:36:18 2015 by ROOT version 5.34/34
// from TChain fTreeQA/
//////////////////////////////////////////////////////////

#ifndef NA61DataEvent_h
#define NA61DataEvent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class NA61DataEvent {
public :

    bool isGoodTrack(int iTrk);
    bool isRefMultTrack(int iTrk);


    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Declaration of leaf types
    Int_t           Event_Id;
    Int_t           Run_Id;
    Int_t           Event_Timestamp;
    Float_t         TPC_track_pT[2000];
    Float_t         TPC_track_eta[2000];
    Float_t         TPC_track_phi[2000];
    Int_t           TPC_track_nClusters_Total[2000];
    Int_t           TPC_track_nClusters_TPCV1[2000];
    Int_t           TPC_track_nClusters_TPCV2[2000];
    Int_t           TPC_track_nClusters_TPCVmain[2000];
    Int_t           TPC_track_nClusters_TPCVgap[2000];
    Int_t           TPC_track_nClustersPotential_Total[2000];
    Int_t           TPC_track_nClustersPotential_TPCV1[2000];
    Int_t           TPC_track_nClustersPotential_TPCV2[2000];
    Int_t           TPC_track_nClustersPotential_TPCVmain[2000];
    Int_t           TPC_track_nClustersPotential_TPCVgap[2000];
    Int_t           TPC_track_nClustersFit_Total[2000];
    Int_t           TPC_track_nClustersFit_TPCV1[2000];
    Int_t           TPC_track_nClustersFit_TPCV2[2000];
    Int_t           TPC_track_nClustersFit_TPCVmain[2000];
    Int_t           TPC_track_nClustersFit_TPCVgap[2000];
    Int_t           TPC_track_nClustersdEdX_Total[2000];
    Int_t           TPC_track_nClustersdEdX_TPCV1[2000];
    Int_t           TPC_track_nClustersdEdX_TPCV2[2000];
    Int_t           TPC_track_nClustersdEdX_TPCVmain[2000];
    Int_t           TPC_track_nClustersdEdX_TPCVgap[2000];
    Float_t         TPC_track_EnergyClusters_Total[2000];
    Float_t         TPC_track_EnergyClusters_TPCV1[2000];
    Float_t         TPC_track_EnergyClusters_TPCV2[2000];
    Float_t         TPC_track_EnergyClusters_TPCVmain[2000];
    Float_t         TPC_track_EnergyClusters_TPCVgap[2000];
    Float_t         TPC_track_DCAtoVertex_X[2000];
    Float_t         TPC_track_DCAtoVertex_Y[2000];
    Float_t         TPC_track_DCAtoVertex_Z[2000];
    Float_t         TPC_track_chi2[2000];
    Int_t           TPC_track_ndf[2000];
    Int_t           TPC_Multiplicity;
    Int_t           TPC_Multiplicity_all;
    Int_t           TPC_Multiplicity_Clusters_VTPC1_VTPC2;
    Int_t           TPC_Multiplicity_Clusters_All;
    Float_t         TPC_cos1;
    Float_t         TPC_sin1;
    Float_t         TPC_cos2;
    Float_t         TPC_sin2;
    Int_t           PSD_module_Number;
    Int_t           PSD_section_Number;
    Float_t         PSD_section_slice_Energy[10];
    Float_t         PSD_module_X[45];
    Float_t         PSD_module_Y[45];
    Float_t         PSD_module_Z[45];
    Float_t         PSD_module_Energy[45];
    Float_t         PSD_module_Energy_default[45];
    Int_t           PSD_module_number_of_sections[45];
    Float_t         PSD_section_Energy[45][10];
    Int_t           PSD_section_Number_array[45][10];
    Float_t         PSD_Energy;
    Float_t         PSD_1_Energy;
    Float_t         PSD_2_Energy;
    Float_t         PSD_3_Energy;
    Float_t         Vertex_X;
    Float_t         Vertex_Y;
    Float_t         Vertex_Z;
    Bool_t          T1;
    Bool_t          T2;
    Bool_t          T3;
    Bool_t          T4;
    Int_t           BPD_Status[2][3][3];
    Float_t         BPD_Position[2][3][3];
    Float_t         BPD_PositionError[2][3][3];
    Float_t         BPD_Z[2][3][3];
    Float_t         BPD_RMS[2][3][3];
    Float_t         BPD_Maximum[2][3][3];
    Float_t         BPD_Charge[2][3][3];
    Float_t         BPD_SumOfAll[2][3][3];
    Float_t         triggersADC[2][6];
    Bool_t          isTriggers_Simple[2][6];
    Bool_t          isTriggers_Combined[2][4];
    Float_t         Beam_Momentum[2][3];
    Float_t         Beam_Fitted2DLineXZ[2][3];
    Float_t         Beam_Fitted2DLineYZ[2][3];
    Int_t           Beam_Status[2];
    Float_t         WFA_TimeStructure[6][2000];
    Int_t           WFA_NumberOfSignalHits[6];
    Float_t         FitVertexX;
    Float_t         FitVertexY;
    Float_t         FitVertexZ;
    Int_t           FitVertexQ;

    // List of branches
    TBranch        *b_Event_Id;   //!
    TBranch        *b_Run_Id;   //!
    TBranch        *b_Event_Timestamp;   //!
    TBranch        *b_TPC_track_pT;   //!
    TBranch        *b_TPC_track_eta;   //!
    TBranch        *b_TPC_track_phi;   //!
    TBranch        *b_TPC_track_nClusters_Total;   //!
    TBranch        *b_TPC_track_nClusters_TPCV1;   //!
    TBranch        *b_TPC_track_nClusters_TPCV2;   //!
    TBranch        *b_TPC_track_nClusters_TPCVmain;   //!
    TBranch        *b_TPC_track_nClusters_TPCVgap;   //!
    TBranch        *b_TPC_track_nClustersPotential_Total;   //!
    TBranch        *b_TPC_track_nClustersPotential_TPCV1;   //!
    TBranch        *b_TPC_track_nClustersPotential_TPCV2;   //!
    TBranch        *b_TPC_track_nClustersPotential_TPCVmain;   //!
    TBranch        *b_TPC_track_nClustersPotential_TPCVgap;   //!
    TBranch        *b_TPC_track_nClustersFit_Total;   //!
    TBranch        *b_TPC_track_nClustersFit_TPCV1;   //!
    TBranch        *b_TPC_track_nClustersFit_TPCV2;   //!
    TBranch        *b_TPC_track_nClustersFit_TPCVmain;   //!
    TBranch        *b_TPC_track_nClustersFit_TPCVgap;   //!
    TBranch        *b_TPC_track_nClustersdEdX_Total;   //!
    TBranch        *b_TPC_track_nClustersdEdX_TPCV1;   //!
    TBranch        *b_TPC_track_nClustersdEdX_TPCV2;   //!
    TBranch        *b_TPC_track_nClustersdEdX_TPCVmain;   //!
    TBranch        *b_TPC_track_nClustersdEdX_TPCVgap;   //!
    TBranch        *b_TPC_track_EnergyClusters_Total;   //!
    TBranch        *b_TPC_track_EnergyClusters_TPCV1;   //!
    TBranch        *b_TPC_track_EnergyClusters_TPCV2;   //!
    TBranch        *b_TPC_track_EnergyClusters_TPCVmain;   //!
    TBranch        *b_TPC_track_EnergyClusters_TPCVgap;   //!
    TBranch        *b_TPC_track_DCAtoVertex_X;   //!
    TBranch        *b_TPC_track_DCAtoVertex_Y;   //!
    TBranch        *b_TPC_track_DCAtoVertex_Z;   //!
    TBranch        *b_TPC_track_chi2;   //!
    TBranch        *b_TPC_track_ndf;   //!
    TBranch        *b_TPC_Multiplicity;   //!
    TBranch        *b_TPC_Multiplicity_all;   //!
    TBranch        *b_TPC_Multiplicity_Clusters_VTPC1_VTPC2;   //!
    TBranch        *b_TPC_Multiplicity_Clusters_All;   //!
    TBranch        *b_TPC_cos1;   //!
    TBranch        *b_TPC_sin1;   //!
    TBranch        *b_TPC_cos2;   //!
    TBranch        *b_TPC_sin2;   //!
    TBranch        *b_PSD_module_Number;   //!
    TBranch        *b_PSD_section_Number;   //!
    TBranch        *b_PSD_section_slice_Energy;   //!
    TBranch        *b_PSD_module_X;   //!
    TBranch        *b_PSD_module_Y;   //!
    TBranch        *b_PSD_module_Z;   //!
    TBranch        *b_PSD_module_Energy;   //!
    TBranch        *b_PSD_module_Energy_default;   //!
    TBranch        *b_PSD_module_number_of_sections;   //!
    TBranch        *b_PSD_section_Energy;   //!
    TBranch        *b_PSD_section_Number_array;   //!
    TBranch        *b_PSD_Energy;   //!
    TBranch        *b_PSD_1_Energy;   //!
    TBranch        *b_PSD_2_Energy;   //!
    TBranch        *b_PSD_3_Energy;   //!
    TBranch        *b_Vertex_X;   //!
    TBranch        *b_Vertex_Y;   //!
    TBranch        *b_Vertex_Z;   //!
    TBranch        *b_T1;   //!
    TBranch        *b_T2;   //!
    TBranch        *b_T3;   //!
    TBranch        *b_T4;   //!
    TBranch        *b_BPD_Status;   //!
    TBranch        *b_BPD_Position;   //!
    TBranch        *b_BPD_PositionError;   //!
    TBranch        *b_BPD_Z;   //!
    TBranch        *b_BPD_RMS;   //!
    TBranch        *b_BPD_Maximum;   //!
    TBranch        *b_BPD_Charge;   //!
    TBranch        *b_BPD_SumOfAll;   //!
    TBranch        *b_triggersADC;   //!
    TBranch        *b_isTriggers_Simple;   //!
    TBranch        *b_isTriggers_Combined;   //!
    TBranch        *b_Beam_Momentum;   //!
    TBranch        *b_Beam_Fitted2DLineXZ;   //!
    TBranch        *b_Beam_Fitted2DLineYZ;   //!
    TBranch        *b_Beam_Status;   //!
    TBranch        *b_WFA_TimeStructure;   //!
    TBranch        *b_WFA_NumberOfSignalHits;   //!
    TBranch        *b_FitVertexX;   //!
    TBranch        *b_FitVertexY;   //!
    TBranch        *b_FitVertexZ;   //!
    TBranch        *b_FitVertexQ;   //!

    TString fOutFile;

    bool isPSD1 ;
    bool isPSD2 ;
    bool isPSD3 ;
    bool isTPC ;   

    Int_t PSD1id ;
    Int_t PSD2id ;
    Int_t PSD3id ;
    Int_t TPCid ;   

    NA61DataEvent() {};
    NA61DataEvent(TTree *tree);

    void SetData ( TString dir, Int_t RunId );


    virtual ~NA61DataEvent();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NA61DataEvent_cxx
NA61DataEvent::NA61DataEvent(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/lustre/nyx/cbm/users/vblinov/NA61/QAtree/QA_15_12_09/run-023006x001.raw_QA.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/lustre/nyx/cbm/users/vblinov/NA61/QAtree/QA_15_12_09/run-023006x001.raw_QA.root");
      }
      f->GetObject("fTreeQA",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("fTreeQA","");
      chain->Add("/lustre/nyx/cbm/users/vblinov/NA61/QAtree/QA_15_12_09/run-023001*.root");

      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

void NA61DataEvent::SetData ( TString dir, Int_t RunId )
{
    TChain * chain = new TChain("fTreeQA","");
    TString filename = dir + Form ( "run-0%d*.root", RunId );
    chain->Add(filename.Data());
    TTree *tree = chain;
    Init(tree);
    
}


NA61DataEvent::~NA61DataEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NA61DataEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NA61DataEvent::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void NA61DataEvent::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event_Id", &Event_Id, &b_Event_Id);
   fChain->SetBranchAddress("Run_Id", &Run_Id, &b_Run_Id);
   fChain->SetBranchAddress("Event_Timestamp", &Event_Timestamp, &b_Event_Timestamp);
   fChain->SetBranchAddress("TPC_track_pT", TPC_track_pT, &b_TPC_track_pT);
   fChain->SetBranchAddress("TPC_track_eta", TPC_track_eta, &b_TPC_track_eta);
   fChain->SetBranchAddress("TPC_track_phi", TPC_track_phi, &b_TPC_track_phi);
   fChain->SetBranchAddress("TPC_track_nClusters_Total", TPC_track_nClusters_Total, &b_TPC_track_nClusters_Total);
   fChain->SetBranchAddress("TPC_track_nClusters_TPCV1", TPC_track_nClusters_TPCV1, &b_TPC_track_nClusters_TPCV1);
   fChain->SetBranchAddress("TPC_track_nClusters_TPCV2", TPC_track_nClusters_TPCV2, &b_TPC_track_nClusters_TPCV2);
   fChain->SetBranchAddress("TPC_track_nClusters_TPCVmain", TPC_track_nClusters_TPCVmain, &b_TPC_track_nClusters_TPCVmain);
   fChain->SetBranchAddress("TPC_track_nClusters_TPCVgap", TPC_track_nClusters_TPCVgap, &b_TPC_track_nClusters_TPCVgap);
   fChain->SetBranchAddress("TPC_track_nClustersPotential_Total", TPC_track_nClustersPotential_Total, &b_TPC_track_nClustersPotential_Total);
   fChain->SetBranchAddress("TPC_track_nClustersPotential_TPCV1", TPC_track_nClustersPotential_TPCV1, &b_TPC_track_nClustersPotential_TPCV1);
   fChain->SetBranchAddress("TPC_track_nClustersPotential_TPCV2", TPC_track_nClustersPotential_TPCV2, &b_TPC_track_nClustersPotential_TPCV2);
   fChain->SetBranchAddress("TPC_track_nClustersPotential_TPCVmain", TPC_track_nClustersPotential_TPCVmain, &b_TPC_track_nClustersPotential_TPCVmain);
   fChain->SetBranchAddress("TPC_track_nClustersPotential_TPCVgap", TPC_track_nClustersPotential_TPCVgap, &b_TPC_track_nClustersPotential_TPCVgap);
   fChain->SetBranchAddress("TPC_track_nClustersFit_Total", TPC_track_nClustersFit_Total, &b_TPC_track_nClustersFit_Total);
   fChain->SetBranchAddress("TPC_track_nClustersFit_TPCV1", TPC_track_nClustersFit_TPCV1, &b_TPC_track_nClustersFit_TPCV1);
   fChain->SetBranchAddress("TPC_track_nClustersFit_TPCV2", TPC_track_nClustersFit_TPCV2, &b_TPC_track_nClustersFit_TPCV2);
   fChain->SetBranchAddress("TPC_track_nClustersFit_TPCVmain", TPC_track_nClustersFit_TPCVmain, &b_TPC_track_nClustersFit_TPCVmain);
   fChain->SetBranchAddress("TPC_track_nClustersFit_TPCVgap", TPC_track_nClustersFit_TPCVgap, &b_TPC_track_nClustersFit_TPCVgap);
   fChain->SetBranchAddress("TPC_track_nClustersdEdX_Total", TPC_track_nClustersdEdX_Total, &b_TPC_track_nClustersdEdX_Total);
   fChain->SetBranchAddress("TPC_track_nClustersdEdX_TPCV1", TPC_track_nClustersdEdX_TPCV1, &b_TPC_track_nClustersdEdX_TPCV1);
   fChain->SetBranchAddress("TPC_track_nClustersdEdX_TPCV2", TPC_track_nClustersdEdX_TPCV2, &b_TPC_track_nClustersdEdX_TPCV2);
   fChain->SetBranchAddress("TPC_track_nClustersdEdX_TPCVmain", TPC_track_nClustersdEdX_TPCVmain, &b_TPC_track_nClustersdEdX_TPCVmain);
   fChain->SetBranchAddress("TPC_track_nClustersdEdX_TPCVgap", TPC_track_nClustersdEdX_TPCVgap, &b_TPC_track_nClustersdEdX_TPCVgap);
   fChain->SetBranchAddress("TPC_track_EnergyClusters_Total", TPC_track_EnergyClusters_Total, &b_TPC_track_EnergyClusters_Total);
   fChain->SetBranchAddress("TPC_track_EnergyClusters_TPCV1", TPC_track_EnergyClusters_TPCV1, &b_TPC_track_EnergyClusters_TPCV1);
   fChain->SetBranchAddress("TPC_track_EnergyClusters_TPCV2", TPC_track_EnergyClusters_TPCV2, &b_TPC_track_EnergyClusters_TPCV2);
   fChain->SetBranchAddress("TPC_track_EnergyClusters_TPCVmain", TPC_track_EnergyClusters_TPCVmain, &b_TPC_track_EnergyClusters_TPCVmain);
   fChain->SetBranchAddress("TPC_track_EnergyClusters_TPCVgap", TPC_track_EnergyClusters_TPCVgap, &b_TPC_track_EnergyClusters_TPCVgap);
   fChain->SetBranchAddress("TPC_track_DCAtoVertex_X", TPC_track_DCAtoVertex_X, &b_TPC_track_DCAtoVertex_X);
   fChain->SetBranchAddress("TPC_track_DCAtoVertex_Y", TPC_track_DCAtoVertex_Y, &b_TPC_track_DCAtoVertex_Y);
   fChain->SetBranchAddress("TPC_track_DCAtoVertex_Z", TPC_track_DCAtoVertex_Z, &b_TPC_track_DCAtoVertex_Z);
   fChain->SetBranchAddress("TPC_track_chi2", TPC_track_chi2, &b_TPC_track_chi2);
   fChain->SetBranchAddress("TPC_track_ndf", TPC_track_ndf, &b_TPC_track_ndf);
   fChain->SetBranchAddress("TPC_Multiplicity", &TPC_Multiplicity, &b_TPC_Multiplicity);
   fChain->SetBranchAddress("TPC_Multiplicity_all", &TPC_Multiplicity_all, &b_TPC_Multiplicity_all);
   fChain->SetBranchAddress("TPC_Multiplicity_Clusters_VTPC1_VTPC2", &TPC_Multiplicity_Clusters_VTPC1_VTPC2, &b_TPC_Multiplicity_Clusters_VTPC1_VTPC2);
   fChain->SetBranchAddress("TPC_Multiplicity_Clusters_All", &TPC_Multiplicity_Clusters_All, &b_TPC_Multiplicity_Clusters_All);
   fChain->SetBranchAddress("TPC_cos1", &TPC_cos1, &b_TPC_cos1);
   fChain->SetBranchAddress("TPC_sin1", &TPC_sin1, &b_TPC_sin1);
   fChain->SetBranchAddress("TPC_cos2", &TPC_cos2, &b_TPC_cos2);
   fChain->SetBranchAddress("TPC_sin2", &TPC_sin2, &b_TPC_sin2);
   fChain->SetBranchAddress("PSD_module_Number", &PSD_module_Number, &b_PSD_module_Number);
   fChain->SetBranchAddress("PSD_section_Number", &PSD_section_Number, &b_PSD_section_Number);
   fChain->SetBranchAddress("PSD_section_slice_Energy", PSD_section_slice_Energy, &b_PSD_section_slice_Energy);
   fChain->SetBranchAddress("PSD_module_X", PSD_module_X, &b_PSD_module_X);
   fChain->SetBranchAddress("PSD_module_Y", PSD_module_Y, &b_PSD_module_Y);
   fChain->SetBranchAddress("PSD_module_Z", PSD_module_Z, &b_PSD_module_Z);
   fChain->SetBranchAddress("PSD_module_Energy", PSD_module_Energy, &b_PSD_module_Energy);
   fChain->SetBranchAddress("PSD_module_Energy_default", PSD_module_Energy_default, &b_PSD_module_Energy_default);
   fChain->SetBranchAddress("PSD_module_number_of_sections", PSD_module_number_of_sections, &b_PSD_module_number_of_sections);
   fChain->SetBranchAddress("PSD_section_Energy", PSD_section_Energy, &b_PSD_section_Energy);
   fChain->SetBranchAddress("PSD_section_Number_array", PSD_section_Number_array, &b_PSD_section_Number_array);
   fChain->SetBranchAddress("PSD_Energy", &PSD_Energy, &b_PSD_Energy);
   fChain->SetBranchAddress("PSD_1_Energy", &PSD_1_Energy, &b_PSD_1_Energy);
   fChain->SetBranchAddress("PSD_2_Energy", &PSD_2_Energy, &b_PSD_2_Energy);
   fChain->SetBranchAddress("PSD_3_Energy", &PSD_3_Energy, &b_PSD_3_Energy);
   fChain->SetBranchAddress("Vertex_X", &Vertex_X, &b_Vertex_X);
   fChain->SetBranchAddress("Vertex_Y", &Vertex_Y, &b_Vertex_Y);
   fChain->SetBranchAddress("Vertex_Z", &Vertex_Z, &b_Vertex_Z);
   fChain->SetBranchAddress("T1", &T1, &b_T1);
   fChain->SetBranchAddress("T2", &T2, &b_T2);
   fChain->SetBranchAddress("T3", &T3, &b_T3);
   fChain->SetBranchAddress("T4", &T4, &b_T4);
   fChain->SetBranchAddress("BPD_Status", BPD_Status, &b_BPD_Status);
   fChain->SetBranchAddress("BPD_Position", BPD_Position, &b_BPD_Position);
   fChain->SetBranchAddress("BPD_PositionError", BPD_PositionError, &b_BPD_PositionError);
   fChain->SetBranchAddress("BPD_Z", BPD_Z, &b_BPD_Z);
   fChain->SetBranchAddress("BPD_RMS", BPD_RMS, &b_BPD_RMS);
   fChain->SetBranchAddress("BPD_Maximum", BPD_Maximum, &b_BPD_Maximum);
   fChain->SetBranchAddress("BPD_Charge", BPD_Charge, &b_BPD_Charge);
   fChain->SetBranchAddress("BPD_SumOfAll", BPD_SumOfAll, &b_BPD_SumOfAll);
   fChain->SetBranchAddress("triggersADC", triggersADC, &b_triggersADC);
   fChain->SetBranchAddress("isTriggers_Simple", isTriggers_Simple, &b_isTriggers_Simple);
   fChain->SetBranchAddress("isTriggers_Combined", isTriggers_Combined, &b_isTriggers_Combined);
   fChain->SetBranchAddress("Beam_Momentum", Beam_Momentum, &b_Beam_Momentum);
   fChain->SetBranchAddress("Beam_Fitted2DLineXZ", Beam_Fitted2DLineXZ, &b_Beam_Fitted2DLineXZ);
   fChain->SetBranchAddress("Beam_Fitted2DLineYZ", Beam_Fitted2DLineYZ, &b_Beam_Fitted2DLineYZ);
   fChain->SetBranchAddress("Beam_Status", Beam_Status, &b_Beam_Status);
   fChain->SetBranchAddress("WFA_TimeStructure", WFA_TimeStructure, &b_WFA_TimeStructure);
   fChain->SetBranchAddress("WFA_NumberOfSignalHits", WFA_NumberOfSignalHits, &b_WFA_NumberOfSignalHits);
   fChain->SetBranchAddress("FitVertexX", &FitVertexX, &b_FitVertexX);
   fChain->SetBranchAddress("FitVertexY", &FitVertexY, &b_FitVertexY);
   fChain->SetBranchAddress("FitVertexZ", &FitVertexZ, &b_FitVertexZ);
   fChain->SetBranchAddress("FitVertexQ", &FitVertexQ, &b_FitVertexQ);
   Notify();
}

Bool_t NA61DataEvent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NA61DataEvent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NA61DataEvent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
    entry++;
    return 1;
}
#endif // #ifdef NA61DataEvent_cxx
