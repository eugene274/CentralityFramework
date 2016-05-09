#include "DataTreeMaker.h"
#include <iostream>
#include <fstream>
#include <vector>
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;

#include "TDirectory.h"
#include "CbmPsdHit.h"
#include "CbmPsdDigi.h"
#include "CbmMCTrack.h"
#include "CbmStsTrack.h"
#include "CbmTrackMatchNew.h"

#include "FairMCEventHeader.h"
#include "UEvent.h"
// #include "QdetVector.h"
// #include "QdetVectorGetter.h"

#include "CbmL1PFFitter.h"
#include "CbmKFVertex.h"
#include "L1Field.h"

#include "DataTreeEvent.h"
#include "DataTreeTrack.h"

//=================================================================> MAIN <===============================================================
DataTreeMaker::DataTreeMaker()
  : FairTask("DataTreeMaker",1),
  uEvent (new UEvent()),
  fTreeEvents (new TChain("events","events"))/*,
  DTEvent (new DataTreeEvent())*/
  
//     l1 (new CbmL1())
{
      DTEvent = new DataTreeEvent();
}
DataTreeMaker::~DataTreeMaker()
{
    
}
//=================================================================> INIT <===============================================================
InitStatus DataTreeMaker::Init()
{
    fCurEvent = 0;
    FairRootManager* ioman = FairRootManager::Instance();
    fPrimVtx = (CbmVertex*) ioman->GetObject("PrimaryVertex");
    fHeader = (FairMCEventHeader*) ioman->GetObject("MCEventHeader.");
    //BEGIN temp for bad phi data
    fTreeEvents -> Add(sInputFileName);
    fTreeEvents -> SetBranchAddress("event",&uEvent);
    //END temp for bad phi data
//     fHeader = (CbmMCEventHeader*) ioman->GetObject("MCEventHeader.");
    flistPSDhit = (TClonesArray*) ioman->GetObject("PsdHit");
    flistPSDdigit = (TClonesArray*) ioman->GetObject("PsdDigi");
    flistSTSMCtrack = (TClonesArray*) ioman->GetObject("MCTrack");
    flistSTSRECOtrack = (TClonesArray*) ioman->GetObject("StsTrack");
    flistSTStrackMATCH = (TClonesArray*) ioman->GetObject("StsTrackMatch");
    
    DataTreeEvent_Init();
    OutputTree_Init();
//     PSD_Init();
}
//--------------------------------------------------------------------------------------------------
void DataTreeMaker::DataTreeEvent_Init()
{
    for (int i=0;i<nPSD_Modules;i++)
    {
	DTEvent -> AddPSDModule(10);
    }
}
//--------------------------------------------------------------------------------------------------
void DataTreeMaker::OutputTree_Init()
{
    fTreeFile = new TFile(sOutputFileName, "RECREATE");
    fTreeFile -> cd();
    fTreeQA = new TTree("fTreeQA","fTreeQA");
    fTreeQA -> SetMaxTreeSize(9000000);
    
    fTreeQA -> Branch("DTEvent", &DTEvent, 256000, 3);    
}
//--------------------------------------------------------------------------------------------------
void DataTreeMaker::PSD_Init()
{
    ifstream fxypos(sPSDCoordinatesFileName);
    
    for (Int_t i=0; i<nPSD_Modules; i++) // to change if different number of PSD modules: use fxypos.end() -> continue? then inc for the other loops on module?
    {
	if (fxypos.eof()) break;
	fxypos >> PSD_module_X[i] >> PSD_module_Y[i];
	cout << "PSD module " << i << ": " << PSD_module_X[i] << "; " << PSD_module_Y[i] << endl;
    }
    fxypos.close();
    cout << PSD_module_X[0] << endl;
}

//=================================================================> EXEC <===============================================================
void DataTreeMaker::Exec(Option_t* opt)
{
    cout << "Event " << ++fCurEvent << endl;
    Clear_Event();
    
    Read_Event();
    Read_PSD();
    Read_STS();
    DTEvent -> Process();
    fTreeQA -> Fill();
}
//--------------------------------------------------------------------------------------------------
void DataTreeMaker::Clear_Event()
{
    DTEvent -> ClearEvent();
}
//--------------------------------------------------------------------------------------------------
void DataTreeMaker::Read_Event()
{
    if (!fHeader) cout << "No fHeader!" << endl;
    else
    {
	DTEvent -> SetRPAngle(fHeader -> GetRotZ());
	DTEvent -> SetImpactParameter(fHeader -> GetB());
    }
    //BEGIN temp for bad phi data
    fTreeEvents->GetEntry(fCurEvent-1);
    DTEvent -> SetRPAngle(uEvent->GetPhi());//used only for DCM-QGSM
//     cout << "RPAngle = " << DTEvent -> GetRPAngle() << endl;
    //END temp for bad phi data
}
//--------------------------------------------------------------------------------------------------
void DataTreeMaker::Read_PSD()
{
    ifstream fxypos(sPSDCoordinatesFileName);
    
    double fX = 0.;
    double fY = 0.;
    
    for (int i=0; i<nPSD_Modules; i++) // to change if different number of PSD modules: use fxypos.end() -> continue? then inc for the other loops on module?
    {
	if (fxypos.eof()) break;
	fxypos >> fX >> fY;
	DTEvent -> GetPSDModule(i) -> SetPosition(fX, fY, nUndefinedValue);
// 	cout << "[" << i << "]: x = " << DTEvent -> GetPSDModule(i) -> GetPositionComponent(0) << "; y = " << DTEvent -> GetPSDModule(i) -> GetPositionComponent(1) << endl;
    }
    fxypos.close();
    
    int nPSDdigits = flistPSDdigit->GetEntriesFast();
    CbmPsdDigi* digit = NULL;
    for (int i=0; i<nPSDdigits; i++)
    {
	digit  = (CbmPsdDigi*) flistPSDdigit -> At(i);
	if (!digit) continue;
	DTEvent -> GetPSDModule(digit->GetModuleID()-1) -> GetSection(digit->GetSectionID()-1) -> AddEnergy(digit->GetEdep());
    }
}
//--------------------------------------------------------------------------------------------------
void DataTreeMaker::Read_STS()
{
    CbmStsTrack* track;
    CbmTrackMatchNew* match;
    CbmMCTrack* mctrack;
	
    Double_t p, px, py, pz, energy, mass, pT, phi, eta, MC_px, MC_py, MC_pz, MC_p, MC_pT, MC_phi, MC_eta;
    Int_t type, trackID;
    
    const FairTrackParam *trackParam;
    TVector3 momRec;

    Int_t nSTStracks = flistSTSRECOtrack->GetEntries();
//     cout << "nSTStracks = " << nSTStracks << endl;
    //cout << "evenPlane::STSRECOtransverseMomMeth: # STS reco tracks = " << nSTStracks << endl; 
    
    // Extrapolation track parameters back to primary vertex
    vector<CbmStsTrack> vRTracks;
    vector<CbmStsTrack> vRTracks_old;
    vRTracks.resize(nSTStracks);   
      
    
    
    //BEGIN fitter
    CbmL1PFFitter fitter;
    vector<float> vChiToPrimVtx;
    CbmKFVertex kfVertex;
    if(fPrimVtx)
    kfVertex = CbmKFVertex(*fPrimVtx);
  
    vector<L1FieldRegion> vField;

    vector<int> Pdg;
    Pdg.resize(nSTStracks); 
    for (Int_t i=0; i<nSTStracks; i++)
    {
	DTEvent->AddTrack();
	vRTracks[i] = *( (CbmStsTrack*) flistSTSRECOtrack->At(i));
	
	if(!track)
	{
	    cout << "ERROR: empty track!";
	    continue;
	}
	track = &vRTracks[i];
	trackParam = track->GetParamFirst();	
	trackParam->Momentum(momRec);
	    
	px = momRec.X();
	py = momRec.Y();
	pz = momRec.Z();
	pT = TMath::Sqrt(px*px+py*py);
	p = TMath::Sqrt(px*px+py*py+pz*pz);
	phi = TMath::ATan2(py,px);//+(TMath::Pi()/2);
	if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
	if (phi>TMath::Pi()) phi-=2*TMath::Pi();
	eta = TMath::Log((p+pz)/(p-pz))/2.;
	
	DTEvent->GetTrack(i)->SetPt(0,pT);
	DTEvent->GetTrack(i)->SetPhi(0,phi);
	DTEvent->GetTrack(i)->SetEta(0,eta);
	
	DTEvent->GetTrack(i)->SetNofHits(0,track->GetNofHits());
	DTEvent->GetTrack(i)->SetFlag(0,track->GetFlag());
	DTEvent->GetTrack(i)->SetChiSq(0,track->GetChiSq());
	DTEvent->GetTrack(i)->SetNDF(0,track->GetNDF());
		
	Pdg[i] = 211;
    }
    
//     double Chi2_test[nSTStracks];
//     int NDF_test[nSTStracks];
//     for (Int_t i=0; i<nSTStracks; i++)
//     {
// 	track = &vRTracks[i];
// 	Chi2_test[i] = track->GetChiSq();
// 	NDF_test[i] = track->GetNDF();
//     }
    
//     cout << "goto fitter" << endl;
    fitter.Fit(vRTracks, Pdg);
//     fitter.Fit(vRTracks, pdg);
    fitter.GetChiToVertex(vRTracks, vField, vChiToPrimVtx, kfVertex, 1.e9f);//tracks array, field, dca over error
    //END fitter
    
    for (Int_t i=0; i<nSTStracks; i++)
    {
	track = &vRTracks[i];
	if(!track)
	{
	    cout << "ERROR: empty track!";
	    continue;
	}
	
// 	match = (CbmTrackMatchNew*) flistSTStrackMATCH->At(i);
// 	trackID = match->GetMCTrackId();
// 	cout << "trackID = " << trackID << endl;
// 	cout << trackID << endl;
// 	if (trackID >= 0)
// 	{
// 	    mctrack = (CbmMCTrack*) flistSTSMCtrack->At(trackID);
// 	    if (!mctrack)
// 	    {
// 		cout << "No mc track found!" << endl;
// 	    }
// 	    else
// 	    {
// // 		type = mctrack->GetPdgCode();
// // 		mass = mctrack->GetMass();	
// 		MC_px = mctrack->GetPx();
// 		MC_py = mctrack->GetPy();
// 		MC_pz = mctrack->GetPz();
// 		MC_p = TMath::Sqrt(MC_px*MC_px + MC_py*MC_py + MC_pz*MC_pz);
// 		MC_pT = TMath::Sqrt(MC_px*MC_px + MC_py*MC_py);
// 		MC_phi = TMath::ATan2(MC_py,MC_px);
// 		if (MC_phi<-TMath::Pi()) MC_phi+=2*TMath::Pi();
// 		if (MC_phi>TMath::Pi()) MC_phi-=2*TMath::Pi();
// 		MC_eta = TMath::Log((MC_p+MC_pz)/(MC_p-MC_pz))/2.;
// 		Sim_TPC_track_pT[i] = MC_pT;
// 		Sim_TPC_track_eta[i] = MC_eta;
// 		Sim_TPC_track_phi[i] = MC_phi;
// 	    }
// 	}
	
	trackParam = track->GetParamFirst();	
	trackParam->Momentum(momRec);
	    
	px = momRec.X();
	py = momRec.Y();
	pz = momRec.Z();
	pT = TMath::Sqrt(px*px+py*py);
	p = TMath::Sqrt(px*px+py*py+pz*pz);
	phi = TMath::ATan2(py,px);//+(TMath::Pi()/2);
	if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
	if (phi>TMath::Pi()) phi-=2*TMath::Pi();
	eta = TMath::Log((p+pz)/(p-pz))/2.;
	
	DTEvent->GetTrack(i)->SetPt(1,pT);
	DTEvent->GetTrack(i)->SetPhi(1,phi);
	DTEvent->GetTrack(i)->SetEta(1,eta);
	
	DTEvent->GetTrack(i)->SetNofHits(1,track->GetNofHits());
	DTEvent->GetTrack(i)->SetFlag(1,track->GetFlag());
	DTEvent->GetTrack(i)->SetChiSq(1,track->GetChiSq());
	DTEvent->GetTrack(i)->SetNDF(1,track->GetNDF());
	
// 	cout << pT << " " << NofHits << " " << Chi2 << endl;
// 	cout << i << " Chi2: after " << track->GetChiSq() << " and before " << Chi2_test[i] << "; Chi2/NDF after " << track->GetChiSq()/track->GetNDF() << " and before " << Chi2_test[i]/NDF_test[i] << endl;
    }
}

//================================================================> FINISH <==============================================================
void DataTreeMaker::Finish()
{
    cout << "DataTreeMaker::Finish" << endl;
    fTreeQA -> Write();
    fTreeFile -> Write();
    fTreeFile -> Close();
}

ClassImp(DataTreeMaker)





// //====================================================================================
// //====================================================================================
// void DataTreeMaker::Get_STS_Track_Data()
//     {
// 	CbmStsTrack* track;
// 	CbmTrackMatch* match;
// 	CbmMCTrack* mctrack;
// 	
// 	Double_t p, px, py, pz, energy, mass, pT, phi, eta, MC_px, MC_py, MC_pz, MC_p, MC_pT, MC_phi, MC_eta;
// 	Int_t type, trackID;
// 	
// 	const FairTrackParam *trackParam;
// 	TVector3 momRec;
// 
// 	Int_t nSTStracks = flistSTSRECOtrack->GetEntriesFast();
// 	//cout << "evenPlane::STSRECOtransverseMomMeth: # STS reco tracks = " << nSTStracks << endl; 
// 	
// 	// Extrapolation track parameters back to primary vertex
// 	vector<CbmStsTrack> vRTracks;
// 	vRTracks.resize(nSTStracks);    
// 	  
// 	fSTS_Qx_1 = 0;
// 	fSTS_Qy_1 = 0;
// 	fSTS_Qx_2 = 0;
// 	fSTS_Qy_2 = 0;
// 	Float_t fSTS_Multiplicity = 0.;
// 	
// 	//BEGIN fitter
// 	CbmL1PFFitter fitter;
// 	vector<int> Pdg;
// 	Pdg.resize(nSTStracks); 
// 	for (Int_t i=0; i<nSTStracks; i++)
// 	{
// 	    vRTracks[i] = *( (CbmStsTrack*) flistSTSRECOtrack->At(i));
// 	    Pdg[i] = 0;
// 	}
// 	cout << "goto fitter" << endl;
// 	fitter.Fit(vRTracks, Pdg);
// 	//END fitter
// 	
// 	for (Int_t i=0; i<nSTStracks; i++)
// 	{
// 	    vRTracks[i] = *( (CbmStsTrack*) flistSTSRECOtrack->At(i));  
// 	    track = &vRTracks[i];
// 	    
// 	    match = (CbmTrackMatch*) flistSTStrackMATCH->At(i);
// 	    trackID = match->GetMCTrackId();
// 	    
// 	    //if (trackID < 0) continue; // ghost track from combinatorix (no corresponding MC tracks)
// 	    if (trackID >= 0)
// 	    {
// 		mctrack = (CbmMCTrack*) flistSTSMCtrack->At(trackID); 
// 		//cout << "trackID = " << trackID << endl;
// 		if (!mctrack)
// 		{
// 		    cout << "No mc track found!" << endl;
// 		}
// 		else
// 		{
// 		    type = mctrack->GetPdgCode();
// 		    mass = mctrack->GetMass();	
// 		    MC_px = mctrack->GetPx();
// 		    MC_py = mctrack->GetPy();
// 		    MC_pz = mctrack->GetPz();
// 		    MC_p = TMath::Sqrt(MC_px*MC_px + MC_py*MC_py + MC_pz*MC_pz);
// 		    MC_pT = TMath::Sqrt(MC_px*MC_px + MC_py*MC_py);
// 		    MC_phi = TMath::ATan2(MC_py,MC_px);
// 		    if (MC_phi<-TMath::Pi()) MC_phi+=2*TMath::Pi();
// 		    if (MC_phi>TMath::Pi()) MC_phi-=2*TMath::Pi();
// 		    MC_eta = TMath::Log((MC_p+MC_pz)/(MC_p-MC_pz))/2.;
// 		}
// 	    }
// 	    trackParam = track->GetParamFirst();	
// 	    trackParam->Momentum(momRec);
// 	    
// 	    px = momRec.X();
// 	    py = momRec.Y();
// 	    pz = momRec.Z();
// 	    pT = TMath::Sqrt(px*px+py*py);
// 	    p = TMath::Sqrt(px*px+py*py+pz*pz);
// 	    phi = TMath::ATan2(py,px);//+(TMath::Pi()/2);
// 	    if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
// 	    if (phi>TMath::Pi()) phi-=2*TMath::Pi();
// 	    eta = TMath::Log((p+pz)/(p-pz))/2.;
// 	    
// 	    fSTS_Qx_1 += pT*TMath::Cos(phi);
// 	    fSTS_Qy_1 += pT*TMath::Sin(phi);
// 	    fSTS_Qx_2 += pT*TMath::Cos(2*phi);
// 	    fSTS_Qy_2 += pT*TMath::Sin(2*phi);
// 	    
// 	    //CUTS
// // 	    if (TMath::Abs(MC_eta - eta) > 0.01) continue;
// //	    if (TMath::Abs(MC_pT - pT) > 0.01) continue;
// // 	    if (TMath::Abs(MC_phi - phi) > 0.01) continue;
// 	    //	    if ((phi < -3.14/2+0.5 and phi > -3.14/2-0.5) or (phi > 3.14/2-0.5 and phi < 3.14/2+0.5)) continue;
// 	    
// 	    //cout << "p = " << p << "; pz = " << pz << "; underlog = " << (p+pz)/(p-pz) << "; afterlog = " << TMath::Log((p+pz)/(p-pz)) << "; alterlog = " << log((p+pz)/(p-pz)) << "; eta = " << eta << endl;
// 	    //cout << "STS Track №" << i << ": px = " << px << "; py = " << py << "; pT = " << pT << "; Phi = " << phi << endl;
// 	    if (fTempdetVectorGetter->Apply_STS_Cuts(pT, eta))
// 	    {
// 		hSTS_Phi->Fill(phi);
// 		hSTS_pT->Fill(pT);
// 		hSTS_Eta->Fill(eta);
// 		hSTS_pT_vs_Eta->Fill(pT, eta);
// 		hSTS_pT_vs_Phi->Fill(pT, phi);
// 		hSTS_Eta_vs_Phi->Fill(eta, phi);
// 		hSTS_NofMvdHits->Fill(track->GetNofMvdHits());
// 		hSTS_pT_MC_vs_pT_RECO->Fill(MC_pT,pT);
// 		hSTS_phi_MC_vs_phi_RECO->Fill(MC_phi,phi);
// 		hSTS_eta_MC_vs_eta_RECO->Fill(MC_eta,eta);
// 		hSTS_diff_eta_MC_eta_RECO->Fill(MC_eta-eta);
// 		fSTS_Multiplicity++;
// 	    }
// 	    hSTS_nocuts_Phi->Fill(phi);
// 	    hSTS_nocuts_pT->Fill(pT);
// 	    hSTS_nocuts_Eta->Fill(eta);
// 	    hSTS_nocuts_pT_vs_Eta->Fill(pT, eta);
// 	    hSTS_nocuts_pT_vs_Phi->Fill(pT, phi);
// 	    hSTS_nocuts_Eta_vs_Phi->Fill(eta, phi);
// 	    hSTS_nocuts_NofMvdHits->Fill(track->GetNofMvdHits());
// 	    hSTS_nocuts_pT_MC_vs_pT_RECO->Fill(MC_pT,pT);
// 	    hSTS_nocuts_phi_MC_vs_phi_RECO->Fill(MC_phi,phi);
// 	    hSTS_nocuts_eta_MC_vs_eta_RECO->Fill(MC_eta,eta);
// 	    hSTS_nocuts_diff_eta_MC_eta_RECO->Fill(MC_eta-eta);
// 	
// //	    uVector.SetXYZM(px, py, pz, pT);
//             //cout << "uVector: px, py, pz, pT = " << px << " " << py << " " << pz << " " << pT << endl;
// //	    QdetVector_STS_plus->AddChannel(i, uVector, 0, 0, 0);
// //	    type = track->GetPdgCode();
// 
// // 	    if (type > 1000000000)
// // 	    {
// // 		mass = TMath::Floor( ( type - 1000000000 - 10000 * ( TMath::Floor( ( type - 1000000000 ) / 10000 ) ) ) / 10 );
// // 		p = track->GetP();
// // 		px = track -> GetPx();
// // 		py = track -> GetPy();
// // 		pz = track -> GetPz();
// // 		energy = TMath::Sqrt(mass*0.93827203*mass*0.93827203 + p*p);
// // 	    }
//     
// // 	    //uVector.SetPxPyPzE(px, py, pz, energy);
// //	    QdetVector *fTempDetVector = (QdetVector*) fDetVector.ConstructedAt(fDetVectorIdx);
// // 	    fTempDetVector->AdduVector(uVector, 1);
// 	}
// 	fEvent_Multiplicity_STS = nSTStracks;
// 	fTree_Multiplicity_STS = nSTStracks;
// 	hSTS_nocuts_Multiplicity->Fill(fEvent_Multiplicity_STS);	
// 	hSTS_Multiplicity->Fill(fSTS_Multiplicity);
//     }    
// //====================================================================================
// //====================================================================================
// //====================================================================================
// //====================================================================================
// void DataTreeMaker::Get_MC_STS_Track_Data()
//     {
// 	CbmMCTrack* mctrack;
// 	
// 	Double_t p, px, py, pz, energy, mass, pT, phi, eta;
// 	Int_t type, trackID;
// 	
// 	const FairTrackParam *trackParam;
// 	TVector3 momRec;
// 
// 	Int_t nMC_STStracks = flistSTSMCtrack->GetEntriesFast();
// 	//cout << "evenPlane::STSRECOtransverseMomMeth: # STS reco tracks = " << nMC_STStracks << endl; 
// 	
// 	// Extrapolation track parameters back to primary vertex
// 	vector<CbmMCTrack> vRTracks;
// 	vRTracks.resize(nMC_STStracks);
// 	Float_t fSTS_Multiplicity = 0.; 
// 	
// 	fEvent_Multiplicity_MC_STS = 0;
// 	for (Int_t i=0; i<nMC_STStracks; i++)
// 	{
// 	    vRTracks[i] = *( (CbmMCTrack*) flistSTSMCtrack->At(i));  
// 	    mctrack = &vRTracks[i];
// 	    
// 	    //============================================================
// 	    //cuts
// 	    if ( ! mctrack ) 
// 	    {
// 		cout << "evenPlane::STSMCtransverseMomMeth: no track pointer!" << endl;
// 		continue;
// 	    }
// 
// 	    // acceptance
// 	    if ( mctrack->AccSTS() > 0 ) hMC_STS_AccSTS->Fill(mctrack->AccSTS());
// 	    if ( mctrack->AccSTS() < 9 ) continue;
// 	    //cuts
// 	    //============================================================
// 	    fEvent_Multiplicity_MC_STS++;
// 	    
// 	    
// 	    px = mctrack->GetPx();
// 	    py = mctrack->GetPy();
// 	    pz = mctrack->GetPz();
// 	    pT = TMath::Sqrt(px*px+py*py);
// 	    p = TMath::Sqrt(px*px+py*py+pz*pz);
// 	    phi = TMath::ATan2(py,px);//+(TMath::Pi()/2);
// 	    if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
// 	    if (phi>TMath::Pi()) phi-=2*TMath::Pi();
// 	    eta = TMath::Log((p+pz)/(p-pz))/2.;
// 	    //cout << "p = " << p << "; pz = " << pz << "; underlog = " << (p+pz)/(p-pz) << "; afterlog = " << TMath::Log((p+pz)/(p-pz)) << "; alterlog = " << log((p+pz)/(p-pz)) << "; eta = " << eta << endl;
// 	    //cout << "STS Track №" << i << ": px = " << px << "; py = " << py << "; pT = " << pT << "; Phi = " << phi << endl;
// 	    if (fTempdetVectorGetter->Apply_STS_Cuts(pT, eta))
// 	    {
// 		hMC_STS_Phi->Fill(phi);
// 		hMC_STS_pT->Fill(pT);
// 		hMC_STS_Eta->Fill(eta);
// 		hMC_STS_pT_vs_Eta->Fill(pT, eta);
// 		hMC_STS_pT_vs_Phi->Fill(pT, phi);
// 		hMC_STS_Eta_vs_Phi->Fill(eta, phi);
// 		fSTS_Multiplicity++;
// 	    }
// 	    hMC_STS_nocuts_Phi->Fill(phi);
// 	    hMC_STS_nocuts_pT->Fill(pT);
// 	    hMC_STS_nocuts_Eta->Fill(eta);
// 	    hMC_STS_nocuts_pT_vs_Eta->Fill(pT, eta);
// 	    hMC_STS_nocuts_pT_vs_Phi->Fill(pT, phi);
// 	    hMC_STS_nocuts_Eta_vs_Phi->Fill(eta, phi);		
// 	    
// 	    
// // 	    type = mctrack->GetPdgCode();
// // 	    mass = mctrack->GetMass();	
// // 	    if (type < 1000000000)
// // 	    {
// // 		charge = track->GetCharge();
// // 		y = track->GetRapidity();
// // 	    }
// // 	    else
// // 	    {
// // 		//pdg = 1000000000 + 10*1000*z + 10*a + i;
// // 		charge = TMath::Floor( ( type - 1000000000 ) / 10000 );
// // 		mass = TMath::Floor( ( type - 1000000000 -  10000 * charge ) / 10 );
// // 		p = track->GetP();
// // 		energy = TMath::Sqrt(mass*0.93827203*mass*0.93827203 + p*p);
// // 		pz = track->GetPz();
// // 		y = 0.5 * TMath::Log( ( energy + pz ) / ( energy - pz ) );      // "pz" convention
// // 	    }
// // 
// // 	    if (charge == 0.) continue; // safety condition
// // 	    
// // 	    if (type > 1000000000)
// // 	    {
// // 		mass = TMath::Floor( ( type - 1000000000 - 10000 * ( TMath::Floor( ( type - 1000000000 ) / 10000 ) ) ) / 10 );
// // 		p = track->GetP();
// // 		px = track -> GetPx();
// // 		py = track -> GetPy();
// // 		pz = track -> GetPz();
// // 		energy = TMath::Sqrt(mass*0.93827203*mass*0.93827203 + p*p);
// // 	    }
// 	    mass = mctrack->GetMass();	
// 	    energy = TMath::Sqrt(mass*0.93827203*mass*0.93827203 + px*px + py*py + pz*pz);
// 	    uVector.SetXYZM(px, py, pz, energy);
// 	    QdetVector_MC_STS_plus->AddChannel(i, uVector, mass, 0, 0);
//     	}
// 	
// 	hMC_STS_nocuts_Multiplicity->Fill(fEvent_Multiplicity_MC_STS);
// 	hMC_STS_Multiplicity->Fill(fSTS_Multiplicity);
//     }    
// //====================================================================================
// //====================================================================================
// //====================================================================================
// //====================================================================================
// void DataTreeMaker::Fill_QA_Histograms()
//     {
// 	hSTS_Multiplicity_vs_PSD_Energy->Fill(fEvent_Multiplicity_STS, fEvent_Energy_PSD);
// 	hSTS_Multiplicity_vs_PSD_1_Energy->Fill(fEvent_Multiplicity_STS/200, fEvent_Energy_PSD_1/80);
// 	hSTS_Multiplicity_vs_PSD_2_Energy->Fill(fEvent_Multiplicity_STS/200, fEvent_Energy_PSD_2/20);
// 	hSTS_Multiplicity_vs_PSD_3_Energy->Fill(fEvent_Multiplicity_STS/200, fEvent_Energy_PSD_3/10);
// // 	hPSD_v1_E->Divide(hPSD_v_E_Entries);
// // 	hPSD_v2_E->Divide(hPSD_v_E_Entries);
// // 	hPSD_v1_R_E->Divide(hPSD_v_R_E_Entries);
// // 	hPSD_v2_R_E->Divide(hPSD_v_R_E_Entries);
// // 	hSTS_v1_M->Divide(hSTS_v_M_Entries);
// // 	hSTS_v2_M->Divide(hSTS_v_M_Entries);
// // 	hSTS_v1_pT_M->Divide(hSTS_v_pT_M_Entries);
// // 	hSTS_v2_pT_M->Divide(hSTS_v_pT_M_Entries);
// 	
// // 	for (Int_t i=0; i<100; i++)
// // 	{
// // 	    if (hPSD_v_E_Entries->GetBinContent(i) > 0)
// // 	    {
// // 		hPSD_v1_E -> SetBinContent(i,(hPSD_v1_E -> GetBinContent(i))/(hPSD_v_E_Entries->GetBinContent(i)));
// // 		hPSD_v2_E -> SetBinContent(i,(hPSD_v2_E -> GetBinContent(i))/(hPSD_v_E_Entries->GetBinContent(i)));
// // 		hSTS_v1_M -> SetBinContent(i,(hSTS_v1_M -> GetBinContent(i))/(hSTS_v_M_Entries->GetBinContent(i)));
// // 		hSTS_v2_M -> SetBinContent(i,(hSTS_v2_M -> GetBinContent(i))/(hSTS_v_M_Entries->GetBinContent(i)));
// // 	    }
// // 	    else
// // 	    {
// // 		hPSD_v1_E -> SetBinContent(i,0);
// // 	    }
// // 	    for (Int_t j=0; j<100; j++)
// // 	    {
// // 		if (hSTS_v_pT_M_Entries->GetBinContent(i,j) > 0)
// // 		{
// // 		    hPSD_v1_R_E -> SetBinContent(i,(hPSD_v2_R_E -> GetBinContent(i,j))/(hPSD_v_R_E_Entries->GetBinContent(i,j)));
// // 		    hPSD_v2_R_E -> SetBinContent(i,(hPSD_v2_R_E -> GetBinContent(i,j))/(hPSD_v_R_E_Entries->GetBinContent(i,j)));
// // 		    hSTS_v1_pT_M -> SetBinContent(i,(hSTS_v1_pT_M -> GetBinContent(i,j))/(hPSD_v_R_E_Entries->GetBinContent(i,j)));
// // 		    hSTS_v2_pT_M -> SetBinContent(i,(hSTS_v2_pT_M -> GetBinContent(i,j))/(hPSD_v_R_E_Entries->GetBinContent(i,j)));
// // 		}
// // 		else
// // 		{
// // 		    hPSD_v2_E -> SetBinContent(i,0);
// // 		}
// // 	    }
// // 	}
// // 	for (Int_t i=0; i<100; i++)
// // 	{
// // 	    if (hSTS_v_M_Entries->GetBinContent(i) > 0)
// // 	    {
// // 		hSTS_v1_M -> SetBinContent(i,(hSTS_v1_M -> GetBinContent(i))/(hSTS_v_M_Entries->GetBinContent(i)));
// // 	    }
// // 	    else
// // 	    {
// // 		hSTS_v1_M -> SetBinContent(i,0);
// // 	    }
// // 	    if (hSTS_v_M_Entries->GetBinContent(i) > 0)
// // 	    {
// // 		hSTS_v2_M -> SetBinContent(i,(hSTS_v2_M -> GetBinContent(i))/(hSTS_v_M_Entries->GetBinContent(i)));
// // 	    }
// // 	    else
// // 	    {
// // 		hSTS_v2_M -> SetBinContent(i,0);
// // 	    }
// // 	    
// // 	}
//     }
// //====================================================================================
// //====================================================================================
// //====================================================================================
// //====================================================================================
// void DataTreeMaker::Init_Tree()
//     {
// // 	fTree->Branch("QdetVector_PSD_plus", &QdetVector_PSD_plus, "QdetVector_PSD_plus");
// 	fTree->Branch("QdetVector_PSD_plus", "QdetVector", &QdetVector_PSD_plus);
// 	fTree->Branch("QdetVector_PSD_minus", "QdetVector", &QdetVector_PSD_minus);
// 	fTree->Branch("QdetVector_STS_plus", "QdetVector", &QdetVector_STS_plus);
// 	fTree->Branch("QdetVector_STS_minus", "QdetVector", &QdetVector_STS_minus);
// 	fTree->Branch("QdetVector_MC_STS_plus", "QdetVector", &QdetVector_MC_STS_plus);
// 	fTree->Branch("QdetVector_MC_STS_minus", "QdetVector", &QdetVector_MC_STS_minus);
// // 	fTree->Branch("", "QdetVector", &);
// 	fTree->Branch("fTree_Energy_PSD", &fTree_Energy_PSD, "fTree_Energy_PSD/D");
// 	fTree->Branch("fTree_Multiplicity_STS", &fTree_Multiplicity_STS, "fTree_Multiplicity_STS/D");
// 	fTree->Branch("fTree_Multiplicity_MC_STS", &fTree_Multiplicity_MC_STS, "fTree_Multiplicity_MC_STS/D");
// 	fTree->Branch("fTree_Event_Id", &fTree_Event_Id, "fTree_Event_Id/I"); 
//     }
// //====================================================================================
// //====================================================================================
// //====================================================================================
// //====================================================================================
// void DataTreeMaker::Fill_Tree()
//     {
// 	fTree->Fill();
//     }
// //====================================================================================
// //====================================================================================
// //====================================================================================
// //====================================================================================
// void DataTreeMaker::Write_Tree()
//     {
//         TDirectory *curr_dir = gDirectory->CurrentDirectory();
//         fOutputFile->cd();
// 	fTree->Write();
//         curr_dir->cd();
//     }
// //====================================================================================
// //====================================================================================
// //====================================================================================
// //====================================================================================
// void DataTreeMaker::Check_Axis_Ranges()
//     {
// // 	if (hPSD_Energy->GetEntries() != hPSD_Energy->GetEffectiveEntries())
// // 	{
// // 	  cout << "hPSD_Energy has the wrong X axis edges. EffectiveEntries = " << hPSD_Energy->GetEffectiveEntries() << "/" << hPSD_Energy->GetEntries() << endl;
// // 	  for (Int_t i = 1; i<10000; i++)
// // 	  {
// // 	    hPSD_Energy->SetAxisRange(0.,i);
// // 	    if (hPSD_Energy->GetEntries() == hPSD_Energy->GetEffectiveEntries())
// // 	    {
// // 	      cout << i << " is a good range! EffectiveEntries = " << hPSD_Energy->GetEffectiveEntries() << "/" << hPSD_Energy->GetEntries() << endl;	      
// // 	      break; 
// // 	    }
// // 	    cout << i << " is not a good range. EffectiveEntries = " << hPSD_Energy->GetEffectiveEntries() << "/" << hPSD_Energy->GetEntries() << endl;	
// // 	  }
// // 	}
//     }
// //====================================================================================
// //====================================================================================
// //====================================================================================
// //====================================================================================
// void DataTreeMaker::Write_QA_Histograms()
//     {
//       // TO DO: Average Phi with weight in PSD as a function of total energy in the whole amount of events.
//       // TO DO: Energy distribution in the whole amount of events.
//       
//       // at Finish stage
// // 	for (Int_t i=0; i<fDetVectorIdx; i++)
// // 	{
// // 	    //QdetVector *fTempDetVector = (QdetVector*) fDetVector.ConstructedAt(fDetVectorIdx);
// // 	    
// // 	    for (Int_t j=0; j<fTempDetVector->GetNumberOfuVectors(); j++)
// // 	    {
// // 		TLorentzVector *fTempLorentzVector = (TLorentzVector*) fTempDetVector->GetuVector(j);
// // 		Int_t TempuVectorId = fTempDetVector->GetuVectorId(j);
// // 		if (TempuVectorId == 0) hPSD_Energy->Fill(fTempLorentzVector->Px(), fTempLorentzVector->Py(), fTempLorentzVector->E());
// // 		else hSTS_Energy->Fill(fTempLorentzVector->Px(), fTempLorentzVector->Py(), fTempLorentzVector->E());
// // 	      
// // 	    }
// // 	}
// 	cout << "fEventIdx = " << fEventIdx << endl;
// // 	for (Int_t i=1; i<finc_mod; i++)
// // 	{
// // 	  fedep_mod_av[i] = fedep_mod_av[i]/fEventIdx;
// // 	  hPSD_mod_Ave_Energy->Fill(fedep_mod_av[i]);
// // 	}
// // 	hPSD_mod_Ave_Energy->Write();
// // 	
// 	TClonesArray &hHistos = *hPSD_mod_Energy;
// 	for (Int_t i=1; i<finc_mod; i++)
// 	{
// 	  TH1F *hTempPSD_mod_Energy = (TH1F*)hHistos.At(i);
// 	  List_hPSD_mod_Energy->Add(hTempPSD_mod_Energy);
// 	  //hTempPSD_mod_Energy->Write();
// 	  //cout << "Module " << i << "; Number of sections = " << finc_sec[i] << endl;
// 	  for (Int_t j=1; j<finc_sec[i]; j++)
// 	  {
// 	    //cout << "Module " << i << "; Section = " << j << "; Energy = " << fedep_sec_av[i][j] << endl;
// 	    hPSD_ModSecEnergy->Fill(i-1, j, fedep_sec_av[i][j]);
// 	  }
// 	}
// // 	List_hPSD_mod_Energy->Write("hPSD_mod_Energy", TObject::kSingleKey);
// // 	hPSD_Energy->Write();
// // 	hSTS_Multiplicity->Write();
// // 	hPSD_Phi->Write();
// // 	hSTS_Phi->Write();
// // 	hSTS_Eta->Write();
// // 	hSTS_pT->Write();
// // 	hPSD_Phi_Energy->Write();
// // 	hPSD_cPhi->Write();
// // 	hPSD_cPhi_Energy->Write();
// // 	hPSD_XYEnergy->Write();
// // 	hPSD_ModSecEnergy->Write();
// // 	hSTS_Multiplicity_vs_PSD_Energy->Write();
// // 	prfx_hSTS_Multiplicity_vs_PSD_Energy=hSTS_Multiplicity_vs_PSD_Energy->ProfileX("prfx_hSTS_Multiplicity_vs_PSD_Energy",1,100);
// // 	hPSD_geometry->Write();
// 	prfx_hSTS_Multiplicity_vs_PSD_Energy=hSTS_Multiplicity_vs_PSD_Energy->ProfileX("prfx_hSTS_Multiplicity_vs_PSD_Energy",1,100);
// 
// 	Int_t fMaxBin = hSTS_Multiplicity->GetMaximumBin();
// 	cout << "Max bin in STS MUltiplicity = " << fMaxBin << endl;
// 	
// // 	List_hFlow->Add(hFlowQx[0]);
// // 	List_hFlow->Add(hFlowQx[1]);
// // 	List_hFlow->Add(hFlowQx[2]);
// // // 	List_hFlow->Add(hFlowQy[0]);
// // // 	List_hFlow->Add(hFlowQy[1]);
// // // 	List_hFlow->Add(hFlowQy[2]);
// // 	List_hFlow->Add(hFlowCosQx);
// // 	List_hFlow->Add(hFlowQxQx[0]);
// // 	List_hFlow->Add(hFlowQxQx[1]);
// // 	List_hFlow->Add(hFlowQxQx[2]);
// // 	hFlowCosQx->Write();
// // 	hFlowQxQx[0]->Write();
// // 	hFlowQxQx[1]->Write();
// // 	hFlowQxQx[2]->Write();
// // 	
// //  	for (Int_t s=0;s<2;s++)
// // 	{
// // 	    for (Int_t v=0;v<5;v++)
// // 	    {
// // 		for (Int_t a=0;a<5;a++)
// // 		{
// // 		    List_hFlowQ->Add(hFlow_Qx_n1[v][a][s]);
// // 		    List_hFlowQ->Add(hFlow_Qy_n1[v][a][s]);
// // 		    List_hFlowQ->Add(hFlow_Qx_n2[v][a][s]);
// // 		    List_hFlowQ->Add(hFlow_Qy_n2[v][a][s]);
// // 		    for (Int_t b=0;b<5;b++)
// // 		    {
// // 			List_hFlowQQ->Add(hFlow_QxQx_n1[v][a][b][s]);
// // 			List_hFlowQQ->Add(hFlow_QyQy_n1[v][a][b][s]);
// // 			List_hFlowQQ->Add(hFlow_QyQx_n1[v][a][b][s]);
// // 			List_hFlowQQ->Add(hFlow_QxQx_n2[v][a][b][s]);
// // 			List_hFlowQQ->Add(hFlow_QyQy_n2[v][a][b][s]);
// // 			List_hFlowQQ->Add(hFlow_QyQx_n2[v][a][b][s]);
// // 		    }
// // 		}
// // 	    }
// // 	}
// 	
// 	
// 	
// 	
//         
// 	List_hPSD->Add(hPSD_Energy);
// 	List_hPSD->Add(hPSD_1_Energy);
// 	List_hPSD->Add(hPSD_2_Energy);
// 	List_hPSD->Add(hPSD_3_Energy);
// 	List_hPSD->Add(hPSD_1_Energy_vs_PSD_2_Energy);
// 	List_hPSD->Add(hPSD_Phi);
// 	List_hPSD->Add(hPSD_Phi_Energy);
// 	List_hPSD->Add(hPSD_XYEnergy);
// 	List_hPSD->Add(hPSD_X_Energy);
// 	List_hPSD->Add(hPSD_Y_Energy);
// 	List_hPSD->Add(hPSD_ModSecEnergy);
// 	List_hPSD->Add(hPSD_geometry);
// 	List_hPSD->Add(List_hPSD_mod_Energy);
// 	
// 	
// 	List_hPSD->Add(hPSD_v1_R_E);
// 	List_hPSD->Add(hPSD_v2_R_E);
// 	List_hPSD->Add(hPSD_v_R_E_Entries);
// 	List_hPSD->Add(hPSD_v1_E);
// 	List_hPSD->Add(hPSD_v2_E);
// 	List_hPSD->Add(hPSD_v_E_Entries);
// 	List_hPSD->Add(hPSD_EPAngle_1);
// 	List_hPSD->Add(hPSD_EPAngle_2);
// 	
// 	List_hSTS->Add(hSTS_Multiplicity);
// 	List_hSTS->Add(hSTS_Phi);
// 	List_hSTS->Add(hSTS_Eta);
// 	List_hSTS->Add(hSTS_pT);
// 	List_hSTS->Add(hSTS_pT_vs_Eta);
// 	List_hSTS->Add(hSTS_pT_vs_Phi);
// 	List_hSTS->Add(hSTS_Eta_vs_Phi);
// 	List_hSTS->Add(hSTS_NofMvdHits);
// 	
// 	List_hSTS->Add(hSTS_nocuts_Multiplicity);
// 	List_hSTS->Add(hSTS_nocuts_Phi);
// 	List_hSTS->Add(hSTS_nocuts_Eta);
// 	List_hSTS->Add(hSTS_nocuts_pT);
// 	List_hSTS->Add(hSTS_nocuts_pT_vs_Eta);
// 	List_hSTS->Add(hSTS_nocuts_pT_vs_Phi);
// 	List_hSTS->Add(hSTS_nocuts_Eta_vs_Phi);
// 	List_hSTS->Add(hSTS_nocuts_NofMvdHits);
// 	
// 	List_hSTS->Add(hSTS_v1_pT_M);
// 	List_hSTS->Add(hSTS_v2_pT_M);
// 	List_hSTS->Add(hSTS_v_pT_M_Entries);
// 	List_hSTS->Add(hSTS_v1_M);
// 	List_hSTS->Add(hSTS_v2_M);
// 	List_hSTS->Add(hSTS_v_M_Entries);
// 	List_hSTS->Add(hSTS_EPAngle_1);
// 	List_hSTS->Add(hSTS_EPAngle_2);
// 	
// 	List_hMC_STS->Add(hMC_STS_Multiplicity);
// 	List_hMC_STS->Add(hMC_STS_Phi);
// 	List_hMC_STS->Add(hMC_STS_Eta);
// 	List_hMC_STS->Add(hMC_STS_pT);
// 	List_hMC_STS->Add(hMC_STS_pT_vs_Eta);
// 	List_hMC_STS->Add(hMC_STS_pT_vs_Phi);
// 	List_hMC_STS->Add(hMC_STS_Eta_vs_Phi);
// 	List_hMC_STS->Add(hMC_STS_AccSTS);
// 
// 	List_hMC_STS->Add(hMC_STS_nocuts_Multiplicity);
// 	List_hMC_STS->Add(hMC_STS_nocuts_Phi);
// 	List_hMC_STS->Add(hMC_STS_nocuts_Eta);
// 	List_hMC_STS->Add(hMC_STS_nocuts_pT);
// 	List_hMC_STS->Add(hMC_STS_nocuts_pT_vs_Eta);
// 	List_hMC_STS->Add(hMC_STS_nocuts_pT_vs_Phi);
// 	List_hMC_STS->Add(hMC_STS_nocuts_Eta_vs_Phi);
// 	List_hMC_STS->Add(hMC_STS_nocuts_AccSTS);
// 	
// 	List_hSTSvsPSD->Add(hSTS_Multiplicity_vs_PSD_Energy);
// 	List_hSTSvsPSD->Add(hSTS_Multiplicity_vs_PSD_1_Energy);
// 	List_hSTSvsPSD->Add(hSTS_Multiplicity_vs_PSD_2_Energy);
// 	List_hSTSvsPSD->Add(hSTS_Multiplicity_vs_PSD_3_Energy);
// 	List_hSTSvsPSD->Add(prfx_hSTS_Multiplicity_vs_PSD_Energy);
//         
// 	List_hSTS_MCvsRECO->Add(hSTS_pT_MC_vs_pT_RECO);
// 	List_hSTS_MCvsRECO->Add(hSTS_phi_MC_vs_phi_RECO);
// 	List_hSTS_MCvsRECO->Add(hSTS_eta_MC_vs_eta_RECO);
// 	List_hSTS_MCvsRECO->Add(hSTS_diff_eta_MC_eta_RECO);
// 	
// 	List_hSTS_MCvsRECO->Add(hSTS_nocuts_pT_MC_vs_pT_RECO);
// 	List_hSTS_MCvsRECO->Add(hSTS_nocuts_phi_MC_vs_phi_RECO);
// 	List_hSTS_MCvsRECO->Add(hSTS_nocuts_eta_MC_vs_eta_RECO);
// 	List_hSTS_MCvsRECO->Add(hSTS_nocuts_diff_eta_MC_eta_RECO);
// 	
// 	
// 	TDirectory *curr_dir = gDirectory->CurrentDirectory();
//         fOutputFile->cd();
// 	
// 	
// 	hPSD_Energy_mod_Energy->Divide(hPSD_Energy_mod_Energy_denominator);
//         hPSD_Energy_mod_Energy_denominator->Write();
//         hPSD_Energy_mod_Energy->Write();
// 	
// //         gDirectory->pwd();
// //         cout << gDirectory->pwd() << endl;
//         
// 	
// 	
// 	
// 	hRPAngle->Write();
// 	List_hPSD_mod_Energy->Write("hPSD_mod_Energy", TObject::kSingleKey);
// 	List_hSTS->Write("hSTS", TObject::kSingleKey);
// 	List_hMC_STS->Write("hMC_STS", TObject::kSingleKey);
// 	List_hPSD->Write("hPSD", TObject::kSingleKey);
// 	List_hSTSvsPSD->Write("hSTSvsPSD", TObject::kSingleKey);
// 	List_hSTS_MCvsRECO->Write("hSTS_MCvsRECO", TObject::kSingleKey);
// 	
// // 	List_hFlow->Write("hFlowlist", TObject::kSingleKey);
// // 	
// // 	List_hFlowQ->Write("hFlowQ", TObject::kSingleKey);
// // 	List_hFlowQQ->Write("hFlowQQ", TObject::kSingleKey);
// // 	List_hFlowpQ->Write("hFlowpQ", TObject::kSingleKey);
// // 	TCanvas* fCanvas= new TCanvas("c1");
// // 	fCanvas->cd();
// // 	hSTS_Multiplicity->Draw();
// // 	fCanvas->Update();
// // 	hMC_STS_Multiplicity->SetLineColor(kRed);
// // 	hMC_STS_Multiplicity->Draw("SAME");
// // 	fCanvas->Update();
// 	
// 	//prfx_hSTS_Multiplicity_vs_PSD_Energy->Write();
// 	//hSTS_Energy->Write();
//         curr_dir->cd();
// // 	TCanvas *c1 = new TCanvas("c1");
// // 	c1->Divide(2,2);
// // 	c1->cd(1);
// // 	hSTS_Multiplicity_vs_PSD_1_Energy->Draw("COL");
// // 	c1->cd(2);
// // 	hSTS_Multiplicity_vs_PSD_2_Energy->Draw("COL");
// // 	c1->cd(3);
// // 	hSTS_Multiplicity_vs_PSD_3_Energy->Draw("COL");
// // 	c1->cd(4);
// // 	hPSD_1_Energy_vs_PSD_2_Energy->Draw("COL");
// // 	c1->SaveAs("./STS_PSD_correlation.png");
// 	cout << "DataTreeMaker::Write_QA_Histograms(): DONE" << endl;
//     }    
// //====================================================================================
// //====================================================================================
// //====================================================================================
// //====================================================================================
//     QdetVector* DataTreeMaker::GetQdetVector(Int_t detVectorSource)
//     {
//       switch (detVectorSource)
//       {
// 	case 0 : // PSD+
// 	  return QdetVector_PSD_plus;
// 	  break;
// 	case 1 : // PSD-
// 	  return QdetVector_PSD_minus;
// 	  break;
// 	case 2 : // STS+
// 	  return QdetVector_STS_plus;
// 	  break;
// 	case 3 : // STS-
// 	  return QdetVector_STS_minus;
// 	  break;
// 	case 4 : // MC STS+
// 	  return QdetVector_MC_STS_plus;
// 	  break;
// 	case 5 : // MC STS-
// 	  return QdetVector_MC_STS_minus;
// 	  break;
// 	default :
// 	  cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
// 	  cout << "ERROR in DataTreeMaker::GetQdetVector : WRONG title for QdetVector! Returned fQdetVector_default ." << endl;
// 	  cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
// 	  return fQdetVector_default;
// 	  break;
//       }
//     }
// //====================================================================================
// //====================================================================================
// //====================================================================================
// //====================================================================================
//     void DataTreeMaker::Clear()
//     {
// 	
//     }
// //====================================================================================
// //====================================================================================
// //====================================================================================
// //====================================================================================
// DataTreeMaker::~DataTreeMaker()
//     {
// 
//     }
// //====================================================================================
// //====================================================================================
// //====================================================================================
// //====================================================================================
// 
// // void DataTreeMaker::SetFlow_Qx(Int_t i, Int_t k, Float_t val)
// // {
// //     fFlowQx[i] = val;
// // //     cout << "ifRecentering = " << ifRecentering << endl;
// // //     if (!ifRecentering)
// // //     {
// // // 	if (k == 0)
// // // 	{
// // // 	    fFlowQx_n1_raw[i] = val;
// // // 	}
// // // 	else
// // // 	{
// // // 	    fFlowQx_n2_raw[i] = val;
// // // 	}
// // //     }
// // //     else
// // //     {
// // // 	if (k == 0)
// // // 	{
// // // 	    fFlowQx_n1_rec[i] = val;
// // // 	}
// // // 	else
// // // 	{
// // // 	    fFlowQx_n2_rec[i] = val;
// // // 	}
// // //     }
// // //     fFlow_Qx[i][k] = val;
// // }
// // //====================================================================================
// // //====================================================================================
// // //====================================================================================
// // //====================================================================================
// // //====================================================================================
// // //====================================================================================
// // 
// // void DataTreeMaker::SetFlow_Qy(Int_t i, Int_t k, Float_t val)
// // {
// //     fFlowQy[i] = val;
// // //     cout << "ifRecentering = " << ifRecentering << endl;
// // //     if (!ifRecentering)
// // //     {
// // // 	if (k == 0)
// // // 	{
// // // 	    fFlowQy_n1_raw[i] = val;
// // // 	}
// // // 	else
// // // 	{
// // // 	    fFlowQy_n2_raw[i] = val;
// // // 	}
// // //     }
// // //     else
// // //     {
// // // 	if (k == 0)
// // // 	{
// // // 	    fFlowQy_n1_rec[i] = val;
// // // 	}
// // // 	else
// // // 	{
// // // 	    fFlowQy_n2_rec[i] = val;
// // // 	}
// // //     }
// // //     fFlow_Qy[i][k] = val;
// // }
// // //====================================================================================
// // //====================================================================================
// // 
// // 
// // 
// // void DataTreeMaker::GetFlowNew()
// // {
// //     CbmStsTrack* track;
// //     CbmTrackMatch* match;
// //     CbmMCTrack* mctrack;
// //     
// //     Double_t p, px, py, pz, energy, mass, pT, phi, eta, MC_px, MC_py, MC_pz, MC_p, MC_pT, MC_phi, MC_eta;
// //     Int_t type, trackID;
// //     
// //     const FairTrackParam *trackParam;
// //     TVector3 momRec;
// // 
// //     Int_t nSTStracks = flistSTSRECOtrack->GetEntriesFast();
// //     //cout << "evenPlane::STSRECOtransverseMomMeth: # STS reco tracks = " << nSTStracks << endl; 
// //     
// //     // Extrapolation track parameters back to primary vertex
// //     vector<CbmStsTrack> vRTracks;
// //     vRTracks.resize(nSTStracks); 
// //     if (nSTStracks > 0 and nSTStracks <10000)
// //     {
// // 	for (Int_t i=0; i<nSTStracks; i++)
// // 	{
// // 	    vRTracks[i] = *( (CbmStsTrack*) flistSTSRECOtrack->At(i));  
// // 	    track = &vRTracks[i];
// // 	    
// // 	    match = (CbmTrackMatch*) flistSTStrackMATCH->At(i);
// // 	    trackID = match->GetMCTrackId();
// // 	    
// // 	    if (trackID < 0) continue; // ghost track from combinatorix (no corresponding MC tracks)
// // 	    
// // 	    mctrack = (CbmMCTrack*) flistSTSMCtrack->At(trackID); 
// // 	    type = mctrack->GetPdgCode();
// // 	    mass = mctrack->GetMass();	
// // 	    MC_px = mctrack->GetPx();
// // 	    MC_py = mctrack->GetPy();
// // 	    MC_pz = mctrack->GetPz();
// // 	    MC_p = TMath::Sqrt(MC_px*MC_px + MC_py*MC_py + MC_pz*MC_pz);
// // 	    MC_pT = TMath::Sqrt(MC_px*MC_px + MC_py*MC_py);
// // 	    MC_phi = TMath::ATan2(MC_py,MC_px);
// // 	    if (MC_phi<-TMath::Pi()) MC_phi+=2*TMath::Pi();
// // 	    if (MC_phi>TMath::Pi()) MC_phi-=2*TMath::Pi();
// // 	    MC_eta = TMath::Log((MC_p+MC_pz)/(MC_p-MC_pz))/2.;
// // 	    
// // 	    trackParam = track->GetParamFirst();	
// // 	    trackParam->Momentum(momRec);
// // 	    
// // 	    px = momRec.X();
// // 	    py = momRec.Y();
// // 	    pz = momRec.Z();
// // 	    pT = TMath::Sqrt(px*px+py*py);
// // 	    p = TMath::Sqrt(px*px+py*py+pz*pz);
// // 	    phi = TMath::ATan2(py,px);//+(TMath::Pi()/2);
// // 	    if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
// // 	    if (phi>TMath::Pi()) phi-=2*TMath::Pi();
// // 	    eta = TMath::Log((p+pz)/(p-pz))/2.;
// // 	    
// // 	    for (Int_t j=0; j<3; j++)
// // 	    {
// // 		hFlowQx[j]->Fill(pT, fFlowQx[j]);
// // 	    }
// // 	    //cout << "pT = " << pT << endl;
// // 	    //cout << "fFlowQx[0] = " << fFlowQx[0] << "; fFlowQx[1] = " << fFlowQx[1] << "; fFlowQx[2] = " << fFlowQx[2] << endl; 
// // 	    hFlowCosQx->Fill(pT, TMath::Cos(2*phi)*fFlowQx[0]);
// // 	    hFlowQxQx[0]->Fill(pT, fFlowQx[0]*fFlowQx[1]);
// // 	    hFlowQxQx[1]->Fill(pT, fFlowQx[0]*fFlowQx[2]);
// // 	    hFlowQxQx[2]->Fill(pT, fFlowQx[1]*fFlowQx[2]);    
// // 	    
// // 	    //cout << "fFlowQx[0]*fFlowQx[1] = " << fFlowQx[0]*fFlowQx[1] << "; fFlowQx[0]*fFlowQx[2]" << fFlowQx[0]*fFlowQx[2] <<  "; fFlowQx[1]*fFlowQx[2]" << fFlowQx[1]*fFlowQx[2] << endl;
// // 	}
// //     }
// // }
// // 
// // 
// // //====================================================================================
// // //====================================================================================
// // //====================================================================================
// // //====================================================================================
// // void DataTreeMaker::CalculateQ()
// // {
// //     for (Int_t v=0; v<5; v++)
// //     {
// // 	for (Int_t s=0; s<2; s++)
// // 	{
// // 	    for (Int_t a=0; a<5; a++)
// // 	    {
// // 		hFlow_Qx_n1[v][a][s] -> Fill(fFlowVar[v],fFlow_Qx_n1[v][a][s]);
// // 		hFlow_Qy_n1[v][a][s] -> Fill(fFlowVar[v],fFlow_Qy_n1[v][a][s]);
// // 		hFlow_Qx_n2[v][a][s] -> Fill(fFlowVar[v],fFlow_Qx_n2[v][a][s]);
// // 		hFlow_Qy_n2[v][a][s] -> Fill(fFlowVar[v],fFlow_Qy_n2[v][a][s]);
// // 		for (Int_t b=0; b<5; b++)
// // 		{
// // 		    hFlow_QxQx_n1[v][a][b][s] -> Fill(fFlowVar[v],fFlow_Qx_n1[v][a][s]*fFlow_Qx_n1[v][b][s]);
// // 		    hFlow_QyQy_n1[v][a][b][s] -> Fill(fFlowVar[v],fFlow_Qy_n1[v][a][s]*fFlow_Qy_n1[v][b][s]);
// // 		    hFlow_QyQx_n1[v][a][b][s] -> Fill(fFlowVar[v],fFlow_Qy_n1[v][a][s]*fFlow_Qx_n1[v][b][s]);
// // 		    hFlow_QxQx_n2[v][a][b][s] -> Fill(fFlowVar[v],fFlow_Qx_n2[v][a][s]*fFlow_Qx_n2[v][b][s]);
// // 		    hFlow_QyQy_n2[v][a][b][s] -> Fill(fFlowVar[v],fFlow_Qy_n2[v][a][s]*fFlow_Qy_n2[v][b][s]);
// // 		    hFlow_QyQx_n2[v][a][b][s] -> Fill(fFlowVar[v],fFlow_Qy_n2[v][a][s]*fFlow_Qx_n2[v][b][s]);
// // 		}
// // 	    }
// // 	}
// //     }
// // }
// // //====================================================================================
// // //====================================================================================
// // //====================================================================================
// // //====================================================================================
// // 
// // void DataTreeMaker::SetFlow_Q_n(Int_t x, Int_t n, Int_t v, Int_t a, Int_t s , Float_t val)
// // {
// //     if (x==0)
// //     {
// // 	if (n==1)
// // 	{
// // 	    //cout << "fFlow_Qx_n1[" << v << "][" << a << "][" << s << "] = " << val << endl;
// // 	    fFlow_Qx_n1[v][a][s] = val;
// // 	}
// // 	if (n==2)
// // 	{
// // 	    //cout << "fFlow_Qx_n2[" << v << "][" << a << "][" << s << "] = " << val << endl;
// // 	    fFlow_Qx_n2[v][a][s] = val;
// // 	}
// //     }
// //     if (x==1)
// //     {
// // 	if (n==1)
// // 	{
// // 	    //cout << "fFlow_Qy_n1[" << v << "][" << a << "][" << s << "] = " << val << endl;
// // 	    fFlow_Qy_n1[v][a][s] = val;
// // 	}
// // 	if (n==2)
// // 	{
// // 	    //cout << "fFlow_Qy_n2[" << v << "][" << a << "][" << s << "] = " << val << endl;
// // 	    fFlow_Qy_n2[v][a][s] = val;
// // 	}
// //     }
// // }
// 
// //====================================================================================
// //====================================================================================
