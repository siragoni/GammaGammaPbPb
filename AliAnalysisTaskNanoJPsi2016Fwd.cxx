/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// c++ headers
#include <iostream>
#include <fstream>

// root headers
#include <TMath.h>
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TFile.h>
#include <TF2.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraph2D.h>
#include <TStopwatch.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>
#include <TLatex.h>
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TObjString.h"
#include "TList.h"
#include "TChain.h"


// aliroot headers
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMuonTrackCuts.h"
#include "AliMCEvent.h"
#include <AliMCParticle.h>
#include <AliAODMCParticle.h>

// my headers
#include "AliAnalysisTaskNanoJPsi2016Fwd.h"



class AliAnalysisTaskNanoJPsi2016Fwd;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskNanoJPsi2016Fwd) // classimp: necessary for root


AliAnalysisTaskNanoJPsi2016Fwd::AliAnalysisTaskNanoJPsi2016Fwd() : AliAnalysisTaskSE(),
  gMuonMass(0.105658), fPeriod(0), fIsMC(0),
  fAOD(0), fOutputList(0),fCounterH(0),
  fMuonTrackCounterH(0),
  fTriggerCounterFwdH(0),
  fGoodTrkFwdH(0),
  fMuonTrackCuts(0x0), fAnaTree(0), fAnaTreeMC(0),
  fRunNum(0), fL0inputs(0), fTracklets(0), fAnaType(-1),
  fGoodPosTrk(0), fGoodNegTrk(0),
  fZNAfired(0), fZNCfired(), fZNCEnergy(0), fZNAEnergy(0), fZPCEnergy(0), fZPAEnergy(0),
  fV0ADecision(-10), fV0CDecision(-10), fV0ATime(0), fV0CTime(0), fADADecision(-10), fADCDecision(-10), fADATime(0), fADCTime(0),
  fIR1Map(0), fIR2Map(0), fTrkTrkPt(0), fTrkTrkPhi(0), fTrkTrkY(0), fTrkTrkM(0), fTrkPt1(0), fTrkPt2(0),
  fTrkEta1(0), fTrkEta2(0), fTrkPhi1(0), fTrkPhi2(0), fTrkQ1(0), fTrkQ2(0), fTrkRabs1(0), fTrkRabs2(0),
  fMCTrkTrkPt(0), fMCTrkTrkPhi(0), fMCTrkTrkY(0), fMCTrkTrkM(0), fMCTrkPt1(0), fMCTrkPt2(0),
  fMCTrkEta1(0), fMCTrkEta2(0), fMCTrkPhi1(0), fMCTrkPhi2(0), fMCTrkQ1(0), fMCTrkQ2(0),
  fV0Hits{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  fV0TotalNCells(0),
  fMCCosThetaHE(0),fMCCosThetaCS(0),fMCPhiHE(0),fMCPhiCS(0),fMCTildePhiHEpos(0),fMCTildePhiHEneg(0),fMCTildePhiCSpos(0),
  fMCTildePhiCSneg(0),fCosThetaHE(0),fCosThetaCS(0),fPhiHE(0),fPhiCS(0),fTildePhiHEpos(0),fTildePhiHEneg(0),
  fTildePhiCSpos(0),fTildePhiCSneg(0)

{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}

//_____________________________________________________________________________
AliAnalysisTaskNanoJPsi2016Fwd::AliAnalysisTaskNanoJPsi2016Fwd(const char* name) : AliAnalysisTaskSE(name),
 gMuonMass(0.105658),fPeriod(0), fIsMC(0),
 fAOD(0), fOutputList(0),fCounterH(0),
  fMuonTrackCounterH(0),
  fTriggerCounterFwdH(0),
  fGoodTrkFwdH(0),
  fMuonTrackCuts(0x0), fAnaTree(0), fAnaTreeMC(0),
  fRunNum(0), fL0inputs(0), fTracklets(0), fAnaType(-1),
  fGoodPosTrk(0), fGoodNegTrk(0),
  fZNAfired(0), fZNCfired(), fZNCEnergy(0), fZNAEnergy(0), fZPCEnergy(0), fZPAEnergy(0),
  fV0ADecision(-10), fV0CDecision(-10), fV0ATime(0), fV0CTime(0),fADADecision(-10), fADCDecision(-10), fADATime(0), fADCTime(0),
  fIR1Map(0), fIR2Map(0), fTrkTrkPt(0), fTrkTrkPhi(0), fTrkTrkY(0), fTrkTrkM(0), fTrkPt1(0), fTrkPt2(0),
  fTrkEta1(0), fTrkEta2(0), fTrkPhi1(0), fTrkPhi2(0), fTrkQ1(0), fTrkQ2(0), fTrkRabs1(0), fTrkRabs2(0),
  fMCTrkTrkPt(0), fMCTrkTrkPhi(0), fMCTrkTrkY(0), fMCTrkTrkM(0), fMCTrkPt1(0), fMCTrkPt2(0),
  fMCTrkEta1(0), fMCTrkEta2(0), fMCTrkPhi1(0), fMCTrkPhi2(0), fMCTrkQ1(0), fMCTrkQ2(0),
  fV0Hits{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  fV0TotalNCells(0),
  fMCCosThetaHE(0),fMCCosThetaCS(0),fMCPhiHE(0),fMCPhiCS(0),fMCTildePhiHEpos(0),fMCTildePhiHEneg(0),fMCTildePhiCSpos(0),
  fMCTildePhiCSneg(0),fCosThetaHE(0),fCosThetaCS(0),fPhiHE(0),fPhiCS(0),fTildePhiHEpos(0),fTildePhiHEneg(0),
  fTildePhiCSpos(0),fTildePhiCSneg(0)
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TTree::Class());
    DefineOutput(2, TList::Class());
    DefineOutput(3, TTree::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskNanoJPsi2016Fwd::~AliAnalysisTaskNanoJPsi2016Fwd()
{
    // destructor
    // liberate all allocated memory
    if(fOutputList) {delete fOutputList;}
    if(fMuonTrackCuts) {delete fMuonTrackCuts;}
    if(fAnaTree) {delete fAnaTree;}
    if(fCounterH) {delete fCounterH;}
    if(fMuonTrackCounterH) {delete fMuonTrackCounterH;}
    if(fTriggerCounterFwdH) {delete fTriggerCounterFwdH;}
    if(fGoodTrkFwdH) {delete fGoodTrkFwdH;}

}
//_____________________________________________________________________________
void AliAnalysisTaskNanoJPsi2016Fwd::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)

  ////////////////////////////////////////
  //muon track cuts
  ////////////////////////////////////////
  fMuonTrackCuts = new AliMuonTrackCuts("StdMuonCuts", "StdMuonCuts");

  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuEta|AliMuonTrackCuts::kMuThetaAbs|AliMuonTrackCuts::kMuPdca|AliMuonTrackCuts::kMuMatchLpt);
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  fMuonTrackCuts->Print("mask");
  fMuonTrackCuts->SetIsMC();

  ////////////////////////////////////////
  //output tree
  ////////////////////////////////////////
  fAnaTree = new TTree("fOutputTree", "fOutputTree");
  fAnaTree ->Branch("fRunNum", &fRunNum, "fRunNum/I");
  // fAnaTree ->Branch("fL0inputs", &fL0inputs, "fL0inputs/i");
  // fAnaTree ->Branch("fAnaType", &fAnaType, "fAnaType/I");
  // fAnaTree ->Branch("fGoodPosTrk", &fGoodPosTrk, "fGoodPosTrk/I");
  // fAnaTree ->Branch("fGoodNegTrk", &fGoodNegTrk, "fGoodNegTrk/I");
  fAnaTree ->Branch("fTrkTrkPt", &fTrkTrkPt, "fTrkTrkPt/D");
  // fAnaTree ->Branch("fTrkTrkPhi", &fTrkTrkPhi, "fTrkTrkPhi/D");
  fAnaTree ->Branch("fTrkTrkY", &fTrkTrkY, "fTrkTrkY/D");
  fAnaTree ->Branch("fTrkTrkM", &fTrkTrkM, "fTrkTrkM/D");




  fAnaTree ->Branch("fV0ADecision",   &fV0ADecision,   "fV0ADecision/I");
  fAnaTree ->Branch("fV0CDecision",   &fV0CDecision,   "fV0CDecision/I");
  fAnaTree ->Branch("fADADecision",   &fADADecision,   "fADADecision/I");
  fAnaTree ->Branch("fADCDecision",   &fADCDecision,   "fADCDecision/I");
  fAnaTree ->Branch("fV0TotalNCells", &fV0TotalNCells, "fV0TotalNCells/I");



  fAnaTree ->Branch("fTrkPt1", &fTrkPt1, "fTrkPt1/D");
  fAnaTree ->Branch("fTrkPt2", &fTrkPt2, "fTrkPt2/D");
  fAnaTree ->Branch("fTrkEta1", &fTrkEta1, "fTrkEta1/D");
  fAnaTree ->Branch("fTrkEta2", &fTrkEta2, "fTrkEta2/D");
  fAnaTree ->Branch("fTrkPhi1", &fTrkPhi1, "fTrkPhi1/D");
  fAnaTree ->Branch("fTrkPhi2", &fTrkPhi2, "fTrkPhi2/D");
  fAnaTree ->Branch("fTrkQ1", &fTrkQ1, "fTrkQ1/D");
  fAnaTree ->Branch("fTrkQ2", &fTrkQ2, "fTrkQ2/D");
  // fAnaTree ->Branch("fTrkRabs1", &fTrkRabs1, "fTrkRabs1/D");
  // fAnaTree ->Branch("fTrkRabs2", &fTrkRabs2, "fTrkRabs2/D");
  fAnaTree ->Branch("fCosThetaHE", &fCosThetaHE, "fCosThetaHE/D");
  fAnaTree ->Branch("fCosThetaCS", &fCosThetaCS, "fCosThetaCS/D");
  fAnaTree ->Branch("fPhiHE", &fPhiHE, "fPhiHE/D");
  fAnaTree ->Branch("fPhiCS", &fPhiCS, "fPhiCS/D");

   // post data
  PostData(1, fAnaTree);

  ////////////////////////////////////////
  //output histograms
  ////////////////////////////////////////

  fOutputList = new TList();          // this is a list which will contain all  histograms
  fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
  //  counter for events passing each cut
  fCounterH = new TH1F("fCounterH", "fCounterH", 25, -0.5, 24.5);
  fCounterH->GetXaxis()->SetBinLabel(0, "Total events");
  fCounterH->GetXaxis()->SetBinLabel(1, "AOD events");
  fCounterH->GetXaxis()->SetBinLabel(2, "Right trigger");
  fCounterH->GetXaxis()->SetBinLabel(3, "Right run");
  fCounterH->GetXaxis()->SetBinLabel(4, "At least one trk");
  fCounterH->GetXaxis()->SetBinLabel(5, "Two good muons");
  fCounterH->GetXaxis()->SetBinLabel(6, "OS muons");
  fCounterH->GetXaxis()->SetBinLabel(7, "All cuts");
  fOutputList->Add(fCounterH);
  //  counter for tracks passing each cut
  fMuonTrackCounterH = new TH1F("fMuonTrackCounterH", "fMuonTrackCounterH", 10, -0.5, 9.5);
  fOutputList->Add(fMuonTrackCounterH);

  // post data
  PostData(2, fOutputList);

  ////////////////////////////////////////
  //MC tree at generator level
  ////////////////////////////////////////

  fAnaTreeMC = new TTree("fOutputTreeMC", "fOutputTreeMC");
  // // // fAnaTreeMC ->Branch("fRunNum", &fRunNum, "fRunNum/I");
  // fAnaTreeMC ->Branch("fMCTrkTrkPt", &fMCTrkTrkPt, "fMCTrkTrkPt/D");
  // // // fAnaTreeMC ->Branch("fMCTrkTrkPhi", &fMCTrkTrkPhi, "fMCTrkTrkPhi/D");
  // fAnaTreeMC ->Branch("fMCTrkTrkY", &fMCTrkTrkY, "fMCTrkTrkY/D");
  // fAnaTreeMC ->Branch("fMCTrkTrkM", &fMCTrkTrkM, "fMCTrkTrkM/D");
  //
  // fAnaTreeMC ->Branch("fMCCosThetaHE", &fMCCosThetaHE, "fMCCosThetaHE/D");
  // fAnaTreeMC ->Branch("fMCCosThetaCS", &fMCCosThetaCS, "fMCCosThetaCS/D");
  // fAnaTreeMC ->Branch("fMCPhiHE", &fMCPhiHE, "fMCPhiHE/D");
  // fAnaTreeMC ->Branch("fMCPhiCS", &fMCPhiCS, "fMCPhiCS/D");
  // fAnaTreeMC ->Branch("fMCTildePhiHEpos", &fMCTildePhiHEpos, "fMCTildePhiHEpos/D");
  // fAnaTreeMC ->Branch("fMCTildePhiHEneg", &fMCTildePhiHEneg, "fMCTildePhiHEneg/D");
  // fAnaTreeMC ->Branch("fMCTildePhiCSpos", &fMCTildePhiCSpos, "fMCTildePhiCSpos/D");
  // fAnaTreeMC ->Branch("fMCTildePhiCSneg", &fMCTildePhiCSneg, "fMCTildePhiCSneg/D");
  // // // fAnaTreeMC ->Branch("fMCTrkPt1", &fMCTrkPt1, "fMCTrkPt1/D");
  // // // fAnaTreeMC ->Branch("fMCTrkPt2", &fMCTrkPt2, "fMCTrkPt2/D");
  // // // fAnaTreeMC ->Branch("fMCTrkEta1", &fMCTrkEta1, "fMCTrkEta1/D");
  // // // fAnaTreeMC ->Branch("fMCTrkEta2", &fMCTrkEta2, "fMCTrkEta2/D");
  // // // fAnaTreeMC ->Branch("fMCTrkPhi1", &fMCTrkPhi1, "fMCTrkPhi1/D");
  // // // fAnaTreeMC ->Branch("fMCTrkPhi2", &fMCTrkPhi2, "fMCTrkPhi2/D");
  // // // fAnaTreeMC ->Branch("fMCTrkQ1", &fMCTrkQ1, "fMCTrkQ1/I");
  // // // fAnaTreeMC ->Branch("fMCTrkQ2", &fMCTrkQ2, "fMCTrkQ2/I");

  // post data
  PostData(3, fAnaTreeMC);

}
//_____________________________________________________________________________
void AliAnalysisTaskNanoJPsi2016Fwd::NotifyRun()
{
  /// Set run number for cuts
 fMuonTrackCuts->SetRun(fInputHandler);
}

//_____________________________________________________________________________
void AliAnalysisTaskNanoJPsi2016Fwd::SetPeriod(Int_t period)
{
  // period = 0 => 2016 s, = 1 => 2016 r
  fPeriod = period;
}

//_____________________________________________________________________________
void AliAnalysisTaskNanoJPsi2016Fwd::SetMC(Bool_t flag)
{
  // set if MC file
  fIsMC = flag;
}
//_____________________________________________________________________________
void AliAnalysisTaskNanoJPsi2016Fwd::TrkTrkKine(AliAODTrack *Track1, AliAODTrack *Track2, Double_t TrkMass)
{
  // --  positive track
  TLorentzVector PosLV;
  PosLV.SetPtEtaPhiM(Track1->Pt(), Track1->Eta(), Track1->Phi(), TrkMass);
  // --  negative track
  TLorentzVector NegLV;
  NegLV.SetPtEtaPhiM(Track2->Pt(), Track2->Eta(), Track2->Phi(), TrkMass);

  // vector of Trk+Trk
  TLorentzVector TrkTrk = NegLV+PosLV;

  // set tree variables
  fTrkTrkPt = TrkTrk.Pt();
  fTrkTrkPhi = TrkTrk.Phi();
  fTrkTrkY = TrkTrk.Rapidity();
  fTrkTrkM = TrkTrk.M();
  fTrkPt1 = Track1->Pt();
  fTrkPt2 = Track2->Pt();
  fTrkEta1 = Track1->Eta();
  fTrkEta2 = Track2->Eta();
  fTrkPhi1 = Track1->Phi();
  fTrkPhi2 = Track2->Phi();
  fTrkQ1 = Track1->Charge();
  fTrkQ2 = Track2->Charge();
  fTrkRabs1 = Track1->GetRAtAbsorberEnd();
  fTrkRabs2 = Track2->GetRAtAbsorberEnd();

  fCosThetaHE   = CosThetaHelicityFrame(PosLV,NegLV,TrkTrk);
  fCosThetaCS   = CosThetaCollinsSoper( PosLV,NegLV,TrkTrk);
  fPhiHE   = CosPhiCollinsSoper(   PosLV,NegLV,TrkTrk);
  fPhiCS  = CosPhiHelicityFrame( PosLV,NegLV,TrkTrk );
  fTildePhiHEpos  = fPhiHE - 0.25 * TMath::Pi() ;
  fTildePhiHEneg  = fPhiHE - 0.75 * TMath::Pi() ;
  fTildePhiCSpos  = fPhiCS - 0.25 * TMath::Pi() ;
  fTildePhiCSneg  = fPhiCS - 0.75 * TMath::Pi() ;
  if( fTildePhiHEpos < 0. ) {
    fTildePhiHEpos += 2. * TMath::Pi();
  }
  if( fTildePhiHEneg < 0. ) {
    fTildePhiHEneg += 2. * TMath::Pi();
  }
  if( fTildePhiCSpos < 0. ) {
    fTildePhiCSpos += 2. * TMath::Pi();
  }
  if( fTildePhiCSneg < 0. ) {
    fTildePhiCSneg += 2. * TMath::Pi();
  }



}

//_____________________________________________________________________________

void AliAnalysisTaskNanoJPsi2016Fwd::CheckTrigger(Bool_t *isTriggered)
// checks if event is triggered according to period and analysis type
{
  // Initialise: 0 = fwd, 1 = cent, 2 = semi-fwd
  isTriggered[0] = isTriggered[1] = isTriggered[2] = kFALSE;

  // read trigger info
  TString trigger = fAOD->GetFiredTriggerClasses();

  // forward analysis
  // in 2016 r : CMUP14-B-NOPF-MUFAST = 0MSL *0VBA *0UBA
  // in 2016 s : CMUP23-B-NOPF-MUFAST = 0MUL *0UBC *0UGC *0VBA *0VGA *0SH2 *0VC5
  // central analysis
  // in 2016 r CCUP20-B-NOPF-CENTNOTRD = !VBA & !VGA & !VBC & !UBA & !UGC & !SH2 & STG & OMU
  // in 2016 s CCUP22-B-SPD2-CENTNOTRD = !VBA & !VGA & !VBC & !UBC & !UGC & !SH2 & STG & OMU
  // semi-forward analysis
  // in 2016 r CMUP15-B-NOPF-ALLNOTRD = *0VBA *0UBA *0VC5 0SMB *0SH2 0MSL
  // in 2016 s CMUP22-B-NOPF-ALLNOTRD = *0UBC *0UGC *0VBA *0VGA *0SH2 *0VC5 0MSL 0SMB
  if (fPeriod == 0) {
    if (trigger.Contains("CMUP11-B-NOPF-MUFAST") ||
	          trigger.Contains("CMUP26-B-NOPF-MUFAST") ||
	          trigger.Contains("CMUP6-B-NOPF-MUFAST")  ||
            trigger.Contains("CMUP10-B-NOPF-MUFAST") ||
            trigger.Contains("CMUP13-B-NOPF-MUFAST") ) isTriggered[0] = kTRUE;
    if (trigger.Contains("CMUP11-B-NOPF-MUFAST") ||
	          trigger.Contains("CMUP26-B-NOPF-MUFAST") ||
	          trigger.Contains("CMUP6-B-NOPF-MUFAST")  ||
            trigger.Contains("CMUP10-B-NOPF-MUFAST") ||
            trigger.Contains("CMUP13-B-NOPF-MUFAST") ) isTriggered[1] = kTRUE;
    if (trigger.Contains("CMUP11-B-NOPF-MUFAST") ||
	          trigger.Contains("CMUP26-B-NOPF-MUFAST") ||
	          trigger.Contains("CMUP6-B-NOPF-MUFAST")  ||
            trigger.Contains("CMUP10-B-NOPF-MUFAST") ||
            trigger.Contains("CMUP13-B-NOPF-MUFAST") ) isTriggered[2] = kTRUE;
  } else if (fPeriod == 1) {
    if (trigger.Contains("CMUP11-B-NOPF-MUFAST") ||
	          trigger.Contains("CMUP26-B-NOPF-MUFAST") ||
	          trigger.Contains("CMUP6-B-NOPF-MUFAST")  ||
            trigger.Contains("CMUP10-B-NOPF-MUFAST") ||
            trigger.Contains("CMUP13-B-NOPF-MUFAST") ) isTriggered[0] = kTRUE;
    if (trigger.Contains("CMUP11-B-NOPF-MUFAST") ||
	          trigger.Contains("CMUP26-B-NOPF-MUFAST") ||
	          trigger.Contains("CMUP6-B-NOPF-MUFAST")  ||
            trigger.Contains("CMUP10-B-NOPF-MUFAST") ||
            trigger.Contains("CMUP13-B-NOPF-MUFAST") ) isTriggered[1] = kTRUE;
    if (trigger.Contains("CMUP11-B-NOPF-MUFAST") ||
	          trigger.Contains("CMUP26-B-NOPF-MUFAST") ||
	          trigger.Contains("CMUP6-B-NOPF-MUFAST")  ||
            trigger.Contains("CMUP10-B-NOPF-MUFAST") ||
            trigger.Contains("CMUP13-B-NOPF-MUFAST") ) isTriggered[2] = kTRUE;
  }
}



//_____________________________________________________________________________

Bool_t AliAnalysisTaskNanoJPsi2016Fwd::GoodMUONTrack(Int_t iTrack)
// selects good MUON tracks
{
  fMuonTrackCounterH->Fill(0);
  AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
  if(!track) return kFALSE;
  fMuonTrackCounterH->Fill(1);
  if(!track->IsMuonTrack())  return kFALSE;
  fMuonTrackCounterH->Fill(2);
  if(!fMuonTrackCuts->IsSelected(track))  return kFALSE;
  fMuonTrackCounterH->Fill(3);
  if(track->GetRAtAbsorberEnd()>89.5)  return kFALSE;
  fMuonTrackCounterH->Fill(4);
  if(track->GetRAtAbsorberEnd()<17.5)  return kFALSE;
  fMuonTrackCounterH->Fill(5);
  return kTRUE;
}


//_____________________________________________________________________________
void AliAnalysisTaskNanoJPsi2016Fwd::UserExec(Option_t *)
{


  ////////////////////////////////////////////
  // general info of the event
  ////////////////////////////////////////////
  Int_t iSelectionCounter = 0; // no selection applied yet
  // How many events
  fCounterH->Fill(iSelectionCounter); // entering UserExec
  iSelectionCounter++;

  fCosThetaHE   = -999;
  fCosThetaCS   = -999;
  fPhiHE        = -999;
  fPhiCS        = -999;
  fTildePhiHEpos  = -999;
  fTildePhiHEneg  = -999;
  fTildePhiCSpos  = -999;
  fTildePhiCSneg  = -999;
  fTrkTrkPt = -10;
  fTrkTrkPhi = -10;
  fTrkTrkY = -10;
  fTrkTrkM = -10;
  fTrkEta1 = -10;
  fTrkEta2 = -10;
  fTrkPhi1 = -10;
  fTrkPhi2 = -10;
  fTrkQ1 = -10;
  fTrkQ2 = -10;


  // get AOD event
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) {
    // fAnaTree->Fill();
    PostData(1, fAnaTree);
    PostData(2, fOutputList);
    return;
  }
  fCounterH->Fill(iSelectionCounter); // AOD event found (1)
  iSelectionCounter++;


  // cout << "AAAAAAAAA" << endl;

  // get the run number
  fRunNum = fAOD ->GetRunNumber();

  // is the right trigger?
  Bool_t isTriggered[3];
  isTriggered[0]=isTriggered[1]=isTriggered[2]=kFALSE;
  fIsMC = 0;
  if (fIsMC) {
    isTriggered[0] = kTRUE;
    GetMCInfo();
  } else {
    CheckTrigger(isTriggered);
  }
  if (!isTriggered[0]) {
    // fAnaTree->Fill();
    PostData(1, fAnaTree);
    PostData(2, fOutputList);
    return;
  }
  fCounterH->Fill(iSelectionCounter); // right trigger found (2)
  iSelectionCounter++;
  // trigger inputs
  fL0inputs = fAOD->GetHeader()->GetL0TriggerInputs();
  // cout << "BBBBBBB" << endl;
  Int_t listOfGoodRunNumbersLHC18q[] = { 295585, 295586, 295587, 295588, 295589, 295612,
                                         295615, 295665, 295666, 295667, 295668, 295671,
                                         295673, 295675, 295676, 295677, 295714, 295716,
                                         295717, 295718, 295719, 295723, 295725, 295753,
                                         295754, 295755, 295758, 295759, 295762, 295763,
                                         295786, 295788, 295791, 295816, 295818, 295819,
                                         295822, 295825, 295826, 295829, 295831, 295854,
                                         295855, 295856, 295859, 295860, 295861, 295863,
                                         295881, 295908, 295909, 295910, 295913, 295936,
                                         295937, 295941, 295942, 295943, 295945, 295947,
                                         296061, 296062, 296063, 296065, 296066, 296068,
                                         296123, 296128, 296132, 296133, 296134, 296135,
                                         296142, 296143, 296191, 296192, 296194, 296195,
                                         296196, 296197, 296198, 296241, 296242, 296243,
                                         296244, 296246, 296247, 296269, 296270, 296273,
                                         296279, 296280, 296303, 296304, 296307, 296309,
                                         296312, /*296376,*/ 296377, 296378, 296379, 296380,
                                         296381, 296383, 296414, 296419, 296420, 296423,
                                         296424, 296433, 296472, 296509, 296510, 296511,
                                         296514, 296516, 296547, 296548, 296549, 296550,
                                         296551, 296552, 296553, 296615, 296616, 296618,
                                         296619, 296622, 296623 };
  Int_t listOfGoodRunNumbersLHC18r[] = { 296690, 296691, 296694, 296749, 296750, 296781,
                                         296784, 296785, 296786, 296787, 296791, 296793,
                                         296794, 296799, 296836, 296838, 296839, 296848,
                                         296849, 296850, 296851, 296852, 296890, 296894,
                                         296899, 296900, 296903, 296930, 296931, 296932,
                                         296934, 296935, 296938, 296941, 296966, 296967,
                                         296968, 296969, 296971, 296975, 296976, /*296977,*/
                                         296979, 297029, 297031, 297035, 297085, 297117,
                                         297118, 297119, 297123, 297124, 297128, 297129,
                                         297132, 297133, 297193, 297194, 297196, 297218,
                                         297219, 297221, 297222, 297278, 297310, 297312,
                                         297315, 297317, 297363, 297366, 297367, 297372,
                                         297379, 297380, 297405, 297408, 297413, 297414,
                                         297415, 297441, 297442, 297446, 297450, 297451,
                                         297452, 297479, 297481, 297483, 297512, 297537,
                                         297540, 297541, 297542, 297544, 297558, 297588,
                                         297590, 297595/*, 297623, 297624*/ };
  /* - This good run number list has been taken from the analysis
     - note of Kay's talk for DIS 2017, see:
     - https://alice-notes.web.cern.ch/system/files/notes/analysis/596/2017-Feb-08-analysis_note-2017-Feb-08-analysis-note.pdf
     -
   */
  Int_t listOfGoodRunNumbersLHC15o[] = { /*244918,*/ 244980, 244982, 244983, 245064, 245066, 245068, 245145, 245146, 245151,
                                         245152, 245231, 245232, 245233, 245253, 245259, 245343, 245345, 245346, 245347,
                                         245353, 245401, 245407, 245409, 245410, 245446, 245450, 245496, 245501, 245504,
                                         245505, 245507, 245535, 245540, 245542, 245543, 245554, 245683, 245692, 245700,
                                         245705, 245729, 245731, 245738, 245752, 245759, 245766, 245775, 245785, 245793,
                                         245829, 245831, 245833, 245949, 245952, 245954, 245963, 245996, 246001, 246003,
                                         246012, 246036, 246037, 246042, 246048, 246049, 246053, 246087, 246089, 246113,
                                         246115, 246148, 246151, 246152, 246153, 246178, 246181, 246182, 246217, 246220,
                                         246222, 246225, 246272, 246275, 246276, 246390, 246391, 246392, 246424, 246428,
                                         246431, 246433, 246434, 246487, 246488, 246493, 246495, 246675, 246676, 246750,
                                         246751, 246755, 246757, 246758, 246759, 246760, 246763, 246765, 246804, 246805,
                                         246806, 246807, 246808, 246809, 246844, 246845, 246846, 246847, 246851, 246855,
                                         246859, 246864, 246865, 246867, 246871, 246930, 246937, 246942, 246945, 246948,
                                         246949, 246980, 246982, 246984, 246989, 246991, 246994
                                       };
  Bool_t checkIfGoodRun = kFALSE;
  for( Int_t iRunLHC18q = 0; iRunLHC18q < 128; iRunLHC18q++){
  // for( Int_t iRunLHC18q = 0; iRunLHC18q < 129; iRunLHC18q++){
  // for( Int_t iRunLHC18q = 0; iRunLHC18q < 125; iRunLHC18q++){
    if( fRunNum == listOfGoodRunNumbersLHC18q[iRunLHC18q] ) checkIfGoodRun = kTRUE;
  }
  for( Int_t iRunLHC18r = 0; iRunLHC18r <  97; iRunLHC18r++){
  // for( Int_t iRunLHC18r = 0; iRunLHC18r <  98; iRunLHC18r++){
  // for( Int_t iRunLHC18r = 0; iRunLHC18r <  82; iRunLHC18r++){
    if( fRunNum == listOfGoodRunNumbersLHC18r[iRunLHC18r] ) checkIfGoodRun = kTRUE;
  }
  for( Int_t iRunLHC15o = 0; iRunLHC15o < 136/*137*/; iRunLHC15o++){
  // for( Int_t iRunLHC15o = 0; iRunLHC15o < 134; iRunLHC15o++){
    if( fRunNum == listOfGoodRunNumbersLHC15o[iRunLHC15o] ) checkIfGoodRun = kTRUE;
  }
  if(checkIfGoodRun != 1) {
    // fAnaTree->Fill();
    PostData(1, fAnaTree);
    PostData(2, fOutputList);
    return;
  }
  fCounterH->Fill(iSelectionCounter); // right run (3)
  iSelectionCounter++;


  // if (isTriggered[0]) fTriggerCounterFwdH->Fill(fRunNum);

  ////////////////////////////////////////////
  //  find events with two good tracks
  ////////////////////////////////////////////

  //are there tracks at all?
  Int_t nTracks(fAOD->GetNumberOfTracks());
  if(nTracks<1) {
    // fAnaTree->Fill();
    PostData(1, fAnaTree);
    PostData(2, fOutputList);
    return;
  }
  fCounterH->Fill(iSelectionCounter); // At least one track (1)
  iSelectionCounter++;

  // loop over tracks and select good muons
  Int_t nGoodPosTrk = 0;
  Int_t nGoodNegTrk = 0;
  Int_t *idxPosTrk = new Int_t[nTracks];
  Int_t *idxNegTrk = new Int_t[nTracks];
  Int_t nGoodMUON = 0;

  for(Int_t iTrack(0); iTrack < nTracks; iTrack++) {
    // initialise
    Bool_t isGoodMUON = kFALSE;

    // check if it is muon, if not check if it is central
    isGoodMUON = GoodMUONTrack(iTrack);

    // if valid track
    if (isGoodMUON) {
      nGoodMUON++; // update counter
      AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
      // set charge, increase counter and store indices
      if (track->Charge() > 0) {
	idxPosTrk[nGoodPosTrk] = iTrack;
	nGoodPosTrk++;
      } else  if (track->Charge() < 0) {
	idxNegTrk[nGoodNegTrk] = iTrack;
	nGoodNegTrk++;
      }
    }
  }

  // fill in number of good tracks
  fGoodPosTrk = nGoodPosTrk;
  fGoodNegTrk = nGoodNegTrk;
  // if (isTriggered[0]) fGoodTrkFwdH->Fill(nGoodPosTrk,nGoodNegTrk);

  // cout << "CCCCCCCC" << endl;


  ////////////////////////////////////////////
  // two track analysis
  ////////////////////////////////////////////

  // set type of analysis
  fAnaType = -1;
  if (isTriggered[0] && (nGoodMUON == 2)) fAnaType = 0;

  // check valid analysis type
  if ( fAnaType == -1 ) {
    // fAnaTree->Fill();
    PostData(1, fAnaTree);
    PostData(2, fOutputList);
    delete [] idxNegTrk;
    delete [] idxPosTrk;
    return;
  }
  fCounterH->Fill(iSelectionCounter); // valid analysis (5)
  iSelectionCounter++;

  // compute trk+trk kinematics
  AliAODTrack *Track1 = NULL;
  AliAODTrack *Track2 = NULL;
  if (nGoodPosTrk == 1 && nGoodNegTrk == 1) { // exactly one positive and one negative track
    fCounterH->Fill(iSelectionCounter);
    Track1 = static_cast<AliAODTrack*>(fAOD->GetTrack(idxPosTrk[0]));
    Track2 = static_cast<AliAODTrack*>(fAOD->GetTrack(idxNegTrk[0]));
  } else if (nGoodPosTrk == 2 && nGoodNegTrk == 0) { //two positive  tracks
    Track1 = static_cast<AliAODTrack*>(fAOD->GetTrack(idxPosTrk[0]));
    Track2 = static_cast<AliAODTrack*>(fAOD->GetTrack(idxPosTrk[1]));
  } else if (nGoodPosTrk == 0 && nGoodNegTrk == 2) { //two positive  tracks
    Track1 = static_cast<AliAODTrack*>(fAOD->GetTrack(idxNegTrk[0]));
    Track2 = static_cast<AliAODTrack*>(fAOD->GetTrack(idxNegTrk[1]));
  }
  iSelectionCounter++; // os muons (6)




  fCosThetaHE   = -999;
  fCosThetaCS   = -999;
  fPhiHE        = -999;
  fPhiCS        = -999;
  fTildePhiHEpos  = -999;
  fTildePhiHEneg  = -999;
  fTildePhiCSpos  = -999;
  fTildePhiCSneg  = -999;

  TrkTrkKine(Track1,Track2,gMuonMass);

  // clean up
  delete [] idxNegTrk;
  delete [] idxPosTrk;












  // V0
  AliVVZERO *dataVZERO = dynamic_cast<AliVVZERO*>(fAOD->GetVZEROData());
  if(!dataVZERO) {
    fAnaTree->Fill();
    PostData(1, fAnaTree);
    PostData(2, fOutputList);
    return;
  }
  fV0ADecision = dataVZERO->GetV0ADecision();
  fV0CDecision = dataVZERO->GetV0CDecision();
  // check AD
  AliVAD *dataAD = dynamic_cast<AliVAD*>(fAOD->GetADData());
  if(!dataAD){
    fAnaTree->Fill();
    PostData(1, fAnaTree);
    PostData(2, fOutputList);
    return;
  }
  fADADecision = dataAD->GetADADecision();
  fADCDecision = dataAD->GetADCDecision();






  fV0TotalNCells = 0;
  for(Int_t iV0Hits = 0; iV0Hits < 64; iV0Hits++) {
        fV0Hits[iV0Hits] = dataVZERO->GetBBFlag(iV0Hits);
        if(fV0Hits[iV0Hits] == kTRUE) {
              fV0TotalNCells += 1;
        }
  }





  if( (fV0TotalNCells > 2) ||
      (fADADecision != 0)  ||
      (fADCDecision != 0)  ||
      (fV0ADecision != 0)  ||
      !(fV0CDecision == 0 || fV0CDecision == 1) ){
    fAnaTree->Fill();
    PostData(1, fAnaTree);
    PostData(2, fOutputList);
    return;
  }
  fCounterH->Fill(iSelectionCounter); // after all cuts (7)
  iSelectionCounter++;




  // fill the tree
  fAnaTree->Fill();

  // post the data
  PostData(1, fAnaTree);
  PostData(2, fOutputList);

}
//_____________________________________________________________________________
/* - The following are code snippets adapted from the AliAODDimuon class.
   - The problem is that that class was adapted specifically for the
   - inclusive people's analysis, hence it is not fit for the UPC...
   -
 */
Double_t AliAnalysisTaskNanoJPsi2016Fwd::CosThetaCollinsSoper( TLorentzVector muonPositive,
                                                            TLorentzVector muonNegative,
                                                            TLorentzVector possibleJPsi )
{
  /* - This function computes the Collins-Soper cos(theta) for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  /* - Determine the z axis for the CS angle.
     -
   */
  TVector3 zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  /* - Determine the CS angle (angle between mu+ and the z axis defined above)
     -
   */
  Double_t CosThetaCS = zaxisCS.Dot((pMu1Dimu.Vect()).Unit());
  return   CosThetaCS;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskNanoJPsi2016Fwd::CosThetaHelicityFrame( TLorentzVector muonPositive,
                                                             TLorentzVector muonNegative,
                                                             TLorentzVector possibleJPsi )
{
  /* - This function computes the Helicity cos(theta) for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the calculation of the polarization angle
  // (i.e. the direction of the dimuon in the CM system)
  //
  TVector3 zaxis = (possibleJPsi.Vect()).Unit();
  /* - Determine the He angle (angle between mu+ and the z axis defined above)
     -
   */
  Double_t CosThetaHE = zaxis.Dot((pMu1Dimu.Vect()).Unit());
  return   CosThetaHE;

}
//_____________________________________________________________________________
Double_t AliAnalysisTaskNanoJPsi2016Fwd::CosPhiCollinsSoper( TLorentzVector muonPositive,
                                                          TLorentzVector muonNegative,
                                                          TLorentzVector possibleJPsi )
{
  /* - This function computes the Collins-Soper PHI for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  /* - Determine the z axis for the CS angle.
     -
   */
  TVector3 zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  //
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  //
  TVector3 yaxisCS=(((pProjDimu.Vect()).Unit()).Cross((pTargDimu.Vect()).Unit())).Unit();
  TVector3 xaxisCS=(yaxisCS.Cross(zaxisCS)).Unit();

  Double_t phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxisCS),((pMu1Dimu.Vect()).Dot(xaxisCS)));
  return   phi;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskNanoJPsi2016Fwd::CosPhiHelicityFrame(  TLorentzVector muonPositive,
                                                            TLorentzVector muonNegative,
                                                            TLorentzVector possibleJPsi )
{
  /* - This function computes the helicity phi for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
  */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the calculation of the polarization angle
  // (i.e. the direction of the dimuon in the CM system)
  //
  TVector3 zaxis = (possibleJPsi.Vect()).Unit();
  TVector3 yaxis = ((pProjDimu.Vect()).Cross(pTargDimu.Vect())).Unit();
  TVector3 xaxis = (yaxis.Cross(zaxis)).Unit();
  //
  // --- Calculation of the azimuthal angle (Helicity)
  //
  Double_t phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxis),(pMu1Dimu.Vect()).Dot(xaxis));
  return   phi;
}
//_____________________________________________________________________________
void AliAnalysisTaskNanoJPsi2016Fwd::Terminate(Option_t *)
{
    cout << endl;
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
