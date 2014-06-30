#include "GenVisualizer.h"
#include "Utilities.h"

#include "TStyle.h"
#include "TMath.h"
#include "TVector2.h"
#include "TLatex.h"
#include "TColor.h"

namespace susy {

  GenVisualizer::GenVisualizer() :
    plotsDir("."),
    canvas_(new TCanvas("GenVisualizer", "GenVisualizer", 500, 500)),
    frame_(new TH2F("frame", "", 100, -3., 3., 100, -TMath::Pi(), TMath::Pi())),
    wMarker_(0., 0., 27),
    connector_(),
    jetCone_(0., 0., 0.5),
    puJetCone_(0., 0., 0.5),
    genJetCone_(0., 0., 0.5),
    metLine_(-3., 0., 3., 0.),
    sumPtLine_(-3., 0., 3., 0.),
    legend_(0.8, 0.9, 0.9, 1.),
    maxPt_(150.),
    nColors_(255),
    col0_(0),
    showGenJet_(true),
    showRecoJet_(true)
  {
    frame_->SetDirectory(0);
  
    gStyle->SetNumberContours(nColors_);
    TColor::InitializeColors();
    const Int_t nRGBs = 5;
    Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    col0_ = TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nColors_);

    frame_->SetStats(0);
    frame_->SetEntries(1);
    frame_->SetBinContent(1, 1, -1.);
    frame_->SetMinimum(0.);
    frame_->SetMaximum(maxPt_);

    marker_[kPhoton].SetMarkerStyle(3);
    marker_[kMuon].SetMarkerStyle(2);
    marker_[kElectron].SetMarkerStyle(4);
    marker_[kHadron].SetMarkerStyle(7);
    marker_[kNeutrino].SetMarkerStyle(22);
    marker_[kLSP].SetMarkerStyle(21);

    wMarker_.SetMarkerSize(2);
    connector_.SetLineStyle(kDashed);
    jetCone_.SetFillStyle(0);
    puJetCone_.SetFillStyle(0);
    puJetCone_.SetLineStyle(kDotted);
    genJetCone_.SetFillStyle(0);
    genJetCone_.SetLineStyle(9);
    genJetCone_.SetLineWidth(2);
    metLine_.SetLineWidth(2);
    sumPtLine_.SetLineStyle(kDotted);

    legend_.SetFillStyle(0);
    legend_.SetBorderSize(0);
    legend_.AddEntry(&marker_[kPhoton], "#gamma", "p");
    legend_.AddEntry(&marker_[kMuon], "#mu", "p");
    legend_.AddEntry(&marker_[kElectron], "e", "p");
    legend_.AddEntry(&marker_[kHadron], "h", "p");
    legend_.AddEntry(&marker_[kNeutrino], "#nu", "p");
  }

  GenVisualizer::~GenVisualizer()
  {
    delete canvas_;
    delete frame_;
  }

  int
  GenVisualizer::getPtColor_(double _pt)
  {
    int iCol = _pt / maxPt_ * nColors_;
    if(iCol >= nColors_) iCol = nColors_ - 1;
    return col0_ + iCol;
  }

  unsigned
  GenVisualizer::formGenJets_(int* _nConst, double* _pt, double* _eta, double* _phi)
  {
    bool used[256];
    int nConst[256];
    double pt[256];
    double eta[256];
    double phi[256];

    unsigned nJets(0);

    for(unsigned iG = 0; iG != genSize_; ++iG){
      if(genStatus_[iG] != 1)
        used[iG] = true;
      else{
        switch(TMath::Abs(genPdgId_[iG])){
        case 12:
        case 14:
        case 16:
        case 1000022:
          used[iG] = true;
          break;
        default:
          used[iG] = false;
          pt[iG] = genPt_[iG];
          eta[iG] = genEta_[iG];
          phi[iG] = genPhi_[iG];
          break;
        }

        nConst[iG] = 1;
      }
    }

    while(true){
      double di[256];
      double dmin = TMath::Infinity();
      int minPair[2] = {-1, -1};
      for(unsigned iG = 0; iG != genSize_; ++iG){
        if(used[iG]) continue;

        di[iG] = 1. / pt[iG] / pt[iG];
        if(di[iG] < dmin){
          minPair[0] = iG;
          minPair[1] = -1;
          dmin = di[iG];
        }
        
        for(unsigned iG2 = 0; iG2 != iG; ++iG2){
          if(used[iG2]) continue;

          double dEta = eta[iG] - eta[iG2];
          double dPhi = TVector2::Phi_mpi_pi(phi[iG] - phi[iG2]);
          double dij = TMath::Min(di[iG], di[iG2]) * (dEta * dEta + dPhi * dPhi) / 0.5 / 0.5;
          
          if(dij < dmin){
            minPair[0] = iG;
            minPair[1] = iG2;
            dmin = dij;
          }
        }
      }

      if(minPair[0] == -1) break;

      if(minPair[1] == -1){
        unsigned iG = minPair[0];

        used[iG] = true;

        _nConst[nJets] = nConst[iG];
        _pt[nJets] = pt[iG];
        _eta[nJets] = eta[iG];
        _phi[nJets] = phi[iG];

        ++nJets;
      }
      else{
        unsigned iG = minPair[0];
        unsigned iG2 = minPair[1];

        used[iG2] = true;

        nConst[iG] += nConst[iG2];

        double sumPt = pt[iG] + pt[iG2];
        eta[iG] = (eta[iG] * pt[iG] + eta[iG2] * pt[iG2]) / sumPt;
        phi[iG] = TVector2::Phi_mpi_pi(phi[iG] + TVector2::Phi_mpi_pi(phi[iG2] - phi[iG]) * pt[iG] / sumPt);
        pt[iG] = sumPt;
      }
    }

    return nJets;
  }

  void
  GenVisualizer::setAddress(TTree& _eventTree)
  {
    _eventTree.SetBranchAddress("eventNumber", &eventNumber_);
    _eventTree.SetBranchAddress("gen.size", &genSize_);
    _eventTree.SetBranchAddress("gen.pt", genPt_);
    _eventTree.SetBranchAddress("gen.eta", genEta_);
    _eventTree.SetBranchAddress("gen.phi", genPhi_);
    _eventTree.SetBranchAddress("gen.status", genStatus_);
    _eventTree.SetBranchAddress("gen.pdgId", genPdgId_);
    _eventTree.SetBranchAddress("gen.motherIndex", genMotherIndex_);
    _eventTree.SetBranchAddress("met", &met_);
    _eventTree.SetBranchAddress("metPhi", &metPhi_);

    _eventTree.SetBranchAddress("jet.size", &jetSize_);
    _eventTree.SetBranchAddress("jet.pt", jetPt_);
    _eventTree.SetBranchAddress("jet.eta", jetEta_);
    _eventTree.SetBranchAddress("jet.phi", jetPhi_);
    _eventTree.SetBranchAddress("jet.passPUJetIdLoose", jetPassPUJetIdLoose_);
    _eventTree.SetBranchAddress("jet.isLoose", jetIsLoose_);
  }

  void
  GenVisualizer::process(TString const& _plotName/* = ""*/)
  {
    frame_->SetTitle(TString::Format("Event %d", eventNumber_));
    frame_->Draw("colz");
    
    for(unsigned iG = 0; iG != genSize_; ++iG){
      if(TMath::Abs(genPdgId_[iG]) == 24){
        wMarker_.SetX(genEta_[iG]);
        wMarker_.SetY(genPhi_[iG]);
        wMarker_.Draw();
      }
      
      if(genStatus_[iG] != 1) continue;

      unsigned pType = nPType;

      switch(TMath::Abs(genPdgId_[iG])){
      case 22:
        pType = kPhoton;
        break;
      case 13:
        pType = kMuon;
        break;
      case 11:
        pType = kElectron;
        break;
      case 12:
      case 14:
      case 16:
        pType = kNeutrino;
        break;
      case 1000022:
        pType = kLSP;
        break;
      default:
        pType = kHadron;
        break;
      }

      marker_[pType].SetMarkerColor(getPtColor_(genPt_[iG]));

      marker_[pType].DrawMarker(genEta_[iG], genPhi_[iG]);

      if(pType != kHadron){
        short idx = genMotherIndex_[iG];
        while(idx != -1 && TMath::Abs(genPdgId_[idx]) != 24) idx = genMotherIndex_[idx];
        if(idx != -1) connector_.DrawLine(wMarker_.GetX(), wMarker_.GetY(), genEta_[iG], genPhi_[iG]);
      }
    }

    if(showGenJet_){
      int nConst[256];
      double pt[256];
      double eta[256];
      double phi[256];

      unsigned nJets(formGenJets_(nConst, pt, eta, phi));

      for(unsigned iJ(0); iJ != nJets; ++iJ){
        if(nConst[iJ] > 1 && pt[iJ] > 20.){
          genJetCone_.SetLineColor(getPtColor_(pt[iJ]));
          genJetCone_.DrawEllipse(eta[iJ], phi[iJ], 0.5, 0.5, 0., 360., 0., "");
        }
      }
    }

    if(showRecoJet_){
      for(unsigned iJ(0); iJ != jetSize_; ++iJ){
        if(!jetIsLoose_[iJ]) continue;

        if(jetPassPUJetIdLoose_[iJ]){
          jetCone_.SetLineColor(getPtColor_(jetPt_[iJ]));
          jetCone_.DrawEllipse(jetEta_[iJ], jetPhi_[iJ], 0.5, 0.5, 0., 360., 0., "");
        }
        else{
          puJetCone_.SetLineColor(getPtColor_(jetPt_[iJ]));
          puJetCone_.DrawEllipse(jetEta_[iJ], jetPhi_[iJ], 0.5, 0.5, 0., 360., 0., "");
        }
      }
    }

    metLine_.SetLineColor(getPtColor_(met_));
    metLine_.SetY1(metPhi_);
    metLine_.SetY2(metPhi_);
    metLine_.Draw();

    sumPtLine_.SetY1(TVector2::Phi_mpi_pi(metPhi_ + TMath::Pi()));
    sumPtLine_.SetY2(TVector2::Phi_mpi_pi(metPhi_ + TMath::Pi()));
    sumPtLine_.Draw();

    if(_plotName != "") print(_plotName);
  }

  void
  GenVisualizer::print(TString const& _plotName)
  {
    legend_.Draw();
    canvas_->Print(plotsDir + "/" + _plotName);
  }

  void
  GenVisualizer::zoomIn(double _eta, double _phi, double _r)
  {
    frame_->GetXaxis()->SetRangeUser(_eta - _r, _eta + _r);
    frame_->GetYaxis()->SetRangeUser(_phi - _r, _phi + _r);
  }

  void
  GenVisualizer::zoomOut()
  {
    frame_->GetXaxis()->SetRangeUser(-3., 3.);
    frame_->GetYaxis()->SetRangeUser(-TMath::Pi(), TMath::Pi());
  }

  void
  GenVisualizer::labelParticles()
  {
    TLatex label;
    label.SetTextSize(0.04);

    double etaMin(frame_->GetXaxis()->GetBinLowEdge(frame_->GetXaxis()->GetFirst()));
    double etaMax(frame_->GetXaxis()->GetBinUpEdge(frame_->GetXaxis()->GetLast()));
    double phiMin(frame_->GetYaxis()->GetBinLowEdge(frame_->GetYaxis()->GetFirst()));
    double phiMax(frame_->GetYaxis()->GetBinUpEdge(frame_->GetYaxis()->GetLast()));

    double etaOffset((etaMax - etaMin) * 0.03);
    double phiOffset((phiMax - phiMin) * 0.03);

    for(unsigned iG = 0; iG != genSize_; ++iG){
      if(genStatus_[iG] != 1) continue;
      if(genEta_[iG] < etaMin || genEta_[iG] > etaMax) continue;
      if(genPhi_[iG] < phiMin || genPhi_[iG] > phiMax) continue;

      label.DrawLatex(genEta_[iG] + etaOffset, genPhi_[iG] + phiOffset, TString::Format("#splitline{%s}{%.0f}", particleName(genPdgId_[iG], true).Data(), genPt_[iG]));
    }
  }

  void
  GenVisualizer::setMaxPt(double _maxPt)
  {
    maxPt_ = _maxPt;
    frame_->SetMaximum(maxPt_);
  }

}
