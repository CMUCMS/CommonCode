/* Partially auto-generated source file - edit where indicated */
/* Add necessary inclusions below */
#include "ObjectVars.h"
#include "Utilities.h"
#include "ObjectSelector.h"
#include "TFile.h"
#include <stdexcept>
#include <limits>


namespace susy {

  PhotonVars::PhotonVars() :
    pt(0.),
    eta(0.),
    phi(0.),
    px(0.),
    py(0.),
    pz(0.),
    energy(0.),
    hOverE(0.),
    sigmaIetaIeta(0.),
    sigmaIphiIphi(0.),
    etaWidth(0.),
    phiWidth(0.),
    r9(0.),
    r5(0.),
    trackerIso(0.),
    ecalIso(0.),
    hcalIso(0.),
    chargedHadronIso(0.),
    neutralHadronIso(0.),
    photonIso(0.),
    caloX(0.),
    caloY(0.),
    caloZ(0.),
    iSubdet(0),
    superClusterIndex(0),
    nPixelSeeds(0),
    nClusters(0),
    hasMatchedElectron(false),
    electronVetoBit(false),
    looseElectronVetoBit(false),
    isLoose(false),
    isMedium(false),
    isTight(false),
    isLoosePix(false),
    isMediumPix(false),
    isTightPix(false),
    isLooseLV(false),
    isMediumLV(false),
    isTightLV(false)
  {
  }

  ElectronVars::ElectronVars() :
    pt(0.),
    eta(0.),
    phi(0.),
    px(0.),
    py(0.),
    pz(0.),
    energy(0.),
    combRelSubdetIso(0.),
    combRelIso(0.),
    deltaEta(0.),
    deltaPhi(0.),
    sigmaIetaIeta(0.),
    sigmaIphiIphi(0.),
    r9(0.),
    r5(0.),
    etaWidth(0.),
    phiWidth(0.),
    hOverE(0.),
    d0(0.),
    dz(0.),
    epDiff(0.),
    vtxFitProb(0.),
    dCot(0.),
    dist(0.),
    caloX(0.),
    caloY(0.),
    caloZ(0.),
    iSubdet(0),
    superClusterIndex(0),
    nClusters(0),
    nPixelHits(0),
    nMissingHits(0),
    passConversionVeto(false),
    isVeto(false),
    isLoose(false),
    isMedium(false),
    isTight(false)
  {
  }

  MuonVars::MuonVars() :
    pt(0.),
    eta(0.),
    phi(0.),
    px(0.),
    py(0.),
    pz(0.),
    energy(0.),
    normChi2(0.),
    dxy(0.),
    dz(0.),
    combRelSubdetIso(0.),
    combRelIso(0.),
    iSubdet(0),
    nMatchedStations(0),
    nLayersWithMmt(0),
    nValidMuonHits(0),
    nValidPixelHits(0),
    isGlobalMuon(false),
    isPFMuon(false),
    hasInnerTrack(false),
    hasGlobalTrack(false),
    hasBestTrack(false),
    isLoose(false),
    isTight(false)
  {
  }

  JetVars::JetVars() :
    pt(0.),
    eta(0.),
    phi(0.),
    px(0.),
    py(0.),
    pz(0.),
    energy(0.),
    jecScale(0.),
    chFraction(0.),
    nhFraction(0.),
    ceFraction(0.),
    neFraction(0.),
    iSubdet(0),
    nConstituents(0),
    nCharged(0),
    isLoose(false)
  {
  }

  VertexVars::VertexVars() :
    x(0.),
    y(0.),
    z(0.),
    rho(0.),
    sumPt2(0.),
    chi2(0.),
    ndof(0.),
    isGood(false)
  {
  }

/* START USER-DEFINED IMPLEMENTATION (DO NOT MODIFY THIS LINE) */

  void
  photonEffectiveAreas(double _eta, double* _effA)
  {
    double& effATrk(_effA[0]);
    double& effAEcal(_effA[1]);
    double& effAHcal(_effA[2]);

    if(_eta < etaGapBegin){
      effATrk = 0.167;
      effAEcal = 0.183;
      effAHcal = 0.062;
    }
    else{
      effATrk = 0.032;
      effAEcal = 0.090;
      effAHcal = 0.180;
    }

    double& effACH(_effA[3]);
    double& effANH(_effA[4]);
    double& effAPh(_effA[5]);

    // CutBasedPhotonID2012
    if(_eta < 1.){
      effACH = 0.012;
      effANH = 0.03;
      effAPh = 0.148;
    }
    else if(_eta < 1.479){
      effACH = 0.010;
      effANH = 0.057;
      effAPh = 0.13;
    }
    else if(_eta < 2.){
      effACH = 0.014;
      effANH = 0.039;
      effAPh = 0.112;
    }
    else if(_eta < 2.2){
      effACH = 0.012;
      effANH = 0.015;
      effAPh = 0.216;
    }
    else if(_eta < 2.3){
      effACH = 0.016;
      effANH = 0.024;
      effAPh = 0.262;
    }
    else if(_eta < 2.4){
      effACH = 0.02;
      effANH = 0.039;
      effAPh = 0.26;
    }
    else{
      effACH = 0.012;
      effANH = 0.072;
      effAPh = 0.266;
    }
  }

  void
  PhotonVars::set(Photon const& _ph, Event const& _event)
  {
    pt = _ph.momentum.Pt();
    eta = _ph.momentum.Eta();
    phi = _ph.momentum.Phi();
    px = _ph.momentum.X();
    py = _ph.momentum.Y();
    pz = _ph.momentum.Z();
    energy = _ph.momentum.E();

    double absEta(std::abs(_ph.caloPosition.Eta()));
    if(absEta < etaGapBegin) iSubdet = 0;
    else if(absEta < etaGapEnd) iSubdet = -1;
    else if(absEta < etaMax) iSubdet = 1;
    else iSubdet = -1;

    hOverE = _ph.hadTowOverEm;

    double effA[6]; // tracker, ecal, hcal, chargedHadron, neutralHadron, photon
    photonEffectiveAreas(absEta, effA);

    // CutBasedPhotonID2012#Effective_Areas_for_rho_correcti
    // "For 52X: double_kt6PFJets_rho_RECO"
    double rho(_event.rho);

    trackerIso = _ph.trkSumPtHollowConeDR04 - rho * effA[0] - 0.001 * pt;

    ecalIso = _ph.ecalRecHitSumEtConeDR04 - rho * effA[1] - 0.006 * pt;

    hcalIso = _ph.hcalDepth1TowerSumEtConeDR03 + _ph.hcalDepth2TowerSumEtConeDR03 - rho * effA[2] - 0.0025 * pt;

    sigmaIetaIeta = _ph.sigmaIetaIeta;

    sigmaIphiIphi = _ph.sigmaIphiIphi;

    r9 = _ph.r9;

    if(!_ph.superCluster)
      throw std::runtime_error("PhotonWithNoCluster");

    SuperCluster const& sc(*_ph.superCluster);

    superClusterIndex = _ph.superClusterIndex;

    caloX = _ph.caloPosition.X();
    caloY = _ph.caloPosition.Y();
    caloZ = _ph.caloPosition.Z();

    r5 = _ph.e1x5 / sc.energy;

    etaWidth = sc.etaWidth;

    phiWidth = sc.phiWidth;

    chargedHadronIso = _ph.chargedHadronIso - rho * effA[3];

    neutralHadronIso = _ph.neutralHadronIso - rho * effA[4] - 0.04 * pt;

    photonIso = _ph.photonIso - rho * effA[5] - 0.005 * pt;

    nPixelSeeds = _ph.nPixelSeeds;

    nClusters = sc.basicClusterIndices.size();

    electronVetoBit = _ph.passelectronveto;

    // Reproducing hasMatchedPromptElectron() implementation in RecoEgamma/EgammaTools/src/ConversionTools.cc
    // but allowing the electron to miss 1 hit
    typename ElectronCollectionMap::const_iterator electronsSrc(_event.electrons.find("gsfElectrons"));
    if(electronsSrc == _event.electrons.end())
      throw std::runtime_error("GsfElectrons not in event");

    hasMatchedElectron = false;
    ElectronCollection const& electrons(electronsSrc->second);
    unsigned nE(electrons.size());
    unsigned iE(0);
    for(; iE != nE; ++iE){
      Electron const& el(electrons[iE]);
      if(el.superClusterIndex != _ph.superClusterIndex) continue;
      hasMatchedElectron = true;
      if(el.nMissingHits > 1 || !el.passConversionVeto) iE = nE; // this is not a prompt electron
      break;
    }
    looseElectronVetoBit = (iE == nE);

    isLoose = ObjectSelector::isGoodPhoton(*this, PhLoose12);
    isMedium = ObjectSelector::isGoodPhoton(*this, PhMedium12);
    isTight = ObjectSelector::isGoodPhoton(*this, PhTight12);
    isLoosePix = ObjectSelector::isGoodPhoton(*this, PhLoose12Pix);
    isMediumPix = ObjectSelector::isGoodPhoton(*this, PhMedium12Pix);
    isTightPix = ObjectSelector::isGoodPhoton(*this, PhTight12Pix);
    isLooseLV = ObjectSelector::isGoodPhoton(*this, PhLoose12LV);
    isMediumLV = ObjectSelector::isGoodPhoton(*this, PhMedium12LV);
    isTightLV = ObjectSelector::isGoodPhoton(*this, PhTight12LV);
  }

  /*static*/
  void
  PhotonVars::setBranchStatus(TTree& _tree)
  {
    _tree.SetBranchStatus("photons_photons*", 1);
    _tree.SetBranchStatus("electrons_gsfElectrons*", 1);
    _tree.SetBranchStatus("superClusters*", 1);
    _tree.SetBranchStatus("rho", 1);
  }

  void
  electronEffectiveAreas(double _eta, double &_effA)
  {
    //EgammaEARhoCorrection
    if(_eta < 1.)
      _effA = 0.13;
    else if(_eta < 1.479)
      _effA = 0.14;
    else if(_eta < 2.)
      _effA = 0.07;
    else if(_eta < 2.2)
      _effA = 0.09;
    else if(_eta < 2.3)
      _effA = 0.11;
    else if(_eta < 2.4)
      _effA = 0.11;
    else
      _effA = 0.14;
  }

  void
  ElectronVars::set(Electron const& _el, Event const& _event)
  {
    pt = _el.momentum.Pt();
    eta = _el.momentum.Eta();
    phi = _el.momentum.Phi();
    px = _el.momentum.X();
    py = _el.momentum.Y();
    pz = _el.momentum.Z();
    energy = _el.momentum.E();

    if(!_el.superCluster)
      throw std::runtime_error("ElectronWithNoCluster");

    SuperCluster const& sc(*_el.superCluster);

    superClusterIndex = _el.superClusterIndex;

    double absEta(std::abs(sc.position.Eta()));

    if(absEta < etaGapBegin) iSubdet = 0;
    else if(absEta < etaGapEnd) iSubdet = -1;
    else if(absEta < etaMax) iSubdet = 1;
    else iSubdet = -1;

    if(!_el.gsfTrack)
      throw std::runtime_error("ElectronWithNoTrack");

    // TODO: Not ideal to do this for every single Electron and Muon
    unsigned nV(_event.vertices.size());
    unsigned iV(0);
    while(iV != nV && !ObjectSelector::isGoodVertex(_event.vertices[iV])) ++iV;
    if(iV == nV)
      throw std::runtime_error("Event with no good vertex");

    Vertex const& primVtx(_event.vertices[iV]);

    // EgammaEARhoCorrection #Rho for 2012-Effective Areas
    // "In this case the rho double_kt6PFJets_rho_RECO already saved in the event (since CMSSW_5XY) needs to be used."
    double rho(_event.rho);

    double effA(0.);
    electronEffectiveAreas(absEta, effA);

    combRelSubdetIso = (std::max(0., _el.dr03EcalRecHitSumEt - 1.) + _el.dr03HcalDepth1TowerSumEt + _el.dr03HcalDepth2TowerSumEt + _el.dr03TkSumPt) / pt;

    combRelIso = (_el.chargedHadronIso + std::max(0., _el.neutralHadronIso + _el.photonIso - rho * effA)) / pt;

    deltaEta = std::abs(_el.deltaEtaSuperClusterTrackAtVtx);

    deltaPhi = std::abs(_el.deltaPhiSuperClusterTrackAtVtx);

    sigmaIetaIeta = _el.sigmaIetaIeta;

    sigmaIphiIphi = _el.sigmaIphiIphi;

    r9 = _el.r9;

    r5 = _el.e1x5 / sc.energy;

    etaWidth = sc.etaWidth;

    phiWidth = sc.phiWidth;

    hOverE = _el.hcalOverEcalBc;

    caloX = sc.position.X();
    caloY = sc.position.Y();
    caloZ = sc.position.Z();

    d0 = std::abs(_el.gsfTrack->dxy(primVtx.position));

    dz = std::abs(_el.gsfTrack->dz(primVtx.position));

    epDiff = std::abs(1. / _el.ecalEnergy - 1. / (_el.ecalEnergy / _el.eSuperClusterOverP));

    dCot = std::abs(_el.convDcot);

    dist = std::abs(_el.convDist);

    nClusters = sc.basicClusterIndices.size();

    nPixelHits = _el.gsfTrack->numberOfValidPixelHits;

    nMissingHits = _el.nMissingHits;

    passConversionVeto = _el.passConversionVeto;

    isVeto = ObjectSelector::isGoodElectron(*this, ElVeto12);
    isLoose = ObjectSelector::isGoodElectron(*this, ElLoose12);
    isMedium = ObjectSelector::isGoodElectron(*this, ElMedium12);
    isTight = ObjectSelector::isGoodElectron(*this, ElTight12);
  }

  /*static*/
  void
  ElectronVars::setBranchStatus(TTree& _tree)
  {
    _tree.SetBranchStatus("electrons_gsfElectrons*", 1);
    _tree.SetBranchStatus("tracks*", 1);
    _tree.SetBranchStatus("vertices*", 1);
    _tree.SetBranchStatus("superClusters*", 1);
    _tree.SetBranchStatus("rho", 1);
  }

  void
  MuonVars::set(Muon const& _mu, Event const& _event)
  {
    pt = _mu.momentum.Pt();
    eta = _mu.momentum.Eta();
    phi = _mu.momentum.Phi();
    px = _mu.momentum.X();
    py = _mu.momentum.Y();
    pz = _mu.momentum.Z();
    energy = _mu.momentum.E();

    double absEta(std::abs(eta));

    if(absEta < etaGapBegin) iSubdet = 0;
    else if(absEta < etaGapEnd) iSubdet = -1;
    else if(absEta < etaMax) iSubdet = 1;
    else iSubdet = -1;

    isGlobalMuon = _mu.isGlobalMuon();

    isPFMuon = _mu.isPFMuon();

    hasInnerTrack = _mu.innerTrack;

    hasGlobalTrack = _mu.globalTrack;

    Track const* bestTrack(pt > 200. ? _mu.highPtBestTrack : _mu.bestTrack);
    hasBestTrack = (bestTrack != 0);

    // TODO: Not ideal to do this for every single Electron and Muon
    unsigned nV(_event.vertices.size());
    unsigned iV(0);
    while(iV != nV && !ObjectSelector::isGoodVertex(_event.vertices[iV])) ++iV;
    if(iV == nV)
      throw std::runtime_error("Event with no good vertex");

    Vertex const& primVtx(_event.vertices[iV]);

    nMatchedStations = _mu.nMatchedStations;

    nLayersWithMmt = _mu.nPixelLayersWithMeasurement + _mu.nStripLayersWithMeasurement;

    if(hasGlobalTrack == 1) normChi2 = _mu.globalTrack->normChi2();
    else normChi2 = -1.;

    if(hasGlobalTrack == 1) nValidMuonHits = _mu.globalTrack->numberOfValidMuonHits;
    else nValidMuonHits = 0;

    if(hasBestTrack == 1) dxy = std::abs(bestTrack->dxy(primVtx.position));
    else if(hasInnerTrack) dxy = std::abs(_mu.innerTrack->dxy(primVtx.position));
    else dxy = -1.;

    if(hasBestTrack == 1) dz = std::abs(bestTrack->dz(primVtx.position));
    else if(hasInnerTrack == 1) dz = std::abs(_mu.innerTrack->dz(primVtx.position));
    else dz = -1.;

    if(hasInnerTrack == 1) nValidPixelHits = _mu.innerTrack->numberOfValidPixelHits;
    else nValidPixelHits = 0;

    combRelSubdetIso = (_mu.ecalIsoR03 + _mu.hcalIsoR03 + _mu.trackIsoR03) / pt; 

    combRelIso = (_mu.sumChargedHadronPt04 + std::max(0., _mu.sumNeutralHadronEt04 + _mu.sumPhotonEt04 - 0.5 * _mu.sumPUPt04)) / pt;

    isLoose = ObjectSelector::isGoodMuon(*this, MuLoose12);
    isTight = ObjectSelector::isGoodMuon(*this, MuTight12);
  }

  /*static*/
  void
  MuonVars::setBranchStatus(TTree& _tree)
  {
    _tree.SetBranchStatus("muons_muons*", 1);
    _tree.SetBranchStatus("tracks*", 1);
    _tree.SetBranchStatus("vertices*", 1);
  }

  void
  JetVars::set(PFJet const& _jet, Event const& _event)
  {
    jecScale = _jet.jecScaleFactors.find("L1FastL2L3")->second;

    TLorentzVector corrP(_jet.momentum * jecScale);

    pt = corrP.Pt();
    if(pt > 0){
      eta = corrP.Eta();
      phi = corrP.Phi();
    }
    else{
      eta = (corrP.Z() > 0. ? 1. : -1.) * std::numeric_limits<float>::max();
      phi = 0.;
    }
    px = corrP.X();
    py = corrP.Y();
    pz = corrP.Z();
    energy = corrP.E();

    double absEta(std::abs(eta));

    if(absEta < etaGapBegin) iSubdet = 0;
    else if(absEta < etaGapEnd) iSubdet = -1;
    else if(absEta < etaMax) iSubdet = 1;
    else iSubdet = -1;

    chFraction = _jet.chargedHadronEnergy / energy;

    nhFraction = _jet.neutralHadronEnergy / energy;

    ceFraction = _jet.chargedEmEnergy / energy;

    neFraction = _jet.neutralEmEnergy / energy;

    nConstituents = _jet.nConstituents;

    nCharged = _jet.chargedMultiplicity;

    isLoose = ObjectSelector::isGoodJet(*this, JtLoose);
  }

  /*static*/
  void
  JetVars::setBranchStatus(TTree& _tree)
  {
    _tree.SetBranchStatus("pfJets_ak5*", 1);
    _tree.SetBranchStatus("pfJets_ak5chs*", 0);
  }

  void
  VertexVars::set(Vertex const& _vtx)
  {
    x = _vtx.position.X();
    y = _vtx.position.Y();
    z = _vtx.position.Z();
    rho = _vtx.position.Perp();

    sumPt2 = _vtx.sumPt2;
    chi2 = _vtx.chi2;
    ndof = _vtx.ndof;

    isGood = ObjectSelector::isGoodVertex(*this);
  }

  /*static*/
  void
  VertexVars::setBranchStatus(TTree& _tree)
  {
    _tree.SetBranchStatus("vertices*", 1);
  }

/* END USER-DEFINED IMPLEMENTATION (DO NOT MODIFY THIS LINE) */

}
