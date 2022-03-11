/* Partially auto-generated source file - edit where indicated */
/* Add necessary inclusions below */
#include "ObjectVars.h"
#include "Utilities.h"
#include "ObjectSelector.h"
#include "TFile.h"
#include <limits>
#ifndef STANDALONE
#include "SusyEvent.h"
#endif

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

  void
  PhotonVars::setBranches(TTree& _tree)
  {
    _tree.Branch("photon.pt", &pt, "pt/F");
    _tree.Branch("photon.eta", &eta, "eta/F");
    _tree.Branch("photon.phi", &phi, "phi/F");
    _tree.Branch("photon.px", &px, "px/F");
    _tree.Branch("photon.py", &py, "py/F");
    _tree.Branch("photon.pz", &pz, "pz/F");
    _tree.Branch("photon.energy", &energy, "energy/F");
    _tree.Branch("photon.hOverE", &hOverE, "hOverE/F");
    _tree.Branch("photon.sigmaIetaIeta", &sigmaIetaIeta, "sigmaIetaIeta/F");
    _tree.Branch("photon.sigmaIphiIphi", &sigmaIphiIphi, "sigmaIphiIphi/F");
    _tree.Branch("photon.etaWidth", &etaWidth, "etaWidth/F");
    _tree.Branch("photon.phiWidth", &phiWidth, "phiWidth/F");
    _tree.Branch("photon.r9", &r9, "r9/F");
    _tree.Branch("photon.r5", &r5, "r5/F");
    _tree.Branch("photon.trackerIso", &trackerIso, "trackerIso/F");
    _tree.Branch("photon.ecalIso", &ecalIso, "ecalIso/F");
    _tree.Branch("photon.hcalIso", &hcalIso, "hcalIso/F");
    _tree.Branch("photon.chargedHadronIso", &chargedHadronIso, "chargedHadronIso/F");
    _tree.Branch("photon.neutralHadronIso", &neutralHadronIso, "neutralHadronIso/F");
    _tree.Branch("photon.photonIso", &photonIso, "photonIso/F");
    _tree.Branch("photon.caloX", &caloX, "caloX/F");
    _tree.Branch("photon.caloY", &caloY, "caloY/F");
    _tree.Branch("photon.caloZ", &caloZ, "caloZ/F");
    _tree.Branch("photon.iSubdet", &iSubdet, "iSubdet/S");
    _tree.Branch("photon.superClusterIndex", &superClusterIndex, "superClusterIndex/S");
    _tree.Branch("photon.nPixelSeeds", &nPixelSeeds, "nPixelSeeds/b");
    _tree.Branch("photon.nClusters", &nClusters, "nClusters/b");
    _tree.Branch("photon.hasMatchedElectron", &hasMatchedElectron, "hasMatchedElectron/O");
    _tree.Branch("photon.electronVetoBit", &electronVetoBit, "electronVetoBit/O");
    _tree.Branch("photon.looseElectronVetoBit", &looseElectronVetoBit, "looseElectronVetoBit/O");
    _tree.Branch("photon.isLoose", &isLoose, "isLoose/O");
    _tree.Branch("photon.isMedium", &isMedium, "isMedium/O");
    _tree.Branch("photon.isTight", &isTight, "isTight/O");
    _tree.Branch("photon.isLoosePix", &isLoosePix, "isLoosePix/O");
    _tree.Branch("photon.isMediumPix", &isMediumPix, "isMediumPix/O");
    _tree.Branch("photon.isTightPix", &isTightPix, "isTightPix/O");
    _tree.Branch("photon.isLooseLV", &isLooseLV, "isLooseLV/O");
    _tree.Branch("photon.isMediumLV", &isMediumLV, "isMediumLV/O");
    _tree.Branch("photon.isTightLV", &isTightLV, "isTightLV/O");
  }

  void
  PhotonVars::setAddress(TTree& _tree)
  {
    std::vector<TString> notFound;
    if(_tree.GetBranchStatus("photon.pt")) _tree.SetBranchAddress("photon.pt", &pt);
    else if(!_tree.GetBranch("photon.pt")) notFound.push_back("photon.pt");
    if(_tree.GetBranchStatus("photon.eta")) _tree.SetBranchAddress("photon.eta", &eta);
    else if(!_tree.GetBranch("photon.eta")) notFound.push_back("photon.eta");
    if(_tree.GetBranchStatus("photon.phi")) _tree.SetBranchAddress("photon.phi", &phi);
    else if(!_tree.GetBranch("photon.phi")) notFound.push_back("photon.phi");
    if(_tree.GetBranchStatus("photon.px")) _tree.SetBranchAddress("photon.px", &px);
    else if(!_tree.GetBranch("photon.px")) notFound.push_back("photon.px");
    if(_tree.GetBranchStatus("photon.py")) _tree.SetBranchAddress("photon.py", &py);
    else if(!_tree.GetBranch("photon.py")) notFound.push_back("photon.py");
    if(_tree.GetBranchStatus("photon.pz")) _tree.SetBranchAddress("photon.pz", &pz);
    else if(!_tree.GetBranch("photon.pz")) notFound.push_back("photon.pz");
    if(_tree.GetBranchStatus("photon.energy")) _tree.SetBranchAddress("photon.energy", &energy);
    else if(!_tree.GetBranch("photon.energy")) notFound.push_back("photon.energy");
    if(_tree.GetBranchStatus("photon.hOverE")) _tree.SetBranchAddress("photon.hOverE", &hOverE);
    else if(!_tree.GetBranch("photon.hOverE")) notFound.push_back("photon.hOverE");
    if(_tree.GetBranchStatus("photon.sigmaIetaIeta")) _tree.SetBranchAddress("photon.sigmaIetaIeta", &sigmaIetaIeta);
    else if(!_tree.GetBranch("photon.sigmaIetaIeta")) notFound.push_back("photon.sigmaIetaIeta");
    if(_tree.GetBranchStatus("photon.sigmaIphiIphi")) _tree.SetBranchAddress("photon.sigmaIphiIphi", &sigmaIphiIphi);
    else if(!_tree.GetBranch("photon.sigmaIphiIphi")) notFound.push_back("photon.sigmaIphiIphi");
    if(_tree.GetBranchStatus("photon.etaWidth")) _tree.SetBranchAddress("photon.etaWidth", &etaWidth);
    else if(!_tree.GetBranch("photon.etaWidth")) notFound.push_back("photon.etaWidth");
    if(_tree.GetBranchStatus("photon.phiWidth")) _tree.SetBranchAddress("photon.phiWidth", &phiWidth);
    else if(!_tree.GetBranch("photon.phiWidth")) notFound.push_back("photon.phiWidth");
    if(_tree.GetBranchStatus("photon.r9")) _tree.SetBranchAddress("photon.r9", &r9);
    else if(!_tree.GetBranch("photon.r9")) notFound.push_back("photon.r9");
    if(_tree.GetBranchStatus("photon.r5")) _tree.SetBranchAddress("photon.r5", &r5);
    else if(!_tree.GetBranch("photon.r5")) notFound.push_back("photon.r5");
    if(_tree.GetBranchStatus("photon.trackerIso")) _tree.SetBranchAddress("photon.trackerIso", &trackerIso);
    else if(!_tree.GetBranch("photon.trackerIso")) notFound.push_back("photon.trackerIso");
    if(_tree.GetBranchStatus("photon.ecalIso")) _tree.SetBranchAddress("photon.ecalIso", &ecalIso);
    else if(!_tree.GetBranch("photon.ecalIso")) notFound.push_back("photon.ecalIso");
    if(_tree.GetBranchStatus("photon.hcalIso")) _tree.SetBranchAddress("photon.hcalIso", &hcalIso);
    else if(!_tree.GetBranch("photon.hcalIso")) notFound.push_back("photon.hcalIso");
    if(_tree.GetBranchStatus("photon.chargedHadronIso")) _tree.SetBranchAddress("photon.chargedHadronIso", &chargedHadronIso);
    else if(!_tree.GetBranch("photon.chargedHadronIso")) notFound.push_back("photon.chargedHadronIso");
    if(_tree.GetBranchStatus("photon.neutralHadronIso")) _tree.SetBranchAddress("photon.neutralHadronIso", &neutralHadronIso);
    else if(!_tree.GetBranch("photon.neutralHadronIso")) notFound.push_back("photon.neutralHadronIso");
    if(_tree.GetBranchStatus("photon.photonIso")) _tree.SetBranchAddress("photon.photonIso", &photonIso);
    else if(!_tree.GetBranch("photon.photonIso")) notFound.push_back("photon.photonIso");
    if(_tree.GetBranchStatus("photon.caloX")) _tree.SetBranchAddress("photon.caloX", &caloX);
    else if(!_tree.GetBranch("photon.caloX")) notFound.push_back("photon.caloX");
    if(_tree.GetBranchStatus("photon.caloY")) _tree.SetBranchAddress("photon.caloY", &caloY);
    else if(!_tree.GetBranch("photon.caloY")) notFound.push_back("photon.caloY");
    if(_tree.GetBranchStatus("photon.caloZ")) _tree.SetBranchAddress("photon.caloZ", &caloZ);
    else if(!_tree.GetBranch("photon.caloZ")) notFound.push_back("photon.caloZ");
    if(_tree.GetBranchStatus("photon.iSubdet")) _tree.SetBranchAddress("photon.iSubdet", &iSubdet);
    else if(!_tree.GetBranch("photon.iSubdet")) notFound.push_back("photon.iSubdet");
    if(_tree.GetBranchStatus("photon.superClusterIndex")) _tree.SetBranchAddress("photon.superClusterIndex", &superClusterIndex);
    else if(!_tree.GetBranch("photon.superClusterIndex")) notFound.push_back("photon.superClusterIndex");
    if(_tree.GetBranchStatus("photon.nPixelSeeds")) _tree.SetBranchAddress("photon.nPixelSeeds", &nPixelSeeds);
    else if(!_tree.GetBranch("photon.nPixelSeeds")) notFound.push_back("photon.nPixelSeeds");
    if(_tree.GetBranchStatus("photon.nClusters")) _tree.SetBranchAddress("photon.nClusters", &nClusters);
    else if(!_tree.GetBranch("photon.nClusters")) notFound.push_back("photon.nClusters");
    if(_tree.GetBranchStatus("photon.hasMatchedElectron")) _tree.SetBranchAddress("photon.hasMatchedElectron", &hasMatchedElectron);
    else if(!_tree.GetBranch("photon.hasMatchedElectron")) notFound.push_back("photon.hasMatchedElectron");
    if(_tree.GetBranchStatus("photon.electronVetoBit")) _tree.SetBranchAddress("photon.electronVetoBit", &electronVetoBit);
    else if(!_tree.GetBranch("photon.electronVetoBit")) notFound.push_back("photon.electronVetoBit");
    if(_tree.GetBranchStatus("photon.looseElectronVetoBit")) _tree.SetBranchAddress("photon.looseElectronVetoBit", &looseElectronVetoBit);
    else if(!_tree.GetBranch("photon.looseElectronVetoBit")) notFound.push_back("photon.looseElectronVetoBit");
    if(_tree.GetBranchStatus("photon.isLoose")) _tree.SetBranchAddress("photon.isLoose", &isLoose);
    else if(!_tree.GetBranch("photon.isLoose")) notFound.push_back("photon.isLoose");
    if(_tree.GetBranchStatus("photon.isMedium")) _tree.SetBranchAddress("photon.isMedium", &isMedium);
    else if(!_tree.GetBranch("photon.isMedium")) notFound.push_back("photon.isMedium");
    if(_tree.GetBranchStatus("photon.isTight")) _tree.SetBranchAddress("photon.isTight", &isTight);
    else if(!_tree.GetBranch("photon.isTight")) notFound.push_back("photon.isTight");
    if(_tree.GetBranchStatus("photon.isLoosePix")) _tree.SetBranchAddress("photon.isLoosePix", &isLoosePix);
    else if(!_tree.GetBranch("photon.isLoosePix")) notFound.push_back("photon.isLoosePix");
    if(_tree.GetBranchStatus("photon.isMediumPix")) _tree.SetBranchAddress("photon.isMediumPix", &isMediumPix);
    else if(!_tree.GetBranch("photon.isMediumPix")) notFound.push_back("photon.isMediumPix");
    if(_tree.GetBranchStatus("photon.isTightPix")) _tree.SetBranchAddress("photon.isTightPix", &isTightPix);
    else if(!_tree.GetBranch("photon.isTightPix")) notFound.push_back("photon.isTightPix");
    if(_tree.GetBranchStatus("photon.isLooseLV")) _tree.SetBranchAddress("photon.isLooseLV", &isLooseLV);
    else if(!_tree.GetBranch("photon.isLooseLV")) notFound.push_back("photon.isLooseLV");
    if(_tree.GetBranchStatus("photon.isMediumLV")) _tree.SetBranchAddress("photon.isMediumLV", &isMediumLV);
    else if(!_tree.GetBranch("photon.isMediumLV")) notFound.push_back("photon.isMediumLV");
    if(_tree.GetBranchStatus("photon.isTightLV")) _tree.SetBranchAddress("photon.isTightLV", &isTightLV);
    else if(!_tree.GetBranch("photon.isTightLV")) notFound.push_back("photon.isTightLV");

    for(unsigned iN(0); iN != notFound.size(); ++iN)
      std::cerr << "Branch " << notFound[iN] << " not found in input" << std::endl;
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
    isEcalDriven(false),
    isVeto(false),
    isLoose(false),
    isMedium(false),
    isTight(false)
  {
  }

  void
  ElectronVars::setBranches(TTree& _tree)
  {
    _tree.Branch("electron.pt", &pt, "pt/F");
    _tree.Branch("electron.eta", &eta, "eta/F");
    _tree.Branch("electron.phi", &phi, "phi/F");
    _tree.Branch("electron.px", &px, "px/F");
    _tree.Branch("electron.py", &py, "py/F");
    _tree.Branch("electron.pz", &pz, "pz/F");
    _tree.Branch("electron.energy", &energy, "energy/F");
    _tree.Branch("electron.combRelSubdetIso", &combRelSubdetIso, "combRelSubdetIso/F");
    _tree.Branch("electron.combRelIso", &combRelIso, "combRelIso/F");
    _tree.Branch("electron.deltaEta", &deltaEta, "deltaEta/F");
    _tree.Branch("electron.deltaPhi", &deltaPhi, "deltaPhi/F");
    _tree.Branch("electron.sigmaIetaIeta", &sigmaIetaIeta, "sigmaIetaIeta/F");
    _tree.Branch("electron.sigmaIphiIphi", &sigmaIphiIphi, "sigmaIphiIphi/F");
    _tree.Branch("electron.r9", &r9, "r9/F");
    _tree.Branch("electron.r5", &r5, "r5/F");
    _tree.Branch("electron.etaWidth", &etaWidth, "etaWidth/F");
    _tree.Branch("electron.phiWidth", &phiWidth, "phiWidth/F");
    _tree.Branch("electron.hOverE", &hOverE, "hOverE/F");
    _tree.Branch("electron.d0", &d0, "d0/F");
    _tree.Branch("electron.dz", &dz, "dz/F");
    _tree.Branch("electron.epDiff", &epDiff, "epDiff/F");
    _tree.Branch("electron.vtxFitProb", &vtxFitProb, "vtxFitProb/F");
    _tree.Branch("electron.dCot", &dCot, "dCot/F");
    _tree.Branch("electron.dist", &dist, "dist/F");
    _tree.Branch("electron.caloX", &caloX, "caloX/F");
    _tree.Branch("electron.caloY", &caloY, "caloY/F");
    _tree.Branch("electron.caloZ", &caloZ, "caloZ/F");
    _tree.Branch("electron.iSubdet", &iSubdet, "iSubdet/S");
    _tree.Branch("electron.superClusterIndex", &superClusterIndex, "superClusterIndex/S");
    _tree.Branch("electron.nClusters", &nClusters, "nClusters/b");
    _tree.Branch("electron.nPixelHits", &nPixelHits, "nPixelHits/b");
    _tree.Branch("electron.nMissingHits", &nMissingHits, "nMissingHits/b");
    _tree.Branch("electron.passConversionVeto", &passConversionVeto, "passConversionVeto/O");
    _tree.Branch("electron.isEcalDriven", &isEcalDriven, "isEcalDriven/O");
    _tree.Branch("electron.isVeto", &isVeto, "isVeto/O");
    _tree.Branch("electron.isLoose", &isLoose, "isLoose/O");
    _tree.Branch("electron.isMedium", &isMedium, "isMedium/O");
    _tree.Branch("electron.isTight", &isTight, "isTight/O");
  }

  void
  ElectronVars::setAddress(TTree& _tree)
  {
    std::vector<TString> notFound;
    if(_tree.GetBranchStatus("electron.pt")) _tree.SetBranchAddress("electron.pt", &pt);
    else if(!_tree.GetBranch("electron.pt")) notFound.push_back("electron.pt");
    if(_tree.GetBranchStatus("electron.eta")) _tree.SetBranchAddress("electron.eta", &eta);
    else if(!_tree.GetBranch("electron.eta")) notFound.push_back("electron.eta");
    if(_tree.GetBranchStatus("electron.phi")) _tree.SetBranchAddress("electron.phi", &phi);
    else if(!_tree.GetBranch("electron.phi")) notFound.push_back("electron.phi");
    if(_tree.GetBranchStatus("electron.px")) _tree.SetBranchAddress("electron.px", &px);
    else if(!_tree.GetBranch("electron.px")) notFound.push_back("electron.px");
    if(_tree.GetBranchStatus("electron.py")) _tree.SetBranchAddress("electron.py", &py);
    else if(!_tree.GetBranch("electron.py")) notFound.push_back("electron.py");
    if(_tree.GetBranchStatus("electron.pz")) _tree.SetBranchAddress("electron.pz", &pz);
    else if(!_tree.GetBranch("electron.pz")) notFound.push_back("electron.pz");
    if(_tree.GetBranchStatus("electron.energy")) _tree.SetBranchAddress("electron.energy", &energy);
    else if(!_tree.GetBranch("electron.energy")) notFound.push_back("electron.energy");
    if(_tree.GetBranchStatus("electron.combRelSubdetIso")) _tree.SetBranchAddress("electron.combRelSubdetIso", &combRelSubdetIso);
    else if(!_tree.GetBranch("electron.combRelSubdetIso")) notFound.push_back("electron.combRelSubdetIso");
    if(_tree.GetBranchStatus("electron.combRelIso")) _tree.SetBranchAddress("electron.combRelIso", &combRelIso);
    else if(!_tree.GetBranch("electron.combRelIso")) notFound.push_back("electron.combRelIso");
    if(_tree.GetBranchStatus("electron.deltaEta")) _tree.SetBranchAddress("electron.deltaEta", &deltaEta);
    else if(!_tree.GetBranch("electron.deltaEta")) notFound.push_back("electron.deltaEta");
    if(_tree.GetBranchStatus("electron.deltaPhi")) _tree.SetBranchAddress("electron.deltaPhi", &deltaPhi);
    else if(!_tree.GetBranch("electron.deltaPhi")) notFound.push_back("electron.deltaPhi");
    if(_tree.GetBranchStatus("electron.sigmaIetaIeta")) _tree.SetBranchAddress("electron.sigmaIetaIeta", &sigmaIetaIeta);
    else if(!_tree.GetBranch("electron.sigmaIetaIeta")) notFound.push_back("electron.sigmaIetaIeta");
    if(_tree.GetBranchStatus("electron.sigmaIphiIphi")) _tree.SetBranchAddress("electron.sigmaIphiIphi", &sigmaIphiIphi);
    else if(!_tree.GetBranch("electron.sigmaIphiIphi")) notFound.push_back("electron.sigmaIphiIphi");
    if(_tree.GetBranchStatus("electron.r9")) _tree.SetBranchAddress("electron.r9", &r9);
    else if(!_tree.GetBranch("electron.r9")) notFound.push_back("electron.r9");
    if(_tree.GetBranchStatus("electron.r5")) _tree.SetBranchAddress("electron.r5", &r5);
    else if(!_tree.GetBranch("electron.r5")) notFound.push_back("electron.r5");
    if(_tree.GetBranchStatus("electron.etaWidth")) _tree.SetBranchAddress("electron.etaWidth", &etaWidth);
    else if(!_tree.GetBranch("electron.etaWidth")) notFound.push_back("electron.etaWidth");
    if(_tree.GetBranchStatus("electron.phiWidth")) _tree.SetBranchAddress("electron.phiWidth", &phiWidth);
    else if(!_tree.GetBranch("electron.phiWidth")) notFound.push_back("electron.phiWidth");
    if(_tree.GetBranchStatus("electron.hOverE")) _tree.SetBranchAddress("electron.hOverE", &hOverE);
    else if(!_tree.GetBranch("electron.hOverE")) notFound.push_back("electron.hOverE");
    if(_tree.GetBranchStatus("electron.d0")) _tree.SetBranchAddress("electron.d0", &d0);
    else if(!_tree.GetBranch("electron.d0")) notFound.push_back("electron.d0");
    if(_tree.GetBranchStatus("electron.dz")) _tree.SetBranchAddress("electron.dz", &dz);
    else if(!_tree.GetBranch("electron.dz")) notFound.push_back("electron.dz");
    if(_tree.GetBranchStatus("electron.epDiff")) _tree.SetBranchAddress("electron.epDiff", &epDiff);
    else if(!_tree.GetBranch("electron.epDiff")) notFound.push_back("electron.epDiff");
    if(_tree.GetBranchStatus("electron.vtxFitProb")) _tree.SetBranchAddress("electron.vtxFitProb", &vtxFitProb);
    else if(!_tree.GetBranch("electron.vtxFitProb")) notFound.push_back("electron.vtxFitProb");
    if(_tree.GetBranchStatus("electron.dCot")) _tree.SetBranchAddress("electron.dCot", &dCot);
    else if(!_tree.GetBranch("electron.dCot")) notFound.push_back("electron.dCot");
    if(_tree.GetBranchStatus("electron.dist")) _tree.SetBranchAddress("electron.dist", &dist);
    else if(!_tree.GetBranch("electron.dist")) notFound.push_back("electron.dist");
    if(_tree.GetBranchStatus("electron.caloX")) _tree.SetBranchAddress("electron.caloX", &caloX);
    else if(!_tree.GetBranch("electron.caloX")) notFound.push_back("electron.caloX");
    if(_tree.GetBranchStatus("electron.caloY")) _tree.SetBranchAddress("electron.caloY", &caloY);
    else if(!_tree.GetBranch("electron.caloY")) notFound.push_back("electron.caloY");
    if(_tree.GetBranchStatus("electron.caloZ")) _tree.SetBranchAddress("electron.caloZ", &caloZ);
    else if(!_tree.GetBranch("electron.caloZ")) notFound.push_back("electron.caloZ");
    if(_tree.GetBranchStatus("electron.iSubdet")) _tree.SetBranchAddress("electron.iSubdet", &iSubdet);
    else if(!_tree.GetBranch("electron.iSubdet")) notFound.push_back("electron.iSubdet");
    if(_tree.GetBranchStatus("electron.superClusterIndex")) _tree.SetBranchAddress("electron.superClusterIndex", &superClusterIndex);
    else if(!_tree.GetBranch("electron.superClusterIndex")) notFound.push_back("electron.superClusterIndex");
    if(_tree.GetBranchStatus("electron.nClusters")) _tree.SetBranchAddress("electron.nClusters", &nClusters);
    else if(!_tree.GetBranch("electron.nClusters")) notFound.push_back("electron.nClusters");
    if(_tree.GetBranchStatus("electron.nPixelHits")) _tree.SetBranchAddress("electron.nPixelHits", &nPixelHits);
    else if(!_tree.GetBranch("electron.nPixelHits")) notFound.push_back("electron.nPixelHits");
    if(_tree.GetBranchStatus("electron.nMissingHits")) _tree.SetBranchAddress("electron.nMissingHits", &nMissingHits);
    else if(!_tree.GetBranch("electron.nMissingHits")) notFound.push_back("electron.nMissingHits");
    if(_tree.GetBranchStatus("electron.passConversionVeto")) _tree.SetBranchAddress("electron.passConversionVeto", &passConversionVeto);
    else if(!_tree.GetBranch("electron.passConversionVeto")) notFound.push_back("electron.passConversionVeto");
    if(_tree.GetBranchStatus("electron.isEcalDriven")) _tree.SetBranchAddress("electron.isEcalDriven", &isEcalDriven);
    else if(!_tree.GetBranch("electron.isEcalDriven")) notFound.push_back("electron.isEcalDriven");
    if(_tree.GetBranchStatus("electron.isVeto")) _tree.SetBranchAddress("electron.isVeto", &isVeto);
    else if(!_tree.GetBranch("electron.isVeto")) notFound.push_back("electron.isVeto");
    if(_tree.GetBranchStatus("electron.isLoose")) _tree.SetBranchAddress("electron.isLoose", &isLoose);
    else if(!_tree.GetBranch("electron.isLoose")) notFound.push_back("electron.isLoose");
    if(_tree.GetBranchStatus("electron.isMedium")) _tree.SetBranchAddress("electron.isMedium", &isMedium);
    else if(!_tree.GetBranch("electron.isMedium")) notFound.push_back("electron.isMedium");
    if(_tree.GetBranchStatus("electron.isTight")) _tree.SetBranchAddress("electron.isTight", &isTight);
    else if(!_tree.GetBranch("electron.isTight")) notFound.push_back("electron.isTight");

    for(unsigned iN(0); iN != notFound.size(); ++iN)
      std::cerr << "Branch " << notFound[iN] << " not found in input" << std::endl;
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
    dpTpT(0.),
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

  void
  MuonVars::setBranches(TTree& _tree)
  {
    _tree.Branch("muon.pt", &pt, "pt/F");
    _tree.Branch("muon.eta", &eta, "eta/F");
    _tree.Branch("muon.phi", &phi, "phi/F");
    _tree.Branch("muon.px", &px, "px/F");
    _tree.Branch("muon.py", &py, "py/F");
    _tree.Branch("muon.pz", &pz, "pz/F");
    _tree.Branch("muon.energy", &energy, "energy/F");
    _tree.Branch("muon.normChi2", &normChi2, "normChi2/F");
    _tree.Branch("muon.dxy", &dxy, "dxy/F");
    _tree.Branch("muon.dz", &dz, "dz/F");
    _tree.Branch("muon.combRelSubdetIso", &combRelSubdetIso, "combRelSubdetIso/F");
    _tree.Branch("muon.combRelIso", &combRelIso, "combRelIso/F");
    _tree.Branch("muon.iSubdet", &iSubdet, "iSubdet/S");
    _tree.Branch("muon.dpTpT", &dpTpT, "dpTpT/F");
    _tree.Branch("muon.nMatchedStations", &nMatchedStations, "nMatchedStations/b");
    _tree.Branch("muon.nLayersWithMmt", &nLayersWithMmt, "nLayersWithMmt/b");
    _tree.Branch("muon.nValidMuonHits", &nValidMuonHits, "nValidMuonHits/b");
    _tree.Branch("muon.nValidPixelHits", &nValidPixelHits, "nValidPixelHits/b");
    _tree.Branch("muon.isGlobalMuon", &isGlobalMuon, "isGlobalMuon/O");
    _tree.Branch("muon.isPFMuon", &isPFMuon, "isPFMuon/O");
    _tree.Branch("muon.hasInnerTrack", &hasInnerTrack, "hasInnerTrack/O");
    _tree.Branch("muon.hasGlobalTrack", &hasGlobalTrack, "hasGlobalTrack/O");
    _tree.Branch("muon.hasBestTrack", &hasBestTrack, "hasBestTrack/O");
    _tree.Branch("muon.isLoose", &isLoose, "isLoose/O");
    _tree.Branch("muon.isTight", &isTight, "isTight/O");
  }

  void
  MuonVars::setAddress(TTree& _tree)
  {
    std::vector<TString> notFound;
    if(_tree.GetBranchStatus("muon.pt")) _tree.SetBranchAddress("muon.pt", &pt);
    else if(!_tree.GetBranch("muon.pt")) notFound.push_back("muon.pt");
    if(_tree.GetBranchStatus("muon.eta")) _tree.SetBranchAddress("muon.eta", &eta);
    else if(!_tree.GetBranch("muon.eta")) notFound.push_back("muon.eta");
    if(_tree.GetBranchStatus("muon.phi")) _tree.SetBranchAddress("muon.phi", &phi);
    else if(!_tree.GetBranch("muon.phi")) notFound.push_back("muon.phi");
    if(_tree.GetBranchStatus("muon.px")) _tree.SetBranchAddress("muon.px", &px);
    else if(!_tree.GetBranch("muon.px")) notFound.push_back("muon.px");
    if(_tree.GetBranchStatus("muon.py")) _tree.SetBranchAddress("muon.py", &py);
    else if(!_tree.GetBranch("muon.py")) notFound.push_back("muon.py");
    if(_tree.GetBranchStatus("muon.pz")) _tree.SetBranchAddress("muon.pz", &pz);
    else if(!_tree.GetBranch("muon.pz")) notFound.push_back("muon.pz");
    if(_tree.GetBranchStatus("muon.energy")) _tree.SetBranchAddress("muon.energy", &energy);
    else if(!_tree.GetBranch("muon.energy")) notFound.push_back("muon.energy");
    if(_tree.GetBranchStatus("muon.normChi2")) _tree.SetBranchAddress("muon.normChi2", &normChi2);
    else if(!_tree.GetBranch("muon.normChi2")) notFound.push_back("muon.normChi2");
    if(_tree.GetBranchStatus("muon.dxy")) _tree.SetBranchAddress("muon.dxy", &dxy);
    else if(!_tree.GetBranch("muon.dxy")) notFound.push_back("muon.dxy");
    if(_tree.GetBranchStatus("muon.dz")) _tree.SetBranchAddress("muon.dz", &dz);
    else if(!_tree.GetBranch("muon.dz")) notFound.push_back("muon.dz");
    if(_tree.GetBranchStatus("muon.combRelSubdetIso")) _tree.SetBranchAddress("muon.combRelSubdetIso", &combRelSubdetIso);
    else if(!_tree.GetBranch("muon.combRelSubdetIso")) notFound.push_back("muon.combRelSubdetIso");
    if(_tree.GetBranchStatus("muon.combRelIso")) _tree.SetBranchAddress("muon.combRelIso", &combRelIso);
    else if(!_tree.GetBranch("muon.combRelIso")) notFound.push_back("muon.combRelIso");
    if(_tree.GetBranchStatus("muon.iSubdet")) _tree.SetBranchAddress("muon.iSubdet", &iSubdet);
    else if(!_tree.GetBranch("muon.iSubdet")) notFound.push_back("muon.iSubdet");
    if(_tree.GetBranchStatus("muon.dpTpT")) _tree.SetBranchAddress("muon.dpTpT", &dpTpT);
    else if(!_tree.GetBranch("muon.dpTpT")) notFound.push_back("muon.dpTpT");
    if(_tree.GetBranchStatus("muon.nMatchedStations")) _tree.SetBranchAddress("muon.nMatchedStations", &nMatchedStations);
    else if(!_tree.GetBranch("muon.nMatchedStations")) notFound.push_back("muon.nMatchedStations");
    if(_tree.GetBranchStatus("muon.nLayersWithMmt")) _tree.SetBranchAddress("muon.nLayersWithMmt", &nLayersWithMmt);
    else if(!_tree.GetBranch("muon.nLayersWithMmt")) notFound.push_back("muon.nLayersWithMmt");
    if(_tree.GetBranchStatus("muon.nValidMuonHits")) _tree.SetBranchAddress("muon.nValidMuonHits", &nValidMuonHits);
    else if(!_tree.GetBranch("muon.nValidMuonHits")) notFound.push_back("muon.nValidMuonHits");
    if(_tree.GetBranchStatus("muon.nValidPixelHits")) _tree.SetBranchAddress("muon.nValidPixelHits", &nValidPixelHits);
    else if(!_tree.GetBranch("muon.nValidPixelHits")) notFound.push_back("muon.nValidPixelHits");
    if(_tree.GetBranchStatus("muon.isGlobalMuon")) _tree.SetBranchAddress("muon.isGlobalMuon", &isGlobalMuon);
    else if(!_tree.GetBranch("muon.isGlobalMuon")) notFound.push_back("muon.isGlobalMuon");
    if(_tree.GetBranchStatus("muon.isPFMuon")) _tree.SetBranchAddress("muon.isPFMuon", &isPFMuon);
    else if(!_tree.GetBranch("muon.isPFMuon")) notFound.push_back("muon.isPFMuon");
    if(_tree.GetBranchStatus("muon.hasInnerTrack")) _tree.SetBranchAddress("muon.hasInnerTrack", &hasInnerTrack);
    else if(!_tree.GetBranch("muon.hasInnerTrack")) notFound.push_back("muon.hasInnerTrack");
    if(_tree.GetBranchStatus("muon.hasGlobalTrack")) _tree.SetBranchAddress("muon.hasGlobalTrack", &hasGlobalTrack);
    else if(!_tree.GetBranch("muon.hasGlobalTrack")) notFound.push_back("muon.hasGlobalTrack");
    if(_tree.GetBranchStatus("muon.hasBestTrack")) _tree.SetBranchAddress("muon.hasBestTrack", &hasBestTrack);
    else if(!_tree.GetBranch("muon.hasBestTrack")) notFound.push_back("muon.hasBestTrack");
    if(_tree.GetBranchStatus("muon.isLoose")) _tree.SetBranchAddress("muon.isLoose", &isLoose);
    else if(!_tree.GetBranch("muon.isLoose")) notFound.push_back("muon.isLoose");
    if(_tree.GetBranchStatus("muon.isTight")) _tree.SetBranchAddress("muon.isTight", &isTight);
    else if(!_tree.GetBranch("muon.isTight")) notFound.push_back("muon.isTight");

    for(unsigned iN(0); iN != notFound.size(); ++iN)
      std::cerr << "Branch " << notFound[iN] << " not found in input" << std::endl;
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
    jecUncert(0.),
    chFraction(0.),
    nhFraction(0.),
    ceFraction(0.),
    neFraction(0.),
    quarkLikelihood(0.),
    gluonLikelihood(0.),
    iSubdet(0),
    algoFlavor(0),
    physFlavor(0),
    nConstituents(0),
    nCharged(0),
    passPUJetIdLoose(false),
    isLoose(false)
  {
  }

  void
  JetVars::setBranches(TTree& _tree)
  {
    _tree.Branch("jet.pt", &pt, "pt/F");
    _tree.Branch("jet.eta", &eta, "eta/F");
    _tree.Branch("jet.phi", &phi, "phi/F");
    _tree.Branch("jet.px", &px, "px/F");
    _tree.Branch("jet.py", &py, "py/F");
    _tree.Branch("jet.pz", &pz, "pz/F");
    _tree.Branch("jet.energy", &energy, "energy/F");
    _tree.Branch("jet.jecScale", &jecScale, "jecScale/F");
    _tree.Branch("jet.jecUncert", &jecUncert, "jecUncert/F");
    _tree.Branch("jet.chFraction", &chFraction, "chFraction/F");
    _tree.Branch("jet.nhFraction", &nhFraction, "nhFraction/F");
    _tree.Branch("jet.ceFraction", &ceFraction, "ceFraction/F");
    _tree.Branch("jet.neFraction", &neFraction, "neFraction/F");
    _tree.Branch("jet.quarkLikelihood", &quarkLikelihood, "quarkLikelihood/F");
    _tree.Branch("jet.gluonLikelihood", &gluonLikelihood, "gluonLikelihood/F");
    _tree.Branch("jet.iSubdet", &iSubdet, "iSubdet/S");
    _tree.Branch("jet.algoFlavor", &algoFlavor, "algoFlavor/S");
    _tree.Branch("jet.physFlavor", &physFlavor, "physFlavor/S");
    _tree.Branch("jet.nConstituents", &nConstituents, "nConstituents/b");
    _tree.Branch("jet.nCharged", &nCharged, "nCharged/b");
    _tree.Branch("jet.passPUJetIdLoose", &passPUJetIdLoose, "passPUJetIdLoose/O");
    _tree.Branch("jet.isLoose", &isLoose, "isLoose/O");
  }

  void
  JetVars::setAddress(TTree& _tree)
  {
    std::vector<TString> notFound;
    if(_tree.GetBranchStatus("jet.pt")) _tree.SetBranchAddress("jet.pt", &pt);
    else if(!_tree.GetBranch("jet.pt")) notFound.push_back("jet.pt");
    if(_tree.GetBranchStatus("jet.eta")) _tree.SetBranchAddress("jet.eta", &eta);
    else if(!_tree.GetBranch("jet.eta")) notFound.push_back("jet.eta");
    if(_tree.GetBranchStatus("jet.phi")) _tree.SetBranchAddress("jet.phi", &phi);
    else if(!_tree.GetBranch("jet.phi")) notFound.push_back("jet.phi");
    if(_tree.GetBranchStatus("jet.px")) _tree.SetBranchAddress("jet.px", &px);
    else if(!_tree.GetBranch("jet.px")) notFound.push_back("jet.px");
    if(_tree.GetBranchStatus("jet.py")) _tree.SetBranchAddress("jet.py", &py);
    else if(!_tree.GetBranch("jet.py")) notFound.push_back("jet.py");
    if(_tree.GetBranchStatus("jet.pz")) _tree.SetBranchAddress("jet.pz", &pz);
    else if(!_tree.GetBranch("jet.pz")) notFound.push_back("jet.pz");
    if(_tree.GetBranchStatus("jet.energy")) _tree.SetBranchAddress("jet.energy", &energy);
    else if(!_tree.GetBranch("jet.energy")) notFound.push_back("jet.energy");
    if(_tree.GetBranchStatus("jet.jecScale")) _tree.SetBranchAddress("jet.jecScale", &jecScale);
    else if(!_tree.GetBranch("jet.jecScale")) notFound.push_back("jet.jecScale");
    if(_tree.GetBranchStatus("jet.jecUncert")) _tree.SetBranchAddress("jet.jecUncert", &jecUncert);
    else if(!_tree.GetBranch("jet.jecUncert")) notFound.push_back("jet.jecUncert");
    if(_tree.GetBranchStatus("jet.chFraction")) _tree.SetBranchAddress("jet.chFraction", &chFraction);
    else if(!_tree.GetBranch("jet.chFraction")) notFound.push_back("jet.chFraction");
    if(_tree.GetBranchStatus("jet.nhFraction")) _tree.SetBranchAddress("jet.nhFraction", &nhFraction);
    else if(!_tree.GetBranch("jet.nhFraction")) notFound.push_back("jet.nhFraction");
    if(_tree.GetBranchStatus("jet.ceFraction")) _tree.SetBranchAddress("jet.ceFraction", &ceFraction);
    else if(!_tree.GetBranch("jet.ceFraction")) notFound.push_back("jet.ceFraction");
    if(_tree.GetBranchStatus("jet.neFraction")) _tree.SetBranchAddress("jet.neFraction", &neFraction);
    else if(!_tree.GetBranch("jet.neFraction")) notFound.push_back("jet.neFraction");
    if(_tree.GetBranchStatus("jet.quarkLikelihood")) _tree.SetBranchAddress("jet.quarkLikelihood", &quarkLikelihood);
    else if(!_tree.GetBranch("jet.quarkLikelihood")) notFound.push_back("jet.quarkLikelihood");
    if(_tree.GetBranchStatus("jet.gluonLikelihood")) _tree.SetBranchAddress("jet.gluonLikelihood", &gluonLikelihood);
    else if(!_tree.GetBranch("jet.gluonLikelihood")) notFound.push_back("jet.gluonLikelihood");
    if(_tree.GetBranchStatus("jet.iSubdet")) _tree.SetBranchAddress("jet.iSubdet", &iSubdet);
    else if(!_tree.GetBranch("jet.iSubdet")) notFound.push_back("jet.iSubdet");
    if(_tree.GetBranchStatus("jet.algoFlavor")) _tree.SetBranchAddress("jet.algoFlavor", &algoFlavor);
    else if(!_tree.GetBranch("jet.algoFlavor")) notFound.push_back("jet.algoFlavor");
    if(_tree.GetBranchStatus("jet.physFlavor")) _tree.SetBranchAddress("jet.physFlavor", &physFlavor);
    else if(!_tree.GetBranch("jet.physFlavor")) notFound.push_back("jet.physFlavor");
    if(_tree.GetBranchStatus("jet.nConstituents")) _tree.SetBranchAddress("jet.nConstituents", &nConstituents);
    else if(!_tree.GetBranch("jet.nConstituents")) notFound.push_back("jet.nConstituents");
    if(_tree.GetBranchStatus("jet.nCharged")) _tree.SetBranchAddress("jet.nCharged", &nCharged);
    else if(!_tree.GetBranch("jet.nCharged")) notFound.push_back("jet.nCharged");
    if(_tree.GetBranchStatus("jet.passPUJetIdLoose")) _tree.SetBranchAddress("jet.passPUJetIdLoose", &passPUJetIdLoose);
    else if(!_tree.GetBranch("jet.passPUJetIdLoose")) notFound.push_back("jet.passPUJetIdLoose");
    if(_tree.GetBranchStatus("jet.isLoose")) _tree.SetBranchAddress("jet.isLoose", &isLoose);
    else if(!_tree.GetBranch("jet.isLoose")) notFound.push_back("jet.isLoose");

    for(unsigned iN(0); iN != notFound.size(); ++iN)
      std::cerr << "Branch " << notFound[iN] << " not found in input" << std::endl;
  }

  VertexVars::VertexVars() :
    x(0.),
    y(0.),
    z(0.),
    rho(0.),
    sumPt2(0.),
    chi2(0.),
    ndof(0.),
    nTracks(0),
    isGood(false)
  {
  }

  void
  VertexVars::setBranches(TTree& _tree)
  {
    _tree.Branch("vertex.x", &x, "x/F");
    _tree.Branch("vertex.y", &y, "y/F");
    _tree.Branch("vertex.z", &z, "z/F");
    _tree.Branch("vertex.rho", &rho, "rho/F");
    _tree.Branch("vertex.sumPt2", &sumPt2, "sumPt2/F");
    _tree.Branch("vertex.chi2", &chi2, "chi2/F");
    _tree.Branch("vertex.ndof", &ndof, "ndof/F");
    _tree.Branch("vertex.nTracks", &nTracks, "nTracks/s");
    _tree.Branch("vertex.isGood", &isGood, "isGood/O");
  }

  void
  VertexVars::setAddress(TTree& _tree)
  {
    std::vector<TString> notFound;
    if(_tree.GetBranchStatus("vertex.x")) _tree.SetBranchAddress("vertex.x", &x);
    else if(!_tree.GetBranch("vertex.x")) notFound.push_back("vertex.x");
    if(_tree.GetBranchStatus("vertex.y")) _tree.SetBranchAddress("vertex.y", &y);
    else if(!_tree.GetBranch("vertex.y")) notFound.push_back("vertex.y");
    if(_tree.GetBranchStatus("vertex.z")) _tree.SetBranchAddress("vertex.z", &z);
    else if(!_tree.GetBranch("vertex.z")) notFound.push_back("vertex.z");
    if(_tree.GetBranchStatus("vertex.rho")) _tree.SetBranchAddress("vertex.rho", &rho);
    else if(!_tree.GetBranch("vertex.rho")) notFound.push_back("vertex.rho");
    if(_tree.GetBranchStatus("vertex.sumPt2")) _tree.SetBranchAddress("vertex.sumPt2", &sumPt2);
    else if(!_tree.GetBranch("vertex.sumPt2")) notFound.push_back("vertex.sumPt2");
    if(_tree.GetBranchStatus("vertex.chi2")) _tree.SetBranchAddress("vertex.chi2", &chi2);
    else if(!_tree.GetBranch("vertex.chi2")) notFound.push_back("vertex.chi2");
    if(_tree.GetBranchStatus("vertex.ndof")) _tree.SetBranchAddress("vertex.ndof", &ndof);
    else if(!_tree.GetBranch("vertex.ndof")) notFound.push_back("vertex.ndof");
    if(_tree.GetBranchStatus("vertex.nTracks")) _tree.SetBranchAddress("vertex.nTracks", &nTracks);
    else if(!_tree.GetBranch("vertex.nTracks")) notFound.push_back("vertex.nTracks");
    if(_tree.GetBranchStatus("vertex.isGood")) _tree.SetBranchAddress("vertex.isGood", &isGood);
    else if(!_tree.GetBranch("vertex.isGood")) notFound.push_back("vertex.isGood");

    for(unsigned iN(0); iN != notFound.size(); ++iN)
      std::cerr << "Branch " << notFound[iN] << " not found in input" << std::endl;
  }


/* START USER-DEFINED IMPLEMENTATION (DO NOT MODIFY THIS LINE) */

#ifdef STANDALONE
  void
  PhotonVars::set(Photon const&, Event const&)
  {
  }
#else
  void
  photonEffectiveAreas(double _eta, double* _effA)
  {
    double& effACH(_effA[0]);
    double& effANH(_effA[1]);
    double& effAPh(_effA[2]);

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

    trackerIso = _ph.trkSumPtHollowConeDR03 - 0.002 * pt;

    ecalIso = _ph.ecalRecHitSumEtConeDR03 - 0.012 * pt;

    hcalIso = _ph.hcalDepth1TowerSumEtConeDR03 + _ph.hcalDepth2TowerSumEtConeDR03 - 0.005 * pt;

    sigmaIetaIeta = _ph.sigmaIetaIeta;

    sigmaIphiIphi = _ph.sigmaIphiIphi;

    r9 = _ph.r9;

    if(!_ph.superCluster)
      throw Exception(Exception::kObjectAnomaly, "Photon with no cluster");

    SuperCluster const& sc(*_ph.superCluster);

    superClusterIndex = _ph.superClusterIndex;

    caloX = _ph.caloPosition.X();
    caloY = _ph.caloPosition.Y();
    caloZ = _ph.caloPosition.Z();

    r5 = _ph.e1x5 / sc.energy;

    etaWidth = sc.etaWidth;

    phiWidth = sc.phiWidth;

    double effA[3]; // chargedHadron, neutralHadron, photon
    photonEffectiveAreas(absEta, effA);
 
    double rho(_event.rho);

    chargedHadronIso = _ph.chargedHadronIso - rho * effA[0];

    neutralHadronIso = _ph.neutralHadronIso - rho * effA[1] - 0.04 * pt;

    photonIso = _ph.photonIso - rho * effA[2] - 0.005 * pt;

    nPixelSeeds = _ph.nPixelSeeds;

    nClusters = sc.basicClusterIndices.size();

    electronVetoBit = _ph.passelectronveto;

    // Reproducing hasMatchedPromptElectron() implementation in RecoEgamma/EgammaTools/src/ConversionTools.cc
    // but allowing the electron to miss 1 hit
    typename ElectronCollectionMap::const_iterator electronsSrc(_event.electrons.find("gsfElectrons"));
    if(electronsSrc == _event.electrons.end())
      throw Exception(Exception::kFormatError, "GsfElectrons not in event");

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
#endif

  /*static*/
  void
  PhotonVars::setBranchStatus(TTree& _tree)
  {
    _tree.SetBranchStatus("photons_photons*", 1);
    _tree.SetBranchStatus("electrons_gsfElectrons*", 1);
    _tree.SetBranchStatus("superClusters*", 1);
    _tree.SetBranchStatus("rho", 1);
  }

#ifdef STANDALONE
  void
  ElectronVars::set(Electron const&, Event const&)
  {
  }
#else
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
      throw Exception(Exception::kObjectAnomaly, "Electron with no cluster");

    SuperCluster const& sc(*_el.superCluster);

    superClusterIndex = _el.superClusterIndex;

    double absEta(std::abs(sc.position.Eta()));

    if(absEta < etaGapBegin) iSubdet = 0;
    else if(absEta < etaGapEnd) iSubdet = -1;
    else if(absEta < etaMax) iSubdet = 1;
    else iSubdet = -1;

    if(!_el.gsfTrack)
      throw Exception(Exception::kObjectAnomaly, "Electron with no track");

    // TODO: Not ideal to do this for every single Electron and Muon
    unsigned nV(_event.vertices.size());
    unsigned iV(0);
    while(iV != nV && !ObjectSelector::isGoodVertex(_event.vertices[iV])) ++iV;
    if(iV == nV)
      throw Exception(Exception::kEventAnomaly, "Event with no good vertex");

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

    isEcalDriven = _el.ecalDriven();

    passConversionVeto = _el.passConversionVeto;

    isVeto = ObjectSelector::isGoodElectron(*this, ElVeto12);
    isLoose = ObjectSelector::isGoodElectron(*this, ElLoose12);
    isMedium = ObjectSelector::isGoodElectron(*this, ElMedium12);
    isTight = ObjectSelector::isGoodElectron(*this, ElTight12);
  }
#endif

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

#ifdef STANDALONE
  void
  MuonVars::set(Muon const&, Event const&)
  {
  }
#else
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

    if(absEta < 1.2) iSubdet = 0;
    else if(absEta < 2.4) iSubdet = 1;
    else iSubdet = -1;

    isGlobalMuon = _mu.isGlobalMuon();

    isPFMuon = _mu.isPFMuon();

    hasInnerTrack = (_mu.innerTrack != 0);

    hasGlobalTrack = (_mu.globalTrack != 0);

    Track const* bestTrack(pt > 200. ? _mu.highPtBestTrack : _mu.bestTrack);
    hasBestTrack = (bestTrack != 0);

    // TODO: Not ideal to do this for every single Electron and Muon
    unsigned nV(_event.vertices.size());
    unsigned iV(0);
    while(iV != nV && !ObjectSelector::isGoodVertex(_event.vertices[iV])) ++iV;
    if(iV == nV)
      throw Exception(Exception::kEventAnomaly, "Event with no good vertex");

    Vertex const& primVtx(_event.vertices[iV]);

    nMatchedStations = _mu.nMatchedStations;

    nLayersWithMmt = _mu.nPixelLayersWithMeasurement + _mu.nStripLayersWithMeasurement;

    if(hasGlobalTrack) normChi2 = _mu.globalTrack->normChi2();
    else normChi2 = -1.;

    if(hasGlobalTrack) nValidMuonHits = _mu.globalTrack->numberOfValidMuonHits;
    else nValidMuonHits = 0;

    if(hasBestTrack) dxy = std::abs(bestTrack->dxy(primVtx.position));
    else if(hasInnerTrack) dxy = std::abs(_mu.innerTrack->dxy(primVtx.position));
    else dxy = -1.;

    if(hasBestTrack) dz = std::abs(bestTrack->dz(primVtx.position));
    else if(hasInnerTrack == 1) dz = std::abs(_mu.innerTrack->dz(primVtx.position));
    else dz = -1.;

    if(hasInnerTrack) nValidPixelHits = _mu.innerTrack->numberOfValidPixelHits;
    else nValidPixelHits = 0;

    combRelSubdetIso = (_mu.ecalIsoR03 + _mu.hcalIsoR03 + _mu.trackIsoR03) / pt; 

    combRelIso = (_mu.sumChargedHadronPt04 + std::max(0., _mu.sumNeutralHadronEt04 + _mu.sumPhotonEt04 - 0.5 * _mu.sumPUPt04)) / pt;

    if(hasBestTrack) dpTpT = bestTrack->ptError / pt;
    else dpTpT = -1.;

    isLoose = ObjectSelector::isGoodMuon(*this, MuLoose12);
    isTight = ObjectSelector::isGoodMuon(*this, MuTight12);
  }
#endif  

  /*static*/
  void
  MuonVars::setBranchStatus(TTree& _tree)
  {
    _tree.SetBranchStatus("muons_muons*", 1);
    _tree.SetBranchStatus("tracks*", 1);
    _tree.SetBranchStatus("vertices*", 1);
  }

#ifdef STANDALONE
  void
  JetVars::set(PFJet const&, Event const&)
  {
  }
#else  
  void
  JetVars::set(PFJet const& _jet, Event const& _event)
  {
    jecScale = _jet.jecScaleFactors.find("L1FastL2L3")->second;
    jecUncert = _jet.jecUncertainty;

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

    double rawE(_jet.momentum.E());

    chFraction = _jet.chargedHadronEnergy / rawE;

    nhFraction = _jet.neutralHadronEnergy / rawE;

    ceFraction = _jet.chargedEmEnergy / rawE;

    neFraction = _jet.neutralEmEnergy / rawE;

    nConstituents = _jet.nConstituents;

    nCharged = _jet.chargedMultiplicity;

    algoFlavor = _jet.algDefFlavour;

    physFlavor = _jet.phyDefFlavour;

    passPUJetIdLoose = _jet.passPuJetIdLoose(kPUJetIdFull);

    quarkLikelihood = _jet.qgDiscriminators[kQuarkLikelihood];

    gluonLikelihood = _jet.qgDiscriminators[kGluonMLP];

    isLoose = ObjectSelector::isGoodJet(*this, JtLoose);
  }
#endif

  /*static*/
  void
  JetVars::setBranchStatus(TTree& _tree)
  {
    _tree.SetBranchStatus("pfJets_ak5*", 1);
    _tree.SetBranchStatus("pfJets_ak5chs*", 0);
  }

#ifdef STANDALONE
  void
  VertexVars::set(Vertex const&)
  {
  }
#else
  void
  VertexVars::set(Vertex const& _vtx)
  {
    x = _vtx.position.X();
    y = _vtx.position.Y();
    z = _vtx.position.Z();
    rho = _vtx.position.Perp();

    nTracks = _vtx.tracksSize;
    sumPt2 = _vtx.sumPt2;
    chi2 = _vtx.chi2;
    ndof = _vtx.ndof;

    isGood = ObjectSelector::isGoodVertex(*this);
  }
#endif
  
  /*static*/
  void
  VertexVars::setBranchStatus(TTree& _tree)
  {
    _tree.SetBranchStatus("vertices*", 1);
  }
/* END USER-DEFINED IMPLEMENTATION (DO NOT MODIFY THIS LINE) */

}
