/* Auto-generated source file */
#include "ObjectTree.h"
#include "Utilities.h"

#include "TFile.h"

#include <iostream>
#ifndef STANDALONE
#include "SusyEvent.h"
#endif

namespace susy {

  void
  PhotonVarsArray::setBranches(TTree& _tree)
  {
    _tree.Branch("photon.size", &size, "photon.size/i");
    _tree.Branch("photon.pt", pt, "pt[photon.size]/F");
    _tree.Branch("photon.eta", eta, "eta[photon.size]/F");
    _tree.Branch("photon.phi", phi, "phi[photon.size]/F");
    _tree.Branch("photon.px", px, "px[photon.size]/F");
    _tree.Branch("photon.py", py, "py[photon.size]/F");
    _tree.Branch("photon.pz", pz, "pz[photon.size]/F");
    _tree.Branch("photon.energy", energy, "energy[photon.size]/F");
    _tree.Branch("photon.hOverE", hOverE, "hOverE[photon.size]/F");
    _tree.Branch("photon.sigmaIetaIeta", sigmaIetaIeta, "sigmaIetaIeta[photon.size]/F");
    _tree.Branch("photon.sigmaIphiIphi", sigmaIphiIphi, "sigmaIphiIphi[photon.size]/F");
    _tree.Branch("photon.etaWidth", etaWidth, "etaWidth[photon.size]/F");
    _tree.Branch("photon.phiWidth", phiWidth, "phiWidth[photon.size]/F");
    _tree.Branch("photon.r9", r9, "r9[photon.size]/F");
    _tree.Branch("photon.r5", r5, "r5[photon.size]/F");
    _tree.Branch("photon.trackerIso", trackerIso, "trackerIso[photon.size]/F");
    _tree.Branch("photon.ecalIso", ecalIso, "ecalIso[photon.size]/F");
    _tree.Branch("photon.hcalIso", hcalIso, "hcalIso[photon.size]/F");
    _tree.Branch("photon.chargedHadronIso", chargedHadronIso, "chargedHadronIso[photon.size]/F");
    _tree.Branch("photon.neutralHadronIso", neutralHadronIso, "neutralHadronIso[photon.size]/F");
    _tree.Branch("photon.photonIso", photonIso, "photonIso[photon.size]/F");
    _tree.Branch("photon.caloX", caloX, "caloX[photon.size]/F");
    _tree.Branch("photon.caloY", caloY, "caloY[photon.size]/F");
    _tree.Branch("photon.caloZ", caloZ, "caloZ[photon.size]/F");
    _tree.Branch("photon.iSubdet", iSubdet, "iSubdet[photon.size]/S");
    _tree.Branch("photon.superClusterIndex", superClusterIndex, "superClusterIndex[photon.size]/S");
    _tree.Branch("photon.nPixelSeeds", nPixelSeeds, "nPixelSeeds[photon.size]/b");
    _tree.Branch("photon.nClusters", nClusters, "nClusters[photon.size]/b");
    _tree.Branch("photon.hasMatchedElectron", hasMatchedElectron, "hasMatchedElectron[photon.size]/O");
    _tree.Branch("photon.electronVetoBit", electronVetoBit, "electronVetoBit[photon.size]/O");
    _tree.Branch("photon.looseElectronVetoBit", looseElectronVetoBit, "looseElectronVetoBit[photon.size]/O");
    _tree.Branch("photon.isLoose", isLoose, "isLoose[photon.size]/O");
    _tree.Branch("photon.isMedium", isMedium, "isMedium[photon.size]/O");
    _tree.Branch("photon.isTight", isTight, "isTight[photon.size]/O");
    _tree.Branch("photon.isLoosePix", isLoosePix, "isLoosePix[photon.size]/O");
    _tree.Branch("photon.isMediumPix", isMediumPix, "isMediumPix[photon.size]/O");
    _tree.Branch("photon.isTightPix", isTightPix, "isTightPix[photon.size]/O");
    _tree.Branch("photon.isLooseLV", isLooseLV, "isLooseLV[photon.size]/O");
    _tree.Branch("photon.isMediumLV", isMediumLV, "isMediumLV[photon.size]/O");
    _tree.Branch("photon.isTightLV", isTightLV, "isTightLV[photon.size]/O");
  }

  void
  PhotonVarsArray::setAddress(TTree& _tree)
  {
    std::vector<TString> notFound;
    _tree.SetBranchAddress("photon.size", &size);
    if(_tree.GetBranch("photon.pt")) _tree.SetBranchAddress("photon.pt", pt);
    else notFound.push_back("photon.pt");
    if(_tree.GetBranch("photon.eta")) _tree.SetBranchAddress("photon.eta", eta);
    else notFound.push_back("photon.eta");
    if(_tree.GetBranch("photon.phi")) _tree.SetBranchAddress("photon.phi", phi);
    else notFound.push_back("photon.phi");
    if(_tree.GetBranch("photon.px")) _tree.SetBranchAddress("photon.px", px);
    else notFound.push_back("photon.px");
    if(_tree.GetBranch("photon.py")) _tree.SetBranchAddress("photon.py", py);
    else notFound.push_back("photon.py");
    if(_tree.GetBranch("photon.pz")) _tree.SetBranchAddress("photon.pz", pz);
    else notFound.push_back("photon.pz");
    if(_tree.GetBranch("photon.energy")) _tree.SetBranchAddress("photon.energy", energy);
    else notFound.push_back("photon.energy");
    if(_tree.GetBranch("photon.hOverE")) _tree.SetBranchAddress("photon.hOverE", hOverE);
    else notFound.push_back("photon.hOverE");
    if(_tree.GetBranch("photon.sigmaIetaIeta")) _tree.SetBranchAddress("photon.sigmaIetaIeta", sigmaIetaIeta);
    else notFound.push_back("photon.sigmaIetaIeta");
    if(_tree.GetBranch("photon.sigmaIphiIphi")) _tree.SetBranchAddress("photon.sigmaIphiIphi", sigmaIphiIphi);
    else notFound.push_back("photon.sigmaIphiIphi");
    if(_tree.GetBranch("photon.etaWidth")) _tree.SetBranchAddress("photon.etaWidth", etaWidth);
    else notFound.push_back("photon.etaWidth");
    if(_tree.GetBranch("photon.phiWidth")) _tree.SetBranchAddress("photon.phiWidth", phiWidth);
    else notFound.push_back("photon.phiWidth");
    if(_tree.GetBranch("photon.r9")) _tree.SetBranchAddress("photon.r9", r9);
    else notFound.push_back("photon.r9");
    if(_tree.GetBranch("photon.r5")) _tree.SetBranchAddress("photon.r5", r5);
    else notFound.push_back("photon.r5");
    if(_tree.GetBranch("photon.trackerIso")) _tree.SetBranchAddress("photon.trackerIso", trackerIso);
    else notFound.push_back("photon.trackerIso");
    if(_tree.GetBranch("photon.ecalIso")) _tree.SetBranchAddress("photon.ecalIso", ecalIso);
    else notFound.push_back("photon.ecalIso");
    if(_tree.GetBranch("photon.hcalIso")) _tree.SetBranchAddress("photon.hcalIso", hcalIso);
    else notFound.push_back("photon.hcalIso");
    if(_tree.GetBranch("photon.chargedHadronIso")) _tree.SetBranchAddress("photon.chargedHadronIso", chargedHadronIso);
    else notFound.push_back("photon.chargedHadronIso");
    if(_tree.GetBranch("photon.neutralHadronIso")) _tree.SetBranchAddress("photon.neutralHadronIso", neutralHadronIso);
    else notFound.push_back("photon.neutralHadronIso");
    if(_tree.GetBranch("photon.photonIso")) _tree.SetBranchAddress("photon.photonIso", photonIso);
    else notFound.push_back("photon.photonIso");
    if(_tree.GetBranch("photon.caloX")) _tree.SetBranchAddress("photon.caloX", caloX);
    else notFound.push_back("photon.caloX");
    if(_tree.GetBranch("photon.caloY")) _tree.SetBranchAddress("photon.caloY", caloY);
    else notFound.push_back("photon.caloY");
    if(_tree.GetBranch("photon.caloZ")) _tree.SetBranchAddress("photon.caloZ", caloZ);
    else notFound.push_back("photon.caloZ");
    if(_tree.GetBranch("photon.iSubdet")) _tree.SetBranchAddress("photon.iSubdet", iSubdet);
    else notFound.push_back("photon.iSubdet");
    if(_tree.GetBranch("photon.superClusterIndex")) _tree.SetBranchAddress("photon.superClusterIndex", superClusterIndex);
    else notFound.push_back("photon.superClusterIndex");
    if(_tree.GetBranch("photon.nPixelSeeds")) _tree.SetBranchAddress("photon.nPixelSeeds", nPixelSeeds);
    else notFound.push_back("photon.nPixelSeeds");
    if(_tree.GetBranch("photon.nClusters")) _tree.SetBranchAddress("photon.nClusters", nClusters);
    else notFound.push_back("photon.nClusters");
    if(_tree.GetBranch("photon.hasMatchedElectron")) _tree.SetBranchAddress("photon.hasMatchedElectron", hasMatchedElectron);
    else notFound.push_back("photon.hasMatchedElectron");
    if(_tree.GetBranch("photon.electronVetoBit")) _tree.SetBranchAddress("photon.electronVetoBit", electronVetoBit);
    else notFound.push_back("photon.electronVetoBit");
    if(_tree.GetBranch("photon.looseElectronVetoBit")) _tree.SetBranchAddress("photon.looseElectronVetoBit", looseElectronVetoBit);
    else notFound.push_back("photon.looseElectronVetoBit");
    if(_tree.GetBranch("photon.isLoose")) _tree.SetBranchAddress("photon.isLoose", isLoose);
    else notFound.push_back("photon.isLoose");
    if(_tree.GetBranch("photon.isMedium")) _tree.SetBranchAddress("photon.isMedium", isMedium);
    else notFound.push_back("photon.isMedium");
    if(_tree.GetBranch("photon.isTight")) _tree.SetBranchAddress("photon.isTight", isTight);
    else notFound.push_back("photon.isTight");
    if(_tree.GetBranch("photon.isLoosePix")) _tree.SetBranchAddress("photon.isLoosePix", isLoosePix);
    else notFound.push_back("photon.isLoosePix");
    if(_tree.GetBranch("photon.isMediumPix")) _tree.SetBranchAddress("photon.isMediumPix", isMediumPix);
    else notFound.push_back("photon.isMediumPix");
    if(_tree.GetBranch("photon.isTightPix")) _tree.SetBranchAddress("photon.isTightPix", isTightPix);
    else notFound.push_back("photon.isTightPix");
    if(_tree.GetBranch("photon.isLooseLV")) _tree.SetBranchAddress("photon.isLooseLV", isLooseLV);
    else notFound.push_back("photon.isLooseLV");
    if(_tree.GetBranch("photon.isMediumLV")) _tree.SetBranchAddress("photon.isMediumLV", isMediumLV);
    else notFound.push_back("photon.isMediumLV");
    if(_tree.GetBranch("photon.isTightLV")) _tree.SetBranchAddress("photon.isTightLV", isTightLV);
    else notFound.push_back("photon.isTightLV");
    
    for(unsigned iN(0); iN != notFound.size(); ++iN)
      std::cerr << "Branch " << notFound[iN] << " not found in input" << std::endl;
  }

  void
  PhotonVarsArray::push_back(PhotonVars const& _vars)
  {
    if(size == NMAX - 1)
      throw Exception(Exception::kEventAnomaly, "Too many Photons");

    pt[size] = _vars.pt;
    eta[size] = _vars.eta;
    phi[size] = _vars.phi;
    px[size] = _vars.px;
    py[size] = _vars.py;
    pz[size] = _vars.pz;
    energy[size] = _vars.energy;
    hOverE[size] = _vars.hOverE;
    sigmaIetaIeta[size] = _vars.sigmaIetaIeta;
    sigmaIphiIphi[size] = _vars.sigmaIphiIphi;
    etaWidth[size] = _vars.etaWidth;
    phiWidth[size] = _vars.phiWidth;
    r9[size] = _vars.r9;
    r5[size] = _vars.r5;
    trackerIso[size] = _vars.trackerIso;
    ecalIso[size] = _vars.ecalIso;
    hcalIso[size] = _vars.hcalIso;
    chargedHadronIso[size] = _vars.chargedHadronIso;
    neutralHadronIso[size] = _vars.neutralHadronIso;
    photonIso[size] = _vars.photonIso;
    caloX[size] = _vars.caloX;
    caloY[size] = _vars.caloY;
    caloZ[size] = _vars.caloZ;
    iSubdet[size] = _vars.iSubdet;
    superClusterIndex[size] = _vars.superClusterIndex;
    nPixelSeeds[size] = _vars.nPixelSeeds;
    nClusters[size] = _vars.nClusters;
    hasMatchedElectron[size] = _vars.hasMatchedElectron;
    electronVetoBit[size] = _vars.electronVetoBit;
    looseElectronVetoBit[size] = _vars.looseElectronVetoBit;
    isLoose[size] = _vars.isLoose;
    isMedium[size] = _vars.isMedium;
    isTight[size] = _vars.isTight;
    isLoosePix[size] = _vars.isLoosePix;
    isMediumPix[size] = _vars.isMediumPix;
    isTightPix[size] = _vars.isTightPix;
    isLooseLV[size] = _vars.isLooseLV;
    isMediumLV[size] = _vars.isMediumLV;
    isTightLV[size] = _vars.isTightLV;
    ++size;
  }

  PhotonVars
  PhotonVarsArray::at(unsigned _pos) const
  {
    if(_pos >= size)
      throw std::range_error("PhotonVars out-of-bounds");
      
    PhotonVars vars;

    vars.pt = pt[_pos];
    vars.eta = eta[_pos];
    vars.phi = phi[_pos];
    vars.px = px[_pos];
    vars.py = py[_pos];
    vars.pz = pz[_pos];
    vars.energy = energy[_pos];
    vars.hOverE = hOverE[_pos];
    vars.sigmaIetaIeta = sigmaIetaIeta[_pos];
    vars.sigmaIphiIphi = sigmaIphiIphi[_pos];
    vars.etaWidth = etaWidth[_pos];
    vars.phiWidth = phiWidth[_pos];
    vars.r9 = r9[_pos];
    vars.r5 = r5[_pos];
    vars.trackerIso = trackerIso[_pos];
    vars.ecalIso = ecalIso[_pos];
    vars.hcalIso = hcalIso[_pos];
    vars.chargedHadronIso = chargedHadronIso[_pos];
    vars.neutralHadronIso = neutralHadronIso[_pos];
    vars.photonIso = photonIso[_pos];
    vars.caloX = caloX[_pos];
    vars.caloY = caloY[_pos];
    vars.caloZ = caloZ[_pos];
    vars.iSubdet = iSubdet[_pos];
    vars.superClusterIndex = superClusterIndex[_pos];
    vars.nPixelSeeds = nPixelSeeds[_pos];
    vars.nClusters = nClusters[_pos];
    vars.hasMatchedElectron = hasMatchedElectron[_pos];
    vars.electronVetoBit = electronVetoBit[_pos];
    vars.looseElectronVetoBit = looseElectronVetoBit[_pos];
    vars.isLoose = isLoose[_pos];
    vars.isMedium = isMedium[_pos];
    vars.isTight = isTight[_pos];
    vars.isLoosePix = isLoosePix[_pos];
    vars.isMediumPix = isMediumPix[_pos];
    vars.isTightPix = isTightPix[_pos];
    vars.isLooseLV = isLooseLV[_pos];
    vars.isMediumLV = isMediumLV[_pos];
    vars.isTightLV = isTightLV[_pos];
    return vars;
  }

  void
  ElectronVarsArray::setBranches(TTree& _tree)
  {
    _tree.Branch("electron.size", &size, "electron.size/i");
    _tree.Branch("electron.pt", pt, "pt[electron.size]/F");
    _tree.Branch("electron.eta", eta, "eta[electron.size]/F");
    _tree.Branch("electron.phi", phi, "phi[electron.size]/F");
    _tree.Branch("electron.px", px, "px[electron.size]/F");
    _tree.Branch("electron.py", py, "py[electron.size]/F");
    _tree.Branch("electron.pz", pz, "pz[electron.size]/F");
    _tree.Branch("electron.energy", energy, "energy[electron.size]/F");
    _tree.Branch("electron.combRelSubdetIso", combRelSubdetIso, "combRelSubdetIso[electron.size]/F");
    _tree.Branch("electron.combRelIso", combRelIso, "combRelIso[electron.size]/F");
    _tree.Branch("electron.deltaEta", deltaEta, "deltaEta[electron.size]/F");
    _tree.Branch("electron.deltaPhi", deltaPhi, "deltaPhi[electron.size]/F");
    _tree.Branch("electron.sigmaIetaIeta", sigmaIetaIeta, "sigmaIetaIeta[electron.size]/F");
    _tree.Branch("electron.sigmaIphiIphi", sigmaIphiIphi, "sigmaIphiIphi[electron.size]/F");
    _tree.Branch("electron.r9", r9, "r9[electron.size]/F");
    _tree.Branch("electron.r5", r5, "r5[electron.size]/F");
    _tree.Branch("electron.etaWidth", etaWidth, "etaWidth[electron.size]/F");
    _tree.Branch("electron.phiWidth", phiWidth, "phiWidth[electron.size]/F");
    _tree.Branch("electron.hOverE", hOverE, "hOverE[electron.size]/F");
    _tree.Branch("electron.d0", d0, "d0[electron.size]/F");
    _tree.Branch("electron.dz", dz, "dz[electron.size]/F");
    _tree.Branch("electron.epDiff", epDiff, "epDiff[electron.size]/F");
    _tree.Branch("electron.vtxFitProb", vtxFitProb, "vtxFitProb[electron.size]/F");
    _tree.Branch("electron.dCot", dCot, "dCot[electron.size]/F");
    _tree.Branch("electron.dist", dist, "dist[electron.size]/F");
    _tree.Branch("electron.caloX", caloX, "caloX[electron.size]/F");
    _tree.Branch("electron.caloY", caloY, "caloY[electron.size]/F");
    _tree.Branch("electron.caloZ", caloZ, "caloZ[electron.size]/F");
    _tree.Branch("electron.iSubdet", iSubdet, "iSubdet[electron.size]/S");
    _tree.Branch("electron.superClusterIndex", superClusterIndex, "superClusterIndex[electron.size]/S");
    _tree.Branch("electron.nClusters", nClusters, "nClusters[electron.size]/b");
    _tree.Branch("electron.nPixelHits", nPixelHits, "nPixelHits[electron.size]/b");
    _tree.Branch("electron.nMissingHits", nMissingHits, "nMissingHits[electron.size]/b");
    _tree.Branch("electron.passConversionVeto", passConversionVeto, "passConversionVeto[electron.size]/O");
    _tree.Branch("electron.isVeto", isVeto, "isVeto[electron.size]/O");
    _tree.Branch("electron.isLoose", isLoose, "isLoose[electron.size]/O");
    _tree.Branch("electron.isMedium", isMedium, "isMedium[electron.size]/O");
    _tree.Branch("electron.isTight", isTight, "isTight[electron.size]/O");
  }

  void
  ElectronVarsArray::setAddress(TTree& _tree)
  {
    std::vector<TString> notFound;
    _tree.SetBranchAddress("electron.size", &size);
    if(_tree.GetBranch("electron.pt")) _tree.SetBranchAddress("electron.pt", pt);
    else notFound.push_back("electron.pt");
    if(_tree.GetBranch("electron.eta")) _tree.SetBranchAddress("electron.eta", eta);
    else notFound.push_back("electron.eta");
    if(_tree.GetBranch("electron.phi")) _tree.SetBranchAddress("electron.phi", phi);
    else notFound.push_back("electron.phi");
    if(_tree.GetBranch("electron.px")) _tree.SetBranchAddress("electron.px", px);
    else notFound.push_back("electron.px");
    if(_tree.GetBranch("electron.py")) _tree.SetBranchAddress("electron.py", py);
    else notFound.push_back("electron.py");
    if(_tree.GetBranch("electron.pz")) _tree.SetBranchAddress("electron.pz", pz);
    else notFound.push_back("electron.pz");
    if(_tree.GetBranch("electron.energy")) _tree.SetBranchAddress("electron.energy", energy);
    else notFound.push_back("electron.energy");
    if(_tree.GetBranch("electron.combRelSubdetIso")) _tree.SetBranchAddress("electron.combRelSubdetIso", combRelSubdetIso);
    else notFound.push_back("electron.combRelSubdetIso");
    if(_tree.GetBranch("electron.combRelIso")) _tree.SetBranchAddress("electron.combRelIso", combRelIso);
    else notFound.push_back("electron.combRelIso");
    if(_tree.GetBranch("electron.deltaEta")) _tree.SetBranchAddress("electron.deltaEta", deltaEta);
    else notFound.push_back("electron.deltaEta");
    if(_tree.GetBranch("electron.deltaPhi")) _tree.SetBranchAddress("electron.deltaPhi", deltaPhi);
    else notFound.push_back("electron.deltaPhi");
    if(_tree.GetBranch("electron.sigmaIetaIeta")) _tree.SetBranchAddress("electron.sigmaIetaIeta", sigmaIetaIeta);
    else notFound.push_back("electron.sigmaIetaIeta");
    if(_tree.GetBranch("electron.sigmaIphiIphi")) _tree.SetBranchAddress("electron.sigmaIphiIphi", sigmaIphiIphi);
    else notFound.push_back("electron.sigmaIphiIphi");
    if(_tree.GetBranch("electron.r9")) _tree.SetBranchAddress("electron.r9", r9);
    else notFound.push_back("electron.r9");
    if(_tree.GetBranch("electron.r5")) _tree.SetBranchAddress("electron.r5", r5);
    else notFound.push_back("electron.r5");
    if(_tree.GetBranch("electron.etaWidth")) _tree.SetBranchAddress("electron.etaWidth", etaWidth);
    else notFound.push_back("electron.etaWidth");
    if(_tree.GetBranch("electron.phiWidth")) _tree.SetBranchAddress("electron.phiWidth", phiWidth);
    else notFound.push_back("electron.phiWidth");
    if(_tree.GetBranch("electron.hOverE")) _tree.SetBranchAddress("electron.hOverE", hOverE);
    else notFound.push_back("electron.hOverE");
    if(_tree.GetBranch("electron.d0")) _tree.SetBranchAddress("electron.d0", d0);
    else notFound.push_back("electron.d0");
    if(_tree.GetBranch("electron.dz")) _tree.SetBranchAddress("electron.dz", dz);
    else notFound.push_back("electron.dz");
    if(_tree.GetBranch("electron.epDiff")) _tree.SetBranchAddress("electron.epDiff", epDiff);
    else notFound.push_back("electron.epDiff");
    if(_tree.GetBranch("electron.vtxFitProb")) _tree.SetBranchAddress("electron.vtxFitProb", vtxFitProb);
    else notFound.push_back("electron.vtxFitProb");
    if(_tree.GetBranch("electron.dCot")) _tree.SetBranchAddress("electron.dCot", dCot);
    else notFound.push_back("electron.dCot");
    if(_tree.GetBranch("electron.dist")) _tree.SetBranchAddress("electron.dist", dist);
    else notFound.push_back("electron.dist");
    if(_tree.GetBranch("electron.caloX")) _tree.SetBranchAddress("electron.caloX", caloX);
    else notFound.push_back("electron.caloX");
    if(_tree.GetBranch("electron.caloY")) _tree.SetBranchAddress("electron.caloY", caloY);
    else notFound.push_back("electron.caloY");
    if(_tree.GetBranch("electron.caloZ")) _tree.SetBranchAddress("electron.caloZ", caloZ);
    else notFound.push_back("electron.caloZ");
    if(_tree.GetBranch("electron.iSubdet")) _tree.SetBranchAddress("electron.iSubdet", iSubdet);
    else notFound.push_back("electron.iSubdet");
    if(_tree.GetBranch("electron.superClusterIndex")) _tree.SetBranchAddress("electron.superClusterIndex", superClusterIndex);
    else notFound.push_back("electron.superClusterIndex");
    if(_tree.GetBranch("electron.nClusters")) _tree.SetBranchAddress("electron.nClusters", nClusters);
    else notFound.push_back("electron.nClusters");
    if(_tree.GetBranch("electron.nPixelHits")) _tree.SetBranchAddress("electron.nPixelHits", nPixelHits);
    else notFound.push_back("electron.nPixelHits");
    if(_tree.GetBranch("electron.nMissingHits")) _tree.SetBranchAddress("electron.nMissingHits", nMissingHits);
    else notFound.push_back("electron.nMissingHits");
    if(_tree.GetBranch("electron.passConversionVeto")) _tree.SetBranchAddress("electron.passConversionVeto", passConversionVeto);
    else notFound.push_back("electron.passConversionVeto");
    if(_tree.GetBranch("electron.isVeto")) _tree.SetBranchAddress("electron.isVeto", isVeto);
    else notFound.push_back("electron.isVeto");
    if(_tree.GetBranch("electron.isLoose")) _tree.SetBranchAddress("electron.isLoose", isLoose);
    else notFound.push_back("electron.isLoose");
    if(_tree.GetBranch("electron.isMedium")) _tree.SetBranchAddress("electron.isMedium", isMedium);
    else notFound.push_back("electron.isMedium");
    if(_tree.GetBranch("electron.isTight")) _tree.SetBranchAddress("electron.isTight", isTight);
    else notFound.push_back("electron.isTight");
    
    for(unsigned iN(0); iN != notFound.size(); ++iN)
      std::cerr << "Branch " << notFound[iN] << " not found in input" << std::endl;
  }

  void
  ElectronVarsArray::push_back(ElectronVars const& _vars)
  {
    if(size == NMAX - 1)
      throw Exception(Exception::kEventAnomaly, "Too many Electrons");

    pt[size] = _vars.pt;
    eta[size] = _vars.eta;
    phi[size] = _vars.phi;
    px[size] = _vars.px;
    py[size] = _vars.py;
    pz[size] = _vars.pz;
    energy[size] = _vars.energy;
    combRelSubdetIso[size] = _vars.combRelSubdetIso;
    combRelIso[size] = _vars.combRelIso;
    deltaEta[size] = _vars.deltaEta;
    deltaPhi[size] = _vars.deltaPhi;
    sigmaIetaIeta[size] = _vars.sigmaIetaIeta;
    sigmaIphiIphi[size] = _vars.sigmaIphiIphi;
    r9[size] = _vars.r9;
    r5[size] = _vars.r5;
    etaWidth[size] = _vars.etaWidth;
    phiWidth[size] = _vars.phiWidth;
    hOverE[size] = _vars.hOverE;
    d0[size] = _vars.d0;
    dz[size] = _vars.dz;
    epDiff[size] = _vars.epDiff;
    vtxFitProb[size] = _vars.vtxFitProb;
    dCot[size] = _vars.dCot;
    dist[size] = _vars.dist;
    caloX[size] = _vars.caloX;
    caloY[size] = _vars.caloY;
    caloZ[size] = _vars.caloZ;
    iSubdet[size] = _vars.iSubdet;
    superClusterIndex[size] = _vars.superClusterIndex;
    nClusters[size] = _vars.nClusters;
    nPixelHits[size] = _vars.nPixelHits;
    nMissingHits[size] = _vars.nMissingHits;
    passConversionVeto[size] = _vars.passConversionVeto;
    isVeto[size] = _vars.isVeto;
    isLoose[size] = _vars.isLoose;
    isMedium[size] = _vars.isMedium;
    isTight[size] = _vars.isTight;
    ++size;
  }

  ElectronVars
  ElectronVarsArray::at(unsigned _pos) const
  {
    if(_pos >= size)
      throw std::range_error("ElectronVars out-of-bounds");
      
    ElectronVars vars;

    vars.pt = pt[_pos];
    vars.eta = eta[_pos];
    vars.phi = phi[_pos];
    vars.px = px[_pos];
    vars.py = py[_pos];
    vars.pz = pz[_pos];
    vars.energy = energy[_pos];
    vars.combRelSubdetIso = combRelSubdetIso[_pos];
    vars.combRelIso = combRelIso[_pos];
    vars.deltaEta = deltaEta[_pos];
    vars.deltaPhi = deltaPhi[_pos];
    vars.sigmaIetaIeta = sigmaIetaIeta[_pos];
    vars.sigmaIphiIphi = sigmaIphiIphi[_pos];
    vars.r9 = r9[_pos];
    vars.r5 = r5[_pos];
    vars.etaWidth = etaWidth[_pos];
    vars.phiWidth = phiWidth[_pos];
    vars.hOverE = hOverE[_pos];
    vars.d0 = d0[_pos];
    vars.dz = dz[_pos];
    vars.epDiff = epDiff[_pos];
    vars.vtxFitProb = vtxFitProb[_pos];
    vars.dCot = dCot[_pos];
    vars.dist = dist[_pos];
    vars.caloX = caloX[_pos];
    vars.caloY = caloY[_pos];
    vars.caloZ = caloZ[_pos];
    vars.iSubdet = iSubdet[_pos];
    vars.superClusterIndex = superClusterIndex[_pos];
    vars.nClusters = nClusters[_pos];
    vars.nPixelHits = nPixelHits[_pos];
    vars.nMissingHits = nMissingHits[_pos];
    vars.passConversionVeto = passConversionVeto[_pos];
    vars.isVeto = isVeto[_pos];
    vars.isLoose = isLoose[_pos];
    vars.isMedium = isMedium[_pos];
    vars.isTight = isTight[_pos];
    return vars;
  }

  void
  MuonVarsArray::setBranches(TTree& _tree)
  {
    _tree.Branch("muon.size", &size, "muon.size/i");
    _tree.Branch("muon.pt", pt, "pt[muon.size]/F");
    _tree.Branch("muon.eta", eta, "eta[muon.size]/F");
    _tree.Branch("muon.phi", phi, "phi[muon.size]/F");
    _tree.Branch("muon.px", px, "px[muon.size]/F");
    _tree.Branch("muon.py", py, "py[muon.size]/F");
    _tree.Branch("muon.pz", pz, "pz[muon.size]/F");
    _tree.Branch("muon.energy", energy, "energy[muon.size]/F");
    _tree.Branch("muon.normChi2", normChi2, "normChi2[muon.size]/F");
    _tree.Branch("muon.dxy", dxy, "dxy[muon.size]/F");
    _tree.Branch("muon.dz", dz, "dz[muon.size]/F");
    _tree.Branch("muon.combRelSubdetIso", combRelSubdetIso, "combRelSubdetIso[muon.size]/F");
    _tree.Branch("muon.combRelIso", combRelIso, "combRelIso[muon.size]/F");
    _tree.Branch("muon.iSubdet", iSubdet, "iSubdet[muon.size]/S");
    _tree.Branch("muon.dpTpT", dpTpT, "dpTpT[muon.size]/F");
    _tree.Branch("muon.nMatchedStations", nMatchedStations, "nMatchedStations[muon.size]/b");
    _tree.Branch("muon.nLayersWithMmt", nLayersWithMmt, "nLayersWithMmt[muon.size]/b");
    _tree.Branch("muon.nValidMuonHits", nValidMuonHits, "nValidMuonHits[muon.size]/b");
    _tree.Branch("muon.nValidPixelHits", nValidPixelHits, "nValidPixelHits[muon.size]/b");
    _tree.Branch("muon.isGlobalMuon", isGlobalMuon, "isGlobalMuon[muon.size]/O");
    _tree.Branch("muon.isPFMuon", isPFMuon, "isPFMuon[muon.size]/O");
    _tree.Branch("muon.hasInnerTrack", hasInnerTrack, "hasInnerTrack[muon.size]/O");
    _tree.Branch("muon.hasGlobalTrack", hasGlobalTrack, "hasGlobalTrack[muon.size]/O");
    _tree.Branch("muon.hasBestTrack", hasBestTrack, "hasBestTrack[muon.size]/O");
    _tree.Branch("muon.isLoose", isLoose, "isLoose[muon.size]/O");
    _tree.Branch("muon.isTight", isTight, "isTight[muon.size]/O");
  }

  void
  MuonVarsArray::setAddress(TTree& _tree)
  {
    std::vector<TString> notFound;
    _tree.SetBranchAddress("muon.size", &size);
    if(_tree.GetBranch("muon.pt")) _tree.SetBranchAddress("muon.pt", pt);
    else notFound.push_back("muon.pt");
    if(_tree.GetBranch("muon.eta")) _tree.SetBranchAddress("muon.eta", eta);
    else notFound.push_back("muon.eta");
    if(_tree.GetBranch("muon.phi")) _tree.SetBranchAddress("muon.phi", phi);
    else notFound.push_back("muon.phi");
    if(_tree.GetBranch("muon.px")) _tree.SetBranchAddress("muon.px", px);
    else notFound.push_back("muon.px");
    if(_tree.GetBranch("muon.py")) _tree.SetBranchAddress("muon.py", py);
    else notFound.push_back("muon.py");
    if(_tree.GetBranch("muon.pz")) _tree.SetBranchAddress("muon.pz", pz);
    else notFound.push_back("muon.pz");
    if(_tree.GetBranch("muon.energy")) _tree.SetBranchAddress("muon.energy", energy);
    else notFound.push_back("muon.energy");
    if(_tree.GetBranch("muon.normChi2")) _tree.SetBranchAddress("muon.normChi2", normChi2);
    else notFound.push_back("muon.normChi2");
    if(_tree.GetBranch("muon.dxy")) _tree.SetBranchAddress("muon.dxy", dxy);
    else notFound.push_back("muon.dxy");
    if(_tree.GetBranch("muon.dz")) _tree.SetBranchAddress("muon.dz", dz);
    else notFound.push_back("muon.dz");
    if(_tree.GetBranch("muon.combRelSubdetIso")) _tree.SetBranchAddress("muon.combRelSubdetIso", combRelSubdetIso);
    else notFound.push_back("muon.combRelSubdetIso");
    if(_tree.GetBranch("muon.combRelIso")) _tree.SetBranchAddress("muon.combRelIso", combRelIso);
    else notFound.push_back("muon.combRelIso");
    if(_tree.GetBranch("muon.iSubdet")) _tree.SetBranchAddress("muon.iSubdet", iSubdet);
    else notFound.push_back("muon.iSubdet");
    if(_tree.GetBranch("muon.dpTpT")) _tree.SetBranchAddress("muon.dpTpT", dpTpT);
    else notFound.push_back("muon.dpTpT");
    if(_tree.GetBranch("muon.nMatchedStations")) _tree.SetBranchAddress("muon.nMatchedStations", nMatchedStations);
    else notFound.push_back("muon.nMatchedStations");
    if(_tree.GetBranch("muon.nLayersWithMmt")) _tree.SetBranchAddress("muon.nLayersWithMmt", nLayersWithMmt);
    else notFound.push_back("muon.nLayersWithMmt");
    if(_tree.GetBranch("muon.nValidMuonHits")) _tree.SetBranchAddress("muon.nValidMuonHits", nValidMuonHits);
    else notFound.push_back("muon.nValidMuonHits");
    if(_tree.GetBranch("muon.nValidPixelHits")) _tree.SetBranchAddress("muon.nValidPixelHits", nValidPixelHits);
    else notFound.push_back("muon.nValidPixelHits");
    if(_tree.GetBranch("muon.isGlobalMuon")) _tree.SetBranchAddress("muon.isGlobalMuon", isGlobalMuon);
    else notFound.push_back("muon.isGlobalMuon");
    if(_tree.GetBranch("muon.isPFMuon")) _tree.SetBranchAddress("muon.isPFMuon", isPFMuon);
    else notFound.push_back("muon.isPFMuon");
    if(_tree.GetBranch("muon.hasInnerTrack")) _tree.SetBranchAddress("muon.hasInnerTrack", hasInnerTrack);
    else notFound.push_back("muon.hasInnerTrack");
    if(_tree.GetBranch("muon.hasGlobalTrack")) _tree.SetBranchAddress("muon.hasGlobalTrack", hasGlobalTrack);
    else notFound.push_back("muon.hasGlobalTrack");
    if(_tree.GetBranch("muon.hasBestTrack")) _tree.SetBranchAddress("muon.hasBestTrack", hasBestTrack);
    else notFound.push_back("muon.hasBestTrack");
    if(_tree.GetBranch("muon.isLoose")) _tree.SetBranchAddress("muon.isLoose", isLoose);
    else notFound.push_back("muon.isLoose");
    if(_tree.GetBranch("muon.isTight")) _tree.SetBranchAddress("muon.isTight", isTight);
    else notFound.push_back("muon.isTight");
    
    for(unsigned iN(0); iN != notFound.size(); ++iN)
      std::cerr << "Branch " << notFound[iN] << " not found in input" << std::endl;
  }

  void
  MuonVarsArray::push_back(MuonVars const& _vars)
  {
    if(size == NMAX - 1)
      throw Exception(Exception::kEventAnomaly, "Too many Muons");

    pt[size] = _vars.pt;
    eta[size] = _vars.eta;
    phi[size] = _vars.phi;
    px[size] = _vars.px;
    py[size] = _vars.py;
    pz[size] = _vars.pz;
    energy[size] = _vars.energy;
    normChi2[size] = _vars.normChi2;
    dxy[size] = _vars.dxy;
    dz[size] = _vars.dz;
    combRelSubdetIso[size] = _vars.combRelSubdetIso;
    combRelIso[size] = _vars.combRelIso;
    iSubdet[size] = _vars.iSubdet;
    dpTpT[size] = _vars.dpTpT;
    nMatchedStations[size] = _vars.nMatchedStations;
    nLayersWithMmt[size] = _vars.nLayersWithMmt;
    nValidMuonHits[size] = _vars.nValidMuonHits;
    nValidPixelHits[size] = _vars.nValidPixelHits;
    isGlobalMuon[size] = _vars.isGlobalMuon;
    isPFMuon[size] = _vars.isPFMuon;
    hasInnerTrack[size] = _vars.hasInnerTrack;
    hasGlobalTrack[size] = _vars.hasGlobalTrack;
    hasBestTrack[size] = _vars.hasBestTrack;
    isLoose[size] = _vars.isLoose;
    isTight[size] = _vars.isTight;
    ++size;
  }

  MuonVars
  MuonVarsArray::at(unsigned _pos) const
  {
    if(_pos >= size)
      throw std::range_error("MuonVars out-of-bounds");
      
    MuonVars vars;

    vars.pt = pt[_pos];
    vars.eta = eta[_pos];
    vars.phi = phi[_pos];
    vars.px = px[_pos];
    vars.py = py[_pos];
    vars.pz = pz[_pos];
    vars.energy = energy[_pos];
    vars.normChi2 = normChi2[_pos];
    vars.dxy = dxy[_pos];
    vars.dz = dz[_pos];
    vars.combRelSubdetIso = combRelSubdetIso[_pos];
    vars.combRelIso = combRelIso[_pos];
    vars.iSubdet = iSubdet[_pos];
    vars.dpTpT = dpTpT[_pos];
    vars.nMatchedStations = nMatchedStations[_pos];
    vars.nLayersWithMmt = nLayersWithMmt[_pos];
    vars.nValidMuonHits = nValidMuonHits[_pos];
    vars.nValidPixelHits = nValidPixelHits[_pos];
    vars.isGlobalMuon = isGlobalMuon[_pos];
    vars.isPFMuon = isPFMuon[_pos];
    vars.hasInnerTrack = hasInnerTrack[_pos];
    vars.hasGlobalTrack = hasGlobalTrack[_pos];
    vars.hasBestTrack = hasBestTrack[_pos];
    vars.isLoose = isLoose[_pos];
    vars.isTight = isTight[_pos];
    return vars;
  }

  void
  JetVarsArray::setBranches(TTree& _tree)
  {
    _tree.Branch("jet.size", &size, "jet.size/i");
    _tree.Branch("jet.pt", pt, "pt[jet.size]/F");
    _tree.Branch("jet.eta", eta, "eta[jet.size]/F");
    _tree.Branch("jet.phi", phi, "phi[jet.size]/F");
    _tree.Branch("jet.px", px, "px[jet.size]/F");
    _tree.Branch("jet.py", py, "py[jet.size]/F");
    _tree.Branch("jet.pz", pz, "pz[jet.size]/F");
    _tree.Branch("jet.energy", energy, "energy[jet.size]/F");
    _tree.Branch("jet.jecScale", jecScale, "jecScale[jet.size]/F");
    _tree.Branch("jet.chFraction", chFraction, "chFraction[jet.size]/F");
    _tree.Branch("jet.nhFraction", nhFraction, "nhFraction[jet.size]/F");
    _tree.Branch("jet.ceFraction", ceFraction, "ceFraction[jet.size]/F");
    _tree.Branch("jet.neFraction", neFraction, "neFraction[jet.size]/F");
    _tree.Branch("jet.iSubdet", iSubdet, "iSubdet[jet.size]/S");
    _tree.Branch("jet.nConstituents", nConstituents, "nConstituents[jet.size]/b");
    _tree.Branch("jet.nCharged", nCharged, "nCharged[jet.size]/b");
    _tree.Branch("jet.passPUJetIdLoose", passPUJetIdLoose, "passPUJetIdLoose[jet.size]/O");
    _tree.Branch("jet.isLoose", isLoose, "isLoose[jet.size]/O");
  }

  void
  JetVarsArray::setAddress(TTree& _tree)
  {
    std::vector<TString> notFound;
    _tree.SetBranchAddress("jet.size", &size);
    if(_tree.GetBranch("jet.pt")) _tree.SetBranchAddress("jet.pt", pt);
    else notFound.push_back("jet.pt");
    if(_tree.GetBranch("jet.eta")) _tree.SetBranchAddress("jet.eta", eta);
    else notFound.push_back("jet.eta");
    if(_tree.GetBranch("jet.phi")) _tree.SetBranchAddress("jet.phi", phi);
    else notFound.push_back("jet.phi");
    if(_tree.GetBranch("jet.px")) _tree.SetBranchAddress("jet.px", px);
    else notFound.push_back("jet.px");
    if(_tree.GetBranch("jet.py")) _tree.SetBranchAddress("jet.py", py);
    else notFound.push_back("jet.py");
    if(_tree.GetBranch("jet.pz")) _tree.SetBranchAddress("jet.pz", pz);
    else notFound.push_back("jet.pz");
    if(_tree.GetBranch("jet.energy")) _tree.SetBranchAddress("jet.energy", energy);
    else notFound.push_back("jet.energy");
    if(_tree.GetBranch("jet.jecScale")) _tree.SetBranchAddress("jet.jecScale", jecScale);
    else notFound.push_back("jet.jecScale");
    if(_tree.GetBranch("jet.chFraction")) _tree.SetBranchAddress("jet.chFraction", chFraction);
    else notFound.push_back("jet.chFraction");
    if(_tree.GetBranch("jet.nhFraction")) _tree.SetBranchAddress("jet.nhFraction", nhFraction);
    else notFound.push_back("jet.nhFraction");
    if(_tree.GetBranch("jet.ceFraction")) _tree.SetBranchAddress("jet.ceFraction", ceFraction);
    else notFound.push_back("jet.ceFraction");
    if(_tree.GetBranch("jet.neFraction")) _tree.SetBranchAddress("jet.neFraction", neFraction);
    else notFound.push_back("jet.neFraction");
    if(_tree.GetBranch("jet.iSubdet")) _tree.SetBranchAddress("jet.iSubdet", iSubdet);
    else notFound.push_back("jet.iSubdet");
    if(_tree.GetBranch("jet.nConstituents")) _tree.SetBranchAddress("jet.nConstituents", nConstituents);
    else notFound.push_back("jet.nConstituents");
    if(_tree.GetBranch("jet.nCharged")) _tree.SetBranchAddress("jet.nCharged", nCharged);
    else notFound.push_back("jet.nCharged");
    if(_tree.GetBranch("jet.passPUJetIdLoose")) _tree.SetBranchAddress("jet.passPUJetIdLoose", passPUJetIdLoose);
    else notFound.push_back("jet.passPUJetIdLoose");
    if(_tree.GetBranch("jet.isLoose")) _tree.SetBranchAddress("jet.isLoose", isLoose);
    else notFound.push_back("jet.isLoose");
    
    for(unsigned iN(0); iN != notFound.size(); ++iN)
      std::cerr << "Branch " << notFound[iN] << " not found in input" << std::endl;
  }

  void
  JetVarsArray::push_back(JetVars const& _vars)
  {
    if(size == NMAX - 1)
      throw Exception(Exception::kEventAnomaly, "Too many Jets");

    pt[size] = _vars.pt;
    eta[size] = _vars.eta;
    phi[size] = _vars.phi;
    px[size] = _vars.px;
    py[size] = _vars.py;
    pz[size] = _vars.pz;
    energy[size] = _vars.energy;
    jecScale[size] = _vars.jecScale;
    chFraction[size] = _vars.chFraction;
    nhFraction[size] = _vars.nhFraction;
    ceFraction[size] = _vars.ceFraction;
    neFraction[size] = _vars.neFraction;
    iSubdet[size] = _vars.iSubdet;
    nConstituents[size] = _vars.nConstituents;
    nCharged[size] = _vars.nCharged;
    passPUJetIdLoose[size] = _vars.passPUJetIdLoose;
    isLoose[size] = _vars.isLoose;
    ++size;
  }

  JetVars
  JetVarsArray::at(unsigned _pos) const
  {
    if(_pos >= size)
      throw std::range_error("JetVars out-of-bounds");
      
    JetVars vars;

    vars.pt = pt[_pos];
    vars.eta = eta[_pos];
    vars.phi = phi[_pos];
    vars.px = px[_pos];
    vars.py = py[_pos];
    vars.pz = pz[_pos];
    vars.energy = energy[_pos];
    vars.jecScale = jecScale[_pos];
    vars.chFraction = chFraction[_pos];
    vars.nhFraction = nhFraction[_pos];
    vars.ceFraction = ceFraction[_pos];
    vars.neFraction = neFraction[_pos];
    vars.iSubdet = iSubdet[_pos];
    vars.nConstituents = nConstituents[_pos];
    vars.nCharged = nCharged[_pos];
    vars.passPUJetIdLoose = passPUJetIdLoose[_pos];
    vars.isLoose = isLoose[_pos];
    return vars;
  }

  void
  VertexVarsArray::setBranches(TTree& _tree)
  {
    _tree.Branch("vertex.size", &size, "vertex.size/i");
    _tree.Branch("vertex.x", x, "x[vertex.size]/F");
    _tree.Branch("vertex.y", y, "y[vertex.size]/F");
    _tree.Branch("vertex.z", z, "z[vertex.size]/F");
    _tree.Branch("vertex.rho", rho, "rho[vertex.size]/F");
    _tree.Branch("vertex.sumPt2", sumPt2, "sumPt2[vertex.size]/F");
    _tree.Branch("vertex.chi2", chi2, "chi2[vertex.size]/F");
    _tree.Branch("vertex.ndof", ndof, "ndof[vertex.size]/F");
    _tree.Branch("vertex.nTracks", nTracks, "nTracks[vertex.size]/s");
    _tree.Branch("vertex.isGood", isGood, "isGood[vertex.size]/O");
  }

  void
  VertexVarsArray::setAddress(TTree& _tree)
  {
    std::vector<TString> notFound;
    _tree.SetBranchAddress("vertex.size", &size);
    if(_tree.GetBranch("vertex.x")) _tree.SetBranchAddress("vertex.x", x);
    else notFound.push_back("vertex.x");
    if(_tree.GetBranch("vertex.y")) _tree.SetBranchAddress("vertex.y", y);
    else notFound.push_back("vertex.y");
    if(_tree.GetBranch("vertex.z")) _tree.SetBranchAddress("vertex.z", z);
    else notFound.push_back("vertex.z");
    if(_tree.GetBranch("vertex.rho")) _tree.SetBranchAddress("vertex.rho", rho);
    else notFound.push_back("vertex.rho");
    if(_tree.GetBranch("vertex.sumPt2")) _tree.SetBranchAddress("vertex.sumPt2", sumPt2);
    else notFound.push_back("vertex.sumPt2");
    if(_tree.GetBranch("vertex.chi2")) _tree.SetBranchAddress("vertex.chi2", chi2);
    else notFound.push_back("vertex.chi2");
    if(_tree.GetBranch("vertex.ndof")) _tree.SetBranchAddress("vertex.ndof", ndof);
    else notFound.push_back("vertex.ndof");
    if(_tree.GetBranch("vertex.nTracks")) _tree.SetBranchAddress("vertex.nTracks", nTracks);
    else notFound.push_back("vertex.nTracks");
    if(_tree.GetBranch("vertex.isGood")) _tree.SetBranchAddress("vertex.isGood", isGood);
    else notFound.push_back("vertex.isGood");
    
    for(unsigned iN(0); iN != notFound.size(); ++iN)
      std::cerr << "Branch " << notFound[iN] << " not found in input" << std::endl;
  }

  void
  VertexVarsArray::push_back(VertexVars const& _vars)
  {
    if(size == NMAX - 1)
      throw Exception(Exception::kEventAnomaly, "Too many Vertexs");

    x[size] = _vars.x;
    y[size] = _vars.y;
    z[size] = _vars.z;
    rho[size] = _vars.rho;
    sumPt2[size] = _vars.sumPt2;
    chi2[size] = _vars.chi2;
    ndof[size] = _vars.ndof;
    nTracks[size] = _vars.nTracks;
    isGood[size] = _vars.isGood;
    ++size;
  }

  VertexVars
  VertexVarsArray::at(unsigned _pos) const
  {
    if(_pos >= size)
      throw std::range_error("VertexVars out-of-bounds");
      
    VertexVars vars;

    vars.x = x[_pos];
    vars.y = y[_pos];
    vars.z = z[_pos];
    vars.rho = rho[_pos];
    vars.sumPt2 = sumPt2[_pos];
    vars.chi2 = chi2[_pos];
    vars.ndof = ndof[_pos];
    vars.nTracks = nTracks[_pos];
    vars.isGood = isGood[_pos];
    return vars;
  }


  ObjectTree::ObjectTree() :
    photonArray_(),
    electronArray_(),
    muonArray_(),
    jetArray_(),
    vertexArray_(),
    runNumber_(0),
    lumiNumber_(0),
    eventNumber_(0),
    output_(0),
    ownOutput_(false)
  {
  }

  ObjectTree::~ObjectTree()
  {
    if(ownOutput_ && output_){
      TFile* outFile(output_->GetCurrentFile());
      outFile->cd();
      output_->Write();
      delete outFile;
    }
  }

  void
  ObjectTree::setOutput(TString const& _fileName, bool _setPhoton/* = true*/, bool _setElectron/* = true*/, bool _setMuon/* = true*/, bool _setJet/* = true*/, bool _setVertex/* = true*/)
  {
    ownOutput_ = true;

    TFile::Open(_fileName, "recreate");
    output_ = new TTree("objectVars", "Object ID variables");

    setBranches_(_setPhoton, _setElectron, _setMuon, _setJet, _setVertex);
  }

  void
  ObjectTree::setOutput(TTree& _tree, bool _setPhoton/* = true*/, bool _setElectron/* = true*/, bool _setMuon/* = true*/, bool _setJet/* = true*/, bool _setVertex/* = true*/)
  {
    output_ = &_tree;

    setBranches_(_setPhoton, _setElectron, _setMuon, _setJet, _setVertex);
  }

  /*static*/
  void
  ObjectTree::setBranchStatus(TTree& _input, bool _setPhoton/* = true*/, bool _setElectron/* = true*/, bool _setMuon/* = true*/, bool _setJet/* = true*/, bool _setVertex/* = true*/)
  {
    _input.SetBranchStatus("runNumber", 1);
    _input.SetBranchStatus("luminosityBlockNumber", 1);
    _input.SetBranchStatus("eventNumber", 1);

    if(_setPhoton) PhotonVars::setBranchStatus(_input);
    if(_setElectron) ElectronVars::setBranchStatus(_input);
    if(_setMuon) MuonVars::setBranchStatus(_input);
    if(_setJet) JetVars::setBranchStatus(_input);
    if(_setVertex) VertexVars::setBranchStatus(_input);
  }

#ifdef STANDALONE
  void
  ObjectTree::initEvent(Event const&)
  {
    runNumber_ = 0;
    lumiNumber_ = 0;
    eventNumber_ = 0;
#else
  void
  ObjectTree::initEvent(Event const& _event)
  {
    runNumber_ = _event.runNumber;
    lumiNumber_ = _event.luminosityBlockNumber;
    eventNumber_ = _event.eventNumber;
#endif
    photonArray_.clear();
    electronArray_.clear();
    muonArray_.clear();
    jetArray_.clear();
    vertexArray_.clear();
  }

  void
  ObjectTree::setBranches_(bool _setPhoton, bool _setElectron, bool _setMuon, bool _setJet, bool _setVertex)
  {
    output_->Branch("runNumber", &runNumber_, "runNumber/i");
    output_->Branch("lumiNumber", &lumiNumber_, "lumiNumber/i");
    output_->Branch("eventNumber", &eventNumber_, "eventNumber/i");

    if(_setPhoton) photonArray_.setBranches(*output_);
    if(_setElectron) electronArray_.setBranches(*output_);
    if(_setMuon) muonArray_.setBranches(*output_);
    if(_setJet) jetArray_.setBranches(*output_);
    if(_setVertex) vertexArray_.setBranches(*output_);
  }
}
