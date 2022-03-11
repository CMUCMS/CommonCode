#ifndef ObjectVars_h
#define ObjectVars_h

#include "TTree.h"

namespace susy {

  class Event;
  class Photon;
  class Electron;
  class Muon;
  class PFJet;
  class Vertex;

  class PhotonVars {
  public:
    PhotonVars();
    PhotonVars(Photon const& _ph, Event const& _ev) { set(_ph, _ev); }
    void setBranches(TTree&);
    void setAddress(TTree&);
    void set(Photon const&, Event const&);
    static void setBranchStatus(TTree&);

    float pt;
    float eta;
    float phi;
    float px;
    float py;
    float pz;
    float energy;
    float hOverE;
    float sigmaIetaIeta;
    float sigmaIphiIphi;
    float etaWidth;
    float phiWidth;
    float r9;
    float r5;
    float trackerIso;
    float ecalIso;
    float hcalIso;
    float chargedHadronIso;
    float neutralHadronIso;
    float photonIso;
    float caloX;
    float caloY;
    float caloZ;
    short iSubdet;
    short superClusterIndex;
    unsigned char nPixelSeeds;
    unsigned char nClusters;
    bool hasMatchedElectron;
    bool electronVetoBit;
    bool looseElectronVetoBit;
    bool isLoose;
    bool isMedium;
    bool isTight;
    bool isLoosePix;
    bool isMediumPix;
    bool isTightPix;
    bool isLooseLV;
    bool isMediumLV;
    bool isTightLV;
  };

  class ElectronVars {
  public:
    ElectronVars();
    ElectronVars(Electron const& _el, Event const& _ev) { set(_el, _ev); }
    void setBranches(TTree&);
    void setAddress(TTree&);
    void set(Electron const&, Event const&);
    static void setBranchStatus(TTree&);

    float pt;
    float eta;
    float phi;
    float px;
    float py;
    float pz;
    float energy;
    float combRelSubdetIso;
    float combRelIso;
    float deltaEta;
    float deltaPhi;
    float sigmaIetaIeta;
    float sigmaIphiIphi;
    float r9;
    float r5;
    float etaWidth;
    float phiWidth;
    float hOverE;
    float d0;
    float dz;
    float epDiff;
    float vtxFitProb;
    float dCot;
    float dist;
    float caloX;
    float caloY;
    float caloZ;
    short iSubdet;
    short superClusterIndex;
    unsigned char nClusters;
    unsigned char nPixelHits;
    unsigned char nMissingHits;
    bool passConversionVeto;
    bool isEcalDriven;
    bool isVeto;
    bool isLoose;
    bool isMedium;
    bool isTight;
  };

  class MuonVars {
  public:
    MuonVars();
    MuonVars(Muon const& _mu, Event const& _ev) { set(_mu, _ev); }
    void setBranches(TTree&);
    void setAddress(TTree&);
    void set(Muon const&, Event const&);
    static void setBranchStatus(TTree&);

    float pt;
    float eta;
    float phi;
    float px;
    float py;
    float pz;
    float energy;
    float normChi2;
    float dxy;
    float dz;
    float combRelSubdetIso;
    float combRelIso;
    short iSubdet;
    float dpTpT;
    unsigned char nMatchedStations;
    unsigned char nLayersWithMmt;
    unsigned char nValidMuonHits;
    unsigned char nValidPixelHits;
    bool isGlobalMuon;
    bool isPFMuon;
    bool hasInnerTrack;
    bool hasGlobalTrack;
    bool hasBestTrack;
    bool isLoose;
    bool isTight;
  };

  class JetVars {
  public:
    JetVars();
    JetVars(PFJet const& _jt, Event const& _ev) { set(_jt, _ev); }
    void setBranches(TTree&);
    void setAddress(TTree&);
    void set(PFJet const&, Event const&);
    static void setBranchStatus(TTree&);

    float pt;
    float eta;
    float phi;
    float px;
    float py;
    float pz;
    float energy;
    float jecScale;
    float jecUncert;
    float chFraction;
    float nhFraction;
    float ceFraction;
    float neFraction;
    float quarkLikelihood;
    float gluonLikelihood;
    short iSubdet;
    short algoFlavor;
    short physFlavor;
    unsigned char nConstituents;
    unsigned char nCharged;
    bool passPUJetIdLoose;
    bool isLoose;
  };

  class VertexVars {
  public:
    VertexVars();
    VertexVars(Vertex const& _vtx) { set(_vtx); }
    void setBranches(TTree&);
    void setAddress(TTree&);
    void set(Vertex const&);
    static void setBranchStatus(TTree&);

    float x;
    float y;
    float z;
    float rho;
    float sumPt2;
    float chi2;
    float ndof;
    unsigned short nTracks;
    bool isGood;
  };

}

#endif
