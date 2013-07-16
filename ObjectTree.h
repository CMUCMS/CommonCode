/* Auto-generated header file */
#ifndef ObjectTree_h
#define ObjectTree_h

#include "ObjectVars.h"

#include "TTree.h"
#include "TString.h"

namespace susy {

  unsigned const NMAX(512);

  class PhotonVarsArray {
  public:
    PhotonVarsArray() {}
    ~PhotonVarsArray() {}
    void setBranches(TTree&);
    void setAddress(TTree&);
    void push_back(PhotonVars const&);
    void clear() { size = 0; }
    PhotonVars at(unsigned) const;

    unsigned size;

    float pt[NMAX];
    float eta[NMAX];
    float phi[NMAX];
    float px[NMAX];
    float py[NMAX];
    float pz[NMAX];
    float energy[NMAX];
    float hOverE[NMAX];
    float sigmaIetaIeta[NMAX];
    float sigmaIphiIphi[NMAX];
    float etaWidth[NMAX];
    float phiWidth[NMAX];
    float r9[NMAX];
    float r5[NMAX];
    float trackerIso[NMAX];
    float ecalIso[NMAX];
    float hcalIso[NMAX];
    float chargedHadronIso[NMAX];
    float neutralHadronIso[NMAX];
    float photonIso[NMAX];
    float caloX[NMAX];
    float caloY[NMAX];
    float caloZ[NMAX];
    short iSubdet[NMAX];
    short superClusterIndex[NMAX];
    unsigned char nPixelSeeds[NMAX];
    unsigned char nClusters[NMAX];
    bool hasMatchedElectron[NMAX];
    bool electronVetoBit[NMAX];
    bool looseElectronVetoBit[NMAX];
    bool isLoose[NMAX];
    bool isMedium[NMAX];
    bool isTight[NMAX];
    bool isLoosePix[NMAX];
    bool isMediumPix[NMAX];
    bool isTightPix[NMAX];
    bool isLooseLV[NMAX];
    bool isMediumLV[NMAX];
    bool isTightLV[NMAX];
  };

  class ElectronVarsArray {
  public:
    ElectronVarsArray() {}
    ~ElectronVarsArray() {}
    void setBranches(TTree&);
    void setAddress(TTree&);
    void push_back(ElectronVars const&);
    void clear() { size = 0; }
    ElectronVars at(unsigned) const;

    unsigned size;

    float pt[NMAX];
    float eta[NMAX];
    float phi[NMAX];
    float px[NMAX];
    float py[NMAX];
    float pz[NMAX];
    float energy[NMAX];
    float combRelSubdetIso[NMAX];
    float combRelIso[NMAX];
    float deltaEta[NMAX];
    float deltaPhi[NMAX];
    float sigmaIetaIeta[NMAX];
    float sigmaIphiIphi[NMAX];
    float r9[NMAX];
    float r5[NMAX];
    float etaWidth[NMAX];
    float phiWidth[NMAX];
    float hOverE[NMAX];
    float d0[NMAX];
    float dz[NMAX];
    float epDiff[NMAX];
    float vtxFitProb[NMAX];
    float dCot[NMAX];
    float dist[NMAX];
    float caloX[NMAX];
    float caloY[NMAX];
    float caloZ[NMAX];
    short iSubdet[NMAX];
    short superClusterIndex[NMAX];
    unsigned char nClusters[NMAX];
    unsigned char nPixelHits[NMAX];
    unsigned char nMissingHits[NMAX];
    bool passConversionVeto[NMAX];
    bool isVeto[NMAX];
    bool isLoose[NMAX];
    bool isMedium[NMAX];
    bool isTight[NMAX];
  };

  class MuonVarsArray {
  public:
    MuonVarsArray() {}
    ~MuonVarsArray() {}
    void setBranches(TTree&);
    void setAddress(TTree&);
    void push_back(MuonVars const&);
    void clear() { size = 0; }
    MuonVars at(unsigned) const;

    unsigned size;

    float pt[NMAX];
    float eta[NMAX];
    float phi[NMAX];
    float px[NMAX];
    float py[NMAX];
    float pz[NMAX];
    float energy[NMAX];
    float normChi2[NMAX];
    float dxy[NMAX];
    float dz[NMAX];
    float combRelSubdetIso[NMAX];
    float combRelIso[NMAX];
    short iSubdet[NMAX];
    unsigned char nMatchedStations[NMAX];
    unsigned char nLayersWithMmt[NMAX];
    unsigned char nValidMuonHits[NMAX];
    unsigned char nValidPixelHits[NMAX];
    bool isGlobalMuon[NMAX];
    bool isPFMuon[NMAX];
    bool hasInnerTrack[NMAX];
    bool hasGlobalTrack[NMAX];
    bool hasBestTrack[NMAX];
    bool isLoose[NMAX];
    bool isTight[NMAX];
  };

  class JetVarsArray {
  public:
    JetVarsArray() {}
    ~JetVarsArray() {}
    void setBranches(TTree&);
    void setAddress(TTree&);
    void push_back(JetVars const&);
    void clear() { size = 0; }
    JetVars at(unsigned) const;

    unsigned size;

    float pt[NMAX];
    float eta[NMAX];
    float phi[NMAX];
    float px[NMAX];
    float py[NMAX];
    float pz[NMAX];
    float energy[NMAX];
    float jecScale[NMAX];
    float chFraction[NMAX];
    float nhFraction[NMAX];
    float ceFraction[NMAX];
    float neFraction[NMAX];
    short iSubdet[NMAX];
    unsigned char nConstituents[NMAX];
    unsigned char nCharged[NMAX];
    bool isLoose[NMAX];
  };

  class VertexVarsArray {
  public:
    VertexVarsArray() {}
    ~VertexVarsArray() {}
    void setBranches(TTree&);
    void setAddress(TTree&);
    void push_back(VertexVars const&);
    void clear() { size = 0; }
    VertexVars at(unsigned) const;

    unsigned size;

    float x[NMAX];
    float y[NMAX];
    float z[NMAX];
    float rho[NMAX];
    float sumPt2[NMAX];
    float chi2[NMAX];
    float ndof[NMAX];
    bool isGood[NMAX];
  };

  class ObjectTree {
  public:    
    ObjectTree();
    ~ObjectTree();

    void setOutput(TString const&, bool = true, bool = true, bool = true, bool = true, bool = true);
    void setOutput(TTree&, bool = true, bool = true, bool = true, bool = true, bool = true);
    static void setBranchStatus(TTree&, bool = true, bool = true, bool = true, bool = true, bool = true);
    void initEvent(Event const&);
    void fill() { output_->Fill(); }
    void save(PhotonVars const& _vars) { photonArray_.push_back(_vars); }
    void save(ElectronVars const& _vars) { electronArray_.push_back(_vars); }
    void save(MuonVars const& _vars) { muonArray_.push_back(_vars); }
    void save(JetVars const& _vars) { jetArray_.push_back(_vars); }
    void save(VertexVars const& _vars) { vertexArray_.push_back(_vars); }
    unsigned getPhotonSize() const { return photonArray_.size; }
    unsigned getElectronSize() const { return electronArray_.size; }
    unsigned getMuonSize() const { return muonArray_.size; }
    unsigned getJetSize() const { return jetArray_.size; }
    unsigned getVertexSize() const { return vertexArray_.size; }
    PhotonVarsArray const& getPhotonArray() const { return photonArray_; }
    ElectronVarsArray const& getElectronArray() const { return electronArray_; }
    MuonVarsArray const& getMuonArray() const { return muonArray_; }
    JetVarsArray const& getJetArray() const { return jetArray_; }
    VertexVarsArray const& getVertexArray() const { return vertexArray_; }
  private:
    void setBranches_(bool, bool, bool, bool, bool);
    PhotonVarsArray photonArray_;
    ElectronVarsArray electronArray_;
    MuonVarsArray muonArray_;
    JetVarsArray jetArray_;
    VertexVarsArray vertexArray_;
    unsigned runNumber_;
    unsigned lumiNumber_;
    unsigned eventNumber_;

    TTree* output_;
    bool ownOutput_;
  };

}

#endif
