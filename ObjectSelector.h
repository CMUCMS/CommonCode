#ifndef ObjectSelectors_h
#define ObjectSelectors_h

#include <bitset>

#include "SusyEvent.h"
#include "ObjectVars.h"

namespace susy {
  enum PhotonId {
    PhLoose12,
    PhMedium12,
    PhTight12,
    PhLoose12Pix,
    PhMedium12Pix,
    PhTight12Pix,
    PhLoose12LV,
    PhMedium12LV,
    PhTight12LV,
    nPhotonId
  };
  enum ElectronId {
    ElVeto12,
    ElLoose12,
    ElMedium12,
    ElTight12,
    nElectronId
  };
  enum MuonId {
    MuLoose12,
    MuSoft12,
    MuTight12,
    nMuonId
  };
  enum JetId {
    JtLoose,
    nJetId
  };

  enum PhotonCriteria {
    PhFiducial,
    PhElectronVeto,
    PhHOverE,
    PhSigmaIetaIeta,
    PhChargedHadronIso,
    PhNeutralHadronIso,
    PhPhotonIso,
    nPhotonCriteria
  };
  enum ElectronCriteria {
    ElFiducial,
    ElCombIso,
    ElDeltaEta,
    ElDeltaPhi,
    ElSigmaIetaIeta,
    ElHOverE,
    ElD0,
    ElDZ,
    ElEPDiff,
    ElMissingHits,
    ElConvVeto,
    nElectronCriteria
  };
  enum MuonCriteria {
    MuFiducial,
    MuGlobalOrTrackerMuon,
    MuGlobalMuon,
    MuPFMuon,
    MuMatchedStations,
    MuLayersWithMmt,
    MuNormChi2,
    MuValidMuonHits,
    MuDxy,
    MuDz,
    MuValidPixelHits,
    MuCombIso,
    nMuonCriteria
  };
  enum JetCriteria {
    JtFiducial,
    JtCHFraction,
    JtNHFraction,
    JtCEFraction,
    JtNEFraction,
    JtNConstituents,
    JtChargedMultiplicity,
    nJetCriteria
  };

  class ObjectSelector {
  public:
    ObjectSelector();
    ~ObjectSelector();

    static ObjectSelector const* getSelector() { return singleton_; }

    static bool isGoodPhoton(PhotonVars const&, PhotonId, std::bitset<nPhotonCriteria>* = 0);
    static bool isGoodElectron(ElectronVars const&, ElectronId, std::bitset<nElectronCriteria>* = 0);
    static bool isGoodMuon(MuonVars const&, MuonId, std::bitset<nMuonCriteria>* = 0);
    static bool isGoodJet(JetVars const&, JetId, std::bitset<nJetCriteria>* = 0);
    static bool isGoodVertex(VertexVars const&);

    static std::bitset<nPhotonCriteria> phReferences[nPhotonId];
    static std::bitset<nElectronCriteria> elReferences[nElectronId];
    static std::bitset<nMuonCriteria> muReferences[nMuonId];
    static std::bitset<nJetCriteria> jtReferences[nJetId];

  private:
    static ObjectSelector const* singleton_;
  };

}

#endif
