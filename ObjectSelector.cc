#include "ObjectSelector.h"

#include "Utilities.h"

#include <cmath>
#include <stdexcept>

// Selector instantiated at the bottom of the file

namespace susy {

  std::bitset<nPhotonCriteria> ObjectSelector::phReferences[nPhotonId];
  std::bitset<nElectronCriteria> ObjectSelector::elReferences[nElectronId];
  std::bitset<nMuonCriteria> ObjectSelector::muReferences[nMuonId];
  std::bitset<nJetCriteria> ObjectSelector::jtReferences[nJetId];
  ObjectSelector const* ObjectSelector::singleton_(0);

  ObjectSelector::ObjectSelector()
  {
    if(getSelector())
      throw std::logic_error("ObjectSelector duplicated");

    singleton_ = this;

    // Reference bit masks - customize whenever a new ID is defined

    for(unsigned iId(0); iId != nPhotonId; ++iId)
      phReferences[iId].set();

    for(unsigned iId(0); iId != nElectronId; ++iId)
      elReferences[iId].set();

    muReferences[MuLoose12].reset();
    unsigned muLooseBits[] = {MuFiducial, MuGlobalOrTrackerMuon, MuPFMuon};
    for(unsigned iB(0); iB != sizeof(muLooseBits) / sizeof(unsigned); ++iB)
      muReferences[MuLoose12][muLooseBits[iB]] = true;
    muReferences[MuSoft12].reset();
    muReferences[MuTight12].reset();
    unsigned muTightBits[] = {MuFiducial, MuGlobalMuon, MuPFMuon, MuMatchedStations, MuLayersWithMmt,
                              MuNormChi2, MuValidMuonHits, MuDxy, MuDz, MuValidPixelHits, MuCombIso, MudpTpT};
    for(unsigned iB(0); iB != sizeof(muTightBits) / sizeof(unsigned); ++iB)
      muReferences[MuTight12][muTightBits[iB]] = true;

    for(unsigned iId(0); iId != nJetId; ++iId)
      jtReferences[iId].set();
  }

  ObjectSelector::~ObjectSelector()
  {
    singleton_ = 0;
  }

  /*static*/
  bool
  ObjectSelector::isGoodPhoton(PhotonVars const& _ph, PhotonId _idType, std::bitset<nPhotonCriteria>* _results/* = 0*/)
  {
    std::bitset<nPhotonCriteria> idResults;

    if(_idType == PhLoose12 || _idType == PhMedium12 || _idType == PhTight12 ||
       _idType == PhLoose12Pix || _idType == PhMedium12Pix || _idType == PhTight12Pix ||
       _idType == PhLoose12LV || _idType == PhMedium12LV || _idType == PhTight12LV){

      double maxSigmaIetaIeta[2][3] =
        {{0.012, 0.011, 0.011},
         {0.034, 0.033, 0.031}};
      double maxChargedHadronIso[2][3] =
        {{2.6, 1.5, 0.7},
         {2.3, 1.2, 0.5}};
      double maxNeutralHadronIso[2][3] =
        {{3.5, 1.0, 0.4},
         {2.9, 1.5, 1.5}};
      double maxPhotonIso[2][3] =
        {{1.3, 0.7, 0.5},
         {std::numeric_limits<double>::max(), 1.0, 1.0}};

      int iId((_idType - PhLoose12) % 3);
      bool usePixelSeed((_idType - PhLoose12) / 3 == 1);
      bool useLooseVeto((_idType - PhLoose12) / 3 == 2);

      idResults[PhFiducial] = _ph.iSubdet != -1;
      if(usePixelSeed) idResults[PhElectronVeto] = _ph.nPixelSeeds == 0;
      else if(useLooseVeto) idResults[PhElectronVeto] = _ph.looseElectronVetoBit;
      else idResults[PhElectronVeto] = _ph.electronVetoBit;
      idResults[PhHOverE] = _ph.hOverE < 0.05;
      if(idResults[PhFiducial]){
        idResults[PhSigmaIetaIeta] = _ph.sigmaIetaIeta < maxSigmaIetaIeta[_ph.iSubdet][iId];
        idResults[PhChargedHadronIso] = _ph.chargedHadronIso < maxChargedHadronIso[_ph.iSubdet][iId];
        idResults[PhNeutralHadronIso] = _ph.neutralHadronIso < maxNeutralHadronIso[_ph.iSubdet][iId];
        idResults[PhPhotonIso] = _ph.photonIso < maxPhotonIso[_ph.iSubdet][iId];
      }

    }

    if(_results) *_results = idResults;

    return idResults == phReferences[_idType];
  }

  /*static*/
  bool
  ObjectSelector::isGoodElectron(ElectronVars const& _el, ElectronId _idType, std::bitset<nElectronCriteria>* _results/* = 0*/)
  {
    std::bitset<nElectronCriteria> idResults;

    if(_idType == ElVeto12 || _idType == ElLoose12 || _idType == ElMedium12 || _idType == ElTight12){

      double maxDeltaEta[2][4] =
        {{0.007, 0.007, 0.004, 0.004},
         {0.01, 0.009, 0.007, 0.005}};
      double maxDeltaPhi[2][4] = 
        {{0.8, 0.15, 0.06, 0.03},
         {0.7, 0.1, 0.03, 0.02}};
      double maxSigmaIetaIeta[2][4] =
        {{0.01, 0.01, 0.01, 0.01},
         {0.03, 0.03, 0.03, 0.03}};
      double maxHOverE[2][4] =
        {{0.15, 0.12, 0.12, 0.12},
         {std::numeric_limits<double>::max(), 0.1, 0.1, 0.1}};
      double maxD0[4] =
         {0.04, 0.02, 0.02, 0.02};
      double maxDZ[4] =
         {0.2, 0.2, 0.1, 0.1};
      double maxEP[4] =
        {std::numeric_limits<double>::max(), 0.05, 0.05, 0.05};
      double maxIso[4] =
        {0.15, 0.15, 0.15, 0.1};
      int maxMissingHits[4] = 
        {std::numeric_limits<int>::max(), 1, 1, 0};

      int iId(_idType - ElVeto12);

      idResults[ElFiducial] = _el.iSubdet != -1;
      idResults[ElCombIso] = _el.combRelIso < maxIso[iId];
      if(idResults[ElFiducial]){
        idResults[ElDeltaEta] = _el.deltaEta < maxDeltaEta[_el.iSubdet][iId];
        idResults[ElDeltaPhi] = _el.deltaPhi < maxDeltaPhi[_el.iSubdet][iId];
        idResults[ElSigmaIetaIeta] = _el.sigmaIetaIeta < maxSigmaIetaIeta[_el.iSubdet][iId];
        idResults[ElHOverE] = _el.hOverE < maxHOverE[_el.iSubdet][iId];
      }
      idResults[ElD0] = _el.d0 < maxD0[iId];
      idResults[ElDZ] = _el.dz < maxDZ[iId];
      idResults[ElEPDiff] = _el.epDiff < maxEP[iId];
      idResults[ElMissingHits] = _el.nMissingHits <= maxMissingHits[iId];
      idResults[ElConvVeto] = _el.passConversionVeto;

    }

    if(_results) *_results = idResults;

    return idResults == elReferences[_idType];
  }

  /*static*/
  bool
  ObjectSelector::isGoodMuon(MuonVars const& _mu, MuonId _idType, std::bitset<nMuonCriteria>* _results/* = 0*/)
  {
    std::bitset<nMuonCriteria> idResults;

    if(_idType == MuLoose12){

      idResults[MuFiducial] = _mu.iSubdet != -1;
      idResults[MuPFMuon] = _mu.isPFMuon;
      idResults[MuGlobalOrTrackerMuon] = _mu.hasInnerTrack || _mu.hasGlobalTrack;

    }
    else if(_idType == MuTight12){

      idResults[MuFiducial] = _mu.iSubdet != -1;
      idResults[MuGlobalMuon] = _mu.isGlobalMuon;
      idResults[MuPFMuon] = _mu.pt > 200. || _mu.isPFMuon;
      idResults[MuMatchedStations] = _mu.nMatchedStations > 1;
      //idResults[MuLayersWithMmt] = _mu.nLayersWithMmt > (_mu.pt < 200. ? 5 : 8);
      idResults[MuLayersWithMmt] = _mu.pt > 200. || ( _mu.nLayersWithMmt >  5 ) ;
      idResults[MuNormChi2] = _mu.pt > 200. || (_mu.normChi2 > 0. && _mu.normChi2 < 10.);
      idResults[MuValidMuonHits] = _mu.nValidMuonHits > 0;
      idResults[MuDxy] = _mu.dxy > 0. && _mu.dxy < 0.2;
      idResults[MuDz] = _mu.dz > 0. && _mu.dz < 0.5;
      idResults[MuValidPixelHits] = _mu.nValidPixelHits > 0;
      idResults[MuCombIso] = _mu.combRelIso < 0.12;
      idResults[MudpTpT] = _mu.pt < 200. || ( _mu.dpTpT < 0.3) ;

    }

    if(_results) *_results = idResults;

    return idResults == muReferences[_idType];
  }

  /*static*/
  bool
  ObjectSelector::isGoodJet(JetVars const& _jet, JetId _idType, std::bitset<nJetCriteria>* _results/* = 0*/)
  {
    std::bitset<nJetCriteria> idResults;

    if(_idType == JtLoose){

      idResults[JtFiducial] = std::abs(_jet.eta) < 2.6;
      idResults[JtCHFraction] = _jet.chFraction > 0.;
      idResults[JtNHFraction] = _jet.nhFraction < 0.99;
      idResults[JtCEFraction] = _jet.ceFraction < 0.99;
      idResults[JtNEFraction] = _jet.neFraction < 0.99;
      idResults[JtNConstituents] = _jet.nConstituents > 1;
      idResults[JtChargedMultiplicity] = _jet.nCharged > 0;

    }

    if(_results) *_results = idResults;

    return idResults == jtReferences[_idType];
  }

  /*static*/
  bool
  ObjectSelector::isGoodVertex(VertexVars const& _vertex)
  {
    return _vertex.ndof >= 4 && std::abs(_vertex.z) < 24. && _vertex.rho < 2.;
  }

  ObjectSelector selector;
}
