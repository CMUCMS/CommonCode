#ifndef SimpleEventProducer_h
#define SimpleEventProducer_h

#include "ObjectTree.h"
#include "ObjectSelector.h"

#include "SusyEvent.h"
#include "SusyTriggerEvent.h"

#include <vector>
#include <map>

#include "TTree.h"
#include "TString.h"
#include "TH1.h"

namespace susy {

  unsigned const NMAXGEN(2048);
  unsigned const NMAXPF(256);

  class SimpleEventProducer {
  public:
    struct EventVars {
      unsigned runNumber;
      unsigned lumiNumber;
      unsigned eventNumber;
      float met;
      float metPhi;
      float ht;
      float mht;
      float mhtPhi;
      float bsx;
      float bsy;
      float bsz;
      float genMet;
      float genMetPhi;
      float mtElectron;
      float mtElectronPhoton;
      float mtMuon;
      float mtMuonPhoton;
      float enuMomentum;
      float munuMomentum;
      float rho;
      unsigned gen_size;
      unsigned short gen_status[NMAXGEN];
      short gen_charge[NMAXGEN];
      short gen_motherIndex[NMAXGEN];
      int gen_pdgId[NMAXGEN];
      float gen_vx[NMAXGEN];
      float gen_vy[NMAXGEN];
      float gen_vz[NMAXGEN];
      float gen_pt[NMAXGEN];
      float gen_eta[NMAXGEN];
      float gen_phi[NMAXGEN];
      float gen_mass[NMAXGEN];
      float gen_px[NMAXGEN];
      float gen_py[NMAXGEN];
      float gen_pz[NMAXGEN];
      float gen_energy[NMAXGEN];
      unsigned pf_size;
      short pf_charge[NMAXPF];
      bool pf_isPU[NMAXPF];
      short pf_pdgId[NMAXPF];
      float pf_vx[NMAXPF];
      float pf_vy[NMAXPF];
      float pf_vz[NMAXPF];
      float pf_pt[NMAXPF];
      float pf_eta[NMAXPF];
      float pf_phi[NMAXPF];
      float pf_mass[NMAXPF];
      float pf_px[NMAXPF];
      float pf_py[NMAXPF];
      float pf_pz[NMAXPF];
      float pf_energy[NMAXPF];
      bool passMetFilters;
      float puWeight;
      std::map<TString, bool> hltBits;
      std::map<TString, bool> hltFilterBits;
      std::map<TString, float> gridParams;

      void bookBranches(TTree&, bool, bool);
      void setAddress(TTree&);
    };

    struct AdditionalObjVars {
      float ph_dRGen[NMAX];
      float ph_genIso[NMAX];
      int ph_nearestGen[NMAX];
      float ph_dRJet[NMAX];
      float ph_dRNextJet[NMAX];
      float ph_dRPF[NMAX];
      short ph_nearestPF[NMAX];
      bool ph_pfIsPU[NMAX];
      float el_dRGen[NMAX];
      float el_genIso[NMAX];
      int el_nearestGen[NMAX];
      float el_dRJet[NMAX];
      float el_dRNextJet[NMAX];
      float el_dRPhoton[NMAX];
      float el_dRNextPhoton[NMAX];
      float el_dRPF[NMAX];
      short el_nearestPF[NMAX];
      bool el_pfIsPU[NMAX];
      float mu_dRGen[NMAX];
      float mu_genIso[NMAX];
      int mu_nearestGen[NMAX];
      float mu_dRJet[NMAX];
      float mu_dRNextJet[NMAX];
      float mu_dRPhoton[NMAX];
      float mu_dRNextPhoton[NMAX];
      float mu_dRPF[NMAX];
      short mu_nearestPF[NMAX];
      bool mu_pfIsPU[NMAX];
      float jt_dRGen[NMAX];
      float jt_genSumPt[NMAX];
      int jt_nearestGen[NMAX];

      std::map<TString, bool*> ph_matchHLTObj;
      std::map<TString, bool*> el_matchHLTObj;
      std::map<TString, bool*> mu_matchHLTObj;

      void setHLTObjFilters(unsigned, std::map<TString, TriggerObjectCollection> const&);
      void bookBranches(TTree&, bool, bool = true, bool = true, bool = true, bool = true);
      void setAddress(TTree&);

      ~AdditionalObjVars();
    };

    SimpleEventProducer();
    ~SimpleEventProducer();

    void initialize(TTree*, TTree*, TTree*, TH1 const*);

    void extractTriggerObjects(TriggerEvent&);
    void produce(Event const&);

    void addHLTPath(TString const& _path) { eventVars_.hltBits[_path] = false; }
    void addHLTEventFilter(TString const& _filter) { eventVars_.hltFilterBits[_filter] = false; }
    void addHLTPhotonFilter(TString const& _filter) { photonHLTObjects_[_filter]; }
    void addHLTElectronFilter(TString const& _filter) { electronHLTObjects_[_filter]; }
    void addHLTMuonFilter(TString const& _filter) { muonHLTObjects_[_filter]; }
    void addGridParam(TString const& _param) { eventVars_.gridParams[_param] = 0.; }

    void addPreselected(TTree&, std::vector<unsigned> const*, std::vector<unsigned> const*, std::vector<unsigned> const*, std::vector<unsigned> const*, std::vector<unsigned> const*);

    void setSavePF(bool _val) { savePF_ = _val; }

    void setPhotonId(PhotonId _id) { photonId_ = _id; }
    void setElectronId(ElectronId _id) { electronId_ = _id; }
    void setMuonId(MuonId _id) { muonId_ = _id; }
    void setJetId(JetId _id) { jetId_ = _id; }

    std::vector<const Photon*> sortPhotons(PhotonCollection const&, std::vector<unsigned>* = 0) const;
    std::vector<const Electron*> sortElectrons(ElectronCollection const&, std::vector<unsigned>* = 0) const;
    std::vector<const Muon*> sortMuons(MuonCollection const&, std::vector<unsigned>* = 0) const;
    std::vector<const PFJet*> sortJets(PFJetCollection const&, std::vector<unsigned>* = 0) const;

    ObjectTree const* getSelectedObjects() const { return selectedObjects_; }
    ObjectTree const* getAllObjects() const { return allObjects_; }
    ObjectTree const* getPreselectedObjects(unsigned _iPre) const { return preselectedObjects_.at(_iPre); }
    AdditionalObjVars const* getSelectedAdd() const { return selectedAdd_; }
    AdditionalObjVars const* getAllAdd() const { return allAdd_; }
    AdditionalObjVars const* getPreselectedAdd(unsigned _iPre) const { return preselectedAdd_.at(_iPre); }

    bool isMC() const { return puWeights_ != 0; }

  private:
    EventVars eventVars_;
    AdditionalObjVars* selectedAdd_;
    AdditionalObjVars* allAdd_;
    ObjectTree* selectedObjects_;
    ObjectTree* allObjects_;

    std::vector<AdditionalObjVars*> preselectedAdd_;
    std::vector<ObjectTree*> preselectedObjects_;

    bool isRealData_;
    bool savePF_;

    PhotonId photonId_;
    ElectronId electronId_;
    MuonId muonId_;
    JetId jetId_;

    std::vector<std::vector<unsigned> const*> photonPreselection_;
    std::vector<std::vector<unsigned> const*> electronPreselection_;
    std::vector<std::vector<unsigned> const*> muonPreselection_;
    std::vector<std::vector<unsigned> const*> jetPreselection_;
    std::vector<std::vector<unsigned> const*> vtxPreselection_;

    std::map<TString, TriggerObjectCollection> photonHLTObjects_;
    std::map<TString, TriggerObjectCollection> electronHLTObjects_;
    std::map<TString, TriggerObjectCollection> muonHLTObjects_;

    TH1* puWeights_;
  };

}

#endif
