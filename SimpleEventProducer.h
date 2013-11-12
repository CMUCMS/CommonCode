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

namespace susy {

  class SimpleEventProducer {
  public:
    struct EventVars {
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
      unsigned gen_size;
      unsigned short gen_status[NMAX];
      short gen_charge[NMAX];
      short gen_motherIndex[NMAX];
      int gen_pdgId[NMAX];
      float gen_vx[NMAX];
      float gen_vy[NMAX];
      float gen_vz[NMAX];
      float gen_pt[NMAX];
      float gen_eta[NMAX];
      float gen_phi[NMAX];
      float gen_mass[NMAX];
      float gen_px[NMAX];
      float gen_py[NMAX];
      float gen_pz[NMAX];
      float gen_energy[NMAX];
      unsigned pf_size;
      short pf_charge[NMAX];
      bool pf_isPU[NMAX];
      short pf_pdgId[NMAX];
      float pf_vx[NMAX];
      float pf_vy[NMAX];
      float pf_vz[NMAX];
      float pf_pt[NMAX];
      float pf_eta[NMAX];
      float pf_phi[NMAX];
      float pf_mass[NMAX];
      float pf_px[NMAX];
      float pf_py[NMAX];
      float pf_pz[NMAX];
      float pf_energy[NMAX];
      bool passMetFilters;
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

      void bookBranches(TTree&, bool, bool = true, bool = true, bool = true, bool = true);
      void setAddress(TTree&);

      ~AdditionalObjVars();
    };

    SimpleEventProducer();
    ~SimpleEventProducer();

    void initialize(TTree*, TTree*, TTree*, bool);

    void extractTriggerObjects(TriggerEvent&);
    void produce(Event const&);

    void addHLTPath(TString const&);
    void addHLTEventFilter(TString const&);
    void addHLTPhotonFilter(TString const&);
    void addHLTElectronFilter(TString const&);
    void addHLTMuonFilter(TString const&);
    void addGridParam(TString const&);

    void addPreselected(TTree&, bool, std::vector<unsigned> const*, std::vector<unsigned> const*, std::vector<unsigned> const*, std::vector<unsigned> const*, std::vector<unsigned> const*);

    void setSavePF(bool _val) { savePF_ = _val; }

    void setPhotonId(PhotonId _id) { photonId_ = _id; }
    void setElectronId(ElectronId _id) { electronId_ = _id; }
    void setMuonId(MuonId _id) { muonId_ = _id; }
    void setJetId(JetId _id) { jetId_ = _id; }

    ObjectTree const& getSelectedObjects() const { return selectedObjects_; }
    ObjectTree const& getAllObjects() const { return allObjects_; }
    ObjectTree const& getPreselectedObjects(unsigned _iPre) const { return *preselectedObjects_.at(_iPre); }
    AdditionalObjVars const& getSelectedAdd() const { return selectedAdd_; }
    AdditionalObjVars const& getAllAdd() const { return allAdd_; }
    AdditionalObjVars const& getPreselectedAdd(unsigned _iPre) const { return *preselectedAdd_.at(_iPre); }

  private:
    EventVars eventVars_;
    AdditionalObjVars selectedAdd_;
    AdditionalObjVars allAdd_;
    ObjectTree selectedObjects_;
    ObjectTree allObjects_;

    std::vector<AdditionalObjVars*> preselectedAdd_;
    std::vector<ObjectTree*> preselectedObjects_;

    bool saveSelected_;
    bool saveAll_;
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
  };

}

#endif
