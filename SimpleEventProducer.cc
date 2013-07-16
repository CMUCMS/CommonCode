#include "SimpleEventProducer.h"
#include "Utilities.h"
#include "ObjectSelector.h"

#include "TMath.h"

#include <iostream>

float const mW(80.39);

namespace susy {

  SimpleEventProducer::SimpleEventProducer() :
    eventVars_(),
    selectedAdd_(),
    allAdd_(),
    selectedObjects_(),
    allObjects_(),
    preselectedAdd_(0),
    preselectedObjects_(0),
    saveSelected_(false),
    saveAll_(false),
    savePF_(false),
    photonId_(PhLoose12),
    electronId_(ElMedium12),
    muonId_(MuTight12),
    jetId_(JtLoose),
    photonPreselection_(0),
    electronPreselection_(0),
    muonPreselection_(0),
    jetPreselection_(0),
    vtxPreselection_(0)
  {
  }

  SimpleEventProducer::~SimpleEventProducer()
  {
    for(unsigned iPre(0); iPre != preselectedObjects_.size(); ++iPre){
      delete preselectedAdd_[iPre];
      delete preselectedObjects_[iPre];
    }
  }

  void
  SimpleEventProducer::initialize(TTree* _evtTree, TTree* _selectedObjTree, TTree* _allObjTree, bool _isRealData)
  {
    eventVars_.bookBranches(*_evtTree, _isRealData, savePF_);
    _evtTree->SetAutoSave(10000000);

    if(_selectedObjTree){
      saveSelected_ = true;

      selectedObjects_.setOutput(*_selectedObjTree);
      selectedAdd_.bookBranches(*_selectedObjTree, _isRealData);
      _selectedObjTree->SetAutoSave(10000000);
    }

    if(_allObjTree){
      saveAll_ = true;

      allObjects_.setOutput(*_allObjTree);
      allAdd_.bookBranches(*_allObjTree, _isRealData);
      _allObjTree->SetAutoSave(10000000);
    }
  }

  void
  SimpleEventProducer::setHLTPaths(std::vector<TString> const& _paths)
  {
    eventVars_.hltBits.clear();
    for(unsigned iP(0); iP != _paths.size(); ++iP)
      eventVars_.hltBits[_paths[iP]] = false;
  }

  void
  SimpleEventProducer::setGridParams(std::vector<TString> const& _params)
  {
    eventVars_.gridParams.clear();
    for(unsigned iG(0); iG != _params.size(); ++iG)
      eventVars_.gridParams[_params[iG]] = 0.;
  }

  void
  SimpleEventProducer::addPreselected(TTree& _tree, bool _isRealData, std::vector<unsigned> const* _photon, std::vector<unsigned> const* _electron, std::vector<unsigned> const* _muon, std::vector<unsigned> const* _jet, std::vector<unsigned> const* _vtx)
  {
    bool addPhoton(_photon != 0);
    bool addElectron(_electron != 0);
    bool addMuon(_muon != 0);
    bool addJet(_jet != 0);
    bool addVertex(_vtx != 0);

    preselectedObjects_.push_back(new ObjectTree);
    preselectedAdd_.push_back(new AdditionalObjVars);

    preselectedObjects_.back()->setOutput(_tree, addPhoton, addElectron, addMuon, addJet, addVertex);
    preselectedAdd_.back()->bookBranches(_tree, _isRealData, addPhoton, addElectron, addMuon, addJet);
    _tree.SetAutoSave(10000000);

    photonPreselection_.push_back(_photon);
    electronPreselection_.push_back(_electron);
    muonPreselection_.push_back(_muon);
    jetPreselection_.push_back(_jet);
    vtxPreselection_.push_back(_vtx);
  }

  void
  SimpleEventProducer::produce(Event const& _event)
  {
    using namespace susy;

    /// TRIGGER BITS

    for(std::map<TString, bool>::iterator bItr(eventVars_.hltBits.begin()); bItr != eventVars_.hltBits.end(); ++bItr)
      bItr->second = _event.hltMap.pass(bItr->first + "_v*");

    if(saveAll_) allObjects_.initEvent(_event);
    if(saveSelected_) selectedObjects_.initEvent(_event);

    unsigned nPre(preselectedObjects_.size());
    for(unsigned iPre(0); iPre != nPre; ++iPre)
      preselectedObjects_[iPre]->initEvent(_event);

    PhotonCollection const& photonsSource(_event.photons.find("photons")->second);
    ElectronCollection const& electronsSource(_event.electrons.find("gsfElectrons")->second);
    MuonCollection const& muonsSource(_event.muons.find("muons")->second);
    PFJetCollection const& jetsSource(_event.pfJets.find("ak5")->second);

    unsigned nP(photonsSource.size());
    unsigned nE(electronsSource.size());
    unsigned nM(muonsSource.size());
    unsigned nJ(jetsSource.size());

    std::vector<Photon const*> photons(nP, 0);
    std::vector<Electron const*> electrons(nE, 0);
    std::vector<Muon const*> muons(nM, 0);
    std::vector<PFJet const*> jets(nJ, 0);
    std::vector<unsigned> phIndices(nP);
    std::vector<unsigned> elIndices(nE);
    std::vector<unsigned> muIndices(nM);
    std::vector<unsigned> jtIndices(nJ);

    for(unsigned iP(0); iP != nP; ++iP){
      photons[iP] = &photonsSource[iP];
      phIndices[iP] = iP;
    }
    if(!_event.isRealData){
      // FIX FOR 52X FASTSIM BUG (DUPLICATE ELECTRONS)
      // FIX FOR NAN-MOMENTUM ELECTRONS
      electrons.clear();
      elIndices.clear();
      std::set<short> superClusterIndices;
      std::set<short> trackIndices;
      for(unsigned iE(0); iE != nE; ++iE){
        Electron const& el(electronsSource[iE]);
        if(superClusterIndices.find(el.superClusterIndex) != superClusterIndices.end() ||
           trackIndices.find(el.gsfTrackIndex) != trackIndices.end())
          continue;

        if(el.momentum.X() != el.momentum.X()) continue;

        superClusterIndices.insert(el.superClusterIndex);
        trackIndices.insert(el.gsfTrackIndex);

        electrons.push_back(&el);
        elIndices.push_back(iE);
      }

      nE = electrons.size();
    }
    else{
      for(unsigned iE(0); iE != nE; ++iE){
        electrons[iE] = &electronsSource[iE];
        elIndices[iE] = iE;
      }
    }
    for(unsigned iM(0); iM != nM; ++iM){
      muons[iM] = &muonsSource[iM];
      muIndices[iM] = iM;
    }
    for(unsigned iJ(0); iJ != nJ; ++iJ){
      jets[iJ] = &jetsSource[iJ];
      jtIndices[iJ] = iJ;
    }

    sortByPt(photons, &phIndices);
    sortByPt(electrons, &elIndices);
    sortByPt(muons, &muIndices);
    sortByPt(jets, &jtIndices);


    /// GEN PARTICLES

    std::vector<Particle const*> fsParticles;

    if(!_event.isRealData){
      TVector2 genMetV;

      eventVars_.gen_size = _event.genParticles.size();
      if(eventVars_.gen_size >= NMAX)
        throw std::runtime_error("Too many GenParticles");

      for(unsigned iG(0); iG != eventVars_.gen_size; ++iG){
        Particle const& particle(_event.genParticles[iG]);

        eventVars_.gen_status[iG] = particle.status;
        eventVars_.gen_charge[iG] = particle.charge;
        eventVars_.gen_motherIndex[iG] = particle.motherIndex;
        eventVars_.gen_pdgId[iG] = particle.pdgId;
        eventVars_.gen_vx[iG] = particle.vertex.X();
        eventVars_.gen_vy[iG] = particle.vertex.Y();
        eventVars_.gen_vz[iG] = particle.vertex.Z();
        eventVars_.gen_pt[iG] = particle.momentum.Pt();
        if(eventVars_.gen_pt[iG] > 0.001){
          eventVars_.gen_eta[iG] = particle.momentum.Eta();
          eventVars_.gen_phi[iG] = particle.momentum.Phi();
        }
        else{
          eventVars_.gen_eta[iG] = 0.;
          eventVars_.gen_phi[iG] = 0.;
        }
        eventVars_.gen_mass[iG] = particle.momentum.M();
        eventVars_.gen_px[iG] = particle.momentum.X();
        eventVars_.gen_py[iG] = particle.momentum.Y();
        eventVars_.gen_pz[iG] = particle.momentum.Z();
        eventVars_.gen_energy[iG] = particle.momentum.E();

        if(particle.status != 1) continue;
       
        if((std::abs(particle.pdgId) / 100) % 10 == 0 && particle.pdgId != 22 && particle.charge == 0)
          genMetV += TVector2(particle.momentum.X(), particle.momentum.Y());

        if(particle.momentum.Pt() > 2.)
          fsParticles.push_back(&particle);
      }

      eventVars_.genMet = genMetV.Mod();
      eventVars_.genMetPhi = TVector2::Phi_mpi_pi(genMetV.Phi());
    }
    unsigned nG(fsParticles.size());


    /// PF PARTICLES

    std::map<std::pair<double, double>, PFParticle const*> uniqueParticles; // bug fix for tag cms533v0 / cms538v1

    unsigned nPF(_event.pfParticles.size());
    for(unsigned iP(0); iP != nPF; ++iP){
      PFParticle const& particle(_event.pfParticles[iP]);
      double pt(particle.momentum.Pt());
      if(std::abs(particle.momentum.Eta()) > etaMax) continue;

      uniqueParticles[std::pair<double, double>(pt, particle.momentum.Eta())] = &particle;
    }

    std::vector<PFParticle const*> pfParticles;
    std::map<std::pair<double, double>, PFParticle const*>::iterator uEnd(uniqueParticles.end());
    for(std::map<std::pair<double, double>, PFParticle const*>::iterator uItr(uniqueParticles.begin()); uItr != uEnd; ++uItr)
      pfParticles.push_back(uItr->second);

    nPF = pfParticles.size();

    if(savePF_){
      unsigned nSavedPF(0);
      for(unsigned iPF(0); iPF != nPF; ++iPF){
        PFParticle const& particle(*pfParticles[iPF]);

        if(particle.momentum.Pt() < 3.) continue;

        if(nSavedPF >= NMAX)
          throw std::runtime_error("Too many PFParticles");

        eventVars_.pf_charge[nSavedPF] = particle.charge;
        eventVars_.pf_isPU[nSavedPF] = particle.isPU;
        eventVars_.pf_pdgId[nSavedPF] = particle.pdgId;
        eventVars_.pf_vx[nSavedPF] = particle.vertex.X();
        eventVars_.pf_vy[nSavedPF] = particle.vertex.Y();
        eventVars_.pf_vz[nSavedPF] = particle.vertex.Z();
        eventVars_.pf_pt[nSavedPF] = particle.momentum.Pt();
        if(eventVars_.pf_pt[nSavedPF] > 0.001){
          eventVars_.pf_eta[nSavedPF] = particle.momentum.Eta();
          eventVars_.pf_phi[nSavedPF] = particle.momentum.Phi();
        }
        else{
          eventVars_.pf_eta[nSavedPF] = 0.;
          eventVars_.pf_phi[nSavedPF] = 0.;
        }
        eventVars_.pf_mass[nSavedPF] = particle.momentum.M();
        eventVars_.pf_px[nSavedPF] = particle.momentum.X();
        eventVars_.pf_py[nSavedPF] = particle.momentum.Y();
        eventVars_.pf_pz[nSavedPF] = particle.momentum.Z();
        eventVars_.pf_energy[nSavedPF] = particle.momentum.E();

        ++nSavedPF;
      }
      eventVars_.pf_size = nSavedPF;
    }

    /// JETS - SELECTION ONLY

    std::vector<bool> isGoodJet(nJ, false);
    std::vector<bool> isMatchedJet(nJ, false);

    for(unsigned iJ(0); iJ != nJ; ++iJ){
      PFJet const& jet(*jets[iJ]);

      JetVars vars(jet, _event);

      isGoodJet[iJ] = ObjectSelector::isGoodJet(vars, jetId_);
    }


    /// PHOTONS

    std::vector<bool> isGoodPhoton(nP, false);
    TLorentzVector const* leadingPhoton(0);

    for(unsigned iP(0); iP != nP; ++iP){
      Photon const& photon(*photons[iP]);

      PhotonVars vars(photon, _event);

      isGoodPhoton[iP] = ObjectSelector::isGoodPhoton(vars, photonId_);

      bool processObject(saveAll_ || isGoodPhoton[iP]);
      std::vector<bool> isPreselected(nPre, false);
      for(unsigned iPre(0); iPre != nPre; ++iPre){
        if(!photonPreselection_[iPre]) continue;
        isPreselected[iPre] = std::find(photonPreselection_[iPre]->begin(), photonPreselection_[iPre]->end(), phIndices[iP]) != photonPreselection_[iPre]->end();
        processObject = processObject || isPreselected[iPre];
      }

      if(!processObject) continue;

      double dRGenMin(0.1);
      unsigned iMatchedGen(-1);
      double genIso(0.);
      for(unsigned iG(0); iG != nG; ++iG){
        Particle const& particle(*fsParticles[iG]);
        double dRGen(particle.momentum.DeltaR(photon.momentum));
        if(dRGen < dRGenMin){
          dRGenMin = dRGen;
          iMatchedGen = iG;
        }
        if(dRGen < 0.3)
          genIso += particle.momentum.Pt();
      }
      int genPdgId(0);
      if(iMatchedGen != unsigned(-1)){
        genIso -= fsParticles[iMatchedGen]->momentum.Pt();
        genPdgId = fsParticles[iMatchedGen]->pdgId;
      }
      else
        dRGenMin = -1.;

      double dRJetMin(-1.);
      unsigned iMatchedJet(-1);
      for(unsigned iJ(0); iJ != nJ; ++iJ){
        if(!isGoodJet[iJ]) continue;
        double dRJet(jets[iJ]->momentum.DeltaR(photon.momentum));
        if(dRJetMin < 0. || dRJet < dRJetMin){
          dRJetMin = dRJet;
          if(dRJet < 0.3) iMatchedJet = iJ;
        }
      }

      double dRPFMin(0.1);
      unsigned iMatchedPF(-1);
      for(unsigned iPF(0); iPF != nPF; ++iPF){
        double dRPF(pfParticles[iPF]->momentum.DeltaR(photon.momentum));
        if(dRPF < dRPFMin){
          dRPFMin = dRPF;
          iMatchedPF = iPF;
        }
      }
      short pfPdgId(0);
      bool isPU(false);
      if(iMatchedPF != unsigned(-1)){
        pfPdgId = pfParticles[iMatchedPF]->pdgId;
        isPU = pfParticles[iMatchedPF]->isPU;
      }
      else
        dRPFMin = -1.;

      if(isGoodPhoton[iP]){
        if(iMatchedJet != unsigned(-1)) isMatchedJet[iMatchedJet] = true;
        if(!leadingPhoton) leadingPhoton = &photon.momentum;
      }

      if(saveSelected_ && isGoodPhoton[iP]){
        unsigned iSel(selectedObjects_.getPhotonSize());
        selectedObjects_.save(vars);
        selectedAdd_.ph_dRGen[iSel] = dRGenMin;
        selectedAdd_.ph_genIso[iSel] = genIso;
        selectedAdd_.ph_nearestGen[iSel] = genPdgId;
        selectedAdd_.ph_dRJet[iSel] = dRJetMin;
        selectedAdd_.ph_dRPF[iSel] = dRPFMin;
        selectedAdd_.ph_nearestPF[iSel] = pfPdgId;
        selectedAdd_.ph_pfIsPU[iSel] = isPU;
      }
      
      if(saveAll_){
        allObjects_.save(vars);
        allAdd_.ph_dRGen[iP] = dRGenMin;
        allAdd_.ph_genIso[iP] = genIso;
        allAdd_.ph_nearestGen[iP] = genPdgId;
        allAdd_.ph_dRJet[iP] = dRJetMin;
        allAdd_.ph_dRPF[iP] = dRPFMin;
        allAdd_.ph_nearestPF[iP] = pfPdgId;
        allAdd_.ph_pfIsPU[iP] = isPU;
      }

      for(unsigned iPre(0); iPre != nPre; ++iPre){
        if(!isPreselected[iPre]) continue;
        unsigned iSel(preselectedObjects_[iPre]->getPhotonSize());
        preselectedObjects_[iPre]->save(vars);
        preselectedAdd_[iPre]->ph_dRGen[iSel] = dRGenMin;
        preselectedAdd_[iPre]->ph_genIso[iSel] = genIso;
        preselectedAdd_[iPre]->ph_nearestGen[iSel] = genPdgId;
        preselectedAdd_[iPre]->ph_dRJet[iSel] = dRJetMin;
        preselectedAdd_[iPre]->ph_dRPF[iSel] = dRPFMin;
        preselectedAdd_[iPre]->ph_nearestPF[iSel] = pfPdgId;
        preselectedAdd_[iPre]->ph_pfIsPU[iSel] = isPU;
      }
    }


    /// ELECTRONS

    std::vector<bool> isGoodElectron(nE, false);
    TLorentzVector const* leadingElectron(0);

    for(unsigned iE(0); iE != nE; ++iE){
      Electron const& electron(*electrons[iE]);

      ElectronVars vars(electron, _event);

      isGoodElectron[iE] = ObjectSelector::isGoodElectron(vars, electronId_);

      bool processObject(saveAll_ || isGoodElectron[iE]);
      std::vector<bool> isPreselected(nPre, false);
      for(unsigned iPre(0); iPre != nPre; ++iPre){
        if(!electronPreselection_[iPre]) continue;
        isPreselected[iPre] = std::find(electronPreselection_[iPre]->begin(), electronPreselection_[iPre]->end(), elIndices[iE]) != electronPreselection_[iPre]->end();
        processObject = processObject || isPreselected[iPre];
      }

      if(!processObject) continue;

      double dRGenMin(0.1);
      unsigned iMatchedGen(-1);
      double genIso(0.);
      for(unsigned iG(0); iG != nG; ++iG){
        Particle const& particle(*fsParticles[iG]);
        double dRGen(particle.momentum.DeltaR(electron.momentum));
        if(dRGen < dRGenMin){
          dRGenMin = dRGen;
          iMatchedGen = iG;
        }
        if(dRGen < 0.3)
          genIso += particle.momentum.Pt();
      }
      int genPdgId(0);
      if(iMatchedGen != unsigned(-1)){
        genIso -= fsParticles[iMatchedGen]->momentum.Pt();
        genPdgId = fsParticles[iMatchedGen]->pdgId;
      }
      else
        dRGenMin = -1.;

      double dRPhotonMin(-1.);
      unsigned iMatchedPhoton(-1);
      for(unsigned iP(0); iP != nP; ++iP){
        double dRPhoton(photons[iP]->momentum.DeltaR(electron.momentum));
        if(dRPhotonMin < 0. || dRPhoton < dRPhotonMin){
          dRPhotonMin = dRPhoton;
          if(dRPhoton < 0.2) iMatchedPhoton = iP;
        }
      }

      double dRNextPhotonMin(-1.);
      for(unsigned iP(0); iP != nP; ++iP){
        if(!isGoodPhoton[iP] || iP == iMatchedPhoton) continue;
        double dRPhoton(photons[iP]->momentum.DeltaR(electron.momentum));
        if(dRNextPhotonMin < 0. || dRPhoton < dRNextPhotonMin) dRNextPhotonMin = dRPhoton;
      }

      double dRJetMin(-1.);
      unsigned iMatchedJet(-1);
      for(unsigned iJ(0); iJ != nJ; ++iJ){
        if(!isGoodJet[iJ]) continue;
        double dRJet(jets[iJ]->momentum.DeltaR(electron.momentum));
        if(dRJetMin < 0. || dRJet < dRJetMin){
          dRJetMin = dRJet;
          if(dRJet < 0.3) iMatchedJet = iJ;
        }
      }

      double dRPFMin(0.1);
      unsigned iMatchedPF(-1);
      for(unsigned iPF(0); iPF != nPF; ++iPF){
        double dRPF(pfParticles[iPF]->momentum.DeltaR(electron.momentum));
        if(dRPF < dRPFMin){
          dRPFMin = dRPF;
          iMatchedPF = iPF;
        }
      }
      short pfPdgId(0);
      bool isPU(false);
      if(iMatchedPF != unsigned(-1)){
        pfPdgId = pfParticles[iMatchedPF]->pdgId;
        isPU = pfParticles[iMatchedPF]->isPU;
      }
      else
        dRPFMin = -1.;

      if(isGoodElectron[iE]){
        if(iMatchedJet != unsigned(-1)) isMatchedJet[iMatchedJet] = true;
        if(!leadingElectron) leadingElectron = &electron.momentum;
      }

      if(saveSelected_ && isGoodElectron[iE]){
        unsigned iSel(selectedObjects_.getElectronSize());
        selectedObjects_.save(vars);
        selectedAdd_.el_dRGen[iSel] = dRGenMin;
        selectedAdd_.el_genIso[iSel] = genIso;
        selectedAdd_.el_nearestGen[iSel] = genPdgId;
        selectedAdd_.el_dRPhoton[iSel] = dRPhotonMin;
        selectedAdd_.el_dRNextPhoton[iSel] = dRNextPhotonMin;
        selectedAdd_.el_dRJet[iSel] = dRJetMin;
        selectedAdd_.el_dRPF[iSel] = dRPFMin;
        selectedAdd_.el_nearestPF[iSel] = pfPdgId;
        selectedAdd_.el_pfIsPU[iSel] = isPU;
      }

      if(saveAll_){
        allObjects_.save(vars);
        allAdd_.el_dRGen[iE] = dRGenMin;
        allAdd_.el_genIso[iE] = genIso;
        allAdd_.el_nearestGen[iE] = genPdgId;
        allAdd_.el_dRPhoton[iE] = dRPhotonMin;
        allAdd_.el_dRNextPhoton[iE] = dRNextPhotonMin;
        allAdd_.el_dRJet[iE] = dRJetMin;
        allAdd_.el_dRPF[iE] = dRPFMin;
        allAdd_.el_nearestPF[iE] = pfPdgId;
        allAdd_.el_pfIsPU[iE] = isPU;
      }

      for(unsigned iPre(0); iPre != nPre; ++iPre){
        if(!isPreselected[iPre]) continue;
        unsigned iSel(preselectedObjects_[iPre]->getElectronSize());
        preselectedObjects_[iPre]->save(vars);
        preselectedAdd_[iPre]->el_dRGen[iSel] = dRGenMin;
        preselectedAdd_[iPre]->el_genIso[iSel] = genIso;
        preselectedAdd_[iPre]->el_nearestGen[iSel] = genPdgId;
        preselectedAdd_[iPre]->el_dRPhoton[iSel] = dRPhotonMin;
        preselectedAdd_[iPre]->el_dRNextPhoton[iSel] = dRNextPhotonMin;
        preselectedAdd_[iPre]->el_dRJet[iSel] = dRJetMin;
        preselectedAdd_[iPre]->el_dRPF[iSel] = dRPFMin;
        preselectedAdd_[iPre]->el_nearestPF[iSel] = pfPdgId;
        preselectedAdd_[iPre]->el_pfIsPU[iSel] = isPU;
      }
    }


    /// MUONS

    std::vector<bool> isGoodMuon(nM, false);
    TLorentzVector const* leadingMuon(0);

    for(unsigned iM(0); iM != nM; ++iM){
      Muon const& muon(*muons[iM]);

      MuonVars vars(muon, _event);

      isGoodMuon[iM] = ObjectSelector::isGoodMuon(vars, muonId_);

      bool processObject(saveAll_ || isGoodMuon[iM]);
      std::vector<bool> isPreselected(nPre, false);
      for(unsigned iPre(0); iPre != nPre; ++iPre){
        if(!muonPreselection_[iPre]) continue;
        isPreselected[iPre] = std::find(muonPreselection_[iPre]->begin(), muonPreselection_[iPre]->end(), muIndices[iM]) != muonPreselection_[iPre]->end();
        processObject = processObject || isPreselected[iPre];
      }

      if(!processObject) continue;

      double dRGenMin(0.1);
      unsigned iMatchedGen(-1);
      double genIso(0.);
      for(unsigned iG(0); iG != nG; ++iG){
        Particle const& particle(*fsParticles[iG]);
        double dRGen(particle.momentum.DeltaR(muon.momentum));
        if(dRGen < dRGenMin){
          dRGenMin = dRGen;
          iMatchedGen = iG;
        }
        if(dRGen < 0.3)
          genIso += particle.momentum.Pt();
      }
      int genPdgId(0);
      if(iMatchedGen != unsigned(-1)){
        genIso -= fsParticles[iMatchedGen]->momentum.Pt();
        genPdgId = fsParticles[iMatchedGen]->pdgId;
      }
      else
        dRGenMin = -1.;

      double dRPhotonMin(-1.);
      unsigned iMatchedPhoton(-1);
      for(unsigned iP(0); iP != nP; ++iP){
        double dRPhoton(photons[iP]->momentum.DeltaR(muon.momentum));
        if(dRPhotonMin < 0. || dRPhoton < dRPhotonMin){
          dRPhotonMin = dRPhoton;
          if(dRPhoton < 0.2) iMatchedPhoton = iP;
        }
      }

      double dRNextPhotonMin(-1.);
      for(unsigned iP(0); iP != nP; ++iP){
        if(!isGoodPhoton[iP] || iP == iMatchedPhoton) continue;
        double dRPhoton(photons[iP]->momentum.DeltaR(muon.momentum));
        if(dRNextPhotonMin < 0. || dRPhoton < dRNextPhotonMin) dRNextPhotonMin = dRPhoton;
      }

      double dRJetMin(-1.);
      unsigned iMatchedJet(-1);
      for(unsigned iJ(0); iJ != nJ; ++iJ){
        if(!isGoodJet[iJ]) continue;
        double dRJet(jets[iJ]->momentum.DeltaR(muon.momentum));
        if(dRJetMin < 0. || dRJet < dRJetMin){
          dRJetMin = dRJet;
          if(dRJet < 0.3) iMatchedJet = iJ;
        }
      }

      double dRPFMin(0.1);
      unsigned iMatchedPF(-1);
      for(unsigned iPF(0); iPF != nPF; ++iPF){
        double dRPF(pfParticles[iPF]->momentum.DeltaR(muon.momentum));
        if(dRPF < dRPFMin){
          dRPFMin = dRPF;
          iMatchedPF = iPF;
        }
      }
      short pfPdgId(0);
      bool isPU(false);
      if(iMatchedPF != unsigned(-1)){
        pfPdgId = pfParticles[iMatchedPF]->pdgId;
        isPU = pfParticles[iMatchedPF]->isPU;
      }
      else
        dRPFMin = -1.;

      if(isGoodMuon[iM]){
        if(iMatchedJet != unsigned(-1)) isMatchedJet[iMatchedJet] = true;
        if(!leadingMuon) leadingMuon = &muon.momentum;
      }

      if(saveSelected_ && isGoodMuon[iM]){
        unsigned iSel(selectedObjects_.getMuonSize());
        selectedObjects_.save(vars);
        selectedAdd_.mu_dRGen[iSel] = dRGenMin;
        selectedAdd_.mu_genIso[iSel] = genIso;
        selectedAdd_.mu_nearestGen[iSel] = genPdgId;
        selectedAdd_.mu_dRPhoton[iSel] = dRPhotonMin;
        selectedAdd_.mu_dRNextPhoton[iSel] = dRNextPhotonMin;
        selectedAdd_.mu_dRJet[iSel] = dRJetMin;
        selectedAdd_.mu_dRPF[iSel] = dRPFMin;
        selectedAdd_.mu_nearestPF[iSel] = pfPdgId;
        selectedAdd_.mu_pfIsPU[iSel] = isPU;
      }

      if(saveAll_){
        allObjects_.save(vars);
        allAdd_.mu_dRGen[iM] = dRGenMin;
        allAdd_.mu_genIso[iM] = genIso;
        allAdd_.mu_nearestGen[iM] = genPdgId;
        allAdd_.mu_dRPhoton[iM] = dRPhotonMin;
        allAdd_.mu_dRNextPhoton[iM] = dRNextPhotonMin;
        allAdd_.mu_dRJet[iM] = dRJetMin;
        allAdd_.mu_dRPF[iM] = dRPFMin;
        allAdd_.mu_nearestPF[iM] = pfPdgId;
        allAdd_.mu_pfIsPU[iM] = isPU;
      }

      for(unsigned iPre(0); iPre != nPre; ++iPre){
        if(!isPreselected[iPre]) continue;
        unsigned iSel(preselectedObjects_[iPre]->getMuonSize());
        preselectedObjects_[iPre]->save(vars);
        preselectedAdd_[iPre]->mu_dRGen[iSel] = dRGenMin;
        preselectedAdd_[iPre]->mu_genIso[iSel] = genIso;
        preselectedAdd_[iPre]->mu_nearestGen[iSel] = genPdgId;
        preselectedAdd_[iPre]->mu_dRPhoton[iSel] = dRPhotonMin;
        preselectedAdd_[iPre]->mu_dRNextPhoton[iSel] = dRNextPhotonMin;
        preselectedAdd_[iPre]->mu_dRJet[iSel] = dRJetMin;
        preselectedAdd_[iPre]->mu_dRPF[iSel] = dRPFMin;
        preselectedAdd_[iPre]->mu_nearestPF[iSel] = pfPdgId;
        preselectedAdd_[iPre]->mu_pfIsPU[iSel] = isPU;
      }
    }


    /// NEAREST UNMATCHED JET

    std::vector<unsigned> iSel((saveSelected_ ? 1 : 0) + nPre, 0); // counter for selected & preselected
    for(unsigned iP(0); iP != nP; ++iP){
      bool processObject(saveAll_ || isGoodPhoton[iP]);
      std::vector<bool> isPreselected(nPre, false);
      for(unsigned iPre(0); iPre != nPre; ++iPre){
        if(!photonPreselection_[iPre]) continue;
        isPreselected[iPre] = std::find(photonPreselection_[iPre]->begin(), photonPreselection_[iPre]->end(), phIndices[iP]) != photonPreselection_[iPre]->end();
        processObject = processObject || isPreselected[iPre];
      }

      if(!processObject) continue;

      Photon const& photon(*photons[iP]);

      double dRMin(-1.);
      for(unsigned iJ(0); iJ != nJ; ++iJ){
        if(!isGoodJet[iJ] || isMatchedJet[iJ]) continue;
        double dRJet(jets[iJ]->momentum.DeltaR(photon.momentum));
        if(dRMin < 0. || dRJet < dRMin) dRMin = dRJet;
      }

      if(saveSelected_ && isGoodPhoton[iP]){
        selectedAdd_.ph_dRNextJet[iSel[0]] = dRMin;
        ++iSel[0];
      }

      if(saveAll_) allAdd_.ph_dRNextJet[iP] = dRMin;

      for(unsigned iPre(0); iPre != nPre; ++iPre){
        if(!isPreselected[iPre]) continue;
        preselectedAdd_[iPre]->ph_dRNextJet[iSel[iPre + 1]] = dRMin;
        ++iSel[iPre + 1];
      }
    }

    iSel.assign(iSel.size(), 0);
    for(unsigned iE(0); iE != nE; ++iE){
      bool processObject(saveAll_ || isGoodElectron[iE]);
      std::vector<bool> isPreselected(nPre, false);
      for(unsigned iPre(0); iPre != nPre; ++iPre){
        if(!electronPreselection_[iPre]) continue;
        isPreselected[iPre] = std::find(electronPreselection_[iPre]->begin(), electronPreselection_[iPre]->end(), elIndices[iE]) != electronPreselection_[iPre]->end();
        processObject = processObject || isPreselected[iPre];
      }

      if(!processObject) continue;

      Electron const& electron(*electrons[iE]);

      double dRMin(-1.);
      for(unsigned iJ(0); iJ != nJ; ++iJ){
        if(!isGoodJet[iJ] || isMatchedJet[iJ]) continue;
        double dRJet(jets[iJ]->momentum.DeltaR(electron.momentum));
        if(dRMin < 0. || dRJet < dRMin) dRMin = dRJet;
      }

      if(saveSelected_ && isGoodElectron[iE]){
        selectedAdd_.el_dRNextJet[iSel[0]] = dRMin;
        ++iSel[0];
      }

      if(saveAll_) allAdd_.el_dRNextJet[iE] = dRMin;

      for(unsigned iPre(0); iPre != nPre; ++iPre){
        if(!isPreselected[iPre]) continue;
        preselectedAdd_[iPre]->el_dRNextJet[iSel[iPre + 1]] = dRMin;
        ++iSel[iPre + 1];
      }
    }

    iSel.assign(iSel.size(), 0);
    for(unsigned iM(0); iM != nM; ++iM){
      bool processObject(saveAll_ || isGoodMuon[iM]);
      std::vector<bool> isPreselected(nPre, false);
      for(unsigned iPre(0); iPre != nPre; ++iPre){
        if(!muonPreselection_[iPre]) continue;
        isPreselected[iPre] = std::find(muonPreselection_[iPre]->begin(), muonPreselection_[iPre]->end(), muIndices[iM]) != muonPreselection_[iPre]->end();
        processObject = processObject || isPreselected[iPre];
      }

      if(!processObject) continue;

      Muon const& muon(*muons[iM]);

      double dRMin(-1.);
      for(unsigned iJ(0); iJ != nJ; ++iJ){
        if(!isGoodJet[iJ] || isMatchedJet[iJ]) continue;
        double dRJet(jets[iJ]->momentum.DeltaR(muon.momentum));
        if(dRMin < 0. || dRJet < dRMin) dRMin = dRJet;
      }

      if(saveSelected_ && isGoodMuon[iM]){
        selectedAdd_.mu_dRNextJet[iSel[0]] = dRMin;
        ++iSel[0];
      }

      if(saveAll_) allAdd_.mu_dRNextJet[iM] = dRMin;

      for(unsigned iPre(0); iPre != nPre; ++iPre){
        if(!isPreselected[iPre]) continue;
        preselectedAdd_[iPre]->mu_dRNextJet[iSel[iPre + 1]] = dRMin;
        ++iSel[iPre + 1];
      }
    }


    /// JETS, HT, & MHT

    eventVars_.ht = 0.;
    TVector2 mhtV;

    for(unsigned iJ(0); iJ != nJ; ++iJ){
      PFJet const& jet(*jets[iJ]);

      JetVars vars(jet, _event);

      if(isGoodJet[iJ] && !isMatchedJet[iJ] && std::abs(vars.eta) < 3. && vars.pt > 10.){
        eventVars_.ht += vars.pt;
        mhtV += TVector2(jet.momentum.X(), jet.momentum.Y());
      }

      bool processObject(saveAll_ || (isGoodJet[iJ] && !isMatchedJet[iJ]));
      std::vector<bool> isPreselected(nPre, false);
      for(unsigned iPre(0); iPre != nPre; ++iPre){
        if(!jetPreselection_[iPre]) continue;
        isPreselected[iPre] = std::find(jetPreselection_[iPre]->begin(), jetPreselection_[iPre]->end(), jtIndices[iJ]) != jetPreselection_[iPre]->end();
        processObject = processObject || isPreselected[iPre];
      }

      if(!processObject) continue;

      double dRGenMin(-1.);
      Particle const* nearest(0);
      for(unsigned iG(0); iG != nG; ++iG){
        Particle const& particle(*fsParticles[iG]);

        double dRGen(particle.momentum.DeltaR(jet.momentum));

        if(dRGenMin < 0. || dRGen < dRGenMin){
          dRGenMin = dRGen;
          nearest = &particle;
        } 
      }

      int nearestGen(0);
      TLorentzVector genJet;

      if(nearest && dRGenMin < 0.6){
        Particle const* particle(nearest);
        while(particle->motherIndex != -1){
          Particle const& mother(_event.genParticles[particle->motherIndex]);

          unsigned motherId(std::abs(mother.pdgId));

          if(motherId == 23 || motherId == 24 || motherId == 25 || mother.momentum.Pt() < 0.001 || mother.momentum.DeltaR(nearest->momentum) > 0.6){
            nearestGen = particle->pdgId;

            for(unsigned iG(0); iG != nG; ++iG){
              Particle const& fs(*fsParticles[iG]);

              if(&fs == particle || fs.momentum.DeltaR(particle->momentum) < 0.6)
                genJet += fs.momentum;
            }

            break;
          }

          particle = &mother;
        }
      }

      if(saveSelected_ && isGoodJet[iJ] && !isMatchedJet[iJ]){
        unsigned iSel(selectedObjects_.getJetSize());
        selectedObjects_.save(vars);
        selectedAdd_.jt_nearestGen[iSel] = nearestGen;
        if(genJet.Pt() > 0.001) selectedAdd_.jt_dRGen[iSel] = genJet.DeltaR(jet.momentum);
        else selectedAdd_.jt_dRGen[iSel] = -1.;
        selectedAdd_.jt_genSumPt[iSel] = genJet.Pt();
      }

      if(saveAll_){
        allObjects_.save(vars);
        allAdd_.jt_nearestGen[iJ] = nearestGen;
        if(genJet.Pt() > 0.001) allAdd_.jt_dRGen[iJ] = genJet.DeltaR(jet.momentum);
        else allAdd_.jt_dRGen[iJ] = -1.;
        allAdd_.jt_genSumPt[iJ] = genJet.Pt();
      }

      for(unsigned iPre(0); iPre != nPre; ++iPre){
        if(!isPreselected[iPre]) continue;
        unsigned iSel(preselectedObjects_[iPre]->getJetSize());
        preselectedObjects_[iPre]->save(vars);
        preselectedAdd_[iPre]->jt_nearestGen[iSel] = nearestGen;
        if(genJet.Pt() > 0.001) preselectedAdd_[iPre]->jt_dRGen[iSel] = genJet.DeltaR(jet.momentum);
        else preselectedAdd_[iPre]->jt_dRGen[iSel] = -1.;
        preselectedAdd_[iPre]->jt_genSumPt[iSel] = genJet.Pt();
      }

    }

    eventVars_.mht = mhtV.Mod();
    eventVars_.mhtPhi = TVector2::Phi_mpi_pi(TMath::Pi() - mhtV.Phi());


    /// VERTICES

    unsigned nV(_event.vertices.size());
    for(unsigned iV(0); iV != nV; ++iV){
      Vertex const& vertex(_event.vertices[iV]);

      VertexVars vars(vertex);

      if(saveSelected_ && vars.isGood)
        selectedObjects_.save(vars);

      if(saveAll_)
        allObjects_.save(vars);

      for(unsigned iPre(0); iPre != nPre; ++iPre){
        if(!vtxPreselection_[iPre]) continue;
        if(std::find(vtxPreselection_[iPre]->begin(), vtxPreselection_[iPre]->end(), iV) != vtxPreselection_[iPre]->end())
          preselectedObjects_[iPre]->save(vars);
      }
    }


    /// EVENT VARIABLES

    TVector2 const& metV(_event.metMap.find("pfType01CorrectedMet")->second.mEt);

    eventVars_.met = metV.Mod();
    eventVars_.metPhi = TVector2::Phi_mpi_pi(metV.Phi());
    eventVars_.mtElectron = -1.;
    eventVars_.mtElectronPhoton = -1.;
    eventVars_.mtMuon = -1.;
    eventVars_.mtMuonPhoton = -1.;
    eventVars_.enuMomentum = -1.;
    eventVars_.munuMomentum = -1.;

    if(leadingElectron){
      double etSum(metV.Mod() + leadingElectron->Et());
      TVector2 ptl(leadingElectron->X(), leadingElectron->Y());
      double ptSum2((metV + ptl).Mod2());

      eventVars_.mtElectron = std::sqrt(etSum * etSum - ptSum2);

      double A((mW * mW / 2. + (metV * ptl)) / leadingElectron->P());
      TVector3 pe(leadingElectron->Vect().Unit());
      if(A * A > metV.Mod2() * pe.Perp2()){
        double mez((A * pe.Z() + std::sqrt(A * A - metV.Mod2() * pe.Perp2())) / pe.Perp2());
        eventVars_.enuMomentum = std::sqrt(ptSum2 + (leadingElectron->Z() + mez) * (leadingElectron->Z() + mez));
      }
    }

    if(leadingMuon){
      double etSum(metV.Mod() + leadingMuon->Et());
      TVector2 ptl(leadingMuon->X(), leadingMuon->Y());
      double ptSum2((metV + ptl).Mod2());

      eventVars_.mtMuon = std::sqrt(etSum * etSum - ptSum2);

      double A((mW * mW / 2. + (metV * ptl)) / leadingMuon->P());
      TVector3 pmu(leadingMuon->Vect().Unit());
      if(A * A > metV.Mod2() * pmu.Perp2()){
        double mez((A * pmu.Z() + std::sqrt(A * A - metV.Mod2() * pmu.Perp2())) / pmu.Perp2());
        eventVars_.munuMomentum = std::sqrt(ptSum2 + (leadingMuon->Z() + mez) * (leadingMuon->Z() + mez));
      }
    }

    if(leadingPhoton && leadingElectron){
      double etSum(metV.Mod() + (*leadingElectron + *leadingPhoton).Et());
      double ptSum2((metV + TVector2(leadingElectron->X() + leadingPhoton->X(), leadingElectron->Y() + leadingPhoton->Y())).Mod());
      
      eventVars_.mtElectronPhoton = std::sqrt(etSum * etSum - ptSum2);
    }

    if(leadingPhoton && leadingMuon){
      double etSum(metV.Mod() + (*leadingMuon + *leadingPhoton).Et());
      double ptSum2((metV + TVector2(leadingMuon->X() + leadingPhoton->X(), leadingMuon->Y() + leadingPhoton->Y())).Mod());

      eventVars_.mtMuonPhoton = std::sqrt(etSum * etSum - ptSum2);
    }

    for(std::map<TString, float>::iterator pItr(eventVars_.gridParams.begin()); pItr != eventVars_.gridParams.end(); ++pItr){
      std::map<TString, float>::const_iterator sItr(_event.gridParams.find(pItr->first));
      pItr->second = sItr != _event.gridParams.end() ? sItr->second : 0.;
    }
  }

  void
  SimpleEventProducer::EventVars::bookBranches(TTree& _tree, bool _isRealData, bool _savePF)
  {
    _tree.Branch("met", &met, "met/F");
    _tree.Branch("metPhi", &metPhi, "metPhi/F");
    _tree.Branch("mht", &mht, "mht/F");
    _tree.Branch("mhtPhi", &mhtPhi, "mhtPhi/F");
    _tree.Branch("mtElectron", &mtElectron, "mtElectron/F");
    _tree.Branch("mtElectronPhoton", &mtElectronPhoton, "mtElectronPhoton/F");
    _tree.Branch("mtMuon", &mtMuon, "mtMuon/F");
    _tree.Branch("mtMuonPhoton", &mtMuonPhoton, "mtMuonPhoton/F");
    _tree.Branch("enuMomentum", &enuMomentum, "enuMomentum/F");
    _tree.Branch("munuMomentum", &munuMomentum, "munuMomentum/F");
    for(std::map<TString, bool>::iterator bItr(hltBits.begin()); bItr != hltBits.end(); ++bItr)
      _tree.Branch(bItr->first, &bItr->second, bItr->first + "/O");
    for(std::map<TString, float>::iterator pItr(gridParams.begin()); pItr != gridParams.end(); ++pItr)
      _tree.Branch(pItr->first, &pItr->second, pItr->first + "/F");

    if(_savePF){
      _tree.Branch("pf.size", &pf_size, "pf.size/i");
      _tree.Branch("pf.charge", pf_charge, "pf_charge[pf.size]/S");
      _tree.Branch("pf.isPU", pf_isPU, "pf_isPU[pf.size]/O");
      _tree.Branch("pf.pdgId", pf_pdgId, "pf_pdgId[pf.size]/S");
      _tree.Branch("pf.vx", pf_vx, "pf_vx[pf.size]/F");
      _tree.Branch("pf.vy", pf_vy, "pf_vy[pf.size]/F");
      _tree.Branch("pf.vz", pf_vz, "pf_vz[pf.size]/F");
      _tree.Branch("pf.pt", pf_pt, "pf_pt[pf.size]/F");
      _tree.Branch("pf.eta", pf_eta, "pf_eta[pf.size]/F");
      _tree.Branch("pf.phi", pf_phi, "pf_phi[pf.size]/F");
      _tree.Branch("pf.mass", pf_mass, "pf_mass[pf.size]/F");
      _tree.Branch("pf.px", pf_px, "pf_px[pf.size]/F");
      _tree.Branch("pf.py", pf_py, "pf_py[pf.size]/F");
      _tree.Branch("pf.pz", pf_pz, "pf_pz[pf.size]/F");
      _tree.Branch("pf.energy", pf_energy, "pf_energy[pf.size]/F");
    }

    if(!_isRealData){
      _tree.Branch("genMet", &genMet, "genMet/F");
      _tree.Branch("genMetPhi", &genMetPhi, "genMetPhi/F");

      _tree.Branch("gen.size", &gen_size, "gen.size/i");
      _tree.Branch("gen.status", gen_status, "status[gen.size]/s");
      _tree.Branch("gen.charge", gen_charge, "charge[gen.size]/S");
      _tree.Branch("gen.motherIndex", gen_motherIndex, "motherIndex[gen.size]/S");
      _tree.Branch("gen.pdgId", gen_pdgId, "pdgId[gen.size]/I");
      _tree.Branch("gen.vx", gen_vx, "vx[gen.size]/F");
      _tree.Branch("gen.vy", gen_vy, "vy[gen.size]/F");
      _tree.Branch("gen.vz", gen_vz, "vz[gen.size]/F");
      _tree.Branch("gen.pt", gen_pt, "pt[gen.size]/F");
      _tree.Branch("gen.eta", gen_eta, "eta[gen.size]/F");
      _tree.Branch("gen.phi", gen_phi, "phi[gen.size]/F");
      _tree.Branch("gen.mass", gen_mass, "mass[gen.size]/F");
      _tree.Branch("gen.px", gen_px, "px[gen.size]/F");
      _tree.Branch("gen.py", gen_py, "py[gen.size]/F");
      _tree.Branch("gen.pz", gen_pz, "pz[gen.size]/F");
      _tree.Branch("gen.energy", gen_energy, "energy[gen.size]/F");
    }
  }

  void
  SimpleEventProducer::AdditionalObjVars::bookBranches(TTree& _tree, bool _isRealData, bool _bookPhoton/* = true*/, bool _bookElectron/* = true*/, bool _bookMuon/* = true*/, bool _bookJet/* = true*/)
  {
    if(_bookPhoton){
      if(!_isRealData){
        _tree.Branch("photon.dRGen", ph_dRGen, "dRGen[photon.size]/F");
        _tree.Branch("photon.genIso", ph_genIso, "genIso[photon.size]/F");
        _tree.Branch("photon.nearestGen", ph_nearestGen, "nearestGen[photon.size]/I");
      }
      _tree.Branch("photon.dRJet", ph_dRJet, "dRJet[photon.size]/F");
      _tree.Branch("photon.dRNextJet", ph_dRNextJet, "dRNextJet[photon.size]/F");
      _tree.Branch("photon.dRPF", ph_dRPF, "dRPF[photon.size]/F");
      _tree.Branch("photon.nearestPF", ph_nearestPF, "nearestPF[photon.size]/S");
      _tree.Branch("photon.pfIsPU", ph_pfIsPU, "pfIsPU[photon.size]/O");
    }
    if(_bookElectron){
      if(!_isRealData){
        _tree.Branch("electron.dRGen", el_dRGen, "dRGen[electron.size]/F");
        _tree.Branch("electron.genIso", el_genIso, "genIso[electron.size]/F");
        _tree.Branch("electron.nearestGen", el_nearestGen, "nearestGen[electron.size]/I");
      }
      _tree.Branch("electron.dRJet", el_dRJet, "dRJet[electron.size]/F");
      _tree.Branch("electron.dRNextJet", el_dRNextJet, "dRNextJet[electron.size]/F");
      _tree.Branch("electron.dRPhoton", el_dRPhoton, "dRPhoton[electron.size]/F");
      _tree.Branch("electron.dRNextPhoton", el_dRNextPhoton, "dRNextPhoton[electron.size]/F");
      _tree.Branch("electron.dRPF", el_dRPF, "dRPF[electron.size]/F");
      _tree.Branch("electron.nearestPF", el_nearestPF, "nearestPF[electron.size]/S");
      _tree.Branch("electron.pfIsPU", el_pfIsPU, "pfIsPU[electron.size]/O");
    }
    if(_bookMuon){
      if(!_isRealData){
        _tree.Branch("muon.dRGen", mu_dRGen, "dRGen[muon.size]/F");
        _tree.Branch("muon.genIso", mu_genIso, "genIso[muon.size]/F");
        _tree.Branch("muon.nearestGen", mu_nearestGen, "nearestGen[muon.size]/I");
      }
      _tree.Branch("muon.dRJet", mu_dRJet, "dRJet[muon.size]/F");
      _tree.Branch("muon.dRNextJet", mu_dRNextJet, "dRNextJet[muon.size]/F");
      _tree.Branch("muon.dRPhoton", mu_dRPhoton, "dRPhoton[muon.size]/F");
      _tree.Branch("muon.dRNextPhoton", mu_dRNextPhoton, "dRNextPhoton[muon.size]/F");
      _tree.Branch("muon.dRPF", mu_dRPF, "dRPF[muon.size]/F");
      _tree.Branch("muon.nearestPF", mu_nearestPF, "nearestPF[muon.size]/S");
      _tree.Branch("muon.pfIsPU", mu_pfIsPU, "pfIsPU[muon.size]/O");
    }
    if(_bookJet){
      if(!_isRealData){
        _tree.Branch("jet.dRGen", jt_dRGen, "dRGen[jet.size]/F");
        _tree.Branch("jet.genSumPt", jt_genSumPt, "genSumPt[jet.size]/F");
        _tree.Branch("jet.nearestGen", jt_nearestGen, "nearestGen[jet.size]/I");
      }
    }
  }

}
