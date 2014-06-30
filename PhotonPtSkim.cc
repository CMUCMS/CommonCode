#include "SusyEvent.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

#include <cmath>
#include <fstream>
#include <iostream>

class PtSkimmer {
public:
  PtSkimmer();
  ~PtSkimmer();
  bool initialize(char const*, double = 25.);
  void addInput(char const*);
  bool run();
  void clearInput();
  bool finalize();

private:
  TChain* input_;
  TChain* fullInput_;
  TTree* output_;
  double threshold_;
};

PtSkimmer::PtSkimmer() :
  input_(0),
  fullInput_(0),
  output_(0),
  threshold_(0.)
{
}

PtSkimmer::~PtSkimmer()
{
  delete input_;
  delete fullInput_;
  delete output_;
}

bool
PtSkimmer::initialize(char const* _outputName, double _pt/* = 25.*/)
{
  TString outputName(_outputName);
  if(!outputName.Contains(".root")) outputName += "/susyEvents.root";

  std::cout << "Output name " << outputName << std::endl;

  TFile* outputFile(TFile::Open(outputName, "recreate"));
  if(!outputFile) return false;
  output_ = new TTree("susyTree", "SUSY Event");
  output_->SetAutoSave(10000000);

  threshold_ = _pt;

  return true;
}

void
PtSkimmer::addInput(char const* _source)
{
  if(!input_) input_ = new TChain("susyTree");
  if(!fullInput_) fullInput_ = new TChain("susyTree");

  input_->Add(_source);
  fullInput_->Add(_source);
}

bool
PtSkimmer::run()
{
  input_->SetBranchStatus("*", 0);
  input_->SetBranchStatus("photons_photons*", 1);
  if(input_->GetBranch("genParticles")) input_->SetBranchStatus("genParticles*", 1);

  susy::Event event;
  event.setInput(*input_);

  susy::Event fullEvent;
  fullEvent.setInput(*fullInput_);
  fullEvent.addOutput(*output_);

  long iEntry(0);
  while(event.getEntry(iEntry++) > 0){
    susy::PhotonCollection const& photons(event.photons["photons"]);
    unsigned nP(photons.size());
    susy::ParticleCollection const& genParticles(event.genParticles);
    unsigned nG(genParticles.size());

    unsigned iP(0);
    for(; iP != nP; ++iP){
      susy::Photon const& photon(photons[iP]);
      if(photon.momentum.Pt() < threshold_) continue;
      if(std::abs(photon.caloPosition.Eta()) > susy::etaGapBegin) continue;

      unsigned iG(0);
      for(; iG != nG; ++iG){
        susy::Particle const& part(genParticles[iG]);
        if(part.status != 1) continue;
        if(std::abs(part.pdgId) != 11) continue;
        if(part.momentum.DeltaR(photon.momentum) > 0.1) continue;

        short idx(part.motherIndex);
        while(idx != -1 && genParticles[idx].pdgId != 23 && std::abs(genParticles[idx].pdgId) != 24) idx = genParticles[idx].motherIndex;
        if(idx != -1) break;
      }
      if(iG == nG) break;
    }

    if(iP == nP) continue;

    fullEvent.getEntry(iEntry - 1);
    output_->Fill();
  }

  fullEvent.releaseTrees();

  return true;
}

void
PtSkimmer::clearInput()
{
  delete input_;
  input_ = 0;
  delete fullInput_;
  fullInput_ = 0;
}

bool
PtSkimmer::finalize()
{
  TFile* file(output_->GetCurrentFile());
  file->cd();
  output_->Write();
  delete file;
  output_ = 0;

  return true;
}

void
PhotonPtSkim(TString const& _source, TString const& _output, double _pt = 25.)
{
  PtSkimmer skim;
  if(!skim.initialize(_output, _pt)) return;
  skim.addInput(_source);
  if(!skim.run()) return;
  skim.clearInput();
  skim.finalize();
}

void
PhotonPtSkim(TString _output, double _pt, TString _list)
{
  std::ifstream list(_list);
  if(!list.is_open()) return;

  PtSkimmer skim;
  if(!skim.initialize(_output, _pt)) return;
  TString source;
  while(list.good()){
    list >> source;
    if(!list.good()) break;
    skim.addInput(source);
  }
  if(!skim.run()) return;
  skim.clearInput();
  skim.finalize();
}
