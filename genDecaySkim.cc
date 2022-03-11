#include "SusyEvent.h"
#include "Toolset/GenTreeViewer/test/GenDecayFilterRA3.cc"

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TString.h"
#include "TList.h"
#include "TROOT.h"
#include "TSystem.h"

#include <iostream>

class GenDecaySkim {
public:
  GenDecaySkim();
  ~GenDecaySkim();
  bool initialize(char const*, char const*);
  bool run(char const*);
  bool finalize();
  void setThrow(bool _val) { throw_ = _val; }
private:
  TTree* output_;
  GenDecayFilterRA3* filter_;
  bool throw_;
};

GenDecaySkim::GenDecaySkim() :
  output_(0),
  filter_(0),
  throw_(false)
{
  gROOT->LoadMacro(TString(gSystem->Getenv("CMSSW_BASE")) + "/src/Toolset/GenTreeViewer/test/GenDecayFilterRA3.cc+");
}

GenDecaySkim::~GenDecaySkim()
{
  if(output_) delete output_->GetCurrentFile();
  delete filter_;
}

bool
GenDecaySkim::initialize(char const* _outputDir, char const* _decayChains)
{
  TString outputName(_outputDir);
  if(!outputName.Contains(".root")) outputName += "/susyEvents.root";

  TFile* outputFile(TFile::Open(outputName, "recreate"));
  if(!outputFile || outputFile->IsZombie()){
    std::cout << "Output " << outputName << " not opened" << std::endl;
    delete outputFile;
    if(throw_) throw std::runtime_error("IOError");
    else return false;
  }

  outputFile->cd();
  output_ = new TTree("susyTree", "SUSY Events");
  output_->SetAutoSave(10000000);

  filter_ = new GenDecayFilterRA3(_decayChains);

  return true;
}

bool
GenDecaySkim::run(char const* _sourcePath)
{
  TChain input("susyTree");
  input.Add(_sourcePath);

  if(input.GetEntries() == 0) return false;

  input.SetBranchStatus("*", 0);
  input.SetBranchStatus("genParticles*", 1);

  susy::Event event;
  event.setInput(input);

  TChain fullInput("susyTree");
  fullInput.Add(_sourcePath);

  susy::Event fullEvent;
  fullEvent.setInput(fullInput);
  fullEvent.addOutput(*output_);

  long iEntry(0);
  while(event.getEntry(iEntry++) > 0){
    if(iEntry % 10000 == 1) std::cout << iEntry << std::endl;

    if(!filter_->pass(event)) continue;

    fullEvent.getEntry(iEntry - 1);
    output_->Fill();
  }

  event.releaseTrees();
  fullEvent.releaseTrees();

  return true;
}

bool
GenDecaySkim::finalize()
{
  if(output_){
    TFile* outputFile(output_->GetCurrentFile());
    outputFile->cd();
    outputFile->Write();
    delete outputFile;
  }
  output_ = 0;

  delete filter_;
  filter_ = 0;

  return true;
}

void
genDecaySkim(TString const& _decayChains, TString const& _sourcePath, TString const& _outputName)
{
  GenDecaySkim skim;
  if(!skim.initialize(_outputName, _decayChains)) return;

  TObjArray* sourcePaths(_sourcePath.Tokenize(" "));
  for(int iS(0); iS != sourcePaths->GetEntries(); ++iS)
    if(!skim.run(sourcePaths->At(iS)->GetName())) return;
  delete sourcePaths;

  skim.finalize();
}

void
genDecaySkim(TString const& _decayChains, TObjArray* _urls, TObjArray* _outputNames)
{
  GenDecaySkim skim;
  if(!skim.initialize(_outputNames->At(0)->GetName(), _decayChains)) return;
  for(int iS(0); iS != _urls->GetEntries(); ++iS)
    if(!skim.run(_urls->At(iS)->GetName())) return;
  skim.finalize();
}
