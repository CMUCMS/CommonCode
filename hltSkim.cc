#include "SusyEvent.h"
#include "SusyTriggerEvent.h"

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TSystem.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <cstring>

class HLTSkimmer {
public:
  HLTSkimmer();
  ~HLTSkimmer();
  bool initialize(char const*, char const*);
  void addInput(char const*, char const* = "");
  bool run();
  void clearInput();
  bool finalize();

private:
  TChain* input_;
  TChain* fullInput_;
  TChain* triggerInput_;
  TTree* output_;
  susy::TriggerEvent* triggerOutput_;

  std::vector<TString> paths_;
  bool fillTriggerEvent_;
  TString triggerOutputName_;
};

HLTSkimmer::HLTSkimmer() :
  input_(0),
  fullInput_(0),
  triggerInput_(0),
  output_(0),
  triggerOutput_(0),
  paths_(),
  fillTriggerEvent_(false),
  triggerOutputName_("")
{
}

HLTSkimmer::~HLTSkimmer()
{
  delete input_;
  delete fullInput_;
  delete triggerInput_;
  delete output_;
  delete triggerOutput_;
}

bool
HLTSkimmer::initialize(char const* _outputDir, char const* _pathList)
{
  TString outputName(_outputDir);
  if(!outputName.Contains(".root")) outputName += "/susyEvents.root";

  TFile* outputFile(TFile::Open(outputName, "recreate"));
  if(!outputFile) return false;
  output_ = new TTree("susyTree", "SUSY Event");
  output_->SetAutoSave(10000000);

  std::cout << "Output name " << outputName;

  if(outputName.Contains("susyEvents")){
    triggerOutputName_ = outputName;
    triggerOutputName_.ReplaceAll("susyEvents", "susyTriggers");
    triggerOutput_ = new susy::TriggerEvent;
    triggerOutput_->bookTrees(triggerOutputName_);

    std::cout << " " << triggerOutputName_;
  }

  std::cout << std::endl;

  TString list(_pathList);
  TObjArray* paths(list.Tokenize(" "));
  for(int iP(0); iP != paths->GetEntries(); ++iP)
    paths_.push_back(TString(paths->At(iP)->GetName()) + "_v*");
  delete paths;

  return paths_.size() != 0;
}

void
HLTSkimmer::addInput(char const* _source, char const* _triggerSource/* = ""*/)
{
  if(!input_) input_ = new TChain("susyTree");
  if(!fullInput_) fullInput_ = new TChain("susyTree");

  input_->Add(_source);
  fullInput_->Add(_source);

  if(triggerOutput_ && _triggerSource && _triggerSource[0] != '\0'){
    if(!triggerInput_) triggerInput_ = new TChain("triggerEventTree");
    triggerInput_->Add(_triggerSource);
    fillTriggerEvent_ = true;
  }
}

bool
HLTSkimmer::run()
{
  input_->SetBranchStatus("*", 0);
  input_->SetBranchStatus("hlt*", 1);

  susy::Event event;
  event.setInput(*input_);

  susy::Event fullEvent;
  fullEvent.setInput(*fullInput_);
  fullEvent.addOutput(*output_);

  susy::TriggerEvent triggerEvent;
  if(fillTriggerEvent_){
    if(!triggerEvent.bindTree(fullInput_, triggerInput_))
      fillTriggerEvent_ = false;
  }

  unsigned nP(paths_.size());

  long iEntry(0);
  while(event.getEntry(iEntry++) > 0){
    unsigned iP(0);
    for(; iP != nP; ++iP)
      if(event.hltMap.pass(paths_[iP])) break;
    if(iP == nP) continue;

    fullEvent.getEntry(iEntry - 1);
    if(fillTriggerEvent_) triggerOutput_->copyEvent(triggerEvent);
    output_->Fill();
  }

  fullEvent.releaseTrees();

  return true;
}

void
HLTSkimmer::clearInput()
{
  delete input_;
  input_ = 0;
  delete fullInput_;
  fullInput_ = 0;
  delete triggerInput_;
  triggerInput_ = 0;
}

bool
HLTSkimmer::finalize()
{
  TFile* file(output_->GetCurrentFile());
  file->cd();
  output_->Write();
  delete file;
  output_ = 0;

  if(fillTriggerEvent_) triggerOutput_->write();
  delete triggerOutput_;
  triggerOutput_ = 0;
  if(!fillTriggerEvent_ && triggerOutputName_.Length() != 0)
    gSystem->Unlink(triggerOutputName_);

  fillTriggerEvent_ = false;

  return true;
}

void
hltSkim(TString const& _source, TString const& _triggerSource, TString const& _output, TString const& _pathList)
{
  HLTSkimmer skim;
  if(!skim.initialize(_output, _pathList)) return;
  skim.addInput(_source, _triggerSource);
  if(!skim.run()) return;
  skim.clearInput();
  skim.finalize();
}
