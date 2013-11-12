#include "SimpleEventProducer.h"
#include "Utilities.h"
#include "HLT.h"

#include <iostream>
#include <stdexcept>

#include "TFile.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"

class SimpleTreeProducer {
public:
  SimpleTreeProducer();
  ~SimpleTreeProducer() {}
  bool initialize(char const*, bool = true, bool = true, bool = false, bool = false);
  bool run(char const*);
  bool finalize();
  void setThrow(bool _val) { throw_ = _val; }
private:
  TFile* outputFile_;
  TTree* evtTree_;
  TTree* selectedObjTree_;
  TTree* allObjTree_;
  susy::SimpleEventProducer eventProducer_;
  bool producerInitialized_;
  bool fillTriggerEvent_;
  bool throw_;
};

SimpleTreeProducer::SimpleTreeProducer() :
  outputFile_(0),
  evtTree_(0),
  selectedObjTree_(0),
  allObjTree_(0),
  eventProducer_(),
  producerInitialized_(false),
  fillTriggerEvent_(false),
  throw_(false)
{
}

bool
SimpleTreeProducer::initialize(char const* _outputDir, bool _fillSelected, bool _fillAll, bool _fillPF, bool _fillTriggerEvent)
{
  /* DEFINE OUTPUT TREES */

  TString outputName(_outputDir);
  if(!outputName.Contains(".root")) outputName += "/simpleTree.root";

  outputFile_ = TFile::Open(outputName, "recreate");
  if(!outputFile_ || outputFile_->IsZombie()){
    std::cerr << "Output " << outputName << " not opened" << std::endl;
    delete outputFile_;
    if(throw_) throw std::runtime_error("IOError");
    else return false;
  }

  outputFile_->cd();
  evtTree_ = new TTree("eventVars", "Event variables");
  evtTree_->SetAutoSave(10000000);
  if(_fillSelected){
    selectedObjTree_ = new TTree("selectedObjects", "Selected objects");
    selectedObjTree_->SetAutoSave(10000000);
  }
  if(_fillAll){
    allObjTree_ = new TTree("allObjects", "All objects");
    allObjTree_->SetAutoSave(10000000);
  }

  fillTriggerEvent_ = _fillTriggerEvent;

  /* SETUP EVENT PRODUCER */

  eventProducer_.setSavePF(_fillPF);

  for(unsigned iT(0); iT != nHLTPaths; ++iT)
    eventProducer_.addHLTPath(hltPaths[iT]);

  if(fillTriggerEvent_){
    for(unsigned iF(0); iF != nHLTEventFilters; ++iF)
      eventProducer_.addHLTEventFilter(hltEventFilters[iF]);
    for(unsigned iF(0); iF != nHLTPhotonFilters; ++iF)
      eventProducer_.addHLTPhotonFilter(hltPhotonFilters[iF]);
    for(unsigned iF(0); iF != nHLTElectronFilters; ++iF)
      eventProducer_.addHLTElectronFilter(hltElectronFilters[iF]);
    for(unsigned iF(0); iF != nHLTMuonFilters; ++iF)
      eventProducer_.addHLTMuonFilter(hltMuonFilters[iF]);
  }

  eventProducer_.setPhotonId(susy::PhLoose12);
  eventProducer_.setElectronId(susy::ElMedium12);
  eventProducer_.setMuonId(susy::MuTight12);
  eventProducer_.setJetId(susy::JtLoose);

  return true;
}

bool
SimpleTreeProducer::run(char const* _sourcePath)
{
  TObjArray* sourcePaths(TString(_sourcePath).Tokenize(","));
  if(sourcePaths->GetEntries() == 0){
    delete sourcePaths;
    std::cerr << "No source name provided" << std::endl;
    if(throw_) throw std::invalid_argument("Input");
    else return false;
  }
  TString susyTreePath(sourcePaths->At(0)->GetName());
  TString susyTriggerPath;
  if(fillTriggerEvent_ && sourcePaths->GetEntries() == 2) susyTriggerPath = sourcePaths->At(1)->GetName();
  delete sourcePaths;

  TChain input("susyTree");
  input.Add(susyTreePath);

  if(input.GetEntries() == 0) return false;

  TChain triggerInput("triggerEventTree");
  if(susyTriggerPath.Length() != 0) triggerInput.Add(susyTriggerPath);

  TBranch* bIsRealData(input.GetBranch("isRealData"));
  TLeaf* lIsRealData(0);
  if(!bIsRealData || !(lIsRealData = bIsRealData->GetLeaf("isRealData"))){
    std::cerr << "Input does not contain the branch isRealData" << std::endl;
    if(throw_) throw std::invalid_argument("Input");
    else return false;
  }
  bIsRealData->GetEntry(0);
  bool isRealData(lIsRealData->GetValue(0) != 0);

  /* DISABLE UNUSED INPUT BRANCHES TO SPEED UP THE PROCESSING */

  input.SetBranchStatus("*", 0);
  input.SetBranchStatus("runNumber", 1);
  input.SetBranchStatus("luminosityBlockNumber", 1);
  input.SetBranchStatus("eventNumber", 1);
  input.SetBranchStatus("metFilter*", 1);
  input.SetBranchStatus("hlt*", 1);
  input.SetBranchStatus("pfParticles*", 1);
  input.SetBranchStatus("met_pfType01CorrectedMet*", 1);
  if(!isRealData){
    input.SetBranchStatus("genParticles*", 1);
    input.SetBranchStatus("met_genMetTrue*", 1);
    input.SetBranchStatus("gridParams*", 1);
  }
  input.SetBranchStatus("beamSpot*", 1);
  susy::ObjectTree::setBranchStatus(input); // set status = 1 for photon-, electron-, muon-, jet-, and vertex-related branches

  /* SET INPUT BRANCH ADDRESS TO EVENT OBJECT */

  susy::Event event;
  susy::TriggerEvent triggerEvent;

  if(fillTriggerEvent_){
    if(susyTriggerPath.Length() != 0)
      triggerEvent.bindTree(&input, &triggerInput);
    else
      triggerEvent.bindTree(&input, "susyEvents", "susyTriggers");
  }

  event.setInput(input);

  if(event.getEntry(0) <= 0){
    std::cerr << "Input is empty or corrupted" << std::endl;
    delete outputFile_;
    if(throw_) throw std::runtime_error("IOError");
    else return false;
  }

  /* ONE-TIME CONFIGURATION OF THE PRODUCER */
  if(!producerInitialized_){
    /* DEFINE LIST OF MC PARAMETERS TO INCLUDE */

    if(!isRealData){
      TObjArray* leaves(input.GetListOfLeaves());
      for(int iL(0); iL != leaves->GetEntries(); ++iL){
        TString name(leaves->At(iL)->GetName());
        if(name.Index("gridParams_") == 0) eventProducer_.addGridParam(name);
      }
    }

    /* INITIALIZE EVENT PRODUCER */
    eventProducer_.initialize(evtTree_, selectedObjTree_, allObjTree_, isRealData);
    producerInitialized_ = true;
  }

  /* EVENT LOOP */

  bool failed(false);
  long iEntry(0);
  int nRead(0);

  while((nRead = event.getEntry(iEntry++)) != 0){
    if(nRead < 0){
      std::cerr << "Input corrupted at entry " << iEntry - 1 << std::endl;
      if(throw_) throw std::runtime_error("IOError");
      else return false;
    }

    if(iEntry % 10000 == 1) std::cout << "Processing event " << iEntry - 1 << "..." << std::endl;

    try{
      if(fillTriggerEvent_) eventProducer_.extractTriggerObjects(triggerEvent);
      eventProducer_.produce(event);
    }
    catch(std::exception& e){
      std::cerr << "Exception caught:" << std::endl;
      std::cerr << e.what() << std::endl;
      std::cerr << "Run " << event.runNumber << ", Lumi " << event.luminosityBlockNumber << ", Event " << event.eventNumber << " in " << std::endl;
      std::cerr << input.GetCurrentFile()->GetName() << std::endl;

      susy::Exception* susyExcept(dynamic_cast<susy::Exception*>(&e));

      if(susyExcept){
        std::cerr << "This was an exception of category " << susyExcept->categoryName() << std::endl;

        switch(susyExcept->category){
        case susy::Exception::kEventAnomaly:
          std::cerr << "Skipping event.." << std::endl;
          continue;
        case susy::Exception::kObjectAnomaly:
        case susy::Exception::kIOError:
        case susy::Exception::kFormatError:
        default:
          break;
        }
      }

      failed = true;
      break;
    }

    evtTree_->Fill();
    if(selectedObjTree_) selectedObjTree_->Fill();
    if(allObjTree_) allObjTree_->Fill();
  }

  triggerEvent.reset();
  event.releaseTrees();

  std::cout << "Processed " << iEntry - 1 << " Events." << std::endl;

  if(failed){
    if(throw_) throw std::runtime_error("produceSimpleTree");
    else return false;
  }

  return true;
}

bool
SimpleTreeProducer::finalize()
{
  /* CLEANUP & FINALZE */
  outputFile_->cd();
  outputFile_->Write();
  delete outputFile_;

  return true;
}

void
produceSimpleTree(TString const& _sourceName, TString const& _outputName, bool _fillSelected = true, bool _fillAll = true, bool _fillPF = false, bool _fillTriggerEvent = false)
{
  SimpleTreeProducer producer;
  producer.setThrow(true);
  producer.initialize(_outputName, _fillSelected, _fillAll, _fillPF, _fillTriggerEvent) &&
    producer.run(_sourceName) &&
    producer.finalize();
}

void
produceSimpleTree(bool _fillSelected, bool _fillAll, bool _fillPF, bool _fillTriggerEvent, TObjArray* _urls, TObjArray* _outputName)
{
  SimpleTreeProducer producer;
  producer.setThrow(true);
  if(!producer.initialize(_outputName->At(0)->GetName(), _fillSelected, _fillAll, _fillPF, _fillTriggerEvent)) return;
  for(int iS(0); iS != _urls->GetEntries(); ++iS)
    if(!producer.run(_urls->At(iS)->GetName())) return;
  producer.finalize();
}

void
produceSimpleTree(TObjArray* _urls, TObjArray* _outputName)
{
  produceSimpleTree(true, true, true, false, _urls, _outputName);
}

void
produceSimpleTree(bool _fillSelected, bool _fillAll, bool _fillPF, bool _fillTriggerEvent, TString const& _dataset, TObjArray* _inputLines, TObjArray* _outputDir)
{
  // For "externalList" mode of dcmu job scheduler
  // The list must have the grid point name in the (N+1)n-th rows (N=files per point, n=0,1,...) and the file names in the folloing N rows
  // example: TChiwg.list
  //   TChiwg_1000_1000
  //   TChiwg_1000_1000_10_4_2zc.root
  //   ...
  // Then the submit command is
  // submitDCMUJobs.py -n 51 -J TChiwg -x "$PWD/TChiwg.list" -s 'output_directory' -z /store/RA3Ntuples/... produceSimpleTree.cc

  TObjArray outputName;
  outputName.SetOwner(true);

  TString gridPoint(_inputLines->At(0)->GetName());

  TString outputFileName(_outputDir->At(0)->GetName());
  outputFileName += "/" + gridPoint + ".root";
  outputName.Add(new TObjString(outputFileName));

  TObjArray urls;
  urls.SetOwner(true);

  for(int iL(1); iL != _inputLines->GetEntries(); ++iL)
    urls.Add(new TObjString(_dataset + "/" + _inputLines->At(iL)->GetName()));

  produceSimpleTree(_fillSelected, _fillAll, _fillPF, _fillTriggerEvent, &urls, &outputName);
}


void
produceSimpleTree(TString const& _dataset, TObjArray* _inputLines, TObjArray* _outputDir)
{
  produceSimpleTree(true, true, true, false, _dataset, _inputLines, _outputDir);
}
