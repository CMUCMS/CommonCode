#include "SimpleEventProducer.h"
#include "Utilities.h"
#include "HLT.h"

#include <iostream>
#include <stdexcept>
#include <cstring>

#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TSystem.h"

class SimpleTreeProducer {
public:
  SimpleTreeProducer();
  ~SimpleTreeProducer();
  bool initialize(char const*, char const* = "", bool = true, bool = true, bool = false, bool = false);
  void addInput(char const*, char const* = "");
  bool run();
  void clearInput();
  bool finalize();
  void setThrow(bool _val) { throw_ = _val; }
protected:
  bool initInput_(TTree&, susy::Event&, TTree&, susy::TriggerEvent&);
  bool processError_(std::exception&, susy::Event const&, TTree const&);

  TChain input_;
  TChain triggerInput_;
  
  TFile* outputFile_;
  TTree* evtTree_;
  TTree* selectedObjTree_;
  TTree* allObjTree_;
  susy::SimpleEventProducer eventProducer_;
  bool fillTriggerEvent_;
  bool throw_;
};

SimpleTreeProducer::SimpleTreeProducer() :
  input_("susyTree"),
  triggerInput_("triggerEventTree"),
  outputFile_(0),
  evtTree_(0),
  selectedObjTree_(0),
  allObjTree_(0),
  eventProducer_(),
  fillTriggerEvent_(false),
  throw_(false)
{
}

SimpleTreeProducer::~SimpleTreeProducer()
{
  delete outputFile_;
}

bool
SimpleTreeProducer::initialize(char const* _outputDir, char const* _puScenario/* = ""*/, bool _fillSelected/* = true*/, bool _fillAll/* = true*/, bool _fillPF/* = false*/, bool _fillTriggerEvent/* = false*/)
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

  TH1* puWeights(0);

  TString puScenario(_puScenario);
  if(puScenario.Length() != 0){
    std::cout << "Running in MC mode with PU scenario " << puScenario << std::endl;

    /* DEFINE LIST OF MC PARAMETERS TO INCLUDE */
    eventProducer_.addGridParam("ptHat");

    if(puScenario != "None"){
      TFile* puWeightSource(TFile::Open(TString(gSystem->DirName(__FILE__)) + "/puReweighting_2012.root"));
      if(!puWeightSource || puWeightSource->IsZombie() || !(puWeights = dynamic_cast<TH1*>(puWeightSource->Get("weight" + puScenario)))){
        std::cerr << "PU weights for scenario " << puScenario << " not found" << std::endl;
        if(throw_) throw std::invalid_argument("PU weights");
        else return false;
      }
    }
  }
  else{
    std::cout << "Running in real-data mode as no PU scenario was specified." << std::endl;
    std::cout << "If this is a MC sample, check the dataset name and rerun the function with a valid scenario name (e.g. S10, PU50)." << std::endl;
  }

  /* INITIALIZE EVENT PRODUCER */
  eventProducer_.initialize(evtTree_, selectedObjTree_, allObjTree_, puWeights);

  return true;
}

void
SimpleTreeProducer::addInput(char const* _eventSource, char const* _triggerSource/* = ""*/)
{
  input_.Add(_eventSource);
  if(fillTriggerEvent_ && std::strlen(_triggerSource) != 0){
    triggerInput_.Add(_triggerSource);
  }
}

bool
SimpleTreeProducer::run()
{
  susy::Event event;
  susy::TriggerEvent triggerEvent;

  if(!initInput_(input_, event, triggerInput_, triggerEvent)){
    std::cerr << "Input incompatible" << std::endl;
    if(throw_) throw std::runtime_error("");
    else return false;
  }

  /* EVENT LOOP */

  long iEntry(0);
  int nRead(0);

  while((nRead = event.getEntry(iEntry++)) != 0){
    try{
      if(nRead < 0)
        throw susy::Exception(susy::Exception::kIOError, "Corrupt input");

      if(iEntry % 10000 == 1) std::cout << "Processing event " << iEntry - 1 << "..." << std::endl;

      if(fillTriggerEvent_) eventProducer_.extractTriggerObjects(triggerEvent);
      eventProducer_.produce(event);

      if(evtTree_->Fill() < 0)
        throw susy::Exception(susy::Exception::kIOError, "eventVars");
      if(selectedObjTree_ && selectedObjTree_->Fill() < 0)
        throw susy::Exception(susy::Exception::kIOError, "selectedObjects");
      if(allObjTree_ && allObjTree_->Fill() < 0)
        throw susy::Exception(susy::Exception::kIOError, "allObjects");
    }
    catch(std::exception& e){
      if(processError_(e, event, input_)) continue;
      if(throw_) throw;
      else return false;
    }
  }

  std::cout << "Processed " << iEntry - 1 << " Events." << std::endl;

  return true;
}

void
SimpleTreeProducer::clearInput()
{
  input_.Reset();
  triggerInput_.Reset();
}

bool
SimpleTreeProducer::finalize()
{
  /* CLEANUP & FINALZE */
  outputFile_->cd();
  outputFile_->Write();
  delete outputFile_;

  outputFile_ = 0;
  evtTree_ = 0;
  selectedObjTree_ = 0;
  allObjTree_ = 0;

  return true;
}

bool
SimpleTreeProducer::initInput_(TTree& _input, susy::Event& _event, TTree& _triggerInput, susy::TriggerEvent& _triggerEvent)
{
  /* DISABLE UNUSED INPUT BRANCHES TO SPEED UP THE PROCESSING */

  _input.SetBranchStatus("*", 0);
  _input.SetBranchStatus("runNumber", 1);
  _input.SetBranchStatus("luminosityBlockNumber", 1);
  _input.SetBranchStatus("eventNumber", 1);
  _input.SetBranchStatus("metFilter*", 1);
  _input.SetBranchStatus("hlt*", 1);
  _input.SetBranchStatus("pfParticles*", 1);
  _input.SetBranchStatus("met_pfType1CorrectedMet*", 1);
  _input.SetBranchStatus("beamSpot*", 1);
  if(_input.GetBranch("pu")) _input.SetBranchStatus("pu*", 1);
  if(_input.GetBranch("genParticles")) _input.SetBranchStatus("genParticles*", 1);
  if(_input.GetBranch("met_genMetTrue.")) _input.SetBranchStatus("met_genMetTrue*", 1);
  if(_input.GetBranch("gridParams_ptHat")) _input.SetBranchStatus("gridParams*", 1);
  susy::ObjectTree::setBranchStatus(_input); // set status = 1 for photon-, electron-, muon-, jet-, and vertex-related branches

  if(fillTriggerEvent_ && !_triggerEvent.bindTree(&_input, &_triggerInput)) return false;

  _event.setInput(_input);

  return true;
}

bool
SimpleTreeProducer::processError_(std::exception& _ex, susy::Event const& _event, TTree const& _input)
{
  std::cerr << "Exception caught:" << std::endl;
  std::cerr << _ex.what() << std::endl;
  std::cerr << "Run " << _event.runNumber << ", Lumi " << _event.luminosityBlockNumber << ", Event " << _event.eventNumber << " in " << std::endl;
  std::cerr << "File " << _input.GetCurrentFile()->GetName() << " Entry " << _input.GetReadEntry() << std::endl;

  susy::Exception* susyExcept(dynamic_cast<susy::Exception*>(&_ex));

  if(susyExcept){
    std::cerr << "This was an exception of category " << susyExcept->categoryName() << std::endl;

    switch(susyExcept->category){
    case susy::Exception::kEventAnomaly:
    case susy::Exception::kObjectAnomaly:
      std::cerr << "Skipping event.." << std::endl;
      return true;
    case susy::Exception::kIOError:
    case susy::Exception::kFormatError:
    default:
      break;
    }
  }

  return false;
}

void
produceSimpleTree(TString const& _sourceName, TString const& _outputName, char const* _puScenario = "", bool _fillSelected = true, bool _fillAll = true, bool _fillPF = false, bool _fillTriggerEvent = false)
{
  SimpleTreeProducer producer;
  producer.setThrow(true);
  producer.initialize(_outputName, _puScenario, _fillSelected, _fillAll, _fillPF, _fillTriggerEvent);

  TObjArray* sourcePaths(_sourceName.Tokenize(","));
  if(sourcePaths->GetEntries() > 1)
    producer.addInput(sourcePaths->At(0)->GetName(), sourcePaths->At(1)->GetName());
  else
    producer.addInput(sourcePaths->At(0)->GetName());
  delete sourcePaths;

  producer.run();
  producer.clearInput();
  producer.finalize();
}

void
produceSimpleTree(char const* _puScenario, bool _fillSelected, bool _fillAll, bool _fillPF, bool _fillTriggerEvent, TObjArray* _urls, TObjArray* _outputName)
{
  SimpleTreeProducer producer;
  producer.setThrow(true);
  producer.initialize(_outputName->At(0)->GetName(), _puScenario, _fillSelected, _fillAll, _fillPF, _fillTriggerEvent);
  for(int iS(0); iS != _urls->GetEntries(); ++iS){
    TObjArray* sourcePaths(static_cast<TObjString*>(_urls->At(iS))->GetString().Tokenize(","));
    if(sourcePaths->GetEntries() > 1)
      producer.addInput(sourcePaths->At(0)->GetName(), sourcePaths->At(1)->GetName());
    else
      producer.addInput(sourcePaths->At(0)->GetName());
  }
  producer.run();
  producer.clearInput();
  producer.finalize();
}

void
produceSimpleTree(char const* _puScenario, TObjArray* _urls, TObjArray* _outputName)
{
  produceSimpleTree(_puScenario, true, true, true, false, _urls, _outputName);
}

void
produceSimpleTree(char const* _puScenario, bool _fillSelected, bool _fillAll, bool _fillPF, bool _fillTriggerEvent, TString const& _dataset, TObjArray* _inputLines, TObjArray* _outputDir)
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

  produceSimpleTree(_puScenario, _fillSelected, _fillAll, _fillPF, _fillTriggerEvent, &urls, &outputName);
}

void
produceSimpleTree(char const* _puScenario, TString const& _dataset, TObjArray* _inputLines, TObjArray* _outputDir)
{
  produceSimpleTree(_puScenario, true, true, true, false, _dataset, _inputLines, _outputDir);
}
