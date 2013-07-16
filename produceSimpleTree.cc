#include "SimpleEventProducer.h"

#include <iostream>
#include <stdexcept>

#include "TFile.h"
#include "TChain.h"

bool fillSelected(true);
bool fillAll(true);

void
produceSimpleTree(TTree& _input, TString const& _outputName, long _maxEvents = -1)
{
  if(_input.GetEntries() == 0) return;

  /* DEFINE OUTPUT TREES */

  TFile* outputFile(TFile::Open(_outputName, "recreate"));
  if(!outputFile || outputFile->IsZombie()){
    std::cerr << "Output " << _outputName << " not opened" << std::endl;
    delete outputFile;
    throw std::runtime_error("IOError");
  }

  TTree* evtTree(new TTree("eventVars", "Event variables"));
  TTree* selectedObjTree(fillSelected ? new TTree("selectedObjects", "Selected objects") : 0);
  TTree* allObjTree(fillAll ? new TTree("allObjects", "All objects") : 0);

  evtTree->SetAutoSave(10000000);
  selectedObjTree->SetAutoSave(10000000);
  allObjTree->SetAutoSave(10000000);

  /* DISABLE UNUSED INPUT BRANCHES TO SPEED UP THE PROCESSING */

  _input.SetBranchStatus("*", 0);
  _input.SetBranchStatus("runNumber", 1);
  _input.SetBranchStatus("luminosityBlockNumber", 1);
  _input.SetBranchStatus("eventNumber", 1);
  _input.SetBranchStatus("metFilter*", 1);
  _input.SetBranchStatus("hlt*", 1);
  _input.SetBranchStatus("genParticles*", 1);
  _input.SetBranchStatus("pfParticles*", 1);
  _input.SetBranchStatus("met_pfType01CorrectedMet*", 1);
  _input.SetBranchStatus("gridParams*", 1);
  susy::ObjectTree::setBranchStatus(_input); // set status = 1 for photon-, electron-, muon-, jet-, and vertex-related branches

  /* SET INPUT BRANCH ADDRESS TO EVENT OBJECT */

  susy::Event* event(new susy::Event);
  event->setInput(_input);

  if(event->getEntry(0) <= 0){
    std::cerr << "Input is empty or corrupted" << std::endl;
    delete outputFile;
    delete event;
    throw std::runtime_error("IOError");
  }

  /* CREATE EVENT PRODUCER */

  susy::SimpleEventProducer producer;

  /* DEFINE LIST OF HLT PATHS TO INCLUDE */

  TString hltPathList[] = {
    "HLT_IsoMu24_eta2p1",
    "HLT_IsoMu24",
    "HLT_Ele27_WP80",
    "HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50",
    "HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50",
    "HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60",
    "HLT_Photon26_R9Id85_IR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60",
    "HLT_Mu22_Photon22_CaloIdL",
    "HLT_Photon70_CaloIdXL_PFHT400",
    "HLT_Photon70_CaloIdXL_PFHT500",
    "HLT_Photon70_CaloIdXL_PFNoPUHT400",
    "HLT_Photon70_CaloIdXL_PFNoPUHT500",
    "HLT_Photon70_CaloIdXL_PFMET100",
    "HLT_Photon60_CaloIdL_HT300",
    "HLT_Photon60_CaloIdL_MHT70",
    "HLT_Photon135",
    "HLT_Photon150"
  };
  producer.setHLTPaths(std::vector<TString>(hltPathList, hltPathList + sizeof(hltPathList) / sizeof(TString)));

  /* DEFINE LIST OF MC PARAMETERS TO INCLUDE */

  TString gridParamList[] = {
    "pthat"
  };
  producer.setGridParams(std::vector<TString>(gridParamList, gridParamList + sizeof(gridParamList) / sizeof(TString)));

  /* DEFINE GOOD OBJECTS */
  producer.setPhotonId(susy::PhLoose12);
  producer.setElectronId(susy::ElMedium12);
  producer.setMuonId(susy::MuTight12);
  producer.setJetId(susy::JtLoose);

  /* DO NOT SAVE PF PARTICLES (WHICH IS THE DEFAULT) */
  producer.setSavePF(true);

  /* INITIALIZE EVENT PRODUCER */
  producer.initialize(evtTree, selectedObjTree, allObjTree, event->isRealData);

  /* EVENT LOOP */

  bool failed(false);
  long iEntry(0);

  try{

    while(iEntry != _maxEvents){
      int nRead(event->getEntry(iEntry++));
      if(nRead == 0) break;
      else if(nRead < 0){
        std::cerr << "Input corrupted at entry " << iEntry - 1 << std::endl;
        throw std::runtime_error("IOError");
      }

      if(iEntry % 10000 == 1) std::cout << "Processing event " << iEntry - 1 << "..." << std::endl;

      if(!event->passMetFilters()) continue;
   
      producer.produce(*event);

      evtTree->Fill();
      if(fillSelected) selectedObjTree->Fill();
      if(fillAll) allObjTree->Fill();
    }

  }
  catch(std::exception& e){
    std::cerr << "Exception caught:" << std::endl;
    std::cerr << e.what() << std::endl;
    std::cerr << "Run " << event->runNumber << ", Lumi " << event->luminosityBlockNumber << ", Event " << event->eventNumber << " in " << std::endl;
    std::cerr << _input.GetCurrentFile()->GetName() << std::endl;

    failed = true;
  }

  std::cout << "Processed " << iEntry << " Events." << std::endl;

  /* CLEANUP & FINALZE */

  delete event;

  if(!failed){
    outputFile->cd();
    outputFile->Write();
  }
  delete outputFile;

  if(failed)
    throw std::runtime_error("");
}

void
produceSimpleTree(TObjArray* _urls, TObjArray* _outputName)
{
  TChain input("susyTree");
  input.AddFileInfoList(_urls);

  produceSimpleTree(input, _outputName->At(0)->GetName());
}

void
produceSimpleTree(TString const& _dataset, TObjArray* _inputLines, TObjArray* _outputDir)
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

  produceSimpleTree(&urls, &outputName);
}
