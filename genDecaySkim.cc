#include "SusyEvent.h"
#include "Toolset/GenTreeViewer/test/GenDecayFilterRA3.cc"

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"

#include <iostream>

void
genDecaySkim(TString const& _decayChains, TTree& _input, TString const& _outputName)
{
  gROOT->LoadMacro(TString(gSystem->Getenv("CMSSW_BASE")) + "/src/Toolset/GenTreeViewer/test/GenDecayFilterRA3.cc+");

  TFile* outputFile(TFile::Open(_outputName, "recreate"));
  if(!outputFile || outputFile->IsZombie())
    throw std::runtime_error("IOError");

  TChain fullInput("susyTree");

  if(_input.InheritsFrom(TChain::Class())){
    TChain& chain(static_cast<TChain&>(_input));
    TObjArray* fileElems(chain.GetListOfFiles());
    for(int iE(0); iE != fileElems->GetEntries(); ++iE)
      fullInput.Add(fileElems->At(iE)->GetTitle());
  }
  else
    fullInput.Add(_input.GetCurrentFile()->GetName());

  susy::Event* event(new susy::Event);
  event->setInput(_input);

  susy::Event* fullEvent(new susy::Event);
  fullEvent->setInput(fullInput);

  outputFile->cd();
  TTree* output(new TTree("susyTree", "SUSY Events"));
  output->SetAutoSave(10000000);
  fullEvent->addOutput(*output);

  GenDecayFilterRA3 filter(_decayChains);

  long iEntry(0);
  while(event->getEntry(iEntry++) > 0){
    if(iEntry % 10000 == 1) std::cout << iEntry << std::endl;

    if(!filter.pass(*event)) continue;

    fullEvent->getEntry(iEntry - 1);
    output->Fill();
  }

  delete event;
  delete fullEvent;

  outputFile->cd();
  output->Write();
  delete outputFile;
}

void
genDecaySkim(TString const& _decayChains, TObjArray* _urls, TObjArray* _outputNames)
{
  TChain input("susyTree");
  input.AddFileInfoList(_urls);
  input.SetBranchStatus("*", 0);
  input.SetBranchStatus("genParticles*", 1);

  TString outputName(_outputNames->At(0)->GetName());

  genDecaySkim(_decayChains, input, outputName);
}
