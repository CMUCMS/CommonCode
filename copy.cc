#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <iostream>
#include <map>
#include <set>

#include "SusyEvent.h"
#include "SusyTriggerEvent.h"

/* 
   Copy susyEvents & susyTriggers over selected lumi ranges
*/

void
copy(TString const& _sourceName, TString const& _outputName, TString const& _lumisToRun, bool _veto = false, long _start = 0, long _nEvents = -1)
{
  if(!_outputName.Contains("susyEvents")){
    std::cerr << "Output name must contain susyEvents" << std::endl;
    return;
  }

  std::map<unsigned, std::set<unsigned> > runList;

  TObjArray* lumiRanges(_lumisToRun.Tokenize(","));
  for(int iR(0); iR != lumiRanges->GetEntries(); ++iR){
    TString range(lumiRanges->At(iR)->GetName());
    if(range.Contains("-")){
      TString begin(range(0, range.Index("-")));
      TString end(range(range.Index("-") + 1, range.Length()));
      if(begin(0, begin.Index(":")) != end(0, end.Index(":"))) continue;
      unsigned run(TString(begin(0, begin.Index(":"))).Atoi());
      unsigned lumiFirst(TString(begin(begin.Index(":") + 1, begin.Length())).Atoi());
      unsigned lumiLast(TString(end(end.Index(":") + 1, end.Length())).Atoi());
      for(unsigned lumi(lumiFirst); lumi != lumiLast + 1; ++lumi)
        runList[run].insert(lumi);
    }
    else
      runList[TString(range(0, range.Index(":"))).Atoi()].insert(TString(range(range.Index(":") + 1, range.Length())).Atoi());
  }
  delete lumiRanges;

  for(std::map<unsigned, std::set<unsigned> >::iterator vItr(runList.begin()); vItr != runList.end(); ++vItr){
    for(std::set<unsigned>::iterator lItr(vItr->second.begin()); lItr != vItr->second.end(); ++lItr)
      std::cout << vItr->first << ":" << *lItr << std::endl;
  }

  TChain input("susyTree");
  input.Add(_sourceName);

  susy::Event event;
  susy::TriggerEvent triggerEvent;

  triggerEvent.bindTree(&input, "susyEvents", "susyTriggers");
  event.setInput(input);

  TFile* outputFile(TFile::Open(_outputName, "recreate"));
  TTree* output(new TTree("susyTree", "SUSY Event"));
  output->SetAutoSave(10000000);
  susy::TriggerEvent triggerOutput;
  TString trigOutputName(_outputName);
  triggerOutput.bookTrees(trigOutputName.ReplaceAll("susyEvents", "susyTriggers"));

  event.addOutput(*output);

  long iEntry(_start);
  long iEnd(_start + _nEvents);
  while(iEntry != iEnd && event.getEntry(iEntry++) > 0){
    bool match(runList.find(event.runNumber) != runList.end() &&
               runList[event.runNumber].find(event.luminosityBlockNumber) != runList[event.runNumber].end());
    if(_veto && match) continue;
    if(!_veto && !match) continue;

    triggerOutput.copyEvent(triggerEvent);
    output->Fill();
  }

  event.releaseTrees();
  triggerEvent.reset();

  triggerOutput.write();
  triggerOutput.reset();

  outputFile->cd();
  output->Write();
  delete outputFile;
}
