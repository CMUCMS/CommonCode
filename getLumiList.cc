#include "SusyEvent.h"

#include "TObjArray.h"
#include "TChain.h"
#include "TSystem.h"

#include <map>
#include <set>
#include <fstream>
#include <iostream>

void
getLumiList(TString const& _sourceName, TString const& _outputName)
{
  TChain input("susyTree");
  input.Add(_sourceName);
  input.SetBranchStatus("*", 0);
  input.SetBranchStatus("runNumber", 1);
  input.SetBranchStatus("luminosityBlockNumber", 1);

  susy::Event* event(new susy::Event);

  event->setInput(input);

  std::map<int, std::set<int> > lumiList;
  unsigned currentRun(0);
  unsigned currentLumi(0);
  std::set<int>* currentSet(0);

  long iEntry(0);
  while(event->getEntry(iEntry++) > 0){
    if(currentRun != event->runNumber){
      currentRun = event->runNumber;
      currentSet = &lumiList[currentRun];
    }
    else if(currentLumi == event->luminosityBlockNumber) continue;

    currentLumi = event->luminosityBlockNumber;

    currentSet->insert(currentLumi);
  }

  std::ofstream output(_outputName.Data());
  output << "{";
  std::map<int, std::set<int> >::iterator rItr(lumiList.begin());
  while(rItr != lumiList.end()){
    std::set<int>::iterator lItr(rItr->second.begin());
    std::set<int>::iterator lEnd(rItr->second.end());
    if(lItr == lEnd) continue;

    output << "\"" << rItr->first << "\": [[";
    int current(*lItr);
    output << current << ", ";
    ++lItr;
    for(; lItr != lEnd; ++lItr){
      if(*lItr == current + 1) ++current;
      else{
        output << current << "], [";
        current = *lItr;
        output << current << ", ";
      }
    }
    output << current << "]]";

    if(++rItr != lumiList.end())
      output << ", ";
    else
      break;
  }

  output << "}";

  output.close();

  delete event;
}
