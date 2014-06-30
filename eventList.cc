#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"

#include <map>
#include <set>
#include <iostream>
#include <fstream>

/*
  Dump event numbers in given susyEvents file
*/

void
eventList(TString const& _sourceDir, TString const& _outputName, unsigned _nJobs = 1, unsigned _iJob = 0)
{
  std::set<TString> fileNames;
  void* dirp(gSystem->OpenDirectory(_sourceDir));
  TString fileName;
  while((fileName = gSystem->GetDirEntry(dirp)) != ""){
    if(fileName == "." || fileName == "..") continue;
    if(!fileName.Contains("susyEvents") || !fileName.Contains(".root")) continue;
    fileNames.insert(_sourceDir + "/" + fileName);
  }

  TChain tree("susyTree");
  unsigned iF(0);
  for(std::set<TString>::iterator fItr(fileNames.begin()); fItr != fileNames.end(); ++fItr, ++iF)
    if(iF % _nJobs == _iJob) tree.Add(*fItr);

  std::cout << tree.GetNtrees() << std::endl;

  tree.SetBranchStatus("*", 0);
  tree.SetBranchStatus("runNumber", 1);
  tree.SetBranchStatus("luminosityBlockNumber", 1);
  tree.SetBranchStatus("eventNumber", 1);

  unsigned runNumber;
  unsigned lumiNumber;
  unsigned eventNumber;

  tree.SetBranchAddress("runNumber", &runNumber);
  tree.SetBranchAddress("luminosityBlockNumber", &lumiNumber);
  tree.SetBranchAddress("eventNumber", &eventNumber);

  std::map<unsigned, std::map<unsigned, std::set<unsigned> > > theList;
  std::map<unsigned, std::map<unsigned, std::vector<unsigned> > > duplicates;

  long iEntry(0);
  while(tree.GetEntry(iEntry++) > 0){
    if(!theList[runNumber][lumiNumber].insert(eventNumber).second){
      std::cout << "Duplicate event " << runNumber << ":" << lumiNumber << ":" << eventNumber << " found in " << tree.GetCurrentFile()->GetName() << std::endl;
      duplicates[runNumber][lumiNumber].push_back(eventNumber);
    }
  }

  std::cout << "====================" << std::endl;

  for(std::map<unsigned, std::map<unsigned, std::vector<unsigned> > >::iterator rItr(duplicates.begin()); rItr != duplicates.end(); ++rItr){
    runNumber = rItr->first;
    for(std::map<unsigned, std::vector<unsigned> >::iterator lItr(rItr->second.begin()); lItr != rItr->second.end(); ++lItr){
      lumiNumber = lItr->first;
      std::set<unsigned>& lumiSet(theList[runNumber][lumiNumber]);
      std::cout << runNumber << ":" << lumiNumber << " has " << lumiSet.size() << " events. ";

      if(lumiSet.size() == lItr->second.size())
        std::cout << "Fully duplicated." << std::endl;
      else{
        std::cout << "Doubly duplicated:" << std::endl;
        std::set<unsigned> uniqify;
        for(unsigned iE(0); iE != lItr->second.size(); ++iE)
          if(!uniqify.insert(lItr->second[iE]).second) std::cout << " " << lItr->second[iE];
      }
    }
  }

  std::ofstream dump(_outputName.Data());
  for(std::map<unsigned, std::map<unsigned, std::set<unsigned> > >::iterator rItr(theList.begin()); rItr != theList.end(); ++rItr){
    dump << "[" << rItr->first << "]" << std::endl;
    for(std::map<unsigned, std::set<unsigned> >::iterator lItr(rItr->second.begin()); lItr != rItr->second.end(); ++lItr){
      dump << lItr->first << ":" << std::endl;
      for(std::set<unsigned>::iterator eItr(lItr->second.begin()); eItr != lItr->second.end(); ++eItr)
        dump << " " << *eItr << std::endl;
    }
    dump << std::endl;
  }
  dump.close();
}
