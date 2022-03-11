#include "TChain.h"
#include "TString.h"

#include <map>
#include <set>
#include <iostream>

void
compareEventByEvent(TString const& _source1, TString const& _source2 = "")
{
  unsigned runNumber;
  unsigned eventNumber;

  TChain input1("susyTree");
  input1.Add(_source1);
  input1.SetBranchStatus("*", 0);
  input1.SetBranchStatus("runNumber", 1);
  input1.SetBranchStatus("eventNumber", 1);

  input1.SetBranchAddress("runNumber", &runNumber);
  input1.SetBranchAddress("eventNumber", &eventNumber);

  std::map<unsigned, std::set<unsigned> > list1;

  long iEntry(0);
  while(input1.GetEntry(iEntry++) > 0){
    if(!list1[runNumber].insert(eventNumber).second)
      std::cout << "Duplicate in 1: " << runNumber << " " << eventNumber << std::endl;
  }

  if(_source2 != ""){
    TChain input2("susyTree");

    input2.Add(_source2);
    input2.SetBranchStatus("*", 0);
    input2.SetBranchStatus("runNumber", 1);
    input2.SetBranchStatus("eventNumber", 1);

    input2.SetBranchAddress("runNumber", &runNumber);
    input2.SetBranchAddress("eventNumber", &eventNumber);

    std::map<unsigned, std::set<unsigned> > list2;

    iEntry = 0;
    while(input2.GetEntry(iEntry++) > 0){
      if(!list2[runNumber].insert(eventNumber).second)
        std::cout << "Duplicate in 2: " << runNumber << " " << eventNumber << std::endl;
    }

    std::set<unsigned> runs;
    for(std::map<unsigned, std::set<unsigned> >::iterator rItr(list1.begin()); rItr != list1.end(); ++rItr)
      runs.insert(rItr->first);
    for(std::map<unsigned, std::set<unsigned> >::iterator rItr(list2.begin()); rItr != list2.end(); ++rItr)
      runs.insert(rItr->first);

    for(std::set<unsigned>::iterator rItr(runs.begin()); rItr != runs.end(); ++rItr){
      std::set<unsigned>& evt1(list1[*rItr]);
      std::set<unsigned>& evt2(list2[*rItr]);
      for(std::set<unsigned>::iterator eItr(evt1.begin()); eItr != evt1.end(); ++eItr){
        if(evt2.find(*eItr) == evt2.end())
          std::cout << "2 Does not contain " << *rItr << " " << *eItr << std::endl;
      }
      for(std::set<unsigned>::iterator eItr(evt2.begin()); eItr != evt2.end(); ++eItr){
        if(evt1.find(*eItr) == evt1.end())
          std::cout << "1 Does not contain " << *rItr << " " << *eItr << std::endl;
      }
    }
  }
}  
