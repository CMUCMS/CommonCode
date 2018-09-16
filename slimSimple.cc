#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjArray.h"

#include <iostream>

void
slimSimple(TString const& _sourcePath, TString const& _outputPath)
{
  TChain obj("allObjects");
  obj.Add(_sourcePath);
  TChain event("eventVars");
  event.Add(_sourcePath);

  TFile* outputFile(TFile::Open(_outputPath, "recreate"));
  TTree* output(obj.CloneTree(0));
  output->SetName("events");

  TObjArray* branches(output->GetListOfBranches());
  TObjArray* leaves(output->GetListOfLeaves());
  int nB(branches->GetEntriesFast());
  for(int iB(0); iB != nB; ++iB){
    TBranch* branch(static_cast<TBranch*>(branches->UncheckedAt(iB)));
    TString bname(branch->GetName());
    if(bname.BeginsWith("photon.")){
      if(bname.EndsWith(".size") ||
         bname.EndsWith(".px") ||
         bname.EndsWith(".py") ||
         bname.EndsWith(".pz") ||
         bname.EndsWith(".energy") ||
         bname.EndsWith(".pt") ||
         bname.EndsWith(".eta") ||
         bname.EndsWith(".phi") ||
         bname.EndsWith(".isTight")) continue;
    }
    else if(bname.BeginsWith("electron.")){
      if(bname.EndsWith(".size") ||
         bname.EndsWith(".px") ||
         bname.EndsWith(".py") ||
         bname.EndsWith(".pz") ||
         bname.EndsWith(".energy") ||
         bname.EndsWith(".pt") ||
         bname.EndsWith(".eta") ||
         bname.EndsWith(".phi") ||
         bname.EndsWith(".isTight")) continue;
    }
    else if(bname.BeginsWith("muon.")){
      if(bname.EndsWith(".size") ||
         bname.EndsWith(".px") ||
         bname.EndsWith(".py") ||
         bname.EndsWith(".pz") ||
         bname.EndsWith(".energy") ||
         bname.EndsWith(".pt") ||
         bname.EndsWith(".eta") ||
         bname.EndsWith(".phi") ||
         bname.EndsWith(".isTight")) continue;
    }
    else if(bname.BeginsWith("jet.")){
      if(bname.EndsWith(".size") ||
         bname.EndsWith(".px") ||
         bname.EndsWith(".py") ||
         bname.EndsWith(".pz") ||
         bname.EndsWith(".energy") ||
         bname.EndsWith(".pt") ||
         bname.EndsWith(".eta") ||
         bname.EndsWith(".phi") ||
         bname.EndsWith(".passPUJetIdLoose")) continue;
    }
    
    branches->RemoveAt(iB);
    delete branch;
    branches->Compress();
    --nB;
    --iB;
  }
  leaves->Compress();

  unsigned runNumber;
  unsigned eventNumber;
  bool isSoft[256];
  event.SetBranchStatus("*", 0);
  event.SetBranchStatus("runNumber", 1);
  event.SetBranchStatus("eventNumber", 1);
  event.SetBranchAddress("runNumber", &runNumber);
  event.SetBranchAddress("eventNumber", &eventNumber);
  output->Branch("runNumber", &runNumber, "runNumber/i");
  output->Branch("eventNumber", &eventNumber, "eventNumber/i");
  output->Branch("muon.isSoft", isSoft, "isSoft[muon.size]/O");

  unsigned size;
  unsigned char nLayersWithMmt[256];
  unsigned char nValidPixelHits[256];
  float normChi2[256];
  float dxy[256];
  float dz[256];
  obj.SetBranchAddress("muon.nLayersWithMmt", nLayersWithMmt);
  obj.SetBranchAddress("muon.nValidPixelHits", nValidPixelHits);
  obj.SetBranchAddress("muon.normChi2", normChi2);
  obj.SetBranchAddress("muon.dxy", dxy);
  obj.SetBranchAddress("muon.dz", dz);

  obj.SetBranchAddress("muon.size", &size);

  long iEntry(0);
  while(obj.GetEntry(iEntry++) > 0){
    if(size < 2) continue;
    event.GetEntry(iEntry - 1);

    for(unsigned iM(0); iM != size; ++iM)
      isSoft[iM] = nLayersWithMmt[iM] > 5 && nValidPixelHits[iM] > 0 && normChi2[iM] < 10. && dxy[iM] < 0.3 && dz[iM] < 20.;

    output->Fill();
  }

  outputFile->cd();
  output->Write();
  delete outputFile;
}
