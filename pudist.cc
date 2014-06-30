#include "TChain.h"
#include "TH1F.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "ObjectTree.h"

#include <cmath>
#include <iostream>

/*
  Plot PU distributions of Z->mumu events
*/

void
pudist(TString const& _sourceName, TString const& _outputName, bool _MC)
{
  TChain objTree("allObjects");
  TChain evtTree("eventVars");
  objTree.Add(_sourceName);
  if(_MC) evtTree.Add(_sourceName);

  objTree.SetBranchStatus("*", 0);
  objTree.SetBranchStatus("muon.size", 1);
  objTree.SetBranchStatus("muon.isTight", 1);
  objTree.SetBranchStatus("muon.px", 1);
  objTree.SetBranchStatus("muon.py", 1);
  objTree.SetBranchStatus("muon.pz", 1);
  objTree.SetBranchStatus("muon.energy", 1);
  objTree.SetBranchStatus("vertex.size", 1);
  objTree.SetBranchStatus("vertex.isGood", 1);
  if(_MC){
    evtTree.SetBranchStatus("*", 0);
    evtTree.SetBranchStatus("puWeight", 1);
  }

  susy::MuonVarsArray muons;
  susy::VertexVarsArray vertices;
  float puWeight(1.);

  muons.setAddress(objTree);
  vertices.setAddress(objTree);
  if(_MC) evtTree.SetBranchAddress("puWeight", &puWeight);

  TFile* outputFile(TFile::Open(_outputName, "recreate"));
  TH1F* h_raw(new TH1F("h_raw", "N_{vtx} distribution;N_{vtx}", 50, 0., 50.));
  TH1F* h_weighted(0);
  if(_MC){
    h_weighted = new TH1F("h_weighted", "Weighted N_{vtx} distribution;N_{vtx}", 50, 0., 50.);
    h_weighted->Sumw2();
  }

  long iEntry(0);
  while(objTree.GetEntry(iEntry++) > 0){
    if(iEntry % 100000 == 1) (std::cout << "\r" << iEntry).flush();

    if(_MC) evtTree.GetEntry(iEntry - 1);

    unsigned iM1(0);
    for(; iM1 != muons.size; ++iM1){
      if(!muons.isTight[iM1]) continue;
      TLorentzVector p1(muons.px[iM1], muons.py[iM1], muons.pz[iM1], muons.energy[iM1]);
      unsigned iM2(0);
      for(; iM2 != muons.size; ++iM2){
        if(iM2 == iM1) continue;
        if(!muons.isTight[iM2]) continue;
        double mass(std::sqrt((muons.energy[iM1] + muons.energy[iM2]) * (muons.energy[iM1] + muons.energy[iM2]) -
                              (muons.px[iM1] + muons.px[iM2]) * (muons.px[iM1] + muons.px[iM2]) -
                              (muons.py[iM1] + muons.py[iM2]) * (muons.py[iM1] + muons.py[iM2]) -
                              (muons.pz[iM1] + muons.pz[iM2]) * (muons.pz[iM1] + muons.pz[iM2])));
        if(mass > 80. && mass < 100.) break;
      }
      if(iM2 != muons.size) break;
    }
    if(iM1 == muons.size) continue;

    double nVtx(0.);
    for(unsigned iV(0); iV != vertices.size; ++iV)
      if(vertices.isGood[iV]) nVtx += 1.;

    h_raw->Fill(nVtx);
    if(_MC) h_weighted->Fill(nVtx, puWeight);
  }

  std::cout << std::endl;

  outputFile->cd();
  outputFile->Write();
  delete outputFile;
}
