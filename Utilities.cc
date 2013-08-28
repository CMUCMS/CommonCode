#include "Utilities.h"

#include "ObjectVars.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include "TArrayI.h"
#include "TPRegexp.h"
#include "TObjArray.h"

namespace susy {

  double
  deltaR2(double _eta1, double _phi1, double _eta2, double _phi2)
  {
    double dEta(_eta1 - _eta2);
    double dPhi(TVector2::Phi_mpi_pi(_phi1 - _phi2));
    return dEta * dEta + dPhi * dPhi;
  }

  double
  deltaR(double _eta1, double _phi1, double _eta2, double _phi2)
  {
    return std::sqrt(deltaR2(_eta1, _phi1, _eta2, _phi2));
  }

  TString
  particleName(int _pdgId, bool _charged/* = false*/)
  {
    TString name;
    int charge;

    switch(std::abs(_pdgId)){
    case 11:
      name = "e";
      charge = -1;
      break;
    case 12:
      name = "#nu_{e}";
      charge = 0;
      break;
    case 13:
      name = "#mu";
      charge = -1;
      break;
    case 14:
      name = "#nu_{#mu}";
      charge = 0;
      break;
    case 15:
      name = "#tau";
      charge = -1;
      break;
    case 16:
      name = "#nu_{#tau}";
      charge = 0;
      break;
    case 22:
      name = "#gamma";
      charge = 0;
      break;
    case 111:
      name = "#pi^{0}";
      charge = 0;
      break;
    case 130:
      name = "K^{0}_{L}";
      charge = 0;
      break;
    case 211:
      name = "#pi";
      charge = 1;
      break;
    case 221:
      name = "#eta";
      charge = 0;
      break;
    case 310:
      name = "K^{0}_{S}";
      charge = 0;
      break;
    case 321:
      name = "K";
      charge = 1;
      break;
    case 2112:
      name = "n";
      charge = 0;
      break;
    case 2212:
      name = "p";
      charge = 1;
      break;
    case 3112:
      name = "#Sigma^{-}";
      charge = 0;
      break;
    case 3122:
      name = "#Lambda";
      charge = 0;
      break;
    case 3222:
      name = "#Sigma^{+}";
      charge = 0;
      break;
    case 3312:
      name = "#Xi";
      charge = 1;
      break;
    case 3322:
      name = "#Xi^{0}";
      charge = 0;
      break;
    case 3334:
      name = "#Omega";
      charge = -1;
      break;
    default:
      name = TString::Format("%d", _pdgId);
      charge = 0;
      break;
    }

    if(_charged){
      if(charge == 0){
        if(_pdgId < 0) name = "#bar{" + name + "}";
      }
      else if(charge == 1 || charge == -1){
        if(charge * _pdgId > 0) name += "^{+}";
        else name += "^{-}";
      }
    }
    else if(charge == 1 || charge == -1) name += "^{#pm}";

    return name;
  }

  GoodLumis::GoodLumis() :
    list_(),
    isGood_(true),
    run_(0),
    lumi_(0)
  {}

  bool
  GoodLumis::parseJSON(TString const& _fileName)
  {
    if(_fileName == "") return true;

    if(_fileName.Contains(",")){
      TObjArray* names(_fileName.Tokenize(","));
      names->SetOwner(true);
      bool result(true);
      for(int i(0); i != names->GetEntries(); ++i)
        result &= parseJSON(names->At(i)->GetName());
      delete names;
      return result;
    }

    std::ifstream inputFile(_fileName);
    if(!inputFile.is_open()){
      std::cerr << "Cannot open JSON file " << _fileName << std::endl;
      return false;
    }

    std::string line;
    TString jsonText;
    std::getline(inputFile, line);
    do{
      jsonText += line;
      std::getline(inputFile, line);
    }while(inputFile.good());
    inputFile.close();

    TPRegexp runBlockPat("\"([0-9]+)\":[ ]*\\[((?:\\[[0-9]+,[ ]*[0-9]+\\](?:,[ ]*|))+)\\]");
    TPRegexp lumiBlockPat("\\[([0-9]+),[ ]*([0-9]+)\\]");

    TArrayI positions(2);
    positions[1] = 0;
    while(runBlockPat.Match(jsonText, "g", positions[1], 10, &positions) == 3){
      TString runBlock(jsonText(positions[0], positions[1] - positions[0]));
      TString lumiPart(jsonText(positions[4], positions[5] - positions[4]));

      int run(TString(jsonText(positions[2], positions[3] - positions[2])).Atoi());
      std::set<int>& lumis(list_[run]);

      TArrayI lumiPos(2);
      lumiPos[1] = 0;
      while(lumiBlockPat.Match(lumiPart, "g", lumiPos[1], 10, &lumiPos) == 3){
        TString lumiBlock(lumiPart(lumiPos[0], lumiPos[1] - lumiPos[0]));
        int begin(TString(lumiPart(lumiPos[2], lumiPos[3] - lumiPos[2])).Atoi());
        int end(TString(lumiPart(lumiPos[4], lumiPos[5] - lumiPos[4])).Atoi());
        for(int lumi(begin); lumi <= end; ++lumi)
          lumis.insert(lumi);
      }
    }

    return true;
  }

  bool
  GoodLumis::isGoodLumi(unsigned _run, unsigned _lumi) const
  {
    if(list_.size() == 0) return true;

    if(_run == run_ && _lumi == lumi_) return isGood_;
    run_ = _run;
    lumi_ = _lumi;

    std::map<int, std::set<int> >::const_iterator rItr(list_.find(run_));
    if(rItr == list_.end()) return (isGood_ = false);
    std::set<int> const& lumis(rItr->second);
    return (isGood_ = (lumis.find(lumi_) != lumis.end()));
  }

}
