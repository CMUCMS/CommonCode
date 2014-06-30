#ifndef Utilities_h
#define Utilities_h

#include "ObjectTree.h"
#include "SimpleEventProducer.h"

#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include "TString.h"
#include "TVector2.h"
#include "TLorentzVector.h"

namespace susy {

  double deltaR2(double, double, double, double);

  template<class T1, class T2> double deltaR2(T1 const& _obj1, T2 const& _obj2)
  {
    return deltaR2(_obj1.eta, _obj1.phi, _obj2.eta, _obj2.phi);
  }

  double deltaR(double, double, double, double);

  template<class T1, class T2> double deltaR(T1 const& _obj1, T2 const _obj2)
  {
    return deltaR(_obj1.eta, _obj1.phi, _obj2.eta, _obj2.phi);
  }

  template<class T1, class T2> double diObjectMass(T1 const& _o1, T2 const& _o2)
  {
    TLorentzVector p1(_o1.px, _o1.py, _o1.pz, _o1.energy);
    TLorentzVector p2(_o2.px, _o2.py, _o2.pz, _o2.energy);
    return (p1 + p2).M();
  }

  template<class T1, class T2, class T3> double triObjectMass(T1 const& _o1, T2 const& _o2, T3 const& _o3)
  {
    TLorentzVector p1(_o1.px, _o1.py, _o1.pz, _o1.energy);
    TLorentzVector p2(_o2.px, _o2.py, _o2.pz, _o2.energy);
    TLorentzVector p3(_o3.px, _o3.py, _o3.pz, _o3.energy);
    return (p1 + p2 + p3).M();
  }

  TString particleName(int, bool = false);

  void genMatch(SimpleEventProducer::EventVars const&, PhotonVarsArray const*, ElectronVarsArray const*, MuonVarsArray const*, unsigned*, unsigned*, unsigned*, double* = 0);

  class GoodLumis {
  public:
    GoodLumis();
    ~GoodLumis(){};

    bool parseJSON(TString const&);
    bool isGoodLumi(unsigned, unsigned) const;

  private:
    std::map<int, std::set<int> > list_;
    mutable bool isGood_;
    mutable unsigned run_;
    mutable unsigned lumi_; 
  };

  class Exception : public std::runtime_error {
  public:
    enum Categories {
      kObjectAnomaly,
      kEventAnomaly,
      kIOError,
      kFormatError,
      nCategories
    };

    explicit Exception(Categories _cat, std::string const& _what) :
      runtime_error(_what),
      category(_cat)
    {}

    std::string categoryName(){
      switch(category){
      case kObjectAnomaly:
        return "ObjectAnomaly";
      case kEventAnomaly:
        return "EventAnomaly";
      case kIOError:
        return "IOError";
      case kFormatError:
        return "FormatError";
      default:
        return "Unknown";
      }
    }

    Categories category;
  };

  template<class T> void sortByPt(std::vector<T const*>& _objects, std::vector<unsigned>* _indices = 0){
    typedef std::map<double, T const*> Sorter;
    typedef typename Sorter::value_type SorterValueType;
    typedef typename Sorter::reverse_iterator SorterIter;

    unsigned nO(_objects.size());

    if(_indices && _indices->size() != nO)
      throw std::range_error("susy::sortByPt: objects and indices size do not match");

    Sorter sortedObjects;
    std::map<double, unsigned> sortedIndices;

    for(unsigned iO(0); iO != nO; ++iO){
      double pt(_objects[iO]->momentum.Pt());
      if(pt != pt)
        throw Exception(Exception::kObjectAnomaly, TString::Format("Pt of %dth object is NaN", iO).Data());

      while(!sortedObjects.insert(SorterValueType(pt, _objects[iO])).second) pt += 1.e-5;
      if(_indices) sortedIndices[pt] = (*_indices)[iO];
    }

    _objects.clear();
    if(_indices) _indices->clear();

    std::map<double, unsigned>::reverse_iterator iItr;
    if(_indices) iItr = sortedIndices.rbegin();

    SorterIter oEnd(sortedObjects.rend());
    for(SorterIter oItr(sortedObjects.rbegin()); oItr != oEnd; ++oItr){
      _objects.push_back(oItr->second);

      if(_indices){
        _indices->push_back(iItr->second);
        ++iItr;
      }
    }
  }
    
}

#endif
