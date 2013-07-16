#ifndef Utilities_h
#define Utilities_h

#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include "SusyEvent.h"

#include "TString.h"
#include "TVector2.h"

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

  TString particleName(int, bool = false);

  unsigned countVertices(susy::VertexCollection const&, std::vector<unsigned>* = 0);

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

  template<class T> void sortByPt(std::vector<T const*>& _objects, std::vector<unsigned>* _indices = 0){
    typedef std::map<double, T const*> Sorter;
    typedef typename Sorter::value_type SorterValueType;
    typedef typename Sorter::reverse_iterator SorterIter;

    unsigned nO(_objects.size());

    if(_indices && _indices->size() != nO)
      throw std::runtime_error("susy::sortByPt: objects and indices size do not match");

    Sorter sortedObjects;
    std::map<double, unsigned> sortedIndices;

    for(unsigned iO(0); iO != nO; ++iO){
      double pt(_objects[iO]->momentum.Pt());
      if(pt != pt)
        throw std::runtime_error(TString::Format("Pt of %dth object is NaN", iO).Data());

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
