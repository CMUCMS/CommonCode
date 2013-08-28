#include "PFParticleBugFix.h"

#include <map>
#include <utility>
#include <cmath>

namespace susy {

  std::vector<const PFParticle*>
  cleanPFParticles(PFParticleCollection const& _input)
  {
    std::vector<const PFParticle*> output;

    std::map<std::pair<double, double>, const PFParticle*> cleanedMap;

    unsigned nPF(_input.size());
    for(unsigned iP(0); iP != nPF; ++iP){
      PFParticle const& particle(_input[iP]);
      double pt(particle.momentum.Pt());
      if(pt < 3.) continue;
      if(std::abs(particle.momentum.Eta()) > etaMax) continue;

      cleanedMap[std::pair<double, double>(pt, particle.momentum.Eta())] = &particle;
    }

    std::map<std::pair<double, double>, const PFParticle*>::iterator mEnd(cleanedMap.end());
    for(std::map<std::pair<double, double>, const PFParticle*>::iterator mItr(cleanedMap.begin()); mItr != mEnd; ++mItr)
      output.push_back(mItr->second);

    return output;
  }

}
