#ifndef PFParticleBugFix_h
#define PFParticleBugFix_h

#include "SusyEvent.h"

#include <vector>

namespace susy {
  std::vector<const PFParticle*> cleanPFParticles(PFParticleCollection const&);
}

#endif
