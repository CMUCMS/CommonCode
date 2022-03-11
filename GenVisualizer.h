#ifndef GenVisualizer_h
#define GenVisualizer_h

#include "TCanvas.h"
#include "TH2F.h"
#include "TTree.h"
#include "TString.h"
#include "TMarker.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TLegend.h"

namespace susy {

  class GenVisualizer {
  public:
    GenVisualizer();
    ~GenVisualizer();

    void setAddress(TTree&);
    void process(TString const& = "");
    void print(TString const&);
    void zoomIn(double, double, double);
    void zoomOut();
    void labelParticles();
    void setMaxPt(double);
    void showGenJet(bool _val) { showGenJet_ = _val; }
    void showRecoJet(bool _val) { showRecoJet_ = _val; }

    TString plotsDir;
    double maxPt;

  private:
    int getPtColor_(double);
    unsigned formGenJets_(int*, double*, double*, double*);

    unsigned eventNumber_;
    unsigned genSize_;
    float genPt_[256];
    float genEta_[256];
    float genPhi_[256];
    unsigned short genStatus_[256];
    int genPdgId_[256];
    short genMotherIndex_[256];
    float met_;
    float metPhi_;
    unsigned jetSize_;
    float jetPt_[256];
    float jetEta_[256];
    float jetPhi_[256];
    bool jetPassPUJetIdLoose_[256];
    bool jetIsLoose_[256];

    enum PType {
      kPhoton,
      kMuon,
      kElectron,
      kHadron,
      kNeutrino,
      kLSP,
      nPType
    };
  
    TCanvas* canvas_;
    TH2F* frame_;
    TMarker marker_[nPType];
    TMarker wMarker_;
    TLine connector_;
    TEllipse jetCone_;
    TEllipse puJetCone_;
    TEllipse genJetCone_;
    TLine metLine_;
    TLine sumPtLine_;
    TLegend legend_;

    double maxPt_;
    int const nColors_;
    int col0_;
    bool showGenJet_;
    bool showRecoJet_;
  };

}

#endif
