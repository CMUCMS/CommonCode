#ifndef HLT_h
#define HLT_h

#include "TString.h"

TString hltPaths[] = {
  "HLT_IsoMu24_eta2p1",
  "HLT_IsoMu24",
  "HLT_Ele27_WP80",
  "HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50",
  "HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50",
  "HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60",
  "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60",
  "HLT_Mu22_Photon22_CaloIdL",
  "HLT_Photon70_CaloIdXL_PFHT400",
  "HLT_Photon70_CaloIdXL_PFNoPUHT400",
  "HLT_Photon70_CaloIdXL_PFHT500",
  "HLT_Photon70_CaloIdXL_PFNoPUHT500",
  "HLT_Photon70_CaloIdXL_PFMET100",
  "HLT_Photon60_CaloIdL_HT300",
  "HLT_Photon60_CaloIdL_MHT70",
  "HLT_Photon135",
  "HLT_Photon150"
};

unsigned const nHLTPaths(sizeof(hltPaths) / sizeof(TString));

// Ph36_OR_Ph22_OR does not save the 'OR' legs - has been verified by Hgg that taking the
// OR offline of IdIso leg and R9Id leg is sufficient
TString hltPhotonFilters[] = {
  "hltL1sL1DoubleEG137", // Photon26_Photon18 seed
  "hltL1sL1SingleEG22", // Photon36_Photon22 seed
  "hltL1sL1SingleEG24", // Photon60/70 seed
  "hltL1sL1SingleEG30", // Photon135/150 seed
  "hltL1sL1Mu3p5EG12", // Mu22_Photon22_CaloIdL seed
  "hltEG36CaloId10Iso50HcalIsoLastFilter", // Ph36_IdIso_Ph22_IdIso first leg up to hcal iso
  "hltEG36R9Id85LastFilter", // Ph36_R9Id_Ph22_R9Id first leg
  "hltEG22CaloId10Iso50TrackIsoDoubleLastFilterUnseeded", // Ph36_IdIso_Ph22_IdIso v>=2 both legs
  "hltEG22R9Id85LastFilterUnseeded", // Ph36_R9Id_Ph22_R9Id both legs
  "hltEG26CaloId10Iso50HcalIsoLastFilter", // Ph26_IdIso_Ph18_IdIso first leg up to hcal iso
  "hltEG18CaloId10Iso50TrackIsoDoubleLastFilterUnseeded", // Ph26_IdIso_Ph18_IdIso v>=2 both legs
  "hltPhoton26CaloId10Iso50Photon18CaloId10Iso50Mass60EgammaCombMassLastFilter", // Ph26_IdIso_Ph18_IdIso mass 60
  "hltMu22Photon22CaloIdLHEFilter", // Mu22_Photon22_CaloIdL
  "hltEG60CaloIdLHEFilter", // Ph60
  "hltEG70CaloIdXLHEFilter", // Ph70
  "hltPhoton135HEFilter", // Ph135
  "hltPhoton150HEFilter" // Ph150
};

unsigned const nHLTPhotonFilters(sizeof(hltPhotonFilters) / sizeof(TString));

TString hltMuonFilters[] = {
  "hltL1sL1Mu3p5EG12", // Mu22_Photon22_CaloIdL seed
  "hltL1sMu16Eta2p1", // IsoMu24_eta2p1 seed
  "hltL1sMu16", // IsoMu24 seed
  "hltL1Mu3p5EG12L3Filtered22", // Mu22_Photon22_CaloIdL
  "hltL2fL1sMu16Eta2p1L1f0L2Filtered16Q", // IsoMu24_eta2p1 first step
  "hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered24Q", // IsoMu24_eta2p1 second step
  "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10", // IsoMu24_eta2p1 v<=12 last step
  "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15", // IsoMu24_eta2p1 v>=13 last step
  "hltL2fL1sMu16L1f0L2Filtered16Q", // IsoMu24 first step
  "hltL3fL1sMu16L1f0L2f16QL3Filtered24Q", // IsoMu24 second step
  "hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15" // IsoMu24 last step
};

unsigned const nHLTMuonFilters(sizeof(hltMuonFilters) / sizeof(TString));

TString hltElectronFilters[] = {
  "hltL1sL1DoubleEG137", // Photon26_Photon18 seed
  "hltL1sL1SingleEG22", // Ele27_WP80 v<=? Photon36_Photon22 seed
  "hltL1sL1SingleEG20ORL1SingleEG22", // Ele27_WP80 v>=? seed
  "hltEle27WP80TrackIsoFilter", // Ele27_WP80
  "hltEG36CaloId10Iso50HcalIsoLastFilter", // Ph36_IdIso_Ph22_IdIso first leg up to hcal iso
  "hltEG36R9Id85LastFilter", // Ph36_R9Id_Ph22_R9Id first leg
  "hltEG22CaloId10Iso50TrackIsoDoubleLastFilterUnseeded", // Ph36_IdIso_Ph22_IdIso both legs
  "hltEG22R9Id85LastFilterUnseeded", // Ph36_R9Id_Ph22_R9Id both legs
  "hltEG26CaloId10Iso50HcalIsoLastFilter", // Ph26_IdIso_Ph18_IdIso first leg up to hcal iso
  "hltEG18CaloId10Iso50TrackIsoDoubleLastFilterUnseeded", // Ph26_IdIso_Ph18_IdIso both legs
  "hltPhoton26CaloId10Iso50Photon18CaloId10Iso50Mass60EgammaCombMassLastFilter" // Ph26_IdIso_Ph18_IdIso mass 60
};

unsigned const nHLTElectronFilters(sizeof(hltElectronFilters) / sizeof(TString));

TString hltEventFilters[] = {
  "hltMET80",
  "hltPFHT400",
  "hltPFHT500",
  "hltPFMET100"
};

unsigned const nHLTEventFilters(sizeof(hltEventFilters) / sizeof(TString));

#endif
