import FWCore.ParameterSet.Config as cms

from RecoTauTag.RecoTau.PFRecoTauQualityCuts_cfi import PFTauQualityCuts

pfRecoTauDiscriminationByLorentzNet = cms.EDProducer("PFRecoTauDiscriminationByLorentzNet2",
  srcPFTaus = cms.InputTag('hpsPFTauProducer'),
  srcPFCandidates = cms.InputTag('particleFlow'),
  addPFCandidatesOutsideJets = cms.bool(True),
  isolationConeSize = cms.double(0.5),
  qualityCuts = PFTauQualityCuts,
  applySignalQualityCuts = cms.bool(False),
  maxChargedPFCand_dZ = cms.double(0.2),
  maxPFCandsPerJet = cms.uint32(100),
  model = cms.FileInPath("TallinnTauTag/LorentzNet/data/LorentzNet2_tau_2022Dec11.onnx"),
  usePdgMass = cms.bool(False),
  verbosity = cms.int32(-1)
)
