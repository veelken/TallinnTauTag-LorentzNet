import FWCore.ParameterSet.Config as cms

process = cms.Process("producePFJetTaggingNtuple")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/store/mc/Run3Winter22DRPremix/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/AODSIM/122X_mcRun3_2021_realistic_v9_ext2-v2/40003/72e1ed1e-0069-4068-a5ae-9532514b7bc6.root'
    )
)

inputFilePath = '/local/veelken/PFMCTruthMatcher/EDMtuples/2022Oct27/'
inputFile_regex = r"producePFMCTruthEDMtuple_[0-9]+.root"
outputFilePath = None
outputFileName = "producePFRecoTauDiscriminationByLorentzNet.root"
is_signal = True

##inputFilePath = None
##inputFileNames = $inputFileNames
##outputFilePath = "$outputFilePath"
##outputFileName = "$outputFileName"
##is_signal = $is_signal

#--------------------------------------------------------------------------------
# set input files
if inputFilePath:
    from TallinnTauTag.Tools.tools.getInputFileNames import getInputFileNames
    print("Searching for input files in path = '%s'" % inputFilePath)
    inputFileNames = getInputFileNames(inputFilePath, inputFile_regex)
    print("Found %i input files." % len(inputFileNames))
    ##process.source.fileNames = cms.untracked.vstring(inputFileNames)
else:
    print("Processing %i input files: %s" % (len(inputFileNames), inputFileNames))
    ##process.source.fileNames = cms.untracked.vstring(inputFileNames)
outputFileName_full = None
if outputFilePath:
    import os
    outputFileName_full = os.path.join(outputFilePath, outputFileName)
else:
    outputFileName_full = outputFileName
print("Writing output to file: %s" % outputFileName_full)
#--------------------------------------------------------------------------------

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2023_realistic', '')

process.analysisSequence = cms.Sequence()

if is_signal:
  process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
  process.analysisSequence += process.tauGenJets

  process.load("PhysicsTools.JetMCAlgos.TauGenJetsDecayModeSelectorAllHadrons_cfi")
  process.analysisSequence += process.tauGenJetsSelectorAllHadrons

  process.selectedGenHadTaus = cms.EDFilter("GenJetSelector",
    src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
    cut = cms.string('pt > 10. & abs(eta) < 2.5'),
    filter = cms.bool(False)
  )
  process.analysisSequence += process.selectedGenHadTaus

process.selectedAK4PFJets = cms.EDFilter("PFJetSelector",
  src = cms.InputTag('ak4PFJets'),
  cut = cms.string('pt > 18. & abs(eta) < 2.5'),
  filter = cms.bool(False)
)
process.analysisSequence += process.selectedAK4PFJets

pfJetCollection = 'selectedAK4PFJets'
if is_signal:
  process.genMatchedAK4PFJets = cms.EDFilter("PFJetAntiOverlapSelector",
    src = cms.InputTag('selectedAK4PFJets'),
    srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
    dRmin = cms.double(0.5),
    invert = cms.bool(True),
    filter = cms.bool(False)
  )
  process.analysisSequence += process.genMatchedAK4PFJets
  pfJetCollection = 'genMatchedAK4PFJets'

process.selectedPFTaus = cms.EDFilter("PFTauAntiOverlapSelector",
  src = cms.InputTag('hpsPFTauProducer'),
  srcNotToBeFiltered = cms.VInputTag(pfJetCollection),
  dRmin = cms.double(0.3),
  invert = cms.bool(True),
  filter = cms.bool(False)
)
process.analysisSequence += process.selectedPFTaus

# CV: load HepPDT::ParticleData in order to access particle masses given their pdgId
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

process.load("TallinnTauTag/LorentzNet/PFRecoTauDiscriminationByLorentzNet_cfi")
process.pfRecoTauDiscriminationByLorentzNet.srcPFTaus = cms.InputTag('selectedPFTaus')
process.pfRecoTauDiscriminationByLorentzNet.verbosity = cms.int32(1)
process.analysisSequence += process.pfRecoTauDiscriminationByLorentzNet

process.p = cms.Path(process.analysisSequence)

process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        'keep *',
        'drop *_hltScouting*Packer_*_*'
    ),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string(outputFileName_full)
)

process.q = cms.EndPath(process.output)

##dump_file = open('dump.py','w')
##dump_file.write(process.dumpPython())
