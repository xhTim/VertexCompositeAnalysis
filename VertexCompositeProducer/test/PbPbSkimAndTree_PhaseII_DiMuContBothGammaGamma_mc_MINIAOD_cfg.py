import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9

process = cms.Process('ANASKIM',Phase2C11I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D77Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('root://eoscms.cern.ch///store/group/phys_heavyions/huangxi/Reco/STARlight_CohJpsi2MuMu_0n0n_PbPb5TeV_Reco_v2/STARlight_CohJpsi2MuMu_0n0n_PbPb5TeV_GenFilter/STARlight_CohJpsi2MuMu_0n0n_PbPb5TeV_Reco_v2/230711_222047/0000/step3_inMINIAODSIM_100.root'),
   inputCommands=cms.untracked.vstring('keep *', 'drop *_hiEvtPlane_*_*')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

# Add the VertexComposite producer
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cff")
process.generalDiMuCandidatesWrongSign = process.generalDiMuCandidates.clone(isWrongSign = cms.bool(True))
from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import doPATMuons
doPATMuons(process, True)

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.colEvtSel = cms.Sequence(process.hfCoincFilter2Th4 * process.primaryVertexFilterAA )

# place holder for changeToMiniAOD
# do not perform selections
process.eventFilter_HM = cms.Sequence(
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

# Define the analysis steps
process.dimurereco_step = cms.Path(process.eventFilter_HM * process.patMuonSequence * process.generalDiMuCandidates)
process.dimurerecowrongsign_step = cms.Path(process.eventFilter_HM * process.patMuonSequence * process.generalDiMuCandidatesWrongSign)

# Add the VertexComposite tree
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_tree_cff")

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('dimuana_mc.root'))
process.p = cms.EndPath(process.dimucontana_mc * process.dimucontana_wrongsign_mc)
process.dimucontana_mc.isGammaGamma = cms.untracked.bool(True)
process.dimucontana_mc.doMuonNtuple = cms.untracked.bool(True)
process.dimucontana_wrongsign_mc.isGammaGamma = cms.untracked.bool(True)
process.dimucontana_wrongsign_mc.doMuonNtuple = cms.untracked.bool(True)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.dimurereco_step,
    process.dimurerecowrongsign_step,
    process.p
)

# Add the event selection filters
# CaloTower never available in MINIAOD
# the flags are useless
process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_hfCoincFilter2Th4 = cms.Path(process.eventFilter_HM * process.hfCoincFilter2Th4)
process.Flag_primaryVertexFilter = cms.Path(process.eventFilter_HM * process.primaryVertexFilterAA)
eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter2Th4 , process.Flag_primaryVertexFilter ]

# MINIAOD
from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD
changeToMiniAOD(process)
from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"genParticles","prunedGenParticles")
process = MassReplaceInputTag(process,"offlinePrimaryVerticesRecovery","unpackedTracksAndVertices")
