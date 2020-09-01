import FWCore.ParameterSet.Config as cms

analysisSkimContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *',
      'keep *_offlinePrimaryVertices_*_*',
      'keep *_*_vertices_*',
      'keep *_offlineBeamSpot_*_*',
      'keep *_TriggerResults_*_*',
      'keep *_hltTriggerSummaryAOD_*_*',
      'keep *_gtStage2Digis_*_*',
      'keep recoCentrality_*_*_*',
      'keep *_*_centralityBin_*'
      'keep *_QWzdcreco_*_*',
      'keep recoEvtPlanes_*_*_*',
      'keep LumiProducerFromBrilcalc_*_*_*',
      'keep LumiScalerss_*_*_*',
      'keep *_externalLHEProducer_*_*',
      'keep GenEventInfoProduct_*_*_*',
      'keep recoGenParticles_*_*_*',
      'keep patGenericParticles_*_*_*',
      'keep *_*_nTracks_*',
      'keep LumiInfo_*_*_*',
      'keep *_onlineMetaDataDigis_*_*',
      )
    )