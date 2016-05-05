import FWCore.ParameterSet.Config as cms

DimuonFilter = cms.EDFilter('DimuonFilter' ,
                            muonInputTag = cms.InputTag("muons"),
                            csvFileName = cms.string("dimuon.csv"),
                            minMuonPt = cms.double(15.0),
                            maxMuonEta = cms.double(2.4),
                            invariantMassMin = cms.double(0.3),
                            invariantMassMax = cms.double(300.0)
                            )
