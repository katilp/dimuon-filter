import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("opendata")

goodJSON = 'Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt'
myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')

import FWCore.Utilities.FileUtils as FileUtils
from FWCore.MessageLogger.MessageLogger_cfi import *

doubleMuFiles = FileUtils.loadListFromFile('CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10000_file_index.txt')
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(*doubleMuFiles)
                            )

process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
process.source.lumisToProcess.extend(myLumis)

process.load('DimuonFilter.DimuonFilter.DimuonFilter_cfi')

process.DimuonFilter.csvFileName = cms.string('DoubleMuRun2011A.csv')
process.DimuonFilter.minMuonPt = cms.double(15.0)
process.DimuonFilter.maxMuonEta = cms.double(2.4)
process.DimuonFilter.invariantMassMin = cms.double(0.3)
process.DimuonFilter.invariantMassMax = cms.double(300.0)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.mypath = cms.Path(process.DimuonFilter)
process.schedule = cms.Schedule(process.mypath)
