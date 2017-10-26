import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("outputAnalyzerDiMuHighMass")

process.load("FWCore.MessageService.MessageLogger_cfi")

# geometry and global tag:

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('PhysicsTools.PatAlgos.patSequences_cff')
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '76X_dataRun2_v15', '')
#**************************************************************

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [ 
'root://cms-xrd-global.cern.ch//store/user/malcles/DoubleMuon/DoubleMuon-Run2017C-PromptReco-NewMini-v1/171018_142815/0000/miniAOD-prod_PAT_my_1.root',
'root://cms-xrd-global.cern.ch//store/user/malcles/DoubleMuon/DoubleMuon-Run2017C-PromptReco-NewMini-v1/171018_142815/0000/miniAOD-prod_PAT_my_10.root',
'root://cms-xrd-global.cern.ch//store/user/malcles/DoubleMuon/DoubleMuon-Run2017C-PromptReco-NewMini-v1/171018_142815/0000/miniAOD-prod_PAT_my_100.root'

#       'file:/afs/cern.ch/user/m/malcles/MonScratch/private/VertexCheck/CMSSW_9_2_1/src/miniAOD-prod_PAT_my.root'
]);
#readFiles.extend( [ 
#'/store/user/malcles/DoubleMuon/DoubleMuon-Run2017C-PromptReco-v2/171013_194312/0000/miniAOD-prod_PAT_DATA25ns_AllVerticesCollFiltered_1.root',
#'/store/user/malcles/DoubleMuon/DoubleMuon-Run2017C-PromptReco-v2/171013_194312/0000/miniAOD-prod_PAT_DATA25ns_AllVerticesCollFiltered_10.root',
#'/store/user/malcles/DoubleMuon/DoubleMuon-Run2017C-PromptReco-v2/171013_194312/0000/miniAOD-prod_PAT_DATA25ns_AllVerticesCollFiltered_100.root',
#'/store/user/malcles/DoubleMuon/DoubleMuon-Run2017C-PromptReco-v2/171013_194312/0000/miniAOD-prod_PAT_DATA25ns_AllVerticesCollFiltered_101.root',
#'/store/user/malcles/DoubleMuon/DoubleMuon-Run2017C-PromptReco-v2/171013_194312/0000/miniAOD-prod_PAT_DATA25ns_AllVerticesCollFiltered_102.root',
#'/store/user/malcles/DoubleMuon/DoubleMuon-Run2017C-PromptReco-v2/171013_194312/0000/miniAOD-prod_PAT_DATA25ns_AllVerticesCollFiltered_103.root',
#'/store/user/malcles/DoubleMuon/DoubleMuon-Run2017C-PromptReco-v2/171013_194312/0000/miniAOD-prod_PAT_DATA25ns_AllVerticesCollFiltered_104.root',
#'/store/user/malcles/DoubleMuon/DoubleMuon-Run2017C-PromptReco-v2/171013_194312/0000/miniAOD-prod_PAT_DATA25ns_AllVerticesCollFiltered_105.root',
#'/store/user/malcles/DoubleMuon/DoubleMuon-Run2017C-PromptReco-v2/171013_194312/0000/miniAOD-prod_PAT_DATA25ns_AllVerticesCollFiltered_106.root'
#       ]);
secFiles.extend( [
               ] )



#**************************************************************

process.load("flashgg/MicroAOD/flashggMicroAODSequenceNoMuData_cff")


process.myPatMuons = cms.EDFilter("PATMuonSelector",
#    src = cms.InputTag("selectedPatMuons"),
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string('pt>10.0')
)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("VtxTreeDATAHighMass.root")
)

process.commissioning = cms.EDAnalyzer('outputAnalyzerDiMuHighMass',
                                       VertexTag=cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
                                       patMuonTag =  cms.InputTag('myPatMuons'),
                                       VertexCandidateMapTag=cms.InputTag('flashggVertexMapUnique'),
                                       VertexCandidateMapTagNoMu=cms.InputTag('flashggVertexMapUniqueNoMu'),
                                       genEventInfoProduct=cms.InputTag('generator'),
                                       VertexSelectorName=cms.string("FlashggLegacyVertexSelector"),
                                       vertexIdMVAweightfile = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxId_SL_2016.xml"),
                                       vertexProbMVAweightfile = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxProb_SL_2016.xml"),
                                       bits                      = cms.InputTag("TriggerResults::HLT"),
#                                       pathName=cms.vstring("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"), 
#                                       pathName=cms.vstring("HLT_Mu30_TkMu11"), 
                                       pathName=cms.vstring("HLT_IsoMu24"), 
                                       useSingleLeg            = cms.bool(True),
                                       nVtxSaveInfo            = cms.untracked.uint32(99),
                                       trackHighPurity         = cms.bool(False),
                                       pureGeomConvMatching    = cms.bool(True),
                                       dRexclude               = cms.double(0.05),
                                       sigma1Pix               = cms.double(0.0125255),
                                       sigma1Tib               = cms.double(0.716301),
                                       sigma1Tob               = cms.double(3.17615),
                                       sigma1PixFwd            = cms.double(0.0581667),
                                       sigma1Tid               = cms.double(0.38521),
                                       sigma1Tec               = cms.double(1.67937),
                                       sigma2Pix               = cms.double(0.0298574),
                                       sigma2Tib               = cms.double(0.414393),
                                       sigma2Tob               = cms.double(1.06805),
                                       sigma2PixFwd            = cms.double(0.180419),
                                       sigma2Tid               = cms.double(0.494722),
                                       sigma2Tec               = cms.double(1.21941),
                                       singlelegsigma1Pix      = cms.double(0.0178107),
                                       singlelegsigma1Tib      = cms.double(1.3188),
                                       singlelegsigma1Tob      = cms.double(2.23662),
                                       singlelegsigma1PixFwd   = cms.double(0.152157),
                                       singlelegsigma1Tid      = cms.double(0.702755),
                                       singlelegsigma1Tec      = cms.double(2.46599),
                                       singlelegsigma2Pix      = cms.double(0.0935307),
                                       singlelegsigma2Tib      = cms.double(0.756568),
                                       singlelegsigma2Tob      = cms.double(0.62143),
                                       singlelegsigma2PixFwd   = cms.double(0.577081),
                                       singlelegsigma2Tid      = cms.double(0.892751),
                                       singlelegsigma2Tec      = cms.double(1.56638)       

)
#**************************************************************

process.p = cms.Path(process.flashggMicroAODSequence*process.myPatMuons*process.commissioning)
