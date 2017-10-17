from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

#config.General.requestName = 'VtxHighMass8020-ZDATA25ns_RunB'
#config.General.requestName = 'VtxHighMass8020-ZDATA25ns_RunCBis'
#config.General.requestName = 'VtxHighMass8020-ZDATA25ns_RunDBis'
#config.General.requestName = 'VtxHighMass8020-ZDATA25ns_RunEBis'
#config.General.requestName = 'VtxHighMass8020-ZDATA25ns_RunF'
#config.General.requestName = 'VtxHighMass8020-ZDATA25ns_RunG'
#config.General.requestName = 'VtxHighMass8020-ZDATA25ns_RunH'# this is v3
config.General.requestName = 'VtxHighMass8020-ZDATA25ns_RunHv2'



config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True 

config.JobType.psetName = 'Validation/python/outputAnalyzerHighMassDataChangeHLT.py'


#config.Data.inputDataset = '/DoubleMuon/malcles-DoubleMuon-Run2016B-23Sep2016-v3-AllVerticesCollFiltered-cc71a0e7222190a450bea2ed0a33ec4b/USER'
#config.Data.inputDataset = '/DoubleMuon/malcles-DoubleMuon-Run2016C-23Sep2016-v1-AllVerticesCollFiltered-cc71a0e7222190a450bea2ed0a33ec4b/USER'
#config.Data.inputDataset = '/DoubleMuon/malcles-DoubleMuon-Run2016D-23Sep2016-v1-AllVerticesCollFiltered-cc71a0e7222190a450bea2ed0a33ec4b/USER'
#config.Data.inputDataset = '/DoubleMuon/malcles-DoubleMuon-Run2016E-23Sep2016-v1-AllVerticesCollFiltered-cc71a0e7222190a450bea2ed0a33ec4b/USER'
#config.Data.inputDataset = '/DoubleMuon/malcles-DoubleMuon-Run2016FPrim-23Sep2016-v1-AllVerticesCollFiltered-cc71a0e7222190a450bea2ed0a33ec4b/USER'
#config.Data.inputDataset = '/DoubleMuon/malcles-DoubleMuon-Run2016G-23Sep2016-v1-AllVerticesCollFiltered-cc71a0e7222190a450bea2ed0a33ec4b/USER'
#config.Data.inputDataset = '/DoubleMuon/malcles-DoubleMuon-Run2016H-PromptReco-v3-a2e3f53203df44f386bf54297150fd64/USER'
config.Data.inputDataset = '/DoubleMuon/malcles-DoubleMuon-Run2016H-PromptReco-v2-a2e3f53203df44f386bf54297150fd64/USER'

config.Data.inputDBS = 'phys03'

config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 1000

config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())

config.Data.publication = False
#config.Data.outputDatasetTag =  'VtxHighMass8020-ZDATA25ns_RunB'
#config.Data.outputDatasetTag =  'VtxHighMass8020-ZDATA25ns_RunCBis'
#config.Data.outputDatasetTag =  'VtxHighMass8020-ZDATA25ns_RunDBis'
#config.Data.outputDatasetTag =  'VtxHighMass8020-ZDATA25ns_RunEBis'
#config.Data.outputDatasetTag =  'VtxHighMass8020-ZDATA25ns_RunFBis'
#config.Data.outputDatasetTag =  'VtxHighMass8020-ZDATA25ns_RunG'
#config.Data.outputDatasetTag =  'VtxHighMass8020-ZDATA25ns_RunH'
config.Data.outputDatasetTag =  'VtxHighMass8020-ZDATA25ns_RunHv2'

config.Site.storageSite = 'T2_FR_GRIF_IRFU'
