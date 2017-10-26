from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

#config.General.requestName = 'VtxHighMass928-ZDATA25ns_Run2017Cv1'
#config.General.requestName = 'VtxHighMass928-ZDATA25ns_Run2017Cv2'
#config.General.requestName = 'VtxHighMass928-ZDATA25ns_Run2017Bv1'
config.General.requestName = 'VtxHighMass928-ZDATA25ns_Run2017Bv2'





config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True 

config.JobType.psetName = 'Validation/python/outputAnalyzerHighMassDataChangeHLT.py'

#config.Data.inputDataset = '/DoubleMuon/malcles-DoubleMuon-Run2017C-PromptReco-NewMini-v1-db2c70837ea49ff752a17c3954ca82e7/USER'
#config.Data.inputDataset = '/DoubleMuon/malcles-DoubleMuon-Run2017C-PromptReco-NewMini-v2-db2c70837ea49ff752a17c3954ca82e7/USER'
#config.Data.inputDataset = '/DoubleMuon/malcles-DoubleMuon-Run2017B-PromptReco-NewMini-v1-db2c70837ea49ff752a17c3954ca82e7/USER'
config.Data.inputDataset = '/DoubleMuon/malcles-DoubleMuon-Run2017B-PromptReco-NewMini-v2-db2c70837ea49ff752a17c3954ca82e7/USER'


config.Data.inputDBS = 'phys03'

config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 1000

config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())

config.Data.publication = False

#config.Data.outputDatasetTag = 'VtxHighMass928-ZDATA25ns_Run2017Cv1'
#config.Data.outputDatasetTag = 'VtxHighMass928-ZDATA25ns_Run2017Cv2'
#config.Data.outputDatasetTag = 'VtxHighMass928-ZDATA25ns_Run2017Bv1'
config.Data.outputDatasetTag = 'VtxHighMass928-ZDATA25ns_Run2017Bv2'



config.Site.storageSite = 'T2_FR_GRIF_IRFU'
