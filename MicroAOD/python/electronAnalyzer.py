import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
from RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi import *
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *


process = cms.Process("FLASHggMicroAOD")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'PLS170_V7AN1::All'
#process.GlobalTag.globaltag = 'auto:run2_mc'
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

# Fix because auto:run2_mc points to MCRUN2_74_V9::All
current_gt = process.GlobalTag.globaltag.value()
if current_gt.count("::All"):
    new_gt = current_gt.replace("::All","")
    print 'Removing "::All" from GlobalTag by hand for condDBv2: was %s, now %s' % (current_gt,new_gt)
    process.GlobalTag.globaltag = new_gt


# 2012 data
#process.GlobalTag = GlobalTag(process.GlobalTag, 'GR_R_74_V8A::All')
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
#        "/store/relval/CMSSW_7_4_0_pre9/DoubleElectron/MINIAOD/GR_R_74_V8A_RelVal_zEl2012D-v1/00000/5A04EF0A-29D4-E411-BB12-003048FFCC2C.root",
#        "/store/relval/CMSSW_7_4_0_pre9/DoubleElectron/MINIAOD/GR_R_74_V8A_RelVal_zEl2012D-v1/00000/A6D4F50A-29D4-E411-98A2-002618943838.root",
#        "/store/relval/CMSSW_7_4_0_pre9/DoubleElectron/MINIAOD/GR_R_74_V8A_RelVal_zEl2012D-v1/00000/DE947AB8-FED3-E411-995F-0025905AA9F0.root"
#        ))

# PHYS14 Files
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
#"/store/mc/Phys14DR/DYToMuMu_M-50_Tune4C_13TeV-pythia8/MINIAODSIM/PU40bx25_tsg_castor_PHYS14_25_V1-v2/00000/622CAFBA-BD9A-E411-BE11-002481E14FFC.root",
#"/store/mc/Phys14DR/DYToMuMu_M-50_Tune4C_13TeV-pythia8/MINIAODSIM/PU40bx25_tsg_castor_PHYS14_25_V1-v2/00000/FA4B46B9-8E9A-E411-A899-002590A3C954.root",
#"/store/mc/Phys14DR/DYToMuMu_M-50_Tune4C_13TeV-pythia8/MINIAODSIM/PU40bx25_tsg_castor_PHYS14_25_V1-v2/10000/8607F88E-F799-E411-A180-0025B3E063F0.root",
#"/store/mc/Phys14DR/DYToMuMu_M-50_Tune4C_13TeV-pythia8/MINIAODSIM/PU40bx25_tsg_castor_PHYS14_25_V1-v2/10000/F620A7C9-F799-E411-8DEF-002590A371AC.root"
#))
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/mc/Phys14DR/GluGluToHToGG_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/3C2EFAB1-B16F-E411-AB34-7845C4FC39FB.root"))
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/mc/Phys14DR/VBF_HToGG_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/4A8E0BD1-026C-E411-8760-00266CFFA418.root"))
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/mc/Phys14DR/WH_ZH_HToGG_M-125_13TeV_pythia6/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/24B70163-5769-E411-93CA-002590200A28.root"))
### process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/mc/Phys14DR/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/101611CC-026E-E411-B8D7-00266CFFBF88.root"))
## process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/mc/Phys14DR/GJets_HT-100to200_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00D67F78-2873-E411-B3BB-0025907DC9C0.root"))

# 740 RelVal
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/relval/CMSSW_7_4_0_pre9_ROOT6/RelValH130GGgluonfusion_13/MINIAODSIM/MCRUN2_74_V7-v1/00000/0A35F6D-DAD1-E411-A8CC-0026189438CC.root"))
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/relval/CMSSW_7_4_0_pre9/RelValH130GGgluonfusion_13/MINIAODSIM/PU25ns_MCRUN2_74_V7-v1/00000/5ABC049C-4CD4-E411-B28A-0025905A613C.root",
#                                                                         "/store/relval/CMSSW_7_4_0_pre9/RelValH130GGgluonfusion_13/MINIAODSIM/PU25ns_MCRUN2_74_V7-v1/00000/C65FAFAA-4CD4-E411-9026-0025905A607E.root"))

# Spring15
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/mc/RunIISpring15DR74/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/70000/0232BC3C-01FF-E411-8779-0025907B4FC2.root"))

#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/mc/RunIISpring15DR74/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/4E503483-EA32-E511-AB07-02163E013542.root"))
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/mc/RunIISpring15DR74/ttHJetToGG_M120_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/60000/32B6B3FC-320C-E511-8DCD-02163E00F3F2.root"))

process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
#'/store/mc/RunIISpring15DR74/VHToGG_M130_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/20000/B6D836A2-2135-E511-A2B0-02163E01413E.root',
'/store/mc/RunIISpring15DR74/ZH_HToGG_ZToAll_M120_13TeV_powheg_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/40000/D67B907E-7B32-E511-9A8F-02163E013576.root',
'/store/mc/RunIISpring15DR74/ZH_HToGG_ZToAll_M120_13TeV_powheg_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/40000/DAD86C5C-4B33-E511-8DA2-02163E00EA96.root',
'/store/mc/RunIISpring15DR74/ZH_HToGG_ZToAll_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/026DD591-E034-E511-8C26-02163E011DFF.root',
'/store/mc/RunIISpring15DR74/ZH_HToGG_ZToAll_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/70EE122A-E434-E511-8705-02163E011DFF.root',
'/store/mc/RunIISpring15DR74/ZH_HToGG_ZToAll_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/84675C89-E034-E511-8623-02163E0146D1.root',
'/store/mc/RunIISpring15DR74/ZH_HToGG_ZToAll_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/C2222F8B-E034-E511-879F-02163E0134CC.root',
'/store/mc/RunIISpring15DR74/ZH_HToGG_ZToAll_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/80000/FCFE5523-A635-E511-B5A6-02163E0137FD.root'
#'/store/mc/RunIISpring15MiniAODv2/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/4CB72B87-C26D-E511-A19A-002590D9D9F0.root',
#'/store/mc/RunIISpring15MiniAODv2/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/00F6D2A1-C26D-E511-BE07-00259029E84C.root',
#'/store/mc/RunIISpring15MiniAODv2/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/0ADFE98D-C26D-E511-87F7-002590D9D98E.root',
#'/store/mc/RunIISpring15MiniAODv2/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/6A5DD488-C26D-E511-9708-002590D9D984.root',
#'/store/mc/RunIISpring15MiniAODv2/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/6C4EEC22-C36D-E511-8CCF-002590AC4C74.root',
#'/store/mc/RunIISpring15MiniAODv2/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/9CBC9A81-C26D-E511-8E89-001E67A3FC1D.root',
#'/store/mc/RunIISpring15MiniAODv2/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/CCDF67C0-C26D-E511-B74C-002590791D60.root',
#'/store/mc/RunIISpring15MiniAODv2/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/FC281F86-C26D-E511-9981-001E67397D73.root',
#'/store/mc/RunIISpring15MiniAODv2/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/FE4C1A28-406E-E511-8B11-00259019A41E.root'
#'/store/mc/RunIISpring15DR74/ZH_HToGG_ZToAll_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/06D038C8-0037-E511-8506-02163E0141A3.root',
#'/store/mc/RunIISpring15DR74/ZH_HToGG_ZToAll_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/00E9A77-0037-E511-A086-02163E013645.root',
#'/store/mc/RunIISpring15DR74/ZH_HToGG_ZToAll_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/A316AF3-FF36-E511-80C5-02163E0141A3.root',
#'/store/mc/RunIISpring15DR74/ZH_HToGG_ZToAll_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/EF3A5DB-FF36-E511-B809-02163E011B74.root',
#'/store/mc/RunIISpring15DR74/ZH_HToGG_ZToAll_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/42AC077-0037-E511-AB7D-02163E014225.root',
#'/store/mc/RunIISpring15DR74/ZH_HToGG_ZToAll_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/8BA068E-4938-E511-90C7-02163E013645.root'
#'/store/mc/RunIISpring15DR74/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/40000/38AC53EB-9235-E511-894A-0CC47A13CBEA.root',
#'/store/mc/RunIISpring15DR74/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/40000/50BA4696-8335-E511-B6FE-00259073E4E8.root',
#'/store/mc/RunIISpring15DR74/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/40000/546AAFD8-8D35-E511-BE78-002590D9D84A.root',
#'/store/mc/RunIISpring15DR74/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/40000/5830BB3D-8A35-E511-A74B-842B2B2B0CFE.root',
#'/store/mc/RunIISpring15DR74/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/40000/5ABB4006-8235-E511-B127-0026182FD7A3.root',
#'/store/mc/RunIISpring15DR74/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/40000/8AA7C13B-F035-E511-A392-001F2908BE52.root',
#'/store/mc/RunIISpring15DR74/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/40000/8E98D947-8C35-E511-8303-008CFA111354.root',
#'/store/mc/RunIISpring15DR74/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/40000/9A99CA89-F135-E511-B262-001E673969D2.root',
#'/store/mc/RunIISpring15DR74/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/40000/9CF606AA-EF35-E511-AC7C-001EC9B22258.root',
#'/store/mc/RunIISpring15DR74/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/40000/A4FA76A0-8735-E511-B261-001517FB21BC.root'
))
#process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring("file:myMicroAODOutputFile.root"))
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("store/mc/RunIISpring15DR74/VHToGG_M130_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/20000/B6D836A2-2135-E511-A2B0-02163E01413E.root"))

process.MessageLogger.cerr.threshold = 'ERROR' # can't get suppressWarning to work: disable all warnings for now
# process.MessageLogger.suppressWarning.extend(['SimpleMemoryCheck','MemoryCheck']) # this would have been better...

# Uncomment the following if you notice you have a memory leak
# This is a lightweight tool to digg further
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#                                        ignoreTotal = cms.untracked.int32(1),
#                                        monitorPssAndPrivate = cms.untracked.bool(True)
#


# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 

dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
		 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff']


#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.ElectronAnalyzer = cms.EDAnalyzer('ElectronValidator',
		verbose = cms.untracked.bool(False),
		electronTag = cms.InputTag('slimmedElectrons'),
		vertexTag = cms.InputTag('offlineSlimmedPrimaryVertices'),
		convTag   = cms.InputTag('reducedEgamma','reducedConversions'),
		beamSpotTag = cms.InputTag( "offlineBeamSpot" ),
		genTag = cms.InputTag('packedGenParticles'),
		reducedEBRecHitCollection = cms.InputTag('reducedEcalRecHitsEB'),
		reducedEERecHitCollection = cms.InputTag('reducedEcalRecHitsEE'),
		rhoFixedGridCollection = cms.InputTag('fixedGridRhoAll'),
		mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),

  		eleMVAMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                eleMVATightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),	     
   		eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
                eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
                eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
                eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51"),
		effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt")

)

process.TFileService = cms.Service("TFileService", fileName = cms.string("histo.root") )
#process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('myElectronOutputFile.root'))

#process.p = cms.Path(process.electronMVAValueMapProducer*process.ElectronAnalyzer)
process.p = cms.Path(process.egmGsfElectronIDSequence*process.ElectronAnalyzer)
#process.e = cms.EndPath(process.out)
