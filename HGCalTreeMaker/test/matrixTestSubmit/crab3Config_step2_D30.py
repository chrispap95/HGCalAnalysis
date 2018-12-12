from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = ''
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'step2_DIGI_L1_L1TrackTrigger_DIGI2RAW_HLT.py'
config.JobType.allowUndistributedCMSSW = False
#config.JobType.outputFiles=['step2.root']

config.section_("Data")
# MC example
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader/'
config.Data.inputDataset = '/TTbar_14TeV/bcaraway-CMSSW_10_4_0_pre2_D30_Step1_v3-bf22a9fcb0c0b8ca16c87498a27c5fad/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#KH config.Data.totalUnits = 100
# MC example ends
# Data example
#config.Data.inputDataset = '/ZeroBias/Run2018A-v1/RAW'
#config.Data.inputDataset = '/HLTPhysics/Run2018A-v1/RAW'
#config.Data.runRange = '315361-315690'
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 10
# Data example ends

config.Data.publication = True
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/' # Parameter Data.publishDbsUrl has been renamed to Data.publishDBS
config.Data.outputDatasetTag = 'CMSSW_10_4_0_pre2_D30_Step2_v4' # <== Check!!!

config.Data.outLFNDirBase = '/store/user/bcaraway/crab_outputs'  # Data.outLFN has been renamed to Data.outLFNDirBase
config.Data.ignoreLocality = True

config.section_("Site")
config.Site.storageSite = 'T3_US_Baylor'
#KH (this whitelisting below is not really necessary. we can use any T2/T3 for running jobs. we can still send output to Baylor)
#KH config.Site.whitelist = ['T3_US_Baylor']
config.Site.whitelist = ['T3_US_Baylor','T2_US_*']
