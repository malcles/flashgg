# NB this command is specific to the configuration on lxplus and is not gaurenteed elsewhere
queue="2nd"
useAAA=0
version="xgboostDiphoBDT2017"
fggRunJobs.py --load tth_sig_jobs_2017.json -d tth_sig_jobs_$version -x cmsRun workspaceStd.py maxEvents=-1 -n 500 -q $queue -D -P useAAA=$useAAA doHTXS=False doStage1=True doFiducial=False tthTagsOnly=False
