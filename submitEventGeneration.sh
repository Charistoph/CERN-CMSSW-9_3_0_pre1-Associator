#! /bin/bash
#BSUB -q 1nh
#echo "Application started"
cd $WORKDIR
#cd /tmp/cbernkop
cd /afs/cern.ch/user/c/cbernkop/ComparingTracks/CMSSW_9_3_0_pre1/src/
eval `scramv1 runtime -sh`
cd $WORKDIR
echo "1 complete"
#cd /tmp/cbernkop
cmsRun /afs/cern.ch/user/c/cbernkop/ComparingTracks/CMSSW_9_3_0_pre1/src/10001.0_SingleElectronPt10+SingleElectronPt10_pythia8_2017_GenSimFull+DigiFull_2017+RecoFull_2017+ALCAFull_2017+HARVESTFull_2017/SingleElectronPt10_pythia8_cfi_GEN_SIM.py
cmsRun /afs/cern.ch/user/c/cbernkop/ComparingTracks/CMSSW_9_3_0_pre1/src/10001.0_SingleElectronPt10+SingleElectronPt10_pythia8_2017_GenSimFull+DigiFull_2017+RecoFull_2017+ALCAFull_2017+HARVESTFull_2017/step2_DIGI_L1_DIGI2RAW_HLT.py
echo "2 complete"
output="step1_new_`date +%y%m%d%l%M`.root"
cp step1.root /afs/cern.ch/user/c/cbernkop/ComparingTracks/CMSSW_9_3_0_pre1/src/10001.0_SingleElectronPt10+SingleElectronPt10_pythia8_2017_GenSimFull+DigiFull_2017+RecoFull_2017+ALCAFull_2017+HARVESTFull_2017/$output
output="step2_new_`date +%y%m%d%l%M`.root"
cp step2.root /afs/cern.ch/user/c/cbernkop/ComparingTracks/CMSSW_9_3_0_pre1/src/10001.0_SingleElectronPt10+SingleElectronPt10_pythia8_2017_GenSimFull+DigiFull_2017+RecoFull_2017+ALCAFull_2017+HARVESTFull_2017/$output
echo "3 complete"
