#! /bin/bash
#BSUB -q 1nh
#
echo "Application started"
cd $WORKDIR
#cd /tmp/cbernkop
cd /afs/cern.ch/user/c/cbernkop/ComparingTracks/CMSSW_9_3_0_pre1/src/
eval `scramv1 runtime -sh`
cd $WORKDIR
echo "1 complete: cmsenv"
#cd /tmp/cbernkop
cmsRun /afs/cern.ch/user/c/cbernkop/ComparingTracks/CMSSW_9_3_0_pre1/src/10001.0_SingleElectronPt10+SingleElectronPt10_pythia8_2017_GenSimFull+DigiFull_2017+RecoFull_2017+ALCAFull_2017+HARVESTFull_2017/SingleElectronPt10_pythia8_cfi_GEN_SIM.py
echo "2 complete: cmsRun 1"
cmsRun /afs/cern.ch/user/c/cbernkop/ComparingTracks/CMSSW_9_3_0_pre1/src/10001.0_SingleElectronPt10+SingleElectronPt10_pythia8_2017_GenSimFull+DigiFull_2017+RecoFull_2017+ALCAFull_2017+HARVESTFull_2017/step2_DIGI_L1_DIGI2RAW_HLT.py
echo "3 complete: cmsRun 2"
#output=`date +%y%m%d%l%M`
#mkdir /afs/cern.ch/user/c/cbernkop/ComparingTracks/CMSSW_9_3_0_pre1/src/10001.0_SingleElectronPt10+SingleElectronPt10_pythia8_2017_GenSimFull+DigiFull_2017+RecoFull_2017+ALCAFull_2017+HARVESTFull_2017/$output
#mkdir /afs/cern.ch/work/c/cbernkop/condor_output/1
echo "4 complete: mkdir"
#output="step1_new_`date +%y%m%d%l%M`.root"
#cp step1.root /afs/cern.ch/user/c/cbernkop/ComparingTracks/CMSSW_9_3_0_pre1/src/10001.0_SingleElectronPt10+SingleElectronPt10_pythia8_2017_GenSimFull+DigiFull_2017+RecoFull_2017+ALCAFull_2017+HARVESTFull_2017/$output/step1.root
cp step1.root /afs/cern.ch/work/c/cbernkop/condor_output/1/step1.root
echo "5 complete: cp 1"
#output="step2_new_`date +%y%m%d%l%M`.root"
cp step2.root /afs/cern.ch/work/c/cbernkop/condor_output/1/step2.root
echo "6 complete: cp 6"

