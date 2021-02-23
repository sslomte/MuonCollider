# MadGraph5_aMC@NLO and DELPHES
## Monte Carlo simulation with MadGraph5
### MG5 installation
```bash
wget https://cms-project-generators.web.cern.ch/cms-project-generators/MG5_aMC_v2.7.2.tar.gz
tar xvf MG5_aMC_v2.7.2.tar.gz
```
### Install pythia8 on MG5
```bash
cd MG5_aMC_v2_7_2
./bin/mg5_aMC
```

```MG5
install pythia8
exit
```
### Simulation
```bash
pushd /cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/cmssw/CMSSW_11_1_6; cmsenv;popd
cd MG5_aMC_v2_7_2
./bin/mg5_aMC
```
Here is a example for generating the signal for diHiggs channel, similar for background or other channel 
```MG5
import model sm
generate mu+ mu- > vm vm~ h h
output mucol-hpair_sig
launch
```
Then you will see a menu for MG5, type 1 to open up the pythia8 for format transformation to .hepmc
```MG5
set mh 125.0
set wh 0.004
set ebeam1 1500.
set ebeam2 1500.
set nevents 1000
set iseed 1823211
```
The process will take about 5-60 min, suggest working on vncserver or add following commnad to ssh config:
```config
    ServerAliveInterval 240
    ServerAliveCountMax 2
```
Don't forget to gunzip the .hepmc.tar.gz for preparation for delphes simulation:
```bash 
cd <your output directory>/Events/run_01/
gunzip tag1_pythia8.hepmc.tar.gz
```
## Detector simulation on Delphes
### Installation of Delphes
```bash
git clone git://github.com/delphes/delphes.git Delphes
cd Delphes

