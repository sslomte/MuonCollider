#!/bin/bash

source /opt/ilcsoft/muonc/init_ilcsoft.sh

git clone https://github.com/MuonColliderSoft/MuonCutil.git

cd MuonCutil/SoftCheck

GEO="/opt/ilcsoft/muonc/detector-simulation/geometries/MuColl_v1/MuColl_v1.xml"

ddsim --compactFile ${GEO} --steeringFile sim_steer.py &> $HOME/sim$$.out

Marlin --InitDD4hep_mod4.DD4hepXMLFile=${GEO} reco_steer.xml &> $HOME/reco$$.out

