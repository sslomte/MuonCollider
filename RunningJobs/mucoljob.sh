#!/bin/bash

mkdir jobHome
mv ilcsoft.sh jobHome

singularity exec -B jobHome:$HOME  /cvmfs/muoncoll.infn.it/sw/singularity/MuonColl_v02-05-MC.sif /bin/bash $HOME/ilcsoft.sh $*


mv jobHome/*.out .

### clean dir

rm -rf jobHome

