#!/bin/bash
 WORKDIR="/afs/cern.ch/user/n/nbower/public/CMU/Stealth_SUSY/ST_plotting/CMSSW_8_0_24/src/STEALTH_PLAY"
 export LSB_JOB_REPORT_MAIL=N
    export X509_USER_PROXY=/afs/cern.ch/user/n/nbower/private/x509up_u94443
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    cd /afs/cern.ch/user/n/nbower/public/CMU/Stealth_SUSY/ST_plotting/CMSSW_8_0_24/src/STEALTH_PLAY/ && eval `scramv1 runtime -sh`
    cd ${WORKDIR}
    python /afs/cern.ch/user/n/nbower/public/CMU/Stealth_SUSY/ST_plotting/CMSSW_8_0_24/src/STEALTH_PLAY/fakephoSel.py -l ${1} -m ${2} -f ${3} -p ${4} 