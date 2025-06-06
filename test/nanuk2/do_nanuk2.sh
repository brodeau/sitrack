#!/bin/bash

. ./conf.bash

EXE1="${SITRACK_DIR}/tools/generate_idealized_seeding.py"; echo; echo "${EXE1}"
EXE2="python3 -u ${SITRACK_DIR}/si3_part_tracker.py"     ; echo; echo "${EXE2}"

echo $NDATE1
echo


#exit

#ICOARSEN

fout="./nc/sitrack_seeding_nemoTsi3_${NDATE1}_HSS${iHSS}.nc"

if [ ! -f ${fout} ]; then
    CMD="${EXE1} -d ${LDATE1} -m ${FNMM}  -i ${FSI3IN} -v ${NM_ICECONC} -k 0 -f ${FFSM} -S ${iHSS}"
    echo; echo " * Launching:"; echo "    ${CMD}"; echo
    ${CMD}
    echo
fi


# Callin `sitrack` to build trajectories:
CMD="${EXE2} -i ${FSI3IN} -m ${FNMM} -s ${fout} -N ${NEMO_CONF} -c ${NM_ICECONC} -R ${RESKM} -p ${FRQ_PLOT_H} -V siconc-t,sivolu-t"
echo
echo " *** About to launch:"; echo "     ${CMD}"; echo
${CMD}
echo
