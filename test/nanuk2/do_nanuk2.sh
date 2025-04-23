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
    CMD="${EXE1} -d ${LDATE1} -m ${FNMM}  -i ${FSI3IN} -v ${NM_ICECONC} -k 0 -f ${FFSM}  -S 2"
    echo; echo " * Launching:"; echo "    ${CMD}"; echo
    ${CMD}
    echo
fi

exit


# Actually that the ice tracker that should look inside the nc file to get date 1 and 2:
CMD="${EXE2} -i ${FSI3IN} -m ${FNMM} -s ${fout}" ; # with nc file for init seed...
echo
echo " *** About to launch:"; echo "     ${CMD}"; echo
${CMD}
