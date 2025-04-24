#!/bin/bash

SITRACK_DIR="${HOME}/DEV/sitrack"


# SI3/NANUQ output file for work with
NEMO_CONF="NANUK2" ; NZ=75
YEAR="2016"
DATE1="${YEAR}0101"; DATE2="${YEAR}0313"
NEMO_EXP="00"
SBDIR="00000001-00003504_dev_BBM" 
RESKM="24km"
NM_ICECONC="siconc-t"

# sitrack:
iHSS=2 
DT_BIN="72"
FRQ_PLOT_H=12

# Where is the data saved depending on host:
host=`hostname | cut -d '.' -f2`
case ${host} in
    "frazilo")
        DATA_DIR="/data/laurent"
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac

#---------------------------------------------------------------------------


FSI3IN="${NEMO_CONF}-${NEMO_EXP}_1h_${DATE1}_${DATE2}_icemod.nc4"

FSI3IN="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}-${NEMO_EXP}-S/${SBDIR}/${FSI3IN}"

FNMM="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}.L${NZ}-I/mesh_mask_${NEMO_CONF}_L${NZ}_4.2.2.nc"

FFSM="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}.L${NZ}-I/mask_${NEMO_CONF}_rgps.nc"

echo
echo "FSI3IN => ${FSI3IN}"
echo "FNMM => ${FNMM}"
echo "FFSM => ${FFSM}"
echo

mkdir -p ./figs ./npz

YYYY=`echo ${DATE1} | cut -c1-4`
MM=`echo ${DATE1} | cut -c5-6`
DD=`echo ${DATE1} | cut -c7-8`
NDATE1="${YYYY}-${MM}-${DD}"
LDATE1="${NDATE1}_00:00:00"

