#!/bin/bash

YEAR="2016"

SITRACK_DIR="${HOME}/DEV/sitrack"

NEMO_CONF="NANUK2" ; NZ=75

cxtra=""

export DATE1="${YEAR}0101"

# /data/gcm_setup/NANUK4/BBM00$ NANUK4_ICE-BBM00_1h_19970101_19970331_icemod.nc4
NEMO_EXP="00"  ; SBDIR="00000001-00003504_dev_BBM" ; export DATE2="${YEAR}0313"

export iHSS=2 ; RESKM=24

NJPAR=1 ; # number of jobs we can launch in //

FREQ_AN_DAYS=3 ; # frequency in days of the deformation analysis...

xtra_sfx=""

LCOARSEN="0"

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        #/MEDIA/data/NANUK4/BBM00/NANUK4_ICE-BBM00_1h_19970101_19970331_icemod_LIGHT480.nc4
        export DATA_DIR="/MEDIA/data"
        export iHSS=1 ; RESKM=12
        LCOARSEN="80"
        #        
        Ndays=6
        xtra_sfx="_LIGHT480"
        #
        ;;
    "mcp-oceannext-01")
        export DATA_DIR="/data/gcm_setup"
        #
        iHSS=8
        NJPAR=4       
        #
        ;;
    "frazilo")
        export DATA_DIR="/data/laurent"
        Ndays=31
        #
        NJPAR=30
        #NJPAR=1
        #
        #SI3IN="${NEMO_CONF}_ICE-${NEMO_EXP}_1h_${YEAR}0101_${YEAR}0205_icemod.nc4" ; # 1 month !!!
        #
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac


FSI3IN="${NEMO_CONF}-${NEMO_EXP}_1h_${DATE1}_${DATE2}_icemod${xtra_sfx}.nc4"

export FSI3IN="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}-${NEMO_EXP}-S/${SBDIR}${cxtra}/${FSI3IN}"

export FNMM="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}.L${NZ}-I/mesh_mask_${NEMO_CONF}_L${NZ}_4.2.2.nc"

export FFSM="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}.L${NZ}-I/mask_${NEMO_CONF}_rgps.nc"

echo
echo "FSI3IN => ${FSI3IN}"
echo "FNMM => ${FNMM}"
echo "FFSM => ${FFSM}"
echo

mkdir -p ./figs ./npz

YYYY=`echo ${DATE1} | cut -c1-4`
MM=`echo ${DATE1} | cut -c5-6`
DD=`echo ${DATE1} | cut -c7-8`
export NDATE1="${YYYY}-${MM}-${DD}"
export LDATE1="${NDATE1}_00:00:00"

