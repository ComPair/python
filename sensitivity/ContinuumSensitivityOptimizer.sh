#! /bin/bash

#Author: Carolyn Kierans (mainly), Henrike Fleischhack

help() {
  echo ""
  echo "Automated script to paralellize SensitivityOptimizer";
  echo "";
  echo "";
}

CMD=( "$@" )

THREADS=1
NICELEVEL=0

#Got that one from Carolyn
CFG="TotalBackground.extracted_mimrec_open.cfg"

#Set paths to continuum source and background sims (and Geometry file) here.
SourceFILE_PATH="/data/slag2/hfleisc1/amego_sims/simfiles/base/background/ContinuumPointSource/"
SOURCEFILES=("FarFieldPointSource_Continuum_0.080MeV_Cos1.0.p.tra.gz" "FarFieldPointSource_Continuum_0.800MeV_Cos1.0.p.tra.gz" "FarFieldPointSource_Continuum_8.000MeV_Cos1.0.p.tra.gz" )

BkgFILE="/data/slag2/hfleisc1/amego_sims/simfiles/base/background/FullBackground.tra.gz"

GEOMETRY="/Home/eud/hfleisc1/amego_software/ComPair/Geometry/AMEGO_Midex/AmegoBase.geo.setup"

#set output directory here
ODIR="$PWD/Sensitivity/base/"
mkdir -p $ODIR

#Energy bins, in keV. Picked such that width of the bin = bin center
ENERGY_Low=( 180 250  350  500  650 1000 1750 2500 3500 5000)
ENERGY_High=(530 750 1050 1500 1950 3000 5250 7500 10500 15000)

# Find the default number of threads
if [[ ${OSTYPE} == darwin* ]]; then
  THREADS=`sysctl -n hw.logicalcpu_max`
elif [[ ${OSTYPE} == linux* ]]; then
  THREADS=`grep processor /proc/cpuinfo | wc -l`
fi

# Check if revan exists
if (`test -f ${MEGAlib}/bin/SensitivityOptimizer`); then
  echo " "
  echo "ERROR: The SensitivityOptimizer executable does not exist. Try to (re-)compile MEGAlib."
  echo " "
  exit 1;
fi


echo " "
echo "Launching mSensitivityOptimizer"
echo " "
echo "Number of threads to use:  ${THREADS}"
echo "Nice level:                ${NICELEVEL}"
if [[ ${GEOMETRY} != "" ]]; then
  echo "Geometry:                ${GEOMETRY}"
fi
echo "Revan configuration file:  ${CFG}"
echo "Source files:             ${SourceFILE_1}"
echo "Bkg file: 		${BkgFILE}"
echo "Energies: 		${ENERGY_Low[@]}"

# Now run

#We produce "signal" sims in three different energy ranges to get good statistics
#(1e8 events total per energy range).
#Set which file to use, energy range and number of simulated signal events here.
for i in ${!ENERGY_Low[@]}; do
  SOURCE="${SourceFILE_PATH}"
  if (( "${ENERGY_Low[$i]}" < 800)); then
    SOURCE+="${SOURCEFILES[0]}"
    PLOW=80
    PHIGH=10000
    EVENTS=100000000
  elif (( "${ENERGY_Low[$i]}" < 8000)); then
    SOURCE+="${SOURCEFILES[1]}"
    PLOW=800
    PHIGH=100000
    EVENTS=100000000
  else
    SOURCE+="${SOURCEFILES[2]}"
    PLOW=8000
    PHIGH=1000000
    EVENTS=100000000
  fi

  if (( "${ENERGY_Low[$i]}" < 1750 )); then
      NARM=8
  else
      NARM=12
  fi

#Adjust the desired observing time (-t 94610000) and simulated bg observing time (-b ... 7200) below.

  echo "Launching SensitivityOptimizer for ${ENERGY_Low[$i]}-${ENERGY_High[$i]}"
  CMD="source ${MEGALIB}/bin/source-megalib.sh; SensitivityOptimizer -n $ODIR/Sensitivity_Continuum_"
  CMD+="${ENERGY_Low[$i]}-"
  CMD+="${ENERGY_High[$i]}_"
  CMD+="${PLOW}"
  CMD+="keV_R1_bash -t 94610000"
  CMD+=" -g ${GEOMETRY}"
  CMD+=" -k ${SOURCE} ${EVENTS} 70685.8 2 ${PLOW} ${PHIGH} 0 0"
  CMD+=" -b ${BkgFILE} 7200"
  CMD+=" -c ${CFG}"
  CMD+=" --contegy ${ENERGY_Low[$i]} ${ENERGY_High[$i]}"
  CMD+=" --ptheta 0 0 1 --pphi 0 0 1 --csl 2 2 10 --tsl 1 2 20 --ehc 90 90 1 --arm 1 8 $NARM --phi 50 180 7"
  CMD+=" --cqf 100 100 1 --fdi 1 4 3 --spd 20 80 7 --pop 10 30 3 --idp 500 700 3"
#  CMD+=" --ptheta 0 0 1 --pphi 0 0 1 --csl 2 2 10 --tsl 1 2 20 --ehc 90 90 1 --arm 1 15 29 --phi 0 180 37"
#  CMD+=" --cqf 100 100 1 --fdi 1 9 19 --spd 10 100 10 --pop 10 60 11 --idp 100 1000 10"
# comment in the last two lines (and comment out the lines before) to open up the cut space to be optimized

  nohup nice -n ${NICELEVEL} bash -c "${CMD}" > /dev/null &
  echo ${CMD}
  sleep 0.5
done


# We always wait until all runs have finished
echo "Waiting till all runs have finished..."
wait


exit 0;
