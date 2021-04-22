GEO=/Home/eud/hfleisc1/amego_software/ComPair/Geometry/AMEGO_Midex/TradeStudies/Tracker/BasePixelTracker/AmegoBase.geo.setup
#GEO=/Home/eud/hfleisc1/amego_software/ComPair/Geometry/AMEGO_Midex/AmegoBase.geo.setup

SIMDIR=/data/slag2/hfleisc1/amego_sims/simfiles/

INDIR=$SIMDIR/pixel/BG_above_0.01MeV/
CWDIR=${PWD}

ODIR=${CWDIR}/pixel_bg/
mkdir -p ${ODIR}

cd $INDIR

for BG in HadronicBackground HadronicBackgroundDecay TrappedHadronicBackground TrappedHadronicBackgroundDecay LeptonicBackground PhotonicBackground TrappedLeptonicBackground
do
    INFILE=${BG}/${BG}.p.tra.gz
    if [[ ! -e $INFILE ]]
    then
	INFILE=${BG}/${BG}.p1.tra.gz
    fi
    ls $INFILE
    ln -s $INFILE .
    for TYPE in all P TC UC C 
    do
	CFG=${CWDIR}/mimrec_AMEGOX_PlotBkgSpectrum_${TYPE}.cfg
	mimrec -f `basename $INFILE` -g $GEO -c $CFG -s -o ${ODIR}/${BG}.spectrum_${TYPE}.root -n >& ${ODIR}/logSpectrum_${BG}_${TYPE}.txt 
    done 
done

cd $CWDIR
