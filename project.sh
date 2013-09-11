#/bin/sh
set -x
SGL_SPK_DIR="/Users/harsha/GRBAnalysis/SingleSpikeData/" 
FULLGRB_DIR="/Users/harsha/GRBAnalysis/FullGRBData/MisalignedJet/Jitter/UniversalGamma/LC2/Sharp/Threshold/Case3/"

mkdir -p $FULLGRB_DIR


EMISSION=Jitter
#inputs are Synchrotron,Formula,Jitter

GAMMA_PROF=Universal
#inputs are TopHat, Lorentzian, Universal,Magnetic(which is just tophat, but with diff comoving spectrum)

LT_CRV_NUM=2

LT_CRV_TYPE=sharp 
#'sharp' or 'thick' are inputs, do not put in quotes.

DILATION_FAC=1.0  
#Dilating the initial single spike to be longer/shorter, inputs are 1.0,2.0,3.0,10.0,20.0,30.0,100.0,1000.0

VWG_ANG=2  
#1 - head on, 2,3,4 varying degrees of misalignment

WIN_LEFT=5e2 
#energy window of fit - left bound

WIN_RIGHT=1e5 
#energy window of fit - right bound

THRSHLD_FRAC=1e-3 
# Threshold value - 1e-1,1e-2,1e-3

ERR_FRAC=0.001; 
# Fraction of Ef(E) included as error in fitting program

#--------------------------------------------------------------------------------------------------------
#DO NOT MODIFY BELOW THIS POINT

LT_CRV_FILE=FluxAmplitudes_${LT_CRV_NUM}_${LT_CRV_TYPE}.txt
N_SPIKES=`cat ${LT_CRV_FILE} | wc -l`

SGL_SPK_FILE=${SGL_SPK_DIR}Angle${VWG_ANG}_onespike_${EMISSION}_${GAMMA_PROF}.txt

FULL_GRB_FILE=${FULLGRB_DIR}FakeGRB_${LT_CRV_TYPE}.txt

FULL_BOLFLUX_FILE=${FULLGRB_DIR}BolFlux_${LT_CRV_TYPE}.txt

BIN_BOL_FILE=${FULLGRB_DIR}BinBolFlux_${LT_CRV_TYPE}_${THRSHLD_FRAC}.txt

BIN_GRB_FILE=${FULLGRB_DIR}BinGRB_${LT_CRV_TYPE}_${THRSHLD_FRAC}.txt

BIN_FIT_TABLE=${FULLGRB_DIR}BinFit_${LT_CRV_TYPE}_${THRSHLD_FRAC}.txt


#Execution

#./FakeGRB $LT_CRV_FILE $N_SPIKES $DILATION_FAC $SGL_SPK_FILE $FULL_GRB_FILE 

./BinData $FULL_GRB_FILE $WIN_LEFT $WIN_RIGHT $FULL_BOLFLUX_FILE $THRSHLD_FRAC $BIN_BOL_FILE $BIN_GRB_FILE

./binmodelfit $BIN_GRB_FILE $BIN_BOL_FILE $BIN_FIT_TABLE $WIN_LEFT $WIN_RIGHT $ERR_FRAC

#./PlotScripts/./plot_A.sm ${BIN_FIT_TABLE} 0 ${WRK_DIR}/FEpvsa.ps 
#./PlotScripts/./plot_A.sm ${BIN_FIT_TABLE} 0 ${WRK_DIR}/ScatterPlot_Magnetic.ps
#./PlotScripts/./plot_B.sm ${BIN_FIT_TABLE} 0 ${WRK_DIR}/Histplots_Magnetic.ps
#./PlotScripts/./plot_C.sm ${BIN_FIT_TABLE} 20

#./PlotScripts/./plot_A.sm ${BIN_FIT_TABLE} 10 
#./PlotScripts/./plot_A_hist.sm ${BIN_FIT_TABLE} 60 
#./PlotScripts/./alpha_hist.sm ${BIN_FIT_TABLE} 60 
#./PlotScripts/./plot_B.sm ${BIN_FIT_TABLE} 10 
#./PlotScripts/./plot_C.sm ${BIN_FIT_TABLE} 10
