#!/bin/sh

SCRIPT=$(readlink -f $0)
SCRIPTPATH=$(dirname "$SCRIPT")

if [ `hostname | grep "nafhh"` ]; then
    echo "Running on the NAF - so let's submit our jobs, job output will be stored in batch_output/..."
    mkdir -p "batch_output"
    
    #define wait function - on the NAF, we just qsub. So let's wait until jobs have been submitted
    w() {
        wait
    }
    
    isNAF=1
    
    if grep -q "slc6" <<< "$SCRAM_ARCH"; then
        echo "Running at SL6 scram architecture. Submitting jobs to SL6 nodes."
        TA="qsub -l distro=sld6 -@ $SCRIPTPATH/qsubParams.txt $SCRIPTPATH/../tele"
	SA="qsub -l distro=sld6 -@ $SCRIPTPATH/qsubParams.txt $SCRIPTPATH/../scope"
	QA="qsub -l distro=sld6 -@ $SCRIPTPATH/qsubParams.txt $SCRIPTPATH/../quad"
    else
        TA="qsub -@ $SCRIPTPATH/qsubParams.txt $SCRIPTPATH/../tele"
	SA="qsub -@ $SCRIPTPATH/qsubParams.txt $SCRIPTPATH/../scope"
	QA="qsub -@ $SCRIPTPATH/qsubParams.txt $SCRIPTPATH/../quad"
    fi
else
    w() {
        while [ `ps ax | grep -E 'tele|scope|quad' | wc -l` -gt 10 ]; do
            sleep 1;
        done
    }

    isNAF=0
    TA=tele
    SA=scope
    QA=quad
fi
