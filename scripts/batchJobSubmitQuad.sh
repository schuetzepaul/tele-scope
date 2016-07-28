#!/bin/sh

# Default runs file
INPUT_FILE='runlist-quad.dat'

USAGE="Usage: ./`basename $0` [-l evts] [-r arg] [-hv] args"

# Define useage function
usage(){
cat <<EOF

Usage: `basename $0` -l 10000 -c 2076 -r 1440-1450,2590

Purpose of script is to submit jobs to the batch system

OPTIONS:
   -h      Show this message
   -l      Number of events
   -c      Run for conversion file
   -r      Range for run numbers
   -v      Verstion of script
EOF
}

# Parse command line options.
while getopts l:c:r:h:v: OPT; do
    case "$OPT" in
        h)
            echo $USAGE
            exit 0
            ;;
        v)
            echo "`basename $0` version 0.1"
            exit 0
            ;;
        l)
            EVENTS=$OPTARG
            ;;
        c)
            CONVRUN=$OPTARG
            ;;
        r)
	    if [ ! -f $OPTARG ]; then
		STRING=`echo $OPTARG | awk -F "," '{for(i=1; i <= NF; i++) printf "%s ",$i}'`

	    elif [ -s $OPTARG ]; then
		STRING=`cat $OPTARG | awk -F "," '{for(i=1; i <= NF; i++) printf "%s ",$i}'`
	    fi

	    array=(${STRING})

	    for elem in "${array[@]}"
	    do
		if [[ $elem =~ "-" ]]; then
		    element=`echo $elem | awk -F "-" '{for(i=$1; i <= $2; i++) printf "%s ",i}'`
		    RUNS+=($element)
		else
		    RUNS+=($elem)
		fi
	    done
            ;;
        \?)
            # getopts issues an error message
            echo $USAGE >&2
            exit 1
            ;;
    esac
done

# We want at least one non-option argument.
# Remove this block if you don't need it.
#if [ $# -eq 0 ]; then
#    usage  >&2
#    exit 1
#fi

# Set environmental variables to corresponding script path
source $(dirname `readlink -f $0`)/parallelTools.sh

SCRIPT=$(readlink -f $0)
SCRIPTPATH=$(dirname "$SCRIPT")

LIST=""
RUN=""

# Setup regex expressions
regex_run="run ([0-9]*)"
regex_list="list (runlist-quad.dat)"

# Check for conversion run
if [ -z ${CONVRUN+x} ]; then
    regex_convrun=""
else
    regex_convrun="-c ${CONVRUN}"
    echo "Get conversion factors from run ${CONVRUN}"
fi

# Dispatch the jobs to the NAF batch system
for RUN in "${RUNS[@]}"; do

#    echo "Starting Run ${RUN}"
    $QA -N quad-${RUN} -l $EVENTS $regex_convrun $RUN &

done


wait

if [ "$isNAF" = 1 ]; then
    echo "Please check your jobs with qstat -u $USER | grep $MODE"
else
    echo "Processing all nominal samples finished!"
fi
