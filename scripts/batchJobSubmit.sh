#!/bin/sh

# Default runs file
INPUT_FILE='runs.dat'

USAGE="Usage: ./`basename $0` [-m arg required] [-o arg] [-f arg] [-r arg] [-hv] args"

# Define useage function
usage(){
cat <<EOF

Usage: `basename $0` -m tele -f runs_DATE.dat -r 1440-1450,2590

Purpose of script is to submit jobs to the batch system

OPTIONS:
   -h      Show this message
   -m      Mode 'tele' or 'scope' to run on (REQUIRED)
   -f      Input runs_DATE.dat file
   -r      Range for run numbers
   -o      Output file
   -v      Verstion of script
EOF
}

# Parse command line options.
while getopts hvo:f:r:m: OPT; do
    case "$OPT" in
        h)
            echo $USAGE
            exit 0
            ;;
        v)
            echo "`basename $0` version 0.1"
            exit 0
            ;;
        o)
            OUTPUT_FILE=$OPTARG
            ;;
        f)
            INPUT_FILE=$OPTARG
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
	m)
	    MODE=$OPTARG
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
if [ $# -eq 0 ]; then
    usage  >&2
    exit 1
fi

# Set environmental variables to corresponding script path
source $(dirname `readlink -f $0`)/parallelTools.sh

SCRIPT=$(readlink -f $0)
SCRIPTPATH=$(dirname "$SCRIPT")

GEO=""
RUN=""

# Setup regex expressions
regex_run="run ([0-9]*)"
regex_geo="geo (geo_.*.dat)"

# Loop through input file to set a map with key-value pair (ie. run, geo.dat file)
while read LINE
do
    if [[ $LINE =~ "#" ]]; then
        continue
    fi

    if [[ $LINE =~ $regex_run ]]; then RUN=${BASH_REMATCH[1]}; fi
    if [[ $LINE =~ $regex_geo ]]; then GEO=${BASH_REMATCH[1]}; fi

    if [[ $RUN != "" ]] && [[ $GEO != "" ]]; then
	if [[ "${RUNS[@]}" =~ "${RUN}" ]] || [[ ${#RUNS[@]} == 0 ]]; then
	    GEOFILES[$RUN]=$GEO
	else
	    continue
	fi
    fi
done < $INPUT_FILE


# Dispatch the jobs to the NAF batch system
for RUN in "${!GEOFILES[@]}"; do

    if test "$MODE" = "tele"; then
	$TA -g ${GEOFILES[$RUN]} $RUN &

    elif test "$MODE" = "scope"; then
	$SA -f $INPUT_FILE $RUN &
    fi
done


wait

if [ "$isNAF" = 1 ]; then
    echo "Please check your jobs with qstat -u $USER | grep $MODE"
else
    echo "Processing all nominal samples finished!"
fi
