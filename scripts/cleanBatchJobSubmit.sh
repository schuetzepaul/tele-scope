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
	    if [ ! -f $OPTARG ];then
		STRING=`echo $OPTARG | awk -F "," '{for(i=1; i <= NF; i++) printf "%s ",$i}'`
		
	    elif [ -s $OPTARG ];then
		STRING=`cat $OPTARG | awk -F "," '{for(i=1; i <= NF; i++) printf "%s ",$i}'`
	    fi

	    array=(${STRING})

	    for elem in "${array[@]}"
	    do
		if [[ $elem =~ "-" ]] ;then
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


JOBID+=($(find batch_output/${MODE}.o* -maxdepth 1 -size +0 -print | xargs head -n +2 -- | awk '{if(match($0,/==> (.*) <==/,arr))print arr[1];if(match($0,/run (.*)/,arr))print arr[1]}' | xargs -n 2 sh -c 'echo $0 $1'))

for (( i=0; i < ${#JOBID[@]}; i+=2 ));do

	if [[ "${RUNS[@]}" =~ "${JOBID[$i+1]}" ]]; then

	       if test "$MODE" = "tele";then
		   rm run_${JOBID[$i+1]}/{${MODE}${JOBID[$i+1]}.root,{align,hot}_${JOBID[$i+1]}.dat}  ${JOBID[$i]}
		   
	       elif test "$MODE" = "scope";then
		   rm run_${JOBID[$i+1]}/{${MODE}${JOBID[$i+1]}.root,align{MOD,REF,DUT}_${JOBID[$i+1]}.dat}  ${JOBID[$i]}
	       fi
	fi
done

