#!/bin/bash

if [ $# == 1 ]
	then echo "$1 will be executed as soon as all previous jobs have finished."
fi

bla=$(qstat | grep schuep -c)
echo "${bla} jobs are running"

while [ ${bla} != "0" ]
do
    sleep 10
    bla=$(qstat | grep schuep -c)
    if [ ${bla} != "0" ]
	then
	echo "${bla} jobs are still running"
    else
	echo "${bla} jobs are running - great!"
    fi
done

$1
