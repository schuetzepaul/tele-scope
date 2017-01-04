#!/bin/bash

if [ $# != 4 ]
	then echo "Usage: ./extend_runlist.sh runlist from to constantStringBehindRunnumber"
	exit
fi

i=$2
while [ "$i" -le "$3" ]; do 
	echo "$i$4" >> "$1"
	echo "$i$4"
	i=$(($i+1))
done
echo "" >> "$1"

