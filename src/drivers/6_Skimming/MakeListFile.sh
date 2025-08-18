#!/bin/bash


basepath="root://eoscms.cern.ch://eos/cms/store/user/amestvir/NanoAOD/"

for entry in `eos ls $1`; do
    for rootfile in `eos ls $1/$entry`; do
	((counter++))
	 echo File: $counter
     xrdcp $basepath/$1/$entry/$rootfile $localbasepath/$Directory
     ##echo $basepath/$1/$entry/$rootfile $localbasepath/$Directory
    done
done
