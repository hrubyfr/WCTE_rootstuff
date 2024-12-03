#!/bin/bash

#get ntuples
#source copy_tuples.sh

#folder for plots
mkdir plots/run${1}

#analyse
./analyse_tuple ${1}
