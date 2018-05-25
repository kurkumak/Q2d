#!/bin/bash
#PBS -V
#PBS -q q-1sem
#PBS -M giampaolofolena@gmail.com
#PBS -m bae
#PBS -l cput=168:00:00,mem=4000mb,nodes=1:ppn=1

#inside_directory
cd /home/gfolena/Documents/GIT/Q2d
#on lance
gcc Q2d.c -o Q2d${1}_${2}_${3}_${4} -lm
./Q2d${1}_${2}_${3}_${4} $1 $2 $3 $4
rm Q2d${1}_${2}_${3}_${4}
