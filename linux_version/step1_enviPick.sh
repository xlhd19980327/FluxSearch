#!/usr/bin/bash
read -p "XML input file path: " N
read -p "XML output file path: " Grp
mkdir $Grp/positive
mkdir $Grp/negative
Rscript enviPick_linux.R $N $Grp
