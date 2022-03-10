#!/usr/bin/bash
read -p "enviPick results file path: " DIR
read -p "12C file path: " N
read -p "Result file path: " Grp
Rscript fluxSearch_procedure_linux.R $DIR $N $Grp