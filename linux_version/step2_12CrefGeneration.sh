#!/usr/bin/bash
read -p "12C input file path: " N
read -p "12C output file path: " Grp
Rscript C12refGeneration_linux.R $N $Grp
