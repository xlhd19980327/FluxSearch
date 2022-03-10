:: set /p mypath=Code file path:
set /p input=12C input file path:
set /p output=12C output file path:
:: Rscript %mypath%\12CrefGeneration.R %input% %output%
Rscript 12CrefGeneration.R %input% %output%
pause