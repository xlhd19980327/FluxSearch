:: set /p mypath=Code file path:
set /p input=enviPick results files path:
set /p output=12C files path:
set /p result=Result files path:
:: Rscript %mypath%\12CrefGeneration.R %input% %output%
Rscript fluxSearch_procedure.R %input% %output% %result%
pause