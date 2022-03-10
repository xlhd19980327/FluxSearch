:: set /p mypath=Code file path:
set /p input=XML input file path:
set /p output=XML output file path:
md %output%\positive
md %output%\negative
:: Rscript %mypath%\envipick.R %input% %output%
Rscript envipick.R %input% %output%
pause