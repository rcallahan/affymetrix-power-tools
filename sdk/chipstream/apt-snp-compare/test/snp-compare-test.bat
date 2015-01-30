rem
rem check apt-snp-compare.exe
echo testing apt-snp-compare.exe
if exist output.generated.txt del output.generated.txt
apt-snp-compare calls.gold.txt calls.generated.txt 2> output.generated.txt
findstr /c:"Call rate: 100 For SNPs with known type: 100" output.generated.txt
if errorlevel 1 goto error
findstr /c:"Overall Accuracy: 100" output.generated.txt
if errorlevel 1 goto error
findstr /c:"Het Accuracy: 100" output.generated.txt
if errorlevel 1 goto error
findstr /c:"Hom Accuracy: 100" output.generated.txt
if errorlevel 1 goto error
echo apt-snp-compare.exe passes tests
goto end
:error
echo string not found
exit 1
:end
