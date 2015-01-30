@Echo off
REM build for CSV2SQLite Command Line Utility.
set buildDate=%Date:~10,4%%Date:~4,2%%Date:~7,2%%Time:~0,2%
set IS2010Home=C:\Progra~2\Installshield\StandaloneBuild2010
set MSDEV9Home=C:\Progra~2\Microsoft Visual Studio 9.0
set dotNet3FxPath=C:\WINDOWS\Microsoft.NET\Framework\v3.5
if "%XercesHome%" == "" set XercesHome=C:\Xerces-C

set PublishDir=C:\Publish\CSV2SQLite\CLI
set Installsrc=C:\Publish
set BuildDir=C:\Bamboo\xml-data\build-dir\CSV2SQLITE-CLI

set TEMP=C:\Temp
set TMP=C:\Temp

call "%MSDEV9Home%\Common7\Tools\vsvars32.bat"

REM *---------- Update the build environment ------------*
set LIB=%LIB%;%XercesHome%\Build\Win32\VC9\Release
set LIBPATH=%LIBPATH%;%XercesHome%\Build\Win32\VC9\Release
set INCLUDE=%INCLUDE%;%XercesHome%\src

:BUILD
REM *---------- Build Project solutions --------*
cd /d %BuildDir%\affy\sdk\dbtools\apt-annotation-converter
MSBuild apt-annotation-converter.sln /p:Configuration=Release;VCBuildAdditionalOptions="/useenv"

if ERRORLEVEL == 1  set BuildStatus=FAILED
if "%BuildStatus%" == "FAILED" GOTO BUILDFAILED

copy /y apt-annotation-converter.exe.config %PublishDir%
copy /y Release\*.exe %PublishDir%

goto PUBLISH
REM *---------- Build Installation -------*
SET MergeModulePath=%publishDir%\Installation\MergeModules,%IS2010Home%\Modules\i386,C:\PROGRA~2\COMMON~1\MERGEM~1

cd /d "%BuildDir%"\Applications\GenotypingConsole\CSVtoSQLite\Installation
DEL /Q %publishDir%\Installation\Release\Release\LogFiles\*.*
DEL /Q %publishDir%\Installation\Release\Release\Reports\*.*

%IS2010Home%\System\IsCmdBld.exe -p CSVToSQLite.ism -r "Release" -b %publishDir%\Installation -o %MergeModulePath% -t %dotNet3FxPath%

if ERRORLEVEL == 1 set BuildStatus=FAILED
if "%BuildStatus%" == "FAILED" GOTO BUILDFAILED

:PUBLISH
NET USE \\SWStorage\SWBuild /USER:SWStorage\swbuilder let1build
set RemotePubDir=\\swstorage\swbuild\CSV2SQLite\CLI\CLI_0.0.0.%BuildVersion%

xcopy /s /i %PublishDir% "%RemotePubDir%"
REM xcopy /i /s /q %publishDir%\Installation\Release\Release\DiskImages\DISK1 "%RemotePubDir%"
if ERRORLEVEL == 1 set BuildStatus=FAILED

NET USE \\SWStorage\SWBuild /D

if "%BuildStatus%" == "FAILED" GOTO BUILDFAILED
GOTO END

:BUILDFAILED
Echo ***** BUILD FAILED *****
Exit /B 1

:END
