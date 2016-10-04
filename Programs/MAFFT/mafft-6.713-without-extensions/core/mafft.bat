@echo off

set CYGDIR=C:\cygwin
@REM Specify the root directory of Cygwin.
@REM CYGDIR=C:\cygwin if you took the defaults.

set PATHBK=%PATH%
set PATH=/usr/bin/
%CYGDIR%\bin\bash %CYGDIR%\usr\local\bin\mafft %* 
set PATH=%PATHBK%

@REM %CYGDIR%\bin\sleep 10000

