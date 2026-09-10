@echo off
rem Genetic storage framework (C23) -- clang-only build
rem   build.bat          build everything (vivi_test.exe, vivi_demo.exe)
rem   build.bat lib      framework only -> vivi.lib
rem   build.bat test     build and run the test suite
rem   build.bat density  build and run the packing-density report
rem   build.bat asan     build and run the test suite under ASan
rem   build.bat clean    remove outputs
setlocal
set CC=clang
set AR=llvm-ar
set CFLAGS=-std=c23 -O2 -Wall -Wextra -Iinclude
set LIBSRC=src\vivi.c src\dna.c src\genome.c src\chromosome.c src\cell.c src\organism.c

if "%1"=="clean" goto clean
if "%1"=="lib" goto lib
if "%1"=="test" goto test
if "%1"=="density" goto density
if "%1"=="asan" goto asan
goto all

:lib
%CC% %CFLAGS% -c %LIBSRC%
if errorlevel 1 exit /b 1
%AR% rcs vivi.lib vivi.o dna.o genome.o chromosome.o cell.o organism.o
if errorlevel 1 exit /b 1
del *.o
echo vivi.lib built.
exit /b 0

:test
%CC% %CFLAGS% -o vivi_test.exe test\test.c %LIBSRC%
if errorlevel 1 exit /b 1
vivi_test.exe
exit /b %errorlevel%

:density
%CC% %CFLAGS% -o density.exe test\density.c %LIBSRC%
if errorlevel 1 exit /b 1
density.exe
exit /b %errorlevel%

:asan
rem copy the ASan runtime next to the test binary, where the loader finds it
for /f "delims=" %%i in ('%CC% -print-resource-dir') do set RD=%%i
if exist "%RD%\lib\windows\clang_rt.asan_dynamic-x86_64.dll" copy /y "%RD%\lib\windows\clang_rt.asan_dynamic-x86_64.dll" . >nul
%CC% %CFLAGS% -O1 -g -fsanitize=address -o vivi_test_asan.exe test\test.c %LIBSRC%
if errorlevel 1 exit /b 1
vivi_test_asan.exe
exit /b %errorlevel%

:all
%CC% %CFLAGS% -o vivi_test.exe test\test.c %LIBSRC%
if errorlevel 1 exit /b 1
%CC% %CFLAGS% -o vivi_demo.exe demo\demo.c %LIBSRC%
if errorlevel 1 exit /b 1
echo Build OK: vivi_test.exe vivi_demo.exe
exit /b 0

:clean
del vivi_test.exe vivi_demo.exe density.exe vivi.lib vivi_test_asan.exe vivi_test_asan.ilk vivi_test_asan.pdb *.obj 2>nul
del clang_rt.asan_dynamic-x86_64.dll 2>nul
exit /b 0

