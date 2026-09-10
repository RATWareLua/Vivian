@echo off
rem Genetic storage framework (C23) -- clang-only build
rem   build.bat          build everything (vivi_test.exe, vivi_demo.exe)
rem   build.bat lib      framework only -> vivi.lib
rem   build.bat test     build and run the test suite
rem   build.bat density  build and run the packing-density report
rem   build.bat sim      build and run a channel experiment
rem   build.bat research run the full parameter sweeps -> research_*.csv
rem   build.bat xcheck   byte-parity check against framework/luau/xcheck.lua
rem   build.bat asan     build and run the test suite under ASan
rem   build.bat fuzz     libFuzzer targets are POSIX-only (run make fuzz)
rem   build.bat clean    remove outputs
setlocal
set CC=clang
set AR=llvm-ar
set CFLAGS=-std=c23 -O2 -Wall -Wextra -Iinclude
set LIBSRC=src\vivi.c src\dna.c src\genome.c src\chromosome.c src\channel.c src\pool.c src\parity.c src\cell.c src\organism.c

if "%1"=="clean" goto clean
if "%1"=="lib" goto lib
if "%1"=="test" goto test
if "%1"=="density" goto density
if "%1"=="sim" goto sim
if "%1"=="research" goto research
if "%1"=="xcheck" goto xcheck
if "%1"=="asan" goto asan
if "%1"=="fuzz" goto fuzz
goto all

:lib
%CC% %CFLAGS% -c %LIBSRC%
if errorlevel 1 exit /b 1
%AR% rcs vivi.lib vivi.o dna.o genome.o chromosome.o channel.o pool.o parity.o cell.o organism.o
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

:sim
%CC% %CFLAGS% -o vivi_sim.exe test\sim.c %LIBSRC%
if errorlevel 1 exit /b 1
vivi_sim.exe --size 256 --trials 1000 --p-sub 0.0005 --header
vivi_sim.exe --size 256 --trials 1000 --p-sub 0.001
vivi_sim.exe --size 256 --trials 1000 --p-sub 0.002
vivi_sim.exe --size 256 --trials 1000 --p-sub 0.004 --p-drop 0.01
vivi_sim.exe --size 1024 --gene-raw 64 --trials 1000 --p-sub 0.001 --header
vivi_sim.exe --size 1024 --gene-raw 64 --access --trials 1000 --p-sub 0.001 --p-access 0.05 --p-cross 0.01
vivi_sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --header
vivi_sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --replicas 3
vivi_sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --parity 4
vivi_sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --replicas 2 --parity 4
exit /b %errorlevel%

:research
%CC% %CFLAGS% -o vivi_sim.exe test\sim.c %LIBSRC%
if errorlevel 1 exit /b 1
vivi_sim.exe --size 1024 --gene-raw 64 --trials 1000 --p-sub 0 --header > research_whole.csv
vivi_sim.exe --size 1024 --gene-raw 64 --trials 1000 --p-sub 0.0002 >> research_whole.csv
vivi_sim.exe --size 1024 --gene-raw 64 --trials 1000 --p-sub 0.0005 >> research_whole.csv
vivi_sim.exe --size 1024 --gene-raw 64 --trials 1000 --p-sub 0.001 >> research_whole.csv
vivi_sim.exe --size 1024 --gene-raw 64 --trials 1000 --p-sub 0.002 >> research_whole.csv
vivi_sim.exe --size 1024 --gene-raw 64 --trials 1000 --p-sub 0.004 >> research_whole.csv
vivi_sim.exe --size 1024 --gene-raw 64 --access --trials 1000 --p-sub 0.001 --p-access 0 --header > research_access.csv
vivi_sim.exe --size 1024 --gene-raw 64 --access --trials 1000 --p-sub 0.001 --p-access 0.01 --p-cross 0.01 >> research_access.csv
vivi_sim.exe --size 1024 --gene-raw 64 --access --trials 1000 --p-sub 0.001 --p-access 0.05 --p-cross 0.01 >> research_access.csv
vivi_sim.exe --size 1024 --gene-raw 64 --access --trials 1000 --p-sub 0.001 --p-access 0.1 --p-cross 0.01 >> research_access.csv
vivi_sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --header > research_library.csv
vivi_sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --replicas 2 >> research_library.csv
vivi_sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --replicas 3 >> research_library.csv
vivi_sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --replicas 4 >> research_library.csv
vivi_sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --parity 2 >> research_library.csv
vivi_sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --parity 4 >> research_library.csv
vivi_sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --parity 8 >> research_library.csv
vivi_sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 --replicas 2 --parity 4 >> research_library.csv
echo wrote research_whole.csv research_access.csv research_library.csv
exit /b 0

:fuzz
echo libFuzzer targets are POSIX-only: run "make fuzz" on Linux/macOS.
exit /b 1

:xcheck
if not exist ..\..\luau.exe (
  echo luau.exe not found at the repo root: get it from
  echo   https://github.com/luau-lang/luau/releases  -- luau-windows.zip
  exit /b 1
)
%CC% %CFLAGS% -o xcheck_c.exe test\xcheck.c %LIBSRC%
if errorlevel 1 exit /b 1
xcheck_c.exe > xcheck_c.txt
if errorlevel 1 exit /b 1
..\..\luau.exe ..\luau\xcheck.lua > xcheck_lua.txt
if errorlevel 1 exit /b 1
fc /b xcheck_c.txt xcheck_lua.txt >nul
if errorlevel 1 (
  echo BYTE PARITY FAILED: xcheck_c.txt vs xcheck_lua.txt
  exit /b 1
)
echo C^<-^>Lua byte parity OK
exit /b 0

:asan
rem copy the ASan runtime next to the test binary, where the loader finds it
for /f "delims=" %%i in ('%CC% -print-resource-dir') do set RD=%%i
if exist "%RD%\lib\windows\clang_rt.asan_dynamic-x86_64.dll" copy /y "%RD%\lib\windows\clang_rt.asan_dynamic-x86_64.dll" . >nul
%CC% %CFLAGS% -O1 -g -fsanitize=address -DVIVI_TEST_TRACK -o vivi_test_asan.exe test\test.c %LIBSRC%
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
del vivi_test.exe vivi_demo.exe density.exe vivi_sim.exe vivi.lib vivi_test_asan.exe vivi_test_asan.ilk vivi_test_asan.pdb research_whole.csv research_access.csv research_library.csv *.obj 2>nul
del xcheck_c.exe xcheck_c.txt xcheck_lua.txt 2>nul
del clang_rt.asan_dynamic-x86_64.dll 2>nul
exit /b 0

