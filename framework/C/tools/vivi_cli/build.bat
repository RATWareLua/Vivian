@echo off
setlocal
set CC=clang
set CFLAGS=-std=c23 -O2 -Wall -Wextra -I..\..\include
set LIBSRC=..\..\src\vivi.c ..\..\src\dna.c ..\..\src\genome.c ..\..\src\chromosome.c ..\..\src\channel.c ..\..\src\pool.c ..\..\src\parity.c ..\..\src\cell.c ..\..\src\organism.c
%CC% %CFLAGS% -o vivi_cli.exe vivi_cli.c %LIBSRC%
exit /b %errorlevel%
